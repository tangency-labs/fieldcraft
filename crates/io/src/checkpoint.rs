//! Explicit native CPU-fluid restart files, separate from a saved project's display results.
//!
//! The host grants the path. A compact header identifies the original inputs/settings, followed
//! by streamed solver state and its SHA-256. Failed writes never publish a partial replacement.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use fieldcraft_lbm::checkpoint::Store;
use sha2::{Digest, Sha256};

use crate::{Error, Result, fs};

const MAGIC: &[u8; 8] = b"FCLBM001";
const HEADER: usize = 40;

fn invalid(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}
fn io(e: std::io::Error) -> Error {
    Error::Io(format!("CPU checkpoint: {e}"))
}
fn lbm(e: Error) -> fieldcraft_lbm::Error {
    fieldcraft_lbm::Error::Invalid(e.to_string())
}

/// Original step inputs, actual triangulated obstacles and CPU execution settings. Geometry
/// coordinates retain their bits, so rebuilding a different tessellation refuses a restart.
pub fn identity(definition: [u8; 32], settings: &[u8], obstacles: &[fieldcraft_lbm::Obstacle]) -> [u8; 32] {
    let mut h = Sha256::new();
    h.update(MAGIC);
    h.update(env!("CARGO_PKG_VERSION").as_bytes());
    h.update(definition);
    h.update((settings.len() as u64).to_le_bytes());
    h.update(settings);
    h.update((obstacles.len() as u64).to_le_bytes());
    for o in obstacles {
        h.update((o.body.len() as u64).to_le_bytes());
        h.update(o.body.as_bytes());
        h.update((o.surface.points().len() as u64).to_le_bytes());
        for p in o.surface.points() {
            for v in p {
                h.update(v.to_bits().to_le_bytes());
            }
        }
        h.update((o.surface.triangles().len() as u64).to_le_bytes());
        for t in o.surface.triangles() {
            for v in t {
                h.update(v.to_le_bytes());
            }
        }
        h.update((o.surface.tri_face().len() as u64).to_le_bytes());
        for v in o.surface.tri_face() {
            h.update(v.to_le_bytes());
        }
        h.update((o.faces.len() as u64).to_le_bytes());
        for face in &o.faces {
            h.update((face.len() as u64).to_le_bytes());
            h.update(face.as_bytes());
        }
    }
    h.finalize().into()
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    struct Temp(PathBuf);
    impl Temp {
        fn new() -> std::io::Result<Self> {
            static SEQ: AtomicU64 = AtomicU64::new(0);
            let p = std::env::temp_dir().join(format!("fieldcraft-checkpoint-test-{}-{}", std::process::id(), SEQ.fetch_add(1, Ordering::Relaxed)));
            std::fs::create_dir(&p)?;
            Ok(Self(p))
        }
        fn file(&self) -> PathBuf {
            self.0.join("flow.restart")
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            if let Ok(entries) = std::fs::read_dir(&self.0) {
                for e in entries.flatten() {
                    let _ = std::fs::remove_file(e.path());
                }
            }
            let _ = std::fs::remove_dir(&self.0);
        }
    }
    fn save(store: &mut CpuCheckpoint, bytes: &[u8]) -> fieldcraft_lbm::Result<()> {
        store.save(7, &mut |w| w.write_all(bytes).map_err(|e| fieldcraft_lbm::Error::Invalid(e.to_string())))
    }
    fn load(store: &mut CpuCheckpoint, expected: &[u8]) -> fieldcraft_lbm::Result<bool> {
        store.load(&mut |r| {
            let mut b = vec![0; expected.len()];
            r.read_exact(&mut b).map_err(|e| fieldcraft_lbm::Error::Invalid(e.to_string()))?;
            if b != expected {
                return Err(fieldcraft_lbm::Error::Invalid("wrong payload".into()));
            }
            Ok(7)
        })
    }
    #[test]
    fn atomic_failure_keeps_the_last_valid_checkpoint_and_resume_verifies_every_byte() {
        let temp = Temp::new().unwrap();
        let path = temp.file();
        let id = [3; 32];
        let mut store = CpuCheckpoint::new(path.clone(), id, 11, false, 4096).unwrap();
        save(&mut store, b"previous state").unwrap();
        let previous = std::fs::read(&path).unwrap();
        save(&mut store, b"original state").unwrap();
        let original = std::fs::read(&path).unwrap();
        assert_ne!(original, previous, "successful replacement of an existing native file");
        assert!(
            store
                .save(8, &mut |w| {
                    w.write_all(b"partial").unwrap();
                    Err(fieldcraft_lbm::Error::Invalid("interrupted writer".into()))
                })
                .is_err()
        );
        assert_eq!(std::fs::read(&path).unwrap(), original);
        assert_eq!(std::fs::read_dir(&temp.0).unwrap().count(), 1, "no partial temp remains");
        let mut resumed = CpuCheckpoint::new(path.clone(), id, 19, true, 4096).unwrap();
        assert!(load(&mut resumed, b"original state").unwrap());
        assert_eq!(resumed.resumed_from, Some(7));
        assert_eq!(resumed.last_step, Some(7));
        assert!(!load(&mut resumed, b"original state").unwrap(), "load once");
        assert!(CpuCheckpoint::new(path.clone(), [4; 32], 11, true, 4096).is_err(), "changed identity");
        let mut damaged = original.clone();
        damaged[HEADER] ^= 1;
        std::fs::write(&path, &damaged).unwrap();
        let mut damaged_store = CpuCheckpoint::new(path.clone(), id, 11, true, 4096).unwrap();
        assert!(
            damaged_store
                .load(&mut |r| {
                    let mut b = [0; 14];
                    r.read_exact(&mut b).unwrap();
                    Ok(7)
                })
                .is_err(),
            "checksum"
        );
        std::fs::write(&path, &original[..original.len() - 1]).unwrap();
        assert!(load(&mut CpuCheckpoint::new(path.clone(), id, 11, true, 4096).unwrap(), b"original state").is_err(), "truncated footer");
        let mut trailing = original;
        trailing.push(0);
        std::fs::write(&path, trailing).unwrap();
        assert!(load(&mut CpuCheckpoint::new(path, id, 11, true, 4096).unwrap(), b"original state").is_err(), "trailing data");
    }
    #[test]
    fn byte_cap_and_unrelated_existing_outputs_are_refused() {
        let temp = Temp::new().unwrap();
        let path = temp.file();
        std::fs::write(&path, b"unrelated output").unwrap();
        assert!(CpuCheckpoint::new(path.clone(), [0; 32], 11, false, 4096).is_err());
        std::fs::remove_file(&path).unwrap();
        let mut store = CpuCheckpoint::new(path.clone(), [0; 32], 11, false, 80).unwrap();
        assert!(store.save(1, &mut |w| w.write_all(&[0; 81]).map_err(|e| fieldcraft_lbm::Error::Invalid(e.to_string()))).is_err());
        assert!(!path.exists());
        assert_eq!(std::fs::read_dir(&temp.0).unwrap().count(), 0);
    }
}

/// One explicitly granted destination. Checkpoints are kept after success and cancellation;
/// callers report the exact path and step, and users deliberately remove or replace them.
pub struct CpuCheckpoint {
    path: PathBuf,
    identity: [u8; 32],
    every: u64,
    resume: bool,
    max_bytes: u64,
    pub last_step: Option<u64>,
    pub resumed_from: Option<u64>,
    pub checked_solver_bytes: Option<u64>,
}

impl CpuCheckpoint {
    pub fn new(path: PathBuf, identity: [u8; 32], every: u64, resume: bool, max_bytes: u64) -> Result<Self> {
        if cfg!(target_arch = "wasm32") {
            return Err(Error::Refused("CPU checkpoint/resume currently requires a persistent native filesystem".into()));
        }
        if every == 0 {
            return Err(invalid("checkpoint_every must be a positive step count"));
        }
        if max_bytes < HEADER as u64 + 32 {
            return Err(Error::Limit("checkpoint byte bound is too small".into()));
        }
        if fs::is_file(&path) {
            // A new run cannot overwrite an unrelated output or another model's restart.
            fs::read_stream(&path, max_bytes, |r| check_header(r, &identity))?;
        } else if resume {
            return Err(invalid(format!("checkpoint to resume does not exist: {}", path.display())));
        }
        Ok(Self { path, identity, every, resume, max_bytes, last_step: None, resumed_from: None, checked_solver_bytes: None })
    }
    pub fn path(&self) -> &Path {
        &self.path
    }
}

fn check_header(reader: &mut dyn Read, identity: &[u8; 32]) -> Result<()> {
    let mut bytes = [0u8; HEADER];
    reader.read_exact(&mut bytes).map_err(io)?;
    if bytes.get(..8) != Some(MAGIC.as_slice()) {
        return Err(invalid("not a supported CPU-fluid checkpoint (format FCLBM001)"));
    }
    if bytes.get(8..) != Some(identity.as_slice()) {
        return Err(invalid(
            "checkpoint model, geometry, solver settings or thread count differ; use the original inputs or choose a new checkpoint path",
        ));
    }
    Ok(())
}

struct HashWriter<'a> {
    inner: &'a mut dyn Write,
    hash: Sha256,
    bytes: u64,
    max: u64,
}
impl Write for HashWriter<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let next = self.bytes.saturating_add(bytes.len() as u64);
        if next > self.max {
            return Err(std::io::Error::other("checkpoint exceeds its configured byte bound"));
        }
        self.inner.write_all(bytes)?;
        self.hash.update(bytes);
        self.bytes = next;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}
struct HashReader<'a> {
    inner: &'a mut dyn Read,
    hash: Sha256,
}
impl Read for HashReader<'_> {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(bytes)?;
        self.hash.update(bytes.get(..n).unwrap_or_default());
        Ok(n)
    }
}

impl Store for CpuCheckpoint {
    fn checked_memory(&mut self, bytes: u64) {
        self.checked_solver_bytes = Some(bytes);
    }
    fn every(&self) -> u64 {
        self.every
    }
    fn load(&mut self, restore: &mut dyn FnMut(&mut dyn Read) -> fieldcraft_lbm::Result<u64>) -> fieldcraft_lbm::Result<bool> {
        if !self.resume {
            return Ok(false);
        }
        let step = fs::read_stream(&self.path, self.max_bytes, |r| {
            check_header(r, &self.identity)?;
            let (step, hash) = {
                let mut reader = HashReader { inner: r, hash: Sha256::new() };
                let step = restore(&mut reader).map_err(|e| invalid(e.to_string()))?;
                (step, <[u8; 32]>::from(reader.hash.finalize()))
            };
            let mut saved = [0u8; 32];
            r.read_exact(&mut saved).map_err(io)?;
            if saved != hash {
                return Err(invalid("checkpoint state checksum differs; the file is damaged"));
            }
            let mut tail = [0u8; 1];
            if r.read(&mut tail).map_err(io)? != 0 {
                return Err(invalid("checkpoint has trailing bytes"));
            }
            Ok(step)
        })
        .map_err(lbm)?;
        self.last_step = Some(step);
        self.resumed_from = Some(step);
        self.resume = false;
        Ok(true)
    }
    fn save(&mut self, step: u64, save: &mut dyn FnMut(&mut dyn Write) -> fieldcraft_lbm::Result<()>) -> fieldcraft_lbm::Result<()> {
        fs::write_output_stream(&self.path, |w| {
            w.write_all(MAGIC).map_err(io)?;
            w.write_all(&self.identity).map_err(io)?;
            let hash = {
                let mut writer = HashWriter { inner: w, hash: Sha256::new(), bytes: HEADER as u64, max: self.max_bytes - 32 };
                save(&mut writer).map_err(|e| invalid(e.to_string()))?;
                <[u8; 32]>::from(writer.hash.finalize())
            };
            w.write_all(&hash).map_err(io)?;
            Ok(())
        })
        .map_err(lbm)?;
        self.last_step = Some(step);
        Ok(())
    }
}
