//! Files kept in memory: the file system of the browser build.
//!
//! A browser page has no file system it can call synchronously, so on `wasm32` [`crate::fs`]
//! keeps files here. The web host puts files in before a command reads them (a file picked or
//! dropped, an autosave kept in the browser's storage) and takes out what commands wrote
//! ([`MemoryFs::take_writes`]) to hand it out as a download or keep it in the browser's storage.
//! Natively the type is only used by its tests.
//!
//! Follows SolidCraft's `crates/io/src/vfs.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::{Error, Result};

/// Most files kept, and most bytes in all: a page that is handed more refuses the write.
pub const MAX_FILES: usize = 10_000;
pub const MAX_TOTAL_BYTES: u64 = 2 << 30;

/// Files by path (shared, so handing one out copies nothing), their total size, and the paths
/// written since the host last took them.
#[derive(Debug, Default)]
pub struct MemoryFs {
    files: BTreeMap<PathBuf, Arc<[u8]>>,
    total: u64,
    written: Vec<PathBuf>,
}

fn missing(path: &Path) -> Error {
    Error::Io(format!("{}: no such file in this browser page", path.display()))
}

impl MemoryFs {
    fn room_for(&self, path: &Path, len: usize) -> Result<()> {
        let replaced = self.files.get(path).map_or(0, |b| b.len() as u64);
        let total = self.total.saturating_sub(replaced).saturating_add(len as u64);
        if !self.files.contains_key(path) && self.files.len() >= MAX_FILES {
            return Err(Error::Limit(format!("{}: more than {MAX_FILES} files in this browser page", path.display())));
        }
        if total > MAX_TOTAL_BYTES {
            return Err(Error::Limit(format!("{}: files in this browser page would take more than {} MB", path.display(), MAX_TOTAL_BYTES >> 20)));
        }
        Ok(())
    }

    fn put(&mut self, path: &Path, bytes: Arc<[u8]>) -> Result<()> {
        self.room_for(path, bytes.len())?;
        self.total = self.total.saturating_add(bytes.len() as u64);
        if let Some(old) = self.files.insert(path.to_path_buf(), bytes) {
            self.total = self.total.saturating_sub(old.len() as u64);
        }
        Ok(())
    }

    /// Put a file in from the host (not recorded as written).
    pub fn insert(&mut self, path: &Path, bytes: Vec<u8>) -> Result<()> {
        self.put(path, bytes.into())
    }

    /// Write a file for a command: kept, and recorded for [`MemoryFs::take_writes`].
    pub fn write(&mut self, path: &Path, bytes: &[u8]) -> Result<()> {
        self.put(path, bytes.into())?;
        self.written.retain(|p| p != path);
        self.written.push(path.to_path_buf());
        Ok(())
    }

    /// A file of at most `max` bytes.
    pub fn read(&self, path: &Path, max: u64) -> Result<Vec<u8>> {
        let bytes = self.files.get(path).ok_or_else(|| missing(path))?;
        if bytes.len() as u64 > max {
            return Err(Error::Limit(format!("{}: {} MB (at most {} MB)", path.display(), bytes.len() >> 20, max >> 20)));
        }
        Ok(bytes.to_vec())
    }

    pub fn is_file(&self, path: &Path) -> bool {
        self.files.contains_key(path)
    }

    pub fn remove(&mut self, path: &Path) -> Result<()> {
        let old = self.files.remove(path).ok_or_else(|| missing(path))?;
        self.total = self.total.saturating_sub(old.len() as u64);
        Ok(())
    }

    /// Whether any file lies under `dir`, at any depth. Paths order by component, so everything
    /// under `dir` follows it at once (a file at `dir` itself is not under it).
    pub fn holds_any(&self, dir: &Path) -> bool {
        use std::ops::Bound::{Excluded, Unbounded};
        self.files.range::<Path, _>((Excluded(dir), Unbounded)).next().is_some_and(|(p, _)| p.starts_with(dir))
    }

    /// The files directly in `dir`, in path order.
    pub fn list(&self, dir: &Path) -> Vec<PathBuf> {
        // Paths order by component, so everything under `dir` follows it.
        self.files
            .range(dir.to_path_buf()..)
            .map(|(p, _)| p)
            .take_while(|p| p.starts_with(dir))
            .filter(|p| p.parent() == Some(dir))
            .cloned()
            .collect()
    }

    /// The files written since the last call, in the order last written, with their content
    /// (a file removed since is left out).
    pub fn take_writes(&mut self) -> Vec<(PathBuf, Arc<[u8]>)> {
        std::mem::take(&mut self.written).into_iter().filter_map(|p| self.files.get(&p).cloned().map(|b| (p, b))).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn files_are_kept_listed_and_handed_out_once_written() {
        let mut m = MemoryFs::default();
        let (a, b) = (Path::new("/browser/a.fieldcraft"), Path::new("/browser/b.inp"));
        m.insert(a, b"picked".to_vec()).unwrap();
        assert_eq!(m.read(a, 100).unwrap(), b"picked");
        assert!(m.take_writes().is_empty(), "a file put in by the host is not a write");
        m.write(b, b"one").unwrap();
        m.write(a, b"saved").unwrap();
        m.write(b, b"two").unwrap();
        let taken: Vec<(PathBuf, Vec<u8>)> = m.take_writes().into_iter().map(|(p, b)| (p, b.to_vec())).collect();
        assert_eq!(taken, vec![(a.to_path_buf(), b"saved".to_vec()), (b.to_path_buf(), b"two".to_vec())]);
        assert!(m.take_writes().is_empty());
        assert_eq!(m.list(Path::new("/browser")), vec![a.to_path_buf(), b.to_path_buf()]);
        assert!(m.list(Path::new("/")).is_empty(), "only the folder's own files");
        m.insert(Path::new("/browser2/c"), vec![]).unwrap();
        m.insert(Path::new("/browser/sub/d"), vec![]).unwrap();
        assert_eq!(m.list(Path::new("/browser")).len(), 2, "not a sibling folder's files, nor a subfolder's");
        m.write(a, b"again").unwrap();
        m.remove(a).unwrap();
        assert!(m.take_writes().is_empty(), "a removed file is not handed out");
        assert!(!m.is_file(a) && m.remove(a).is_err() && m.read(a, 100).is_err());
    }

    #[test]
    fn a_folder_holds_what_lies_under_it_at_any_depth() {
        let mut m = MemoryFs::default();
        m.insert(Path::new("/case2/x"), vec![]).unwrap();
        assert!(!m.holds_any(Path::new("/case")), "a sibling folder's file is not under it");
        m.insert(Path::new("/case/system/controlDict"), vec![]).unwrap();
        assert!(m.holds_any(Path::new("/case")), "a file in a subfolder");
        assert!(!m.holds_any(Path::new("/case/constant")));
        m.insert(Path::new("/b"), vec![]).unwrap();
        m.insert(Path::new("/b/c"), vec![]).unwrap();
        assert!(m.holds_any(Path::new("/b")), "a file at the folder's own path does not hide what is under it");
    }

    #[test]
    fn reads_and_sizes_are_capped() {
        let mut m = MemoryFs::default();
        let p = Path::new("/x");
        m.insert(p, vec![0; 100]).unwrap();
        assert!(matches!(m.read(p, 99), Err(Error::Limit(_))));
        assert_eq!(m.read(p, 100).map(|b| b.len()), Ok(100));
        for i in 1..MAX_FILES {
            m.insert(&PathBuf::from(format!("/f{i}")), vec![]).unwrap();
        }
        assert!(matches!(m.write(Path::new("/one-too-many"), b""), Err(Error::Limit(_))));
        assert!(m.write(p, b"replacing is fine").is_ok());
        assert_eq!(m.total, 17, "the total follows replacements");
        m.remove(p).unwrap();
        assert_eq!(m.total, 0);
    }
}
