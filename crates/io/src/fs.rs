//! The file system edge: crash-safe writes and capped reads.
//!
//! A save writes the new content to a temporary file next to the target, flushes it to disk and
//! replaces the target with one rename, then flushes the directory. A crash or power cut at any
//! moment leaves the old file or the new one, never a torn mix. The file being replaced can be
//! kept as `<file>.bak` (itself written the same way).
//!
//! A file that is not positively a model of the format this version writes is never replaced by
//! a model without a lasting copy: [`save_model`] first copies it beside it, to
//! `<stem>.v<N>.fieldcraft` for a model of format N, `<stem>.backup.<ext>` for anything else (a
//! name not taken yet; `.bak` is replaced by the next save, this copy never is).
//!
//! Adapted from SolidCraft's `crates/io/src/safe_write.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.
//!
//! In the browser (`wasm32`) the same calls read and write files kept in memory
//! ([`crate::memory`]): the web host puts files in ([`insert`]) and takes out what commands wrote
//! ([`take_writes`]). Everything else in FieldCraft reaches files only through this module.

#[cfg(target_arch = "wasm32")]
pub use browser::*;
#[cfg(not(target_arch = "wasm32"))]
pub use disk::*;

use std::path::{Path, PathBuf};

/// The backup of a file written with `backup: true`.
pub fn backup_path(path: &Path) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".bak");
    PathBuf::from(s)
}

/// Most names tried for a format backup before a save gives up.
const MAX_VERSION_BACKUPS: u32 = 1000;

/// The names a copy of the file at `path` may take before a model replaces it, in order:
/// `<stem>.v<N>.<ext>` for a model of format N, `<stem>.backup.<ext>` for any other file, then
/// the same with `-2`, `-3`, … (`ext` is the file's own extension, `fieldcraft` when it has none).
pub fn version_backup_paths(path: &Path, format: Option<u32>) -> impl Iterator<Item = PathBuf> + '_ {
    let stem = path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "model".into());
    let ext = path.extension().map(|e| e.to_string_lossy().into_owned()).unwrap_or_else(|| crate::EXTENSION.into());
    let tag = format.map_or_else(|| "backup".to_string(), |f| format!("v{f}"));
    (1..=MAX_VERSION_BACKUPS).map(move |k| {
        let suffix = if k == 1 { String::new() } else { format!("-{k}") };
        path.with_file_name(format!("{stem}.{tag}{suffix}.{ext}"))
    })
}

/// What the existing file at `path` is when it is a FieldCraft model or a SolidCraft design,
/// which no export, deck or picture may replace: by its extension (`.fieldcraft`,
/// `.solidcraft`, any case) or by its content (a FieldCraft package with a format tag, a
/// SolidCraft design's `"format": "solidcraft/…"`). `None` when there is no file or it is
/// neither.
pub fn protected_kind(path: &Path) -> Option<&'static str> {
    if !is_file(path) {
        return None;
    }
    let ext = path.extension().map(|e| e.to_string_lossy().to_ascii_lowercase());
    if ext.as_deref() == Some(crate::EXTENSION) || format_of_file(path).is_some() {
        return Some("FieldCraft model");
    }
    let head = read_prefix(path, 4096);
    let text = String::from_utf8_lossy(&head);
    let design = text.trim_start().starts_with('{') && text.contains("\"solidcraft/");
    (ext.as_deref() == Some("solidcraft") || design).then_some("SolidCraft design")
}

/// Write an export, a deck, a picture or any other output that is not the model's own file:
/// atomically, with no backup, but never over a FieldCraft model or SolidCraft design
/// ([`protected_kind`]), which that would destroy with nothing kept.
pub fn write_output(path: &Path, bytes: &[u8]) -> crate::Result<()> {
    if let Some(kind) = protected_kind(path) {
        return Err(crate::Error::Refused(format!("{} is a {kind}; an export never replaces one (choose another file name)", path.display())));
    }
    write_atomic(path, bytes, false)
}

#[cfg(target_arch = "wasm32")]
mod browser {
    use std::cell::RefCell;
    use std::path::{Path, PathBuf};

    use crate::Result;
    use crate::memory::MemoryFs;

    /// Streaming restart files require a persistent native filesystem in this first slice.
    pub fn read_stream<T>(_path: &Path, _max: u64, _read: impl FnOnce(&mut dyn std::io::Read) -> Result<T>) -> Result<T> {
        Err(crate::Error::Refused("streaming CPU checkpoints currently require the native filesystem".into()))
    }

    pub fn write_output_stream<T>(_path: &Path, _write: impl FnOnce(&mut dyn std::io::Write) -> Result<T>) -> Result<T> {
        Err(crate::Error::Refused("streaming CPU checkpoints currently require the native filesystem".into()))
    }

    thread_local! {
        // The browser runs the app on one thread.
        static FILES: RefCell<MemoryFs> = RefCell::new(MemoryFs::default());
    }

    /// Write `bytes` to `path` (a backup is not kept: the host hands each save out as a
    /// download).
    pub fn write_atomic(path: &Path, bytes: &[u8], _backup: bool) -> Result<()> {
        FILES.with(|f| f.borrow_mut().write(path, bytes))
    }

    /// Save a model file. The host hands each save out as a download, which replaces nothing,
    /// so no copy is kept.
    pub fn save_model(path: &Path, bytes: &[u8]) -> Result<Option<PathBuf>> {
        write_atomic(path, bytes, false).map(|()| None)
    }

    /// The format of the model file at `path` (see [`crate::format_of`]).
    pub fn format_of_file(path: &Path) -> Option<u32> {
        FILES.with(|f| f.borrow().read(path, crate::MAX_FILE_BYTES).ok()).and_then(|b| crate::format_of(std::io::Cursor::new(b)))
    }

    /// Read a file of at most `max` bytes.
    pub fn read_capped(path: &Path, max: u64) -> Result<Vec<u8>> {
        FILES.with(|f| f.borrow().read(path, max))
    }

    /// The first `n` bytes of a file (fewer when it is shorter; none when it cannot be read).
    pub fn read_prefix(path: &Path, n: usize) -> Vec<u8> {
        FILES.with(|f| f.borrow().read(path, crate::MAX_FILE_BYTES).map(|b| b.get(..n.min(b.len())).unwrap_or_default().to_vec()).unwrap_or_default())
    }

    pub fn is_file(path: &Path) -> bool {
        FILES.with(|f| f.borrow().is_file(path))
    }

    pub fn remove(path: &Path) -> Result<()> {
        FILES.with(|f| f.borrow_mut().remove(path))
    }

    /// The files directly in `dir`.
    pub fn list(dir: &Path) -> Vec<PathBuf> {
        FILES.with(|f| f.borrow().list(dir))
    }

    /// Files in memory keep no modification time.
    pub fn modified(_path: &Path) -> Option<u64> {
        None
    }

    /// Folders are implied by the files' paths.
    pub fn create_dir_all(_dir: &Path) -> Result<()> {
        Ok(())
    }

    /// The host puts a file in (one picked, dropped or kept in the browser's storage).
    pub fn insert(path: &Path, bytes: Vec<u8>) -> Result<()> {
        FILES.with(|f| f.borrow_mut().insert(path, bytes))
    }

    /// What commands wrote since the last call, for the host to hand out or store.
    pub fn take_writes() -> Vec<(PathBuf, std::sync::Arc<[u8]>)> {
        FILES.with(|f| f.borrow_mut().take_writes())
    }
}

#[cfg(not(target_arch = "wasm32"))]
mod disk {
    use std::fs::{self, OpenOptions};
    use std::io::{Read, Write};
    use std::path::{Path, PathBuf};

    use super::{backup_path, version_backup_paths};
    use crate::{Error, FORMAT_VERSION, Result};

    fn io_error(path: &Path, e: std::io::Error) -> Error {
        Error::Io(format!("{}: {e}", path.display()))
    }

    /// A capped streaming input, avoiding a whole-file allocation. The caller checks its schema,
    /// payload lengths and end marker; metadata is checked before any payload is read.
    pub fn read_stream<T>(path: &Path, max: u64, read: impl FnOnce(&mut dyn Read) -> Result<T>) -> Result<T> {
        // Opening a FIFO can block before an opened-file metadata check on Unix.
        let before = fs::metadata(path).map_err(|e| io_error(path, e))?;
        if !before.is_file() || before.len() > max {
            return Err(Error::Limit(format!("{} is not a regular file of at most {max} bytes", path.display())));
        }
        let file = fs::File::open(path).map_err(|e| io_error(path, e))?;
        let meta = file.metadata().map_err(|e| io_error(path, e))?;
        if !meta.is_file() || meta.len() > max {
            return Err(Error::Limit(format!("{} is not a regular file of at most {max} bytes", path.display())));
        }
        read(&mut std::io::BufReader::with_capacity(64 * 1024, file.take(max)))
    }

    /// Write a non-model output using a bounded streaming writer. A failed callback, flush or
    /// rename leaves the previous target intact and removes only this write's temporary file.
    pub fn write_output_stream<T>(path: &Path, write: impl FnOnce(&mut dyn Write) -> Result<T>) -> Result<T> {
        if fs::metadata(path).is_ok_and(|m| !m.is_file()) {
            return Err(Error::Refused(format!("{} is not a regular output file", path.display())));
        }
        if let Some(kind) = super::protected_kind(path) {
            return Err(Error::Refused(format!("{} is a {kind}; checkpoint output never replaces one", path.display())));
        }
        let tmp = temp_beside(path, "new");
        let result = (|| {
            let mut file = OpenOptions::new().write(true).create_new(true).open(&tmp).map_err(|e| io_error(&tmp, e))?;
            let value = {
                let mut buffered = std::io::BufWriter::with_capacity(64 * 1024, &mut file);
                let value = write(&mut buffered)?;
                buffered.flush().map_err(|e| io_error(&tmp, e))?;
                value
            };
            file.sync_all().map_err(|e| io_error(&tmp, e))?;
            drop(file);
            fs::rename(&tmp, path).map_err(|e| io_error(path, e))?;
            sync_dir(path);
            Ok(value)
        })();
        if result.is_err() {
            let _ = fs::remove_file(&tmp);
        }
        result
    }

    fn temp_beside(path: &Path, tag: &str) -> PathBuf {
        let name = path.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_else(|| "file".into());
        static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let seq = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        path.with_file_name(format!(".{name}.{tag}.{}.{seq}.tmp", std::process::id()))
    }

    fn write_synced(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
        let mut f = OpenOptions::new().write(true).create_new(true).open(path)?;
        f.write_all(bytes)?;
        f.sync_all()
    }

    /// Flush the directory entry (the rename) to disk where the platform allows it.
    fn sync_dir(path: &Path) {
        #[cfg(unix)]
        if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty())
            && let Ok(d) = fs::File::open(dir)
        {
            let _ = d.sync_all();
        }
        #[cfg(not(unix))]
        let _ = path;
    }

    /// Steps of a save, so tests can stop one part way (as a crash would).
    #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
    enum Step {
        TempWritten,
        BackupWritten,
        Replaced,
    }

    fn write_steps(path: &Path, bytes: &[u8], backup: bool, stop: Step) -> std::io::Result<()> {
        let tmp = temp_beside(path, "new");
        if let Err(e) = write_synced(&tmp, bytes) {
            let _ = fs::remove_file(&tmp);
            return Err(e);
        }
        if stop == Step::TempWritten {
            return Ok(());
        }
        if backup && path.is_file() {
            // Copy (not move) the old file, so the target exists at every moment.
            let bak_tmp = temp_beside(path, "bak");
            let r = fs::copy(path, &bak_tmp)
                .and_then(|_| OpenOptions::new().write(true).open(&bak_tmp))
                .and_then(|f| f.sync_all())
                .and_then(|_| fs::rename(&bak_tmp, backup_path(path)));
            if let Err(e) = r {
                let _ = fs::remove_file(&bak_tmp);
                let _ = fs::remove_file(&tmp);
                return Err(e);
            }
        }
        if stop == Step::BackupWritten {
            return Ok(());
        }
        if let Err(e) = fs::rename(&tmp, path) {
            let _ = fs::remove_file(&tmp);
            return Err(e);
        }
        sync_dir(path);
        // Temporary files of saves that crashed part way (a minute old or more) go.
        for t in stale_temps(path) {
            if fs::metadata(&t).and_then(|m| m.modified()).ok().and_then(|m| m.elapsed().ok()).is_some_and(|age| age.as_secs() >= 60) {
                let _ = fs::remove_file(t);
            }
        }
        Ok(())
    }

    /// Write `bytes` to `path` atomically; with `backup`, the file being replaced is kept as
    /// `<path>.bak` first.
    pub fn write_atomic(path: &Path, bytes: &[u8], backup: bool) -> Result<()> {
        write_steps(path, bytes, backup, Step::Replaced).map_err(|e| io_error(path, e))
    }

    /// The format of the model file at `path`, from its manifest alone (see
    /// [`crate::format_of`]); `None` when there is no such file or it is not a FieldCraft model.
    pub fn format_of_file(path: &Path) -> Option<u32> {
        let f = fs::File::open(path).ok().filter(|f| f.metadata().is_ok_and(|m| m.is_file() && m.len() <= crate::MAX_FILE_BYTES))?;
        crate::format_of(std::io::BufReader::new(f))
    }

    /// Save a model file atomically, the old file kept as `<path>.bak`. Unless there is no file
    /// at `path` or it is positively a model of the format this version writes, the file is first
    /// copied to the first free name of [`version_backup_paths`], which no later save replaces:
    /// an older model opened and migrated, a later one, a damaged one, one whose manifest this
    /// version cannot read, or any other file a Save As replaces. The save is refused when that
    /// copy cannot be made. Returns the copy's path.
    pub fn save_model(path: &Path, bytes: &[u8]) -> Result<Option<PathBuf>> {
        let format = format_of_file(path);
        let kept = if path.exists() && format != Some(FORMAT_VERSION) { Some(keep_copy(path, format)?) } else { None };
        write_atomic(path, bytes, true)?;
        Ok(kept)
    }

    /// Copy `path` (a model of `format`, or another file) to the first of its backup names not
    /// taken.
    fn keep_copy(path: &Path, format: Option<u32>) -> Result<PathBuf> {
        let what = format.map_or_else(|| "file".to_string(), |f| format!("format {f} file"));
        let refused =
            |e: std::io::Error| Error::Io(format!("{}: not saved, because a copy of the {what} could not be kept first: {e}", path.display()));
        for target in version_backup_paths(path, format) {
            match copy_new(path, &target) {
                Ok(()) => return Ok(target),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(refused(e)),
            }
        }
        Err(refused(std::io::Error::other("every backup name is taken")))
    }

    /// Copy `from` to `to`, which must not exist yet (never replaced), flushed to disk. A copy
    /// that fails part way is removed.
    fn copy_new(from: &Path, to: &Path) -> std::io::Result<()> {
        let mut out = OpenOptions::new().write(true).create_new(true).open(to)?;
        let copied = fs::File::open(from).and_then(|mut src| std::io::copy(&mut src, &mut out)).and_then(|_| out.sync_all());
        if let Err(e) = copied {
            drop(out);
            let _ = fs::remove_file(to);
            return Err(e);
        }
        sync_dir(to);
        Ok(())
    }

    /// Leftover temporary files of saves that never finished (a crash part way), next to `path`.
    pub fn stale_temps(path: &Path) -> Vec<PathBuf> {
        let Some(name) = path.file_name().map(|n| n.to_string_lossy().to_string()) else { return Vec::new() };
        let dir = path.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
        let prefix = format!(".{name}.");
        fs::read_dir(dir)
            .into_iter()
            .flatten()
            .flatten()
            .map(|e| e.path())
            .filter(|p| p.file_name().map(|n| n.to_string_lossy().to_string()).is_some_and(|n| n.starts_with(&prefix) && n.ends_with(".tmp")))
            .take(10_000)
            .collect()
    }

    /// Read a regular file of at most `max` bytes.
    pub fn read_capped(path: &Path, max: u64) -> Result<Vec<u8>> {
        let meta = fs::metadata(path).map_err(|e| io_error(path, e))?;
        if !meta.is_file() {
            return Err(Error::Io(format!("{}: not a regular file", path.display())));
        }
        if meta.len() > max {
            return Err(Error::Limit(format!("{}: {} MB (at most {} MB)", path.display(), meta.len() >> 20, max >> 20)));
        }
        let f = fs::File::open(path).map_err(|e| io_error(path, e))?;
        let mut bytes = Vec::with_capacity(usize::try_from(meta.len()).unwrap_or(0));
        // The file may grow between the size check and the read.
        f.take(max.saturating_add(1)).read_to_end(&mut bytes).map_err(|e| io_error(path, e))?;
        if bytes.len() as u64 > max {
            return Err(Error::Limit(format!("{}: larger than {} MB", path.display(), max >> 20)));
        }
        Ok(bytes)
    }

    pub fn is_file(path: &Path) -> bool {
        path.is_file()
    }

    /// The first `n` bytes of a file (fewer when it is shorter; none when it cannot be read).
    pub fn read_prefix(path: &Path, n: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(n);
        let _ = fs::File::open(path).and_then(|f| f.take(n as u64).read_to_end(&mut out));
        out
    }

    pub fn remove(path: &Path) -> Result<()> {
        fs::remove_file(path).map_err(|e| io_error(path, e))
    }

    /// The files directly in `dir` (at most a million).
    pub fn list(dir: &Path) -> Vec<PathBuf> {
        fs::read_dir(dir)
            .into_iter()
            .flatten()
            .flatten()
            .take(1_000_000)
            .filter(|e| e.file_type().is_ok_and(|t| t.is_file()))
            .map(|e| e.path())
            .collect()
    }

    pub fn create_dir_all(dir: &Path) -> Result<()> {
        fs::create_dir_all(dir).map_err(|e| io_error(dir, e))
    }

    /// When a regular file was last modified (seconds since 1970), if the platform says.
    pub fn modified(path: &Path) -> Option<u64> {
        let meta = fs::metadata(path).ok().filter(fs::Metadata::is_file)?;
        meta.modified().ok()?.duration_since(std::time::UNIX_EPOCH).ok().map(|d| d.as_secs())
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn dir(tag: &str) -> PathBuf {
            let d = std::env::temp_dir().join(format!("fieldcraft-fs-{tag}-{}", std::process::id()));
            let _ = fs::remove_dir_all(&d);
            fs::create_dir_all(&d).unwrap();
            d
        }

        #[test]
        fn replaces_and_keeps_a_backup() {
            let d = dir("bak");
            let p = d.join("a.fieldcraft");
            write_atomic(&p, b"one", true).unwrap();
            assert!(!backup_path(&p).exists(), "nothing to back up the first time");
            write_atomic(&p, b"two", true).unwrap();
            assert_eq!(fs::read(&p).unwrap(), b"two");
            assert_eq!(fs::read(backup_path(&p)).unwrap(), b"one");
            write_atomic(&p, b"three", false).unwrap();
            assert_eq!(fs::read(backup_path(&p)).unwrap(), b"one");
            assert!(stale_temps(&p).is_empty());
            let _ = fs::remove_dir_all(&d);
        }

        /// A save stopped after any step leaves the old file whole (or the new one).
        #[test]
        fn a_save_stopped_part_way_keeps_the_old_file() {
            let d = dir("stop");
            let p = d.join("a.fieldcraft");
            write_atomic(&p, b"old", true).unwrap();
            for stop in [Step::TempWritten, Step::BackupWritten, Step::Replaced] {
                write_steps(&p, b"new content", true, stop).unwrap();
                let now = fs::read(&p).unwrap();
                if stop < Step::Replaced {
                    assert_eq!(now, b"old", "{stop:?}");
                    assert!(!stale_temps(&p).is_empty(), "the unfinished save's temporary file is left over");
                    for t in stale_temps(&p) {
                        fs::remove_file(t).unwrap();
                    }
                } else {
                    assert_eq!(now, b"new content");
                }
                write_atomic(&p, b"old", false).unwrap();
            }
            let _ = fs::remove_dir_all(&d);
        }

        #[test]
        fn a_failed_write_leaves_the_target() {
            let d = dir("fail");
            let p = d.join("missing-dir").join("a.fieldcraft");
            assert!(write_atomic(&p, b"x", true).is_err());
            assert!(!p.exists());
            let _ = fs::remove_dir_all(&d);
        }

        #[test]
        fn reads_are_capped() {
            let d = dir("read");
            let p = d.join("a.bin");
            fs::write(&p, [7u8; 100]).unwrap();
            assert_eq!(read_capped(&p, 100).unwrap().len(), 100);
            assert!(matches!(read_capped(&p, 99), Err(Error::Limit(_))));
            assert!(read_capped(&d, 1 << 20).is_err());
            assert!(read_capped(&d.join("missing"), 1 << 20).is_err());
            let _ = fs::remove_dir_all(&d);
        }
    }
}
