//! Folders an automation client may read from and write to.
//!
//! An agent driving FieldCraft (the MCP server, docs/mcp.md) can be kept to the folders it was
//! given: every file a command reads must lie inside the read root, every file it writes inside
//! the write root. [`Roots::check`] decides, before the file is touched:
//!
//! - The path is made absolute against the working folder, then its links are followed: a file
//!   (or folder) link that leads out of a root is outside it. For a file that does not exist yet,
//!   its nearest existing folder's links are followed and the rest must be plain names, so `..`
//!   below a folder that does not exist is refused, as is a link that leads nowhere (writing
//!   through it would create its target).
//! - The result is compared with the root by whole folder names: `/work2/x` is not inside
//!   `/work`.
//! - Read and write access are separate; a root left out grants none of its access. Without
//!   either root ([`Roots::default`]) nothing is confined. [`Roots::writes_within`] confines
//!   writes alone (a crash replay may read the model's geometry wherever it is, and writes
//!   nothing outside a scratch folder).
//!
//! A path is checked, then opened: another program that swaps a link in between can still win
//! that race, so the roots keep an agent to the files it was given rather than guarding against
//! other software on the computer.

use std::path::{Component, Path, PathBuf};

use crate::{Error, Result};

/// What a command does with a file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Access {
    Read,
    Write,
}

impl Access {
    fn name(self) -> &'static str {
        match self {
            Access::Read => "read",
            Access::Write => "write",
        }
    }
}

/// The folders file access is confined to; the default confines nothing.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Roots {
    /// `None` while unconfined; otherwise what each of read and write access is granted.
    confined: Option<[Grant; 2]>,
}

/// What one kind of access is granted while confined.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Grant {
    /// Not granted.
    None,
    /// Granted anywhere.
    Anywhere,
    /// Granted inside this resolved folder.
    Within(PathBuf),
}

/// `path` as a user reads it, without the Windows verbatim prefix [`std::fs::canonicalize`] adds
/// to every resolved root and path: `\\?\C:\x` is `C:\x`, `\\?\UNC\server\share\x` is
/// `\\server\share\x`. For display only; the root checks keep the canonical form. Any other
/// verbatim form (a volume GUID, a device) has no plain equivalent and is kept, as is a path
/// without the prefix.
pub fn without_verbatim(path: &Path) -> PathBuf {
    let Some(s) = path.to_str() else { return path.to_path_buf() };
    if let Some(unc) = s.strip_prefix(r"\\?\UNC\") {
        return PathBuf::from(format!(r"\\{unc}"));
    }
    match s.strip_prefix(r"\\?\") {
        Some(rest) if matches!(rest.as_bytes(), [d, b':', b'\\', ..] if d.is_ascii_alphabetic()) => PathBuf::from(rest),
        _ => path.to_path_buf(),
    }
}

fn resolve_root(p: &Path, access: Access) -> Result<PathBuf> {
    let shown = p.display();
    let dir = std::fs::canonicalize(p).map_err(|e| Error::Io(format!("the {} root `{shown}`: {e}", access.name())))?;
    if !dir.is_dir() {
        return Err(Error::Io(format!("the {} root `{shown}` is not a folder", access.name())));
    }
    Ok(dir)
}

impl Roots {
    /// Nothing confined.
    pub fn unconfined() -> Roots {
        Roots::default()
    }

    /// Confine reads to `read` and writes to `write`. A root left out grants none of its access.
    /// Each root must be an existing folder; it is resolved (links followed) once, here.
    pub fn new(read: Option<&Path>, write: Option<&Path>) -> Result<Roots> {
        let grant = |p: Option<&Path>, access: Access| -> Result<Grant> {
            Ok(p.map(|p| resolve_root(p, access)).transpose()?.map_or(Grant::None, Grant::Within))
        };
        Ok(Roots { confined: Some([grant(read, Access::Read)?, grant(write, Access::Write)?]) })
    }

    /// Reads anywhere, writes only inside `write` (an existing folder, resolved here).
    pub fn writes_within(write: &Path) -> Result<Roots> {
        Ok(Roots { confined: Some([Grant::Anywhere, Grant::Within(resolve_root(write, Access::Write)?)]) })
    }

    /// Whether any access is confined.
    pub fn is_confined(&self) -> bool {
        self.confined.is_some()
    }

    /// The root of `access`, when confined to one.
    pub fn root(&self, access: Access) -> Option<&Path> {
        match self.confined.as_ref()?.get(access as usize)? {
            Grant::Within(root) => Some(root),
            Grant::None | Grant::Anywhere => None,
        }
    }

    /// The path to use for `access` to `path`: unchanged while unconfined, else the resolved path
    /// (absolute, links followed), which lies inside the root. Refused with an error naming the
    /// root and the path.
    pub fn check(&self, path: &Path, access: Access) -> Result<PathBuf> {
        let Some(roots) = &self.confined else { return Ok(path.to_path_buf()) };
        let shown = path.display();
        let root = match roots.get(access as usize) {
            Some(Grant::Anywhere) => return Ok(path.to_path_buf()),
            Some(Grant::Within(root)) => root,
            Some(Grant::None) | None => {
                return Err(Error::Refused(format!("automation filesystem access is not granted: {} authority is absent: {shown}", access.name())));
            }
        };
        let resolved = resolve(path).map_err(|why| Error::Refused(format!("automation path rejected: {why}: {shown}")))?;
        if !resolved.starts_with(root) {
            return Err(Error::Refused(format!("automation path rejected: outside the {} root {}: {shown}", access.name(), root.display())));
        }
        Ok(resolved)
    }
}

/// `path` absolute, its links followed: wholly when it exists, else up to its nearest existing
/// folder, below which only plain names may follow.
fn resolve(path: &Path) -> std::result::Result<PathBuf, String> {
    if path.as_os_str().is_empty() {
        return Err("empty path".into());
    }
    #[cfg(windows)]
    if let Some(name) = path.file_stem().and_then(|n| n.to_str()) {
        const DEVICES: [&str; 22] = [
            "CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9", "LPT1", "LPT2", "LPT3", "LPT4",
            "LPT5", "LPT6", "LPT7", "LPT8", "LPT9",
        ];
        if DEVICES.iter().any(|d| d.eq_ignore_ascii_case(name.trim_end())) {
            return Err("a device name".into());
        }
    }
    let absolute = std::path::absolute(path).map_err(|e| e.to_string())?;
    if let Ok(p) = std::fs::canonicalize(&absolute) {
        return Ok(p);
    }
    // A link that leads nowhere: writing through it would create its target.
    if std::fs::symlink_metadata(&absolute).is_ok() {
        return Err("a link that leads nowhere".into());
    }
    let mut rest = Vec::new();
    let mut base = absolute.as_path();
    loop {
        let Some(parent) = base.parent() else { return Err("no existing folder above it".into()) };
        let Some(name) = base.file_name() else { return Err("not a plain file name".into()) };
        rest.push(name.to_os_string());
        base = parent;
        if let Ok(dir) = std::fs::canonicalize(base) {
            if !dir.is_dir() {
                return Err(format!("{} is not a folder", base.display()));
            }
            let mut out = dir;
            for name in rest.iter().rev() {
                out.push(name);
            }
            return Ok(out);
        }
        if std::fs::symlink_metadata(base).is_ok() {
            return Err("a link that leads nowhere".into());
        }
        // `..` or `.` below a folder that does not exist cannot be resolved.
        if !matches!(Path::new(base).components().next_back(), Some(Component::Normal(_))) {
            return Err("`..` below a folder that does not exist".into());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Dir(PathBuf);
    impl Dir {
        fn new(tag: &str) -> Dir {
            let d = std::env::temp_dir().join(format!("fieldcraft-roots-{tag}-{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&d);
            std::fs::create_dir_all(&d).unwrap();
            Dir(std::fs::canonicalize(&d).unwrap())
        }
    }
    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn verbatim_prefix_is_taken_off_for_display() {
        let plain = |s: &str| without_verbatim(Path::new(s));
        assert_eq!(plain(r"\\?\C:\Users\Ada\Dev\x"), Path::new(r"C:\Users\Ada\Dev\x"));
        assert_eq!(plain(r"\\?\UNC\server\share\x"), Path::new(r"\\server\share\x"));
        assert_eq!(plain(r"\\?\Volume{1234}\x"), Path::new(r"\\?\Volume{1234}\x"));
        assert_eq!(plain(r"C:\Users\Ada"), Path::new(r"C:\Users\Ada"));
        assert_eq!(plain("/home/ada"), Path::new("/home/ada"));
    }

    #[test]
    fn unconfined_changes_nothing() {
        let r = Roots::unconfined();
        assert!(!r.is_confined());
        assert_eq!(r.check(Path::new("/etc/passwd"), Access::Write).unwrap(), Path::new("/etc/passwd"));
    }

    #[test]
    fn reads_and_writes_stay_inside_their_roots() {
        let d = Dir::new("inside");
        let (read, write) = (d.0.join("in"), d.0.join("in").join("out"));
        std::fs::create_dir_all(&write).unwrap();
        std::fs::write(read.join("a.solidcraft"), b"x").unwrap();
        let r = Roots::new(Some(&read), Some(&write)).unwrap();
        assert_eq!(r.check(&read.join("a.solidcraft"), Access::Read).unwrap(), read.join("a.solidcraft"));
        assert_eq!(r.check(&write.join("new.inp"), Access::Write).unwrap(), write.join("new.inp"));
        // A new file in a new folder below the root: plain names only.
        assert_eq!(r.check(&write.join("x/y/z.inp"), Access::Write).unwrap(), write.join("x/y/z.inp"));
        assert!(r.check(&write.join("x/../../../z.inp"), Access::Write).is_err());
        // `..` through an existing folder resolves first, then is compared.
        assert!(r.check(&write.join("../a.solidcraft"), Access::Write).is_err());
        assert!(r.check(&write.join("../a.solidcraft"), Access::Read).is_ok());
        let e = r.check(&read.join("x.png"), Access::Write).unwrap_err().to_string();
        assert!(e.contains("outside the write root"), "{e}");
        assert!(r.check(Path::new("/etc/passwd"), Access::Read).is_err());
        // Whole folder names: `in2` is not inside `in`.
        std::fs::create_dir_all(d.0.join("in2")).unwrap();
        assert!(r.check(&d.0.join("in2/a"), Access::Read).is_err());
        assert!(r.check(Path::new(""), Access::Read).is_err());
    }

    #[test]
    fn a_root_left_out_grants_nothing() {
        let d = Dir::new("readonly");
        let r = Roots::new(Some(&d.0), None).unwrap();
        assert!(r.is_confined());
        let e = r.check(&d.0.join("x.inp"), Access::Write).unwrap_err().to_string();
        assert!(e.contains("write authority is absent"), "{e}");
        assert!(r.check(&d.0.join("x.inp"), Access::Read).is_ok());
        assert!(Roots::new(Some(&d.0.join("missing")), None).is_err());
        std::fs::write(d.0.join("file"), b"").unwrap();
        assert!(Roots::new(None, Some(&d.0.join("file"))).is_err(), "a root is a folder");
    }

    #[test]
    fn writes_within_reads_anywhere() {
        let d = Dir::new("writes");
        let r = Roots::writes_within(&d.0).unwrap();
        assert!(r.is_confined());
        assert_eq!(r.check(Path::new("/etc/passwd"), Access::Read).unwrap(), Path::new("/etc/passwd"));
        assert!(r.check(Path::new("/etc/passwd"), Access::Write).is_err());
        assert_eq!(r.check(&d.0.join("x.inp"), Access::Write).unwrap(), d.0.join("x.inp"));
        assert_eq!(r.root(Access::Write), Some(d.0.as_path()));
        assert_eq!(r.root(Access::Read), None);
    }

    #[cfg(unix)]
    #[test]
    fn links_are_followed() {
        let d = Dir::new("links");
        let (root, outside) = (d.0.join("root"), d.0.join("outside"));
        std::fs::create_dir_all(&root).unwrap();
        std::fs::create_dir_all(&outside).unwrap();
        std::fs::write(outside.join("secret"), b"x").unwrap();
        std::os::unix::fs::symlink(outside.join("secret"), root.join("file-link")).unwrap();
        std::os::unix::fs::symlink(&outside, root.join("dir-link")).unwrap();
        std::os::unix::fs::symlink(outside.join("nowhere"), root.join("dangling")).unwrap();
        std::os::unix::fs::symlink(&root, outside.join("back")).unwrap();
        let r = Roots::new(Some(&root), Some(&root)).unwrap();
        for p in [root.join("file-link"), root.join("dir-link/secret"), root.join("dir-link/new.inp")] {
            for a in [Access::Read, Access::Write] {
                assert!(r.check(&p, a).is_err(), "{} {a:?}", p.display());
            }
        }
        assert!(r.check(&root.join("dangling"), Access::Write).unwrap_err().to_string().contains("leads nowhere"));
        // A link from outside that leads in is inside.
        assert_eq!(r.check(&outside.join("back/x.inp"), Access::Write).unwrap(), root.join("x.inp"));
    }
}
