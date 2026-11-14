//! How FieldCraft is installed here, which release asset updates it, and the plan that swaps a
//! verified, staged copy in: pure functions of the facts they are given, so every platform's
//! plan is tested on every platform.
//!
//! | Install | Asset | Swap |
//! |---|---|---|
//! | macOS `FieldCraft.app` (Apple silicon) | `fieldcraft-V-macos-aarch64.dmg` | the bundle, then `open -n` it |
//! | macOS `fieldcraft-cli` (Apple silicon) | `fieldcraft-cli-V-macos-aarch64.zip` | the binary |
//! | Windows folder | `fieldcraft-V-windows-x64.zip` | each file (a running `.exe` can be renamed, not overwritten) |
//! | Linux AppImage | `fieldcraft-V-linux-x86_64.AppImage` | the file |
//! | Linux tar.gz folder | `fieldcraft-V-linux-x86_64.tar.gz` | the folder |
//! | `.deb` / `.rpm` (`/usr/bin`) | `fieldcraft_V_amd64.deb` / `fieldcraft-V.x86_64.rpm` | none: the system's package manager installs it |
//! | Flatpak (`/app/bin`) | none | none: `flatpak update` installs it from the Flatpak repository |
//!
//! Each swap moves the installed item into the staging folder and the new one into its place,
//! both renames within one folder's file system, and undoes the moves made when one fails.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use crate::{Error, Result, Version};

/// The macOS app's bundle identifier (packaging/macos/Info.plist.in).
pub const BUNDLE_ID: &str = "io.github.tangencylabs.fieldcraft";
/// The bundle identifiers a downloaded macOS app may carry and still be FieldCraft: [`BUNDLE_ID`]
/// and the one FieldCraft was signed with up to 0.7.x. A 0.7.x app accepts only its own, so it
/// refuses a release with [`BUNDLE_ID`]: those installs download one release by hand
/// (docs/releasing.md).
pub const ACCEPTED_BUNDLE_IDS: [&str; 2] = [BUNDLE_ID, "io.github.jyviko.fieldcraft"];
/// The prefix of a staging folder, beside what it replaces (`.fieldcraft-update-0.5.9`).
pub const STAGE_PREFIX: &str = ".fieldcraft-update-";

/// How this copy of FieldCraft is installed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Install {
    /// The macOS app bundle (`…/FieldCraft.app`).
    MacApp { bundle: PathBuf },
    /// The macOS command-line tool from its own zip.
    MacCli { exe: PathBuf },
    /// The folder the Windows zip was unpacked into (`fieldcraft.exe`, `fieldcraft-cli.exe`, …).
    Windows { dir: PathBuf },
    /// A Linux AppImage file.
    AppImage { file: PathBuf },
    /// The folder the Linux tar.gz was unpacked into (`bin/`, `share/`).
    TarGz { root: PathBuf },
    /// Installed by the Debian package (`/usr/bin`).
    Deb,
    /// Installed by the RPM package (`/usr/bin`).
    Rpm,
    /// The Flatpak (`/app/bin`, read-only): `flatpak update` installs new releases.
    Flatpak,
    /// Not a release package (a development build, an unknown layout): it checks, but does not
    /// install. The reason, in words.
    Unmanaged(String),
}

/// What [`detect`] looks at.
pub struct Facts<'a> {
    /// The running executable (`std::env::current_exe`).
    pub exe: &'a Path,
    /// `APPIMAGE`: the AppImage file the AppImage runtime started this from.
    pub appimage: Option<&'a Path>,
    /// `std::env::consts::OS`.
    pub os: &'a str,
    /// Whether a path exists.
    pub exists: &'a dyn Fn(&Path) -> bool,
}

/// How this process's program is installed.
pub fn here() -> Install {
    let exe = std::env::current_exe().ok().map(|p| std::fs::canonicalize(&p).unwrap_or(p));
    let appimage = std::env::var_os("APPIMAGE").filter(|a| !a.is_empty()).map(PathBuf::from);
    let Some(exe) = exe else { return Install::Unmanaged("the running program's path is unknown".into()) };
    detect(&Facts { exe: &exe, appimage: appimage.as_deref(), os: std::env::consts::OS, exists: &|p: &Path| p.exists() })
}

/// How the program at `facts.exe` is installed.
pub fn detect(f: &Facts) -> Install {
    let dir = f.exe.parent().unwrap_or(Path::new(""));
    let name = f.exe.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let unmanaged = || Install::Unmanaged(format!("{} is not from a FieldCraft release package (a development build?)", f.exe.display()));
    match f.os {
        "macos" => {
            // …/FieldCraft.app/Contents/MacOS/FieldCraft
            let bundle = dir.parent().and_then(Path::parent);
            if let Some(bundle) = bundle.filter(|b| b.extension().is_some_and(|e| e == "app") && dir.ends_with("Contents/MacOS"))
                && (f.exists)(&bundle.join("Contents/Info.plist"))
            {
                return Install::MacApp { bundle: bundle.to_path_buf() };
            }
            if name == "fieldcraft-cli" && !is_build_dir(dir) {
                return Install::MacCli { exe: f.exe.to_path_buf() };
            }
            unmanaged()
        }
        "windows" => {
            let packaged = ["fieldcraft.exe", "fieldcraft-cli.exe", "THIRD-PARTY-LICENSES.txt"].iter().all(|n| (f.exists)(&dir.join(n)));
            if packaged { Install::Windows { dir: dir.to_path_buf() } } else { unmanaged() }
        }
        "linux" => {
            // The Flatpak sandbox (packaging/linux/flatpak): /app is read-only, and the Flatpak
            // repository on Pages delivers each release.
            if dir == Path::new("/app/bin") && (f.exists)(Path::new("/.flatpak-info")) {
                return Install::Flatpak;
            }
            if let Some(file) = f.appimage {
                return Install::AppImage { file: file.to_path_buf() };
            }
            if dir == Path::new("/usr/bin") {
                return if (f.exists)(Path::new("/etc/debian_version")) { Install::Deb } else { Install::Rpm };
            }
            // <root>/bin/fieldcraft with <root>/share/doc/fieldcraft (packaging/linux/package.sh).
            match dir.parent() {
                Some(root) if dir.ends_with("bin") && (f.exists)(&root.join("share/doc/fieldcraft/THIRD-PARTY-LICENSES.txt")) => {
                    Install::TarGz { root: root.to_path_buf() }
                }
                _ => unmanaged(),
            }
        }
        other => Install::Unmanaged(format!("there are no FieldCraft release packages for {other}")),
    }
}

/// A cargo build folder (`target/debug`, `target/release`).
fn is_build_dir(dir: &Path) -> bool {
    dir.components().any(|c| c.as_os_str().to_str().is_some_and(|s| s == "target" || s.starts_with("target-")))
}

impl Install {
    /// The release asset that updates this install to `v` (release builds are macOS aarch64,
    /// Linux x86_64 and Windows x64; `arch` is `std::env::consts::ARCH`). `None` when it cannot
    /// update itself.
    pub fn asset_name(&self, v: &Version, arch: &str) -> Option<String> {
        let x64 = arch == "x86_64";
        let arm64 = arch == "aarch64";
        Some(match self {
            Install::MacApp { .. } if arm64 => format!("fieldcraft-{v}-macos-aarch64.dmg"),
            Install::MacCli { .. } if arm64 => format!("fieldcraft-cli-{v}-macos-aarch64.zip"),
            // 64-bit Arm Windows runs the x64 build.
            Install::Windows { .. } => format!("fieldcraft-{v}-windows-x64.zip"),
            Install::AppImage { .. } if x64 => format!("fieldcraft-{v}-linux-x86_64.AppImage"),
            Install::TarGz { .. } if x64 => format!("fieldcraft-{v}-linux-x86_64.tar.gz"),
            Install::Deb if x64 => format!("fieldcraft_{v}_amd64.deb"),
            Install::Rpm if x64 => format!("fieldcraft-{v}.x86_64.rpm"),
            _ => return None,
        })
    }

    /// The file or folder an update replaces (`None` for system packages and unmanaged builds).
    pub fn target(&self) -> Option<&Path> {
        match self {
            Install::MacApp { bundle } => Some(bundle),
            Install::MacCli { exe } => Some(exe),
            Install::Windows { dir } => Some(dir),
            Install::AppImage { file } => Some(file),
            Install::TarGz { root } => Some(root),
            Install::Deb | Install::Rpm | Install::Flatpak | Install::Unmanaged(_) => None,
        }
    }

    /// Where a download for it is staged: beside what it replaces (one file system, so the swap
    /// is a rename); `downloads` for a system package.
    pub fn stage_dir(&self, v: &Version, downloads: &Path) -> Option<PathBuf> {
        match self {
            Install::Deb | Install::Rpm => Some(downloads.to_path_buf()),
            other => other.target().and_then(Path::parent).map(|p| p.join(format!("{STAGE_PREFIX}{v}"))),
        }
    }

    /// Why it cannot install updates itself, when it cannot.
    pub fn unsupported(&self) -> Option<String> {
        match self {
            Install::Flatpak => Some("installed with Flatpak: `flatpak update` installs new releases".into()),
            Install::Unmanaged(why) => Some(why.clone()),
            _ => None,
        }
    }

    /// In words, for messages.
    pub fn describe(&self) -> String {
        match self {
            Install::MacApp { bundle } => format!("the app at {}", bundle.display()),
            Install::MacCli { exe } => format!("fieldcraft-cli at {}", exe.display()),
            Install::Windows { dir } => format!("the folder {}", dir.display()),
            Install::AppImage { file } => format!("the AppImage {}", file.display()),
            Install::TarGz { root } => format!("the folder {}", root.display()),
            Install::Deb => "the Debian package".into(),
            Install::Rpm => "the RPM package".into(),
            Install::Flatpak => "the Flatpak".into(),
            Install::Unmanaged(why) => why.clone(),
        }
    }
}

/// A verified release, unpacked and ready to swap in ([`crate::stage::stage`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Staged {
    pub install: Install,
    pub version: Version,
    /// The staging folder (removed once the swap is done); for a system package, the folder the
    /// package was downloaded to.
    pub dir: PathBuf,
    /// The new app bundle, AppImage, binary, unpacked folder or package file.
    pub payload: PathBuf,
    /// For a Windows folder: the names of the files the new folder holds.
    pub entries: Vec<String>,
}

/// One item replaced: `target` (when it exists) moves to `backup`, then `new` moves to `target`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Replace {
    pub target: PathBuf,
    pub new: PathBuf,
    pub backup: PathBuf,
}

/// The swap: replace each item in order, start the new version (when asked), then remove the
/// staging folder.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Plan {
    pub replace: Vec<Replace>,
    /// The program (and its arguments) that starts the new version.
    pub launch: Option<(PathBuf, Vec<OsString>)>,
    pub cleanup: PathBuf,
}

/// How a staged release is finished.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Finish {
    /// Swapped in by FieldCraft after the window closes (or at once by the CLI).
    Swap(Plan),
    /// A system package: the downloaded file, opened with the system's installer or installed with
    /// the command given.
    Package { file: PathBuf, command: String },
}

/// The plan that finishes `s`; `relaunch` starts the new app after the swap (the window's
/// update, not the CLI's).
pub fn plan(s: &Staged, relaunch: bool) -> Result<Finish> {
    let one = |target: &Path, backup_name: &str, launch: Option<(PathBuf, Vec<OsString>)>| {
        Finish::Swap(Plan {
            replace: vec![Replace { target: target.to_path_buf(), new: s.payload.clone(), backup: s.dir.join(backup_name) }],
            launch: if relaunch { launch } else { None },
            cleanup: s.dir.clone(),
        })
    };
    Ok(match &s.install {
        Install::MacApp { bundle } => {
            one(bundle, "previous.app", Some((PathBuf::from("/usr/bin/open"), vec!["-n".into(), bundle.as_os_str().to_owned()])))
        }
        Install::MacCli { exe } => one(exe, "previous-fieldcraft-cli", None),
        Install::AppImage { file } => one(file, "previous.AppImage", Some((file.clone(), Vec::new()))),
        Install::TarGz { root } => one(root, "previous", Some((root.join("bin/fieldcraft"), Vec::new()))),
        Install::Windows { dir } => {
            if s.entries.is_empty() {
                return Err(Error::Install("the staged folder is empty".into()));
            }
            let replace =
                s.entries.iter().map(|n| Replace { target: dir.join(n), new: s.payload.join(n), backup: s.dir.join("previous").join(n) }).collect();
            Finish::Swap(Plan { replace, launch: relaunch.then(|| (dir.join("fieldcraft.exe"), Vec::new())), cleanup: s.dir.clone() })
        }
        Install::Deb => Finish::Package { file: s.payload.clone(), command: format!("sudo apt install {}", shell_quote(&s.payload)) },
        Install::Rpm => Finish::Package { file: s.payload.clone(), command: format!("sudo dnf install {}", shell_quote(&s.payload)) },
        Install::Flatpak | Install::Unmanaged(_) => return Err(Error::NotSupported(s.install.unsupported().unwrap_or_default())),
    })
}

/// `path` as one POSIX shell word.
fn shell_quote(path: &Path) -> String {
    let s = path.to_string_lossy();
    if s.bytes().all(|b| b.is_ascii_alphanumeric() || b"/._-+".contains(&b)) { s.into_owned() } else { format!("'{}'", s.replace('\'', r"'\''")) }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn detect_with(os: &str, exe: &str, appimage: Option<&str>, present: &[&str]) -> Install {
        let present: Vec<PathBuf> = present.iter().map(PathBuf::from).collect();
        let exists = |p: &Path| present.iter().any(|q| q == p);
        detect(&Facts { exe: Path::new(exe), appimage: appimage.map(Path::new), os, exists: &exists })
    }

    #[test]
    fn installs_are_recognised_per_platform() {
        let mac = detect_with(
            "macos",
            "/Applications/FieldCraft.app/Contents/MacOS/FieldCraft",
            None,
            &["/Applications/FieldCraft.app/Contents/Info.plist"],
        );
        assert_eq!(mac, Install::MacApp { bundle: "/Applications/FieldCraft.app".into() });
        assert_eq!(detect_with("macos", "/usr/local/bin/fieldcraft-cli", None, &[]), Install::MacCli { exe: "/usr/local/bin/fieldcraft-cli".into() });
        assert!(matches!(detect_with("macos", "/src/fieldcraft/target/release/fieldcraft-cli", None, &[]), Install::Unmanaged(_)));
        assert!(matches!(detect_with("macos", "/src/fieldcraft/target/release/fieldcraft", None, &[]), Install::Unmanaged(_)));

        let win = ["C:/FieldCraft/fieldcraft.exe", "C:/FieldCraft/fieldcraft-cli.exe", "C:/FieldCraft/THIRD-PARTY-LICENSES.txt"];
        assert_eq!(detect_with("windows", "C:/FieldCraft/fieldcraft.exe", None, &win), Install::Windows { dir: "C:/FieldCraft".into() });
        assert!(matches!(detect_with("windows", "C:/src/target/release/fieldcraft.exe", None, &[]), Install::Unmanaged(_)));

        assert_eq!(
            detect_with("linux", "/tmp/.mount_x/usr/bin/fieldcraft", Some("/home/a/FieldCraft.AppImage"), &[]),
            Install::AppImage { file: "/home/a/FieldCraft.AppImage".into() }
        );
        assert_eq!(detect_with("linux", "/usr/bin/fieldcraft", None, &["/etc/debian_version"]), Install::Deb);
        assert_eq!(detect_with("linux", "/usr/bin/fieldcraft", None, &[]), Install::Rpm);
        assert_eq!(
            detect_with("linux", "/opt/fc/bin/fieldcraft", None, &["/opt/fc/share/doc/fieldcraft/THIRD-PARTY-LICENSES.txt"]),
            Install::TarGz { root: "/opt/fc".into() }
        );
        assert!(matches!(detect_with("linux", "/src/target/debug/fieldcraft", None, &[]), Install::Unmanaged(_)));
        // The Flatpak carries share/doc/fieldcraft/THIRD-PARTY-LICENSES.txt like the tar.gz, but /app is read-only.
        let flatpak = detect_with("linux", "/app/bin/fieldcraft", None, &["/.flatpak-info", "/app/share/doc/fieldcraft/THIRD-PARTY-LICENSES.txt"]);
        assert_eq!(flatpak, Install::Flatpak);
        assert!(flatpak.unsupported().is_some_and(|why| why.contains("flatpak update")));
        assert!(matches!(detect_with("freebsd", "/usr/local/bin/fieldcraft", None, &[]), Install::Unmanaged(_)));
    }

    #[test]
    fn each_install_picks_its_release_asset() {
        let v = Version::parse("0.5.9").unwrap();
        let cases = [
            (Install::MacApp { bundle: "/A/FieldCraft.app".into() }, "aarch64", Some("fieldcraft-0.5.9-macos-aarch64.dmg")),
            (Install::MacCli { exe: "/b/fieldcraft-cli".into() }, "aarch64", Some("fieldcraft-cli-0.5.9-macos-aarch64.zip")),
            // Intel Macs have no release build since 0.8.0.
            (Install::MacApp { bundle: "/A/FieldCraft.app".into() }, "x86_64", None),
            (Install::MacCli { exe: "/b/fieldcraft-cli".into() }, "x86_64", None),
            (Install::Windows { dir: "C:/F".into() }, "x86_64", Some("fieldcraft-0.5.9-windows-x64.zip")),
            (Install::Windows { dir: "C:/F".into() }, "aarch64", Some("fieldcraft-0.5.9-windows-x64.zip")),
            (Install::AppImage { file: "/a.AppImage".into() }, "x86_64", Some("fieldcraft-0.5.9-linux-x86_64.AppImage")),
            (Install::AppImage { file: "/a.AppImage".into() }, "aarch64", None),
            (Install::TarGz { root: "/opt/fc".into() }, "x86_64", Some("fieldcraft-0.5.9-linux-x86_64.tar.gz")),
            (Install::Deb, "x86_64", Some("fieldcraft_0.5.9_amd64.deb")),
            (Install::Rpm, "x86_64", Some("fieldcraft-0.5.9.x86_64.rpm")),
            (Install::Unmanaged("dev".into()), "x86_64", None),
        ];
        for (install, arch, want) in cases {
            assert_eq!(install.asset_name(&v, arch).as_deref(), want, "{install:?} {arch}");
        }
        // The names release.yml checks for (its "Check the set" step).
        let workflow = include_str!("../../../.github/workflows/release.yml");
        for name in [
            "fieldcraft-$VERSION-macos-aarch64.dmg",
            "fieldcraft-cli-$VERSION-macos-aarch64.zip",
            "fieldcraft-$VERSION-linux-x86_64.tar.gz",
            "fieldcraft-$VERSION-linux-x86_64.AppImage",
            "fieldcraft_${VERSION}_amd64.deb",
            "fieldcraft-$VERSION.x86_64.rpm",
            "fieldcraft-$VERSION-windows-x64.zip",
        ] {
            assert!(workflow.contains(name), "{name}");
        }
    }

    fn staged(install: Install, payload: &str, entries: &[&str]) -> Staged {
        let v = Version::parse("0.5.9").unwrap();
        let dir = install.stage_dir(&v, Path::new("/dl")).unwrap();
        Staged { install, version: v, dir, payload: payload.into(), entries: entries.iter().map(|e| e.to_string()).collect() }
    }

    #[test]
    fn each_install_has_its_swap_plan() {
        let mac =
            staged(Install::MacApp { bundle: "/Applications/FieldCraft.app".into() }, "/Applications/.fieldcraft-update-0.5.9/FieldCraft.app", &[]);
        let Finish::Swap(p) = plan(&mac, true).unwrap() else { panic!() };
        assert_eq!(
            p.replace,
            [Replace {
                target: "/Applications/FieldCraft.app".into(),
                new: "/Applications/.fieldcraft-update-0.5.9/FieldCraft.app".into(),
                backup: "/Applications/.fieldcraft-update-0.5.9/previous.app".into()
            }]
        );
        assert_eq!(p.launch, Some(("/usr/bin/open".into(), vec!["-n".into(), "/Applications/FieldCraft.app".into()])));
        assert_eq!(p.cleanup, PathBuf::from("/Applications/.fieldcraft-update-0.5.9"));
        let Finish::Swap(p) = plan(&mac, false).unwrap() else { panic!() };
        assert_eq!(p.launch, None);

        let win = staged(
            Install::Windows { dir: "C:/Tools/FieldCraft".into() },
            "C:/Tools/.fieldcraft-update-0.5.9/new/fieldcraft-0.5.9-windows-x64",
            &["fieldcraft.exe", "fieldcraft-cli.exe", "README.md"],
        );
        let Finish::Swap(p) = plan(&win, true).unwrap() else { panic!() };
        assert_eq!(p.replace.len(), 3);
        assert_eq!(p.replace[0].target, PathBuf::from("C:/Tools/FieldCraft/fieldcraft.exe"));
        assert_eq!(p.replace[0].backup, PathBuf::from("C:/Tools/.fieldcraft-update-0.5.9/previous/fieldcraft.exe"));
        assert_eq!(p.launch.map(|l| l.0), Some(PathBuf::from("C:/Tools/FieldCraft/fieldcraft.exe")));
        assert!(plan(&staged(Install::Windows { dir: "C:/F".into() }, "C:/x", &[]), true).is_err());

        let app_image =
            staged(Install::AppImage { file: "/home/a/Apps/FieldCraft.AppImage".into() }, "/home/a/Apps/.fieldcraft-update-0.5.9/new.AppImage", &[]);
        let Finish::Swap(p) = plan(&app_image, true).unwrap() else { panic!() };
        assert_eq!(
            (p.replace[0].target.clone(), p.launch.map(|l| l.0)),
            ("/home/a/Apps/FieldCraft.AppImage".into(), Some("/home/a/Apps/FieldCraft.AppImage".into()))
        );

        let tar = staged(Install::TarGz { root: "/opt/fieldcraft".into() }, "/opt/.fieldcraft-update-0.5.9/new/fieldcraft-0.5.9-linux-x86_64", &[]);
        let Finish::Swap(p) = plan(&tar, true).unwrap() else { panic!() };
        assert_eq!(
            (p.replace[0].backup.clone(), p.launch.map(|l| l.0)),
            ("/opt/.fieldcraft-update-0.5.9/previous".into(), Some("/opt/fieldcraft/bin/fieldcraft".into()))
        );

        let cli =
            staged(Install::MacCli { exe: "/usr/local/bin/fieldcraft-cli".into() }, "/usr/local/bin/.fieldcraft-update-0.5.9/fieldcraft-cli", &[]);
        let Finish::Swap(p) = plan(&cli, true).unwrap() else { panic!() };
        assert_eq!(p.launch, None, "the CLI is not relaunched");

        // System packages are installed by the system's package manager, never replaced.
        let deb = staged(Install::Deb, "/home/a/.local/share/fieldcraft/updates/fieldcraft_0.5.9_amd64.deb", &[]);
        assert_eq!(
            plan(&deb, true).unwrap(),
            Finish::Package {
                file: "/home/a/.local/share/fieldcraft/updates/fieldcraft_0.5.9_amd64.deb".into(),
                command: "sudo apt install /home/a/.local/share/fieldcraft/updates/fieldcraft_0.5.9_amd64.deb".into()
            }
        );
        let rpm = staged(Install::Rpm, "/home/a b/fieldcraft-0.5.9.x86_64.rpm", &[]);
        let Finish::Package { command, .. } = plan(&rpm, true).unwrap() else { panic!() };
        assert_eq!(command, "sudo dnf install '/home/a b/fieldcraft-0.5.9.x86_64.rpm'");
        let dev = Staged {
            install: Install::Unmanaged("dev".into()),
            version: Version::parse("1.0.0").unwrap(),
            dir: "/x".into(),
            payload: "/x/y".into(),
            entries: vec![],
        };
        assert!(plan(&dev, true).is_err());
    }
}
