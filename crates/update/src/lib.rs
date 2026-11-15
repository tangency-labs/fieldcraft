//! FieldCraft's updates (docs/updates.md).
//!
//! - **Check** ([`net::Client::check`]): ask the GitHub Releases API for the latest published
//!   release (a nightly build also for the published nightlies; [`Check::choose`]) and compare
//!   its [`Version`] with the running one. The desktop app checks at start-up at most
//!   once a day ([`due`]) when the user leaves automatic checks on, and from Help › Check for
//!   Updates; `app.update.check` and `fieldcraft-cli update --check` run the same check.
//! - **Verify** ([`verify`]): a release's `SHA256SUMS.txt` must carry a signature by the release
//!   key ([`PUBLIC_KEY`], Ed25519 in minisign's format), and the downloaded asset must have the
//!   SHA-256 that list gives it. Anything else is refused; an unsigned release can only be
//!   downloaded by hand from its page.
//! - **Install** ([`install`], [`stage`]): how this copy is installed decides the asset and the
//!   swap. The verified release is unpacked beside the install and swapped in by renames
//!   (after the window closes, for the app), or, for `.deb`/`.rpm`, handed to the system's
//!   package manager.
//!
//! Downloading and unpacking are native only: the browser build has the pure parts (versions,
//! releases, plans), and the web app updates through its service worker.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]

pub mod install;
mod release;
pub mod verify;
mod version;

#[cfg(all(any(test, feature = "mock"), not(target_arch = "wasm32")))]
pub mod mock;
#[cfg(not(target_arch = "wasm32"))]
pub mod net;
#[cfg(not(target_arch = "wasm32"))]
pub mod stage;

pub use install::{Finish, Install, Plan, Staged};
pub use release::{Asset, Check, Release};
pub use version::{NIGHTLY, Version};

/// The repository releases come from.
pub const REPO: &str = "tangency-labs/fieldcraft";
/// The GitHub API.
pub const DEFAULT_HOST: &str = "https://api.github.com";
/// Replaces [`DEFAULT_HOST`] (and allows plain `http`), for tests against a local server only.
pub const HOST_ENV: &str = "FIELDCRAFT_UPDATE_HOST";
/// The public half of the release signing key (minisign format, key id FEF410EFEF5BAB21;
/// docs/releasing.md says how it is kept and rotated).
pub const PUBLIC_KEY: &str = "RWQhq1vv7xD0/uAOt8fqERP1hP+1UCM9U9QHTDYTYytWKiuedTkT7Lgr";
/// The checksum list every release carries.
pub const SUMS: &str = "SHA256SUMS.txt";
/// Its signature by the release key.
pub const SIGNATURE: &str = "SHA256SUMS.txt.minisig";
/// The least time between two automatic checks (seconds).
pub const CHECK_INTERVAL: u64 = 24 * 60 * 60;

/// An update failed.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    /// The release server could not be reached or answered with an error.
    #[error("{0}")]
    Network(String),
    /// The release description is not usable.
    #[error("{0}")]
    Release(String),
    /// A signature or checksum does not match: nothing is installed.
    #[error("{0}")]
    Verify(String),
    /// Unpacking or swapping in failed.
    #[error("{0}")]
    Install(String),
    /// This copy cannot install updates itself (a development build, an unknown layout).
    #[error("this copy of FieldCraft cannot update itself: {0}")]
    NotSupported(String),
}

pub type Result<T> = std::result::Result<T, Error>;

/// The running version as text: the version release.yml built this binary as
/// (`FIELDCRAFT_VERSION`, which carries a nightly's suffix), else the crate version.
pub const VERSION: &str = match option_env!("FIELDCRAFT_VERSION") {
    Some(v) if !v.is_empty() => v,
    _ => env!("CARGO_PKG_VERSION"),
};

/// The running version.
pub fn current() -> Version {
    Version::parse(VERSION).unwrap_or(Version { major: 0, minor: 0, patch: 0, pre: Vec::new() })
}

/// Whether an automatic check is due: checks are on and the last one (seconds since the Unix
/// epoch, 0 for never) was at least [`CHECK_INTERVAL`] ago, or the clock went back past it.
pub fn due(enabled: bool, last_checked: u64, now: u64) -> bool {
    enabled && (last_checked == 0 || now < last_checked || now - last_checked >= CHECK_INTERVAL)
}

/// The release server: [`HOST_ENV`] when set, else GitHub.
pub fn host() -> String {
    std::env::var(HOST_ENV)
        .ok()
        .map(|h| h.trim().trim_end_matches('/').to_string())
        .filter(|h| !h.is_empty())
        .unwrap_or_else(|| DEFAULT_HOST.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn automatic_checks_are_daily_and_only_when_on() {
        let day = CHECK_INTERVAL;
        assert!(due(true, 0, 1_000));
        assert!(!due(true, 1_000, 1_000 + day - 1));
        assert!(due(true, 1_000, 1_000 + day));
        assert!(due(true, 5_000, 1_000), "the clock went back");
        assert!(!due(false, 0, 10 * day), "off means no check");
        assert_eq!(current().to_string(), VERSION);
    }
}
