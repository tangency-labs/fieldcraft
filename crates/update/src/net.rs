//! Talking to the release server: the GitHub Releases API for the release description, and its
//! download links for the files. HTTPS only (checked against the operating system's
//! certificates), except against a local test server named by [`crate::HOST_ENV`].

use std::io::{Read, Write};
use std::path::Path;
use std::time::Duration;

use sha2::{Digest, Sha256};

use crate::{Check, Error, REPO, Release, Result, Version};

/// The longest a release description may be (bytes).
const MAX_JSON: u64 = 4 * 1024 * 1024;
/// The largest asset downloaded (bytes; FieldCraft's are well under 200 MB).
pub const MAX_ASSET: u64 = 2 * 1024 * 1024 * 1024;
/// The most entries a release list may have (GitHub's largest page, `per_page=100`).
const MAX_LISTED: usize = 100;
/// How long connecting and the release description may take.
const TIMEOUT: Duration = Duration::from_secs(15);
/// How long a whole download may take.
const DOWNLOAD_TIME: Duration = Duration::from_secs(30 * 60);

/// Which release to ask for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Which {
    /// The latest published, non-draft, non-pre-release release (`releases/latest`).
    Latest,
    /// The release of this version (to reinstall it).
    Version(Version),
}

/// A client for one release server.
pub struct Client {
    host: String,
    agent: ureq::Agent,
    downloads: ureq::Agent,
}

impl Client {
    /// A client for `host` (`https://api.github.com`, or a local test server).
    pub fn new(host: &str) -> Client {
        // rustls's process-wide crypto (ring); another part of the app may have installed it.
        let _ = rustls::crypto::ring::default_provider().install_default();
        let tls = || ureq::tls::TlsConfig::builder().root_certs(ureq::tls::RootCerts::PlatformVerifier).build();
        let agent = ureq::Agent::config_builder().timeout_global(Some(TIMEOUT)).tls_config(tls()).build().into();
        let downloads = ureq::Agent::config_builder()
            .timeout_connect(Some(TIMEOUT))
            .timeout_recv_response(Some(TIMEOUT))
            .timeout_recv_body(Some(DOWNLOAD_TIME))
            .tls_config(tls())
            .build()
            .into();
        Client { host: host.trim_end_matches('/').to_string(), agent, downloads }
    }

    /// The client for [`crate::host`].
    pub fn from_env() -> Client {
        Client::new(&crate::host())
    }

    /// Whether `url` may be fetched: HTTPS, or anything when testing against a local server.
    fn allowed(&self, url: &str) -> Result<()> {
        if url.starts_with("https://") || self.host != crate::DEFAULT_HOST {
            Ok(())
        } else {
            Err(Error::Network(format!("refusing a download that is not HTTPS: {}", url.chars().take(120).collect::<String>())))
        }
    }

    /// The JSON the API answers at `path`; `None` when it answers 404 (no such release).
    fn api(&self, path: &str) -> Result<Option<serde_json::Value>> {
        let url = format!("{}{path}", self.host);
        let reply = self
            .agent
            .get(&url)
            .header("Accept", "application/vnd.github+json")
            .header("User-Agent", &format!("FieldCraft/{}", crate::current()))
            .config()
            .http_status_as_error(false)
            .build()
            .call()
            .map_err(|e| Error::Network(format!("could not reach the release server: {e}")))?;
        let status = reply.status().as_u16();
        let body = reply
            .into_body()
            .with_config()
            .limit(MAX_JSON)
            .read_to_string()
            .map_err(|e| Error::Network(format!("the release server's reply: {e}")))?;
        match status {
            200 => {}
            404 => return Ok(None),
            403 | 429 => return Err(Error::Network("the release server is limiting requests; try again later".into())),
            other => return Err(Error::Network(format!("the release server answered {other}"))),
        }
        serde_json::from_str(&body).map(Some).map_err(|e| Error::Release(format!("the release description is not JSON: {e}")))
    }

    /// The release `which` describes.
    pub fn release(&self, which: &Which) -> Result<Release> {
        let path = match which {
            Which::Latest => format!("/repos/{REPO}/releases/latest"),
            Which::Version(v) => format!("/repos/{REPO}/releases/tags/v{v}"),
        };
        let json = self.api(&path)?.ok_or_else(|| {
            Error::Release(match which {
                Which::Latest => "no FieldCraft release is published yet".into(),
                Which::Version(v) => format!("there is no published release of FieldCraft {v}"),
            })
        })?;
        Release::from_github(&json)
    }

    /// The published nightlies among the newest releases (`GET …/releases`, GitHub's first page:
    /// newest first, and nightlies are published daily, so the newest nightly is on it).
    /// Entries that are not nightlies, or not usable, are left out.
    pub fn nightlies(&self) -> Result<Vec<Release>> {
        let json = self.api(&format!("/repos/{REPO}/releases"))?.unwrap_or_default();
        let list = json.as_array().ok_or_else(|| Error::Release("the release list is not a list".into()))?;
        if list.len() > MAX_LISTED {
            return Err(Error::Release(format!("the release list has {} entries (at most {MAX_LISTED})", list.len())));
        }
        Ok(list.iter().filter_map(|v| Release::from_github(v).ok()).filter(|r| r.version.is_nightly()).collect())
    }

    /// Compare the latest release for the running version with it ([`Check::choose`]).
    pub fn check(&self) -> Result<Check> {
        self.check_for(crate::current())
    }

    /// [`Client::check`] for an install of `current`: a nightly also asks for the nightlies.
    pub fn check_for(&self, current: Version) -> Result<Check> {
        let stable = match self.api(&format!("/repos/{REPO}/releases/latest"))? {
            Some(json) => Some(Release::from_github(&json)?),
            None => None,
        };
        let nightlies = if current.is_nightly() { self.nightlies()? } else { Vec::new() };
        Check::choose(current, stable, nightlies)
    }

    /// A small file (the checksum list, its signature), at most `max` bytes.
    pub fn small(&self, url: &str, max: usize) -> Result<Vec<u8>> {
        self.allowed(url)?;
        let reply = self.agent.get(url).header("User-Agent", "FieldCraft").call().map_err(|e| Error::Network(format!("download failed: {e}")))?;
        reply.into_body().with_config().limit(max as u64).read_to_vec().map_err(|e| Error::Network(format!("download failed: {e}")))
    }

    /// Download `url` to `to`, at most `max` bytes; its SHA-256 (hex). `progress` hears the bytes
    /// so far and the total (0 when not known).
    pub fn download(&self, url: &str, to: &Path, max: u64, progress: &dyn Fn(u64, u64)) -> Result<String> {
        self.allowed(url)?;
        let reply = self.downloads.get(url).header("User-Agent", "FieldCraft").call().map_err(|e| Error::Network(format!("download failed: {e}")))?;
        let total = reply.headers().get("content-length").and_then(|v| v.to_str().ok()).and_then(|v| v.parse::<u64>().ok()).unwrap_or(0);
        if total > max {
            return Err(Error::Network(format!("the download is {total} bytes, more than the {max} accepted")));
        }
        let mut body = reply.into_body();
        let mut reader = body.with_config().limit(max).reader();
        let mut file = std::fs::File::create(to).map_err(|e| Error::Install(format!("{}: {e}", to.display())))?;
        let mut hash = Sha256::new();
        let mut buf = vec![0u8; 256 * 1024];
        let mut done = 0u64;
        loop {
            let n = reader.read(&mut buf).map_err(|e| Error::Network(format!("download failed: {e}")))?;
            let Some(chunk) = buf.get(..n).filter(|c| !c.is_empty()) else { break };
            hash.update(chunk);
            file.write_all(chunk).map_err(|e| Error::Install(format!("{}: {e}", to.display())))?;
            done += n as u64;
            progress(done, total);
        }
        file.sync_all().map_err(|e| Error::Install(format!("{}: {e}", to.display())))?;
        Ok(crate::verify::hex(&hash.finalize()))
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;
    use crate::mock::Server;

    #[test]
    fn the_latest_release_is_asked_for_and_compared() {
        let server = Server::start().unwrap();
        server.release("latest", &json!({"tag_name": "v99.0.0", "html_url": "https://example.invalid/r", "assets": []}));
        let check = Client::new(&server.url()).check().unwrap();
        assert!(check.newer());
        assert_eq!(check.latest.version.to_string(), "99.0.0");
        assert_eq!(server.requests(), 1);

        // Errors are messages, never panics.
        let missing = Client::new(&server.url()).release(&Which::Version(Version::parse("1.2.3").unwrap())).unwrap_err();
        assert!(missing.to_string().contains("no published release of FieldCraft 1.2.3"), "{missing}");
        server.route("/repos/tangency-labs/fieldcraft/releases/latest", 200, b"not json".to_vec());
        assert!(matches!(Client::new(&server.url()).check(), Err(Error::Release(_))));
        server.route("/repos/tangency-labs/fieldcraft/releases/latest", 403, b"{}".to_vec());
        assert!(matches!(Client::new(&server.url()).check(), Err(Error::Network(_))));
        assert!(matches!(Client::new("http://127.0.0.1:1").check(), Err(Error::Network(_))));
    }

    /// `releases/latest` never names a pre-release (GitHub's definition), but the answer comes
    /// from the network: a stable install refuses one rather than offer it.
    #[test]
    fn a_stable_install_is_never_offered_a_nightly() {
        let server = Server::start().unwrap();
        let nightly = json!({"tag_name": "v99.0.0-nightly.20261116.40.1", "prerelease": true, "assets": []});
        server.release("latest", &nightly);
        let offered = Client::new(&server.url()).check().map(|c| c.newer());
        assert!(!matches!(offered, Ok(true)), "a stable install was offered a nightly: {offered:?}");
    }

    #[test]
    fn a_nightly_install_is_offered_the_newest_nightly() {
        let server = Server::start().unwrap();
        server.release("latest", &json!({"tag_name": "v0.8.0", "assets": []}));
        let nightly = |tag: &str| json!({"tag_name": tag, "prerelease": true, "assets": []});
        server.releases(&json!([
            nightly("v0.8.0-nightly.20261117.41.1"),
            {"tag_name": "v0.8.0", "assets": []},
            {"tag_name": "v0.8.0-rc.1", "prerelease": true, "assets": []},
            nightly("v0.8.0-nightly.20261116.40.1"),
            {"tag_name": "not a version"},
        ]));
        let client = Client::new(&server.url());
        let check = client.check_for(Version::parse("0.8.0-nightly.20261116.40.1").unwrap()).unwrap();
        assert!(check.newer());
        assert_eq!(check.latest.version.to_string(), "0.8.0-nightly.20261117.41.1");
        assert_eq!(server.requests(), 2);
        // A stable install asks for the latest stable release only.
        let check = client.check_for(Version::parse("0.7.5").unwrap()).unwrap();
        assert_eq!((check.latest.version.to_string().as_str(), server.requests()), ("0.8.0", 3));

        // Before any stable release, a nightly still finds the nightlies; a hostile list is an error.
        let server = Server::start().unwrap();
        server.releases(&json!([nightly("v0.8.0-nightly.20261117.41.1")]));
        let client = Client::new(&server.url());
        assert!(client.check_for(Version::parse("0.8.0-nightly.20261116.40.1").unwrap()).unwrap().newer());
        server.releases(&json!(vec![nightly("v0.8.0-nightly.20261117.41.1"); MAX_LISTED + 1]));
        assert!(matches!(client.check_for(Version::parse("0.8.0-nightly.20261116.40.1").unwrap()), Err(Error::Release(_))));
        server.releases(&json!({"message": "not a list"}));
        assert!(matches!(client.check_for(Version::parse("0.8.0-nightly.20261116.40.1").unwrap()), Err(Error::Release(_))));
    }

    #[test]
    fn downloads_are_hashed_and_capped() {
        let server = Server::start().unwrap();
        server.route("/f.bin", 200, b"0123456789".to_vec());
        let dir = std::env::temp_dir().join(format!("fieldcraft-update-net-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let to = dir.join("f.bin");
        let c = Client::new(&server.url());
        let sha = c.download(&format!("{}/f.bin", server.url()), &to, 100, &|_, _| {}).unwrap();
        assert_eq!(sha, crate::verify::sha256_hex(b"0123456789"));
        assert!(c.download(&format!("{}/f.bin", server.url()), &to, 5, &|_, _| {}).is_err());
        assert!(c.small(&format!("{}/f.bin", server.url()), 5).is_err());
        // Against GitHub, only HTTPS.
        assert!(Client::new(crate::DEFAULT_HOST).small("http://example.invalid/x", 5).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
