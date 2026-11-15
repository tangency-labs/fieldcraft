//! A published release, as the GitHub Releases API describes it
//! (`GET /repos/{owner}/{repo}/releases/latest`, `…/releases/tags/{tag}`, `…/releases`), and
//! which one an install is offered ([`Check::choose`]). The reply is hostile input: every string
//! and list is capped, and a draft is refused.

use serde_json::{Value, json};

use crate::{Error, Result, SIGNATURE, SUMS, Version};

/// The most assets a release may list (FieldCraft's have nine).
const MAX_ASSETS: usize = 100;
/// The longest asset name, URL or tag kept.
const MAX_NAME: usize = 512;
/// The most of the release notes kept (bytes).
const MAX_NOTES: usize = 64 * 1024;

/// A file of a release.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Asset {
    pub name: String,
    /// Where it downloads from (`browser_download_url`).
    pub url: String,
    /// Its size in bytes, as the release lists it (0 when not given).
    pub size: u64,
}

/// A published (non-draft) release.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Release {
    pub version: Version,
    /// Its tag (`v0.5.8`).
    pub tag: String,
    /// Its page on GitHub: the release notes and every download.
    pub url: String,
    /// The release notes (Markdown, capped).
    pub notes: String,
    /// When it was published (ISO 8601), as GitHub gives it.
    pub published: String,
    pub prerelease: bool,
    pub assets: Vec<Asset>,
}

/// A string field of `v`, capped at `max` bytes (cut at a character boundary).
fn text(v: &Value, key: &str, max: usize) -> String {
    let s = v.get(key).and_then(Value::as_str).unwrap_or_default();
    let mut end = s.len().min(max);
    while !s.is_char_boundary(end) {
        end -= 1;
    }
    s.get(..end).unwrap_or_default().to_string()
}

impl Release {
    /// The release a GitHub Releases API reply describes.
    pub fn from_github(v: &Value) -> Result<Release> {
        if v.get("draft").and_then(Value::as_bool) == Some(true) {
            return Err(Error::Release("the release is a draft".into()));
        }
        let tag = text(v, "tag_name", MAX_NAME);
        let version =
            Version::parse(&tag).ok_or_else(|| Error::Release(format!("tag `{}` is not a version", tag.chars().take(64).collect::<String>())))?;
        let list = v.get("assets").and_then(Value::as_array).map(Vec::as_slice).unwrap_or_default();
        if list.len() > MAX_ASSETS {
            return Err(Error::Release(format!("the release lists {} assets (at most {MAX_ASSETS})", list.len())));
        }
        let assets = list
            .iter()
            .map(|a| Asset {
                name: text(a, "name", MAX_NAME),
                url: text(a, "browser_download_url", MAX_NAME),
                size: a.get("size").and_then(Value::as_u64).unwrap_or(0),
            })
            .filter(|a| !a.name.is_empty() && !a.url.is_empty())
            .collect();
        Ok(Release {
            version,
            tag,
            url: text(v, "html_url", MAX_NAME),
            notes: text(v, "body", MAX_NOTES),
            published: text(v, "published_at", 64),
            prerelease: v.get("prerelease").and_then(Value::as_bool).unwrap_or(false),
            assets,
        })
    }

    /// The asset named `name`.
    pub fn asset(&self, name: &str) -> Option<&Asset> {
        self.assets.iter().find(|a| a.name == name)
    }

    /// Whether the release carries the signed checksum list an install needs
    /// (`SHA256SUMS.txt` and `SHA256SUMS.txt.minisig`).
    pub fn signed(&self) -> bool {
        self.asset(SUMS).is_some() && self.asset(SIGNATURE).is_some()
    }
}

/// What a check found: the running version and the latest release for it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Check {
    pub current: Version,
    pub latest: Release,
}

impl Check {
    /// The release `current` is offered, from the latest stable release (`releases/latest`;
    /// `None` when none is published) and the published nightlies.
    ///
    /// A stable or candidate install is offered the latest stable release only: never a nightly,
    /// and a pre-release in `stable` is refused (GitHub's `releases/latest` never names one, so
    /// the answer is not GitHub's). A nightly install is offered the newest of the latest stable
    /// release and the nightlies, in [`Version::cmp_build`] order.
    pub fn choose(current: Version, stable: Option<Release>, nightlies: Vec<Release>) -> Result<Check> {
        if let Some(r) = stable.as_ref().filter(|r| r.prerelease || r.version.is_prerelease()) {
            return Err(Error::Release(format!("the latest stable release is a pre-release ({})", r.version)));
        }
        let candidates = stable.into_iter().chain(nightlies.into_iter().filter(|r| current.is_nightly() && r.version.is_nightly()));
        let latest = candidates
            .max_by(|a, b| a.version.cmp_build(&b.version))
            .ok_or_else(|| Error::Release("no FieldCraft release is published yet".into()))?;
        Ok(Check { current, latest })
    }

    /// The latest release is newer than the running version.
    pub fn newer(&self) -> bool {
        self.latest.version.cmp_build(&self.current).is_gt()
    }

    /// As `app.update.check` and `fieldcraft-cli update --check` report it; `asset` is the
    /// file this install would download (`None` when it cannot install itself).
    pub fn to_json(&self, asset: Option<&str>) -> Value {
        json!({
            "current": self.current.to_string(),
            "latest": self.latest.version.to_string(),
            "newer": self.newer(),
            "url": self.latest.url,
            "published": self.latest.published,
            "signed": self.latest.signed(),
            "asset": asset,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn github(tag: &str) -> Value {
        json!({
            "tag_name": tag,
            "html_url": format!("https://github.com/tangency-labs/fieldcraft/releases/tag/{tag}"),
            "body": "## What's new\n- updates",
            "draft": false,
            "prerelease": false,
            "published_at": "2026-10-15T04:06:51Z",
            "assets": [
                {"name": "SHA256SUMS.txt", "browser_download_url": "https://example.invalid/SHA256SUMS.txt", "size": 700},
                {"name": "SHA256SUMS.txt.minisig", "browser_download_url": "https://example.invalid/SHA256SUMS.txt.minisig", "size": 300},
                {"name": "", "browser_download_url": "https://example.invalid/x"},
            ],
        })
    }

    #[test]
    fn a_release_reads_from_the_api_reply() {
        let r = Release::from_github(&github("v0.5.9")).unwrap();
        assert_eq!((r.version.to_string().as_str(), r.tag.as_str(), r.assets.len()), ("0.5.9", "v0.5.9", 2));
        assert!(r.signed() && r.notes.starts_with("## What's new"));
        let check = Check { current: Version::parse("0.5.8").unwrap(), latest: r.clone() };
        assert!(check.newer());
        assert_eq!(check.to_json(Some("a.dmg"))["latest"], "0.5.9");
        assert!(!Check { current: Version::parse("0.5.9").unwrap(), latest: r.clone() }.newer());
        assert!(!Check { current: Version::parse("0.6.0-rc.1").unwrap(), latest: r }.newer());

        let mut unsigned = github("v0.5.9");
        unsigned["assets"] = json!([]);
        assert!(!Release::from_github(&unsigned).unwrap().signed());
        let mut draft = github("v0.5.9");
        draft["draft"] = json!(true);
        assert!(Release::from_github(&draft).is_err());
        for hostile in [
            json!(null),
            json!([]),
            json!({"tag_name": 7}),
            json!({"tag_name": "latest"}),
            json!({"tag_name": "v1.0.0", "assets": vec![json!({}); 1000]}),
        ] {
            assert!(Release::from_github(&hostile).is_err(), "{hostile}");
        }
        let mut long = github("v1.0.0");
        long["body"] = json!("é".repeat(MAX_NOTES));
        assert!(Release::from_github(&long).unwrap().notes.len() <= MAX_NOTES);
    }

    fn release(tag: &str) -> Release {
        let mut r = github(tag);
        r["prerelease"] = json!(tag.contains('-'));
        Release::from_github(&r).unwrap()
    }

    fn offered(current: &str, stable: Option<&str>, nightlies: &[&str]) -> Option<String> {
        let check = Check::choose(Version::parse(current).unwrap(), stable.map(release), nightlies.iter().map(|t| release(t)).collect()).unwrap();
        check.newer().then(|| check.latest.version.to_string())
    }

    /// A stable install is never offered a nightly; a nightly install is offered the newest
    /// nightly, or a stable release newer than the one it is numbered from.
    #[test]
    fn each_install_is_offered_its_own_channel() {
        let nightlies = ["v0.8.0-nightly.20261117.41.1", "v0.8.0-nightly.20261116.40.1"];
        assert_eq!(offered("0.7.5", Some("v0.8.0"), &nightlies).as_deref(), Some("0.8.0"));
        assert_eq!(offered("0.8.0", Some("v0.8.0"), &nightlies), None, "a stable install never sees nightlies");
        assert_eq!(offered("0.8.0-rc.1", Some("v0.7.5"), &nightlies), None);
        assert_eq!(offered("0.8.0-nightly.20261116.40.1", Some("v0.8.0"), &nightlies).as_deref(), Some("0.8.0-nightly.20261117.41.1"));
        assert_eq!(offered("0.8.0-nightly.20261117.41.1", Some("v0.8.0"), &nightlies), None, "its own release is not newer");
        assert_eq!(offered("0.8.0-nightly.20261117.41.1", Some("v0.8.1"), &nightlies).as_deref(), Some("0.8.1"));
        assert_eq!(offered("0.8.0-nightly.20261117.41.1", None, &["v0.8.1-nightly.20261201.60.1"]).as_deref(), Some("0.8.1-nightly.20261201.60.1"));
        assert_eq!(offered("0.8.0-nightly.20261116.40.1", Some("v0.7.5"), &["v0.8.0-rc.2"]), None, "a candidate is not a nightly");

        // A pre-release where the latest stable release belongs is refused; no release at all is an error.
        let current = Version::parse("0.8.0").unwrap();
        assert!(Check::choose(current.clone(), Some(release("v0.9.0-nightly.20261201.60.1")), Vec::new()).is_err());
        let mut flagged = github("v0.9.0");
        flagged["prerelease"] = json!(true);
        assert!(Check::choose(current.clone(), Some(Release::from_github(&flagged).unwrap()), Vec::new()).is_err());
        assert!(Check::choose(current, None, vec![release("v0.9.0-nightly.20261201.60.1")]).is_err());
    }
}
