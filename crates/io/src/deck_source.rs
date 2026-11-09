//! Immutable original input files; archival data, never solver state or filesystem authority.
use crate::{Error, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
pub const PART: &str = "deck-source.json";
pub const MAX_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_JSON_BYTES: u64 = (MAX_BYTES as u64) * 6 + 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DeckSource {
    pub root: String,
    pub files: BTreeMap<String, String>,
}
impl DeckSource {
    pub fn validate(&self) -> Result<()> {
        if self.files.is_empty() || self.files.len() > 1024 || !self.files.contains_key(&self.root) {
            return Err(Error::Invalid("deck source needs a root in at most 1024 files".into()));
        }
        let mut bytes = 0usize;
        for (name, text) in &self.files {
            if name.len() > 256 || name.contains(['\\', ':', '\0']) || name.split('/').any(|s| s.is_empty() || s == "." || s == "..") {
                return Err(Error::Invalid("deck source path must stay relative to its export directory".into()));
            }
            bytes = bytes.saturating_add(text.len());
            if bytes > MAX_BYTES {
                return Err(Error::Limit("deck source exceeds 64 MiB".into()));
            }
        }
        Ok(())
    }
}
