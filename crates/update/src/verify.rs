//! What makes a download trustworthy: the release's `SHA256SUMS.txt` is signed with the release
//! key (an Ed25519 key in minisign's format, docs/releasing.md) whose public half is built into
//! the app ([`crate::PUBLIC_KEY`]), and every file installed has the SHA-256 that list gives it.
//! A list that is not signed by that key, or a file that does not match it, is refused.

use std::collections::BTreeMap;

use sha2::{Digest, Sha256};

use crate::{Error, Result};

/// The longest checksum list accepted (bytes; FieldCraft's has nine lines).
pub const MAX_SUMS: usize = 64 * 1024;
/// The longest signature file accepted (bytes; a minisign signature is about 300).
pub const MAX_SIGNATURE: usize = 4 * 1024;

/// Check that `sums` is signed by `public_key` (minisign's base64 public key) with `signature`
/// (a `.minisig` file). Only prehashed signatures (minisign's default since 0.8) are accepted.
pub fn signature(sums: &[u8], signature: &str, public_key: &str) -> Result<()> {
    let key =
        minisign_verify::PublicKey::from_base64(public_key).map_err(|e| Error::Verify(format!("the built-in release key is not usable: {e}")))?;
    let sig = minisign_verify::Signature::decode(signature).map_err(|e| Error::Verify(format!("SHA256SUMS.txt.minisig is not a signature: {e}")))?;
    key.verify(sums, &sig, false).map_err(|e| Error::Verify(format!("SHA256SUMS.txt is not signed by the FieldCraft release key ({e})")))
}

/// The file names and SHA-256 digests (lower-case hex) of a `sha256sum` list
/// (`<64 hex>  <name>` per line; a `*` before the name marks binary mode).
pub fn parse_sums(text: &str) -> Result<BTreeMap<String, String>> {
    if text.len() > MAX_SUMS {
        return Err(Error::Verify("SHA256SUMS.txt is too large".into()));
    }
    let mut sums = BTreeMap::new();
    for line in text.lines().map(str::trim_end).filter(|l| !l.is_empty()) {
        let bad = || Error::Verify(format!("SHA256SUMS.txt: `{}` is not a checksum line", line.chars().take(80).collect::<String>()));
        let (hash, name) = line.split_once(' ').ok_or_else(bad)?;
        let name = name.trim_start_matches(' ').trim_start_matches('*');
        if hash.len() != 64 || !hash.bytes().all(|b| b.is_ascii_hexdigit()) || name.is_empty() {
            return Err(bad());
        }
        if sums.insert(name.to_string(), hash.to_ascii_lowercase()).is_some() {
            return Err(Error::Verify(format!("SHA256SUMS.txt lists `{name}` twice")));
        }
    }
    Ok(sums)
}

/// The SHA-256 of `bytes`, as lower-case hex.
pub fn sha256_hex(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

/// Bytes as lower-case hex.
pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Check that file `name`, whose SHA-256 is `actual`, is the one the list names.
pub fn checksum(sums: &BTreeMap<String, String>, name: &str, actual: &str) -> Result<()> {
    let expected = sums.get(name).ok_or_else(|| Error::Verify(format!("SHA256SUMS.txt does not list {name}")))?;
    if !expected.eq_ignore_ascii_case(actual) {
        return Err(Error::Verify(format!(
            "{name} does not match its checksum (expected {expected}, downloaded {actual}): the download is damaged or was changed"
        )));
    }
    Ok(())
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A fresh key pair (the test's own; never the release key) and its public key in base64.
    pub(crate) fn keypair() -> (minisign::SecretKey, String) {
        let kp = minisign::KeyPair::generate_unencrypted_keypair().unwrap();
        let public = kp.pk.to_base64();
        (kp.sk, public)
    }

    /// `data` signed with `sk`, as a `.minisig` file.
    pub(crate) fn sign(sk: &minisign::SecretKey, data: &[u8]) -> String {
        minisign::sign(None, sk, std::io::Cursor::new(data), Some("timestamp:0\tfile:SHA256SUMS.txt"), None).unwrap().into_string()
    }

    #[test]
    fn checksum_lists_parse() {
        let a = "a".repeat(64);
        let sums =
            parse_sums(&format!("{a}  fieldcraft-0.5.9-macos-aarch64.dmg\n{}  *fieldcraft-0.5.9-windows-x64.zip\n\n", "B".repeat(64))).unwrap();
        assert_eq!(sums.len(), 2);
        assert_eq!(sums["fieldcraft-0.5.9-windows-x64.zip"], "b".repeat(64));
        for bad in
            [a.clone(), format!("{}  x", "a".repeat(63)), format!("{}  x", "g".repeat(64)), format!("{a}  x\n{a}  x"), "x".repeat(MAX_SUMS + 1)]
        {
            assert!(parse_sums(&bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn a_checksum_mismatch_is_refused() {
        let sums = parse_sums(&format!("{}  app.zip\n", sha256_hex(b"the app"))).unwrap();
        assert!(checksum(&sums, "app.zip", &sha256_hex(b"the app")).is_ok());
        let e = checksum(&sums, "app.zip", &sha256_hex(b"the app, changed")).unwrap_err();
        assert!(e.to_string().contains("does not match its checksum"), "{e}");
        assert!(checksum(&sums, "other.zip", &sha256_hex(b"the app")).is_err());
    }

    #[test]
    fn only_a_signature_by_the_key_is_accepted() {
        let (sk, public) = keypair();
        let sums = b"abc  fieldcraft.zip\n";
        let sig = sign(&sk, sums);
        assert!(signature(sums, &sig, &public).is_ok());
        // The list changed after signing.
        assert!(signature(b"abd  fieldcraft.zip\n", &sig, &public).is_err());
        // Signed by another key.
        let (other, _) = keypair();
        assert!(signature(sums, &sign(&other, sums), &public).is_err());
        // No signature, a damaged one, a bad key.
        for bad in ["", "untrusted comment: x\n", &sig.replace('A', "B")] {
            assert!(signature(sums, bad, &public).is_err(), "{bad}");
        }
        assert!(signature(sums, &sig, "RWQ").is_err());
        // The built-in release key is a valid key.
        assert!(minisign_verify::PublicKey::from_base64(crate::PUBLIC_KEY).is_ok());
    }

    /// The built-in key accepts what the release workflow's tool (rsign2) signs with the
    /// release key: tests/data holds a list signed that way. Rotating the key re-signs it.
    #[test]
    fn the_release_key_verifies_the_release_tooling() {
        let sums = include_bytes!("../tests/data/SHA256SUMS.txt");
        let sig = include_str!("../tests/data/SHA256SUMS.txt.minisig");
        assert_eq!(signature(sums, sig, crate::PUBLIC_KEY), Ok(()));
        assert!(signature(b"changed", sig, crate::PUBLIC_KEY).is_err());
    }
}
