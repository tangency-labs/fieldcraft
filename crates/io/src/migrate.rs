//! Format versions and the migration chain of `.fieldcraft` files.
//!
//! The manifest's `format` is `fieldcraft/<N>`, N an integer that grows by one with every change
//! to what a file holds or how ([`FORMAT_VERSION`] is the one this version writes). Parts with an
//! encoding of their own also carry its version: `mesh.bin` starts with `MCMESH01`.
//!
//! A file of any format up to [`FORMAT_VERSION`] is read: its parts are parsed (the JSON parts to
//! [`Value`], `mesh.bin` kept as bytes), then [`migrate`] applies the steps `N → N+1`, one after
//! the other, up to the current format, and only then are the parts deserialised into the
//! model's types, with every check a current file gets. A file of a later format is refused
//! before anything is deserialised ([`version`]), so the user is told to update FieldCraft rather
//! than about a field this version does not know.
//!
//! **Rules for steps.**
//! - A step is a pure function of the parts: no clock, no file system, no settings. The same
//!   file always migrates to the same model.
//! - Each step takes format N exactly as it was released to format N+1 and sets the manifest's
//!   `format` to N+1. It has its own tests.
//! - **A step is never edited once released** (a FieldCraft release has written format N+1).
//!   Files that went through it are already out there, and later steps were written against its
//!   output. A mistake is fixed by a new step and a new format number.
//! - Steps change the parts' content, never the package (the ZIP, the parts' names, the
//!   manifest's `parts` sizes and hashes): those are checked against the bytes on disk before
//!   any step runs.

use serde_json::Value;

use crate::{Error, FORMAT_VERSION, Result, invalid};

/// The format tag's prefix (`fieldcraft/<N>`).
const TAG_PREFIX: &str = "fieldcraft/";

/// Longest generator name quoted in a message.
const MAX_GENERATOR_LEN: usize = 64;

/// A file's parts as read, before they are deserialised: what migration steps work on.
#[derive(Debug, Clone, PartialEq)]
pub struct Parts {
    /// `fieldcraft.json`.
    pub manifest: Value,
    /// `model.json`: the definition.
    pub model: Value,
    /// `mesh.bin`, whose own magic says its encoding's version.
    pub mesh: Vec<u8>,
    /// `journal.json`.
    pub journal: Value,
}

/// One step: format N to N+1.
type Step = fn(&mut Parts) -> Result<()>;

/// `STEPS[i]` takes format `i + 1` to `i + 2`.
const STEPS: [Step; FORMAT_VERSION as usize - 1] = [v1_to_v2];

/// The format version a manifest declares, checked before anything else in it is read: a later
/// format than this version reads is [`Error::Version`] (naming the program that wrote it), a
/// missing or malformed tag [`Error::Invalid`].
pub fn version(manifest: &Value) -> Result<u32> {
    let not_ours = || invalid(format!("not a FieldCraft file: `fieldcraft.json` has no `{TAG_PREFIX}<N>` format tag"));
    let Some(tag) = manifest.get("format").and_then(Value::as_str) else { return not_ours() };
    let Some(number) = tag.strip_prefix(TAG_PREFIX) else { return not_ours() };
    if number.is_empty() || !number.bytes().all(|b| b.is_ascii_digit()) || number.starts_with('0') {
        return not_ours();
    }
    // Too many digits for a u32 is a later format still.
    match number.parse::<u32>() {
        Ok(n) if n <= FORMAT_VERSION => Ok(n),
        _ => Err(Error::Version(newer_message(manifest, number))),
    }
}

/// The message for a file of a later format: who made it, its format and what this version reads.
fn newer_message(manifest: &Value, number: &str) -> String {
    let number: String = number.chars().take(12).collect();
    let generator = manifest
        .get("generator")
        .and_then(Value::as_str)
        // Printable ASCII only: no control, bidirectional or other format characters reach the
        // message (a generator is `FieldCraft <version>`).
        .map(|g| g.chars().filter(|c| c.is_ascii_graphic() || *c == ' ').take(MAX_GENERATOR_LEN).collect::<String>())
        .filter(|g| !g.trim().is_empty());
    let made_by = match generator {
        Some(g) => format!("made by {g}"),
        None => "made by a later FieldCraft".to_string(),
    };
    format!(
        "this model was {made_by} (format {number}); FieldCraft {} reads formats up to {FORMAT_VERSION}. Update FieldCraft to open it.",
        env!("CARGO_PKG_VERSION")
    )
}

/// Bring `parts` from the format their manifest declares to the current one
/// ([`FORMAT_VERSION`]); returns the format they were in.
pub fn migrate(parts: &mut Parts) -> Result<u32> {
    let from = version(&parts.manifest)?;
    for (i, step) in STEPS.iter().enumerate().skip(from as usize - 1) {
        step(parts)?;
        let reached = i as u32 + 2;
        if parts.manifest.get("format").and_then(Value::as_str) != Some(tag(reached).as_str()) {
            return invalid(format!("the migration to format {reached} did not reach it"));
        }
    }
    Ok(from)
}

/// Run the one step that takes `parts`, of format `from`, to `from + 1` (for tests of a step on
/// real files; [`migrate`] runs the chain).
pub fn step(parts: &mut Parts, from: u32) -> Result<()> {
    let run = from.checked_sub(1).and_then(|i| STEPS.get(i as usize)).ok_or_else(|| Error::Version(format!("no step from format {from}")))?;
    run(parts)
}

/// The format tag of version `n`.
pub fn tag(n: u32) -> String {
    format!("{TAG_PREFIX}{n}")
}

fn set_format(parts: &mut Parts, n: u32) -> Result<()> {
    match parts.manifest.as_object_mut() {
        Some(m) => {
            m.insert("format".into(), Value::String(tag(n)));
            Ok(())
        }
        None => invalid("`fieldcraft.json` is not a JSON object"),
    }
}

/// Format 1 → 2 (FieldCraft 0.7). Format 2 brought format versioning and migration; the parts
/// hold the same content, so only the tag changes. Released: never edit.
fn v1_to_v2(parts: &mut Parts) -> Result<()> {
    set_format(parts, 2)
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;
    use crate::FORMAT;

    fn parts(format: &str) -> Parts {
        Parts { manifest: json!({"format": format, "generator": "FieldCraft 0.7.0"}), model: json!({}), mesh: vec![1, 2], journal: json!([]) }
    }

    #[test]
    fn the_tag_and_the_number_agree() {
        assert_eq!(FORMAT, tag(FORMAT_VERSION));
        assert_eq!(STEPS.len() + 1, FORMAT_VERSION as usize);
    }

    #[test]
    fn versions_are_read_from_the_tag() {
        for n in 1..=FORMAT_VERSION {
            assert_eq!(version(&json!({"format": tag(n)})).unwrap(), n);
        }
        for bad in [
            json!({}),
            json!({"format": 2}),
            json!({"format": "fieldcraft/"}),
            json!({"format": "fieldcraft/0"}),
            json!({"format": "fieldcraft/01"}),
            json!({"format": "fieldcraft/-1"}),
            json!({"format": "fieldcraft/2.0"}),
            json!({"format": "solidcraft/2"}),
            json!({"format": " fieldcraft/2"}),
            json!([]),
        ] {
            assert!(matches!(version(&bad), Err(Error::Invalid(_))), "{bad}");
        }
    }

    #[test]
    fn a_later_format_names_its_generator_and_what_this_version_reads() {
        let later = FORMAT_VERSION + 1;
        let e = version(&json!({"format": tag(later), "generator": "FieldCraft 9.1.0"})).unwrap_err();
        let Error::Version(m) = e else { panic!("{e:?}") };
        assert!(m.contains("made by FieldCraft 9.1.0") && m.contains(&format!("(format {later})")), "{m}");
        assert!(m.contains(&format!("reads formats up to {FORMAT_VERSION}")) && m.ends_with("Update FieldCraft to open it."), "{m}");
        // A hostile generator is cut short and loses its control characters; a missing one is
        // a later FieldCraft; a number past u32 is a later format too.
        let e = version(
            &json!({"format": "fieldcraft/99999999999999999999999", "generator": format!("x\u{1b}[2J\u{202e}\u{200b}{}", "y".repeat(10_000))}),
        )
        .unwrap_err();
        let Error::Version(m) = e else { panic!("{e:?}") };
        assert!(m.len() < 400 && m.is_ascii() && !m.contains('\u{1b}'), "{m}");
        let e = version(&json!({"format": tag(later)})).unwrap_err();
        assert!(e.to_string().contains("made by a later FieldCraft"), "{e}");
    }

    #[test]
    fn v1_to_v2_changes_only_the_tag() {
        let mut p = parts("fieldcraft/1");
        let before = p.clone();
        v1_to_v2(&mut p).unwrap();
        assert_eq!(p.manifest, json!({"format": "fieldcraft/2", "generator": "FieldCraft 0.7.0"}));
        assert_eq!((&p.model, &p.mesh, &p.journal), (&before.model, &before.mesh, &before.journal));
        assert!(v1_to_v2(&mut Parts { manifest: json!([]), ..before }).is_err());
    }

    #[test]
    fn the_chain_runs_from_any_version_to_the_current_one() {
        for n in 1..=FORMAT_VERSION {
            let mut p = parts(&tag(n));
            assert_eq!(migrate(&mut p).unwrap(), n);
            assert_eq!(p.manifest["format"], FORMAT);
        }
        for bad in ["fieldcraft/0", "fieldcraft/x", &tag(FORMAT_VERSION + 1)] {
            let mut p = parts(bad);
            assert!(migrate(&mut p).is_err(), "{bad}");
            assert_eq!(p, parts(bad), "a refused file is left as it was");
        }
    }
}
