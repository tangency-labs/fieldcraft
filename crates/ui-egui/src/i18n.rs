//! Interface localization (`en` / `he`) and BCP 47 / legacy `iw` tag normalization.
//!
//! `egui 0.37` and `solidcraft-ui-text` handle visual UAX #9 BiDi line formatting inside
//! `Catalog::shown()`, so `locales/he.json` stores logical Hebrew text.

/// Preference value that follows the operating system or browser language.
pub const SYSTEM: &str = solidcraft_ui_text::SYSTEM;
/// English catalog id (fallback).
pub const ENGLISH: &str = "en";
/// Hebrew catalog id (`עברית`, right-to-left).
pub const HEBREW: &str = "he";

/// A selectable language entry for Settings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LanguageChoice {
    pub pref: &'static str,
    pub label: String,
}

/// Strip POSIX encoding/modifier suffixes, normalize separators to `-`, lowercase, reject `x-*`
/// private-use tags, and map legacy ISO 639-1 `iw` primary subtags to `he`.
pub fn normalize_language_tag(tag: &str) -> String {
    let base = tag.trim().split(['.', '@']).next().unwrap_or_default().replace('_', "-").to_ascii_lowercase();
    let mut parts = base.split('-').filter(|part| !part.is_empty());
    let Some(primary) = parts.next() else {
        return String::new();
    };
    if primary == "x" {
        return String::new();
    }
    let normalized_primary = if primary == "iw" { HEBREW } else { primary };
    let mut out = String::from(normalized_primary);
    for part in parts {
        out.push('-');
        out.push_str(part);
    }
    out
}

/// Normalize and validate a saved or requested language preference (`system`, `en`, `he`, `iw`).
pub fn validate_preference(pref: &str) -> Result<String, String> {
    let trimmed = pref.trim();
    if trimmed.eq_ignore_ascii_case(SYSTEM) {
        return Ok(SYSTEM.into());
    }
    let normalized = normalize_language_tag(trimmed);
    let primary = normalized.split('-').next().unwrap_or_default();
    match primary {
        "en" => Ok(ENGLISH.into()),
        "he" => Ok(HEBREW.into()),
        _ => Err(format!("unknown language {pref:?} (expected system, en, or he)")),
    }
}

/// Resolve a stored/requested preference (`system`, `en`, `he`, `iw`) to an installed catalog id.
pub fn resolve_preference(pref: &str) -> &'static str {
    match validate_preference(pref).as_deref() {
        Ok(SYSTEM) => solidcraft_ui_text::system_language(),
        Ok(HEBREW) => HEBREW,
        _ => ENGLISH,
    }
}

/// Install the catalogs and enter a scoped thread-local locale for `pref`.
pub fn scope(pref: &str) -> solidcraft_ui_text::LanguageScope {
    crate::theme::install_look();
    solidcraft_ui_text::scope(resolve_preference(pref))
}

/// Language options shown in Settings: the system's, then each catalogue's own name.
pub fn selectable_languages() -> Vec<LanguageChoice> {
    [ENGLISH, HEBREW].into_iter().fold(
        vec![LanguageChoice { pref: SYSTEM, label: solidcraft_ui_text::tr("settings.language_system").to_owned() }],
        |mut choices, id| {
            choices.push(LanguageChoice { pref: id, label: solidcraft_ui_text::Locale::for_language(id).display_name().into_owned() });
            choices
        },
    )
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeSet, HashMap};

    use super::*;

    fn placeholder_names(text: &str) -> Vec<String> {
        let mut names = Vec::new();
        let mut rest = text;
        while let Some((_, after)) = rest.split_once('{') {
            let Some((name, tail)) = after.split_once('}') else {
                break;
            };
            names.push(name.to_owned());
            rest = tail;
        }
        names.sort();
        names
    }

    #[test]
    fn legacy_iw_normalization() {
        assert_eq!(normalize_language_tag("iw"), "he");
        assert_eq!(normalize_language_tag("iw-IL"), "he-il");
        assert_eq!(normalize_language_tag("iw_IL.UTF-8"), "he-il");
        assert_eq!(normalize_language_tag("he_IL@euro"), "he-il");
        assert_eq!(normalize_language_tag("x-he"), "");

        assert_eq!(validate_preference("system"), Ok("system".into()));
        assert_eq!(validate_preference("en-US"), Ok("en".into()));
        assert_eq!(validate_preference("he-IL"), Ok("he".into()));
        assert_eq!(validate_preference("iw_IL.UTF-8"), Ok("he".into()));
        assert!(validate_preference("x-he").is_err());
        assert!(validate_preference("xx").is_err());
    }

    #[test]
    fn hebrew_catalog_covers_all_english_and_workspace_keys() {
        let en: HashMap<String, String> = serde_json::from_str(include_str!("../locales/en.json")).expect("valid en.json");
        let he: HashMap<String, String> = serde_json::from_str(include_str!("../locales/he.json")).expect("valid he.json");

        let mut expected_keys: BTreeSet<String> = en.keys().cloned().collect();
        for w in fieldcraft_engine::catalogue::Workspace::ALL {
            let (about, tabs) = crate::toolbar::workspace_ids(w);
            expected_keys.insert(about.to_owned());
            expected_keys.insert(tabs.to_owned());
        }

        let missing: Vec<&String> = expected_keys.iter().filter(|k| !he.contains_key(k.as_str())).collect();
        assert!(missing.is_empty(), "missing Hebrew translations in locales/he.json: {missing:?}");

        for (key, he_val) in &he {
            assert!(!he_val.trim().is_empty(), "{key} has empty Hebrew translation");
            if let Some(en_val) = en.get(key) {
                assert_eq!(placeholder_names(en_val), placeholder_names(he_val), "placeholder mismatch for {key}");
                if !matches!(key.as_str(), "settings.group_app" | "settings.updates_skipped_version" | "links.discord_label" | "links.discord_name") {
                    assert_ne!(en_val, he_val, "{key} was left untranslated in locales/he.json");
                }
            }
        }
    }

    #[test]
    fn scoped_hebrew_translation() {
        crate::theme::install_look();
        assert_eq!(resolve_preference("iw-IL"), "he");
        {
            let _guard = scope("iw");
            assert_eq!(solidcraft_ui_text::language(), "he");
            assert_ne!(solidcraft_ui_text::tr("settings.title"), "Settings");
            assert_ne!(solidcraft_ui_text::tr("settings.language"), "Language");
            assert_ne!(solidcraft_ui_text::tr("workspace.prepare.about"), "Open geometry and mesh it");
            let formatted = solidcraft_ui_text::format("settings.updates_skipped_version", &[("version", "0.8.1")]);
            assert!(formatted.contains("0.8.1") && formatted.contains("FieldCraft"), "got {formatted:?}");
            let labels: Vec<(&str, String)> = selectable_languages().into_iter().map(|c| (c.pref, c.label)).collect();
            assert_eq!(labels.iter().map(|(p, _)| *p).collect::<Vec<_>>(), [SYSTEM, ENGLISH, HEBREW]);
            assert_eq!(labels[1].1, "English");
            assert_eq!(labels[2].1, "עברית");
        }
    }
}
