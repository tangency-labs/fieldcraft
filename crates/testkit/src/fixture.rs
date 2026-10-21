//! The corpus parts' fixtures: each part's accepted results, frozen in
//! `tests/corpus/<part>/fixture.json` (an analytic shape's directory holds only its fixture) by
//! `cargo xtask corpus --bless`, so a mesher change that moves one fails the corpus run and CI's
//! heavy job (docs/meshing/corpus.md says how and when to bless).
//!
//! The format is the oracle's ([`crate::oracle`]): `values` maps a measure to its accepted
//! `value` and the `band` it may move by. Beside it, `stages` gives each stage's outcome, which is
//! compared exactly: `ok` (with the method that meshed the part, for the hex stage), `no: <kind>`
//! (the mesher refused the part, with its error's kind) or `error: <kind>`. Wall times are never
//! frozen.
//!
//! **Platforms.** The tessellation of curved faces differs in its last bits between platforms
//! (the system maths library), and the tet mesher turns that into a different mesh, so a part's
//! results can differ by platform. `platforms` names every platform the fixture was blessed on
//! ([`platform`]) and, for each, what it measured differently from `stages` and `values`: a
//! stage's outcome, a value outside the band (`null` when the platform does not measure it).
//! A run checks against `stages` and `values` with its own platform's differences laid over them.
//! Blessing on a platform the fixture already records rewrites the fixture from that platform
//! alone when its results moved (the other platforms' results are stale then, and blessing there
//! adds them back); blessing on a new platform adds its differences.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::oracle::{FIXTURE_RELATIVE, Frozen};

/// The absolute part of a real measure's band, in the measure's own unit (mm³, degrees, or a
/// ratio): far below any change that matters, it keeps a measure that is zero from having no band.
const ABSOLUTE_FLOOR: f64 = 1e-12; // unit-factors: a floor in the measure's unit
const NOTE: &str = "The corpus part's accepted results, written by `cargo xtask corpus --bless` (docs/meshing/corpus.md). A stage whose outcome changes, or a value outside its band, fails the corpus run and CI's heavy job. `platforms` lists what each platform blessed measured differently.";

/// The platform a run is on, as fixtures name it: `linux-x86_64`, `macos-aarch64`, …
pub fn platform() -> String {
    format!("{}-{}", std::env::consts::OS, std::env::consts::ARCH)
}

/// A number a stage measured.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Metric {
    /// Frozen exactly.
    Count(usize),
    /// Frozen with a band of [`FIXTURE_RELATIVE`] of its size (and [`ABSOLUTE_FLOOR`]).
    Real(f64),
    /// A relative difference (a volume against the exact one), frozen within [`FIXTURE_RELATIVE`]
    /// of one: its own size says nothing about the precision of what it compares.
    Relative(f64),
}

impl Metric {
    pub fn value(self) -> f64 {
        match self {
            Metric::Count(n) => n as f64,
            Metric::Real(x) | Metric::Relative(x) => x,
        }
    }

    /// The band it gets when it is blessed.
    pub fn band(self) -> f64 {
        match self {
            Metric::Count(_) => 0.0,
            Metric::Real(x) => FIXTURE_RELATIVE * x.abs() + ABSOLUTE_FLOOR,
            Metric::Relative(_) => FIXTURE_RELATIVE,
        }
    }

    fn frozen(self) -> Frozen {
        Frozen { value: self.value(), band: self.band() }
    }
}

/// What a corpus run found for one part: each stage's outcome and the measures, keyed
/// `<stage>.<measure>`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Results {
    pub stages: BTreeMap<String, String>,
    pub values: BTreeMap<String, Metric>,
}

/// What one platform measured differently from the fixture's `stages` and `values`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Differences {
    pub stages: BTreeMap<String, String>,
    /// `None`: the platform does not measure it.
    pub values: BTreeMap<String, Option<Frozen>>,
}

/// A part's frozen results.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Fixture {
    pub stages: BTreeMap<String, String>,
    pub values: BTreeMap<String, Frozen>,
    /// The platforms it was blessed on, each with what it measured differently.
    pub platforms: BTreeMap<String, Differences>,
}

/// What a run on one platform is checked against.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Expected {
    pub stages: BTreeMap<String, String>,
    pub values: BTreeMap<String, Frozen>,
}

impl Fixture {
    /// The results `platform` is checked against: the fixture's, with its differences laid over.
    pub fn expected(&self, platform: &str) -> Expected {
        let mut e = Expected { stages: self.stages.clone(), values: self.values.clone() };
        if let Some(d) = self.platforms.get(platform) {
            e.stages.extend(d.stages.iter().map(|(k, v)| (k.clone(), v.clone())));
            for (k, v) in &d.values {
                match v {
                    Some(f) => e.values.insert(k.clone(), *f),
                    None => e.values.remove(k),
                };
            }
        }
        e
    }

    /// The fixture of `results` measured on `platform` alone.
    pub fn of(results: &Results, platform: &str) -> Fixture {
        Fixture {
            stages: results.stages.clone(),
            values: results.values.iter().map(|(k, m)| (k.clone(), m.frozen())).collect(),
            platforms: [(platform.to_string(), Differences::default())].into_iter().collect(),
        }
    }

    /// Bless `results` measured on `platform` (see the module docs): the new fixture, or `None`
    /// when this one already holds them.
    pub fn bless(old: Option<&Fixture>, results: &Results, platform: &str) -> Option<Fixture> {
        let Some(old) = old else { return Some(Fixture::of(results, platform)) };
        if old.platforms.contains_key(platform) || old.platforms.is_empty() {
            let unchanged = old.platforms.contains_key(platform) && compare(results, &old.expected(platform)).is_empty();
            return if unchanged { None } else { Some(Fixture::of(results, platform)) };
        }
        // A new platform: record what it measures differently.
        let mut d = Differences::default();
        for (k, now) in &results.stages {
            if old.stages.get(k) != Some(now) {
                d.stages.insert(k.clone(), now.clone());
            }
        }
        for k in results.values.keys().chain(old.values.keys()).collect::<BTreeSet<_>>() {
            match (results.values.get(k), old.values.get(k)) {
                (Some(m), Some(f)) if (m.value() - f.value).abs() <= f.band => {}
                (Some(m), _) => {
                    d.values.insert(k.clone(), Some(m.frozen()));
                }
                (None, Some(_)) => {
                    d.values.insert(k.clone(), None);
                }
                (None, None) => {}
            }
        }
        let mut new = old.clone();
        new.platforms.insert(platform.to_string(), d);
        Some(new)
    }
}

/// Where the part's fixture lives.
pub fn path(corpus_dir: &Path, part: &str) -> PathBuf {
    corpus_dir.join(part).join("fixture.json")
}

fn stage_map(v: Option<&Value>) -> Result<BTreeMap<String, String>, String> {
    let Some(v) = v else { return Ok(BTreeMap::new()) };
    v.as_object()
        .ok_or("`stages` must be an object")?
        .iter()
        .map(|(k, s)| s.as_str().map(|s| (k.clone(), s.to_string())).ok_or_else(|| format!("stage `{k}` must be a string")))
        .collect()
}

fn parse(v: &Value) -> Result<Fixture, String> {
    let values = crate::oracle::fixture_values(v)?;
    let stages = stage_map(Some(v.get("stages").ok_or("no `stages`")?))?;
    let mut platforms = BTreeMap::new();
    for (name, p) in v.get("platforms").and_then(Value::as_object).into_iter().flatten() {
        let mut d = Differences { stages: stage_map(p.get("stages"))?, values: BTreeMap::new() };
        for (k, f) in p.get("values").and_then(Value::as_object).into_iter().flatten() {
            let frozen = match f {
                Value::Null => None,
                f => crate::oracle::fixture_values(&json!({"values": {k: f}}))?.remove(k),
            };
            d.values.insert(k.clone(), frozen);
        }
        platforms.insert(name.clone(), d);
    }
    Ok(Fixture { stages, values, platforms })
}

/// Read a fixture (`None` when there is none).
pub fn read(path: &Path) -> Result<Option<Fixture>, String> {
    if !path.is_file() {
        return Ok(None);
    }
    let text = crate::corpus::read(path).map_err(|e| e.to_string())?;
    let v: Value = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&v).map(Some).map_err(|e| format!("{}: {e}", path.display()))
}

fn frozen_json(f: &Frozen) -> Value {
    // Counts read as whole numbers.
    if f.band == 0.0 && f.value.fract() == 0.0 && f.value.abs() < 2f64.powi(53) {
        json!({"value": f.value as i64, "band": 0})
    } else {
        json!({"value": f.value, "band": f.band})
    }
}

/// The fixture's text.
pub fn render(f: &Fixture) -> Result<String, String> {
    let values: serde_json::Map<String, Value> = f.values.iter().map(|(k, v)| (k.clone(), frozen_json(v))).collect();
    let platforms: serde_json::Map<String, Value> = f
        .platforms
        .iter()
        .map(|(name, d)| {
            let values: serde_json::Map<String, Value> =
                d.values.iter().map(|(k, v)| (k.clone(), v.as_ref().map_or(Value::Null, frozen_json))).collect();
            (name.clone(), json!({"stages": d.stages, "values": values}))
        })
        .collect();
    let text = serde_json::to_string_pretty(&json!({"note": NOTE, "stages": f.stages, "values": values, "platforms": platforms}))
        .map_err(|e| e.to_string())?;
    Ok(text + "\n")
}

/// How `results` differ from what was `expected`: a stage whose outcome changed, a value outside
/// its band, a value measured now but not frozen or frozen but not measured now. Empty when they
/// match.
pub fn compare(results: &Results, expected: &Expected) -> Vec<String> {
    let mut out = Vec::new();
    for s in results.stages.keys().chain(expected.stages.keys()).collect::<BTreeSet<_>>() {
        let (now, frozen) = (results.stages.get(s), expected.stages.get(s));
        if now != frozen {
            out.push(format!("stage {s}: `{}` in the fixture, `{}` now", frozen.map_or("—", String::as_str), now.map_or("—", String::as_str)));
        }
    }
    for k in results.values.keys().chain(expected.values.keys()).collect::<BTreeSet<_>>() {
        match (results.values.get(k), expected.values.get(k)) {
            (Some(m), Some(f)) if (m.value() - f.value).abs() <= f.band => {}
            (Some(m), Some(f)) => out.push(format!("{k}: moved from {} to {} (band {:e})", f.value, m.value(), f.band)),
            (Some(m), None) => out.push(format!("{k}: {} is not in the fixture", m.value())),
            (None, Some(f)) => out.push(format!("{k}: {} in the fixture, not measured now", f.value)),
            (None, None) => {}
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const MAC: &str = "macos-aarch64";
    const LINUX: &str = "linux-x86_64";

    fn results() -> Results {
        Results {
            stages: [("tet".to_string(), "ok".to_string()), ("hex".to_string(), "no: not supported yet".to_string())].into_iter().collect(),
            values: [("tet.tets".to_string(), Metric::Count(120)), ("tet.min_dihedral".to_string(), Metric::Real(21.5))].into_iter().collect(),
        }
    }

    fn round_trip(f: &Fixture, test: &str) -> Fixture {
        let dir = std::env::temp_dir().join(format!("fieldcraft-fixture-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("fixture.json");
        std::fs::write(&p, render(f).unwrap()).unwrap();
        let back = read(&p).unwrap().unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        back
    }

    #[test]
    fn a_blessed_fixture_matches_its_results() {
        let r = results();
        let f = round_trip(&Fixture::bless(None, &r, MAC).unwrap(), "blessed");
        assert_eq!(compare(&r, &f.expected(MAC)), Vec::<String>::new());
        assert_eq!(Fixture::bless(Some(&f), &r, MAC), None, "blessing the same results again changes nothing");
    }

    /// Counts and outcomes are exact; reals move by their band only; a value that appears or
    /// disappears is a change.
    #[test]
    fn changes_are_found() {
        let f = round_trip(&Fixture::of(&results(), MAC), "changes").expected(MAC);
        let mut r = results();
        r.values.insert("tet.tets".into(), Metric::Count(121));
        r.values.insert("tet.min_dihedral".into(), Metric::Real(21.5 * (1.0 + 1e-10)));
        assert_eq!(compare(&r, &f).len(), 1, "{:?}", compare(&r, &f));
        r.values.insert("tet.min_dihedral".into(), Metric::Real(21.5 * (1.0 + 1e-8)));
        r.stages.insert("hex".into(), "ok: extrude".into());
        r.values.insert("hex.elements".into(), Metric::Count(10));
        r.values.remove("tet.tets");
        let got = compare(&r, &f);
        assert_eq!(got.len(), 4, "{got:?}");
        assert!(got.iter().any(|g| g.starts_with("stage hex: `no: not supported yet` in the fixture, `ok: extrude` now")), "{got:?}");
    }

    /// A second platform's differences are recorded beside the first's, each platform is checked
    /// against its own, and a change blessed on one platform drops the other's stale results.
    #[test]
    fn platforms_keep_their_own_results() {
        let mac = results();
        let mut linux = results();
        linux.values.insert("tet.tets".into(), Metric::Count(130));
        linux.stages.insert("hex".into(), "ok: sweep".into());
        linux.values.insert("hex.elements".into(), Metric::Count(40));
        linux.values.remove("tet.min_dihedral");
        let f = Fixture::bless(None, &mac, MAC).unwrap();
        let f = round_trip(&Fixture::bless(Some(&f), &linux, LINUX).unwrap(), "platforms");
        assert_eq!(f.platforms.len(), 2);
        let d = &f.platforms[LINUX];
        assert_eq!(d.stages.len(), 1);
        assert_eq!(d.values.len(), 3, "{d:?}");
        assert_eq!(d.values["tet.min_dihedral"], None);
        assert!(compare(&mac, &f.expected(MAC)).is_empty());
        assert!(compare(&linux, &f.expected(LINUX)).is_empty());
        assert_eq!(compare(&linux, &f.expected(MAC)).len(), 4, "a stage, a moved, a new and a missing value");
        // Unchanged on either platform: nothing to bless.
        assert_eq!(Fixture::bless(Some(&f), &linux, LINUX), None);
        assert_eq!(Fixture::bless(Some(&f), &mac, MAC), None);
        // A change on the Mac starts over from the Mac's results.
        let mut moved = mac.clone();
        moved.values.insert("tet.tets".into(), Metric::Count(125));
        let g = Fixture::bless(Some(&f), &moved, MAC).unwrap();
        assert_eq!(g.platforms.keys().collect::<Vec<_>>(), vec![MAC]);
        assert!(compare(&moved, &g.expected(MAC)).is_empty());
    }

    #[test]
    fn bands() {
        assert!(Metric::Real(0.0).band() > 0.0);
        assert_eq!(Metric::Count(0).band(), 0.0);
        assert_eq!(Metric::Relative(1e-15).band(), FIXTURE_RELATIVE);
        assert!(platform().contains('-'));
    }

    #[test]
    fn hostile_fixtures_are_an_error() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-fixture-hostile-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("fixture.json");
        for text in [
            "",
            "{",
            r#"{"values": {}}"#,
            r#"{"stages": {"a": 1}, "values": {}}"#,
            r#"{"stages": {}, "values": {"a": {"value": 1}}}"#,
            r#"{"stages": {}, "values": {}, "platforms": {"x": {"values": {"a": {"band": 1}}}}}"#,
            r#"{"stages": {}, "values": {}, "platforms": {"x": {"stages": {"a": null}}}}"#,
        ] {
            std::fs::write(&p, text).unwrap();
            assert!(read(&p).is_err(), "{text}");
        }
        assert_eq!(read(&dir.join("none.json")).unwrap(), None);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
