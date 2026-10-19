//! The verification and validation oracle (docs/verification.md): every case in
//! `tests/oracle/<case>/` is our own model of a problem with a known answer, run through the
//! command engine and compared with its reference values and with its frozen fixture.
//!
//! A case directory holds:
//! - `expected.json`: the case's kind (theory, experiment, benchmark), source, status and its
//!   quantities: how each is measured, the reference value, the tolerance and the reason for it;
//! - `recipe.json`: the model as FieldCraft commands (a JSON array of `{"command", "params"}`,
//!   each step optionally with an `id` that measures can read its output by); absent while a
//!   case is pending;
//! - `design.json` (optional): a SolidCraft script that builds the geometry, opened first;
//! - `fixture.json`: our accepted values, written by `--bless` once the case passes;
//! - `README.md`: the problem in our own words.
//!
//! A case marked `nightly` (too long for every change) runs only with `--nightly` (the scheduled,
//! manual and release runs) or when named with `--case` or `--bless`; otherwise it is checked by
//! its fixture alone. docs/oracle.md always shows a nightly case's fixture, so the page reads the
//! same whether or not the run included it.
//!
//! [`run`] (the `oracle` binary, `cargo xtask oracle`) runs the cases, prints and checks the
//! results and writes `docs/oracle.md` after a full run.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use fieldcraft_engine::Session;
use serde_json::{Value, json};

/// Fixture band: a result may move by this fraction of its size (and a thousandth of the
/// reference's absolute tolerance) before CI calls it a change. Results are deterministic on one
/// platform; the band only absorbs last-digit differences of the system maths library.
pub const FIXTURE_RELATIVE: f64 = 1e-9; // unit-factors: relative to the result
/// The widest band a case may ask for (`fixture_relative`).
const MAX_FIXTURE_RELATIVE: f64 = 1e-3; // unit-factors: relative to the result
/// The share of the reference's absolute tolerance a fixture band adds.
const ABSOLUTE_SHARE: f64 = 1e-3; // unit-factors: a fraction of a tolerance
/// Magnitudes the page writes as plain decimals.
const PLAIN: std::ops::Range<f64> = 1e-3..1e5; // unit-factors: magnitudes of a written number
/// Bounds below this are written in scientific notation.
const SMALL_BOUND: f64 = 1e-3; // unit-factors: a magnitude of a written number
/// An error's floor is this many times its band, on the decade above.
const FLOOR_OVER_BAND: f64 = 1000.0; // unit-factors: a ratio
/// The slack that keeps a decade's own bound on it (1000 · 1e-9 is 1.0000000000000002e-6).
const DECADE_SLACK: f64 = 1e-9; // unit-factors: a fraction of a decade's logarithm

/// The kinds of reference (docs/verification.md).
pub const KINDS: [&str; 3] = ["theory", "benchmark", "experiment"];

/// Options of one run.
#[derive(Debug, Default, Clone)]
pub struct Options {
    /// Run only these cases (all when empty).
    pub cases: Vec<String>,
    /// Write `fixture.json` for these cases when they pass (`all` for every passing case).
    pub bless: Vec<String>,
    /// Skip cases marked `slow` (and `nightly`).
    pub quick: bool,
    /// Run the cases marked `nightly` too.
    pub nightly: bool,
    /// Write docs/oracle.md (only after a run of every case).
    pub write: bool,
}

/// A quantity's tolerance: |ours − reference| ≤ relative·|reference| + absolute.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Tolerance {
    pub relative: f64,
    pub absolute: f64,
}

impl Tolerance {
    fn from_json(v: &Value) -> Result<Tolerance, String> {
        let get = |k: &str| -> Result<f64, String> {
            match v.get(k) {
                None => Ok(0.0),
                Some(x) => x.as_f64().filter(|x| x.is_finite() && *x >= 0.0).ok_or_else(|| format!("tolerance `{k}` must be a number ≥ 0")),
            }
        };
        let t = Tolerance { relative: get("relative")?, absolute: get("absolute")? };
        if t.relative == 0.0 && t.absolute == 0.0 && v.get("exact") != Some(&Value::Bool(true)) {
            return Err("a tolerance gives `relative`, `absolute` or `\"exact\": true`".into());
        }
        Ok(t)
    }

    pub fn allows(&self, reference: f64, ours: f64) -> bool {
        (ours - reference).abs() <= self.relative * reference.abs() + self.absolute
    }

    fn describe(&self) -> String {
        match (self.relative > 0.0, self.absolute > 0.0) {
            (true, false) => fraction(self.relative),
            (false, true) => format!("± {}", number(self.absolute)),
            (true, true) => format!("{} + {}", fraction(self.relative), number(self.absolute)),
            (false, false) => "exact".into(),
        }
    }
}

/// One quantity of a case.
#[derive(Debug, Clone)]
pub struct Quantity {
    pub id: String,
    pub description: String,
    pub measure: Value,
    pub reference: Option<f64>,
    pub tolerance: Tolerance,
    pub reason: String,
    pub known_failure: Option<String>,
    /// The fixture band's relative part, when the quantity is more sensitive to rounding than
    /// [`FIXTURE_RELATIVE`] allows (an ill-conditioned problem); `fixture_reason` says why.
    pub fixture_relative: f64,
    /// Whether the result is frozen in `fixture.json`. A quantity that depends on a mesh the
    /// mesher generates from geometry (`"fixture": false`, with a `fixture_reason`) is not: the
    /// mesh is not bit-identical across platforms, so only its reference tolerance gates it.
    pub frozen: bool,
}

/// One case's `expected.json`.
#[derive(Debug, Clone)]
pub struct Expected {
    pub title: String,
    pub kind: String,
    pub source: String,
    pub pending: Option<String>,
    pub slow: bool,
    /// Runs only nightly (`--nightly`); the page shows its fixture.
    pub nightly: bool,
    pub quantities: Vec<Quantity>,
}

fn text(v: &Value, key: &str) -> Result<String, String> {
    v.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .map(str::to_string)
        .ok_or_else(|| format!("`{key}` must be a non-empty string"))
}

impl Expected {
    pub fn parse(v: &Value) -> Result<Expected, String> {
        let kind = text(v, "kind")?;
        if !KINDS.contains(&kind.as_str()) {
            return Err(format!("`kind` must be one of {KINDS:?}"));
        }
        let pending = match v.get("status").and_then(Value::as_str) {
            None | Some("active") => None,
            Some("pending") => Some(text(v, "requires")?),
            Some(other) => return Err(format!("`status` must be active or pending, not `{other}`")),
        };
        let list = v.get("quantities").and_then(Value::as_array).filter(|q| !q.is_empty()).ok_or("`quantities` must be a non-empty list")?;
        let mut quantities = Vec::with_capacity(list.len());
        for q in list {
            let id = text(q, "id")?;
            let context = |e: String| format!("quantity `{id}`: {e}");
            let reference = match q.get("value") {
                None | Some(Value::Null) if pending.is_some() => None,
                Some(x) => Some(x.as_f64().filter(|x| x.is_finite()).ok_or_else(|| context("`value` must be a finite number".into()))?),
                None => return Err(context("`value` is required".into())),
            };
            let tolerance = match q.get("tolerance") {
                None if pending.is_some() => Tolerance::default(),
                Some(t) => Tolerance::from_json(t).map_err(context)?,
                None => return Err(context("`tolerance` is required".into())),
            };
            let reason = if pending.is_some() {
                q.get("reason").and_then(Value::as_str).unwrap_or_default().to_string()
            } else {
                text(q, "reason").map_err(context)?
            };
            let measure = q.get("measure").cloned().unwrap_or(Value::Null);
            if pending.is_none() && measure.get("pointer").and_then(Value::as_str).is_none() {
                return Err(context("`measure` needs a JSON `pointer` into a step's output or a command's".into()));
            }
            if quantities.iter().any(|x: &Quantity| x.id == id) {
                return Err(format!("quantity `{id}` is listed twice"));
            }
            let frozen = match q.get("fixture") {
                None => true,
                Some(Value::Bool(f)) => *f,
                Some(_) => return Err(context("`fixture` must be true or false".into())),
            };
            if !frozen {
                text(q, "fixture_reason").map_err(|_| context("`\"fixture\": false` needs a `fixture_reason`".into()))?;
                if q.get("fixture_relative").is_some() {
                    return Err(context("`fixture_relative` has no meaning with `\"fixture\": false`".into()));
                }
            }
            let fixture_relative = match q.get("fixture_relative") {
                None => FIXTURE_RELATIVE,
                Some(r) => {
                    let r = r
                        .as_f64()
                        .filter(|r| *r > 0.0 && *r <= MAX_FIXTURE_RELATIVE)
                        .ok_or_else(|| context("`fixture_relative` must be in (0, 1e-3]".into()))?;
                    text(q, "fixture_reason").map_err(|_| context("`fixture_relative` needs a `fixture_reason`".into()))?;
                    r
                }
            };
            quantities.push(Quantity {
                fixture_relative,
                frozen,
                description: text(q, "quantity").map_err(context)?,
                known_failure: q.get("known_failure").and_then(Value::as_str).map(str::to_string),
                id,
                measure,
                reference,
                tolerance,
                reason,
            });
        }
        Ok(Expected {
            title: text(v, "title")?,
            kind,
            source: text(v, "source")?,
            pending,
            slow: v.get("slow").and_then(Value::as_bool).unwrap_or(false),
            nightly: v.get("nightly").and_then(Value::as_bool).unwrap_or(false),
            quantities,
        })
    }
}

/// A frozen value and how far it may move.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frozen {
    pub value: f64,
    pub band: f64,
}

/// The band a fixture value gets when it is blessed. A quantity whose reference is zero (a
/// rigid-body frequency, an equilibrium residual, an invalid-element count) is rounding noise
/// inside its tolerance: its band is that tolerance, as freezing noise would only flag rounding.
pub fn fixture_band(value: f64, reference: f64, tolerance: Tolerance, relative: f64) -> f64 {
    if reference == 0.0 {
        return tolerance.absolute;
    }
    relative * value.abs().max(reference.abs()) + ABSOLUTE_SHARE * tolerance.absolute
}

fn read_json(path: &Path) -> Result<Value, String> {
    let text = crate::corpus::read(path).map_err(|e| e.to_string())?;
    serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))
}

fn read_fixture(path: &Path) -> Result<Option<BTreeMap<String, Frozen>>, String> {
    if !path.is_file() {
        return Ok(None);
    }
    let v = read_json(path)?;
    let values = v.get("values").and_then(Value::as_object).ok_or_else(|| format!("{}: no `values`", path.display()))?;
    values
        .iter()
        .map(|(k, f)| {
            let value = f.get("value").and_then(Value::as_f64);
            let band = f.get("band").and_then(Value::as_f64).filter(|b| *b >= 0.0);
            match (value, band) {
                (Some(value), Some(band)) => Ok((k.clone(), Frozen { value, band })),
                _ => Err(format!("{}: `{k}` needs a `value` and a `band` ≥ 0", path.display())),
            }
        })
        .collect::<Result<_, _>>()
        .map(Some)
}

/// Run a case's recipe and measure its quantities, in order. A failing command fails every
/// quantity.
pub fn run_case(dir: &Path, expected: &Expected) -> Result<Vec<Result<f64, String>>, String> {
    let mut session = Session::default();
    let design = dir.join("design.json");
    if design.is_file() {
        let script = crate::corpus::read(&design).map_err(|e| e.to_string())?;
        let bytes = fieldcraft_geometry::script_design(&script).map_err(|e| format!("design.json: {e}"))?;
        let opened = fieldcraft_geometry::open_design(&bytes).map_err(|e| format!("design.json: {e}"))?;
        if let Some(e) = opened.errors.first() {
            return Err(format!("design.json: feature {}: {}", e.feature, e.message));
        }
        session.design = Some(Arc::new(opened));
    }
    let recipe = read_json(&dir.join("recipe.json"))?;
    let results = session.run_script(&recipe).map_err(|e| format!("recipe.json: {e}"))?;
    let outputs: BTreeMap<String, Value> =
        recipe.as_array().into_iter().flatten().zip(results).filter_map(|(step, out)| Some((step.get("id")?.as_str()?.to_string(), out))).collect();
    Ok(expected.quantities.iter().map(|q| measure(&mut session, &outputs, &q.measure)).collect())
}

/// Read one number: from a recipe step's output (`"from": id`) or a command run now
/// (`"command"`, `"params"`), at a JSON pointer, times `scale`; with `"over"` (a measure of its
/// own), divided by that measure's number.
fn measure(session: &mut Session, outputs: &BTreeMap<String, Value>, m: &Value) -> Result<f64, String> {
    let pointer = m.get("pointer").and_then(Value::as_str).ok_or("no `pointer`")?;
    let fresh;
    let out = match (m.get("from").and_then(Value::as_str), m.get("command").and_then(Value::as_str)) {
        (Some(id), None) => outputs.get(id).ok_or_else(|| format!("no recipe step has id `{id}`"))?,
        (None, Some(command)) => {
            let params = m.get("params").cloned().unwrap_or_else(|| json!({}));
            fresh = session.execute(command, &params).map_err(|e| format!("{command}: {e}"))?;
            &fresh
        }
        _ => return Err("a measure gives either `from` (a recipe step id) or `command`".into()),
    };
    let v = out.pointer(pointer).ok_or_else(|| format!("nothing at `{pointer}`"))?;
    let x = v.as_f64().ok_or_else(|| format!("`{pointer}` is {v}, not a number"))?;
    let mut x = x * m.get("scale").and_then(Value::as_f64).unwrap_or(1.0);
    // A ratio: this value over another measure's (both from the same run).
    if let Some(over) = m.get("over") {
        let y = measure(session, outputs, over)?;
        if y == 0.0 {
            return Err(format!("`{pointer}` is measured over a zero"));
        }
        x /= y;
    }
    if x.is_finite() { Ok(x) } else { Err(format!("`{pointer}` is not finite")) }
}

/// How one quantity came out.
#[derive(Debug, Clone, PartialEq)]
pub enum Verdict {
    Pass,
    Fail,
    /// Fails, and expected.json says why.
    Known,
    /// Marked as a known failure but passes: remove the mark.
    Fixed,
    /// Could not be measured.
    Error(String),
    Pending,
    /// Not run here (a nightly case) and not kept in its fixture.
    Nightly,
}

/// The fixture check of one quantity.
#[derive(Debug, Clone, PartialEq)]
pub enum FixtureCheck {
    Matches,
    Moved {
        frozen: f64,
    },
    Missing,
    /// No fixture is kept (pending, failing or known-failing quantity).
    NotKept,
    Blessed,
    /// A nightly case's fixture, shown without running the case.
    Nightly,
}

/// One row of the report.
#[derive(Debug, Clone)]
pub struct Row {
    pub case: String,
    pub quantity: Quantity,
    pub ours: Option<f64>,
    pub verdict: Verdict,
    pub fixture: FixtureCheck,
}

/// Judge a measured value against the reference and the fixture.
pub fn judge(q: &Quantity, ours: &Result<f64, String>, frozen: Option<&Frozen>) -> (Verdict, FixtureCheck) {
    let (Some(reference), Ok(ours)) = (q.reference, ours) else {
        return match ours {
            Err(e) => (Verdict::Error(e.clone()), FixtureCheck::NotKept),
            Ok(_) => (Verdict::Pending, FixtureCheck::NotKept),
        };
    };
    let pass = q.tolerance.allows(reference, *ours);
    let verdict = match (pass, &q.known_failure) {
        (true, None) => Verdict::Pass,
        (true, Some(_)) => Verdict::Fixed,
        (false, None) => Verdict::Fail,
        (false, Some(_)) => Verdict::Known,
    };
    let fixture = match frozen {
        _ if !q.frozen => FixtureCheck::NotKept,
        Some(f) if (ours - f.value).abs() <= f.band => FixtureCheck::Matches,
        Some(f) => FixtureCheck::Moved { frozen: f.value },
        None if verdict == Verdict::Pass => FixtureCheck::Missing,
        None => FixtureCheck::NotKept,
    };
    (verdict, fixture)
}

/// The case directories under `tests/oracle`, sorted.
pub fn case_dirs(root: &Path) -> Result<Vec<PathBuf>, String> {
    let dir = root.join("tests").join("oracle");
    let mut cases: Vec<PathBuf> = std::fs::read_dir(&dir)
        .map_err(|e| format!("{}: {e}", dir.display()))?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.join("expected.json").is_file())
        .collect();
    cases.sort();
    Ok(cases)
}

fn name_of(dir: &Path) -> String {
    dir.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default()
}

/// Run the oracle. Fails when a quantity fails its reference, a fixture moved or is missing,
/// a known failure now passes, or a case cannot run.
pub fn run(root: &Path, options: &Options) -> Result<(), String> {
    let all = case_dirs(root)?;
    for want in options.cases.iter().chain(options.bless.iter().filter(|b| *b != "all")) {
        if !all.iter().any(|d| name_of(d) == *want) {
            return Err(format!("no case `{want}` in tests/oracle"));
        }
    }
    let selected: Vec<&PathBuf> = all.iter().filter(|d| options.cases.is_empty() || options.cases.contains(&name_of(d))).collect();
    let mut rows = Vec::new();
    // What docs/oracle.md shows: every row as run, but a nightly case's from its fixture.
    let mut page_rows = Vec::new();
    let mut rows_of_skipped = Vec::new();
    let mut cases: Vec<(String, Expected)> = Vec::new();
    let mut problems = Vec::new();
    let mut skipped = Vec::new();
    for dir in selected {
        let case = name_of(dir);
        let expected = Expected::parse(&read_json(&dir.join("expected.json"))?).map_err(|e| format!("{case}/expected.json: {e}"))?;
        if (expected.slow || expected.nightly) && options.quick {
            skipped.push(case);
            continue;
        }
        let fixture_path = dir.join("fixture.json");
        let fixture = read_fixture(&fixture_path)?;
        let named = options.cases.contains(&case) || options.bless.contains(&case);
        let nightly = expected.nightly && expected.pending.is_none();
        if nightly && !(options.nightly || named) {
            eprintln!("oracle: {case} runs nightly (`--nightly`); checked by its fixture");
            let rows = fixture_rows(&case, &expected, fixture.as_ref());
            problems.extend(nightly_problems(&case, &rows));
            page_rows.extend(rows.iter().cloned());
            rows_of_skipped.extend(rows);
            cases.push((case, expected));
            continue;
        }
        let mut case_rows: Vec<Row> = if expected.pending.is_some() {
            let row =
                |q: &Quantity| Row { case: case.clone(), quantity: q.clone(), ours: None, verdict: Verdict::Pending, fixture: FixtureCheck::NotKept };
            expected.quantities.iter().map(row).collect()
        } else {
            let start = Instant::now();
            let measured = run_case(dir, &expected).unwrap_or_else(|e| vec![Err(e); expected.quantities.len()]);
            eprintln!("oracle: {case} ({:.2} s)", start.elapsed().as_secs_f64());
            expected
                .quantities
                .iter()
                .zip(&measured)
                .map(|(q, m)| {
                    let (verdict, fixture) = judge(q, m, fixture.as_ref().and_then(|f| f.get(&q.id)));
                    Row { case: case.clone(), quantity: q.clone(), ours: m.as_ref().ok().copied(), verdict, fixture }
                })
                .collect()
        };
        let bless = options.bless.iter().any(|b| *b == case || b == "all");
        if bless && expected.pending.is_none() {
            let passing = case_rows.iter().all(|r| matches!(r.verdict, Verdict::Pass | Verdict::Known));
            if !passing {
                problems.push(format!("{case}: not blessed, it does not pass against expected.json"));
            } else {
                let mut values = serde_json::Map::new();
                for r in case_rows.iter_mut().filter(|r| r.verdict == Verdict::Pass && r.quantity.frozen) {
                    if let (Some(ours), Some(reference)) = (r.ours, r.quantity.reference) {
                        values.insert(
                            r.quantity.id.clone(),
                            json!({"value": ours, "band": fixture_band(ours, reference, r.quantity.tolerance, r.quantity.fixture_relative)}),
                        );
                        r.fixture = FixtureCheck::Blessed;
                    }
                }
                let text = serde_json::to_string_pretty(&json!({
                    "note": "FieldCraft's accepted results for this case, written by `cargo xtask oracle --bless`. A result outside its band fails CI.",
                    "values": values,
                }))
                .map_err(|e| e.to_string())?;
                std::fs::write(&fixture_path, text + "\n").map_err(|e| format!("{}: {e}", fixture_path.display()))?;
            }
        }
        for r in &case_rows {
            let id = format!("{case}/{}", r.quantity.id);
            match &r.verdict {
                Verdict::Fail => {
                    problems.push(format!("{id}: {} vs reference {} ({})", show(r.ours), show(r.quantity.reference), r.quantity.tolerance.describe()))
                }
                Verdict::Fixed => problems.push(format!("{id}: passes now; remove its `known_failure` mark")),
                Verdict::Error(e) => problems.push(format!("{id}: {e}")),
                _ => {}
            }
            match &r.fixture {
                FixtureCheck::Moved { frozen } => {
                    problems.push(format!("{id}: moved from the fixture {frozen:e} to {:e}", r.ours.unwrap_or(f64::NAN)))
                }
                FixtureCheck::Missing => problems.push(format!("{id}: no fixture; bless the case (`cargo xtask oracle --bless {case}`)")),
                _ => {}
            }
        }
        // A nightly case's page shows its fixture (as blessed by this run, if it was).
        if nightly {
            page_rows.extend(fixture_rows(&case, &expected, read_fixture(&fixture_path)?.as_ref()));
        } else {
            page_rows.extend(case_rows.iter().cloned());
        }
        rows.extend(case_rows);
        cases.push((case, expected));
    }
    print_summary(&rows);
    print_summary(&rows_of_skipped);
    if options.write && options.cases.is_empty() && !options.quick {
        let path = root.join("docs").join("oracle.md");
        std::fs::write(&path, report(&cases, &page_rows)).map_err(|e| format!("{}: {e}", path.display()))?;
        eprintln!("oracle: wrote docs/oracle.md");
    }
    if !skipped.is_empty() {
        eprintln!("oracle: skipped slow cases: {}", skipped.join(", "));
    }
    if problems.is_empty() { Ok(()) } else { Err(format!("{} problem(s):\n  {}", problems.len(), problems.join("\n  "))) }
}

/// What fails a nightly case that was not run: a fixture that no longer meets its reference (an
/// edited expected.json), and a quantity that should be frozen but has no fixture (a case added
/// without being blessed), as for any other case.
fn nightly_problems(case: &str, rows: &[Row]) -> Vec<String> {
    rows.iter()
        .filter_map(|r| {
            let id = &r.quantity.id;
            match (&r.verdict, &r.fixture) {
                (Verdict::Fail, _) => Some(format!("{case}/{id}: its fixture misses the reference; run it (`cargo xtask oracle --case {case}`)")),
                (_, FixtureCheck::Missing) => Some(format!("{case}/{id}: no fixture; bless the case (`cargo xtask oracle --bless {case}`)")),
                _ => None,
            }
        })
        .collect()
}

/// A nightly case's rows from its fixture: a frozen quantity judged by its fixture value, a known
/// failure as known, a frozen quantity without a fixture as missing, anything else as nightly
/// only.
fn fixture_rows(case: &str, expected: &Expected, fixture: Option<&BTreeMap<String, Frozen>>) -> Vec<Row> {
    expected
        .quantities
        .iter()
        .map(|q| {
            let frozen = fixture.and_then(|f| f.get(&q.id)).filter(|_| q.frozen);
            let (ours, verdict, check) = match (frozen, &q.known_failure) {
                (Some(f), _) => (Some(f.value), judge(q, &Ok(f.value), Some(f)).0, FixtureCheck::Nightly),
                (None, Some(_)) => (None, Verdict::Known, FixtureCheck::NotKept),
                (None, None) if q.frozen => (None, Verdict::Nightly, FixtureCheck::Missing),
                (None, None) => (None, Verdict::Nightly, FixtureCheck::NotKept),
            };
            Row { case: case.to_string(), quantity: q.clone(), ours, verdict, fixture: check }
        })
        .collect()
}

fn print_summary(rows: &[Row]) {
    for r in rows {
        if r.verdict == Verdict::Pending {
            continue;
        }
        println!(
            "{:<34} {:<22} ours {:<14} ref {:<14} err {:<12} {:<14} {}",
            r.case,
            r.quantity.id,
            show(r.ours),
            show(r.quantity.reference),
            error_of(r),
            verdict_text(&r.verdict),
            fixture_text(&r.fixture)
        );
    }
}

/// A number with five significant digits. The notation follows the rounded value, so a number a
/// rounding error below a decade (0.000999999999) prints like the decade itself (0.001).
pub fn number(x: f64) -> String {
    if x == 0.0 || !x.is_finite() {
        return if x == 0.0 { "0".into() } else { x.to_string() };
    }
    // Round to five significant digits first; `{:e}` output always parses back.
    let x = format!("{x:.4e}").parse::<f64>().unwrap_or(x);
    let a = x.abs();
    if PLAIN.contains(&a) {
        let digits = 4 - a.log10().floor() as i32;
        let s = format!("{x:.*}", digits.max(0) as usize);
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s }
    } else {
        format!("{x:.4e}")
    }
}

/// A fraction as a percentage with three significant digits.
fn percent(r: f64) -> String {
    let p = r * 100.0;
    if p == 0.0 {
        return "0%".into();
    }
    let digits = (2 - p.abs().log10().floor() as i32).max(0) as usize;
    let s = format!("{p:.digits$}");
    let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s };
    format!("{s}%")
}

/// A tolerance's relative part: a percentage, or a plain fraction when it only admits rounding.
fn fraction(r: f64) -> String {
    if r >= 1e-4 { percent(r) } else { format!("{r:.0e} rel.") }
}

fn show(x: Option<f64>) -> String {
    x.map(number).unwrap_or_else(|| "—".into())
}

/// Below this fraction of its absolute tolerance, a quantity whose reference is zero is shown as
/// zero: it is rounding noise, and its digits change with any reordering of the arithmetic
/// (docs/oracle.md must not).
const NOISE: f64 = 1e-2;

fn is_noise(r: &Row, value: f64) -> bool {
    r.quantity.reference == Some(0.0) && value.abs() < NOISE * r.quantity.tolerance.absolute
}

/// Our value as the report shows it.
fn ours_text(r: &Row) -> String {
    match r.ours {
        Some(v) if is_noise(r, v) => "≈ 0".into(),
        other => show(other),
    }
}

/// docs/oracle.md's result column. A quantity that is not frozen depends on a generated mesh
/// that differs by platform, so the page names that instead of one platform's number (the run's
/// console prints the value).
fn page_ours(r: &Row) -> String {
    if r.quantity.frozen { ours_text(r) } else { "mesh-dependent".into() }
}

/// docs/oracle.md's error column; see [`page_ours`].
fn page_error(r: &Row) -> String {
    if r.quantity.frozen { error_of(r) } else { format!("within {}", r.quantity.tolerance.describe()) }
}

/// The error against the reference: relative where the reference is not zero.
///
/// A frozen result may move by its fixture band (`fixture_relative` of its size) without being a
/// change, so digits of the error at that scale are rounding noise and must not reach
/// docs/oracle.md. The rule: the page's floor is 1000 bands (rounded up to a power of ten; the
/// band's own size, not the result, sets it, so the floor cannot move with the result); an error
/// under the floor is shown as that bound ("< 1e-4%" for the usual 10⁻⁹ band) and a larger one is
/// rounded to a tenth of the floor. A quantity whose reference is zero keeps the "≈ 0" rule of
/// [`is_noise`].
fn error_of(r: &Row) -> String {
    let (Some(ours), Some(reference)) = (r.ours, r.quantity.reference) else { return "—".into() };
    let d = ours - reference;
    if reference == 0.0 {
        return if is_noise(r, d) || d == 0.0 { "≈ 0".into() } else { format!("{d:.1e}") };
    }
    let band = r.quantity.fixture_relative;
    if r.quantity.tolerance.relative > 0.0 {
        match banded(d / reference.abs(), band) {
            Banded::Below(bound) => format!("< {}%", decade(bound * 100.0)),
            Banded::Error(rel) => percent(rel),
        }
    } else {
        match banded(d, band * reference.abs()) {
            Banded::Below(bound) => format!("< {}", decade(bound)),
            Banded::Error(d) => format!("{d:.1e}"),
        }
    }
}

/// An error as [`error_of`] shows it.
#[derive(Debug, PartialEq)]
enum Banded {
    /// Smaller than this bound: below the page's floor.
    Below(f64),
    /// Rounded to the page's resolution.
    Error(f64),
}

/// `error` against the floor of `band` (both in the same unit); see [`error_of`].
fn banded(error: f64, band: f64) -> Banded {
    if !(band > 0.0 && band.is_finite() && error.is_finite()) {
        return Banded::Error(error);
    }
    let exponent = ((FLOOR_OVER_BAND * band).log10() - DECADE_SLACK).ceil() as i32;
    let bound = 10f64.powi(exponent);
    let step = 10f64.powi(exponent - 1);
    if error.abs() < bound { Banded::Below(bound) } else { Banded::Error((error / step).round() * step) }
}

/// A power of ten as the page writes a bound: "1e-7" when small, a plain number otherwise.
fn decade(x: f64) -> String {
    if x < SMALL_BOUND { format!("{x:.0e}") } else { number(x) }
}

fn verdict_text(v: &Verdict) -> &'static str {
    match v {
        Verdict::Pass => "pass",
        Verdict::Fail => "**FAIL**",
        Verdict::Known => "known failure",
        Verdict::Fixed => "**passes: remove known_failure**",
        Verdict::Error(_) => "**ERROR**",
        Verdict::Pending => "pending",
        Verdict::Nightly => "nightly",
    }
}

fn fixture_text(f: &FixtureCheck) -> &'static str {
    match f {
        FixtureCheck::Matches => "matches",
        FixtureCheck::Moved { .. } => "**MOVED**",
        FixtureCheck::Missing => "**missing**",
        FixtureCheck::NotKept => "—",
        FixtureCheck::Blessed => "blessed",
        FixtureCheck::Nightly => "nightly",
    }
}

/// The fixture column of docs/oracle.md. A value blessed in this run matches its new fixture by
/// construction, so the page reads the same whether or not the run blessed anything; otherwise
/// a page written by a blessing run would differ from the one CI regenerates.
fn page_fixture_text(f: &FixtureCheck) -> &'static str {
    match f {
        FixtureCheck::Blessed => fixture_text(&FixtureCheck::Matches),
        other => fixture_text(other),
    }
}

fn cell(s: &str) -> String {
    s.replace('|', "\\|").replace('\n', " ")
}

/// docs/oracle.md.
pub fn report(cases: &[(String, Expected)], rows: &[Row]) -> String {
    let active: Vec<&(String, Expected)> = cases.iter().filter(|(_, e)| e.pending.is_none()).collect();
    let count = |kind: &str| active.iter().filter(|(_, e)| e.kind == kind).count();
    let quantities = rows.iter().filter(|r| r.verdict != Verdict::Pending).count();
    let passing = rows.iter().filter(|r| r.verdict == Verdict::Pass).count();
    let known = rows.iter().filter(|r| r.verdict == Verdict::Known).count();
    let mut md = String::new();
    let _ = write!(
        md,
        "# Verification and validation oracle\n\n\
Generated by `cargo xtask oracle` from [`tests/oracle/`](../tests/oracle/); the policy is [verification.md](verification.md). \
Every case is our own model, built with FieldCraft commands from the problem's published description, and compared with \
the reference value its `expected.json` cites, within the stated tolerance and for the stated reason. Accepted results are \
frozen in each case's `fixture.json`: a later change that moves one beyond its band (10⁻⁹ of its size) fails CI even while \
the reference still passes. Errors are rounded to 100 bands and an error under 1000 bands \
(rounded up to a power of ten) is shown as that bound: digits at the band's scale are rounding noise. No reference value was produced with a commercial product.\n\n\
**{passing} of {quantities} quantities pass** ({known} known failure(s)) in {} active cases: {} theory, {} published \
benchmark, {} experiment; {} pending cases wait for capabilities still being built.\n\n",
        active.len(),
        count("theory"),
        count("benchmark"),
        count("experiment"),
        cases.len() - active.len(),
    );
    for kind in KINDS {
        let in_kind: Vec<&&(String, Expected)> = active.iter().filter(|(_, e)| e.kind == kind).collect();
        if in_kind.is_empty() {
            continue;
        }
        let heading = match kind {
            "theory" => "Theory: closed-form solutions, patch tests and mesh invariants",
            "benchmark" => "Published benchmarks",
            _ => "Experiments",
        };
        let _ = write!(
            md,
            "## {heading}\n\n| Case | Quantity | Reference | FieldCraft | Error | Tolerance | Result | Fixture |\n|---|---|---|---|---|---|---|---|\n"
        );
        for (case, _) in &in_kind {
            for r in rows.iter().filter(|r| &r.case == case) {
                let mut result = verdict_text(&r.verdict).to_string();
                if let Some(k) = &r.quantity.known_failure {
                    result = format!("{result}: {k}");
                }
                if let Verdict::Error(e) = &r.verdict {
                    result = format!("{result}: {e}");
                }
                let _ = writeln!(
                    md,
                    "| [{case}](../tests/oracle/{case}/) | {} | {} | {} | {} | {} | {} | {} |",
                    cell(&r.quantity.description),
                    show(r.quantity.reference),
                    page_ours(r),
                    page_error(r),
                    r.quantity.tolerance.describe(),
                    cell(&result),
                    page_fixture_text(&r.fixture),
                );
            }
        }
        md.push('\n');
    }
    md += "## Sources and tolerances\n\n";
    for (case, e) in &active {
        let _ = write!(md, "- **[{case}](../tests/oracle/{case}/)** ({}): {}. Source: {}.", e.kind, cell(&e.title), cell(&e.source));
        if e.nightly {
            md += " Runs nightly (`--nightly`); the table shows its fixture.";
        }
        let reasons: Vec<String> = e.quantities.iter().map(|q| format!("{}: {}", q.id, cell(&q.reason))).collect();
        let _ = writeln!(md, " Tolerances — {}", reasons.join("; "));
    }
    let pending: Vec<&(String, Expected)> = cases.iter().filter(|(_, e)| e.pending.is_some()).collect();
    if !pending.is_empty() {
        md += "\n## Pending cases\n\nReference values already cited; each case gets its recipe when the capability it needs lands.\n\n| Case | Kind | Problem | Needs | Reference values | Source |\n|---|---|---|---|---|---|\n";
        for (case, e) in pending {
            let values: Vec<String> = e
                .quantities
                .iter()
                .map(|q| format!("{} = {}", cell(&q.description), q.reference.map(number).unwrap_or_else(|| "to be transcribed".into())))
                .collect();
            let _ = writeln!(
                md,
                "| [{case}](../tests/oracle/{case}/) | {} | {} | {} | {} | {} |",
                e.kind,
                cell(&e.title),
                cell(e.pending.as_deref().unwrap_or_default()),
                values.join("; "),
                cell(&e.source)
            );
        }
    }
    md
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quantity(reference: f64, tolerance: Tolerance, known: Option<&str>) -> Quantity {
        Quantity {
            id: "q".into(),
            description: "q".into(),
            measure: json!({"from": "s", "pointer": "/x"}),
            reference: Some(reference),
            tolerance,
            reason: "r".into(),
            known_failure: known.map(str::to_string),
            fixture_relative: FIXTURE_RELATIVE,
            frozen: true,
        }
    }

    /// A quantity of a generated mesh is gated by its tolerance only: a stored fixture is ignored
    /// and blessing skips it.
    #[test]
    fn unfrozen_quantities_skip_the_fixture() {
        let t = Tolerance { relative: 0.03, absolute: 0.0 };
        let q = Quantity { frozen: false, ..quantity(213.0, t, None) };
        let old = Frozen { value: 213.3, band: 2e-7 };
        assert_eq!(judge(&q, &Ok(212.4), Some(&old)), (Verdict::Pass, FixtureCheck::NotKept));
        assert_eq!(judge(&q, &Ok(200.0), Some(&old)).0, Verdict::Fail, "the tolerance still gates it");
    }

    /// A nightly case not run here is judged by its fixture: frozen values pass as they did, a
    /// known failure stays known, and a fixture that misses an edited reference fails.
    #[test]
    fn nightly_cases_show_their_fixtures() {
        let t = Tolerance { relative: 0.01, absolute: 0.0 };
        let expected = Expected {
            title: "t".into(),
            kind: "benchmark".into(),
            source: "s".into(),
            pending: None,
            slow: false,
            nightly: true,
            quantities: vec![
                Quantity { id: "a".into(), ..quantity(2.0, t, None) },
                Quantity { id: "b".into(), ..quantity(2.0, t, Some("why")) },
                Quantity { id: "c".into(), ..quantity(3.0, t, None) },
            ],
        };
        let fixture: BTreeMap<String, Frozen> =
            [("a".to_string(), Frozen { value: 2.01, band: 1e-9 }), ("c".to_string(), Frozen { value: 2.0, band: 1e-9 })].into_iter().collect();
        let rows = fixture_rows("n", &expected, Some(&fixture));
        let got: Vec<(Option<f64>, Verdict, FixtureCheck)> = rows.into_iter().map(|r| (r.ours, r.verdict, r.fixture)).collect();
        assert_eq!(
            got,
            vec![
                (Some(2.01), Verdict::Pass, FixtureCheck::Nightly),
                (None, Verdict::Known, FixtureCheck::NotKept),
                (Some(2.0), Verdict::Fail, FixtureCheck::Nightly),
            ]
        );
        assert_eq!(nightly_problems("n", &fixture_rows("n", &expected, Some(&fixture))).len(), 1, "c's fixture misses its reference");
        // Never blessed: every frozen quantity is missing and fails, as for a case run here; the
        // known failure and an unfrozen quantity do not.
        let mut unblessed = expected.clone();
        unblessed.quantities.push(Quantity { id: "d".into(), frozen: false, ..quantity(1.0, t, None) });
        let none = fixture_rows("n", &unblessed, None);
        let checks: Vec<(Verdict, FixtureCheck)> = none.iter().map(|r| (r.verdict.clone(), r.fixture.clone())).collect();
        assert_eq!(
            checks,
            vec![
                (Verdict::Nightly, FixtureCheck::Missing),
                (Verdict::Known, FixtureCheck::NotKept),
                (Verdict::Nightly, FixtureCheck::Missing),
                (Verdict::Nightly, FixtureCheck::NotKept),
            ]
        );
        let problems = nightly_problems("n", &none);
        assert_eq!(problems.len(), 2, "{problems:?}");
        assert!(problems.iter().all(|p| p.contains("no fixture; bless")), "{problems:?}");
    }

    #[test]
    fn verdicts_and_fixtures() {
        let t = Tolerance { relative: 0.01, absolute: 0.0 };
        let q = quantity(2.0, t, None);
        let frozen = Frozen { value: 2.01, band: fixture_band(2.01, 2.0, t, FIXTURE_RELATIVE) };
        assert_eq!(judge(&q, &Ok(2.01), Some(&frozen)), (Verdict::Pass, FixtureCheck::Matches));
        assert_eq!(judge(&q, &Ok(2.011), Some(&frozen)), (Verdict::Pass, FixtureCheck::Moved { frozen: 2.01 }));
        assert_eq!(judge(&q, &Ok(2.01), None), (Verdict::Pass, FixtureCheck::Missing));
        assert_eq!(judge(&q, &Ok(2.1), None).0, Verdict::Fail);
        assert_eq!(judge(&quantity(2.0, t, Some("why")), &Ok(2.1), None), (Verdict::Known, FixtureCheck::NotKept));
        assert_eq!(judge(&quantity(2.0, t, Some("why")), &Ok(2.0), None).0, Verdict::Fixed);
        assert!(matches!(judge(&q, &Err("boom".into()), None).0, Verdict::Error(_)));
        // Zero references use the absolute tolerance, and the band keeps a floor from it.
        let a = Tolerance { relative: 0.0, absolute: 1e-9 };
        assert!(a.allows(0.0, 5e-10) && !a.allows(0.0, 2e-9));
        assert_eq!(fixture_band(1e-16, 0.0, a, FIXTURE_RELATIVE), 1e-9, "a zero reference keeps its tolerance as the band");
    }

    #[test]
    fn measures_read_pointers_scale_and_divide() {
        let mut session = Session::default();
        let outputs: BTreeMap<String, Value> = [("s".to_string(), json!({"a": 6.0, "b": -2.0, "z": 0.0}))].into_iter().collect();
        let m = |v: Value| measure(&mut Session::default(), &outputs, &v);
        assert_eq!(m(json!({"from": "s", "pointer": "/a", "scale": 0.5})), Ok(3.0));
        assert_eq!(m(json!({"from": "s", "pointer": "/a", "over": {"from": "s", "pointer": "/b", "scale": -1.0}})), Ok(3.0));
        assert!(m(json!({"from": "s", "pointer": "/a", "over": {"from": "s", "pointer": "/z"}})).is_err());
        assert!(m(json!({"from": "s", "pointer": "/a", "over": {"from": "s", "pointer": "/nope"}})).is_err());
        assert!(measure(&mut session, &outputs, &json!({"from": "t", "pointer": "/a"})).is_err());
    }

    #[test]
    fn expected_files_are_checked() {
        let good = json!({"title": "t", "kind": "theory", "source": "s", "quantities": [
            {"id": "a", "quantity": "a", "measure": {"from": "x", "pointer": "/a"}, "value": 1.0, "tolerance": {"relative": 0.1}, "reason": "r"}]});
        assert!(Expected::parse(&good).is_ok());
        let pending = json!({"title": "t", "kind": "benchmark", "source": "s", "status": "pending", "requires": "shells",
            "quantities": [{"id": "a", "quantity": "a", "value": null}]});
        assert!(Expected::parse(&pending).unwrap().pending.is_some());
        for (key, bad) in [("kind", json!("guess")), ("quantities", json!([])), ("status", json!("maybe")), ("title", json!(""))] {
            let mut v = good.clone();
            v[key] = bad;
            assert!(Expected::parse(&v).is_err(), "{key}");
        }
        for (key, bad) in
            [("tolerance", json!({})), ("tolerance", json!({"relative": -1})), ("value", json!("x")), ("reason", json!("")), ("measure", json!({}))]
        {
            let mut v = good.clone();
            v["quantities"][0][key] = bad;
            assert!(Expected::parse(&v).is_err(), "{key}");
        }
    }

    #[test]
    fn numbers_print_with_five_digits() {
        assert_eq!(number(92.7), "92.7");
        assert_eq!(number(-4.0241234), "-4.0241");
        assert_eq!(number(0.30241), "0.30241");
        assert_eq!(number(1.8248e-5), "1.8248e-5");
        assert_eq!(number(210000.0), "2.1000e5");
        assert_eq!(number(0.0), "0");
        assert_eq!(percent(-0.0050317), "-0.503%");
        assert_eq!(percent(0.000132), "0.0132%");
        assert_eq!(percent(0.02), "2%");
        assert_eq!(fraction(1e-9), "1e-9 rel.");
    }

    #[test]
    fn numbers_take_the_notation_of_their_rounded_value() {
        assert_eq!(number(0.000999999999), "0.001");
        assert_eq!(number(0.001), "0.001");
        assert_eq!(number(-0.000999999999), "-0.001");
        assert_eq!(number(99999.99999), "1.0000e5");
        assert_eq!(number(1.0e5), "1.0000e5");
        assert_eq!(number(3.48751e-4), "3.4875e-4");
        assert_eq!(number(1.234567e9), "1.2346e9");
        assert_eq!(number(-2.5e-12), "-2.5000e-12");
    }

    fn row(reference: f64, tolerance: Tolerance, ours: f64) -> Row {
        Row {
            case: "c".into(),
            quantity: quantity(reference, tolerance, None),
            ours: Some(ours),
            verdict: Verdict::Pass,
            fixture: FixtureCheck::Matches,
        }
    }

    /// Errors at the scale of the fixture band are noise: they print as the page's bound.
    #[test]
    fn noise_errors_print_as_a_bound() {
        let rel = Tolerance { relative: 0.001, absolute: 0.0 };
        assert_eq!(error_of(&row(12572.0, rel, 12572.0 * (1.0 + 1.78e-8))), "< 1e-4%");
        assert_eq!(error_of(&row(12572.0, rel, 12572.0 * (1.0 + 1.74e-8))), "< 1e-4%");
        assert_eq!(error_of(&row(1.0, rel, 1.0)), "< 1e-4%");
        let abs = Tolerance { relative: 0.0, absolute: 1e-9 };
        assert_eq!(error_of(&row(-0.011149, abs, -0.011149 - 1.2e-16)), "< 1e-7");
        assert_eq!(error_of(&row(-0.011149, abs, -0.011149 - 1.4e-16)), "< 1e-7");
    }

    #[test]
    fn real_errors_still_print() {
        let rel = Tolerance { relative: 0.02, absolute: 0.0 };
        assert_eq!(error_of(&row(1.0, rel, 0.99183)), "-0.817%");
        assert_eq!(error_of(&row(0.024355, rel, 0.024357)), "0.00821%");
        let abs = Tolerance { relative: 0.0, absolute: 0.1 };
        assert_eq!(error_of(&row(0.21104, abs, 0.21895)), "7.9e-3");
    }

    /// The page must not change when a result moves by a few ulps (summation order, platform maths).
    #[test]
    fn the_page_ignores_a_few_ulps() {
        let rel = Tolerance { relative: 0.01, absolute: 0.0 };
        let abs = Tolerance { relative: 0.0, absolute: 1e-9 };
        let cases = [
            (0.001, rel, 0.000999999999),
            (12572.0, rel, 12572.000223),
            (-0.011149, abs, -0.011149 - 1.2e-16),
            (0.0, abs, 3e-17),
            (1.0, rel, 0.99183),
            (3.4875e-4, rel, 3.4875e-4 * (1.0 - 6.65e-9)),
            (210000.0, rel, 209999.99999),
        ];
        for (reference, tolerance, ours) in cases {
            let base = row(reference, tolerance, ours);
            let page = (number(reference), ours_text(&base), error_of(&base));
            for ulps in [-4i64, -1, 1, 4] {
                let moved = f64::from_bits((ours.to_bits() as i64 + ulps) as u64);
                let r = row(reference, tolerance, moved);
                assert_eq!((number(reference), ours_text(&r), error_of(&r)), page, "{ours} moved by {ulps} ulps");
            }
        }
    }
}
