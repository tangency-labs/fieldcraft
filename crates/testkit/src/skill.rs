//! The meshing skill's test (#482), behind `cargo xtask skill-results`: every corpus part and
//! every sample with a design, meshed by following `skills/fieldcraft-meshing/SKILL.md` step by
//! step with today's commands, through [`Session::execute`] as an agent over MCP would. The
//! results go to docs/meshing/skill-results.md, one row per part.
//!
//! The steps, per part: look (`geometry.bodies`, `geometry.features`, `geometry.check` at the
//! planned size); plan (`mesh.plan`) and apply it (`mesh.plan.apply`, which uses the body's mirror
//! symmetry itself); where the plan fails, the skill's fallbacks in order (`mesh.hex.generate`,
//! then `mesh.tet.generate` tet10), each failure kept as the reason; `mesh.quality`;
//! `mesh.quality.mark` at the skill's limits; `mesh.node.relax` round each marked element; mark
//! again; `mesh.check.duplicates`. The commands that built the final mesh are written as a
//! script per part, so `fieldcraft-cli snapshot --section` can draw its section picture.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use fieldcraft_engine::Session;
use serde_json::{Value, json};

use crate::corpus;

pub const FILE: &str = "docs/meshing/skill-results.md";
/// The section pictures, one per part, relative to the workspace root.
pub const PICTURES: &str = "docs/meshing/skill";
/// The skill's acceptance limits (§6).
const MIN_SCALED_JACOBIAN: f64 = 0.3;
const MAX_EDGE_RATIO: f64 = 10.0;
/// A quality measure that changes by less than this is unchanged.
const ROUND_OFF: f64 = 1e-9; // unit-factors: a change of a quality measure
/// The mesh's volume against the geometry's (§6.1): a larger difference is a missing cell or a
/// mesh too coarse for the curvature.
const MAX_VOLUME_ERROR: f64 = 0.01;
/// The most marked elements relaxed one by one; more than this is a clustered failure (the
/// blocks are wrong, §6.3), which relaxing does not fix.
const MAX_RELAXED: usize = 20;

/// The Raspberry Pi 5 enclosure analysis design (base, lid, heat sink, board, chips, air).
const PI5: &str = include_str!("../../engine/tests/designs/pi5-enclosure.json");
/// The cup of the projection mesher's tests (#418): a Ø40 × 30 cylinder with a 4 mm floor and
/// wall, and a box handle joined to its side.
const CUP: &str = r#"[
    {"command": "solid.cylinder", "params": {"radius": 20, "height": 30, "body_name": "Cup"}},
    {"command": "solid.cylinder", "params": {"base": [0, 0, 4], "radius": 16, "height": 27, "operation": "cut", "targets": ["Cup"]}},
    {"command": "solid.box", "params": {"length": 11, "width": 6, "height": 16, "corner": [19, -3, 8], "operation": "new", "body_name": "Handle"}},
    {"command": "solid.box", "params": {"length": 8, "width": 8, "height": 8, "corner": [18, -4, 12], "operation": "cut", "targets": ["Handle"]}},
    {"command": "solid.combine", "params": {"target": "Cup", "tools": ["Handle"], "operation": "join"}}
]"#;
const PLAIN_CUP: &str = r#"[
    {"command": "solid.cylinder", "params": {"radius": 20, "height": 30, "body_name": "Cup"}},
    {"command": "solid.cylinder", "params": {"base": [0, 0, 4], "radius": 16, "height": 27, "operation": "cut", "targets": ["Cup"]}}
]"#;

/// Why a part fails the skill's limits today, read from its numbers and its section picture.
/// Shown only while the part fails.
const CAUSES: &[(&str, &str)] = &[
    (
        "loft-transition",
        "the lofted sides are swept one to one, and the elements where the section changes shape are skewed (min SJ 0.19). Relaxing the four marked elements did not raise it: the blocks are wrong, not the nodes (§6.3). A loft has no template yet (`unmapped`).",
    ),
    (
        "sharp-wedge",
        "the wedge's sharp edge is swept as hex, so every element tapers to it (min SJ 0.01, edge ratio 150). A column of wedge6 at the edge (§4) raises the SJ but not the ratio: any element with a corner on a 2.9° edge has sides at least 1/tan 2.9° ≈ 20 times its width there, twice the limit. The edge has to be blunted (defeatured) or the limit waived along it.",
    ),
    (
        "analytic-sphere",
        "no structured method takes a sphere (it has no face to sweep from or axis to revolve about), so the plan does not apply and the fallback's sphere template meshes it at the plan's size (2.45 mm, which keeps a sphere's chords to 0.75 % on paper), 1.1 % short of the volume: the template's elements at the cube's eight corners span more of the sphere than the size. A sphere template in the plan (§3) would size its caps; hex20 (`order: 2`) closes it too.",
    ),
    (
        "cross-cylinders",
        "the boundary tessellation of the joined cylinders is invalid (inconsistent orientation at the two singular points), so every mesher refuses the body.",
    ),
];

/// A part the skill is run on.
pub struct Case {
    pub name: String,
    /// Where it comes from and what it is.
    pub about: String,
    /// Its `.solidcraft` bytes.
    pub design: Vec<u8>,
    /// The bodies meshed (none: every body of the design).
    pub bodies: Vec<String>,
    /// The bodies are cells of one part (split along shared faces), joined by merging nodes.
    pub cells: bool,
}

fn script(commands: &str) -> Result<Vec<u8>, String> {
    let wrapped = if commands.trim_start().starts_with('[') { format!("{{\"commands\": {commands}}}") } else { commands.to_string() };
    fieldcraft_geometry::script_design(&wrapped).map_err(|e| e.to_string())
}

/// The corpus parts, then the samples' designs that are not corpus parts.
pub fn cases() -> Result<Vec<Case>, String> {
    let mut out = Vec::new();
    for part in corpus::parts().map_err(|e| e.to_string())? {
        let design = corpus::design_bytes(&part).map_err(|e| e.to_string())?;
        let cells = part.categories.iter().any(|c| c == "partitioned");
        out.push(Case { name: part.name.clone(), about: format!("corpus: {}", part.description), design, bodies: Vec::new(), cells });
    }
    use fieldcraft_engine::sample;
    let samples: [(&str, &str, &str, &[&str]); 6] = [
        ("crank-arm", "sample `crank-arm`: the 170 mm bicycle crank, mirror symmetric across y = 0", sample::CRANK_ARM, &[]),
        ("upright-bracket", "samples `drop-test`, `flow-past`: the upright L-bracket", sample::UPRIGHT_BRACKET, &[]),
        ("pi5-base", "the Raspberry Pi 5 enclosure's base: a shelled box with fillets and standoffs", PI5, &["Base"]),
        ("pi5-lid", "the Raspberry Pi 5 enclosure's lid", PI5, &["Lid"]),
        ("cup", "the cup with a handle of #418", CUP, &[]),
        ("plain-cup", "the cup without its handle (#418)", PLAIN_CUP, &[]),
    ];
    for (name, about, commands, bodies) in samples {
        out.push(Case {
            name: name.into(),
            about: about.into(),
            design: script(commands)?,
            bodies: bodies.iter().map(|b| b.to_string()).collect(),
            cells: false,
        });
    }
    Ok(out)
}

/// One part's result.
#[derive(Default)]
pub struct Outcome {
    pub method: String,
    /// Why each earlier method failed, in order.
    pub failures: Vec<String>,
    pub symmetry: String,
    pub geometry_findings: String,
    pub elements: usize,
    pub kinds: String,
    pub min_sj: Option<f64>,
    pub mean_sj: Option<f64>,
    pub max_edge_ratio: Option<f64>,
    pub volume_error: Option<f64>,
    pub untagged: usize,
    /// Marked before relaxing, then after.
    pub marked: (usize, usize),
    /// The minimum scaled Jacobian and the largest edge ratio before relaxing.
    pub before: (Option<f64>, Option<f64>),
    pub relaxed: usize,
    /// Relaxations undone because they left their region worse.
    pub undone: usize,
    pub cracks: String,
    /// Coincident nodes merged between cells.
    pub merged: usize,
    pub seconds: f64,
    /// The commands that built the final mesh (for the section picture).
    pub replay: Vec<Value>,
    /// The section: axis (0..3) and its coordinate.
    pub section: (usize, f64),
}

/// Run a command, keeping it in `replay` when it succeeds.
fn exec(s: &mut Session, replay: &mut Vec<Value>, id: &str, params: Value) -> Result<Value, String> {
    let r = s.execute(id, &params).map_err(|e| e.to_string())?;
    replay.push(json!({"command": id, "params": params}));
    Ok(r)
}

fn short(e: &str) -> String {
    let one = e.replace(['\n', '|'], " ");
    if one.chars().count() > 240 { format!("{}…", one.chars().take(240).collect::<String>()) } else { one }
}

/// Follow the skill on one part, its design written to `path`.
pub fn run_case(case: &Case, path: &Path) -> Result<Outcome, String> {
    std::fs::write(path, &case.design).map_err(|e| format!("{}: {e}", path.display()))?;
    let start = Instant::now();
    let mut s = Session::default();
    let mut o = Outcome::default();
    let mut replay = Vec::new();
    exec(&mut s, &mut replay, "geometry.open", json!({"path": path.to_string_lossy()}))?;

    // §1 Look.
    let bodies = s.execute("geometry.bodies", &json!({})).map_err(|e| e.to_string())?;
    let all: Vec<&Value> = bodies.as_array().map(|a| a.iter().collect()).unwrap_or_default();
    let targets: Vec<&Value> =
        all.into_iter().filter(|b| case.bodies.is_empty() || case.bodies.iter().any(|n| b["name"].as_str() == Some(n))).collect();
    if targets.is_empty() {
        return Err("no body to mesh".into());
    }
    let names: Vec<String> = targets.iter().filter_map(|b| b["name"].as_str().map(String::from)).collect();
    let volume: f64 = targets.iter().filter_map(|b| b["volume"].as_f64()).sum();
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for b in &targets {
        for k in 0..3 {
            lo[k] = lo[k].min(b["bbox"][0][k].as_f64().unwrap_or(f64::INFINITY));
            hi[k] = hi[k].max(b["bbox"][1][k].as_f64().unwrap_or(f64::NEG_INFINITY));
        }
    }
    // The section across the longest side, through the middle: it cuts the primary features.
    let axis = (0..3).max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b]))).unwrap_or(0);
    o.section = (axis, 0.5 * (lo[axis] + hi[axis]));

    let mut findings = Vec::new();
    let mut methods = Vec::new();
    let mut symmetry = Vec::new();
    for body in &names {
        let _ = s.execute("geometry.features", &json!({"body": body}));
        // §3 Plan: the size it picks is the size the geometry is checked at (§1.2).
        let plan = exec(&mut s, &mut replay, "mesh.plan", json!({"body": body}));
        let size = plan.as_ref().ok().and_then(|p| p["size"].as_f64());
        if let Ok(c) = s.execute("geometry.check", &json!({"body": body, "size": size})) {
            let counts: Vec<String> = c["counts"]
                .as_object()
                .map(|m| m.iter().filter(|(_, n)| n.as_u64().unwrap_or(0) > 0).map(|(k, n)| format!("{n} {k}")).collect())
                .unwrap_or_default();
            if !counts.is_empty() {
                findings.push(if names.len() > 1 { format!("{body}: {}", counts.join(", ")) } else { counts.join(", ") });
            }
        }
        let mut method = None;
        match plan {
            Ok(p) => match exec(
                &mut s,
                &mut replay,
                "mesh.plan.apply",
                json!({"plan": p["name"].as_str().map_or_else(|| format!("Plan-{body}"), String::from)}),
            ) {
                Ok(r) => {
                    let skipped = r["skipped_cuts"].as_array().map_or(0, Vec::len);
                    let mesher = r["mesher"].as_str().or_else(|| p["mesher"].as_str()).unwrap_or("?");
                    method = Some(if skipped > 0 { format!("plan ({mesher}, {skipped} cuts skipped)") } else { format!("plan ({mesher})") });
                    if let Some(planes) = r["symmetry"].as_array().filter(|a| !a.is_empty()) {
                        symmetry.push(format!("{} plane(s)", planes.len()));
                    }
                }
                Err(e) => o.failures.push(format!("{body}: mesh.plan.apply: {}", short(&e))),
            },
            Err(e) => o.failures.push(format!("{body}: mesh.plan: {}", short(&e))),
        }
        // §3.4 The plan does not fit: the structured meshers, then tet10 for what is not worth
        // blocking.
        if method.is_none() {
            match exec(&mut s, &mut replay, "mesh.hex.generate", json!({"body": body, "size": size})) {
                Ok(r) => method = Some(format!("mesh.hex.generate ({})", r["method"].as_str().unwrap_or("auto"))),
                Err(e) => o.failures.push(format!("{body}: mesh.hex.generate: {}", short(&e))),
            }
        }
        if method.is_none() {
            match exec(&mut s, &mut replay, "mesh.hex.project", json!({"body": body, "size": size, "template": "sphere"})) {
                Ok(_) => method = Some("mesh.hex.project (sphere)".into()),
                Err(e) => o.failures.push(format!("{body}: mesh.hex.project (sphere): {}", short(&e))),
            }
        }
        if method.is_none() {
            match exec(&mut s, &mut replay, "mesh.tet.generate", json!({"body": body, "size": size, "order": "quadratic"})) {
                Ok(_) => method = Some("tet10".into()),
                Err(e) => o.failures.push(format!("{body}: mesh.tet.generate: {}", short(&e))),
            }
        }
        methods.push(method.unwrap_or_else(|| "none".into()));
    }
    methods.dedup();
    o.method = methods.join(", ");
    o.symmetry = if symmetry.is_empty() { "–".into() } else { symmetry.join(", ") };
    o.geometry_findings = if findings.is_empty() { "–".into() } else { findings.join("; ") };
    if s.model.elements().is_empty() {
        o.seconds = start.elapsed().as_secs_f64();
        return Ok(o);
    }

    // §6.5 Cells of one part meshed one by one join by merging their coincident nodes (separate
    // parts are tied instead, never merged).
    if case.cells {
        let r = exec(&mut s, &mut replay, "mesh.node.merge", json!({}))?;
        o.merged = r["merged"].as_u64().unwrap_or(0) as usize;
    }

    // §6 Quality, mark, relax.
    let before = s.execute("mesh.quality", &json!({})).map_err(|e| e.to_string())?;
    o.before = (before["metrics"]["scaled_jacobian"]["min"].as_f64(), worst_ratio(&before));
    let mark_params = json!({"limits": {"scaled_jacobian": MIN_SCALED_JACOBIAN, "edge_ratio": MAX_EDGE_RATIO}});
    let marked = exec(&mut s, &mut replay, "mesh.quality.mark", mark_params.clone())?;
    let labels: Vec<u64> = marked["elements"].as_array().map(|a| a.iter().filter_map(|e| e["label"].as_u64()).collect()).unwrap_or_default();
    o.marked.0 = marked["marked"].as_u64().unwrap_or(0) as usize;
    if o.marked.0 > 0 && o.marked.0 <= MAX_RELAXED {
        for label in labels.iter().take(MAX_RELAXED) {
            let params = json!({"around_element": label, "reapply": true});
            let Ok(r) = s.execute("mesh.node.relax", &params) else { continue };
            // §6.3: a relaxation that leaves its region worse is undone.
            let worse = |k: &str, higher_is_better: bool| match (r["before"][k].as_f64(), r["after"][k].as_f64()) {
                (Some(b), Some(a)) => {
                    if higher_is_better {
                        a < b - ROUND_OFF
                    } else {
                        a > b.max(MAX_EDGE_RATIO) + ROUND_OFF
                    }
                }
                _ => false,
            };
            if worse("min_scaled_jacobian", true) || worse("max_edge_ratio", false) {
                s.execute("edit.undo", &json!({})).map_err(|e| e.to_string())?;
                o.undone += 1;
            } else {
                replay.push(json!({"command": "mesh.node.relax", "params": params}));
                o.relaxed += 1;
            }
        }
    }
    let marked = if o.relaxed + o.undone > 0 { exec(&mut s, &mut replay, "mesh.quality.mark", mark_params)? } else { marked };
    o.marked.1 = marked["marked"].as_u64().unwrap_or(0) as usize;

    let q = s.execute("mesh.quality", &json!({})).map_err(|e| e.to_string())?;
    let metric = |m: &str, k: &str| q["metrics"][m][k].as_f64();
    o.min_sj = metric("scaled_jacobian", "min");
    o.mean_sj = metric("scaled_jacobian", "mean");
    o.max_edge_ratio = worst_ratio(&q);
    o.volume_error = q["volume"].as_f64().filter(|_| volume > 0.0).map(|v| (v - volume) / volume);
    o.elements = s.model.elements().len();
    let mut kinds = std::collections::BTreeMap::<String, usize>::new();
    for e in s.model.elements().values() {
        *kinds.entry(format!("{:?}", e.kind).to_lowercase()).or_default() += 1;
    }
    o.kinds = kinds.iter().map(|(k, n)| format!("{n} {k}")).collect::<Vec<_>>().join(", ");
    o.untagged = fieldcraft_engine::cmd::untagged_faces(&s, &names).len();
    o.cracks = match s.execute("mesh.check.duplicates", &json!({})) {
        Ok(d) => {
            let n = |k: &str| d[k].as_array().map_or(0, Vec::len);
            let (coincident, cracks) = (n("coincident"), n("cracks"));
            if coincident + cracks == 0 { "none".into() } else { format!("{coincident} coincident, {cracks} cracks") }
        }
        Err(e) => short(&e.to_string()),
    };
    o.seconds = start.elapsed().as_secs_f64();
    o.replay = replay;
    Ok(o)
}

/// The largest edge ratio of a quality reply (hex, wedge), else the largest aspect ratio (tet,
/// which has no edge ratio).
fn worst_ratio(q: &Value) -> Option<f64> {
    q["metrics"]["edge_ratio"]["max"].as_f64().or_else(|| q["metrics"]["aspect_ratio"]["max"].as_f64())
}

fn num(v: Option<f64>, digits: usize) -> String {
    v.map_or("–".into(), |x| format!("{x:.digits$}"))
}

/// Whether the part meets the skill's limits.
fn verdict(o: &Outcome) -> String {
    let (Some(sj), Some(ratio)) = (o.min_sj, o.max_edge_ratio) else { return "no mesh".into() };
    let mut failed = Vec::new();
    if sj < MIN_SCALED_JACOBIAN {
        failed.push("SJ");
    }
    if ratio > MAX_EDGE_RATIO {
        failed.push("ratio");
    }
    if o.volume_error.is_some_and(|v| v.abs() > MAX_VOLUME_ERROR) {
        failed.push("volume");
    }
    match (failed.is_empty(), o.method.starts_with("plan")) {
        (false, _) => format!("fail ({})", failed.join(", ")),
        (true, true) => "pass".into(),
        (true, false) => "pass (fallback)".into(),
    }
}

/// Run the skill on every case (or the one named `only`), write the replay scripts and the
/// picture list to `work`, and the results to [`FILE`] unless `only` is given or `write` is
/// false. The picture list (`pictures.json`) is what `cargo xtask skill-results` draws with
/// `fieldcraft-cli snapshot`.
pub fn run(root: &Path, work: &Path, only: Option<&str>, write: bool) -> Result<(), String> {
    std::fs::create_dir_all(work).map_err(|e| format!("{}: {e}", work.display()))?;
    let work = &std::fs::canonicalize(work).map_err(|e| format!("{}: {e}", work.display()))?;
    let mut cases = cases()?;
    if let Some(name) = only {
        cases.retain(|c| c.name == name);
        if cases.is_empty() {
            return Err(format!("no part named `{name}`"));
        }
    }
    let mut rows = Vec::new();
    let mut notes = Vec::new();
    let mut pictures = Vec::new();
    for case in &cases {
        eprintln!("skill: {}", case.name);
        let design: PathBuf = work.join(format!("{}.solidcraft", case.name));
        let o = match run_case(case, &design) {
            Ok(o) => o,
            Err(e) => Outcome { method: "none".into(), failures: vec![short(&e)], ..Outcome::default() },
        };
        for f in &o.failures {
            eprintln!("skill: {}: {f}", case.name);
        }
        let verdict = verdict(&o);
        let picture = if o.replay.is_empty() {
            "–".to_string()
        } else {
            let replay = work.join(format!("{}.json", case.name));
            std::fs::write(&replay, Value::Array(o.replay.clone()).to_string()).map_err(|e| format!("{}: {e}", replay.display()))?;
            let png = format!("{PICTURES}/{}.png", case.name);
            let section = format!("{}={:.3}", ["x", "y", "z"][o.section.0], o.section.1);
            pictures.push(json!({"script": replay.to_string_lossy(), "out": png, "section": section}));
            format!("[{section}](skill/{}.png)", case.name)
        };
        rows.push(vec![
            case.name.clone(),
            o.method.clone(),
            o.symmetry.clone(),
            if o.elements > 0 { format!("{} ({})", o.elements, o.kinds) } else { "–".into() },
            num(o.min_sj, 3),
            num(o.mean_sj, 3),
            num(o.max_edge_ratio, 1),
            o.volume_error.map_or("–".into(), |v| format!("{v:.1e}")),
            if o.elements > 0 { o.untagged.to_string() } else { "–".into() },
            if o.elements == 0 {
                "–".into()
            } else if o.relaxed == 0 {
                let why = if o.undone > 0 { format!("{} relaxations undone: worse", o.undone) } else { "not relaxed".into() };
                format!("{} ({why})", o.marked.0)
            } else {
                format!(
                    "{} → {} ({} relaxed, {} undone; min SJ {} → {}, max ratio {} → {})",
                    o.marked.0,
                    o.marked.1,
                    o.relaxed,
                    o.undone,
                    num(o.before.0, 3),
                    num(o.min_sj, 3),
                    num(o.before.1, 1),
                    num(o.max_edge_ratio, 1)
                )
            },
            match (o.cracks.is_empty(), o.merged) {
                (true, _) => "–".into(),
                (false, 0) => o.cracks.clone(),
                (false, m) => format!("{} (after merging {m} nodes)", o.cracks),
            },
            o.geometry_findings.clone(),
            format!("{:.1}", o.seconds),
            verdict.clone(),
            picture,
        ]);
        let cause = CAUSES.iter().find(|(name, _)| *name == case.name).map(|(_, c)| *c).filter(|_| !verdict.starts_with("pass"));
        if !o.failures.is_empty() || cause.is_some() {
            let mut note = format!("- `{}` ({}):", case.name, case.about);
            for f in &o.failures {
                let _ = write!(note, "\n  - {f}");
            }
            if let Some(c) = cause {
                let _ = write!(note, "\n  - **{verdict}**: {c}");
            }
            notes.push(note);
        }
    }
    let header = [
        "part",
        "method",
        "symmetry",
        "elements",
        "min SJ",
        "mean SJ",
        "max edge ratio¹",
        "volume error",
        "untagged faces",
        "marked → after relax",
        "duplicates",
        "geometry.check",
        "time s",
        "verdict",
        "section",
    ];
    let mut table = format!("| {} |\n|{}\n", header.join(" | "), "---|".repeat(header.len()));
    for r in &rows {
        let _ = writeln!(table, "| {} |", r.join(" | "));
    }
    println!("{table}");
    let list = work.join("pictures.json");
    std::fs::write(&list, Value::Array(pictures).to_string()).map_err(|e| format!("{}: {e}", list.display()))?;
    if write && only.is_none() {
        let path = root.join(FILE);
        std::fs::write(&path, document(&table, &notes, &cases)).map_err(|e| format!("{}: {e}", path.display()))?;
        println!("skill: wrote {FILE}");
    }
    Ok(())
}

fn document(table: &str, notes: &[String], cases: &[Case]) -> String {
    let mut d = String::from(
        "# The meshing skill on the corpus and the samples\n\n\
Generated by `cargo xtask skill-results`; do not edit. Each part is meshed by following\n\
[the meshing skill](../../skills/fieldcraft-meshing/SKILL.md) with today's commands, headless, as\n\
an agent over MCP would (`crates/testkit/src/skill.rs`): `geometry.bodies`, `geometry.features`,\n\
`mesh.plan` and `geometry.check` at the plan's size, `mesh.plan.apply` (which meshes one cell of a\n\
mirror-symmetric body and mirrors it), and where the plan fails the skill's fallbacks in order,\n\
`mesh.hex.generate`, `mesh.hex.project` (sphere template) and then `mesh.tet.generate` (tet10);\n\
the cells of a split part are joined with `mesh.node.merge`; then `mesh.quality.mark` at the\n\
skill's limits (scaled Jacobian ≥ 0.3, edge ratio ≤ 10), `mesh.node.relax` round each marked\n\
element when at most 20 are marked (more is a clustered failure: the blocks are wrong, and relaxing\n\
is not the fix), `mesh.quality.mark` again, `mesh.quality` and `mesh.check.duplicates`. Sizes are\n\
the plan's (1/25 of the body's diagonal, 2 layers through a shell).\n\n\
Columns: `method` is what meshed the part (`plan (project|cells|sweep)` is the skill's target,\n\
anything else a fallback); `symmetry` the mirror planes `mesh.plan.apply` used; `volume error` the\n\
mesh's volume against the geometry's; `untagged faces` B-rep faces no element face is tagged\n\
with (a set on them finds nothing); `marked` the elements failing the limits before and after\n\
relaxing, and what relaxing changed; ¹ a tet mesh has no edge ratio: its largest aspect ratio\n\
stands in; `duplicates` coincident nodes and cracks; `geometry.check` its findings at the size;\n\
`verdict` pass when the minimum scaled Jacobian and the largest edge ratio meet the limits and\n\
the volume is within 1 % (`pass (fallback)` when a fallback meshed it).\n\
The section pictures cut each part through the middle of its longest side, marked elements red.\n\
Times are wall clock on the machine that ran it.\n\n## Results\n\n",
    );
    d.push_str(table);
    d.push_str("\n## Why a part fails\n\nEach earlier method's error, in the order the skill tries them, and for a part that fails the\nlimits the cause read from its numbers and section picture.\n\n");
    d.push_str(&notes.join("\n"));
    d.push_str("\n\n## Parts\n\n");
    for c in cases {
        let _ = writeln!(d, "- `{}`: {}", c.name, c.about);
    }
    d
}
