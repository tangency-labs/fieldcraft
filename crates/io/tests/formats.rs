//! Files of every released format open, migrate to what their own writer's reader gave, and
//! save again at the current format (crates/io/tests/formats/v<N>/README.md).

use std::io::Cursor;
use std::path::{Path, PathBuf};

use fieldcraft_io::migrate::{self, Parts};
use fieldcraft_io::{FORMAT, FORMAT_VERSION, ModelFile, format_of, read_versioned, write};
use serde_json::{Value, json};

/// A fixture: the program that wrote it (its folder), the file, and what that program's own
/// reader gave for it (`<name>.read.json`).
struct Fixture {
    writer: String,
    path: PathBuf,
    recorded: Value,
}

impl Fixture {
    fn name(&self) -> String {
        format!("{}/{}", self.writer, self.path.file_name().unwrap().to_string_lossy())
    }
}

/// Every fixture of format `n`, from every writer's folder.
fn fixtures(n: u32) -> Vec<Fixture> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats").join(format!("v{n}"));
    let mut out = Vec::new();
    for writer in std::fs::read_dir(&root).unwrap().map(|e| e.unwrap().path()).filter(|p| p.is_dir()) {
        for path in std::fs::read_dir(&writer).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|x| x == "fieldcraft")) {
            let recorded = serde_json::from_slice(&std::fs::read(path.with_extension("read.json")).unwrap()).unwrap();
            out.push(Fixture { writer: writer.file_name().unwrap().to_string_lossy().into_owned(), path, recorded });
        }
    }
    out.sort_by(|a, b| a.path.cmp(&b.path));
    out
}

/// What the current reader must give for a file of format `from` whose writer's reader gave
/// `recorded`: the steps from `from` to the current format, applied to the recorded JSON.
/// Format 1 → 2 changes nothing. A record from a reader that had no view state (0.4.0, 0.5.0)
/// holds `null`: the model has no hidden items. Format 2 → 3 stores SI ([`to_si`]) and heads the
/// journal with the system its numbers were typed in.
fn expected(mut recorded: Value, from: u32) -> Value {
    assert!((1..=FORMAT_VERSION).contains(&from));
    if recorded["view"].is_null() {
        recorded["view"] = json!({});
    }
    if from <= 2 {
        let declared = recorded["model"]["units"].as_str().unwrap_or("unspecified").to_string();
        let system = if declared == "unspecified" { "mm-t-s-N-MPa".to_string() } else { declared.clone() };
        let model = recorded["model"].as_object_mut().unwrap();
        model.remove("units");
        model.insert("display_units".into(), json!(system));
        let in_mm = millimetre_nodes(&recorded, &declared);
        let built = built_on_geometry(&recorded["geometry"]);
        to_si(&mut recorded["model"], &declared, &in_mm, built);
        fluid_zero(&mut recorded["model"], &declared);
        let commands = recorded["journal"].take();
        recorded["journal"] = json!({"units": system, "model_units": "mm-t-s-N-MPa", "commands": commands});
    }
    recorded
}

/// The nodes (labels) of the mesh bodies a model of metres or inches made from its millimetre
/// design or 3MF meshes, as its writer recorded them: every body's but those `mesh.block` made
/// (in the journal), when the file names a design or meshes.
fn millimetre_nodes(recorded: &Value, declared: &str) -> std::collections::BTreeSet<String> {
    let mut out = std::collections::BTreeSet::new();
    if !matches!(declared, "m-kg-s-N-Pa" | "in-lbf-s-psi") || !built_on_geometry(&recorded["geometry"]) {
        return out;
    }
    let blocks: Vec<&str> = recorded["journal"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(|c| c["command"] == "mesh.block")
        .filter_map(|c| c["params"]["body"].as_str())
        .collect();
    let model = &recorded["model"];
    let entities = model["classification"]["entities"].as_array().cloned().unwrap_or_default();
    for (label, tags) in model["classification"]["elements"].as_object().into_iter().flatten() {
        let Some(body) = tags["region"].as_u64().and_then(|r| entities.get(r as usize)).and_then(|e| e["body"].as_str()) else { continue };
        if blocks.contains(&body) {
            continue;
        }
        for n in model["elements"][label]["nodes"].as_array().into_iter().flatten() {
            out.insert(n.to_string());
        }
    }
    out
}

/// Whether a file was built on a design or 3MF meshes (its geometry as recorded).
fn built_on_geometry(geometry: &Value) -> bool {
    !geometry["design"].is_null() || geometry["meshes"].as_array().is_some_and(|m| !m.is_empty())
}

/// A fluid's heat keeps counting from the zero of the system it was typed in (format 3's
/// `temperature_zero`, written where it is not 0 °C).
fn fluid_zero(model: &mut Value, declared: &str) {
    let zero = match declared {
        "m-kg-s-N-Pa" => 0.0,
        "in-lbf-s-psi" => 459.67 * (5.0 / 9.0),
        _ => return,
    };
    if let Some(domains) = model.get_mut("fluid_domains").and_then(Value::as_object_mut) {
        for d in domains.values_mut() {
            if let Some(t) = d.get_mut("thermal").and_then(Value::as_object_mut) {
                t.insert("temperature_zero".into(), json!(zero));
            }
        }
    }
}

/// The base units of a format 2 system, as (length, mass, time) powers of ten of SI's (the
/// metric ones), and its absolute temperature unit.
fn base(system: &str) -> ([i32; 3], &'static str) {
    match system {
        "m-kg-s-N-Pa" => ([0, 0, 0], "K"),
        "mm-t-s-N-MPa" | "unspecified" => ([-3, 3, 0], "°C"),
        "mm-kg-ms-kN-GPa" => ([-3, 0, -3], "°C"),
        "mm-g-ms-N-MPa" => ([-3, -3, -3], "°C"),
        // Not powers of ten: `to_si` scales it by its own factors.
        "in-lbf-s-psi" => ([0, 0, 0], "°F"),
        other => panic!("no format 1 fixture is in {other}"),
    }
}

/// What a number of a format 2 model (at `path` below the model) is: its dimension as powers
/// of length, mass, time and temperature, whether it is an absolute temperature, and whether it
/// was always typed in millimetres. Written from docs/units.md's table, apart from the step's
/// own map.
fn dimension_at(path: &[&str], model: &Value) -> Option<([i32; 4], bool, bool)> {
    const L: [i32; 4] = [1, 0, 0, 0];
    const PRESSURE: [i32; 4] = [-1, 1, -2, 0];
    let parent = |up: usize| -> &Value {
        path[..path.len() - up]
            .iter()
            .fold(model, |v, k| v.get(*k).or_else(|| k.parse::<usize>().ok().and_then(|i| v.get(i))).unwrap_or(&Value::Null))
    };
    let kind = |up: usize| parent(up).get("type").and_then(Value::as_str).unwrap_or("").to_string();
    let temperature = Some(([0, 0, 0, 1], true, false));
    let d = |dim: [i32; 4]| Some((dim, false, false));
    // Drop list indices at the end: a vector shares its entries' dimension.
    let mut p: Vec<&str> = path.to_vec();
    let mut indices = Vec::new();
    while p.last().is_some_and(|k| k.parse::<usize>().is_ok()) {
        indices.push(p.pop().unwrap().parse::<usize>().unwrap());
    }
    let up = indices.len();
    match p.as_slice() {
        ["nodes", _, "xyz"] | ["reference_points", _, "xyz"] => d(L),
        ["materials", _, "elastic", "young" | "shear"] => d(PRESSURE),
        ["materials", _, "density"] => d([-3, 1, 0, 0]),
        ["materials", _, "fluid", "viscosity"] => d([-1, 1, -1, 0]),
        ["materials", _, "hyperelastic", _] => d(PRESSURE),
        ["materials", _, "plastic", "yield_stress" | "hardening_modulus" | "a" | "b"] => d(PRESSURE),
        ["materials", _, "plastic", "points", _, "stress"] => d(PRESSURE),
        ["materials", _, "plastic", "reference_strain_rate"] => d([0, 0, -1, 0]),
        ["materials", _, "plastic", "room_temperature" | "melt_temperature"] => temperature,
        ["materials", _, "plastic", "specific_heat"] | ["materials", _, "thermal", "specific_heat"] => d([2, 0, -2, -1]),
        ["materials", _, "thermal", "conductivity"] => d([1, 1, -3, -1]),
        ["materials", _, "thermal", "expansion"] => d([0, 0, 0, -1]),
        ["materials", _, "thermal", "reference_temperature"] => temperature,
        ["sections", _, "kind", "thickness"] | ["sections", _, "kind", "profile", _] => d(L),
        ["sections", _, "kind", "area"] => d([2, 0, 0, 0]),
        ["contacts", _, "friction", "decay"] | ["steps", _, "procedure", "rigid_walls", _, "friction", "decay"] => d([-1, 0, 1, 0]),
        ["ties", _, "position_tolerance"] => d(L),
        ["masses", _, "mass"] => d([0, 1, 0, 0]),
        ["masses", _, "inertia"] => d([2, 1, 0, 0]),
        ["springs", name, "stiffness" | "coefficient"] => {
            let torsional = model["springs"][*name]["direction"]["dof"].as_u64().is_some_and(|d| d >= 4);
            let damping = p.last() == Some(&"coefficient");
            d(match (torsional, damping) {
                (false, false) => [0, 1, -2, 0],
                (true, false) => [2, 1, -2, 0],
                (false, true) => [0, 1, -1, 0],
                (true, true) => [2, 1, -1, 0],
            })
        }
        ["boundary_conditions", _, "values"] => (indices.first().is_some_and(|i| *i < 3)).then_some((L, false, false)),
        ["boundary_conditions", _, "value"] => temperature,
        ["loads", _, "vector"] => d(match kind(up + 1).as_str() {
            "force" => [1, 1, -2, 0],
            "moment" => [2, 1, -2, 0],
            "traction" => PRESSURE,
            "body_force" => [-2, 1, -2, 0],
            other => panic!("a {other} load has no vector"),
        }),
        ["loads", _, "magnitude"] => d(PRESSURE),
        ["loads", _, "acceleration"] | ["fluid_domains", _, "acceleration"] | ["fluid_domains", _, "thermal", "gravity"] => d([1, 0, -2, 0]),
        ["loads", _, "angular_velocity"] => d([0, 0, -1, 0]),
        ["loads", _, "point"] => d(L),
        ["loads", _, "flux"] => d([0, 1, -3, 0]),
        ["loads", _, "film"] => d([0, 1, -3, -1]),
        ["loads", _, "power"] => d([-1, 1, -3, 0]),
        ["loads", _, "sink" | "value"] => temperature,
        ["fluid_domains", _, "min" | "max" | "spacing" | "tolerance"] => d(L),
        ["fluid_domains", _, "faces", _, "velocity"] | ["fluid_domains", _, "initial_velocity"] => d([1, 0, -1, 0]),
        ["fluid_domains", _, "faces", _, "pressure"] => d(PRESSURE),
        ["fluid_domains", _, "thermal", "initial_temperature"] => temperature,
        ["fluid_domains", _, "thermal", "expansion"] => d([0, 0, 0, -1]),
        ["fluid_domains", _, "thermal", "faces", _, "value"] => {
            if kind(1) == "flux" {
                d([0, 1, -3, 0])
            } else {
                temperature
            }
        }
        ["steps", _, "procedure", "min_hz" | "max_hz"] => d([0, 0, -1, 0]),
        ["steps", _, "procedure", "end_time"] => d([0, 0, 1, 0]),
        ["steps", _, "procedure", "initial_velocities", _, "velocity"] => d([1, 0, -1, 0]),
        ["steps", _, "procedure", "rigid_walls", _, "point"] => d(L),
        ["steps", _, "procedure", "transient", "end_time" | "increment"] => d([0, 0, 1, 0]),
        ["steps", _, "procedure", "transient", "initial_temperature"] => temperature,
        ["datums", _, "point"] => Some((L, false, true)),
        ["seeds", _, "size"] | ["size_controls", _, "size"] => Some((L, false, true)),
        ["mesh_plans", _, "size" | "tolerance"] | ["mesh_plans", _, "cuts", _, "point"] => Some((L, false, true)),
        ["mesh_plans", _, "rings", _, "center" | "radius" | "frame"] => Some((L, false, true)),
        ["relaxations", _, "cells", _, "at" | "moved"] => d(L),
        ["size_controls", _, "region", "box", "min" | "max"] | ["size_controls", _, "region", "sphere", "center" | "radius"] => {
            Some((L, false, true))
        }
        _ => None,
    }
}

/// The numbers of a format 2 model in SI, by the system it declared.
fn to_si(model: &mut Value, system: &str, in_mm: &std::collections::BTreeSet<String>, built_on_geometry: bool) {
    let snapshot = model.clone();
    let ([lb, mb, tb], temperature) = base(system);
    fn walk(v: &mut Value, path: &mut Vec<String>, f: &mut dyn FnMut(&[&str], &mut Value)) {
        match v {
            Value::Object(m) => m.iter_mut().for_each(|(k, x)| {
                path.push(k.clone());
                walk(x, path, f);
                path.pop();
            }),
            Value::Array(a) => a.iter_mut().enumerate().for_each(|(i, x)| {
                path.push(i.to_string());
                walk(x, path, f);
                path.pop();
            }),
            Value::Number(_) => f(&path.iter().map(String::as_str).collect::<Vec<_>>(), v),
            _ => {}
        }
    }
    walk(model, &mut Vec::new(), &mut |path, v| {
        let Some(([l, m, t, _], absolute, millimetre)) = dimension_at(path, &snapshot) else { return };
        // Lengths typed against the geometry are millimetres in a file built on it (else the
        // declared unit), and so is a node of a mesh made from it.
        let millimetre = (millimetre && built_on_geometry) || (path.first() == Some(&"nodes") && path.get(1).is_some_and(|n| in_mm.contains(*n)));
        let x = v.as_f64().unwrap();
        let si = if temperature == "°F" {
            inch_pound_to_si(x, [l, m, t, path_temperature(path, &snapshot)], absolute, millimetre)
        } else {
            let k = if millimetre { -3 } else { l * lb + m * mb + t * tb };
            let si = scaled(x, k);
            if absolute && temperature == "°C" { celsius_to_kelvin(si) } else { si }
        };
        *v = json!(si);
    });
}

/// Where `a` and `b` first differ, with both values there.
fn first_difference(a: &Value, b: &Value, at: &str) -> String {
    match (a, b) {
        (Value::Object(x), Value::Object(y)) => {
            for k in x.keys().chain(y.keys()) {
                if x.get(k) != y.get(k) {
                    return first_difference(x.get(k).unwrap_or(&Value::Null), y.get(k).unwrap_or(&Value::Null), &format!("{at}/{k}"));
                }
            }
            format!("{at}: (same)")
        }
        (Value::Array(x), Value::Array(y)) if x.len() == y.len() => x
            .iter()
            .zip(y)
            .enumerate()
            .find(|(_, (p, q))| p != q)
            .map_or(format!("{at}: (same)"), |(i, (p, q))| first_difference(p, q, &format!("{at}/{i}"))),
        _ => format!("{at}: {a} vs {b}"),
    }
}

/// The temperature power of the number at `path` (its dimension's fourth entry).
fn path_temperature(path: &[&str], model: &Value) -> i32 {
    dimension_at(path, model).map_or(0, |(d, ..)| d[3])
}

/// A number of an in-lbf-s-psi model in SI: inches (0.0254 m), lbf·s²/in (the pound-force,
/// 4.4482216152605 N, per inch per s²), seconds and °F steps (5/9 K); an absolute temperature
/// (F + 459.67) × 5/9; a millimetre field from millimetres.
fn inch_pound_to_si(x: f64, [l, m, t, th]: [i32; 4], absolute: bool, millimetre: bool) -> f64 {
    const INCH: f64 = 0.0254;
    const MASS: f64 = 4.4482216152605 / INCH;
    if millimetre {
        return scaled(x, -3);
    }
    if absolute {
        let n = (x * 100.0).round() as i64 + 45967;
        // In hundredths: one correctly rounded division of exact integers.
        return (n * 5) as f64 / 900.0;
    }
    let _ = t;
    let factor = INCH.powi(l) * MASS.powi(m) * (5.0f64 / 9.0).powi(th);
    if factor == 1.0 { x } else { x * factor }
}

/// `x × 10^k`: on the decimal digits `x` was written with when there are at most 15 of them (a
/// typed number), else in binary by an exact power of ten.
fn scaled(x: f64, k: i32) -> f64 {
    if x == 0.0 || k == 0 {
        return x;
    }
    let text = format!("{x:e}");
    let (digits, exponent) = text.split_once('e').unwrap();
    let significant = digits.chars().filter(char::is_ascii_digit).count();
    if significant <= 15 {
        return format!("{digits}e{}", exponent.parse::<i32>().unwrap() + k).parse().unwrap();
    }
    if k > 0 { x * 10f64.powi(k) } else { x / 10f64.powi(-k) }
}

/// t °C in kelvin, added in decimal (t is a typed number).
fn celsius_to_kelvin(t: f64) -> f64 {
    // Both have at most two decimals here; hundredths are exact integers.
    let hundredths = (t * 100.0).round() as i64 + 27315;
    format!("{}e-2", hundredths).parse().unwrap()
}

/// The current reader's result in the record's shape, as raw JSON.
fn as_recorded(file: &ModelFile) -> Value {
    json!({
        "model": serde_json::to_value(&file.model).unwrap(),
        "geometry": serde_json::to_value(&file.geometry).unwrap(),
        "journal": {
            "units": serde_json::to_value(&file.journal_header.units).unwrap(),
            "model_units": serde_json::to_value(&file.journal_header.model_units).unwrap(),
            "commands": file.journal.iter().map(|(c, p)| json!({"command": c, "params": p})).collect::<Vec<_>>(),
        },
        "view": serde_json::to_value(&file.view).unwrap(),
    })
}

/// A mesh made from the millimetre design in a model of metres or inches (alone, or beside a
/// block typed in metres) opens with a note saying what the release that wrote it did, and now
/// overlays the design: its bounding box is the design's, in metres; the block stays in metres.
#[test]
fn geometry_meshes_of_other_lengths_migrate_from_millimetres() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats/v1/0.7.0");
    for (name, system, unit, times) in [
        ("design-metres", "m-kg-s-N-Pa", "metres", "1000×"),
        ("design-inches", "in-lbf-s-psi", "inches", "25.4×"),
        ("mixed-metres", "m-kg-s-N-Pa", "metres", "1000×"),
    ] {
        let f = fieldcraft_io::read(&std::fs::read(dir.join(format!("{name}.fieldcraft"))).unwrap()).unwrap();
        let note = f.notes.join(" ");
        for want in [system, &format!("solved those millimetre coordinates as {unit}"), times, "`Case`", "overlays the design", "results will differ"]
        {
            assert!(note.contains(want), "{name}: `{want}` in {note}");
        }
        // The hollow box is 40 × 30 × 22 mm from the origin.
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for label in f.model.body_elements("Case") {
            for n in &f.model.elements()[&label].nodes {
                let p = f.model.node(*n).unwrap().xyz;
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        assert_eq!((lo, hi), ([0.0; 3], [0.04, 0.03, 0.022]), "{name}");
    }
    // The block typed in metres stays in metres: 40 mm long at z = 0.1 m.
    let f = fieldcraft_io::read(&std::fs::read(dir.join("mixed-metres.fieldcraft")).unwrap()).unwrap();
    assert_eq!(f.model.node(1).unwrap().xyz, [0.0, 0.0, 0.1]);
}

#[test]
fn every_format_1_fixture_opens_as_its_writer_read_it() {
    let all = fixtures(1);
    let writers: std::collections::BTreeSet<&str> = all.iter().map(|f| f.writer.as_str()).collect();
    assert_eq!(writers.into_iter().collect::<Vec<_>>(), ["0.4.0", "0.5.0", "0.6.0", "0.7.0", "develop-2eb3515c"]);
    for f in &all {
        let name = f.name();
        let bytes = std::fs::read(&f.path).unwrap();
        assert_eq!(format_of(Cursor::new(&bytes)), Some(1), "{name}");
        let (file, format) = read_versioned(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(format, 1, "{name}");
        let (got, want) = (as_recorded(&file), expected(f.recorded.clone(), 1));
        for key in ["model", "geometry", "journal", "view"] {
            assert!(got[key] == want[key], "{name}: {key} differs from what its writer read at {}", first_difference(&got[key], &want[key], ""));
        }
        assert_eq!(got, want, "{name}");
    }
}

/// The step 1 → 2 on every real format 1 file: only the tag changes; every part's content is
/// the same.
#[test]
fn the_step_from_1_to_2_changes_only_the_tag() {
    for f in fixtures(1) {
        let name = f.name();
        let ps = parts(&std::fs::read(&f.path).unwrap());
        let json = |part: &str| -> Value { serde_json::from_slice(&ps[part]).unwrap() };
        let before = Parts {
            manifest: json("fieldcraft.json"),
            model: json("model.json"),
            mesh: ps["mesh.bin"].clone(),
            journal: json("journal.json"),
            notes: Vec::new(),
        };
        let mut after = before.clone();
        migrate::step(&mut after, 1).unwrap();
        assert_eq!(after.manifest["format"], "fieldcraft/2", "{name}");
        let mut manifest = after.manifest.clone();
        manifest["format"] = before.manifest["format"].clone();
        assert_eq!(
            (&manifest, &after.model, &after.mesh, &after.journal),
            (&before.manifest, &before.model, &before.mesh, &before.journal),
            "{name}"
        );
    }
}

/// The step 2 → 3 on every real file (through 1 → 2): `model.json` comes out in SI as the
/// independent transform of [`to_si`] says, `mesh.bin` as MCMESH02 with every coordinate scaled
/// from the declared length unit and every other byte as it was, the journal's commands untouched
/// under a header naming the system they were typed in, and the tag 3.
#[test]
fn the_step_from_2_to_3_scales_to_si_and_keeps_the_journal() {
    for f in fixtures(1) {
        let name = f.name();
        let ps = parts(&std::fs::read(&f.path).unwrap());
        let json = |part: &str| -> Value { serde_json::from_slice(&ps[part]).unwrap() };
        let mut p = Parts {
            manifest: json("fieldcraft.json"),
            model: json("model.json"),
            mesh: ps["mesh.bin"].clone(),
            journal: json("journal.json"),
            notes: Vec::new(),
        };
        migrate::step(&mut p, 1).unwrap();
        let before = p.clone();
        migrate::step(&mut p, 2).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(p.manifest["format"], "fieldcraft/3", "{name}");
        let mut manifest = p.manifest.clone();
        manifest["format"] = before.manifest["format"].clone();
        assert_eq!(manifest, before.manifest, "{name}: only the tag changes in the manifest");

        let declared = before.model["units"].as_str().unwrap_or("unspecified").to_string();
        let shown = if declared == "unspecified" { "mm-t-s-N-MPa" } else { declared.as_str() };
        let mut want = before.model.clone();
        want.as_object_mut().unwrap().remove("units");
        want["display_units"] = json!(shown);
        let in_mm = millimetre_nodes(&f.recorded, &declared);
        to_si(&mut want, &declared, &in_mm, built_on_geometry(&f.recorded["geometry"]));
        fluid_zero(&mut want, &declared);
        assert!(p.model == want, "{name}: model.json differs at {}", first_difference(&p.model, &want, ""));

        assert_eq!(&p.mesh[..8], b"MCMESH02", "{name}");
        assert_eq!(p.mesh.len(), before.mesh.len(), "{name}");
        let nodes = u64::from_le_bytes(before.mesh[8..16].try_into().unwrap()) as usize;
        let ([length, ..], _) = base(&declared);
        for i in 0..nodes {
            let at = 48 + i * 28;
            assert_eq!(p.mesh[at..at + 4], before.mesh[at..at + 4], "{name}: label");
            let label = u32::from_le_bytes(before.mesh[at..at + 4].try_into().unwrap());
            for c in 0..3 {
                let x = |b: &[u8]| f64::from_le_bytes(b[at + 4 + 8 * c..at + 12 + 8 * c].try_into().unwrap());
                let want = if in_mm.contains(&label.to_string()) {
                    scaled(x(&before.mesh), -3)
                } else if declared == "in-lbf-s-psi" {
                    x(&before.mesh) * 0.0254
                } else {
                    scaled(x(&before.mesh), length)
                };
                assert_eq!(x(&p.mesh).to_bits(), want.to_bits(), "{name}: node {i}");
            }
        }
        assert_eq!(p.mesh[48 + nodes * 28..], before.mesh[48 + nodes * 28..], "{name}: elements and tags");
        assert_eq!(p.journal, json!({"units": shown, "model_units": "mm-t-s-N-MPa", "commands": before.journal}), "{name}");
    }
}

/// Saving a migrated model writes the current format; reading that gives the same model, and
/// saving it again the same bytes.
#[test]
fn migrated_fixtures_save_and_reopen_bit_identical() {
    for f in fixtures(1) {
        let name = f.name();
        let (file, _) = read_versioned(&std::fs::read(&f.path).unwrap()).unwrap();
        let saved = write(&file.model, &file.geometry, &file.journal, &file.journal_header, &file.view).unwrap();
        assert_eq!(format_of(Cursor::new(&saved)), Some(FORMAT_VERSION), "{name}");
        let (again, format) = read_versioned(&saved).unwrap();
        assert_eq!(format, FORMAT_VERSION, "{name}");
        // The same model; the migration's note is the old file's, not the saved one's.
        assert!(again.notes.is_empty(), "{name}");
        assert_eq!(again, fieldcraft_io::ModelFile { notes: Vec::new(), ..file }, "{name}");
        assert_eq!(write(&again.model, &again.geometry, &again.journal, &again.journal_header, &again.view).unwrap(), saved, "{name}: bytes");
        let manifest: Value = serde_json::from_slice(&parts(&saved)["fieldcraft.json"]).unwrap();
        assert_eq!(manifest["format"], FORMAT);
    }
}

fn parts(bytes: &[u8]) -> std::collections::BTreeMap<String, Vec<u8>> {
    use std::io::Read;
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).unwrap();
    (0..zip.len())
        .map(|i| {
            let mut f = zip.by_index(i).unwrap();
            let mut b = Vec::new();
            f.read_to_end(&mut b).unwrap();
            (f.name().unwrap().to_string(), b)
        })
        .collect()
}

/// What `format_of` makes of things that are not a FieldCraft file, or are a later one.
#[test]
fn format_of_reads_only_the_tag() {
    assert_eq!(format_of(Cursor::new(b"not a zip")), None);
    let (file, _) = read_versioned(&std::fs::read(&fixtures(1)[0].path).unwrap()).unwrap();
    let saved = write(&file.model, &file.geometry, &file.journal, &file.journal_header, &file.view).unwrap();
    let mut ps = parts(&saved);
    let mut manifest: Value = serde_json::from_slice(&ps["fieldcraft.json"]).unwrap();
    for (tag, expect) in [("fieldcraft/9", Some(9)), ("fieldcraft/99999999999", Some(u32::MAX)), ("other/1", None)] {
        manifest["format"] = Value::from(tag);
        ps.insert("fieldcraft.json".into(), serde_json::to_vec(&manifest).unwrap());
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for (name, b) in &ps {
            zip.start_file(name.as_str(), zip::write::SimpleFileOptions::default()).unwrap();
            std::io::Write::write_all(&mut zip, b).unwrap();
        }
        let bytes = zip.finish().unwrap().into_inner();
        assert_eq!(format_of(Cursor::new(&bytes)), expect, "{tag}");
    }
}

/// A fresh folder holding a copy of the format 1 fixture `name` as `bracket.fieldcraft`.
fn folder_with_v1(tag: &str, name: &str) -> (PathBuf, PathBuf) {
    let dir = std::env::temp_dir().join(format!("fieldcraft-formats-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let model = dir.join("bracket.fieldcraft");
    std::fs::copy(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats/v1/0.7.0").join(name), &model).unwrap();
    (dir, model)
}

fn listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().into_owned()).collect();
    names.sort();
    names
}

#[cfg(not(target_arch = "wasm32"))]
#[test]
fn saving_over_an_older_format_keeps_a_copy_once() {
    use fieldcraft_io::fs::{backup_path, format_of_file, save_model};
    let (dir, model) = folder_with_v1("backup", "hex-static.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    let (file, _) = read_versioned(&original).unwrap();
    let bytes = write(&file.model, &file.geometry, &file.journal, &file.journal_header, &file.view).unwrap();
    assert_eq!(format_of_file(&model), Some(1));

    let kept = save_model(&model, &bytes).unwrap().expect("a copy of the format 1 file");
    assert_eq!(kept, dir.join("bracket.v1.fieldcraft"));
    assert_eq!(std::fs::read(&kept).unwrap(), original, "the copy is the original, byte for byte");
    assert_eq!(std::fs::read(&model).unwrap(), bytes);
    assert_eq!(format_of_file(&model), Some(FORMAT_VERSION));

    // The next save replaces `.bak` as always, but makes no copy and leaves the first one.
    assert_eq!(save_model(&model, &bytes).unwrap(), None);
    assert_eq!(std::fs::read(&kept).unwrap(), original);
    assert_eq!(std::fs::read(backup_path(&model)).unwrap(), bytes);
    assert_eq!(listing(&dir), ["bracket.fieldcraft", "bracket.fieldcraft.bak", "bracket.v1.fieldcraft"]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// A name already taken (an earlier copy, or a file of the user's) is never replaced: the copy
/// takes the next free name.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn a_format_backup_never_replaces_a_file() {
    use fieldcraft_io::fs::save_model;
    let (dir, model) = folder_with_v1("collision", "every-kind.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    std::fs::write(dir.join("bracket.v1.fieldcraft"), b"the user's own file").unwrap();
    std::fs::write(dir.join("bracket.v1-2.fieldcraft"), b"another").unwrap();
    let (file, _) = read_versioned(&original).unwrap();
    let bytes = write(&file.model, &file.geometry, &file.journal, &file.journal_header, &file.view).unwrap();
    let kept = save_model(&model, &bytes).unwrap().unwrap();
    assert_eq!(kept, dir.join("bracket.v1-3.fieldcraft"));
    assert_eq!(std::fs::read(&kept).unwrap(), original);
    assert_eq!(std::fs::read(dir.join("bracket.v1.fieldcraft")).unwrap(), b"the user's own file");
    assert_eq!(std::fs::read(dir.join("bracket.v1-2.fieldcraft")).unwrap(), b"another");

    // A later format's file is kept the same way before this version writes over it; a file
    // that is not a model, or a current one, needs no copy.
    let later = dir.join("later.fieldcraft");
    let mut ps = parts(&bytes);
    let mut manifest: Value = serde_json::from_slice(&ps["fieldcraft.json"]).unwrap();
    manifest["format"] = Value::from(format!("fieldcraft/{}", FORMAT_VERSION + 1));
    ps.insert("fieldcraft.json".into(), serde_json::to_vec(&manifest).unwrap());
    let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for (name, b) in &ps {
        zip.start_file(name.as_str(), zip::write::SimpleFileOptions::default()).unwrap();
        std::io::Write::write_all(&mut zip, b).unwrap();
    }
    std::fs::write(&later, zip.finish().unwrap().into_inner()).unwrap();
    assert_eq!(save_model(&later, &bytes).unwrap(), Some(dir.join(format!("later.v{}.fieldcraft", FORMAT_VERSION + 1))));
    // Anything that is not positively a current model is kept: a file that is not a model at
    // all (a Save As over it), a package whose manifest cannot be read, a malformed tag.
    let plain = dir.join("notes.txt");
    std::fs::write(&plain, b"the user's notes").unwrap();
    let kept = save_model(&plain, &bytes).unwrap().unwrap();
    assert_eq!((kept.file_name().unwrap().to_str().unwrap(), std::fs::read(&kept).unwrap()), ("notes.backup.txt", b"the user's notes".to_vec()));
    for (name, manifest) in
        [("unreadable", b"{not json".to_vec()), ("malformed", br#"{"format": "fieldcraft/01"}"#.to_vec()), ("untagged", b"{}".to_vec())]
    {
        let mut ps = parts(&bytes);
        ps.insert("fieldcraft.json".into(), manifest);
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for (part, b) in &ps {
            zip.start_file(part.as_str(), zip::write::SimpleFileOptions::default()).unwrap();
            std::io::Write::write_all(&mut zip, b).unwrap();
        }
        let odd = zip.finish().unwrap().into_inner();
        let path = dir.join(format!("{name}.fieldcraft"));
        std::fs::write(&path, &odd).unwrap();
        let kept = save_model(&path, &bytes).unwrap().unwrap_or_else(|| panic!("{name}: no copy"));
        assert_eq!(kept, dir.join(format!("{name}.backup.fieldcraft")));
        assert_eq!(std::fs::read(&kept).unwrap(), odd, "{name}");
    }
    // A current model, or no file at all, needs no copy.
    assert_eq!(save_model(&dir.join("new.fieldcraft"), &bytes).unwrap(), None);
    assert_eq!(save_model(&dir.join("new.fieldcraft"), &bytes).unwrap(), None);
    let _ = std::fs::remove_dir_all(&dir);
}

/// When the copy cannot be made, nothing is saved: the original stays as it was.
#[cfg(unix)]
#[test]
fn no_copy_no_save() {
    use std::os::unix::fs::PermissionsExt;
    let (dir, model) = folder_with_v1("readonly", "hex-static.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o555)).unwrap();
    let r = fieldcraft_io::fs::save_model(&model, b"new bytes");
    std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o755)).unwrap();
    // Running as root writes anyway; then the copy was made and the save went through.
    if let Err(e) = r {
        assert!(e.to_string().contains("could not be kept"), "{e}");
        assert_eq!(std::fs::read(&model).unwrap(), original);
        assert_eq!(listing(&dir), ["bracket.fieldcraft"]);
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// No export, deck or picture replaces a model or a design, whatever its name.
#[cfg(not(target_arch = "wasm32"))]
#[test]
fn outputs_never_replace_a_model_or_a_design() {
    use fieldcraft_io::fs::{protected_kind, write_output};
    let (dir, model) = folder_with_v1("outputs", "hex-static.fieldcraft");
    let original = std::fs::read(&model).unwrap();
    let renamed = dir.join("deck.inp");
    std::fs::copy(&model, &renamed).unwrap();
    let design = dir.join("plate.json");
    std::fs::write(&design, b"{\n  \"format\": \"solidcraft/21\",\n  \"features\": []\n}").unwrap();
    let named_design = dir.join("empty.SolidCraft");
    std::fs::write(&named_design, b"").unwrap();
    let named_model = dir.join("blank.FIELDCRAFT");
    std::fs::write(&named_model, b"").unwrap();
    for (path, kind) in [
        (&model, "FieldCraft model"),
        (&renamed, "FieldCraft model"),
        (&named_model, "FieldCraft model"),
        (&design, "SolidCraft design"),
        (&named_design, "SolidCraft design"),
    ] {
        assert_eq!(protected_kind(path), Some(kind), "{}", path.display());
        let before = std::fs::read(path).unwrap();
        let e = write_output(path, b"*HEADING").unwrap_err();
        assert!(e.to_string().contains(kind) && e.to_string().contains("choose another file name"), "{e}");
        assert_eq!(std::fs::read(path).unwrap(), before);
    }
    assert_eq!(std::fs::read(&model).unwrap(), original);
    // Other files, and new ones, are written as before.
    let notes = dir.join("notes.inp");
    std::fs::write(&notes, b"{\"format\": \"something/1\"}").unwrap();
    assert_eq!(protected_kind(&notes), None);
    write_output(&notes, b"*HEADING").unwrap();
    write_output(&dir.join("new.inp"), b"*HEADING").unwrap();
    assert_eq!(std::fs::read(&notes).unwrap(), b"*HEADING");
    assert_eq!(protected_kind(&dir.join("missing.fieldcraft")), None);
    let _ = std::fs::remove_dir_all(&dir);
}
