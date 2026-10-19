//! The geometry commands and the meshers through the display system: they work in the
//! geometry's millimetres, and take and give everything in the display units.

use serde_json::{Value, json};

use crate::Session;

/// A session with the design the SolidCraft `commands` build open.
fn design(commands: Value) -> (crate::tests::TempDir, Session) {
    let bytes = fieldcraft_geometry::script_design(&json!({"commands": commands}).to_string()).unwrap();
    let (dir, path) = crate::tests::design_file(&bytes);
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": path})).unwrap();
    (dir, s)
}

const INCH: &str = "in-lbf-s-psi";
const SI: &str = "m-kg-s-N-Pa";
const MM: &str = "mm-t-s-N-MPa";

/// A 4 × 2 × 1 inch block (101.6 × 50.8 × 25.4 mm), shown in `system`.
fn block(system: &str) -> (crate::tests::TempDir, Session) {
    let (dir, mut s) = design(json!([{"command": "solid.box", "params": {"length": 101.6, "width": 50.8, "height": 25.4, "body_name": "Block"}}]));
    s.execute("model.units", &json!({"system": system})).unwrap();
    (dir, s)
}

/// One inch (25.4 mm) in `system`'s length unit.
fn inch_in(system: &str) -> f64 {
    match system {
        INCH => 1.0,
        SI => 0.0254,
        _ => 25.4,
    }
}

#[test]
fn a_model_shown_in_inches_or_si_meshes_as_in_millimetres() {
    let mesh = |system: &str, command: &str| {
        let (_dir, mut s) = block(system);
        let h = inch_in(system);
        let r = s.execute(command, &json!({"body": "Block", "size": h / 2.0, "tolerance": h / 100.0})).unwrap();
        (s.model.nodes().clone(), s.model.elements().clone(), r)
    };
    for command in ["mesh.tet.generate", "mesh.hex.generate"] {
        let (nodes, elements, in_mm) = mesh(MM, command);
        assert!(!elements.is_empty(), "{command}");
        for system in [INCH, SI] {
            let (n, e, r) = mesh(system, command);
            assert_eq!((&n, &e), (&nodes, &elements), "{command} in {system}");
            // The reply's volume is the same volume, in the display units.
            let unit = if system == INCH { "in³" } else { "m³" };
            assert_eq!(r["units"]["volume"], unit, "{command}: {r}");
            let v = |r: &Value| r["bodies"][0]["volume"].as_f64().unwrap();
            let expected = if system == INCH { 8.0 } else { 8.0 * 0.0254f64.powi(3) };
            assert!((v(&r) - expected).abs() < 1e-9 * expected, "{command} in {system}: {r}");
            assert!((v(&in_mm) - 8.0 * 25.4f64.powi(3)).abs() < 1e-6, "{in_mm}");
        }
    }
}

#[test]
fn datums_seeds_partitions_and_plans_round_trip_in_inches() {
    let (_dir, mut s) = block(INCH);
    // A datum typed in inches is that point, and its edit form gives the inches back.
    let r = s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 2}})).unwrap();
    assert_eq!((r["point"].clone(), r["units"]["length"].clone()), (json!([2.0, 0.0, 0.0]), json!("in")));
    let stored = *s.model.datum("Mid").unwrap();
    let point = match stored {
        fieldcraft_model::Datum::Plane { point, .. } => point,
        _ => panic!("not a plane"),
    };
    assert!((point[0] - 0.0508).abs() < 1e-15, "{point:?}");
    let (command, params) = crate::cmd::datum_edit_params("Mid", &stored);
    let shown = super::shown_params(command, &params, s.units()).unwrap();
    assert_eq!(shown["by"]["point"], json!([2.0, 0.0, 0.0]));
    s.execute(command, &shown).unwrap();
    assert_eq!(s.model.datum("Mid"), Some(&stored));
    // The partition cuts the block in two halves of 2 inches.
    s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
    let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
    assert_eq!(bodies.as_array().unwrap().len(), 2, "{bodies}");
    for b in bodies.as_array().unwrap() {
        let width = b["bbox"][1][0].as_f64().unwrap() - b["bbox"][0][0].as_f64().unwrap();
        assert!((width - 2.0).abs() < 1e-9, "{b}");
    }
    // A seed's size typed in inches lists as typed.
    let edge = bodies[0]["edges"][0]["name"].as_str().unwrap().to_string();
    let body = bodies[0]["name"].as_str().unwrap().to_string();
    s.execute("mesh.seed.edge", &json!({"edges": [{"body": body, "edge": edge}], "size": 0.25})).unwrap();
    let seeds = s.execute("mesh.seed.list", &json!({})).unwrap();
    assert_eq!(seeds[0]["size"], json!(0.25), "{seeds}");
    // A plan in inches, applied as the object it replied, is the plan.
    let (_dir, mut s) = block(INCH);
    let mut plan = s.execute("mesh.plan", &json!({"body": "Block", "size": 0.5})).unwrap();
    assert_eq!(plan["size"], json!(0.5), "{plan}");
    assert!(plan["notes"].as_array().into_iter().flatten().all(|n| !n.as_str().unwrap_or_default().contains("mm")), "{plan}");
    let saved = s.model.mesh_plan("Plan-Block").unwrap().clone();
    for k in ["name", "decomposition", "units"] {
        plan.as_object_mut().unwrap().remove(k);
    }
    s.execute("mesh.plan.apply", &json!({"plan": plan, "name": "Again", "symmetry": false})).unwrap();
    let again = s.model.mesh_plan("Again").unwrap();
    assert!((again.size - saved.size).abs() <= 1e-15 * saved.size, "{} {}", again.size, saved.size);
    assert_eq!(again.cuts.len(), saved.cuts.len());
}

/// A journal of a model from before SI storage built on a design typed the numbers of its
/// geometry commands in millimetres, whatever its system: replayed, it gives the same model; a
/// datum typed into it now is recorded with its unit, so the journal replays to the model again.
/// Without a design they followed the declared unit, as the stored datums migrate.
#[test]
fn a_journal_from_before_si_storage_reads_its_millimetres_on_a_design() {
    let bytes = fieldcraft_geometry::script_design(
        &json!({"commands": [{"command": "solid.box", "params": {"length": 20, "width": 20, "height": 20, "body_name": "Block"}}]}).to_string(),
    )
    .unwrap();
    let (_dir, path) = crate::tests::design_file(&bytes);
    let point = |p: f64| json!({"command": "geometry.datum.point", "params": {"name": "P", "by": {"type": "coordinates", "point": [p, 0, 0]}}});
    let script = json!({"units": INCH, "model_units": MM, "commands": [{"command": "geometry.open", "params": {"path": path}}, point(10.0)]});
    let mut s = Session::default();
    s.run_script(&script).unwrap();
    assert_eq!(s.model.datum("P"), Some(&fieldcraft_model::Datum::Point { point: [0.01, 0.0, 0.0] }));
    // The form shows the journal's entry in the display units (millimetres here).
    let (c, p) = s.journal.last().cloned().unwrap();
    assert_eq!(s.journal_params_shown(&c, &p)["by"]["point"], json!([10.0, 0.0, 0.0]));
    // Typed now, in the display system, it is not millimetres.
    s.execute("model.units", &json!({"system": INCH})).unwrap();
    s.execute("geometry.datum.point", &json!({"name": "Q", "by": {"type": "coordinates", "point": [1, 0, 0]}})).unwrap();
    assert_eq!(s.model.datum("Q"), Some(&fieldcraft_model::Datum::Point { point: [0.0254, 0.0, 0.0] }));
    let commands: Vec<Value> = s.journal.iter().map(|(c, p)| json!({"command": c, "params": p})).collect();
    let mut replayed = Session::default();
    let h = &s.journal_header;
    replayed.run_script(&json!({"units": h.units, "model_units": h.model_units, "commands": commands})).unwrap();
    assert_eq!(replayed.model.datums(), s.model.datums());

    // Without a design, in the journal's inches.
    let mut free = Session::default();
    free.run_script(&json!({"units": INCH, "model_units": INCH, "commands": [point(10.0)]})).unwrap();
    assert_eq!(free.model.datum("P"), Some(&fieldcraft_model::Datum::Point { point: [0.254, 0.0, 0.0] }));
}

/// Every string a geometry or mesh command gives back (its reply's text, its errors) names
/// lengths in the display units: run in inches, none says millimetres.
#[test]
fn geometry_and_mesh_replies_name_no_fixed_unit() {
    fn strings(v: &Value, out: &mut Vec<String>) {
        match v {
            Value::String(s) => out.push(s.clone()),
            Value::Array(a) => a.iter().for_each(|x| strings(x, out)),
            // A feature's `definition` is the design's own text, as its author wrote it.
            Value::Object(o) => o.iter().filter(|(k, _)| *k != "definition").for_each(|(_, x)| strings(x, out)),
            _ => {}
        }
    }
    let (_dir, mut s) = design(json!([
        {"command": "solid.box", "params": {"length": 60, "width": 40, "height": 6, "body_name": "Plate"}},
        {"command": "solid.cylinder", "params": {"center": [30, 20, 0], "diameter": 10, "height": 6, "operation": "cut", "targets": ["Plate"]}},
    ]));
    s.execute("model.units", &json!({"system": INCH})).unwrap();
    let mut texts = Vec::new();
    let runs = [
        ("geometry.bodies", json!({})),
        ("geometry.features", json!({})),
        ("geometry.check", json!({"size": 0.2})),
        ("geometry.defeature.propose", json!({"body": "Plate", "below": {"hole_diameter": 1}})),
        ("geometry.datum", json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 0.5}})),
        ("mesh.plan", json!({"body": "Plate"})),
        ("mesh.tet.generate", json!({"body": "Plate", "size": 0.2})),
        ("mesh.quality", json!({})),
        ("mesh.check.duplicates", json!({})),
        ("mesh.node.merge", json!({})),
        ("mesh.size.set", json!({"name": "Fine", "size": 0.1, "sphere": {"center": [1, 1, 0], "radius": 0.2}})),
        ("mesh.size.list", json!({})),
        // Refused: their messages too.
        ("mesh.tet.generate", json!({"body": "Plate", "tolerance": 1e-9})),
        ("mesh.hex.generate", json!({"body": "Plate", "size": 1e-9})),
        ("mesh.hex.project", json!({"body": "Plate", "size": 1e-9})),
        ("mesh.hex.octree", json!({"body": "Plate", "size": 1e-9})),
        // The meshers' own refusals: a weld limit, an octree grid limit, a point.
        ("mesh.surface.check", json!({"body": "Plate", "weld_distance": 10})),
        ("mesh.tet.generate", json!({"body": "Plate", "weld_distance": 10})),
        ("mesh.hex.octree", json!({"body": "Plate", "size": 0.01, "max_depth": 1, "max_elements": 10})),
        ("mesh.planar", json!({"polygon": [[0, 0], [2, 0], [1, 1], [2, 2], [0, 2], [1, 1]], "size": 0.1})),
    ];
    for (command, params) in runs {
        match s.execute(command, &params) {
            Ok(r) => strings(&r, &mut texts),
            Err(e) => texts.push(e.to_string()),
        }
    }
    let millimetre = regex_free_word(&texts, "mm");
    assert!(millimetre.is_empty(), "named millimetres in an inch display: {millimetre:#?}");
    // The meshers' refusals name their numbers in inches.
    let said = |what: &str| texts.iter().any(|t| t.contains(what) && t.contains(" in"));
    for what in ["weld distance", "a size of", "coincide at ("] {
        assert!(said(what), "no `{what}` in inches: {texts:#?}");
    }
}

/// A mesher's refusal keeps its numbers as data: shown at the boundary in the display units,
/// whatever the frame it was made in (the volume gate, a point, a size).
#[test]
fn mesher_errors_are_shown_in_the_display_units() {
    use fieldcraft_mesh::message::{length, point, volume};
    let inch = fieldcraft_model::DisplaySystem::new(fieldcraft_model::units::Preset::InLbfS);
    let gate = fieldcraft_mesh::Error::failed(fieldcraft_mesh::message![
        "mesh validity gate: mesh volume ",
        volume(16387.064),
        " differs from the boundary's ",
        volume(32774.128)
    ]);
    let e = crate::EngineError::context("body `Block`", gate).shown(super::geometry(), &inch).to_string();
    assert_eq!(e, "body `Block`: meshing failed: mesh validity gate: mesh volume 1 in³ differs from the boundary's 2 in³");
    let near = fieldcraft_mesh::Error::invalid(fieldcraft_mesh::message!["size ", length(-25.4), " at ", point([25.4, 0.0, 50.8])]);
    assert_eq!(crate::EngineError::Mesh(near).shown(super::geometry(), &inch).to_string(), "invalid input: size -1 in at (1, 0, 2) in");
    // A command in SI shows its SI numbers in the display units too.
    let si = fieldcraft_mesh::Error::invalid(fieldcraft_mesh::message!["at ", point([0.0254, 0.0, 0.0])]);
    assert_eq!(crate::EngineError::Mesh(si).shown(super::si(), &inch).to_string(), "invalid input: at (1, 0, 0) in");
}

/// The texts that hold `word` (or `word²`, `word³`) as a word of its own.
fn regex_free_word(texts: &[String], word: &str) -> Vec<String> {
    texts
        .iter()
        .filter(|t| {
            t.match_indices(word).any(|(i, _)| {
                let before = t[..i].chars().next_back();
                let after = t[i + word.len()..].chars().next();
                !before.is_some_and(char::is_alphanumeric) && !after.is_some_and(|c| c.is_alphanumeric() || c == '-')
            })
        })
        .cloned()
        .collect()
}

/// No string literal of the geometry and mesh commands, the geometry crate or the replies the
/// UI draws for them names millimetres: their numbers go through the display system.
#[test]
fn geometry_and_mesh_sources_name_no_fixed_unit() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let cmd = "engine/src/cmd";
    let mut files: Vec<std::path::PathBuf> = [
        "geometry.rs",
        "datum.rs",
        "defeature.rs",
        "stitch.rs",
        "virtual_faces.rs",
        "mesh.rs",
        "hex.rs",
        "octree.rs",
        "extrude.rs",
        "project.rs",
        "cells.rs",
        "plan.rs",
        "seeds.rs",
        "sizes.rs",
        "mirror.rs",
        "nodes.rs",
        "node_move.rs",
        "relax.rs",
        "symmetry.rs",
        "replace.rs",
        "plan/blocks.rs",
        "plan/decompose.rs",
        "plan/pillow.rs",
        "cells/transition.rs",
    ]
    .iter()
    .map(|f| root.join(cmd).join(f))
    .collect();
    files.extend(["ui-egui/src/reply.rs", "ui-egui/src/feature_graph.rs"].map(|f| root.join(f)));
    // The meshers, whose messages keep their numbers as data.
    let mut dirs = vec![root.join("mesh/src")];
    while let Some(d) = dirs.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                dirs.push(p);
            } else if p.extension().is_some_and(|x| x == "rs") && !p.file_name().is_some_and(|n| n == "tests.rs") {
                files.push(p);
            }
        }
    }
    // The geometry crate, but for its 3MF unit table and the SolidCraft scripts of its samples.
    for e in std::fs::read_dir(root.join("geometry/src")).unwrap() {
        let p = e.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        if p.extension().is_some_and(|x| x == "rs") && !["tagged_mesh.rs", "sample.rs"].contains(&name.as_str()) {
            files.push(p);
        }
    }
    let mut found = Vec::new();
    for f in &files {
        let text = std::fs::read_to_string(f).unwrap();
        // The code before its tests.
        let code = text.split("#[cfg(test)]").next().unwrap_or_default();
        for (n, line) in code.lines().enumerate() {
            if line.trim_start().starts_with("//") {
                continue;
            }
            // A unit written as a `Typed` amount's own ("mm", "mm^3") is the amount's unit, not text.
            let literals: Vec<String> =
                line.split('"').skip(1).step_by(2).filter(|l| !matches!(*l, "mm" | "mm^2" | "mm^3")).map(str::to_string).collect();
            for l in regex_free_word(&literals, "mm") {
                found.push(format!("{}:{}: {l}", f.display(), n + 1));
            }
        }
    }
    assert!(found.is_empty(), "literal millimetres in what geometry and mesh commands say: {found:#?}");
}
