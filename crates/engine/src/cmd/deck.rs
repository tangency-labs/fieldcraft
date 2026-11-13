//! Solver decks: write model data, or read a flat Abaqus deck into a new session.

use std::path::Path;

use fieldcraft_io::roots::Access;
use serde_json::{Value, json};

use super::{opt, req};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const WRITE: CommandSpec = CommandSpec {
    id: "deck.write",
    label: "Write Solver Deck",
    params: r#"{"format": "abaqus|nastran|lsdyna", "path": "part.inp", "steps": ["Load"], "units": "mm-t-s-N-MPa|m-kg-s-N-Pa|mm-kg-ms-kN-GPa|mm-g-ms-N-MPa|in-lbf-s-psi"} — writes the model as a solver deck in one consistent unit system, stated in its header (default: the display system's, ignoring its overrides) (formats: fieldcraft_decks::FORMATS); `steps` (optional) writes only those steps, as a Nastran run needs (static or frequency steps, not both); LS-DYNA writes several steps as *CASE analyses (run it with CASE); fails naming every error of the model check and everything the format cannot express"#,
    form: Some(super::forms::DECK),
    kind: CommandKind::Session,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection],
    run: write,
};

pub const READ: CommandSpec = CommandSpec {
    id: "deck.read",
    label: "Read Abaqus Deck",
    params: r#"{"format":"abaqus","path":"part.inp","units"?:"mm-t-s-N-MPa|m-kg-s-N-Pa|mm-kg-ms-kN-GPa|mm-g-ms-N-MPa|in-lbf-s-psi"} — replaces the model with a flat Abaqus deck's mesh, sets, element surfaces, elastic/density materials and sections; input numbers use this consistent preset (default: the current display preset, ignoring overrides), then store SI; the reply states units, unread keywords with original lines (first 100 and total count), and warnings. No steps, loads, conditions or unknown text are retained; INCLUDE/parts/assemblies and unsupported physical variants are refused. Empty supported data is refused. Success clears old results/design/path/undo and is unsaved; errors leave the whole session intact."#,
    form: Some(super::forms::DECK_READ),
    kind: CommandKind::Session,
    requires: &[],
    run: read,
};

pub const OPENFOAM: CommandSpec = CommandSpec {
    id: "deck.openfoam",
    label: "Write OpenFOAM Case",
    params: r#"{"path": "case_dir", "step"?: "Flow", "cell_factor"?: 2, "processes"?: 1, "max_courant"?: 0.8, "refine_boxes"?: [{"min": [x,y,z], "max": [x,y,z]}], "location"?: [x,y,z]} — writes a fluid step (the only one when `step` is left out) as an OpenFOAM case directory: blockMesh background of spacing × cell_factor (1, 2, 4, 8), snappyHexMesh around the obstacles (one patch per named face, refined back to the spacing, as are `refine_boxes`), 0/U p nut, pimpleFoam with an adaptive time step (max_courant), Smagorinsky LES, decomposeParDict for `processes` ranks and an Allrun script; lengths in the input units; `location` is a point in the fluid (found when left out). Fails naming what OpenFOAM cannot be given yet (heat, porous zones, fans, baffles, free surface, periodic faces, parabolic inlets, non-Newtonian fluids). Writes SI; nothing of the model changes"#,
    form: None,
    kind: CommandKind::Session,
    requires: &[Requirement::FluidDomain, Requirement::Step],
    run: openfoam,
};

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct BoxParam {
    min: [f64; 3],
    max: [f64; 3],
}

fn openfoam(s: &mut Session, p: &Value) -> Result<Value> {
    use fieldcraft_decks::openfoam::{CaseOptions, ObstacleSurface, RefineBox, write_case};
    use fieldcraft_model::units::Quantity;

    let dir = super::file_param(s, p, "path", Access::Write)?;
    let step: Option<String> = opt(p, "step")?;
    let length = |key: &str, v: [f64; 3]| -> Result<[f64; 3]> {
        let mut out = [0.0; 3];
        for (o, x) in out.iter_mut().zip(v) {
            *o = s.input_units().to_si(Quantity::Length, x).map_err(|e| EngineError::invalid(key, e.to_string()))?;
        }
        Ok(out)
    };
    let defaults = CaseOptions::default();
    let mut options = CaseOptions {
        cell_factor: opt(p, "cell_factor")?.unwrap_or(defaults.cell_factor),
        processes: opt(p, "processes")?.unwrap_or(defaults.processes),
        max_courant: opt(p, "max_courant")?.unwrap_or(defaults.max_courant),
        ..defaults
    };
    for b in opt::<Vec<BoxParam>>(p, "refine_boxes")?.unwrap_or_default() {
        options.refine_boxes.push(RefineBox { min: length("refine_boxes", b.min)?, max: length("refine_boxes", b.max)? });
    }
    if let Some(l) = opt::<[f64; 3]>(p, "location")? {
        options.location = Some(length("location", l)?);
    }

    // The obstacles of the step's domain, tessellated as the lattice solver does, in metres.
    let domain = match step.as_deref() {
        Some(name) => s.model.steps().iter().find(|st| st.name == name),
        None => s.model.steps().iter().find(|st| matches!(st.procedure, fieldcraft_model::Procedure::Fluid { .. })),
    }
    .and_then(|st| st.procedure.fluid_domain())
    .and_then(|d| s.model.fluid_domains().get(d));
    let mut meshes = Vec::new();
    if let Some(domain) = domain {
        let tolerance = super::solve::obstacle_tolerance(domain);
        for body in &domain.obstacles {
            let (mesh, faces) = super::solve::obstacle_mesh(s, body, tolerance)?;
            let positions: Vec<[f64; 3]> = mesh.positions.iter().map(|q| crate::units::geometry_point(*q)).collect();
            meshes.push((body.as_str(), positions, mesh.into_owned(), faces));
        }
    }
    let surfaces: Vec<ObstacleSurface<'_>> = meshes
        .iter()
        .map(|(body, positions, mesh, faces)| ObstacleSurface { body, positions, triangles: &mesh.triangles, tri_face: &mesh.tri_face, faces })
        .collect();
    let case = write_case(&s.model, step.as_deref(), &surfaces, &options)?;

    let root = Path::new(&dir);
    let mut bytes = 0;
    for (rel, text) in &case.files {
        let path = root.join(rel);
        if let Some(parent) = path.parent() {
            fieldcraft_io::fs::create_dir_all(parent)?;
        }
        fieldcraft_io::fs::write_output(&path, text.as_bytes())?;
        bytes += text.len();
    }
    Ok(json!({"format": "openfoam", "path": dir, "solver": case.solver, "files": case.files.keys().collect::<Vec<_>>(), "bytes": bytes,
        "background_cells": case.background_cells, "location": case.location, "patches": case.patches, "notes": case.notes,
        "run": format!("cd {dir} && sh Allrun")}))
}

fn read(s: &mut Session, p: &Value) -> Result<Value> {
    let format: &str = req(p, "format")?;
    if format != "abaqus" {
        return Err(EngineError::invalid("format", "only abaqus is supported by deck.read"));
    }
    let units = match opt::<String>(p, "units")? {
        Some(label) => fieldcraft_model::units::Preset::from_label(&label)
            .ok_or_else(|| EngineError::invalid("units", "choose a known consistent unit preset"))?,
        None => s.units().preset(),
    };
    let path = super::file_param(s, p, "path", Access::Read)?;
    let bytes = super::file::read_file(&path, fieldcraft_decks::abaqus::read::MAX_INPUT_BYTES as u64)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::invalid("path", "the input deck is not UTF-8 text"))?;
    let limits = fieldcraft_decks::abaqus::Limits {
        nodes: s.limits.max_import_nodes,
        elements: s.limits.max_import_elements,
        members: s.limits.max_import_members,
    };
    let imported = fieldcraft_decks::abaqus::read_limited(text, units, limits)?;
    let model = &imported.model;
    if model.nodes().is_empty()
        && model.elements().is_empty()
        && model.sets().values().all(fieldcraft_model::Set::is_empty)
        && model.materials().is_empty()
        && model.sections().is_empty()
    {
        return Err(EngineError::invalid(
            "path",
            format!("the deck contains no supported model data ({} unread keywords); the current model is unchanged", imported.unknown_count),
        ));
    }
    let mut kinds = std::collections::BTreeMap::<&str, usize>::new();
    for element in model.elements().values() {
        *kinds.entry(element.kind.name()).or_default() += 1;
    }
    let unknown: Vec<Value> = imported.unknown.iter().map(|u| json!({"keyword":u.keyword,"line":u.line})).collect();
    let reply = json!({"format":"abaqus","path":path,"units":units.label(),"nodes":model.nodes().len(),"elements":kinds,
        "sets":model.sets().len(),"materials":model.materials().len(),"sections":model.sections().len(),"unknown":unknown,
        "unknown_count":imported.unknown_count,"warnings":imported.warnings,"history_imported":false});
    // The new, unsaved document has no journal, but still needs its own dirty revision.
    let opened = Session { model: imported.model, revision: crate::next_revision(), ..Session::default() };
    s.replace_with(opened);
    Ok(reply)
}

fn write(s: &mut Session, p: &Value) -> Result<Value> {
    let name: &str = req(p, "format")?;
    let format = fieldcraft_decks::format(name).ok_or_else(|| {
        let shown: String = name.chars().take(32).collect();
        let known: Vec<&str> = fieldcraft_decks::FORMATS.iter().map(|f| f.name).collect();
        EngineError::Params(format!("unknown deck format `{shown}`; known: {}", known.join(", ")))
    })?;
    let path = super::file_param(s, p, "path", Access::Write)?;
    let steps: Option<Vec<String>> = super::opt(p, "steps")?;
    let units = match p.get("units").and_then(Value::as_str) {
        None => s.units().preset(),
        Some(label) => fieldcraft_model::units::Preset::from_label(label).ok_or_else(|| {
            let known: Vec<&str> = fieldcraft_model::units::Preset::ALL.iter().map(|u| u.label()).collect();
            EngineError::invalid("units", format!("must be one of {}", known.join(", ")))
        })?,
    };
    let deck = format.write_steps(&s.model, steps.as_deref(), units)?;
    fieldcraft_io::fs::write_output(Path::new(&path), deck.text.as_bytes())?;
    Ok(json!({"format": format.name, "path": path, "bytes": deck.text.len(), "units": units.label(), "warnings": deck.warnings}))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sample plate in a box of water flowing along x, written as an OpenFOAM case: a patch
    /// per named face of the plate, lengths in metres, nothing of the model changed.
    #[test]
    fn writes_a_fluid_step_as_an_openfoam_case() {
        use fieldcraft_geometry::sample::{Plate, plate_design};
        let dir = std::env::temp_dir().join(format!("fieldcraft-openfoam-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, plate_design(&Plate::default()).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": design.to_str().unwrap()})).unwrap();
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let [lo, hi] = [0, 1].map(|k| bodies[0]["bbox"][k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>());
        let min: Vec<f64> = lo.iter().map(|v| (v / 5.0).floor() * 5.0 - 20.0).collect();
        let max: Vec<f64> = hi.iter().map(|v| (v / 5.0).ceil() * 5.0 + 20.0).collect();
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-9}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [1000, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0},
                          "z_min": {"type": "slip"}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10}}},
        ]);
        s.run_script(&script).unwrap();
        let journal = s.journal.len();
        let case = dir.join("case");
        let out = s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "processes": 4, "cell_factor": 1})).unwrap();
        assert_eq!(s.journal.len(), journal, "writing a case does not change the model");
        assert_eq!(out["solver"], "pimpleFoam");
        let read = |rel: &str| std::fs::read_to_string(case.join(rel)).unwrap();
        // Lengths in metres: the 5 mm spacing is the background cell.
        let extent: Vec<f64> = (0..3).map(|i| (max[i] - min[i]) / 5.0).collect();
        assert_eq!(out["background_cells"], json!(extent.iter().map(|n| n.round() as u64).collect::<Vec<_>>()), "{out}");
        assert!(read("system/blockMeshDict").contains(&format!("({:e} ", min[0] / 1000.0)), "{}", read("system/blockMeshDict"));
        // A solid and a patch per named face of the plate; m/s and m²/s.
        let faces = bodies[0]["faces"].as_array().unwrap();
        let stl = read("constant/triSurface/Body1.stl");
        assert_eq!(stl.matches("endsolid").count(), faces.len(), "{out}");
        assert_eq!(out["patches"].as_object().unwrap().len(), 6 + faces.len());
        assert!(read("0/U").contains("x_min { type fixedValue; value uniform (1e0 0 0); }"), "{}", read("0/U"));
        assert!(read("constant/transportProperties").contains("nu 1e-6;"));
        assert!(read("system/decomposeParDict").contains("numberOfSubdomains 4;"));
        assert!(read("Allrun").contains("snappyHexMesh -overwrite"));
        // What OpenFOAM cannot be given yet is named, and nothing is written.
        s.execute(
            "model.fluid_domain.edit",
            &json!({"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"], "acceleration": [0, 0, -9810]}),
        )
        .unwrap();
        let other = dir.join("other");
        let e = s.execute("deck.openfoam", &json!({"path": other.to_str().unwrap()})).unwrap_err().to_string();
        assert!(e.contains("acceleration"), "{e}");
        assert!(!other.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn writes_the_defined_model_and_refuses_what_it_cannot() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-deck-write-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bar.inp");
        let path_str = path.to_string_lossy().into_owned();

        let mut s = Session::default();
        crate::cmd::define::tests::defined_bar(&mut s);
        let journal = s.journal.len();
        let out = s.execute("deck.write", &json!({"format": "abaqus", "path": path_str})).unwrap();
        let written = std::fs::read_to_string(&path).unwrap();
        assert_eq!(written, fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MmTonneS).unwrap());
        assert_eq!(out["bytes"], written.len());
        assert_eq!(out["warnings"].as_array().map(Vec::len), Some(1), "`press` is in no step: {out}");
        assert!(written.contains("*ELEMENT, TYPE=C3D8I"), "{written}");
        assert_eq!(s.journal.len(), journal, "writing a deck does not change the model");
        assert!(std::fs::read_dir(&dir).unwrap().count() == 1, "no partial file is left behind");

        // Unknown formats, missing parameters and models that are not ready fail without writing.
        let other = dir.join("other.inp").to_string_lossy().into_owned();
        assert!(s.execute("deck.write", &json!({"format": "nope", "path": other})).unwrap_err().to_string().contains("known: abaqus"));
        assert!(s.execute("deck.write", &json!({"path": other})).is_err());
        assert!(s.execute("deck.write", &json!({"format": "abaqus"})).is_err());
        let mut empty = Session::default();
        assert!(matches!(empty.execute("deck.write", &json!({"format": "abaqus", "path": other})), Err(EngineError::Unavailable(_))));
        assert!(!dir.join("other.inp").exists());
        assert!(matches!(s.execute("deck.write", &json!({"format": "abaqus", "path": dir.to_string_lossy()})), Err(EngineError::File(_))));

        // Nastran: one run is static or frequency steps, so the bar's two steps need two decks.
        let bdf = dir.join("bar.bdf").to_string_lossy().into_owned();
        let e = s.execute("deck.write", &json!({"format": "nastran", "path": bdf})).unwrap_err().to_string();
        assert!(e.contains("not supported by Nastran") && e.contains("separate decks"), "{e}");
        for (step, sol) in [("Bend", "SOL 101"), ("Modes", "SOL 103")] {
            s.execute("deck.write", &json!({"format": "nastran", "path": bdf, "steps": [step]})).unwrap();
            let deck = std::fs::read_to_string(&bdf).unwrap();
            assert!(deck.contains(sol) && deck.contains("RBE2") && deck.ends_with("ENDDATA\n"), "{deck}");
        }
        assert_eq!(s.journal.len(), journal);
        for steps in [json!([]), json!(["Nope"]), json!("Bend"), json!([1, 2]), json!([["Bend"]])] {
            assert!(s.execute("deck.write", &json!({"format": "nastran", "path": bdf, "steps": steps})).is_err(), "{steps}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn writes_lsdyna_keyword_decks() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-deck-lsdyna-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bar.k");
        let mut s = Session::default();
        crate::cmd::define::tests::defined_bar(&mut s);
        let out = s.execute("deck.write", &json!({"format": "lsdyna", "path": path.to_string_lossy()})).unwrap();
        let written = std::fs::read_to_string(&path).unwrap();
        assert_eq!(written, fieldcraft_decks::lsdyna::write(&s.model, fieldcraft_model::units::Preset::MmTonneS).unwrap().text);
        assert_eq!(out["format"], "lsdyna");
        assert!(written.starts_with("*KEYWORD long=y\n") && written.contains("*CASE_BEGIN_2"), "two steps, two cases");
        // One chosen step: a single analysis, no cases.
        s.execute("deck.write", &json!({"format": "lsdyna", "path": path.to_string_lossy(), "steps": ["Modes"]})).unwrap();
        let one = std::fs::read_to_string(&path).unwrap();
        assert!(one.contains("*CONTROL_IMPLICIT_EIGENVALUE") && !one.contains("*CASE_BEGIN") && !one.contains("*CONTROL_IMPLICIT_SOLUTION"), "{one}");
        // What LS-DYNA cannot take is named and nothing is written.
        s.execute("model.load.create", &json!({"name": "heat", "type": "body_force", "set": "all", "vector": [0, 0, 1]})).unwrap();
        s.execute("model.step.activate", &json!({"step": "Bend", "loads": ["heat"]})).unwrap();
        let other = dir.join("heat.k");
        let e = s.execute("deck.write", &json!({"format": "lsdyna", "path": other.to_string_lossy()})).unwrap_err().to_string();
        assert!(e.contains("not supported by LS-DYNA") && e.contains("load `heat`"), "{e}");
        assert!(!other.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[cfg(test)]
#[path = "deck/read_tests.rs"]
mod read_tests;
