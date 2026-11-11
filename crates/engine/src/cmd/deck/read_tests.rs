use super::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);
fn folder() -> PathBuf {
    let path = std::env::temp_dir().join(format!("fieldcraft-deck-read-{}-{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn bar() -> Session {
    let mut s = Session::default();
    super::super::define::tests::tagged_bar(&mut s);
    s.run_script(&json!([
        {"command":"model.units","params":{"system":"mm-t-s-N-MPa"}},
        {"command":"model.set.create","params":{"name":"all","geometry":[{"body":"Bar"}]}},
        {"command":"model.set.create","params":{"name":"root","geometry":[{"body":"Bar","face":"left"}]}},
        {"command":"model.set.create","params":{"name":"end","geometry":[{"body":"Bar","face":"right"}]}},
        {"command":"model.material.create","params":{"name":"Steel","elastic":{"type":"isotropic","young":"210 GPa","poisson":0.3},"density":"7850 kg/m^3"}},
        {"command":"model.section.create","params":{"name":"BarSection","set":"all","material":"Steel","kind":{"type":"solid"},"formulation":{"incompatible_modes":true}}}
    ])).unwrap();
    constrain(&mut s);
    s
}

fn constrain(s: &mut Session) {
    s.run_script(&json!([
        {"command":"model.boundary_condition.create","params":{"name":"clamp","type":"fixed","target":{"set":"root"}}},
        {"command":"model.load.create","params":{"name":"push","type":"force","target":{"set":"end"},"vector":["0 N","0 N","-100 N"]}},
        {"command":"model.step.create","params":{"name":"Bend","procedure":{"type":"static"},"boundary_conditions":["clamp"],"loads":["push"]}}
    ]))
    .unwrap();
}

fn assert_unchanged(s: &Session, before: &Session) {
    assert_eq!(s.model, before.model);
    assert_eq!(s.document(), before.document());
    assert_eq!(s.journal, before.journal);
    assert_eq!(s.results_report(), before.results_report());
    assert_eq!(s.results.len(), before.results.len());
    for (name, (revision, result)) in &before.results {
        let (kept_revision, kept) = s.results.get(name).unwrap();
        assert_eq!(revision, kept_revision);
        assert!(std::sync::Arc::ptr_eq(result, kept));
    }
    assert_eq!(s.undo_depth(), before.undo_depth());
    assert_eq!(s.revision(), before.revision());
    assert_eq!(s.opened(), before.opened());
    assert_eq!(s.is_dirty(), before.is_dirty());
    assert_eq!(s.file_access, before.file_access);
    assert_eq!(s.view, before.view);
    assert_eq!(s.design.is_some(), before.design.is_some());
    if let (Some(kept), Some(original)) = (&s.design, &before.design) {
        assert!(std::sync::Arc::ptr_eq(kept, original));
    }
}

#[test]
fn read_replaces_only_after_success_and_writer_history_solves_the_same() {
    let dir = folder();
    let path = dir.join("original.inp");
    let mut s = bar();
    let original = s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap();
    s.path = Some(dir.join("previous.fieldcraft").to_string_lossy().into_owned());
    let result = s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap();
    assert_eq!(result["nodes"], 12);
    assert_eq!(result["elements"]["hex8"], 2);
    assert_eq!(result["steps"], 1);
    assert_eq!(result["boundary_conditions"], 1);
    assert_eq!(result["loads"], 1);
    assert_eq!(result["units"], "mm-t-s-N-MPa");
    assert_eq!(result["history_imported"], true);
    assert!(!s.model.steps().is_empty() && !s.model.loads().is_empty() && !s.model.boundary_conditions().is_empty());
    assert_eq!(s.undo_depth(), (0, 0));
    assert!(s.journal.is_empty() && s.path.is_none() && s.result("Bend").is_none() && s.is_dirty());
    assert!(s.design.is_none() && s.surfaces.is_empty());
    let resumed = s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    let before = original["max_displacement"]["value"].as_f64().unwrap();
    let after = resumed["max_displacement"]["value"].as_f64().unwrap();
    assert!((before - after).abs() <= before.abs() * 1e-12, "{before} vs {after}");
    assert_eq!(original["max_displacement"]["node"], resumed["max_displacement"]["node"]);
}

#[test]
fn writer_roundtrips_fixed_and_symmetry_conditions_with_shared_nodes() {
    let dir = folder();
    let path = dir.join("shared-boundary.inp");
    let mut s = bar();
    s.execute("model.boundary_condition.create", &json!({"name":"sym","type":"symmetry","target":{"set":"all"},"plane":"x"})).unwrap();
    s.execute("model.step.activate", &json!({"step":"Bend","boundary_conditions":["sym"]})).unwrap();
    let original = s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    s.execute("deck.write", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap();
    assert_eq!(s.model.steps().len(), 1);
    assert_eq!(s.model.steps()[0].boundary_conditions.len(), 2);
    assert!(
        s.model
            .boundary_conditions()
            .values()
            .any(|bc| matches!(bc, fieldcraft_model::BoundaryCondition::Fixed { target: fieldcraft_model::Target::Set(set) } if set == "root"))
    );
    assert!(s.model.boundary_conditions().values().any(|bc| matches!(bc, fieldcraft_model::BoundaryCondition::Symmetry { target: fieldcraft_model::Target::Set(set), plane: fieldcraft_model::Symmetry::X } if set == "all")));
    let resumed = s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    let before = original["max_displacement"]["value"].as_f64().unwrap();
    let after = resumed["max_displacement"]["value"].as_f64().unwrap();
    assert!((before - after).abs() <= before.abs() * 1e-12, "{before} vs {after}");
    assert_eq!(original["max_displacement"]["node"], resumed["max_displacement"]["node"]);
}

#[test]
fn command_built_four_node_total_force_roundtrips_through_multiple_steps() {
    let dir = folder();
    let path = dir.join("multi-step.inp");
    let mut s = bar();
    s.execute("model.step.create", &json!({"name":"Hold","procedure":{"type":"static"},"boundary_conditions":["clamp"],"loads":["push"]})).unwrap();
    s.execute("model.step.create", &json!({"name":"Modes","procedure":{"type":"frequency","modes":2},"boundary_conditions":["clamp"]})).unwrap();
    s.execute("deck.write", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap();
    assert_eq!(s.model.steps().iter().map(|step| step.name.as_str()).collect::<Vec<_>>(), ["Bend", "Hold", "Modes"]);
    assert_eq!(s.model.boundary_conditions().len(), 1);
    assert_eq!(s.model.loads().len(), 1);
    let Some(fieldcraft_model::Load::Force {
        target: fieldcraft_model::Target::Set(set),
        vector,
        distribution: fieldcraft_model::Distribution::PerNode,
    }) = s.model.loads().values().next()
    else {
        panic!("{:?}", s.model.loads());
    };
    assert_eq!(set, "end");
    assert_eq!(*vector, [0.0, 0.0, -25.0]);
    assert_eq!(s.model.set(set).map(|set| set.nodes.len()), Some(4));
    assert_eq!(s.model.steps()[0].loads, s.model.steps()[1].loads);
    assert!(s.execute("solve.static", &json!({"step":"Hold"})).is_ok());
}

#[test]
fn every_read_error_keeps_model_results_path_undo_and_unsaved_work() {
    let dir = folder();
    let path = dir.join("bad.inp");
    let mut s = bar();
    s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    s.execute("file.save_as", &json!({"path":dir.join("saved.fieldcraft"),"results":true})).unwrap();
    s.execute("mesh.node.add", &json!({"label":99,"xyz":["5 mm","0 mm","0 mm"]})).unwrap();
    assert!(s.is_dirty());
    let before = s.clone();
    for text in [
        "",
        "*HEADING\nTitle\n",
        "*UNKNOWN\ndata\n",
        "*NSET,NSET=Empty\n",
        "*NODE\n0,0,0,0",
        "*NODE\n1,0,0,0\n*ELEMENT,TYPE=C3D4\n1,1,2,3,4",
        "*INCLUDE,INPUT=secret.inp",
        "*STEP\n*NODE\n9,1,2,3",
    ] {
        std::fs::write(&path, text).unwrap();
        assert!(s.execute("deck.read", &json!({"format":"abaqus","path":path})).is_err(), "{text}");
        assert_unchanged(&s, &before);
    }
    std::fs::write(&path, [0xff, 0xfe]).unwrap();
    for params in [
        json!({"format":"abaqus","path":path}),
        json!({"format":"nastran","path":path}),
        json!({"format":"abaqus"}),
        json!({"format":"abaqus","path":dir}),
        json!({"format":"abaqus","path":"/dev/zero"}),
        json!({"format":"abaqus","path":"x".repeat(100_000)}),
        json!({"format":"abaqus","path":path,"units":"guess"}),
    ] {
        assert!(s.execute("deck.read", &params).is_err());
        assert_unchanged(&s, &before);
    }
    // The existing undo stack is still usable, not merely the same depth.
    s.execute("edit.undo", &json!({})).unwrap();
    assert!(s.model.node(99).is_none());
}

#[test]
fn cooperative_read_cancellation_preserves_the_whole_previous_session() {
    let dir = folder();
    let path = dir.join("original.inp");
    let mut s = bar();
    s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap();
    s.path = Some(dir.join("prior.fieldcraft").to_string_lossy().into_owned());
    let before = s.clone();
    let progress = fieldcraft_model::progress::Progress::new();
    progress.cancel();
    let error =
        fieldcraft_model::progress::scoped(Some(progress), || s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})))
            .unwrap_err();
    assert!(error.is_cancelled(), "{error}");
    assert_unchanged(&s, &before);
    // Neither cancellation nor its reporter poisons the next ordinary import.
    assert!(s.execute("deck.read", &json!({"format":"abaqus","path":path})).is_ok());
}

#[test]
fn read_roots_default_units_and_unknown_output_are_explicit() {
    let dir = folder();
    let allowed = dir.join("allowed");
    std::fs::create_dir(&allowed).unwrap();
    let path = allowed.join("mesh.inp");
    let outside = dir.join("outside.inp");
    let text = format!("*NODE\n1,1,0,0\n{}", "*IGNORED\n".repeat(120));
    std::fs::write(&path, &text).unwrap();
    std::fs::write(&outside, &text).unwrap();
    let mut s = bar();
    s.execute("model.units", &json!({"system":"in-lbf-s-psi","overrides":{"length":"mm"}})).unwrap();
    s.file_access = fieldcraft_io::roots::Roots::new(Some(&allowed), None).unwrap();
    let before = s.clone();
    assert!(s.execute("deck.read", &json!({"format":"abaqus","path":outside})).is_err());
    assert_unchanged(&s, &before);
    let reply = s.execute("deck.read", &json!({"format":"abaqus","path":path})).unwrap();
    assert_eq!(reply["units"], "in-lbf-s-psi");
    assert_eq!(reply["unknown_count"], 120);
    assert_eq!(reply["unknown"].as_array().unwrap().len(), 100);
    let inch = fieldcraft_model::units::Quantity::Length.parse("1 in", crate::units::si()).unwrap();
    assert_eq!(s.model.node(1).unwrap().xyz, [inch, 0., 0.]);
    assert_eq!(s.file_access, before.file_access);
}

#[test]
fn the_session_import_caps_bound_the_read() {
    let dir = folder();
    let path = dir.join("original.inp");
    let mut s = bar();
    s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap();
    for limits in [
        crate::Limits { max_import_nodes: 11, ..crate::Limits::default() },
        crate::Limits { max_import_elements: 1, ..crate::Limits::default() },
        crate::Limits { max_import_members: 4, ..crate::Limits::default() },
    ] {
        s.limits = limits;
        let before = s.clone();
        let error = s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap_err().to_string();
        assert!(error.contains("limit"), "{error}");
        assert_unchanged(&s, &before);
    }
    s.limits = crate::Limits { max_import_nodes: 12, max_import_elements: 2, ..crate::Limits::default() };
    assert!(s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).is_ok());
    assert_eq!(s.limits.max_import_nodes, 12, "the new document keeps the session's caps");
}
