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
        {"command":"model.section.create","params":{"name":"BarSection","set":"all","material":"Steel","kind":{"type":"solid"},"formulation":{"b_bar":true}}}
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
    assert_eq!(s.deck_source.is_some(), before.deck_source.is_some());
    if let (Some(now), Some(old)) = (&s.deck_source, &before.deck_source) {
        assert!(std::sync::Arc::ptr_eq(now, old));
    }
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
fn read_replaces_only_after_success_and_imported_bar_solves_the_same() {
    let dir = folder();
    let path = dir.join("original.inp");
    let mut s = bar();
    let original = s.execute("solve.static", &json!({"step":"Bend"})).unwrap();
    s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap();
    s.path = Some(dir.join("previous.fieldcraft").to_string_lossy().into_owned());
    let result = s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"})).unwrap();
    assert_eq!(result["nodes"], 12);
    assert_eq!(result["elements"]["hex8"], 2);
    assert_eq!(result["units"], "mm-t-s-N-MPa");
    assert_eq!(result["history_imported"], true);
    assert_eq!(result["analysis_complete"], true);
    assert_eq!(s.model.steps().len(), 1);
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

#[test]
fn real_commands_replay_unknown_text_history_and_include_bytes() {
    let dir = folder();
    std::fs::create_dir_all(dir.join("sub")).unwrap();
    let root = "** untouched header\r\n*INCLUDE,INPUT=sub/mesh.inc\r\n*Future Option,flag=YES\r\nabc, DEF\r\n*STEP,NAME=Load\r\n*STATIC\r\n0.1,1\r\n*CLOAD\r\n1,2,-5\r\n*END STEP";
    let fragment = "*NODE\n1,0,0,0\n2,1,0,0";
    std::fs::write(dir.join("input.inp"), root).unwrap();
    std::fs::write(dir.join("sub/mesh.inc"), fragment).unwrap();
    let mut session = Session::default();
    let reply = session.execute("deck.read", &json!({"format":"abaqus","path":dir.join("input.inp")})).unwrap();
    assert_eq!(reply["source_retained"], true);
    assert_eq!(reply["history_imported"], false);
    assert!(session.model.steps().is_empty());
    assert!(
        reply["unknown"].as_array().unwrap().iter().any(|u| u["keyword"] == "FUTUREOPTION" && u["original"].as_str().unwrap().contains("abc, DEF"))
    );
    let out = dir.join("out");
    std::fs::create_dir(&out).unwrap();
    session.execute("deck.source.write", &json!({"path":out.join("replay.inp")})).unwrap();
    assert_eq!(std::fs::read(out.join("replay.inp")).unwrap(), root.as_bytes());
    assert_eq!(std::fs::read(out.join("sub/mesh.inc")).unwrap(), fragment.as_bytes());
    let model = session.model.clone();
    session.execute("deck.read", &json!({"format":"abaqus","path":out.join("replay.inp")})).unwrap();
    assert_eq!(session.model, model);
    // Replaying originals never misrepresents a changed model as an edited source deck.
    session.execute("mesh.node.add", &json!({"label":3,"xyz":[0,0,1]})).unwrap();
    let reply = session.execute("deck.source.write", &json!({"path":out.join("replay.inp")})).unwrap();
    assert_eq!(reply["model_edits_exported"], false);
    assert_eq!(std::fs::read(out.join("replay.inp")).unwrap(), root.as_bytes());
    session.execute("file.new", &json!({})).unwrap();
    assert!(session.execute("deck.source.write", &json!({"path":out.join("empty.inp")})).is_err());
}

#[test]
fn include_failure_and_output_collision_do_not_replace_session_or_write_root() {
    let dir = folder();
    let path = dir.join("input.inp");
    let mut session = bar();
    let before = session.clone();
    for input in ["*INCLUDE,INPUT=../escape.inc", "*INCLUDE,INPUT=input.inp", "*INCLUDE,INPUT=missing.inc"] {
        std::fs::write(&path, input).unwrap();
        assert!(session.execute("deck.read", &json!({"format":"abaqus","path":path})).is_err());
        assert_unchanged(&session, &before);
    }
    std::fs::write(&path, "*INCLUDE,INPUT=mesh.inc\n").unwrap();
    std::fs::write(dir.join("mesh.inc"), "*NODE\n1,0,0,0\n").unwrap();
    session.execute("deck.read", &json!({"format":"abaqus","path":path})).unwrap();
    let out = dir.join("out");
    std::fs::create_dir(&out).unwrap();
    assert!(session.execute("deck.source.write", &json!({"path":out.join("mesh.inc")})).is_err());
    assert!(!out.join("mesh.inc").exists());
}

#[cfg(unix)]
#[test]
fn includes_cannot_escape_through_symlinks_even_when_session_is_unconfined() {
    let dir = folder();
    let outside = folder();
    std::fs::write(outside.join("nodes.inc"), "*NODE\n1,0,0,0\n").unwrap();
    std::os::unix::fs::symlink(outside.join("nodes.inc"), dir.join("alias.inc")).unwrap();
    std::fs::write(dir.join("root.inp"), "*INCLUDE,INPUT=alias.inc\n").unwrap();
    let mut session = Session::default();
    let error = session.execute("deck.read", &json!({"format":"abaqus","path":dir.join("root.inp")})).unwrap_err().to_string();
    assert!(error.contains("outside") && error.contains("alias.inc"), "{error}");
}

#[test]
fn commands_expand_part_include_and_replay_assembly_without_flattening_source() {
    let dir = folder();
    let part = "*PART,NAME=P\n*NODE,NSET=All\n1,0,0,0\n2,1,0,0\n*END PART\n";
    let root = "*INCLUDE,INPUT=part.inc\n*ASSEMBLY,NAME=A\n*INSTANCE,NAME=I,PART=P\n3,0,0\n*END INSTANCE\n*NSET,NSET=Pick,INSTANCE=I\n2\n*END ASSEMBLY\n*STEP,NAME=Modes\n*FREQUENCY\n6\n*END STEP\n";
    std::fs::write(dir.join("part.inc"), part).unwrap();
    std::fs::write(dir.join("root.inp"), root).unwrap();
    let mut session = Session::default();
    let reply = session.execute("deck.read", &json!({"format":"abaqus","path":dir.join("root.inp"),"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(session.model.node(2).unwrap().xyz, [4.0, 0.0, 0.0]);
    assert_eq!(session.model.set("PICK").unwrap().nodes, vec![2]);
    assert_eq!(reply["instances"][0]["name"], "I");
    assert_eq!(reply["history_imported"], true);
    let out = dir.join("copy");
    std::fs::create_dir(&out).unwrap();
    session.execute("deck.source.write", &json!({"path":out.join("root.inp")})).unwrap();
    assert_eq!(std::fs::read_to_string(out.join("root.inp")).unwrap(), root);
    assert_eq!(std::fs::read_to_string(out.join("part.inc")).unwrap(), part);
    let expected = session.model.clone();
    session.execute("deck.read", &json!({"format":"abaqus","path":out.join("root.inp"),"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(session.model, expected);
    let before = session.clone();
    std::fs::write(dir.join("bad.inp"), "*INCLUDE,INPUT=missing.inc\n").unwrap();
    assert!(session.execute("deck.read", &json!({"format":"abaqus","path":dir.join("bad.inp")})).is_err());
    assert_unchanged(&session, &before);
}

const AXIAL_BEAM: &str = "** Original handwritten FieldCraft import validation; SI: m, kg, s, N, Pa\n*NODE,NSET=All\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=B31,ELSET=Beam\n1,1,2\n*MATERIAL,NAME=Steel\n*ELASTIC\n200000000000,0.3\n*DENSITY\n7800\n*BEAM SECTION,ELSET=Beam,MATERIAL=Steel,SECTION=RECT\n0.1,0.1\n0,0,-1\n*AMPLITUDE,NAME=Ramp,DEFINITION=TABULAR\n0,0,1,2\n*BOUNDARY\n1,ENCASTRE\n*STEP,NAME=Pull,NLGEOM=NO\n*STATIC\n1,1\n*CLOAD,AMPLITUDE=Ramp\n2,1,500\n*OUTPUT,FIELD\n*NODE OUTPUT\nU,RF\n*ELEMENT OUTPUT\nS,E\n*END STEP\n";

#[test]
fn imported_beam_solves_and_saved_source_survives_reopen_without_original_files() {
    let dir = folder();
    let inp = dir.join("beam.inp");
    let saved = dir.join("beam.fieldcraft");
    std::fs::write(&inp, AXIAL_BEAM).unwrap();
    let mut s = Session::default();
    let imported = s.execute("deck.read", &json!({"format":"abaqus","path":inp,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(imported["history_imported"], true, "{imported}");
    assert_eq!(imported["analysis_complete"], true, "{imported}");
    let solved = s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let displacement = solved["max_displacement"]["value"].as_f64().unwrap();
    assert!((displacement - 5e-7).abs() < 5e-17, "{solved}");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    std::fs::remove_file(&inp).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    let original = dir.join("original.inp");
    s.execute("deck.source.write", &json!({"path":original})).unwrap();
    assert_eq!(std::fs::read_to_string(original).unwrap(), AXIAL_BEAM);
    let again = s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    assert_eq!(again["max_displacement"], solved["max_displacement"]);
    let normalized = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normalized})).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    let roundtrip = s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    assert_eq!(roundtrip["max_displacement"], solved["max_displacement"]);
}

#[test]
fn harmonic_beam_stress_is_refused_without_dropping_requested_output() {
    let dir = folder();
    let inp = dir.join("harmonic.inp");
    let source = AXIAL_BEAM
        .replace("*AMPLITUDE,NAME=Ramp,DEFINITION=TABULAR\n0,0,1,2\n", "")
        .replace("*STATIC\n1,1\n", "*STEADY STATE DYNAMICS,DIRECT\n0,0\n")
        .replace("*CLOAD,AMPLITUDE=Ramp\n2,1,500", "*CLOAD,REAL\n2,1,1000")
        + "*STEP,NAME=Matrices\n*MATRIX GENERATE,STIFFNESS,MASS\n*MATRIX OUTPUT,STIFFNESS,MASS\n*END STEP\n";
    std::fs::write(&inp, source).unwrap();
    let mut s = Session::default();
    let reply = s.execute("deck.read", &json!({"format":"abaqus","path":inp,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["imported_steps"], json!(["Pull", "Matrices"]));
    assert_eq!(reply["analysis_complete"], true);
    assert!(matches!(s.model.step("Matrices").unwrap().procedure, fieldcraft_model::Procedure::MatrixGeneration { .. }));
    let error = s.execute("solve.harmonic", &json!({"step":"Pull"})).unwrap_err();
    assert!(error.to_string().contains("harmonic beam section stress"), "{error}");
    assert!(s.result("Pull").is_none());
}

#[test]
fn imported_solid_bar_solves_direct_harmonic_without_rebuilding_the_step() {
    let dir = folder();
    let inp = dir.join("harmonic-solid.inp");
    let mut s = bar();
    s.execute("deck.write", &json!({"format":"abaqus","path":inp})).unwrap();
    let source = std::fs::read_to_string(&inp).unwrap();
    assert!(source.contains("*STATIC\n1.0, 1.0"), "{source}");
    let source = source.replace("*STATIC\n1.0, 1.0", "*STEADY STATE DYNAMICS,DIRECT\n0,0\n10,0");
    std::fs::write(&inp, source).unwrap();
    let reply = s.execute("deck.read", &json!({"format":"abaqus","path":inp,"units":"mm-t-s-N-MPa"})).unwrap();
    assert_eq!(reply["analysis_complete"], true);
    assert_eq!(reply["imported_steps"], json!(["Bend"]));
    let result = s.execute("solve.harmonic", &json!({"step":"Bend"})).unwrap();
    assert_eq!(result["frequencies"].as_array().unwrap().len(), 2);
    assert!(s.result("Bend").is_some());
}

#[test]
fn imported_discrete_oscillator_solves_after_save_reopen_and_normalized_write() {
    // Original analytical oscillator, f = sqrt(200 / 2) / (2 pi).
    let source = "*NODE\n7,0,0,0\n*ELEMENT,TYPE=MASS,ELSET=M\n10,7\n*ELEMENT,TYPE=SPRING1,ELSET=K\n11,7\n*MASS,ELSET=M\n2\n*SPRING,ELSET=K\n1\n200\n*BOUNDARY\n7,2,3\n*STEP,NAME=Modes\n*FREQUENCY\n1\n*END STEP\n";
    let dir = folder();
    let inp = dir.join("oscillator.inp");
    std::fs::write(&inp, source).unwrap();
    let mut s = Session::default();
    let reply = s.execute("deck.read", &json!({"format":"abaqus","path":inp,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["analysis_complete"], true, "{reply}");
    let expected = 10.0 / (2.0 * std::f64::consts::PI);
    let first = s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap();
    let frequency = first["modes"][0]["frequency"].as_f64().unwrap();
    assert!((frequency - expected).abs() < expected * 1e-10, "{first}");
    let saved = dir.join("oscillator.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    std::fs::remove_file(inp).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), source);
    let reopened = s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap();
    assert_eq!(reopened["modes"], first["modes"]);
    let normalized = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normalized})).unwrap();
    s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    let again = s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap();
    assert_eq!(again["modes"], first["modes"]);
}

#[test]
fn imported_random_oscillator_runs_modal_psd_and_survives_reopen() {
    // Original SDOF: m=2 kg, k=200 N/m, zeta=.05, unit force multiplier PSD.
    let source = "*NODE\n7,0,0,0\n*ELEMENT,TYPE=MASS,ELSET=M\n10,7\n*ELEMENT,TYPE=SPRING1,ELSET=K\n11,7\n*MASS,ELSET=M\n2\n*SPRING,ELSET=K\n1\n200\n*BOUNDARY\n7,2,3\n*PSD-DEFINITION,NAME=White,TYPE=FORCE\n1,0,0\n1,0,4\n*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY\n1\n*END STEP\n*STEP,NAME=Random,PERTURBATION\n*RANDOM RESPONSE\n0,4,101,3,1\n*MODAL DAMPING,DEFINITION=MODE NUMBERS\n1,1,0.05\n*CLOAD,LOAD CASE=2\n7,1,1\n*CORRELATION,TYPE=CORRELATED,PSD=White\n2,1\n*END STEP\n";
    let dir = folder();
    let inp = dir.join("random.inp");
    std::fs::write(&inp, source).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":inp,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    let r = s.execute("solve.random", &json!({"step":"Random"})).unwrap();
    assert_eq!(r["frequency_points"], 201);
    let samples = r["probes"][0]["samples"].as_array().unwrap();
    for sample in samples {
        let omega = 2.0 * std::f64::consts::PI * sample["frequency"].as_f64().unwrap();
        let expected = 1.0 / ((200.0 - 2.0 * omega * omega).powi(2) + (2.0 * omega).powi(2));
        let actual = sample["displacement_components_psd"][0].as_f64().unwrap();
        assert!((actual - expected).abs() < expected * 1e-8, "{sample}");
    }
    let saved = dir.join("random.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved})).unwrap();
    std::fs::remove_file(inp).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    let again = s.execute("solve.random", &json!({"step":"Random"})).unwrap();
    assert_eq!(again, r);
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), source);
}

#[test]
fn imported_matrices_export_theory_operators_and_survive_reopen() {
    // Two original axial bars: each E A/L=300, rho A L/3=10; C=2M+.1K, G=.05K.
    // Labels deliberately disagree with physical order; exported triangle uses labels.
    let source = "*NODE\n10,0,0,0\n20,2,0,0\n5,4,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,10,20\n2,20,5\n*MATERIAL,NAME=Steel\n*ELASTIC\n200,0.3\n*DENSITY\n5\n*DAMPING,ALPHA=2,BETA=0.1,STRUCTURAL=0.05\n*SOLID SECTION,ELSET=Bar,MATERIAL=Steel\n3\n*BOUNDARY\n10,1,3\n20,2,3\n5,2,3\n*STEP,NAME=Matrices,PERTURBATION\n*MATRIX GENERATE,STIFFNESS,MASS,VISCOUS DAMPING,STRUCTURAL DAMPING\n*MATRIX OUTPUT,STIFFNESS,MASS,VISCOUS DAMPING,STRUCTURAL DAMPING,FORMAT=MATRIX INPUT\n*END STEP\n";
    let dir = folder();
    let inp = dir.join("matrices.inp");
    std::fs::write(&inp, source).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":inp,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    let manifest = dir.join("operators.json");
    let result = s.execute("solve.matrices", &json!({"step":"Matrices","path":manifest})).unwrap();
    let degrees: std::collections::BTreeSet<(u32, u8)> = serde_json::from_value(result["degrees"].clone()).unwrap();
    assert_eq!(degrees, [(5, 1), (20, 1)].into_iter().collect());
    for (name, diagonal, off) in [("stiffness", 300., -300.), ("mass", 10., 5.), ("viscous_damping", 50., -20.), ("structural_damping", 15., -15.)] {
        let path = result["matrices"][name]["path"].as_str().unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        let rows: Vec<_> = text.lines().filter(|l| !l.starts_with("**")).collect();
        assert_eq!(rows.len(), 3);
        let mut seen = std::collections::BTreeSet::new();
        for line in rows {
            let row: Vec<_> = line.split(',').map(str::trim).collect();
            let left: (u32, u8) = (row[0].parse().unwrap(), row[1].parse().unwrap());
            let right: (u32, u8) = (row[2].parse().unwrap(), row[3].parse().unwrap());
            assert!(left >= right);
            assert!(seen.insert((left, right)));
            let expected = match (left, right) {
                ((5, 1), (5, 1)) => diagonal,
                ((20, 1), (20, 1)) => 2. * diagonal,
                ((20, 1), (5, 1)) => off,
                _ => panic!("unexpected coordinate"),
            };
            assert!((row[4].parse::<f64>().unwrap() - expected).abs() < 1e-10);
        }
    }
    let saved = dir.join("matrices.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved})).unwrap();
    std::fs::remove_file(inp).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    assert_eq!(s.execute("solve.matrices", &json!({"step":"Matrices","path":manifest})).unwrap(), result);
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), source);
}

#[test]
fn named_boundaries_reopen_replay_and_normalized_read_preserve_dofs() {
    use fieldcraft_model::{BoundaryCondition, Model, Target};
    let constraints = |m: &Model| {
        let step = m.step("Load").unwrap();
        let mut dofs = std::collections::BTreeSet::new();
        for name in &step.boundary_conditions {
            let BoundaryCondition::Displacement { target: Target::Set(set), values } = &m.boundary_conditions()[name] else { panic!() };
            for node in &m.set(set).unwrap().nodes {
                for (dof, value) in values.iter().enumerate() {
                    if let Some(value) = value {
                        assert_eq!(*value, 0.);
                        dofs.insert((*node, dof + 1));
                    }
                }
            }
        }
        dofs
    };
    for (kind, dofs) in [
        ("ENCASTRE", vec![1, 2, 3, 4, 5, 6]),
        ("PINNED", vec![1, 2, 3]),
        ("XSYMM", vec![1, 5, 6]),
        ("YSYMM", vec![2, 4, 6]),
        ("ZSYMM", vec![3, 4, 5]),
        ("XASYMM", vec![2, 3, 4]),
        ("YASYMM", vec![1, 3, 5]),
        ("ZASYMM", vec![1, 2, 6]),
    ] {
        let dir = folder();
        let input = dir.join("input.inp");
        let text = format!(
            "** own boundary fixture\n*NODE,NSET=Edge\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=B31,ELSET=Beam\n1,1,2\n*MATERIAL,NAME=M\n*ELASTIC\n1000,.3\n*BEAM SECTION,ELSET=Beam,MATERIAL=M,SECTION=RECT\n1,1\n0,1,0\n*STEP,NAME=Load\n*STATIC\n*BOUNDARY\nEdge,{kind}\n*END STEP\n"
        );
        std::fs::write(&input, &text).unwrap();
        let mut s = Session::default();
        let reply = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(reply["analysis_complete"], true, "{kind}: {reply}");
        let expected = dofs.iter().flat_map(|d| [(1, *d), (2, *d)]).collect();
        assert_eq!(constraints(&s.model), expected);
        let project = dir.join("model.fieldcraft");
        s.execute("file.save_as", &json!({"path":project})).unwrap();
        s.execute("file.new", &json!({})).unwrap();
        s.execute("file.open", &json!({"path":project})).unwrap();
        assert_eq!(constraints(&s.model), expected);
        let replay = dir.join("replay.inp");
        s.execute("deck.source.write", &json!({"path":replay})).unwrap();
        assert_eq!(std::fs::read(&replay).unwrap(), text.as_bytes());
        let normalized = dir.join("normalized.inp");
        s.execute("deck.write", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
        s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(constraints(&s.model), expected);
        std::fs::remove_dir_all(dir).unwrap();
    }
}

#[test]
fn unmapped_constraints_and_materials_survive_reopen_without_runnable_steps() {
    for (keyword, block) in [
        ("TIE", "*TIE,NAME=Joint\nSlave,Master\n"),
        ("COUPLING", "*COUPLING,CONSTRAINT NAME=Joint,REF NODE=1,SURFACE=Slave,ORIENTATION=Unmapped\n*KINEMATIC\n1,6\n"),
        ("PLASTIC", "*MATERIAL,NAME=M\n*ELASTIC\n1000,.3\n*PLASTIC\n10,0\n"),
        ("TRANSFORM", "*TRANSFORM,NSET=All,TYPE=R\n0,1,0,-1,0,0\n"),
    ] {
        let dir = folder();
        let input = dir.join("input.inp");
        let text = format!("*NODE,NSET=All\n1,0,0,0\n2,1,0,0\n{block}*STEP,NAME=Load\n*STATIC\n*BOUNDARY\nAll,XASYMM\n*END STEP\n");
        std::fs::write(&input, &text).unwrap();
        let mut s = Session::default();
        let reply = s.execute("deck.read", &json!({"format":"abaqus","path":input})).unwrap();
        assert_eq!(reply["analysis_complete"], false, "{keyword}: {reply}");
        assert!(s.model.steps().is_empty());
        assert!(
            reply["unknown"]
                .as_array()
                .unwrap()
                .iter()
                .any(|x| x["keyword"] == keyword && x["original"].as_str().unwrap().contains(&format!("*{keyword}")))
        );
        let project = dir.join("model.fieldcraft");
        s.execute("file.save_as", &json!({"path":project})).unwrap();
        s.execute("file.new", &json!({})).unwrap();
        s.execute("file.open", &json!({"path":project})).unwrap();
        assert!(s.model.steps().is_empty());
        let output = dir.join("replay.inp");
        s.execute("deck.source.write", &json!({"path":output})).unwrap();
        assert_eq!(std::fs::read(output).unwrap(), text.as_bytes());
        std::fs::remove_dir_all(dir).unwrap();
    }
}

// Original two-cube axial patch: disjoint interface nodes, unit area, length 2,
// E=1000, nu=0, one newton through a rigid top reference point: u_z=0.002.
fn constraint_input() -> String {
    let mut text = String::from("*NODE\n");
    let corners = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    for cube in 0..2 {
        for (i, p) in corners.iter().enumerate() {
            text += &format!("{},{},{},{}\n", cube * 8 + i + 1, p[0], p[1], p[2] + cube as f64);
        }
    }
    text + "99,.5,.5,2\n*ELEMENT,TYPE=C3D8,ELSET=All\n1,1,2,3,4,5,6,7,8\n2,9,10,11,12,13,14,15,16\n*NSET,NSET=Base\n1,2,3,4\n*NSET,NSET=Top\n13,14,15,16\n*NSET,NSET=RP\n99\n*SURFACE,TYPE=ELEMENT,NAME=Main\n1,S2\n*SURFACE,TYPE=ELEMENT,NAME=Secondary\n2,S1\n*SURFACE,TYPE=NODE,NAME=Pad\nTop,1\n*MATERIAL,NAME=M\n*ELASTIC\n1000,0\n*DENSITY\n1\n*SOLID SECTION,ELSET=All,MATERIAL=M\n*TIE,NAME=Glue,ADJUST=NO,POSITION TOLERANCE=.01\nSecondary,Main\n*COUPLING,CONSTRAINT NAME=Rigid,REF NODE=RP,SURFACE=Pad\n*KINEMATIC\n1,6\n*STEP,NAME=Pull\n*STATIC\n*BOUNDARY\nBase,1,3\n*CLOAD\nRP,3,1\n*END STEP\n"
}
fn native_constraint_patch() -> Session {
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    for (i, xyz) in [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [0, 0, 2],
        [1, 0, 2],
        [1, 1, 2],
        [0, 1, 2],
    ]
    .iter()
    .enumerate()
    {
        s.execute("mesh.node.add", &json!({"label":i+1,"xyz":xyz})).unwrap();
    }
    for (id, nodes) in [(1, vec![1, 2, 3, 4, 5, 6, 7, 8]), (2, vec![9, 10, 11, 12, 13, 14, 15, 16])] {
        s.execute("mesh.element.add", &json!({"label":id,"kind":"hex8","nodes":nodes})).unwrap();
    }
    s.run_script(&json!([
        {"command":"model.set.create","params":{"name":"All","elements":[1,2]}},
        {"command":"model.set.create","params":{"name":"Base","nodes":[1,2,3,4]}},
        {"command":"model.set.create","params":{"name":"Top","nodes":[13,14,15,16]}},
        {"command":"model.set.create","params":{"name":"Main","faces":[{"element":1,"face":2}]}},
        {"command":"model.set.create","params":{"name":"Secondary","faces":[{"element":2,"face":1}]}},
        {"command":"model.material.create","params":{"name":"M","elastic":{"type":"isotropic","young":1000,"poisson":0},"density":1}},
        {"command":"model.section.create","params":{"name":"S","set":"All","material":"M","kind":{"type":"solid"},"formulation":{"b_bar":true}}},
        {"command":"interaction.tie","params":{"name":"Glue","main":"Main","secondary":"Secondary","position_tolerance":0.01}},
        {"command":"model.reference_point.create","params":{"name":"RP","xyz":[0.5,0.5,2],"label":99}},
        {"command":"model.coupling.create","params":{"name":"Rigid","reference_point":"RP","set":"Top","kind":"kinematic"}},
        {"command":"model.boundary_condition.create","params":{"name":"Clamp","type":"displacement","target":{"set":"Base"},"values":[0,0,0,null,null,null]}},
        {"command":"model.load.create","params":{"name":"F","type":"force","target":{"reference_point":"RP"},"vector":[0,0,1]}},
        {"command":"model.step.create","params":{"name":"Pull","procedure":{"type":"static"},"boundary_conditions":["Clamp"],"loads":["F"]}}
    ])).unwrap();
    s
}
#[test]
fn imported_constraints_solve_reopen_replay_and_normalize_like_native_commands() {
    let dir = folder();
    let input = dir.join("input.inp");
    let text = constraint_input();
    std::fs::write(&input, &text).unwrap();
    let mut native = native_constraint_patch();
    let expected = native.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let expected = expected["max_displacement"]["value"].as_f64().unwrap();
    assert!((expected - 0.002).abs() < 1e-12, "native {expected}");
    let mut s = Session::default();
    let reply = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["analysis_complete"], true, "{reply}");
    assert_eq!(s.model.ties().len(), 1);
    assert_eq!(s.model.couplings().len(), 1);
    assert!(s.model.node(99).is_none());
    assert_eq!(s.model.reference_points().values().next().unwrap().label, 99);
    let saved = dir.join("saved.fieldcraft");
    for phase in 0..3 {
        let actual = s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
        let actual = actual["max_displacement"]["value"].as_f64().unwrap();
        assert!((actual - expected).abs() < 1e-12, "phase {phase}: {actual} vs {expected}");
        assert_constraint_results(&native, &s);
        println!("constraint patch phase={phase} native={expected:.15e} imported={actual:.15e} error={:.3e}", (actual - expected).abs());
        if phase == 0 {
            s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
            s.execute("file.new", &json!({})).unwrap();
            s.execute("file.open", &json!({"path":saved})).unwrap();
            let replay = dir.join("replay.inp");
            s.execute("deck.source.write", &json!({"path":replay})).unwrap();
            assert_eq!(std::fs::read(replay).unwrap(), text.as_bytes());
        } else if phase == 1 {
            let normalized = dir.join("normalized.inp");
            s.execute("deck.write", &json!({"path":normalized,"format":"abaqus"})).unwrap();
            let read = s.execute("deck.read", &json!({"path":normalized,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
            assert_eq!(read["analysis_complete"], true, "{read}");
        }
    }
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn imported_constraints_named_rejections_preserve_source_and_bad_references_are_atomic() {
    let dir = folder();
    let input = dir.join("input.inp");
    for (from, to, reason) in [
        ("ADJUST=NO", "ADJUST=YES", "ADJUST"),
        ("*STEP,NAME=Pull", "*ELEMENT,TYPE=MASS,ELSET=PointMass\n100,99\n*MASS,ELSET=PointMass\n1\n*STEP,NAME=Pull", "MASS/SPRING"),
        ("POSITION TOLERANCE=.01", "POSITION TOLERANCE=.01,NO ROTATION", "NOROTATION"),
        ("SURFACE=Pad", "SURFACE=Pad,ORIENTATION=Local", "ORIENTATION"),
        ("SURFACE=Pad", "SURFACE=Pad,INFLUENCE RADIUS=1", "INFLUENCERADIUS"),
        ("*KINEMATIC\n1,6", "*DISTRIBUTING,COUPLING=STRUCTURAL\n1,6", "STRUCTURAL"),
        ("*KINEMATIC\n1,6", "*KINEMATIC\n4,6", "at least one translation"),
        ("*KINEMATIC\n1,6", "*DISTRIBUTING\n1,3", "distributing requires 1..6"),
    ] {
        let text = constraint_input().replace(from, to);
        std::fs::write(&input, &text).unwrap();
        let mut s = Session::default();
        let read = s.execute("deck.read", &json!({"path":input,"format":"abaqus"})).unwrap();
        assert_eq!(read["analysis_complete"], false, "{read}");
        assert!(read["warnings"].to_string().contains(reason), "{reason}: {read}");
        assert!(s.model.steps().is_empty());
        let replay = dir.join(format!("{}.inp", reason.replace('/', "-")));
        s.execute("deck.source.write", &json!({"path":replay})).unwrap();
        assert_eq!(std::fs::read(replay).unwrap(), text.as_bytes());
    }
    let mut s = native_constraint_patch();
    for (from, to) in [("Secondary,Main", "Missing,Main"), ("REF NODE=RP", "REF NODE=100"), ("SURFACE=Pad", "SURFACE=Missing")] {
        std::fs::write(&input, constraint_input().replace(from, to)).unwrap();
        let before = s.clone();
        let err = s.execute("deck.read", &json!({"path":input,"format":"abaqus"})).unwrap_err();
        assert!(err.to_string().contains("missing"), "{err}");
        assert_unchanged(&s, &before);
    }
    std::fs::remove_dir_all(dir).unwrap();
}

fn assert_constraint_results(native: &Session, imported: &Session) {
    fn result(s: &Session) -> &fieldcraft_solver::StaticResult {
        let crate::StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = s.result("Pull").unwrap().as_ref() else {
            panic!("static result");
        };
        r
    }
    let (a, b) = (result(native), result(imported));
    for (x, y) in a.displacement.iter().flatten().zip(b.displacement.iter().flatten()) {
        assert!((x - y).abs() < 1e-12, "displacement {x} vs {y}");
    }
    for (x, y) in a.stress.iter().flatten().zip(b.stress.iter().flatten()) {
        assert!((x - y).abs() < 1e-9, "stress {x} vs {y}");
    }
    for (x, y) in a.reaction_total.iter().zip(b.reaction_total.iter()) {
        assert!((x - y).abs() < 1e-12, "reaction {x} vs {y}");
    }
    assert_eq!(a.displacement.len(), b.displacement.len());
    let du = a.displacement.iter().flatten().zip(b.displacement.iter().flatten()).map(|(x, y)| (x - y).abs()).fold(0.0_f64, f64::max);
    let ds = a.stress.iter().flatten().zip(b.stress.iter().flatten()).map(|(x, y)| (x - y).abs()).fold(0.0_f64, f64::max);
    println!("constraint native comparison max_du={du:.6e} max_ds={ds:.6e} applied={:?} reaction={:?}", b.applied_total, b.reaction_total);

    assert!((b.reaction_total[2] + b.applied_total[2]).abs() < 1e-12);
}
#[test]
fn imported_constraints_instance_references_and_moment_transfer() {
    let dir = folder();
    let input = dir.join("instance.inp");
    let material = "*MATERIAL,NAME=M\n*ELASTIC\n1000,0\n*DENSITY\n1\n";
    let text = constraint_input().replace(material, "");
    let (part, step) = text.split_once("*STEP").unwrap();
    let text = format!(
        "{material}*PART,NAME=P\n{part}*END PART\n*ASSEMBLY,NAME=A\n*INSTANCE,NAME=I,PART=P\n3,0,0\n*END INSTANCE\n*END ASSEMBLY\n*STEP{}",
        step.replace("Base,1,3", "I.Base,1,3").replace("RP,3,1", "I.RP,3,1\nI.RP,5,.2")
    );
    std::fs::write(&input, text).unwrap();
    let mut native = native_constraint_patch();
    native.execute("model.load.create", &json!({"name":"Moment","type":"moment","target":{"reference_point":"RP"},"vector":[0,0.2,0]})).unwrap();
    native
        .execute("model.step.edit", &json!({"name":"Pull","procedure":{"type":"static"},"boundary_conditions":["Clamp"],"loads":["F","Moment"]}))
        .unwrap();
    native.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let mut s = Session::default();
    let reply = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["analysis_complete"], true, "{reply}");
    let rp = s.model.reference_points().values().next().unwrap();
    assert_eq!(rp.xyz, [3.5, 0.5, 2.0]);
    assert_eq!(reply["instances"][0]["name"], "I");
    s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    assert_constraint_results(&native, &s);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn imported_constraints_uniform_distribution_reopens_and_round_trips_with_moment() {
    let dir = folder();
    let input = dir.join("distributed.inp");
    let text = constraint_input().replace("*KINEMATIC", "*DISTRIBUTING,WEIGHTING METHOD=UNIFORM").replace("RP,3,1", "RP,3,1\nRP,5,.2");
    std::fs::write(&input, &text).unwrap();
    let mut native = native_constraint_patch();
    native.execute("model.coupling.edit", &json!({"name":"Rigid","reference_point":"RP","set":"Top","kind":"distributing"})).unwrap();
    native.execute("model.load.create", &json!({"name":"Moment","type":"moment","target":{"reference_point":"RP"},"vector":[0,0.2,0]})).unwrap();
    native
        .execute("model.step.edit", &json!({"name":"Pull","procedure":{"type":"static"},"boundary_conditions":["Clamp"],"loads":["F","Moment"]}))
        .unwrap();
    native.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["analysis_complete"], true, "{read}");
    for phase in 0..3 {
        s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
        assert_constraint_results(&native, &s);
        if phase == 0 {
            let saved = dir.join("saved.fieldcraft");
            s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
            s.execute("file.new", &json!({})).unwrap();
            s.execute("file.open", &json!({"path":saved})).unwrap();
            let replay = dir.join("replay.inp");
            s.execute("deck.source.write", &json!({"path":replay})).unwrap();
            assert_eq!(std::fs::read(replay).unwrap(), text.as_bytes());
        } else if phase == 1 {
            let path = dir.join("normalized.inp");
            s.execute("deck.write", &json!({"path":path,"format":"abaqus"})).unwrap();
            let read = s.execute("deck.read", &json!({"path":path,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
            assert_eq!(read["analysis_complete"], true, "{read}");
        }
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn imported_constraints_split_mixed_node_point_loads_and_reference_boundaries() {
    let dir = folder();
    let input = dir.join("mixed.inp");
    let text = constraint_input()
        .replace("*STEP,NAME=Pull", "*NSET,NSET=Mixed\n13,99\n*STEP,NAME=Pull")
        .replace("Base,1,3", "Base,1,3\nRP,1,1")
        .replace("RP,3,1", "Mixed,3,1");
    std::fs::write(&input, text).unwrap();
    let mut native = native_constraint_patch();
    native.run_script(&json!([
        {"command":"model.set.create","params":{"name":"Corner","nodes":[13]}},
        {"command":"model.load.create","params":{"name":"CornerForce","type":"force","target":{"set":"Corner"},"vector":[0,0,1]}},
        {"command":"model.boundary_condition.create","params":{"name":"Hold","type":"displacement","target":{"reference_point":"RP"},"values":[0,null,null,null,null,null]}},
        {"command":"model.step.edit","params":{"name":"Pull","procedure":{"type":"static"},"boundary_conditions":["Clamp","Hold"],"loads":["F","CornerForce"]}}
    ])).unwrap();
    native.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let mut s = Session::default();
    let reply = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["analysis_complete"], true, "{reply}");
    assert!(s.model.boundary_conditions().values().any(|b| b.target().reference_point().is_some()));
    assert_eq!(s.model.loads().len(), 2);
    s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    assert_constraint_results(&native, &s);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn imported_mass_requirements_are_unapplied_and_survive_source_and_model_round_trips() {
    let dir = folder();
    let input = dir.join("input.inp");
    std::fs::write(&input, constraint_input()).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    let report = read["mass_semantics"].clone();
    assert_eq!(report["applied"], false);
    assert_eq!(report["conditions_verified"], false);
    assert_eq!(report["requirements"][0]["choice"], "lumped");
    assert_eq!(report["requirements"][0]["rule"], "row_sum");
    assert!(report["requirements"][0]["native_gap"].as_str().unwrap().contains("HRZ"));
    let saved = dir.join("model.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    let again = s.execute("deck.read", &json!({"path":replay,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(again["mass_semantics"], report);
    let normalized = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"path":normalized,"format":"abaqus"})).unwrap();
    let again = s.execute("deck.read", &json!({"path":normalized,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    // The native writer emits ELSET separately from ELEMENT. Source-card metadata changes,
    // while type, rule, choice, applicability and the explicit unapplied status must not.
    let mut normalized_report = report;
    normalized_report["requirements"][0]["elset"] = Value::Null;
    assert_eq!(again["mass_semantics"], normalized_report);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn imported_mass_requirement_row_sum_is_not_native_hrz_on_a_tapered_hex() {
    // Original analytic map x=(1+z)*xi, y=eta, z=zeta on [0,1]^3.
    // detJ=1+z; volume=3/2. Integrating each Ni gives bottom 1/6, top 5/24.
    // Integrating Ni^2 gives bottom 5/108, top 7/108. HRZ scales by 27/8.
    let xyz = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [2., 0., 1.], [2., 1., 1.], [0., 1., 1.]];
    let solid = fieldcraft_fem::element::Solid::new(fieldcraft_model::ElementKind::Hex8, &xyz, fieldcraft_fem::element::Formulation::Full).unwrap();
    let c = solid.mass(1.0).unwrap();
    let hrz = solid.lumped_mass(1.0).unwrap();
    for (i, h) in hrz.iter().enumerate() {
        let row: f64 = (0..8).map(|j| c.get(i, j)).sum();
        let expected = if i < 4 { 1.0 / 6.0 } else { 5.0 / 24.0 };
        let expected_hrz = if i < 4 { 5.0 / 32.0 } else { 7.0 / 32.0 };
        assert!((row - expected).abs() < 1e-14);
        assert!((h - expected_hrz).abs() < 1e-14);
        assert!((row - h).abs() > 0.01, "lumped must not be reported as row-sum");
    }
    println!("analytic tapered hex: total mass=1.5 kg; row_sum=[1/6 x4,5/24 x4]; native_HRZ=[5/32 x4,7/32 x4]; max_difference=1/96 kg");
}

#[test]
fn section_mass_policy_is_effective_after_native_reopen_and_deck_round_trip() {
    let dir = folder();
    let input = dir.join("mass.inp");
    let text = "*NODE\n1,0,0,0\n2,1,0,0\n3,1,1,0\n4,0,1,0\n5,0,0,1\n6,2,0,1\n7,2,1,1\n8,0,1,1\n*ELEMENT,TYPE=C3D8,ELSET=Body\n1,1,2,3,4,5,6,7,8\n*NSET,NSET=Base\n1,2,3,4\n*MATERIAL,NAME=M\n*ELASTIC\n1000,0\n*DENSITY\n1\n*SOLID SECTION,ELSET=Body,MATERIAL=M\n** FIELDCRAFT LUMPED MASS, ELSET=Body, STRATEGY=ROW_SUM\n*STEP,NAME=Modes\n*FREQUENCY\n1\n*BOUNDARY\nBase,1,3\n*END STEP\n";
    std::fs::write(&input, text).unwrap();
    let mut s = Session::default();
    let read = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(read["native_lumped_mass_policies"][0]["strategy"], "row_sum");
    let a = s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap();
    let b = s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"})).unwrap();
    assert_eq!(a, b, "default remains consistent");
    let row = s.execute("solve.frequency", &json!({"step":"Modes","mass":"lumped"})).unwrap();
    let expected = row["modes"][0]["frequency"].as_f64().unwrap();
    assert_eq!(row["mass"], "lumped");
    let saved = dir.join("mass.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    assert!(s.result("Modes").is_some(), "frequency result and unchanged format remain readable");
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), text);
    for phase in 0..2 {
        let reply = s.execute("solve.frequency", &json!({"step":"Modes","mass":"lumped"})).unwrap();
        let actual = reply["modes"][0]["frequency"].as_f64().unwrap();
        assert!((actual - expected).abs() < 1e-12);
        println!("native row_sum phase={phase} frequency={actual:.16e} Hz");
        if phase == 0 {
            let path = dir.join("normalized.inp");
            s.execute("deck.write", &json!({"path":path,"format":"abaqus"})).unwrap();
            assert!(std::fs::read_to_string(&path).unwrap().contains("FIELDCRAFT LUMPED MASS"));
            s.execute("deck.read", &json!({"path":path,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
        }
    }
    let (name, section) = s.model.sections().iter().next().unwrap();
    let mut edit = serde_json::to_value(section).unwrap();
    edit["name"] = json!(name);
    edit["formulation"]["lumped_mass"] = json!("hrz");
    s.execute("model.section.edit", &edit).unwrap();
    let hrz = s.execute("solve.frequency", &json!({"step":"Modes","mass":"lumped"})).unwrap();
    let hrz = hrz["modes"][0]["frequency"].as_f64().unwrap();
    assert!((hrz - expected).abs() > expected * 0.01, "the selected mass changes the actual eigenproblem");
    assert!((expected / hrz - (21.0_f64 / 20.0).sqrt()).abs() < 1e-12, "the free top-node masses differ by exactly 21/20");
    println!("mass-only tapered hex: row_sum={expected:.16e} Hz HRZ={hrz:.16e} Hz");
    for bad in [text.replace("STRATEGY=ROW_SUM", "STRATEGY=GUESS"), text.replace("ELSET=Body, STRATEGY", "ELSET=Missing, STRATEGY")] {
        std::fs::write(&input, bad).unwrap();
        let before = s.clone();
        assert!(s.execute("deck.read", &json!({"path":input,"format":"abaqus"})).is_err());
        assert_unchanged(&s, &before);
    }

    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn bbar_two_bricks_reopen_and_round_trip_preserve_response_and_identity() {
    let dir = folder();
    let input = dir.join("bbar.inp");
    let text = r#"*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0,0,1
6,1,0,1
7,1,1,1
8,0,1,1
9,0,0,1
10,1,0,1
11,1,1,1
12,0,1,1
13,0,0,2
14,1,0,2
15,1,1,2
16,0,1,2
*ELEMENT,TYPE=C3D8,ELSET=All
1,1,2,3,4,5,6,7,8
2,9,10,11,12,13,14,15,16
*NSET,NSET=Base
1,2,3,4
*NSET,NSET=Top
13,14,15,16
*SURFACE,TYPE=ELEMENT,NAME=Main
1,S2
*SURFACE,TYPE=ELEMENT,NAME=Secondary
2,S1
*MATERIAL,NAME=M
*ELASTIC
1000,0
*DENSITY
1
*SOLID SECTION,ELSET=All,MATERIAL=M
*TIE,NAME=Glue,ADJUST=NO,POSITION TOLERANCE=.01
Secondary,Main
*STEP,NAME=Pull
*STATIC
*BOUNDARY
Base,1,3
*CLOAD
13,3,.35
14,3,.15
15,3,.15
16,3,.35
*END STEP
"#;
    std::fs::write(&input, text).unwrap();
    let mut s = Session::default();
    s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    let mut reference = None;
    for phase in 0..3 {
        assert!(s.model.sections().values().all(|sec| sec.formulation.b_bar));
        s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
        if let Some(ref baseline) = reference {
            assert_constraint_results(baseline, &s);
        } else {
            reference = Some(s.clone());
        }
        let crate::StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = s.result("Pull").unwrap().as_ref() else { panic!("static") };
        let mut mean = [0.; 3];
        let mut theta = 0.;
        for id in 13..=16 {
            let i = r.nodes.binary_search(&id).unwrap();
            let u = r.displacement[i];
            for j in 0..3 {
                mean[j] += u[j] / 4.;
            }
            let x = if id == 13 || id == 16 { -0.5 } else { 0.5 };
            theta -= x * u[2];
        }
        println!("B-bar two bricks phase={phase}: mean={mean:?}, theta_y={theta:.16e}");
        if phase == 0 {
            let path = dir.join("bbar.fieldcraft");
            s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
            s.execute("file.new", &json!({})).unwrap();
            s.execute("file.open", &json!({"path":path})).unwrap();
            let replay = dir.join("source.inp");
            s.execute("deck.source.write", &json!({"path":replay})).unwrap();
            assert_eq!(std::fs::read_to_string(replay).unwrap(), text);
        }
        if phase == 1 {
            let path = dir.join("normal.inp");
            s.execute("deck.write", &json!({"path":path,"format":"abaqus"})).unwrap();
            s.execute("deck.read", &json!({"path":path,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
        }
    }
    // Native standard stays standard and is refused by the Abaqus writer, rather than relabelled.
    let (name, sec) = s.model.sections().iter().next().unwrap();
    let mut edit = serde_json::to_value(sec).unwrap();
    edit["name"] = json!(name);
    edit["formulation"]["b_bar"] = json!(false);
    s.execute("model.section.edit", &edit).unwrap();
    s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let error = s.execute("deck.write", &json!({"path":dir.join("standard.inp"),"format":"abaqus"})).unwrap_err().to_string();
    assert!(error.contains("standard full Hex8"), "{error}");
    for incompatible in [false, true] {
        let (name, sec) = s.model.sections().iter().next().unwrap();
        let mut section = sec.clone();
        let name = name.clone();
        section.formulation.incompatible_modes = incompatible;
        s.model.replace_section(&name, section).unwrap();
        s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
        let crate::StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = s.result("Pull").unwrap().as_ref() else { panic!("static") };
        let mut mean = [0.; 3];
        let mut theta = 0.;
        for id in 13..=16 {
            let u = r.displacement[r.nodes.binary_search(&id).unwrap()];
            for j in 0..3 {
                mean[j] += u[j] / 4.;
            }
            let x = if id == 13 || id == 16 { -0.5 } else { 0.5 };
            theta -= x * u[2];
        }
        println!(
            "Native two bricks incompatible={incompatible}: mean={mean:?}, theta_y={theta:.16e}; beam theory ux=0.0048, uz=0.002, theta_y=0.0048 (limited rectangular nu=0 case)"
        );
    }

    let before = s.clone();
    std::fs::write(&input, text.replace("TYPE=C3D8,", "TYPE=C3D8I,")).unwrap();
    let error = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap_err().to_string();
    assert!(error.contains("thirteen internal modes"), "{error}");
    assert_unchanged(&s, &before);
    let replay = dir.join("retained-after-refusal.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert!(std::fs::read_to_string(replay).unwrap().contains("TYPE=C3D8"));
    std::fs::remove_dir_all(dir).unwrap();
}

fn property_closure(text: &str, native: &mut Session, expected: f64) {
    let dir = folder();
    let input = dir.join("input.inp");
    std::fs::write(&input, text).unwrap();
    native.execute("solve.static", &json!({"step":"Pull"})).unwrap();
    let mut s = Session::default();
    let reply = s.execute("deck.read", &json!({"path":input,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["analysis_complete"], true, "{reply}");
    for phase in 0..3 {
        let reply = s.execute("solve.static", &json!({"step":"Pull"})).unwrap();
        let u = reply["max_displacement"]["value"].as_f64().unwrap();
        assert!((u - expected).abs() < 1e-12, "phase={phase}, u={u}, expected={expected}");
        assert_constraint_results(native, &s);
        if phase == 0 {
            let saved = dir.join("saved.fieldcraft");
            s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
            s.execute("file.new", &json!({})).unwrap();
            s.execute("file.open", &json!({"path":saved})).unwrap();
            let replay = dir.join("source.inp");
            s.execute("deck.source.write", &json!({"path":replay})).unwrap();
            assert_eq!(std::fs::read(replay).unwrap(), text.as_bytes());
        } else if phase == 1 {
            let normalized = dir.join("normalized.inp");
            s.execute("deck.write", &json!({"path":normalized,"format":"abaqus"})).unwrap();
            let read = s.execute("deck.read", &json!({"path":normalized,"format":"abaqus","units":"m-kg-s-N-Pa"})).unwrap();
            assert_eq!(read["analysis_complete"], true, "{read}");
        }
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn orthotropic_stiffness_rotated_solid_solves_reopens_and_normalizes() {
    let text = "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n*ELEMENT,TYPE=C3D4,ELSET=All\n1,1,2,3,4\n*NSET,NSET=Base\n1,3,4\n*NSET,NSET=Top\n2\n*MATERIAL,NAME=M\n*ELASTIC,TYPE=ORTHOTROPIC\n1000,0,2000,0,0,3000,400,500\n600\n*ORIENTATION,NAME=Axes\n1,0,0,0,1,0\n3,90\n*SOLID SECTION,ELSET=All,MATERIAL=M,ORIENTATION=Axes\n*STEP,NAME=Pull\n*STATIC\n*BOUNDARY\nBase,1,3\n*CLOAD\nTop,1,1\n*END STEP\n";
    let mut native = Session::default();
    native.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    for (i, xyz) in [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]].iter().enumerate() {
        native.execute("mesh.node.add", &json!({"label":i+1,"xyz":xyz})).unwrap();
    }
    native.run_script(&json!([
        {"command":"mesh.element.add","params":{"label":1,"kind":"tet4","nodes":[1,2,3,4]}},
        {"command":"model.set.create","params":{"name":"All","elements":[1]}},
        {"command":"model.set.create","params":{"name":"Base","nodes":[1,3,4]}},
        {"command":"model.set.create","params":{"name":"Top","nodes":[2]}},
        {"command":"model.material.create","params":{"name":"M","elastic":{"type":"engineering_constants","young":[1000,2000,3000],"poisson":[0,0,0],"shear":[400,500,600]}}},
        {"command":"model.section.create","params":{"name":"S","set":"All","material":"M","kind":{"type":"solid"},"orientation":{"frame":{"type":"rectangular","axis1":[0,1,0],"axis2":[-1,0,0]}}}},
        {"command":"model.boundary_condition.create","params":{"name":"Clamp","type":"displacement","target":{"set":"Base"},"values":[0,0,0,null,null,null]}},
        {"command":"model.load.create","params":{"name":"F","type":"force","target":{"set":"Top"},"vector":[1,0,0]}},
        {"command":"model.step.create","params":{"name":"Pull","procedure":{"type":"static"},"boundary_conditions":["Clamp"],"loads":["F"]}}
    ])).unwrap();
    // Unit tetra V=1/6, grad(N2)=ex: Kxx=E_local2/6 after 90 degree rotation.
    property_closure(text, &mut native, 6. / 2000.);
}

#[test]
fn partial_kinematic_translation_with_explicit_reference_restraints_closes() {
    let text = constraint_input().replace("*KINEMATIC\n1,6", "*KINEMATIC\n3,3").replace("Base,1,3\n*CLOAD", "Base,1,3\nRP,1,2\nRP,4,6\n*CLOAD");
    let mut native = native_constraint_patch();
    let mut c = native.model.coupling("Rigid").unwrap().clone();
    c.dofs = fieldcraft_model::Dofs::from_bits(4);
    native.model.replace_coupling("Rigid", c).unwrap();
    native
        .execute(
            "model.boundary_condition.create",
            &json!({"name":"RPClamp","type":"displacement","target":{"reference_point":"RP"},"values":[0,0,null,0,0,0]}),
        )
        .unwrap();
    let mut step = native.model.step("Pull").unwrap().clone();
    step.boundary_conditions.insert("RPClamp".into());
    native.model.replace_step(step).unwrap();
    property_closure(&text, &mut native, 0.002);
}

#[test]
fn kinematic_element_surface_resolves_face_nodes_and_roundtrips() {
    let text = constraint_input().replace("*SURFACE,TYPE=NODE,NAME=Pad\nTop,1", "*SURFACE,TYPE=ELEMENT,NAME=Pad\n2,S2");
    property_closure(&text, &mut native_constraint_patch(), 0.002);
}

#[test]
fn random_continuum_mass_uses_requested_operator_not_cached_modes_and_reopens() {
    // One free tetrahedron DOF: k=E*V=1000/6; C22=rho*V/10=.2,
    // row-sum L22=rho*V/4=.5. These are independent simplex integrals.
    let source = "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n*ELEMENT,TYPE=C3D4,ELSET=Body\n1,1,2,3,4\n*MATERIAL,NAME=M\n*ELASTIC\n1000,0\n*DENSITY\n12\n*SOLID SECTION,ELSET=Body,MATERIAL=M\n** FIELDCRAFT LUMPED MASS,ELSET=Body,STRATEGY=ROW_SUM\n*BOUNDARY\n1,1,3\n3,1,3\n4,1,3\n2,2,3\n*PSD-DEFINITION,NAME=White,TYPE=FORCE\n1,0,0\n1,0,10\n*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY\n1\n*END STEP\n*STEP,NAME=Random,PERTURBATION\n*RANDOM RESPONSE\n0,10,11,3,1\n*MODAL DAMPING,DEFINITION=MODE NUMBERS\n1,1,.05\n*CLOAD,LOAD CASE=2\n2,1,1\n*CORRELATION,TYPE=CORRELATED,PSD=White\n2,1\n*END STEP\n";
    let dir = folder();
    let input = dir.join("input.inp");
    std::fs::write(&input, source).unwrap();
    let mut s = Session::default();
    let imported = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(imported["analysis_complete"], true, "{imported}");
    assert!(s.model.step("Random").is_some(), "{:?}", s.model.steps());
    let grid = json!([0, 1, 2, 3, 4, 5, 6, 8, 10]);
    let mut lumped = Value::Null;
    for (choice, mass) in [("consistent", 0.2_f64), ("lumped", 0.5)] {
        // Deliberately seed an incompatible cached frequency result before the lumped request.
        s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"})).unwrap();
        let out = s.execute("solve.random", &json!({"step":"Random","mass":choice,"nodes":[2],"sampling_frequencies":grid})).unwrap();
        let k = 1000. / 6.;
        let f = (k / mass).sqrt() / (2. * std::f64::consts::PI);
        assert!((out["natural_frequencies"][0].as_f64().unwrap() - f).abs() < 1e-10);
        assert_eq!(out["mass"], choice);
        assert_eq!(out["sampling_source"], "explicit_replay_grid");
        let samples = out["probes"][0]["samples"].as_array().unwrap();
        for sample in samples {
            let omega = 2. * std::f64::consts::PI * sample["frequency"].as_f64().unwrap();
            let c = 0.1 * (k * mass).sqrt();
            let expected = 1. / ((k - omega * omega * mass).powi(2) + (omega * c).powi(2));
            let actual = sample["displacement_components_psd"][0].as_f64().unwrap();
            assert!((actual - expected).abs() < expected * 1e-8);
        }
        if choice == "lumped" {
            lumped = out;
        }
    }
    let saved = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), source);
    let normalized = dir.join("normalized.inp");
    let error = s.execute("deck.write", &json!({"path":normalized,"format":"abaqus"})).unwrap_err().to_string();
    assert!(error.contains("random response is not written yet"), "{error}");
    let again = s.execute("solve.random", &json!({"step":"Random","mass":"lumped","nodes":[2],"sampling_frequencies":grid})).unwrap();
    assert_eq!(again["natural_frequencies"], lumped["natural_frequencies"]);
    assert_eq!(again["node_displacement_rms"], lumped["node_displacement_rms"]);
    for bad in [json!([0, 11]), json!([0, 0]), json!([1])] {
        assert!(s.execute("solve.random", &json!({"step":"Random","mass":"lumped","sampling_frequencies":bad})).is_err());
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn enhanced13_native_static_persists_and_refuses_unvalidated_deck_mapping() {
    let dir = folder();
    let mut s = bar();
    let section = s.model.sections().get("BarSection").unwrap();
    let mut edit = serde_json::to_value(section).unwrap();
    edit["formulation"] = json!({"enhanced_strain13":true});
    edit["name"] = json!("BarSection");
    s.execute("model.section.edit", &edit).unwrap();
    let result = s.execute("solve.static", &json!({"step":"Bend"}));
    assert!(result.is_ok(), "{result:?}");
    let probe = s.execute("solve.probe", &json!({"step":"Bend","set":"end"})).unwrap();
    s.execute("file.save_as", &json!({"path":dir.join("enhanced.fieldcraft"),"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":dir.join("enhanced.fieldcraft")})).unwrap();
    assert!(s.model.sections().get("BarSection").unwrap().formulation.enhanced_strain13);
    assert_eq!(s.execute("solve.probe", &json!({"step":"Bend","set":"end"})).unwrap(), probe);
    let error = s.execute("deck.write", &json!({"path":dir.join("not-equivalent.inp"),"format":"abaqus"})).unwrap_err().to_string();
    assert!(error.contains("enhanced_strain13"), "{error}");
    s.execute("model.step.create", &json!({"name":"Modes","procedure":{"type":"frequency","modes":2}})).unwrap();
    let error = s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap_err().to_string();
    assert!(error.contains("enhanced_strain13"), "{error}");
    std::fs::remove_dir_all(dir).unwrap();
}
