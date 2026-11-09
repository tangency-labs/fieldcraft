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
        ("COUPLING", "*COUPLING,CONSTRAINT NAME=Joint,REF NODE=1,SURFACE=Slave\n*KINEMATIC\n1,6\n"),
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
