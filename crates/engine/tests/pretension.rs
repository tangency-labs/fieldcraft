//! End-to-end independent cut load, contact equilibrium, lock and linear perturbation.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

#[test]
fn bolt_clamp_load_lock_persistence_and_perturbation() {
    let recipe: Value = serde_json::from_str(include_str!("fixtures/pretension.json")).unwrap();
    let mut session = Session::default();
    let results = session.run_script(&recipe).unwrap();
    let outputs: Vec<_> = recipe.as_array().unwrap().iter().zip(&results).collect();
    let pre = outputs.iter().find(|(c, _)| c["command"] == "solve.pretension" && c["params"]["step"] == "Preload").unwrap().1;
    let locked = outputs.iter().find(|(c, _)| c["command"] == "solve.pretension" && c["params"]["step"] == "Locked").unwrap().1;
    let q = pre["adjustment"].as_f64().unwrap();
    // Two equal axial compliances (bolt and members), each L/(EA)=2/1e6.
    let base = outputs.iter().find(|(c, _)| c["command"] == "solve.nonlinear" && c["params"]["step"] == "Preload").unwrap().1;
    let penetration = base["contacts"][0]["penetration"].as_f64().unwrap();
    // The contact algorithm permits gap <= 1e-5 times element thickness.
    // Cut closure includes this measured gap as well as the two elastic strains.
    assert!((0. ..=1e-5).contains(&penetration));
    assert!(((q - penetration) / 0.0004 - 1.).abs() < 0.005, "{pre}; penetration={penetration}");
    assert!((locked["adjustment"].as_f64().unwrap() - q).abs() < 1e-12);
    assert!((locked["section_force"].as_f64().unwrap() / 100. - 1.).abs() < 0.005);
    for (c, r) in &outputs {
        if c["command"] == "solve.nonlinear" {
            let unloaded = format!("step `{}`: nothing loads", c["params"]["step"].as_str().unwrap());
            assert!(r["warnings"].as_array().unwrap().iter().all(|w| !w.as_str().unwrap().contains(&unloaded)));
            let f = r["contacts"][0]["normal_force"].as_f64().unwrap();
            assert!((f / 100. - 1.).abs() < 0.005, "clamped member force: {f}");
        }
        if c["command"] == "solve.frequency" {
            assert_eq!(r["modes"].as_array().unwrap().len(), 6);
        }
    }
    check_frozen_perturbations(&session);
    for step in ["Preload", "Locked"] {
        let fieldcraft_engine::StepResult::Implicit(fieldcraft_solver::StepResult::Static(result)) = session.result(step).unwrap().as_ref() else {
            panic!("equilibrium")
        };
        let bytes = fieldcraft_io::codec::encode(result).unwrap();
        let restored: fieldcraft_solver::StaticResult = fieldcraft_io::codec::decode(&bytes).unwrap();
        assert_eq!(*result, restored);
        assert_eq!(restored.pretension.len(), 1);
    }
    let stored = serde_json::to_vec(&session.model).unwrap();
    let loaded: fieldcraft_model::Model = serde_json::from_slice(&stored).unwrap();
    assert_eq!(session.model, loaded);
    // The lock journals a numerical, unit-bearing adjustment, so replay needs no results.
    let journal: Value = session.journal.iter().map(|(command, params)| json!({"command":command,"params":params})).collect();
    assert!(session.journal.iter().any(|(c, _)| c == "model.pretension.adjust"));
    assert!(session.journal.iter().all(|(c, _)| c != "model.pretension.lock"));
    let mut replay = Session::default();
    replay.run_script(&journal).unwrap();
    assert_eq!(session.model, replay.model);
    let err = fieldcraft_decks::abaqus::write(&session.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(err.contains("disconnected topology"), "{err}");
}

// Check physical increments, not just that a solver returned successfully.
fn check_frozen_perturbations(session: &Session) {
    use fieldcraft_engine::StepResult;
    use fieldcraft_solver::StepResult as Implicit;
    let StepResult::Implicit(Implicit::Frequency(modes)) = session.result("Modes").unwrap().as_ref() else { panic!("modes") };
    let StepResult::Implicit(Implicit::Harmonic(response)) = session.result("Response").unwrap().as_ref() else { panic!("response") };
    let pairs = [(5, 9), (6, 10), (7, 11), (8, 12), (21, 25), (22, 26), (23, 27), (24, 28)];
    let check = |labels: &[u32], values: &[[f64; 3]]| {
        let scale = values.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs())).max(1e-30);
        for (a, b) in pairs {
            let a = labels.iter().position(|n| *n == a).unwrap();
            let b = labels.iter().position(|n| *n == b).unwrap();
            for (a, b) in values[a].iter().zip(&values[b]) {
                assert!((a - b).abs() <= 1e-8 * scale, "frozen interface increment");
            }
        }
    };
    let control = modes.reference_points.iter().position(|(n, _)| n == "CutControl").unwrap();
    for mode in &modes.modes {
        check(&modes.nodes, &mode.shape);
        assert_eq!(mode.point_shape[control][0], 0.0);
    }
    assert_eq!(response.frames.len(), 4);
    assert_eq!(response.stats.equations, modes.stats.equations);
    let control = response.reference_points.iter().position(|(n, _)| n == "CutControl").unwrap();
    for frame in &response.frames {
        assert!(frame.relative_residual < 1e-8, "harmonic residual: {}", frame.relative_residual);
        check(&response.nodes, &frame.real);
        check(&response.nodes, &frame.imaginary);
        assert_eq!(frame.points_real[control][0], 0.0);
        assert_eq!(frame.points_imaginary[control][0], 0.0);
        let scale = frame.input_power.abs().max(frame.dissipated_power.abs()).max(1e-30);
        assert!((frame.input_power - frame.dissipated_power).abs() < 1e-7 * scale);
    }
}

#[test]
fn tied_clamp_keeps_its_mapping_through_lock_and_perturbations() {
    let mut recipe: Value = serde_json::from_str(include_str!("fixtures/pretension.json")).unwrap();
    let contact = recipe.as_array_mut().unwrap().iter_mut().find(|c| c["command"] == "interaction.contact.pair").unwrap();
    *contact = json!({"command":"interaction.tie", "params":{"name":"ClampTie", "main":"member_top", "secondary":"member_bottom"}});
    let mut session = Session::default();
    session.run_script(&recipe).unwrap();
    check_frozen_perturbations(&session);
}

#[test]
fn abaqus_native_chain_rejects_missing_lock_and_late_unbased_frequency() {
    let recipe: Value = serde_json::from_str(include_str!("fixtures/pretension.json")).unwrap();
    let mut session = Session::default();
    session.run_script(&recipe).unwrap();
    let mut locked = session.model.step("Locked").cloned().unwrap();
    locked.pretension.clear();
    session.model.replace_step(locked).unwrap();
    let err = fieldcraft_decks::abaqus::write(&session.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(err.contains("missing nonzero adjustment"), "{err}");

    let mut session = Session::default();
    session.run_script(&recipe).unwrap();
    session
        .run_script(&serde_json::json!([{
            "command":"model.step.create",
            "params":{"name":"LateFreeModes","procedure":{"type":"frequency","modes":2},"boundary_conditions":["Bottom","CapGuide"]}
        }]))
        .unwrap();
    let err = fieldcraft_decks::abaqus::write(&session.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(err.contains("unbased frequency/response follows pretension"), "{err}");
}
