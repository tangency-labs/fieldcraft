//! End-to-end independent cut load, contact equilibrium, lock and linear perturbation.
use fieldcraft_engine::Session;
use fieldcraft_model::{
    BoundaryCondition, Elastic, ElementFace, ElementKind, Formulation, Material, Model, Procedure, ReferencePoint, Section, SectionKind, Set, Step,
    Target,
};
use serde_json::{Value, json};

#[test]
fn native_fixed_releases_the_cut_force_and_retains_ordinary_loads() {
    let recipe: Vec<Value> = serde_json::from_str(include_str!("fixtures/internal-pretension.json")).unwrap();
    let mut session = Session::default();
    session.run_script(&Value::Array(recipe[..recipe.len() - 2].to_vec())).unwrap();
    let mut model = session.model.clone();
    model.add_set("OneNode", Set { nodes: vec![9], ..Set::default() }).unwrap();
    model
        .add_load(
            "Ordinary",
            fieldcraft_model::Load::Force { target: Target::Set("OneNode".into()), vector: [0., 0., 25.], distribution: Default::default() },
        )
        .unwrap();
    let mut preload = model.step("Preload").unwrap().clone();
    preload.loads.insert("Ordinary".into());
    model.replace_step(preload.clone()).unwrap();
    let mut locked = preload;
    locked.name = "Locked".into();
    locked.pretension.insert("Cut".into(), fieldcraft_model::PretensionControl::Fixed);
    model.add_step(locked, None).unwrap();
    let deck = fieldcraft_decks::abaqus::write(&model, fieldcraft_model::units::Preset::MKgS).unwrap();
    let locked = deck.split("** Step Locked\n").nth(1).unwrap();
    assert!(locked.contains("*CLOAD, OP=NEW\n9, 3, 25.0\n"), "{locked}");
    assert!(locked.contains("*BOUNDARY, FIXED\n100, 1, 1\n"));
    assert!(!locked.contains("100, 1, 1000"));
    assert!(!deck.contains("*EQUATION"));
    assert_eq!(&model.element(1).unwrap().nodes[4..], &model.element(2).unwrap().nodes[..4]);
}

#[test]
fn continuous_bolt_fixed_chain_keeps_terminal_contact_and_source_topology() {
    let recipe: Vec<Value> = serde_json::from_str(include_str!("fixtures/pretension.json")).unwrap();
    let label = |n: u64| if (9..=12).contains(&n) { n - 4 } else { n };
    let mut commands = Vec::new();
    for mut c in recipe {
        let command = c["command"].as_str().unwrap().to_owned();
        if command == "mesh.node.add" && (9..=12).contains(&c["params"]["label"].as_u64().unwrap()) {
            continue;
        }
        if matches!(command.as_str(), "mesh.element.add" | "model.set.create")
            && let Some(nodes) = c.get_mut("params").and_then(|p| p.get_mut("nodes")).and_then(Value::as_array_mut)
        {
            for n in nodes {
                *n = json!(label(n.as_u64().unwrap()));
            }
        }
        if command == "model.pretension.create" {
            commands.push(json!({"command":"model.set.create","params":{"name":"BoltFace","faces":[{"element":1,"face":2}]}}));
            c["params"]["face"] = json!("BoltFace");
        }
        if command == "model.pretension.lock" {
            c["command"] = json!("model.pretension.fix");
            c["params"].as_object_mut().unwrap().remove("from");
        }
        // This new current-value chain checks frictionless contact; the original
        // frictional, explicit-adjustment recipe remains an independent regression.
        if command == "interaction.contact.pair" {
            c["params"]["friction"] = json!({"static":0.0});
        }
        commands.push(c);
    }
    let mut session = Session::default();
    session.run_script(&Value::Array(commands)).unwrap();
    assert_eq!(&session.model.element(1).unwrap().nodes[4..], &session.model.element(2).unwrap().nodes[..4]);
    assert!(matches!(session.model.step("Locked").unwrap().pretension["Bolt"], fieldcraft_model::PretensionControl::Fixed));
    use fieldcraft_engine::StepResult;
    use fieldcraft_solver::StepResult as Implicit;
    let StepResult::Implicit(Implicit::Static(pre)) = session.result("Preload").unwrap().as_ref() else { panic!("preload") };
    let StepResult::Implicit(Implicit::Static(locked)) = session.result("Locked").unwrap().as_ref() else { panic!("locked") };
    assert_eq!(pre.pretension[0].adjustment.to_bits(), locked.pretension[0].adjustment.to_bits());
    assert!((locked.pretension[0].section_force / 100. - 1.).abs() < 0.005);
    for result in [pre, locked] {
        let history = result.nonlinear.as_ref().unwrap();
        assert!(history.increments.last().unwrap().converged);
        assert_eq!(history.increments.last().unwrap().load_factor, 1.);
        assert!((result.contacts[0].normal_force / 100. - 1.).abs() < 0.005);
    }
    let StepResult::Implicit(Implicit::Frequency(modes)) = session.result("Modes").unwrap().as_ref() else { panic!("modes") };
    let StepResult::Implicit(Implicit::Harmonic(response)) = session.result("Response").unwrap().as_ref() else { panic!("response") };
    assert_eq!(modes.nodes.len(), 28);
    assert_eq!(modes.nodes, response.nodes);
    let check = |values: &[[f64; 3]]| {
        let scale = values.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs())).max(1e-30);
        let mut mismatch = 0.0_f64;
        for (a, b) in [(21, 25), (22, 26), (23, 27), (24, 28)] {
            let a = modes.nodes.iter().position(|n| *n == a).unwrap();
            let b = modes.nodes.iter().position(|n| *n == b).unwrap();
            for (a, b) in values[a].iter().zip(values[b]) {
                mismatch = mismatch.max((a - b).abs() / scale);
            }
        }
        // Same terminal-mortar trace tolerance as the existing contact regression
        // below: the deformed projection need not join equal reference labels exactly.
        assert!(mismatch <= 1e-8, "terminal contact trace mismatch {mismatch}");
        println!("terminal contact trace mismatch={mismatch:.17e}");
    };
    let control = modes.reference_points.iter().position(|(name, _)| name == "CutControl").unwrap();
    for mode in &modes.modes {
        check(&mode.shape);
        assert_eq!(mode.point_shape[control][0], 0.);
    }
    for frame in &response.frames {
        check(&frame.real);
        check(&frame.imaginary);
        assert_eq!(frame.points_real[control][0], 0.);
        assert_eq!(frame.points_imaginary[control][0], 0.);
        assert!(frame.relative_residual < 1e-8);
    }
    fn roundtrip<T: serde::Serialize + serde::de::DeserializeOwned + PartialEq + std::fmt::Debug>(result: &T) {
        let bytes = fieldcraft_io::codec::encode(result).unwrap();
        let restored: T = fieldcraft_io::codec::decode(&bytes).unwrap();
        assert_eq!(*result, restored);
    }
    roundtrip(pre);
    roundtrip(locked);
    roundtrip(modes);
    roundtrip(response);
    session.execute("model.pretension.force", &json!({"name":"Bolt","step":"Preload","value":101})).unwrap();
    for name in ["Preload", "Locked", "Modes", "Response"] {
        assert!(session.result(name).is_none(), "stale {name}");
    }
}

#[test]
fn continuous_internal_cut_runs_through_commands_and_keeps_native_source_connected() {
    let recipe: Value = serde_json::from_str(include_str!("fixtures/internal-pretension.json")).unwrap();
    let mut session = Session::default();
    let commands = recipe.as_array().unwrap();
    session.run_script(&Value::Array(commands[..commands.len() - 2].to_vec())).unwrap();
    let original = session.model.clone();
    let before = fieldcraft_decks::abaqus::write(&original, fieldcraft_model::units::Preset::MKgS).unwrap();
    let results = session.run_script(&Value::Array(commands[commands.len() - 2..].to_vec())).unwrap();
    assert!((results[1]["adjustment"].as_f64().unwrap() / 0.002 - 1.).abs() < 1e-12);
    assert!((results[1]["section_force"].as_f64().unwrap() / 1000. - 1.).abs() < 1e-12);
    assert_eq!(session.model, original);
    assert_eq!(session.model.nodes().len(), 12);
    assert_eq!(&session.model.element(1).unwrap().nodes[4..], &session.model.element(2).unwrap().nodes[..4]);
    assert_eq!(before, fieldcraft_decks::abaqus::write(&session.model, fieldcraft_model::units::Preset::MKgS).unwrap());
    assert!(!before.contains("*EQUATION"));
    assert!(before.contains("*PRE-TENSION SECTION"));
    let stored = serde_json::to_vec(&session.model).unwrap();
    assert_eq!(session.model, serde_json::from_slice::<Model>(&stored).unwrap());
}

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
    // This nonlinear pretension chain has no verified B-bar nonlinear mapping;
    // standard full Hex8 is rejected explicitly and the solved physics stays intact.
    let standard_error = fieldcraft_decks::abaqus::write(&session.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(standard_error.contains("standard full Hex8"), "{standard_error}");
    let mut fixed = session.model.clone();
    let mut locked_step = fixed.step("Locked").cloned().unwrap();
    locked_step.pretension.insert("Bolt".into(), fieldcraft_model::PretensionControl::Fixed);
    fixed.replace_step(locked_step).unwrap();
    let fixed_error = fieldcraft_decks::abaqus::write(&fixed, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(fixed_error.contains("standard full Hex8"), "{fixed_error}");
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
    let error = fieldcraft_decks::abaqus::write(&session.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(error.contains("standard full Hex8"), "{error}");
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

fn tet10_writer_fixture(disconnected: bool) -> Model {
    let mut m = Model::default();
    // Two positively oriented tetrahedra sharing all six face nodes.
    let base = [
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [0.5, 0., 0.],
        [0.5, 0.5, 0.],
        [0., 0.5, 0.],
        [0., 0., 0.5],
        [0.5, 0., 0.5],
        [0., 0.5, 0.5],
        [0., 0., -1.],
        [0., 0., -0.5],
        [0., 0.5, -0.5],
        [0.5, 0., -0.5],
    ];
    for (i, xyz) in base.iter().enumerate() {
        m.set_node(i as u32 + 1, *xyz).unwrap();
    }
    m.set_element(1, ElementKind::Tet10, (1..=10).collect()).unwrap();
    let mut second_nodes = vec![1, 3, 2, 11, 7, 6, 5, 12, 13, 14];
    if disconnected {
        for node in &mut second_nodes {
            if [1, 2, 3, 5, 6, 7].contains(node) {
                let xyz = m.node(*node).unwrap().xyz;
                *node += 20;
                m.set_node(*node, xyz).unwrap();
            }
        }
    }
    m.set_element(2, ElementKind::Tet10, second_nodes).unwrap();
    m.add_set("First", Set { faces: vec![ElementFace { element: 1, face: 1 }], ..Set::default() }).unwrap();
    m.add_set("Second", Set { faces: vec![ElementFace { element: 2, face: 1 }], ..Set::default() }).unwrap();
    m.add_set("CutFace", Set { faces: vec![ElementFace { element: 1, face: 1 }], ..Set::default() }).unwrap();
    m.add_set("All", Set { elements: vec![1, 2], ..Set::default() }).unwrap();
    m.add_set("Anchors", Set { nodes: vec![4, 11], ..Set::default() }).unwrap();
    m.add_material("Steel", Material { elastic: Some(Elastic::Isotropic { young: 1000., poisson: 0.3 }), density: Some(1.), ..Default::default() })
        .unwrap();
    m.add_section(
        "Solid",
        Section { set: "All".into(), material: "Steel".into(), kind: SectionKind::Solid {}, orientation: None, formulation: Formulation::default() },
    )
    .unwrap();
    m.add_reference_point("Control", ReferencePoint { xyz: [0., 0., 0.], label: 100 }).unwrap();
    m.add_pretension(
        "Cut",
        fieldcraft_model::PretensionSection {
            first: "First".into(),
            second: "Second".into(),
            face: Some("CutFace".into()),
            normal: [0., 0., -1.],
            control: "Control".into(),
        },
    )
    .unwrap();
    m.add_boundary_condition("AnchorsFixed", BoundaryCondition::Fixed { target: Target::Set("Anchors".into()) }).unwrap();
    let mut before = Step::new("Before", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    let mut preload = Step::new("Preload", Procedure::Static {});
    let mut locked = Step::new("Locked", Procedure::Static {});
    let mut after = Step::new("After", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None });
    for s in [&mut before, &mut preload, &mut locked, &mut after] {
        s.boundary_conditions.insert("AnchorsFixed".into());
    }
    preload.pretension.insert("Cut".into(), fieldcraft_model::PretensionControl::Force { value: 1. });
    locked.pretension.insert("Cut".into(), fieldcraft_model::PretensionControl::Adjustment { value: 0.001 });
    after.base_state = Some("Locked".into());
    for s in [before, preload, locked, after] {
        m.add_step(s, None).unwrap();
    }
    m
}

#[test]
fn tet10_writer_hits_disconnected_topology_after_face_resolution() {
    let err = fieldcraft_decks::abaqus::write(&tet10_writer_fixture(true), fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(err.contains("coincident coordinates but no element connectivity"), "{err}");
    let deck = fieldcraft_decks::abaqus::write(&tet10_writer_fixture(false), fieldcraft_model::units::Preset::MKgS)
        .unwrap_or_else(|err| panic!("legal Tet10 writer failed: {err}"));
    assert!(deck.contains("*PRE-TENSION SECTION"));
    assert!(!deck.contains("*EQUATION"));
    assert_eq!(tet10_writer_fixture(false).nodes().len(), 14);
    assert!(deck.contains("*SURFACE, NAME=FC-PRET-Cut, TYPE=ELEMENT"));
    for step in ["Before", "Preload", "Locked", "After"] {
        assert!(deck.contains(step), "missing step {step}");
    }
}
