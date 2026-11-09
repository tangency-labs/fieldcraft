use fieldcraft_engine::Session;
use serde_json::{Value, json};

fn oscillator() -> Session {
    let mut s = Session::default();
    s.run_script(&json!([
        {"command":"model.units","params":{"system":"m-kg-s-N-Pa"}},
        {"command":"mesh.node.add","params":{"label":1,"xyz":[0,0,0]}},
        {"command":"model.set.create","params":{"name":"tip","nodes":[1]}},
        {"command":"model.mass.create","params":{"name":"m","target":{"set":"tip"},"mass":1}},
        {"command":"model.spring.create","params":{"name":"k","type":"spring","first":{"set":"tip"},"direction":{"dof":1},"stiffness":4}},
        {"command":"model.spring.create","params":{"name":"c","type":"dashpot","first":{"set":"tip"},"direction":{"dof":1},"coefficient":0.04}},
        {"command":"model.boundary_condition.create","params":{"name":"guide","type":"displacement","target":{"set":"tip"},"values":[null,0,0,null,null,null]}},
        {"command":"model.load.create","params":{"name":"force","type":"force","target":{"set":"tip"},"vector":[1,0,0]}},
        {"command":"model.step.create","params":{"name":"Sweep","procedure":{"type":"harmonic","frequencies":[0,1.0/std::f64::consts::PI],"method":{"type":"direct"},"load_factors":{"force":{"amplitude":1,"phase":"90 deg"}}},"boundary_conditions":["guide"],"loads":["force"],"output":["displacement"]}}
    ])).unwrap();
    s
}

#[test]
fn command_phasors_units_modal_and_round_trip() {
    let mut s = oscillator();
    let direct = s.execute("solve.harmonic", &json!({"threads":1})).unwrap();
    assert!(direct["frequencies"][1]["relative_residual"].as_f64().unwrap() < 1e-12);
    let p = s.execute("solve.probe", &json!({"step":"Sweep","nodes":[1]})).unwrap();
    let value = &p["nodes"][0]["samples"][1];
    assert!((value["displacement_real"][0].as_f64().unwrap() - 12.5).abs() < 1e-11, "{p}");
    assert!(value["displacement_imaginary"][0].as_f64().unwrap().abs() < 1e-11);
    assert_eq!(p["units"]["displacement_amplitude"], "m");
    s.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa", "overrides":{"angle":"deg"}})).unwrap();
    let shown = s.execute("solve.probe", &json!({"step":"Sweep","nodes":[1],"frame":1})).unwrap();
    assert!((shown["nodes"][0]["samples"][0]["displacement_amplitude"][0].as_f64().unwrap() - 12500.0).abs() < 1e-8);
    s.execute("model.step.edit", &json!({"name":"Sweep", "procedure":{"type":"harmonic", "frequencies":{"start":"0 Hz", "end":format!("{} Hz",1.0/std::f64::consts::PI),"count":2},"method":{"type":"modal","modes":1},"load_factors":{"force":{"phase":"90 deg"}}}})).unwrap();
    let modal = s.execute("solve.harmonic", &json!({"threads":2})).unwrap();
    assert!(modal["frequencies"][1]["relative_residual"].as_f64().unwrap() < 1e-12);
    let p = s.execute("solve.probe", &json!({"step":"Sweep","nodes":[1],"frame":1})).unwrap();
    assert!((p["nodes"][0]["samples"][0]["displacement_amplitude"][0].as_f64().unwrap() - 12500.0).abs() < 1e-8);
    assert!(s.execute("solve.probe", &json!({"step":"Sweep","nodes":[1],"frame":usize::MAX})).is_err());
    let result = s.result("Sweep").unwrap();
    let fieldcraft_engine::StepResult::Implicit(h) = result.as_ref() else { panic!() };
    let stored = fieldcraft_io::results::Solved::Implicit(h.clone());
    fieldcraft_io::results::check(&stored).unwrap();
    let bytes = serde_json::to_vec(h).unwrap();
    let back: fieldcraft_solver::StepResult = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(*h, back);
    use fieldcraft_io::results::{Decoded, Stored, StoredRef, decode, encode};
    let definition = [42; 32];
    let binary = encode(&[("Sweep", StoredRef { definition, result: stored.as_ref() })]).unwrap();
    assert_eq!(decode(&binary).unwrap(), Decoded::Results(vec![Stored { definition, result: stored.clone() }]));
    let mut damaged = h.clone();
    let fieldcraft_solver::StepResult::Harmonic(ref mut harmonic) = damaged else { panic!() };
    harmonic.frames[0].imaginary.clear();
    assert!(fieldcraft_io::results::check(&fieldcraft_io::results::Solved::Implicit(damaged)).is_err());
}

#[test]
fn harmonic_then_random_keeps_frequency_and_probe_indices() {
    let mut s = oscillator();
    let h = s.execute("solve.harmonic", &json!({"threads": 1})).unwrap();
    assert_eq!(h["frequencies"].as_array().unwrap().len(), 2);
    let r = s
        .execute(
            "solve.random",
            &json!({
            "natural_frequency": 1.0, "damping_ratio": 0.05, "base_acceleration_psd": 2.0,
            "frequencies": [0.0, 1.0], "transfer": [1.0, 1.0], "base_psd": [2.0, 2.0],
                    "probe_indices": [0], "step": "Sweep"
                }),
        )
        .unwrap();
    assert_eq!(r["frequencies"], json!([0.0, 1.0]));
    assert_eq!(r["probe_indices"], json!([0]));
    assert_eq!(r["units"], Value::Null);
    assert!(r["displacement_rms"].as_f64().unwrap() > 0.0);
}

#[test]
fn random_projects_input_psd_onto_stored_harmonic_frames() {
    let mut s = oscillator();
    s.execute("solve.harmonic", &json!({"threads": 1})).unwrap();
    let r = s
        .execute(
            "solve.random",
            &json!({
                "step": "Sweep",
                "input_psd": [[0.0, 2.0], [1.0 / std::f64::consts::PI, 4.0]],
            }),
        )
        .unwrap();
    assert_eq!(r["method"], "harmonic_frame_psd");
    assert_eq!(r["frequencies"], json!([0.0, 1.0 / std::f64::consts::PI]));
    assert_eq!(r["input_psd"], json!([2.0, 4.0]));
    assert!(r["node_displacement_rms"][0].as_f64().unwrap() > 0.0);
}

#[test]
fn hostile_harmonic_inputs_rollback_and_decks_name_limits() {
    let mut s = oscillator();
    for frequencies in [json!([]), json!([-1]), json!([1, 1]), json!({"start":0,"end":1,"count":usize::MAX})] {
        let before = s.model.clone();
        assert!(
            s.execute(
                "model.step.edit",
                &json!({"name":"Sweep","procedure":{"type":"harmonic","frequencies":frequencies,"method":{"type":"direct"}}})
            )
            .is_err()
        );
        assert_eq!(s.model, before);
    }
    for format in ["nastran", "lsdyna"] {
        let error = fieldcraft_decks::FORMATS
            .iter()
            .find(|f| f.name == format)
            .unwrap()
            .write(&s.model, fieldcraft_model::units::Preset::MKgS)
            .unwrap_err()
            .to_string();
        assert!(error.contains("harmonic") && error.contains("Sweep"), "{error}");
    }
    let deck = fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MKgS).unwrap();
    assert!(deck.contains("*STEADY STATE DYNAMICS, DIRECT"));
    assert!(deck.contains("*CLOAD, IMAGINARY, OP=NEW"));
    s.execute(
        "model.step.edit",
        &json!({"name":"Sweep", "procedure":{"type":"harmonic", "frequencies":[0], "method":{"type":"modal","modes":1,"residual_vectors":1}}}),
    )
    .unwrap();
    let result = s.execute("solve.harmonic", &json!({})).unwrap();
    assert_eq!(result["residual_vectors"], 0, "a full basis needs no correction vector");
    let error = fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(error.contains("Sweep") && error.contains("residual-vector"), "{error}");
    let before = s.model.clone();
    assert!(
        s.execute(
            "model.step.edit",
            &json!({"name":"Sweep", "procedure":{"type":"harmonic", "frequencies":[0], "method":{"type":"modal","modes":1,"residual_vectors":17}}})
        )
        .is_err()
    );
    assert_eq!(s.model, before);
}

#[test]
fn material_damping_units_and_unsupported_exports_do_not_drop_coefficients() {
    let mut s = oscillator();
    s.run_script(&json!([
        {"command":"model.units","params":{"system":"mm-kg-ms-kN-GPa"}},
        {"command":"mesh.node.add","params":{"label":2,"xyz":[1000,0,0]}},
        {"command":"mesh.element.add","params":{"label":1,"kind":"rod2","nodes":[1,2]}},
        {"command":"model.set.create","params":{"name":"member","elements":[1]}},
        {"command":"model.material.create","params":{"name":"Damped","elastic":{"type":"isotropic","young":1,"poisson":0.3},"density":1,"damping":{"alpha":0.002,"beta":10,"loss_factor":0.03}}},
        {"command":"model.section.create","params":{"name":"Bar","set":"member","material":"Damped","kind":{"type":"rod","area":1}}},
        {"command":"model.step.edit","params":{"name":"Sweep","procedure":{"type":"static"}}}
    ])).unwrap();
    let d = s.model.material("Damped").unwrap().damping.unwrap();
    assert!((d.alpha - 2.0).abs() < 1e-14);
    assert!((d.beta - 0.01).abs() < 1e-14);
    for format in ["nastran", "lsdyna"] {
        let error = fieldcraft_decks::FORMATS
            .iter()
            .find(|f| f.name == format)
            .unwrap()
            .write(&s.model, fieldcraft_model::units::Preset::MKgS)
            .unwrap_err()
            .to_string();
        assert!(error.contains("Damped") && error.contains("damping"), "{error}");
    }
    let deck = fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MKgS).unwrap();
    let damping = deck.lines().find(|line| line.starts_with("*DAMPING,")).unwrap();
    let values: std::collections::BTreeMap<_, _> = damping
        .split(',')
        .skip(1)
        .map(|part| {
            let (name, value) = part.trim().split_once('=').unwrap();
            (name, value.parse::<f64>().unwrap())
        })
        .collect();
    assert_eq!(values.len(), 3);
    for (name, expected) in [("ALPHA", 2.0), ("BETA", 0.01), ("STRUCTURAL", 0.03)] {
        assert!((values[name] - expected).abs() < 1e-14, "{damping}");
    }
}

#[test]
fn harmonic_complex_memory_cap_is_explicit_and_preserves_results_on_refusal() {
    let mut s = oscillator();
    let expected = s.execute("solve.harmonic", &json!({"step":"Sweep"})).unwrap();
    let actual = s.execute("solve.harmonic", &json!({"step":"Sweep","max_factor_entries":3_200_000_000u64})).unwrap();
    assert_eq!(expected, actual);
    let count = s.solved_steps().count();
    let error = s.execute("solve.harmonic", &json!({"step":"Sweep","max_factor_entries":4_000_000_001u64})).unwrap_err().to_string();
    assert!(error.contains("max_factor_entries"), "{error}");
    assert_eq!(s.solved_steps().count(), count);
    let schema = fieldcraft_engine::schema::json_schema(fieldcraft_engine::find_command("solve.harmonic").unwrap().form.unwrap());
    assert_eq!(schema["properties"]["max_factor_entries"]["maximum"], 4_000_000_000u64);
}

#[test]
fn random_base_motion_adds_coherent_acceleration_and_validates_inputs() {
    let mut s = oscillator();
    s.execute("solve.harmonic", &json!({"threads":1})).unwrap();
    let p = json!({"step":"Sweep","nodes":[1],"input_psd":[[0.,2.],[1.,2.]],"base_direction":[1.,0.,0.]});
    let r = s.execute("solve.random", &p).unwrap();
    assert_eq!(r["probes"][0]["samples"][0]["acceleration_components_psd"], json!([2., 0., 0.]));
    assert_eq!(r["acceleration_reference"], "absolute_base_motion");
    let probe = s.execute("solve.probe", &json!({"step":"Sweep","nodes":[1]})).unwrap();
    let sample = &probe["nodes"][0]["samples"][1];
    let w2 = (2. * std::f64::consts::PI * sample["frequency"].as_f64().unwrap()).powi(2);
    let re = sample["displacement_real"][0].as_f64().unwrap();
    let im = sample["displacement_imaginary"][0].as_f64().unwrap();
    let expected = ((1. - w2 * re).powi(2) + (w2 * im).powi(2)) * 2.;
    assert!((r["probes"][0]["samples"][1]["acceleration_components_psd"][0].as_f64().unwrap() - expected).abs() < 1e-12);
    for (key, value) in
        [("nodes", json!([1, 1])), ("nodes", json!([999])), ("base_direction", json!([0, 0, 0])), ("input_psd", json!([[0., 1.], [1., -1.]]))]
    {
        let mut bad = p.clone();
        bad[key] = value;
        assert!(s.execute("solve.random", &bad).is_err());
    }
}

#[test]
fn solid_element_stress_locations_phasors_psd_and_reopen() {
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let xyz = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]];
    for (i, x) in xyz.iter().enumerate() {
        s.execute("mesh.node.add", &json!({"label":i+1,"xyz":x})).unwrap();
    }
    s.run_script(&json!([
        {"command":"mesh.element.add","params":{"label":1,"kind":"hex8","nodes":[1,2,3,4,5,6,7,8]}},
        {"command":"model.set.create","params":{"name":"body","elements":[1],"nodes":[1,2,3,4,5,6,7,8]}},
        {"command":"model.set.create","params":{"name":"root","nodes":[1,4,5,8]}},
        {"command":"model.set.create","params":{"name":"tip","nodes":[2,3,6,7]}},
        {"command":"model.material.create","params":{"name":"m","elastic":{"type":"isotropic","young":100,"poisson":0},"density":1}},
        {"command":"model.section.create","params":{"name":"s","set":"body","material":"m","kind":{"type":"solid"}}},
        {"command":"model.boundary_condition.create","params":{"name":"root","type":"displacement","target":{"set":"root"},"values":[0,null,null,null,null,null]}},
        {"command":"model.boundary_condition.create","params":{"name":"guide","type":"displacement","target":{"set":"body"},"values":[null,0,0,null,null,null]}},
        {"command":"model.load.create","params":{"name":"f","type":"force","target":{"set":"tip"},"vector":[1,0,0],"distribution":"total"}},
        {"command":"model.step.create","params":{"name":"H","procedure":{"type":"harmonic","frequencies":[0,0.1],"method":{"type":"direct"}},"boundary_conditions":["root","guide"],"loads":["f"]}}
    ])).unwrap();
    s.execute("solve.harmonic", &json!({"step":"H","threads":1})).unwrap();
    let params = json!({"step":"H","input_psd":[[0,2],[0.1,2]],"nodes":[2],"elements":[1],"element_phasors":true});
    let out = s.execute("solve.random", &params).unwrap();
    let e = &out["element_stress"];
    assert_eq!(e["locations"].as_array().unwrap().len(), 16);
    for (i, location) in e["locations"].as_array().unwrap().iter().enumerate() {
        assert_eq!(location["element"], 1);
        if i < 8 {
            assert_eq!(location["node"], i + 1);
            assert_eq!(location["xyz"], json!(xyz[i]));
        } else {
            for j in 0..3 {
                let xi = location["natural_coordinates"][j].as_f64().unwrap();
                assert!((xi.abs() - 1. / 3_f64.sqrt()).abs() < 1e-14);
                assert!((location["xyz"][j].as_f64().unwrap() - (1. + xi) / 2.).abs() < 1e-14);
            }
        }
        assert!((e["stress_real"][0][i][0].as_f64().unwrap() - 1.).abs() < 1e-12);
        assert!(e["stress_imaginary"][0][i][0].as_f64().unwrap().abs() < 1e-12);
        let p0 = e["components_psd"][0][i][0].as_f64().unwrap();
        let p1 = e["components_psd"][1][i][0].as_f64().unwrap();
        assert!((p0 - 2.).abs() < 1e-12);
        assert!((e["components_rms"][i][0].as_f64().unwrap().powi(2) - 0.05 * (p0 + p1)).abs() < 1e-12);
    }
    let dir = std::env::temp_dir().join(format!("fieldcraft-element-stress-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("model.fieldcraft");
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(s.execute("solve.random", &params).unwrap(), out);
    assert!(s.execute("solve.random", &json!({"step":"H","input_psd":[[0,2],[0.1,2]],"elements":[1,1]})).is_err());
    // Constitutive patch: two independently labelled overlapping solids share
    // every node, but their moduli differ. A manufactured affine field must
    // retain both stresses, never their shared-node average.
    use fieldcraft_model::{Elastic, ElementKind, Material, Section, SectionKind, Set};
    let mut model = s.model.clone();
    model.set_element(2, ElementKind::Hex8, (1..=8).collect()).unwrap();
    model.add_set("second", Set { elements: vec![2], ..Set::default() }).unwrap();
    model
        .add_material("twice", Material { elastic: Some(Elastic::Isotropic { young: 200., poisson: 0. }), density: Some(1.), ..Material::default() })
        .unwrap();
    model
        .add_section(
            "second",
            Section {
                set: "second".into(),
                material: "twice".into(),
                kind: SectionKind::Solid {},
                orientation: None,
                formulation: Default::default(),
            },
        )
        .unwrap();
    let fieldcraft_engine::StepResult::Implicit(fieldcraft_solver::StepResult::Harmonic(h)) = s.result("H").unwrap().as_ref() else { panic!() };
    let recovery = fieldcraft_solver::harmonic::ElementStressRecovery::new(&model, h, &[1, 2], &fieldcraft_solver::Limits::default()).unwrap();
    let affine = xyz.iter().map(|x| [0.01 * x[0] + 0.02 * x[1] + 0.03 * x[2], 0.04 * x[1] + 0.05 * x[2], 0.06 * x[2]]).collect::<Vec<_>>();
    let values = recovery.recover(&affine).unwrap();
    for (location, stress) in recovery.locations.iter().zip(values) {
        for (a, b) in stress.iter().zip([1., 4., 6., 1., 1.5, 2.5]) {
            assert!((a - b * f64::from(location.element)).abs() < 1e-12, "{location:?}: {stress:?}");
        }
    }
    std::fs::remove_dir_all(dir).unwrap();
}
