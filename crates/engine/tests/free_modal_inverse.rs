//! Free–free rectangular beam, EB beta_L^2=22.373 (first bending mode).
//! Synthetic inverse observations are independently solved at the true FE material.
use fieldcraft_engine::Session;
use fieldcraft_model::*;
use fieldcraft_solver::{Limits, solve_frequency};
use serde_json::{Value, json};

const E: f64 = 1000.0;
const K: f64 = 0.00001;
fn beam(young: f64, spring: Option<f64>) -> Session {
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let m = &mut s.model;
    let n = 20u32;
    for i in 0..=n {
        m.set_node(i + 1, [10.0 * i as f64 / n as f64, 0.0, 0.0]).unwrap();
    }
    for i in 1..=n {
        m.set_element(i, ElementKind::Beam2, vec![i, i + 1]).unwrap();
    }
    m.add_set("all", Set { nodes: (1..=n + 1).collect(), elements: (1..=n).collect(), ..Default::default() }).unwrap();
    m.add_set("tip", Set { nodes: vec![1], ..Default::default() }).unwrap();
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young, poisson: 0.3 }), density: Some(1.0), ..Default::default() }).unwrap();
    m.add_section(
        "beam",
        Section {
            set: "all".into(),
            material: "Mat".into(),
            kind: SectionKind::Beam { profile: BeamProfile::Rectangle { width: 0.1, height: 0.2 }, orientation: [0.0, 1.0, 0.0] },
            orientation: None,
            formulation: Default::default(),
        },
    )
    .unwrap();
    if let Some(k) = spring {
        m.add_spring(
            "soft",
            Spring::Spring { first: Target::Set("tip".into()), second: None, direction: fieldcraft_model::discrete::Direction::Dof(2), stiffness: k },
        )
        .unwrap();
    }
    m.add_step(Step::new("F", Procedure::Frequency { modes: 8, min_hz: None, max_hz: None }), None).unwrap();
    s
}
fn modes(s: &Session) -> fieldcraft_solver::FrequencyResult {
    solve_frequency(&s.model, Some("F"), &Default::default(), &Limits::default()).unwrap()
}
fn parameter() -> Value {
    json!({"type":"material","name":"Mat","constant":"young"})
}
fn joint(f: f64, soft: f64, matched: bool, shape: Value) -> Value {
    let freq =
        if matched { json!({"kind":"matched_frequency","group":"g","measured_mode":0}) } else { json!({"kind":"frequency","step":"F","mode":6}) };
    let mut p = json!({"parameters":[{"parameter":parameter(),"lower_si":800.,"upper_si":1200.,"scale_si":1000.},{"parameter":{"type":"spring_stiffness","name":"soft"},"lower_si":K*0.5,"upper_si":K*1.5,"scale_si":K}],"observations":[{"response":freq,"value_si":f,"standard_deviation_si":f*0.001},{"response":{"kind":"frequency","step":"F","mode":5},"value_si":soft,"standard_deviation_si":soft*0.001}],"max_iterations":64,"relative_increment":1e-4});
    if matched {
        p["modal_observations"] = json!([{"id":"g","step":"F","sensors":(1..=21).map(|i|json!({"node":i,"dof":2})).collect::<Vec<_>>(),"measured_modes":[shape],"min_mac":0.99,"min_assignment_margin":0.001,"normalization":{"kind":"reference","reference_sensor":0}}]);
    }
    p
}
fn check(r: &Value, label: &str) {
    for (i, t) in [E, K].iter().enumerate() {
        let error = (r["estimate_si"][i].as_f64().unwrap() / t - 1.).abs();
        eprintln!("{label} parameter={i} error={error:.12e}");
        assert!(error <= 1e-6, "{r}");
    }
}
#[test]
fn free_beam_three_command_paths_recover_e() {
    let truth = beam(E, None);
    let target = modes(&truth);
    let eb = 22.373 / (std::f64::consts::TAU * 100.) * (E * 0.1_f64.powi(2) / 12.).sqrt();
    let f = target.modes[6].frequency;
    eprintln!(
        "free frequencies {:?}; EB={eb:.12e}, FE={f:.12e}, error={:.12e}",
        target.modes.iter().map(|m| m.frequency).collect::<Vec<_>>(),
        (f / eb - 1.).abs()
    );
    assert!(target.modes[..6].iter().all(|m| m.frequency < f * 1e-5));
    assert!((f / eb - 1.).abs() <= 0.01);
    for uniform in [false, true] {
        let mut initial = beam(E * 0.9, None);
        let p = json!({"step":"F","parameter":parameter(),"lower":800.,"upper":1200.,"observed":[f],"mode_indices":[6],"standard_deviation":[f*0.001],"max_iterations":64,"relative_increment":1e-4,"uniform_young_scaling":uniform});
        let r = initial.execute("solve.inverse.frequency", &p).unwrap();
        let error = (r["estimate"].as_f64().unwrap() / E - 1.).abs();
        eprintln!("scalar uniform={uniform} E error={error:.12e}");
        assert!(error <= 1e-6, "{r}");
        let mut bad = p;
        bad["mode_indices"] = json!([0]);
        let err = initial.execute("solve.inverse.frequency", &bad).unwrap_err().to_string();
        assert!(err.contains("near-zero") && err.contains("index 0"), "{err}");
    }
    let target = modes(&beam(E, Some(K)));
    let f = target.modes[6].frequency;
    eprintln!("soft free frequencies {:?}", target.modes.iter().map(|m| m.frequency).collect::<Vec<_>>());
    assert!(target.modes[..5].iter().all(|m| m.frequency < f * 1e-5));
    assert!((f / eb - 1.).abs() <= 0.01);
    let shape = json!(target.modes[6].shape.iter().map(|u| u[1]).collect::<Vec<_>>());
    for matched in [false, true] {
        let mut initial = beam(E * 0.9, Some(K * 1.1));
        let p = joint(f, target.modes[5].frequency, matched, shape.clone());
        let r = initial.execute("solve.inverse.joint", &p).unwrap();
        check(&r, if matched { "joint MAC" } else { "joint Frequency" });
        if matched {
            assert_eq!(r["modal_pairings"][0]["pairing"]["excluded_near_zero_modes"], json!([0, 1, 2, 3, 4]));
            assert_eq!(r["modal_pairings"][0]["pairing"]["pairs"][0]["computed_mode"], 6);
        }
    }
}

#[test]
fn mac_removal_reports_original_mode_and_keeps_good_fit() {
    let target = modes(&beam(E, Some(K)));
    let good = json!(target.modes[6].shape.iter().map(|u| u[1]).collect::<Vec<_>>());
    let mut p = joint(target.modes[6].frequency, target.modes[5].frequency, true, good.clone());
    // Alternating sensor readings have no low-frequency beam counterpart.
    let bad = json!((0..21).map(|i| if i % 2 == 0 { 1.0 } else { -1.0 }).collect::<Vec<_>>());
    // Put unmatched row first to prove original indices are preserved.
    p["modal_observations"][0]["measured_modes"] = json!([bad, good]);
    p["observations"][0]["response"]["measured_mode"] = json!(1);
    p["observations"]
        .as_array_mut()
        .unwrap()
        .push(json!({"response":{"kind":"matched_frequency","group":"g","measured_mode":0},"value_si":999.,"standard_deviation_si":0.001}));
    let mut initial = beam(E * 0.9, Some(K * 1.1));
    let err = initial.execute("solve.inverse.joint", &p).unwrap_err().to_string();
    assert!(err.contains("MAC"), "{err}");
    p["modal_observations"][0]["allow_mac_fallback"] = json!(true);
    let r = initial.execute("solve.inverse.joint", &p).unwrap();
    check(&r, "MAC removal");
    let report = &r["modal_pairings"][0]["pairing"];
    eprintln!("MAC removal report {report}");
    assert_eq!(report["removed_modes"].as_array().unwrap().len(), 1);
    assert_eq!(report["removed_modes"][0]["measured_mode"], 0);
    assert_eq!(report["removed_modes"][0]["reason"], "below_min_mac");
    assert!(report["removed_modes"][0]["best_mac"].as_f64().unwrap() < 0.99);
    assert_eq!(report["measured_indices"], json!([1]));
    assert_eq!(r["predicted_si"].as_array().unwrap().len(), 2);
}
