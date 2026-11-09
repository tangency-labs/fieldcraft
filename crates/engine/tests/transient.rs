//! End-to-end original rod model: consistent free mass rho*A*L/3=1, K=EA/L=4.
//! This checks command -> FE assembly -> integration -> stress -> file -> probe.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn rod() -> Session {
    let mut s = Session::default();
    let recipe: Value = serde_json::from_str(include_str!("data/transient-rod.json")).unwrap();
    s.run_script(&recipe).unwrap();
    s
}
fn probe(s: &mut Session) -> Value {
    s.execute("solve.probe", &json!({"step":"Pulse","nodes":[2]})).unwrap()
}
#[test]
fn finite_element_transient_round_trip_and_display_units() {
    let mut s = rod();
    s.execute("solve.transient", &json!({"step":"Pulse","threads":1})).unwrap();
    let before = probe(&mut s);
    for sample in before["nodes"][0]["samples"].as_array().unwrap() {
        let t = sample["time"].as_f64().unwrap();
        let u = (1.0 - (2.0 * t).cos()) / 4.0;
        assert!((sample["displacement"][0].as_f64().unwrap() - u).abs() < 2e-7);
        assert!((sample["stress"][0].as_f64().unwrap() - 4.0 * u).abs() < 8e-7);
    }
    assert!((before["nodes"][0]["samples"][0]["acceleration"][0].as_f64().unwrap() - 1.0).abs() < 1e-12);
    let path = std::env::temp_dir().join(format!("fieldcraft-transient-roundtrip-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut restored = Session::default();
    restored.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(probe(&mut restored), before);
    restored.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let after = probe(&mut restored);
    let a = &after["nodes"][0]["samples"][10];
    let b = &before["nodes"][0]["samples"][10];
    assert_eq!(a["time"], 1000.0);
    assert!((a["displacement"][0].as_f64().unwrap() - 1000.0 * b["displacement"][0].as_f64().unwrap()).abs() < 1e-10);
    assert!((a["acceleration"][0].as_f64().unwrap() - 0.001 * b["acceleration"][0].as_f64().unwrap()).abs() < 1e-12);
    let _ = std::fs::remove_file(path);
}
#[test]
fn modal_transient_matches_direct_on_same_finite_element_model() {
    let mut s = rod();
    s.execute("solve.transient", &json!({"threads":1})).unwrap();
    let direct = probe(&mut s);
    s.execute("model.step.edit",&json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.001,"increments":1000,"output_every":100,"method":{"type":"modal","modes":1,"damping_ratios":[0]}}})).unwrap();
    s.execute("solve.run", &json!({})).unwrap();
    let modal = probe(&mut s);
    for (d, m) in direct["nodes"][0]["samples"].as_array().unwrap().iter().zip(modal["nodes"][0]["samples"].as_array().unwrap()) {
        assert!((d["displacement"][0].as_f64().unwrap() - m["displacement"][0].as_f64().unwrap()).abs() < 1e-11);
        assert!((d["stress"][0].as_f64().unwrap() - m["stress"][0].as_f64().unwrap()).abs() < 1e-11);
    }
}
#[test]
fn finite_element_spectrum_participation_and_persistence() {
    let mut s = rod();
    // Consistent full mass gives free row M*r=1.5, not free diagonal Mff=1.
    s.execute("model.step.edit",&json!({"name":"Pulse","loads":[],"procedure":{"type":"response_spectrum","modes":1,"damping_ratio":0.05,"combination":"cqc","direction":[1,0,0],"pseudo_acceleration":[[0.1,2.0],[1.0,2.0]]}})).unwrap();
    s.execute("solve.spectrum", &json!({"threads":1})).unwrap();
    let p = probe(&mut s);
    let sample = &p["nodes"][0]["samples"][0];
    assert!(sample["time"].is_null());
    assert!((sample["displacement"][0].as_f64().unwrap() - 0.75).abs() < 1e-10);
    assert!((sample["stress"][0].as_f64().unwrap() - 3.0).abs() < 1e-10);
    let r = s.result("Pulse").unwrap();
    let fieldcraft_engine::StepResult::Implicit(r) = r.as_ref() else { panic!() };
    let stored = fieldcraft_io::results::Solved::Implicit(r.clone());
    use fieldcraft_io::results::{Decoded, Stored, StoredRef, decode, encode};
    let data = encode(&[("Pulse", StoredRef { definition: [7; 32], result: stored.as_ref() })]).unwrap();
    assert_eq!(decode(&data).unwrap(), Decoded::Results(vec![Stored { definition: [7; 32], result: stored }]));
    let mut bad = r.clone();
    let fieldcraft_solver::StepResult::Dynamics(ref mut d) = bad else { panic!() };
    d.frames[0].time = Some(0.0);
    assert!(fieldcraft_io::results::check(&fieldcraft_io::results::Solved::Implicit(bad)).is_err());
}
#[test]
fn initial_constraints_and_inactive_amplitudes_are_rejected() {
    let mut s = rod();
    assert!(s.execute("model.step.edit",&json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.01,"increments":10,"load_amplitudes":{"missing":{"type":"tabular","points":[[0,1],[1,1]]}}}})).is_err());
    s.execute(
        "model.step.edit",
        &json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.01,"increments":10,"initial_displacement":{"1":[0.1,0,0]}}}),
    )
    .unwrap();
    assert!(s.execute("solve.transient", &json!({})).is_err());
    assert!(s.result("Pulse").is_none());
}

#[test]
fn damping_interpolation_period_units_and_direction_envelopes() {
    let mut s = rod();
    s.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let mut procedure = json!({"type":"response_spectrum","modes":1,"damping_ratio":0.05,
        "combination":"cqc","directions":[{"direction":[0.6,0.8,0],"scale":1},{"direction":[-0.8,0.6,0],"scale":1}],
        "direction_combination":"srss","period_tables":[
            {"damping_ratio":0.02,"points":[["1000 ms","2 m/s^2"],["5000 ms","2 m/s^2"]]},
            {"damping_ratio":0.08,"points":[["1000 ms","4 m/s^2"],["5000 ms","4 m/s^2"]]}]});
    s.execute("model.step.edit", &json!({"name":"Pulse","loads":[],"procedure":procedure})).unwrap();
    s.execute("solve.spectrum", &json!({})).unwrap();
    let p = probe(&mut s);
    assert!((p["nodes"][0]["samples"][0]["displacement"][0].as_f64().unwrap() - 1125.0).abs() < 1e-7);
    procedure["direction_combination"] = json!("percent1004040");
    s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).unwrap();
    s.execute("solve.spectrum", &json!({})).unwrap();
    let p = probe(&mut s);
    assert!((p["nodes"][0]["samples"][0]["displacement"][0].as_f64().unwrap() - 1170.0).abs() < 1e-7);
    // Extrapolation and non-orthogonal direction sets are rejected before solving.
    procedure["damping_ratio"] = json!(0.1);
    assert!(s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).is_err());
    procedure["damping_ratio"] = json!(0.05);
    procedure["directions"][1]["direction"] = json!([0.6, 0.8, 0]);
    assert!(s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).is_err());
}
