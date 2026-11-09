//! Original FE oscillator -> stored dynamics -> exported SI stress -> rainflow/Miner.
//! K=4, M=1, unit step: sigma(t)=1-cos(2t). On [0,1] this is monotonic,
//! hence one residual half cycle, amplitude=sigma(1)/2, N=100/amplitude^3.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

fn rod() -> Session {
    let mut s = Session::default();
    s.run_script(&serde_json::from_str::<Value>(include_str!("data/transient-rod.json")).unwrap()).unwrap();
    s
}

fn probe(s: &mut Session) -> Value {
    s.execute("solve.probe", &json!({"step":"Pulse","nodes":[2]})).unwrap()
}

fn request(history: Value) -> Value {
    json!({"history":history,"source_id":"original FE rod: E4 rho3 A1 L1 unit step, consistent mass",
        "step":"Pulse","nodes":[2],"max_time_step_seconds":0.1,"weights":[1,0,0,0,0,0],
        "curve":{"amplitude_ref":1,"cycles_ref":100,"exponent":3,"endurance_limit":0},
        "correction":{"type":"none"}})
}

#[test]
fn stored_transient_stress_drives_fatigue_without_relabeling_or_reconstruction() {
    let mut s = rod();
    s.execute("solve.transient", &json!({"step":"Pulse","threads":1})).unwrap();
    let history = probe(&mut s);
    assert_eq!(history["units"]["stress"], "Pa");
    assert_eq!(history["units"]["time"], "s");
    let samples = history["nodes"][0]["samples"].as_array().unwrap();
    assert_eq!(samples.len(), 11);
    for sample in samples {
        let t = sample["time"].as_f64().unwrap();
        assert!((sample["stress"][0].as_f64().unwrap() - (1. - (2. * t).cos())).abs() < 8e-7);
    }
    let p = request(history.clone());
    let fatigue = s.execute("solve.fatigue.transient", &p).unwrap();
    assert_eq!(fatigue["source_id"], p["source_id"]);
    let node = &fatigue["nodes"][0];
    assert_eq!(node["cycle_count"], 0.5);
    assert_eq!(node["cycles"].as_array().unwrap().len(), 1);
    let amplitude = samples.last().unwrap()["stress"][0].as_f64().unwrap() / 2.;
    let damage = node["damage"].as_f64().unwrap();
    assert!((damage - 0.5 * amplitude * amplitude * amplitude / 100.).abs() < 1e-14);
    let exact_amplitude = (1. - 2_f64.cos()) / 2.;
    assert!((damage / (0.5 * exact_amplitude * exact_amplitude * exact_amplitude / 100.) - 1.).abs() < 2e-6);

    let fieldcraft_engine::StepResult::Implicit(result) = s.result("Pulse").unwrap().as_ref() else { panic!("implicit dynamics") };
    assert!(matches!(result, fieldcraft_solver::StepResult::Dynamics(_)));
    use fieldcraft_io::results::{Decoded, Solved, Stored, StoredRef, decode, encode};
    let stored = Solved::Implicit(result.clone());
    let bytes = encode(&[("Pulse", StoredRef { definition: [7; 32], result: stored.as_ref() })]).unwrap();
    assert_eq!(&bytes[..8], b"MCRSLT07");
    assert_eq!(decode(&bytes).unwrap(), Decoded::Results(vec![Stored { definition: [7; 32], result: stored }]));
    let path = std::env::temp_dir().join(format!("fieldcraft-integrated-transient-fatigue-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut reopened = Session::default();
    reopened.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(probe(&mut reopened), history);
    // Fatigue is an exported calculation, not a stored dynamics field. Recompute from
    // the actual reopened probe, keeping source identity and all declared SI metadata.
    let p = request(probe(&mut reopened));
    assert_eq!(reopened.execute("solve.fatigue.transient", &p).unwrap(), fatigue);
    reopened.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
    let displayed = request(probe(&mut reopened));
    let error = reopened.execute("solve.fatigue.transient", &displayed).unwrap_err().to_string();
    assert!(error.contains("seconds and Pa"), "{error}");
    // The unchanged SI export is independent of the receiving session's display units.
    assert_eq!(reopened.execute("solve.fatigue.transient", &p).unwrap(), fatigue);
    std::fs::remove_file(path).unwrap();
}

#[test]
fn actual_spectrum_envelope_is_not_a_time_history_for_fatigue() {
    let mut s = rod();
    s.execute(
        "model.step.edit",
        &json!({"name":"Pulse","loads":[],"procedure":{
        "type":"response_spectrum","modes":1,"damping_ratio":0.05,"combination":"cqc",
        "direction":[1,0,0],"pseudo_acceleration":[[0.1,2.0],[1.0,2.0]]}}),
    )
    .unwrap();
    s.execute("solve.spectrum", &json!({"step":"Pulse","threads":1})).unwrap();
    let before = probe(&mut s);
    assert_eq!(before["spectrum"], true);
    assert!(before["frames"][0]["time"].is_null());
    let error = s.execute("solve.fatigue.transient", &request(before.clone())).unwrap_err().to_string();
    assert!(error.contains("response spectra"), "{error}");
    assert_eq!(probe(&mut s), before);
}
