use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn input() -> Value {
    let samples: Vec<_> = [0.0, 2.0, 0.0].iter().enumerate().map(|(i, x)| json!({"time":i as f64,"stress":[x,0,0,0,0,0]})).collect();
    json!({"history":{"step":"Pulse","spectrum":false,"units":{"stress":"Pa","time":"s"},"stress_components":["xx","yy","zz","xy","xz","yz"],"stress_location":"global nodal average; shell top surface","frames":[{"time":0},{"time":1},{"time":2}],"nodes":[{"node":2,"samples":samples}]},"source_id":"original analytic triangle; external mesh A","step":"Pulse","nodes":[2],"max_time_step_seconds":1,"weights":[1,0,0,0,0,0],"curve":{"amplitude_ref":1,"cycles_ref":100,"exponent":3,"endurance_limit":0},"correction":{"type":"none"}})
}
#[test]
fn transient_triangle_counts_residuals_and_keeps_source_identity() {
    let mut s = Session::default();
    let p = input();
    let r = s.execute("solve.fatigue.transient", &p).unwrap();
    assert_eq!(r["source_id"], p["source_id"]);
    assert_eq!(r["nodes"][0]["cycle_count"], 1.0);
    assert!((r["nodes"][0]["damage"].as_f64().unwrap() - 0.01).abs() < 1e-14);
    let cycles = r["nodes"][0]["cycles"].as_array().unwrap();
    assert_eq!(cycles.len(), 2);
    assert_eq!(cycles[0]["count"], 0.5);
    assert_eq!(cycles[0]["mean_pa"], 1.0);
    assert_eq!(cycles.last().unwrap()["cumulative_damage"], r["nodes"][0]["damage"]);
    s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
    assert_eq!(s.execute("solve.fatigue.transient", &p).unwrap(), r);
}
#[test]
fn transient_contract_rejects_spectra_wrong_identity_units_time_and_missing_stress() {
    let changes = [
        ("/history/spectrum", json!(true)),
        ("/history/units/stress", json!("MPa")),
        ("/history/units/time", json!("ms")),
        ("/history/step", json!("Other")),
        ("/history/nodes/0/node", json!(3)),
        ("/history/frames/1/time", json!(0)),
        ("/history/nodes/0/samples/1/time", json!(0.5)),
        ("/history/nodes/0/samples/1/stress", Value::Null),
        ("/history/stress_components", json!(["xx", "yy", "zz", "xy", "yz", "xz"])),
        ("/max_time_step_seconds", json!(0.5)),
        ("/nodes", json!([2, 2])),
    ];
    for (path, value) in changes {
        let mut p = input();
        *p.pointer_mut(path).unwrap() = value;
        assert!(Session::default().execute("solve.fatigue.transient", &p).is_err(), "accepted {path}");
    }
}
