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

fn tensor_input(kind: &str, tensor: [f64; 6]) -> Value {
    let mut p = input();
    p.as_object_mut().unwrap().remove("weights");
    p["criterion"] = json!({"type":kind});
    p["curve"]["exponent"] = json!(2);
    p["history"]["nodes"][0]["samples"] =
        json!([-1., 1., -1.].iter().enumerate().map(|(i, f)| json!({"time":i as f64,"stress":tensor.map(|v|v*f)})).collect::<Vec<_>>());
    p
}
#[test]
fn complete_tensor_criteria_count_scalar_means_and_miner() {
    // Eigenvalues 4,0,-1; rotated tensor has the same invariants.
    for tensor in [[3., 0., 0., 2., 0., 0.], [-0.5, 3.5, 0., 1.5, 0., 0.]] {
        for (kind, amplitude, mean, damage) in [("signed_von_mises", 21_f64.sqrt(), 0., 0.21), ("max_principal", 1.5, 2.5, 0.04)] {
            let mut p = tensor_input(kind, tensor);
            if kind == "max_principal" {
                p["correction"] = json!({"type":"goodman","ultimate":10});
            }
            let r = Session::default().execute("solve.fatigue.transient", &p).unwrap();
            assert!((r["nodes"][0]["damage"].as_f64().unwrap() - damage).abs() < 1e-13);
            assert_eq!(r["nodes"][0]["cycle_count"], 1.0);
            for cycle in r["nodes"][0]["cycles"].as_array().unwrap() {
                assert_eq!(cycle["count"], 0.5);
                assert!((cycle["mean_pa"].as_f64().unwrap() - mean).abs() < 1e-13);
                assert!((cycle["amplitude_pa"].as_f64().unwrap() - amplitude).abs() < 1e-13);
            }
        }
    }
    for (kind, damage) in [("signed_von_mises", 0.04), ("max_principal", 0.01)] {
        let p = tensor_input(kind, [2., 0., 0., 0., 0., 0.]);
        let r = Session::default().execute("solve.fatigue.transient", &p).unwrap();
        assert!((r["nodes"][0]["damage"].as_f64().unwrap() - damage).abs() < 1e-13);
    }
}
#[test]
fn nonlinear_tensor_criteria_refuse_ambiguous_or_incomplete_information() {
    let p = tensor_input("signed_von_mises", [0., 0., 0., 1., 0., 0.]);
    assert!(Session::default().execute("solve.fatigue.transient", &p).is_err());
    for kind in ["signed_von_mises", "max_principal"] {
        let p = tensor_input(kind, [2., 0., 0., 0., 0., 0.]);
        for (path, value) in [
            ("/history/nodes/0/samples/1/stress", json!([0, 2, 0, 0, 0, 0])),
            ("/history/nodes/0/samples/1/stress", json!([2])),
            ("/history/stress_location", json!("local axes")),
        ] {
            let mut bad = p.clone();
            *bad.pointer_mut(path).unwrap() = value;
            assert!(Session::default().execute("solve.fatigue.transient", &bad).is_err());
        }
        let mut bad = p.clone();
        bad["weights"] = json!([1, 0, 0, 0, 0, 0]);
        assert!(Session::default().execute("solve.fatigue.transient", &bad).is_err());
        for request in [json!({"criterion":{"type":kind}}), json!({"stress_projection":{"criterion":{"type":kind}}})] {
            let e = Session::default().execute("solve.random", &request).unwrap_err().to_string();
            if request.get("criterion").is_some() {
                assert!(e.contains("has no parameter `criterion`"), "{e}");
            } else {
                assert!(e.contains("PSD"), "{e}");
            }
        }
    }
}
