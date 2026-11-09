use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn input(values: &[f64]) -> Value {
    json!({"source_id":"synthetic rod element 7","step":"Motion","max_time_step_seconds":1,
        "axial_history":{"step":"Motion","element":7,"location":"constant_rod2_midpoint","stress_measure":"nominal_axial_reference_area","axis":"material_following_axial","reference_material_axis":[1,0,0],"reference_area_m2":2,"stress_unit":"Pa","time_unit":"s","samples":values.iter().enumerate().map(|(i,v)|json!({"time":i,"stress":v})).collect::<Vec<_>>()},
        "curve_stress_measure":"nominal_axial_reference_area","calibration_source":"synthetic nominal S-N and strengths, not measured material",
        "curve":{"amplitude_ref":10,"cycles_ref":1,"exponent":2,"endurance_limit":0},"correction":{"type":"none"}})
}
#[test]
fn nominal_axial_closed_cycles_half_cycles_and_mean_correction() {
    let mut s = Session::default();
    // Outer 0..6 range remains as two halves, nested 2..4 is one full cycle.
    let p = input(&[0., 6., 2., 4., 0.]);
    let r = s.execute("solve.fatigue.transient", &p).unwrap();
    assert_eq!(r["element"], 7);
    assert_eq!(r["persisted"], false);
    let c = r["cycles"].as_array().unwrap();
    assert_eq!(c.len(), 3);
    assert_eq!(c.iter().filter(|c| c["count"] == 1.).count(), 1);
    assert!((r["damage"].as_f64().unwrap() - 0.10).abs() < 1e-13);
    let mut p = input(&[0., 2., 0.]);
    p["correction"] = json!({"type":"goodman","ultimate":4});
    let r = s.execute("solve.fatigue.transient", &p).unwrap();
    assert!((r["damage"].as_f64().unwrap() - (4_f64 / 3. / 10.).powi(2)).abs() < 1e-13);
    for c in r["cycles"].as_array().unwrap() {
        assert_eq!(c["count"], 0.5);
        assert_eq!(c["mean_pa"], 1.);
        assert_eq!(c["amplitude_pa"], 1.);
    }
    assert_eq!(r["stress_measure"], p["curve_stress_measure"]);
    assert_eq!(r["calibration_source"], p["calibration_source"]);
}
#[test]
fn axial_measure_location_time_and_tensor_mixing_are_rejected() {
    let mut s = Session::default();
    for (pointer, value, reason) in [
        ("/axial_history/stress_measure", json!("PK2"), "PK2"),
        ("/curve_stress_measure", json!("Cauchy"), "Cauchy"),
        ("/axial_history/location", json!("global nodal average; shell top surface"), "constant_rod2_midpoint"),
        ("/axial_history/axis", json!("global"), "nonproportional"),
        ("/axial_history/samples/1/time", json!(0), "times"),
        ("/axial_history/stress_unit", json!("MPa"), "Pa"),
        ("/axial_history/reference_area_m2", json!(0), "area"),
        ("/axial_history/reference_material_axis", json!([0, 0, 0]), "axis"),
        ("/calibration_source", json!(""), "calibration_source"),
    ] {
        let mut p = input(&[0., 2., 0.]);
        *p.pointer_mut(pointer).unwrap() = value;
        let e = s.execute("solve.fatigue.transient", &p).unwrap_err().to_string();
        assert!(e.contains(reason), "{e}");
    }
    let mut p = input(&[0., 2., 0.]);
    p["criterion"] = json!({"type":"signed_von_mises"});
    assert!(s.execute("solve.fatigue.transient", &p).unwrap_err().to_string().contains("criterion"));
    let mut p = input(&[0., 2., 0.]);
    p["history"] = json!({});
    assert!(s.execute("solve.fatigue.transient", &p).is_err());
}
