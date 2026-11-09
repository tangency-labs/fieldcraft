use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn params(correction: &str, amplitude: f64) -> Value {
    json!({"curve":{"young_pa":200e9,"strength_pa":1e9,"ductility":0.5,"b":-0.5,"c":-0.5},"cycles":[{"total_strain_amplitude":amplitude,"count":2.5}],"correction":correction})
}
#[test]
fn counted_total_strain_has_closed_form_life_and_miner_damage() {
    for (kind, amplitude) in [("none", 0.00505), ("morrow", 0.005025), ("swt", 0.00505)] {
        let mut p = params(kind, amplitude);
        if kind == "morrow" {
            p["cycles"][0]["mean_stress_pa"] = json!(0.5e9);
        }
        if kind == "swt" {
            p["cycles"][0]["maximum_stress_pa"] = json!(1e7);
        }
        let mut s = Session::default();
        let r = s.execute("solve.fatigue.strain_cycles", &p).unwrap();
        assert!((r["cycles"][0]["life_cycles"].as_f64().unwrap() / 5000.0 - 1.0).abs() < 1e-10);
        assert!((r["damage"].as_f64().unwrap() / 0.0005 - 1.0).abs() < 1e-10);
        s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
        assert_eq!(r, s.execute("solve.fatigue.strain_cycles", &p).unwrap());
    }
}
#[test]
fn counted_strain_requires_corresponding_stress_and_finite_positive_counts() {
    for p in [params("morrow", 0.005), params("swt", 0.005), params("none", 1.0)] {
        assert!(Session::default().execute("solve.fatigue.strain_cycles", &p).is_err());
    }
    for (path, value) in
        [("/cycles/0/count", json!(0)), ("/cycles/0/total_strain_amplitude", json!(-1)), ("/curve/b", json!(0)), ("/correction", json!("goodman"))]
    {
        let mut p = params("none", 0.00505);
        *p.pointer_mut(path).unwrap() = value;
        assert!(Session::default().execute("solve.fatigue.strain_cycles", &p).is_err());
    }
    let r = Session::default().execute("solve.fatigue.strain_cycles", &params("none", 0.0)).unwrap();
    assert!(r["cycles"][0]["life_cycles"].is_null());
    assert_eq!(r["damage"], 0.0);
}
