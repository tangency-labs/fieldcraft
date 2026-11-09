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

#[test]
fn exact_elastic_plastic_limits_and_zero_mean_equivalence() {
    // Analytic power laws at 2N=10000, not a second numerical root solver.
    for (strength, ductility, amplitude, correction, stress_key, stress) in [
        (1e9, 0., 0.00005, "none", "mean_stress_pa", 0.),
        (0., 0.5, 0.005, "none", "mean_stress_pa", 0.),
        (1e9, 0., 0.000025, "morrow", "mean_stress_pa", 0.5e9),
        (1e9, 0., 0.00005, "swt", "maximum_stress_pa", 1e7),
        (1e9, 0.5, 0.00505, "morrow", "mean_stress_pa", 0.),
        (1e9, 0.5, 0.00505, "morrow", "mean_stress_pa", -0.5e9),
        (1e9, 0.5, 0.00505, "swt", "maximum_stress_pa", 1e7),
    ] {
        let mut p = params(correction, amplitude);
        p["curve"]["strength_pa"] = json!(strength);
        p["curve"]["ductility"] = json!(ductility);
        p["cycles"][0][stress_key] = json!(stress);
        let r = Session::default().execute("solve.fatigue.strain_cycles", &p).unwrap();
        assert!((r["cycles"][0]["life_cycles"].as_f64().unwrap() / 5000. - 1.).abs() < 1e-10);
        assert!((r["damage"].as_f64().unwrap() / 0.0005 - 1.).abs() < 1e-10);
    }
    for (strength, ductility, correction) in [(0., 0., "none"), (-1., 0.5, "none"), (0., 0.5, "morrow"), (0., 0.5, "swt")] {
        let mut p = params(correction, 0.005);
        p["curve"]["strength_pa"] = json!(strength);
        p["curve"]["ductility"] = json!(ductility);
        p["cycles"][0]["mean_stress_pa"] = json!(-1.);
        p["cycles"][0]["maximum_stress_pa"] = json!(1.);
        assert!(Session::default().execute("solve.fatigue.strain_cycles", &p).is_err());
    }
}

#[test]
fn public_1020_mixed_terms_match_independent_decimal_reference() {
    // Public 1020 L HR plate row: Union College MER419 Fatigue Review I.
    // E=205 GPa, sigma_f'=895 MPa, epsilon_f'=.41, b=-.12, c=-.51.
    // Independent 60-digit direct-N bisection at total amplitude .005;
    // these are equation references using public coefficients, not measured lives.
    // Source: https://rbb.union.edu/courses/mer419/content/data/Lectures/MER419%20L01%20Fatigue%20Review%20I.pdf
    for (kind, n, damage) in [
        ("none", 5475.357354894572, 0.0004565912027212583),
        ("morrow", 5057.214341758491, 0.0004943432947575446),
        ("swt", 3037.0386417023486, 0.0008231702967726078),
    ] {
        let p = json!({"curve":{"young_pa":205e9,"strength_pa":895e6,"ductility":0.41,"b":-0.12,"c":-0.51},
            "cycles":[{"total_strain_amplitude":0.005,"count":2.5,"mean_stress_pa":100e6,"maximum_stress_pa":400e6}],"correction":kind});
        let r = Session::default().execute("solve.fatigue.strain_cycles", &p).unwrap();
        assert!((r["cycles"][0]["life_cycles"].as_f64().unwrap() / n - 1.).abs() < 1e-10);
        assert!((r["damage"].as_f64().unwrap() / damage - 1.).abs() < 1e-10);
    }
}
