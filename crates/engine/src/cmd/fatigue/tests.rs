//! The fatigue commands against closed forms: a uniform bar reproduces the S-N curve and the mean
//! stress corrections exactly, a two-block history gives Miner's sum, and a plate with a hole
//! gives the life of its stress concentration factor (Kirsch's Kt = 3).

use serde_json::{Value, json};

use crate::Session;
use crate::results::{Field, Pick, field_values};

/// σ_f′ = 1000 MPa, b = −0.1: an illustrative curve in round numbers (not a real material's).
const SIGMA_F: f64 = 1000.0;
const B: f64 = -0.1;

/// Basquin's life N = ½ (σ_ar/σ_f′)^(1/b).
fn basquin(amplitude: f64) -> f64 {
    0.5 * (amplitude / SIGMA_F).powf(1.0 / B)
}

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * b.abs()
}

/// A bar 10 × 1 × 1 mm of ten hex8 elements along X, held on its three symmetry planes (ν = 0)
/// and pulled by 200 N at its end: a uniform σ_xx = 200 MPa (step `Pull`).
fn bar(s: &mut Session) {
    let label = |i: usize, j: usize, k: usize| 1 + i + 11 * (j + 2 * k);
    let mut script = Vec::new();
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..11 {
                script.push(json!({"command": "mesh.node.add", "params": {"label": label(i, j, k), "xyz": [i, j, k]}}));
            }
        }
    }
    for i in 0..10 {
        let nodes = [
            label(i, 0, 0),
            label(i + 1, 0, 0),
            label(i + 1, 1, 0),
            label(i, 1, 0),
            label(i, 0, 1),
            label(i + 1, 0, 1),
            label(i + 1, 1, 1),
            label(i, 1, 1),
        ];
        script.push(json!({"command": "mesh.element.add", "params": {"label": i + 1, "kind": "hex8", "nodes": nodes}}));
    }
    let plane = |f: &dyn Fn(usize, usize, usize) -> bool| -> Vec<usize> {
        let mut v = Vec::new();
        for k in 0..2 {
            for j in 0..2 {
                for i in 0..11 {
                    if f(i, j, k) {
                        v.push(label(i, j, k));
                    }
                }
            }
        }
        v
    };
    script.extend([
        json!({"command": "model.set.create", "params": {"name": "all", "elements": (1..=10).collect::<Vec<_>>()}}),
        json!({"command": "model.set.create", "params": {"name": "x0", "nodes": plane(&|i, _, _| i == 0)}}),
        json!({"command": "model.set.create", "params": {"name": "y0", "nodes": plane(&|_, j, _| j == 0)}}),
        json!({"command": "model.set.create", "params": {"name": "z0", "nodes": plane(&|_, _, k| k == 0)}}),
        json!({"command": "model.set.create", "params": {"name": "end", "nodes": plane(&|i, _, _| i == 10)}}),
        json!({"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 200000, "poisson": 0}, "density": 7.85e-9}}),
        json!({"command": "model.section.create", "params": {"name": "bar", "set": "all", "material": "Steel", "kind": {"type": "solid"}}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "sx", "type": "symmetry", "target": {"set": "x0"}, "plane": "x"}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "sy", "type": "symmetry", "target": {"set": "y0"}, "plane": "y"}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "sz", "type": "symmetry", "target": {"set": "z0"}, "plane": "z"}}),
        json!({"command": "model.load.create", "params": {"name": "pull", "type": "force", "target": {"set": "end"}, "vector": [200, 0, 0]}}),
        json!({"command": "model.step.create", "params": {"name": "Pull", "procedure": {"type": "static"}, "boundary_conditions": ["sx", "sy", "sz"], "loads": ["pull"]}}),
        json!({"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 2}, "boundary_conditions": ["sx", "sy", "sz"]}}),
        json!({"command": "solve.static", "params": {"step": "Pull"}}),
    ]);
    s.run_script(&Value::Array(script)).unwrap();
}

/// The bar's material with the illustrative S-N data (MPa, the session's display units).
fn give_data(s: &mut Session, endurance: Option<f64>) {
    let mut fatigue = json!({"curve": {"type": "basquin", "coefficient": SIGMA_F, "exponent": B}, "ultimate_strength": 800, "yield_strength": 600});
    if let Some(e) = endurance {
        fatigue["endurance_limit"] = json!(e);
    }
    s.execute("model.material.fatigue", &json!({"name": "Steel", "fatigue": fatigue})).unwrap();
}

fn assess(s: &mut Session, loading: Value, mean_stress: &str) {
    s.execute("model.fatigue.set", &json!({"step": "Pull", "loading": loading, "measure": "signed_von_mises", "mean_stress": mean_stress})).unwrap();
}

fn values(s: &Session, field: Field) -> Vec<f64> {
    let r = s.result("Pull").unwrap();
    field_values(&s.model, r, field, Pick::default()).unwrap()
}

/// A uniform bar at constant amplitude: every node's life is Basquin's at the bar's stress, and
/// the mean stress corrections give their formulas' equivalent amplitudes (σ = 200 MPa, S_u = 800,
/// S_y = 600: Goodman σ_a/(1 − σ_m/S_u), Gerber σ_a/(1 − (σ_m/S_u)²), Soderberg σ_a/(1 − σ_m/S_y),
/// Smith-Watson-Topper √((σ_a + σ_m) σ_a)).
#[test]
fn a_uniform_bar_reproduces_the_basquin_curve_and_the_mean_stress_corrections() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, None);
    // The data is stored in SI and read by no solver: the result stays current.
    assert!(s.result("Pull").is_some(), "fatigue data keeps the result current");
    let stored = s.model.material("Steel").and_then(|m| m.fatigue.clone()).unwrap();
    assert_eq!(stored.curve, fieldcraft_model::SnCurve::Basquin { coefficient: 1e9, exponent: B });
    let sigma = 200.0;
    let cases: [(Value, &str, f64); 7] = [
        (json!({"type": "stress_ratio", "ratio": -1}), "none", sigma),
        (json!({"type": "stress_ratio", "ratio": -1}), "goodman", sigma),
        (json!({"type": "stress_ratio", "ratio": 0}), "goodman", 100.0 / (1.0 - 100.0 / 800.0)),
        (json!({"type": "stress_ratio", "ratio": 0}), "gerber", 100.0 / (1.0 - (100.0f64 / 800.0).powi(2))),
        (json!({"type": "stress_ratio", "ratio": 0}), "soderberg", 100.0 / (1.0 - 100.0 / 600.0)),
        (json!({"type": "constant_amplitude", "max": 1.5, "min": -0.5}), "goodman", 200.0 / (1.0 - 100.0 / 800.0)),
        // σ_a = σ_m = 100, σ_max = 200: √(200 · 100).
        (json!({"type": "stress_ratio", "ratio": 0}), "swt", 20000.0f64.sqrt()),
    ];
    for (loading, correction, equivalent) in cases {
        assess(&mut s, loading.clone(), correction);
        assert!(s.result("Pull").is_some(), "an assessment keeps the result current");
        let n = basquin(equivalent);
        for v in values(&s, Field::FatigueLogLife) {
            assert!((v - n.log10()).abs() < 1e-9, "{loading} {correction}: log life {v}, Basquin {}", n.log10());
        }
        for d in values(&s, Field::FatigueDamage) {
            assert!(close(d, 1.0 / n, 1e-9), "{loading} {correction}: {d}");
        }
        let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
        assert_eq!(r["life_unit"], "cycles");
        assert_eq!(r["nodes"]["assessed"], 44, "{r}");
        assert!(close(r["min_life"]["life"].as_f64().unwrap(), n, 1e-9), "{r}");
        assert!(close(r["min_life"]["equivalent_amplitude"].as_f64().unwrap(), equivalent, 1e-9), "{r}");
        assert_eq!(r["units"]["stress"], "MPa");
    }
    // A path plot along the bar reads the same field.
    let path = s
        .execute("results.path", &json!({"step": "Pull", "field": "fatigue_log_life", "points": [[0.5, 0.5, 0.5], [9.5, 0.5, 0.5]], "samples": 5}))
        .unwrap();
    // The last case assessed above is SWT at R = 0: σ_ar = √(200 · 100).
    let last = basquin(20000.0f64.sqrt()).log10();
    let along = path["values"].as_array().unwrap();
    assert!(along.len() == 5 && along.iter().all(|v| (v.as_f64().unwrap() - last).abs() < 1e-9), "{path}");
    // No endurance limit: no safety factor field.
    assert!(field_values(&s.model, s.result("Pull").unwrap(), Field::FatigueSafetyFactor, Pick::default()).is_err());
}

/// Fatigue data goes in and replies come out in the display units, and the stored data is SI:
/// the same bar in Pa and in GPa gives the life of the MPa session.
#[test]
fn fatigue_data_and_replies_follow_the_display_units() {
    for (system, unit, per_mpa) in [("m-kg-s-N-Pa", "Pa", 1e6), ("mm-kg-ms-kN-GPa", "GPa", 1e-3)] {
        let mut s = Session::default();
        bar(&mut s);
        s.execute("model.units", &json!({"system": system})).unwrap();
        let fatigue = json!({
            "curve": {"type": "basquin", "coefficient": SIGMA_F * per_mpa, "exponent": B},
            "endurance_limit": 150.0 * per_mpa, "ultimate_strength": 800.0 * per_mpa, "yield_strength": 600.0 * per_mpa,
        });
        s.execute("model.material.fatigue", &json!({"name": "Steel", "fatigue": fatigue})).unwrap();
        let stored = s.model.material("Steel").and_then(|m| m.fatigue.clone()).unwrap();
        assert!(close(stored.endurance_limit.unwrap(), 1.5e8, 1e-12) && close(stored.ultimate_strength.unwrap(), 8e8, 1e-12), "{system}: {stored:?}");
        assert_eq!(stored.curve, fieldcraft_model::SnCurve::Basquin { coefficient: 1e9, exponent: B }, "{system}");
        assess(&mut s, json!({"type": "stress_ratio", "ratio": -1}), "swt");
        for d in values(&s, Field::FatigueDamage) {
            assert!(close(d, 1.0 / basquin(200.0), 1e-9), "{system}: {d}");
        }
        let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
        assert_eq!(r["units"]["stress"], unit, "{r}");
        for key in ["stress_amplitude", "equivalent_amplitude"] {
            assert!(close(r["min_life"][key].as_f64().unwrap(), 200.0 * per_mpa, 1e-9), "{system} {key}: {r}");
        }
        assert!(close(r["min_life"]["safety_factor"].as_f64().unwrap(), 0.75, 1e-9), "{system}: {r}");
    }
}

/// The safety factor against the endurance limit, and infinite life below it: S_e = 150 MPa,
/// fully reversed σ_a = 200 (n = 0.75, finite life) and σ_a = 100 (n = 1.5, infinite life, no
/// damage, drawn in the top band and left out of the range).
#[test]
fn the_endurance_limit_gives_infinite_life_and_the_safety_factor() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, Some(150.0));
    assess(&mut s, json!({"type": "stress_ratio", "ratio": -1}), "goodman");
    for v in values(&s, Field::FatigueSafetyFactor) {
        assert!(close(v, 0.75, 1e-9), "{v}");
    }
    assert!(values(&s, Field::FatigueLogLife).iter().all(|v| v.is_finite()));
    assess(&mut s, json!({"type": "stress_ratio", "ratio": -1, "scale": 0.5}), "goodman");
    assert!(values(&s, Field::FatigueSafetyFactor).iter().all(|v| close(*v, 1.5, 1e-9)));
    assert!(values(&s, Field::FatigueLogLife).iter().all(|v| *v == f64::INFINITY));
    assert!(values(&s, Field::FatigueDamage).iter().all(|v| *v == 0.0));
    let r = s.execute("results.fatigue", &json!({})).unwrap();
    assert_eq!((r["min_life"].clone(), r["nodes"]["infinite_life"].clone()), (Value::Null, json!(44)), "{r}");
    assert!(close(r["min_safety_factor"]["value"].as_f64().unwrap(), 1.5, 1e-9), "{r}");
    let c = s.execute("results.contour", &json!({"step": "Pull", "field": "fatigue_safety_factor"})).unwrap();
    assert_eq!(c["contour"]["field"], "fatigue_safety_factor", "{c}");
    assert_eq!(c["contour"]["unit"], Value::Null);
    let life = s.execute("results.contour", &json!({"field": "fatigue_log_life"})).unwrap();
    assert_eq!(life["contour"]["range"], Value::Null, "every value infinite: no finite range: {life}");
}

/// Miner's sum of a two-block history, by hand: three cycles at ±1 and five at ±0.5 times the
/// bar's 200 MPa, repeated: damage 3/N(200) + 5/N(100) per block, life its inverse in blocks.
#[test]
fn a_two_block_history_gives_miners_sum() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, None);
    let mut factors = Vec::new();
    for _ in 0..3 {
        factors.extend([1.0, -1.0]);
    }
    for _ in 0..5 {
        factors.extend([0.5, -0.5]);
    }
    assess(&mut s, json!({"type": "history", "factors": factors}), "goodman");
    let hand = 3.0 / basquin(200.0) + 5.0 / basquin(100.0);
    for d in values(&s, Field::FatigueDamage) {
        assert!(close(d, hand, 1e-9), "{d} against {hand}");
    }
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["life_unit"], "blocks");
    assert_eq!(r["cycles"], json!([{"range": 2.0, "mean": 0.0, "count": 3.0}, {"range": 1.0, "mean": 0.0, "count": 5.0}]), "{r}");
    assert_eq!(r["cycles_per_block"], 8.0);
    assert!(close(r["min_life"]["life"].as_f64().unwrap(), 1.0 / hand, 1e-9), "{r}");
    // Counted once from rest instead: 17 points, 16 ranges, half a cycle each in all (ASTM
    // E1049's count of a history that does not repeat).
    let mut once = factors.clone();
    once.insert(0, 0.0);
    assess(&mut s, json!({"type": "history", "factors": once, "counting": "once"}), "goodman");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["cycles_per_block"].as_f64(), Some(8.0), "{r}");
}

/// A plate with a small circular hole in tension, 2D plane stress (the oracle case
/// kirsch-plate-plane-stress's mesh: W = 60a, a quarter by symmetry). Kirsch's Kt = 3 puts
/// σ_xx = 3S at the top of the hole; at S = 100 MPa fully reversed the life there follows from
/// the S-N curve at 300 MPa (Kirsch, Z. VDI 42 (1898) 797–807; Timoshenko & Goodier §35). The FE
/// Kt is within the oracle case's 2 %, which moves log₁₀ N by at most |1/b| log₁₀ 1.02 = 0.086.
#[test]
fn the_life_at_a_hole_follows_its_stress_concentration() {
    let mut s = Session::default();
    let r45 = std::f64::consts::FRAC_1_SQRT_2;
    let script = json!([
        {"command": "mesh.block", "params": {"body": "Plate", "kind": "quad8", "corners": [[1.0, 0, 0], [30.0, 0, 0], [30.0, 30.0, 0], [r45, r45, 0]], "divisions": [40, 16], "bias": [40, 1], "edges": [{"corners": [1, 4], "arc": {"center": [0, 0, 0]}}, {"corners": [1, 2], "name": "y0"}, {"corners": [2, 3], "name": "load"}]}},
        {"command": "mesh.block", "params": {"body": "Plate", "kind": "quad8", "corners": [[r45, r45, 0], [30.0, 30.0, 0], [0, 30.0, 0], [0, 1.0, 0]], "divisions": [40, 16], "bias": [40, 1], "edges": [{"corners": [1, 4], "arc": {"center": [0, 0, 0]}}, {"corners": [3, 4], "name": "x0"}]}},
        {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Plate", "face": "area"}]}},
        {"command": "model.set.create", "params": {"name": "x0", "geometry": [{"body": "Plate", "edge": "x0"}]}},
        {"command": "model.set.create", "params": {"name": "y0", "geometry": [{"body": "Plate", "edge": "y0"}]}},
        {"command": "model.set.create", "params": {"name": "load", "geometry": [{"body": "Plate", "edge": "load"}]}},
        {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 200000.0, "poisson": 0.3},
            "fatigue": {"curve": {"type": "basquin", "coefficient": SIGMA_F, "exponent": B}, "ultimate_strength": 800}}},
        {"command": "model.section.create", "params": {"name": "plate", "set": "all", "material": "Steel", "kind": {"type": "plane_stress", "thickness": 0.1}}},
        {"command": "model.boundary_condition.create", "params": {"name": "sx", "type": "symmetry", "target": {"set": "x0"}, "plane": "x"}},
        {"command": "model.boundary_condition.create", "params": {"name": "sy", "type": "symmetry", "target": {"set": "y0"}, "plane": "y"}},
        {"command": "model.load.create", "params": {"name": "S", "type": "traction", "set": "load", "vector": [100, 0, 0]}},
        {"command": "model.step.create", "params": {"name": "Tension", "procedure": {"type": "static"}, "boundary_conditions": ["sx", "sy"], "loads": ["S"]}},
        {"command": "solve.static", "params": {}},
        {"command": "model.fatigue.set", "params": {"step": "Tension", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "max_abs_principal", "mean_stress": "goodman"}},
    ]);
    s.run_script(&script).unwrap();
    let r = s.execute("results.fatigue", &json!({"step": "Tension"})).unwrap();
    // The least life is at the top of the hole, (0, 1).
    let xyz = r["min_life"]["xyz"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>();
    assert!(xyz[0].abs() < 1e-9 && (xyz[1] - 1.0).abs() < 1e-9, "{r}");
    let amplitude = r["min_life"]["stress_amplitude"].as_f64().unwrap();
    let kt = amplitude / 100.0;
    assert!((kt - 3.0).abs() < 0.02 * 3.0, "Kt {kt}");
    // The node's life is Basquin's at its own stress exactly, and Kirsch's within the Kt band.
    let life = r["min_life"]["life"].as_f64().unwrap();
    assert!(close(life, basquin(amplitude), 1e-12), "{life}");
    let theory = basquin(300.0);
    let log_error = (life.log10() - theory.log10()).abs();
    assert!(log_error < 10.0 * 1.02f64.log10(), "log₁₀ N {} against Kirsch's {} ({log_error})", life.log10(), theory.log10());
    // Far from the hole the stress is the remote one and the life Basquin's at S.
    let far = s.execute("solve.probe", &json!({"step": "Tension", "at": [[30, 15, 0]]})).unwrap();
    let remote = far["nodes"][0]["fatigue"]["stress_amplitude"].as_f64().unwrap();
    assert!((remote - 100.0).abs() < 2.0, "{far}");
    assert_eq!(far["units"]["stress_amplitude"], "MPa", "{far}");
    assert!(far["nodes"][0]["fatigue"]["life"].as_f64().unwrap() > 1e3 * life);
    // The contour and the path plot read the same field.
    let c = s.execute("results.contour", &json!({"step": "Tension", "field": "fatigue_log_life"})).unwrap();
    assert!((c["contour"]["range"]["min"]["value"].as_f64().unwrap() - life.log10()).abs() < 1e-9, "{c}");
}

/// The fields follow the saved definitions: a file saved with the result opens with the same
/// fatigue fields, recomputed (the results encoding is unchanged).
#[test]
fn fatigue_fields_survive_save_and_open() {
    let d = std::env::temp_dir().join(format!("fieldcraft-fatigue-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, Some(150.0));
    assess(&mut s, json!({"type": "history", "factors": [0, 1, -0.5, 0.8, -1, 0]}), "gerber");
    let before = values(&s, Field::FatigueDamage);
    let path = d.join("bar.fieldcraft");
    s.execute("file.save_as", &json!({"path": path})).unwrap();
    let mut t = Session::default();
    t.execute("file.open", &json!({"path": path})).unwrap();
    assert_eq!(t.model, s.model);
    assert!(t.result("Pull").is_some(), "the saved result is current");
    assert_eq!(values(&t, Field::FatigueDamage), before);
    // Undo takes an assessment back like any edit.
    t.execute("model.fatigue.delete", &json!({"step": "Pull"})).unwrap();
    assert!(field_values(&t.model, t.result("Pull").unwrap(), Field::FatigueDamage, Pick::default()).is_err());
    t.execute("edit.undo", &json!({})).unwrap();
    assert_eq!(values(&t, Field::FatigueDamage), before);
    let _ = std::fs::remove_dir_all(&d);
}

/// Hostile and unsupported input is refused by name, leaving the model as it was.
#[test]
fn hostile_fatigue_input_is_refused() {
    let mut s = Session::default();
    bar(&mut s);
    let e = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap_err().to_string();
    assert!(e.contains("model.fatigue.set"), "{e}");
    let e = s.execute("results.contour", &json!({"step": "Pull", "field": "fatigue_damage"})).unwrap_err().to_string();
    assert!(e.contains("model.fatigue.set"), "{e}");
    // A static step that is defined but not solved has no result to assess.
    s.execute(
        "model.step.create",
        &json!({"name": "Unsolved", "procedure": {"type": "static"}, "boundary_conditions": ["sx", "sy", "sz"], "loads": ["pull"]}),
    )
    .unwrap();
    s.execute(
        "model.fatigue.set",
        &json!({"step": "Unsolved", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "swt"}),
    )
    .unwrap();
    assert!(s.execute("results.fatigue", &json!({"step": "Unsolved"})).is_err());
    s.execute("model.fatigue.delete", &json!({"step": "Unsolved"})).unwrap();
    assess(&mut s, json!({"type": "stress_ratio", "ratio": -1}), "soderberg");
    let e = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap_err().to_string();
    assert!(e.contains("fatigue data") && e.contains("model.material.fatigue"), "{e}");
    s.execute("model.material.fatigue", &json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": 900, "exponent": -0.09}}}))
        .unwrap();
    let e = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap_err().to_string();
    assert!(e.contains("Steel") && e.contains("yield strength"), "{e}");
    // The probe says why it has no fatigue, and keeps the rest.
    let p = s.execute("solve.probe", &json!({"step": "Pull", "nodes": [1]})).unwrap();
    assert!(p["fatigue_unavailable"].as_str().unwrap().contains("yield strength") && p["nodes"][0]["stress"].is_array(), "{p}");
    let before = s.model.clone();
    let assessment = |loading: Value| json!({"step": "Pull", "loading": loading, "measure": "signed_von_mises", "mean_stress": "none"});
    for bad in [
        json!({"step": "Modes", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "none"}),
        json!({"step": "Nope", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "none"}),
        json!({"loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "none"}),
        json!({"step": "Pull", "loading": {"type": "stress_ratio", "ratio": -1}, "mean_stress": "none"}),
        json!({"step": "Pull", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "tresca", "mean_stress": "none"}),
        json!({"step": "Pull", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "morrow"}),
        json!({"step": "Pull", "loading": {"type": "constant_amplitude", "max": "NaN", "min": 0}, "measure": "signed_von_mises", "mean_stress": "swt"}),
        json!({"step": "Pull", "loading": {"type": "history", "factors": []}, "measure": "signed_von_mises", "mean_stress": "swt"}),
        json!({"step": "Pull", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "none", "stress_factor": 0}),
        json!({"step": "Pull", "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "none", "stress_factor": "1e400"}),
        json!({"step": 7, "loading": {"type": "stress_ratio", "ratio": -1}, "measure": "signed_von_mises", "mean_stress": "none"}),
        assessment(json!({"type": "stress_ratio", "ratio": 3})),
        assessment(json!({"type": "stress_ratio", "ratio": "half"})),
        assessment(json!({"type": "constant_amplitude", "max": -1, "min": 1})),
        assessment(json!({"type": "constant_amplitude", "max": 1e300, "min": 0})),
        assessment(json!({"type": "history", "factors": [1]})),
        assessment(json!({"type": "history", "factors": vec![0.5; fieldcraft_model::fatigue::MAX_HISTORY + 1]})),
        assessment(json!({"type": "history", "factors": [0, "x"]})),
        assessment(json!({"type": "spectrum", "levels": [1]})),
        assessment(json!(null)),
    ] {
        assert!(s.execute("model.fatigue.set", &bad).is_err(), "{bad}");
    }
    for bad in [
        json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": -1, "exponent": -0.1}}}),
        json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": 900, "exponent": 0.1}}}),
        json!({"name": "Steel", "fatigue": {"curve": {"type": "table", "points": [{"cycles": 1e6, "amplitude": 300}, {"cycles": 1e3, "amplitude": 500}]}}}),
        json!({"name": "Steel", "fatigue": {"curve": {"type": "table", "points": vec![json!({"cycles": 1e3, "amplitude": 500}); 65]}}}),
        json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": 900, "exponent": -0.1}, "endurance_limit": 1000}}),
        json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": 900, "exponent": -0.1}, "yield_strength": 900, "ultimate_strength": 800}}),
        json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": "9e999", "exponent": -0.1}}}),
        json!({"name": "Nope", "fatigue": null}),
        json!({"name": "Steel"}),
    ] {
        assert!(s.execute("model.material.fatigue", &bad).is_err(), "{bad}");
    }
    assert!(s.execute("model.fatigue.delete", &json!({"step": "Modes"})).is_err());
    assert!(s.execute("results.fatigue", &json!({"step": "Modes"})).is_err(), "a step with no result");
    assert_eq!(s.model, before, "refusals change nothing");
    // Not a linear static result: a frequency step's fields have no fatigue.
    s.execute("solve.frequency", &json!({"step": "Modes"})).unwrap();
    assert!(s.execute("results.fatigue", &json!({"step": "Modes"})).is_err());
    // Removing the data, then the assessment.
    s.execute("model.material.fatigue", &json!({"name": "Steel", "fatigue": null})).unwrap();
    assert!(s.model.material("Steel").unwrap().fatigue.is_none());
    s.execute("model.fatigue.delete", &json!({"step": "Pull"})).unwrap();
    assert!(s.model.fatigue_assessments().is_empty());
    assert!(s.execute("model.fatigue.delete", &json!({"step": "Pull"})).is_err());
}

/// The longest history runs (a repeated block of 100 000 factors with tens of thousands of
/// distinct cycles over the bar's 44 nodes) and its reply lists a bounded number of cycles.
#[test]
fn the_longest_history_runs_with_a_bounded_reply() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, None);
    let factors: Vec<f64> = (0..fieldcraft_model::fatigue::MAX_HISTORY).map(|i| if i % 2 == 0 { 1.0 + i as f64 * 1e-6 } else { -1.0 }).collect();
    assess(&mut s, json!({"type": "history", "factors": factors}), "none");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert!(r["distinct_cycles"].as_u64().unwrap() > 1000, "{}", r["distinct_cycles"]);
    assert!(r["cycles"].as_array().unwrap().len() <= 200);
    assert!(r["min_life"]["life"].as_f64().unwrap() > 0.0);
}

/// A step that has a fatigue assessment cannot be edited into one that is not linear static
/// (the saved file would not open again), and nothing changes when it is refused. A model that
/// was saved opens again after every accepted command of the sequence.
#[test]
fn an_assessed_step_cannot_become_a_step_fatigue_cannot_read() {
    let d = std::env::temp_dir().join(format!("fieldcraft-fatigue-edit-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    let count = std::cell::Cell::new(0);
    let reopens = |s: &mut Session, what: &str| {
        count.set(count.get() + 1);
        let path = d.join(format!("m{}.fieldcraft", count.get()));
        s.execute("file.save_as", &json!({"path": path})).unwrap();
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": path})).unwrap_or_else(|e| panic!("after {what} the saved file does not open: {e}"));
        assert_eq!(t.model, s.model, "{what}");
    };
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, None);
    reopens(&mut s, "the material data");
    assess(&mut s, json!({"type": "stress_ratio", "ratio": -1}), "goodman");
    reopens(&mut s, "the assessment");
    let before = s.model.clone();
    for procedure in [json!({"type": "nonlinear_static"}), json!({"type": "frequency", "modes": 2})] {
        let e = s.execute("model.step.edit", &json!({"name": "Pull", "procedure": procedure})).unwrap_err().to_string();
        assert!(e.contains("fatigue assessment") && e.contains("model.fatigue.delete"), "{e}");
        assert_eq!(s.model, before, "a refusal changes nothing");
    }
    // An edit that keeps the procedure linear static is fine; deleting the assessment first
    // allows the change; deleting the step takes its assessment along.
    s.execute("model.step.edit", &json!({"name": "Pull", "procedure": {"type": "static"}})).unwrap();
    reopens(&mut s, "a static step edit");
    s.execute("model.fatigue.delete", &json!({"step": "Pull"})).unwrap();
    s.execute("model.step.edit", &json!({"name": "Pull", "procedure": {"type": "nonlinear_static"}})).unwrap();
    reopens(&mut s, "the procedure change after deleting the assessment");
    let _ = std::fs::remove_dir_all(&d);
}

/// A load with equal maximum and minimum has no cycle but can overload the material: a mean of
/// 1000 MPa against an ultimate strength of 800 is a static failure (life 0, safety factor
/// below 1, counted in the reply), and without an ultimate strength the reply says that no
/// static check was made.
#[test]
fn a_static_overload_is_reported_not_called_safe() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, Some(150.0));
    assess(&mut s, json!({"type": "constant_amplitude", "max": 5, "min": 5}), "goodman");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["nodes"]["static_failures"], 44, "{r}");
    assert_eq!(r["nodes"]["infinite_life"], 0, "{r}");
    assert_eq!(r["min_life"]["life"], 0.0, "{r}");
    assert!(close(r["min_safety_factor"]["value"].as_f64().unwrap(), 0.8, 1e-9), "{r}");
    assert_eq!(r["static_check"], true);
    assert!(values(&s, Field::FatigueLogLife).iter().all(|v| v.is_finite()), "a static failure is drawn at the shortest life");
    // A cycle that peaks past the strength (σ_a = 1000, σ_m = 100) is one, whatever the correction.
    assess(&mut s, json!({"type": "constant_amplitude", "max": 5.5, "min": -4.5}), "none");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["nodes"]["static_failures"], 44, "{r}");
    // No ultimate strength: no static check, and the reply says so.
    s.execute("model.material.fatigue", &json!({"name": "Steel", "fatigue": {"curve": {"type": "basquin", "coefficient": SIGMA_F, "exponent": B}}}))
        .unwrap();
    assess(&mut s, json!({"type": "constant_amplitude", "max": 5, "min": 5}), "swt");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["static_check"], false, "{r}");
    assert!(r["warnings"].to_string().contains("no static check"), "{r}");
}

/// The stress factor multiplies the amplitude and the mean alike: R = 0 at 200 MPa with factor 2
/// reads σ_a = σ_m = 200, Goodman σ_ar = 200 / (1 − 200/800) (amplitude alone gave 228.6).
#[test]
fn the_stress_factor_scales_the_mean_too() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, None);
    s.execute(
        "model.fatigue.set",
        &json!({"step": "Pull", "loading": {"type": "stress_ratio", "ratio": 0}, "measure": "signed_von_mises", "mean_stress": "goodman", "stress_factor": 2}),
    )
    .unwrap();
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert!(close(r["min_life"]["mean_stress"].as_f64().unwrap(), 200.0, 1e-9), "{r}");
    assert!(close(r["min_life"]["equivalent_amplitude"].as_f64().unwrap(), 200.0 / 0.75, 1e-9), "{r}");
}

/// One history, both countings: repeated gives a life in blocks (Basquin's life at 200 MPa for
/// [0, 1, -1, 0]); once gives its damage and no life (no life field, null life in the reply).
#[test]
fn a_history_counted_once_has_damage_and_no_life() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, None);
    let history = json!([0, 1, -1, 0]);
    assess(&mut s, json!({"type": "history", "factors": history, "counting": "repeated"}), "none");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["life_unit"], "blocks");
    assert!(close(r["min_life"]["life"].as_f64().unwrap(), basquin(200.0), 1e-9), "{r}");
    let repeated_damage = r["min_life"]["damage"].as_f64().unwrap();
    assert!(close(repeated_damage, 1.0 / basquin(200.0), 1e-9), "{r}");
    assess(&mut s, json!({"type": "history", "factors": history, "counting": "once"}), "none");
    let r = s.execute("results.fatigue", &json!({"step": "Pull"})).unwrap();
    assert_eq!(r["life_unit"], "none", "{r}");
    assert!(r["life_note"].as_str().unwrap().contains("no life"), "{r}");
    assert!(r["min_life"]["life"].is_null() && r["min_life"]["log10_life"].is_null(), "{r}");
    // Counted once: half a cycle of 200, and two halves of range 1 (amplitude 100) at different means.
    let once_damage = 0.5 / basquin(200.0) + 1.0 / basquin(100.0);
    assert!(close(r["min_life"]["damage"].as_f64().unwrap(), once_damage, 1e-9), "{r}");
    assert!(!r["fields"].to_string().contains("fatigue_log_life"), "{r}");
    let e = s.execute("results.contour", &json!({"step": "Pull", "field": "fatigue_log_life"})).unwrap_err().to_string();
    assert!(e.contains("counted once"), "{e}");
    assert!(field_values(&s.model, s.result("Pull").unwrap(), Field::FatigueLogLife, Pick::default()).is_err());
    assert!(values(&s, Field::FatigueDamage).iter().all(|d| close(*d, once_damage, 1e-9)));
    let p = s.execute("solve.probe", &json!({"step": "Pull", "nodes": [1]})).unwrap();
    assert!(p["nodes"][0]["fatigue"]["life"].is_null(), "{p}");
}

/// A path plot through nodes of infinite life marks them: the values are null (an infinite life
/// has no interpolated number) and `infinite_life` is true there, so it is not read as outside
/// the part (null in both `values` and `infinite_life`).
#[test]
fn a_path_marks_where_the_life_is_infinite() {
    let mut s = Session::default();
    bar(&mut s);
    give_data(&mut s, Some(150.0));
    assess(&mut s, json!({"type": "stress_ratio", "ratio": -1, "scale": 0.5}), "goodman");
    let path = json!({"step": "Pull", "field": "fatigue_log_life", "points": [[0.5, 0.5, 0.5], [9.5, 0.5, 0.5], [12.0, 0.5, 0.5]], "samples": 7});
    let r = s.execute("results.path", &path).unwrap();
    let values = r["values"].as_array().unwrap();
    let marks = r["infinite_life"].as_array().unwrap();
    assert!(values[0].is_null() && marks[0] == json!(true), "{r}");
    assert!(values.last().unwrap().is_null() && marks.last().unwrap().is_null(), "outside the part is null in both: {r}");
    // A finite life is a number and is not marked.
    assess(&mut s, json!({"type": "stress_ratio", "ratio": -1}), "goodman");
    let r = s.execute("results.path", &path).unwrap();
    assert!(r["values"][0].is_number() && r["infinite_life"][0] == json!(false), "{r}");
}
