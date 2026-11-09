use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn setup() -> Session {
    let mut s = Session::default();
    let recipe:Value=serde_json::from_str(r#"[{"command": "model.units", "params": {"system": "m-kg-s-N-Pa"}}, {"command": "mesh.node.add", "params": {"label": 1, "xyz": [0, 0, 0]}}, {"command": "mesh.node.add", "params": {"label": 2, "xyz": [1, 0, 0]}}, {"command": "mesh.element.add", "params": {"label": 1, "kind": "rod2", "nodes": [1, 2]}}, {"command": "model.set.create", "params": {"name": "root", "nodes": [1]}}, {"command": "model.set.create", "params": {"name": "tip", "nodes": [2]}}, {"command": "model.set.create", "params": {"name": "member", "elements": [1]}}, {"command": "model.material.create", "params": {"name": "synthetic", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0.3}, "density": 1}}, {"command": "model.section.create", "params": {"name": "rod", "set": "member", "material": "synthetic", "kind": {"type": "rod", "area": 1}}}, {"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}}, {"command": "model.boundary_condition.create", "params": {"name": "guide", "type": "displacement", "target": {"set": "tip"}, "values": [null, 0, 0, null, null, null]}}, {"command": "model.load.create", "params": {"name": "Aforce", "type": "force", "target": {"set": "tip"}, "vector": [1, 0, 0]}}, {"command": "model.load.create", "params": {"name": "Bforce", "type": "force", "target": {"set": "tip"}, "vector": [2, 0, 0]}}, {"command": "model.step.create", "params": {"name": "A", "procedure": {"type": "static"}, "boundary_conditions": ["clamp", "guide"], "loads": ["Aforce"], "output": ["displacement", "stress"]}}, {"command": "model.step.create", "params": {"name": "B", "procedure": {"type": "static"}, "boundary_conditions": ["clamp", "guide"], "loads": ["Bforce"], "output": ["displacement", "stress"]}}, {"command": "solve.static", "params": {"step": "A", "threads": 1}}, {"command": "solve.static", "params": {"step": "B", "threads": 1}}]"#).unwrap();
    s.run_script(&recipe).unwrap();
    s
}
fn params() -> Value {
    serde_json::from_str(r#"{"channels": [{"step": "A", "coefficients": [-1, 1, -1]}, {"step": "B", "coefficients": [-0.5, 0.5, -0.5]}], "times_seconds": [0, 1, 2], "nodes": [2], "weights": [1, 0, 0, 0, 0, 0], "curve": {"amplitude_ref": 2, "cycles_ref": 100, "exponent": 3, "endurance_limit": 0}, "correction": {"type": "none"}}"#).unwrap()
}
#[test]
fn two_actual_static_channels_sum_before_counting() {
    let mut s = setup();
    let mut p = params();
    let r = s.execute("solve.fatigue.history", &p).unwrap();
    assert!((r["nodes"][0]["damage"].as_f64().unwrap() / 0.01 - 1.0).abs() < 1e-10);
    assert_eq!(r["nodes"][0]["cycle_count"], 1.0);
    p["channels"][1]["coefficients"] = json!([0.5, -0.5, 0.5]);
    let zero = s.execute("solve.fatigue.history", &p).unwrap();
    assert_eq!(zero["nodes"][0]["damage"], 0.0);
    assert_eq!(zero["nodes"][0]["repeated_record_life_seconds"], Value::Null);
}
#[test]
fn bad_time_and_channel_dimensions_are_refused() {
    let mut s = setup();
    for (key, value) in [("times_seconds", json!([0, 0, 2])), ("nodes", json!([2, 2])), ("weights", json!([0, 0, 0, 0, 0, 0]))] {
        let mut p = params();
        p[key] = value;
        assert!(s.execute("solve.fatigue.history", &p).is_err());
    }
    let mut p = params();
    p["channels"][0]["coefficients"] = json!([1, 2]);
    assert!(s.execute("solve.fatigue.history", &p).is_err());
    s.execute(
        "model.boundary_condition.create",
        &json!({"name":"driven","type":"displacement","target":{"set":"root"},"values":[0.1,0,0,null,null,null]}),
    )
    .unwrap();
    s.execute("model.step.edit", &json!({"name":"A","boundary_conditions":["driven","guide"]})).unwrap();
    assert!(s.execute("solve.fatigue.history", &params()).is_err());
}

#[test]
fn full_field_segmented_curve_contours_persistence_and_dependencies() {
    use fieldcraft_engine::results::{Field, Pick, field_values};
    let mut s = setup();
    let mut p = params();
    p["nodes"] = json!("all");
    p["curve"] = json!({"points":[{"cycles":100,"amplitude":4},{"cycles":10000,"amplitude":0.4},{"cycles":10000000,"amplitude":0.04}],"endurance_limit":0,"surface_factor":0.8,"size_factor":0.5,"extrapolation":"end_slopes"});
    let r = s.execute("solve.fatigue.history", &p).unwrap();
    assert_eq!(r["coverage_all"], true);
    // Effective amplitude 2/(.8*.5)=5 Pa, first slope k=2: D=(5/4)^2/100.
    let expected = 0.015625;
    for n in r["nodes"].as_array().unwrap() {
        assert!((n["damage"].as_f64().unwrap() - expected).abs() < 1e-12);
    }
    let c = s.execute("results.contour", &json!({"step":"A","field":"fatigue_damage"})).unwrap();
    assert!((c["contour"]["range"]["max"]["value"].as_f64().unwrap() - expected).abs() < 1e-12);
    s.execute("results.contour", &json!({"step":"A","field":"fatigue_life"})).unwrap();
    assert!(s.execute("results.animate", &json!({"on":true})).is_err());
    let path = std::env::temp_dir().join(format!("fieldcraft-fatigue-{}.fieldcraft", std::process::id()));
    s.execute("file.save", &json!({"path":path.to_str().unwrap()})).unwrap();
    let mut opened = Session::default();
    opened.execute("file.open", &json!({"path":path.to_str().unwrap()})).unwrap();
    let stored = opened.result("A").unwrap();
    let lives = field_values(&opened.model, stored, Field::FatigueLife, Pick::default()).unwrap();
    assert_eq!(lives.len(), 2);
    for v in lives {
        assert!((v - 128.).abs() < 1e-9);
    }
    opened.execute("model.load.edit", &json!({"name":"Bforce","type":"force","target":{"set":"tip"},"vector":[3,0,0]})).unwrap();
    assert!(opened.result("A").is_none(), "second channel must invalidate fatigue on first step");
    let vtu = path.with_extension("vtu");
    s.execute("results.export.vtu", &json!({"step":"A","path":vtu.to_str().unwrap(),"encoding":"ascii"})).unwrap();
    let exported = std::fs::read_to_string(&vtu).unwrap();
    for name in ["fatigue_damage", "fatigue_life", "has_fatigue_damage", "has_fatigue_life"] {
        assert!(exported.contains(name));
    }
    std::fs::remove_file(vtu).unwrap();
    std::fs::remove_file(path).unwrap();
}

#[test]
fn selected_coverage_does_not_fabricate_zero_damage_elsewhere() {
    use fieldcraft_engine::results::{Field, Pick, field_values};
    let mut s = setup();
    s.execute("solve.fatigue.history", &params()).unwrap();
    let values = field_values(&s.model, s.result("A").unwrap(), Field::FatigueDamage, Pick::default()).unwrap();
    assert!(values[0].is_nan());
    assert!((values[1] - 0.01).abs() < 1e-12);
}

#[test]
fn phase_aligned_channels_are_added_as_signed_tensors_not_equivalent_amplitudes() {
    let mut s = setup();
    let mut p = params();
    p["times_seconds"] = json!([0., 0.25, 0.5, 0.75, 1.]);
    p["channels"][0]["coefficients"] = json!([0., 1., 0., -1., 0.]);
    p["channels"][1]["coefficients"] = json!([0.5, 0., -0.5, 0., 0.5]);
    // F/A references are 1 and 2 Pa. Quadrature samples sum to [1,1,-1,-1,1].
    // Plateaus reduce to [1,-1,1]: two half cycles, amplitude 1, mean 0,
    // D=(1/2)^3/100=1/800, repeated-block life=800 seconds.
    for criterion in [json!({"type":"projection","weights":[1,0,0,0,0,0]}), json!({"type":"signed_von_mises"})] {
        p.as_object_mut().unwrap().remove("weights");
        p["criterion"] = criterion;
        let r = s.execute("solve.fatigue.history", &p).unwrap();
        assert!((r["nodes"][0]["damage"].as_f64().unwrap() - 0.00125).abs() < 1e-13);
        assert_eq!(r["nodes"][0]["cycle_count"], 1.);
        assert!((r["nodes"][0]["repeated_record_life_seconds"].as_f64().unwrap() - 800.).abs() < 1e-9);
    }
    // Per-channel clocks/phase shifts are not guessed or silently interpolated.
    p["channels"][1]["times_seconds"] = json!([0.1, 0.35, 0.6, 0.85, 1.1]);
    assert!(s.execute("solve.fatigue.history", &p).is_err());
}
