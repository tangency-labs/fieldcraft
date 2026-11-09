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
