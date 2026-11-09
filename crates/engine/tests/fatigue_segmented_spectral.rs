use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn setup() -> (Session, Value) {
    let mut session = Session::default();
    let recipe:Value=serde_json::from_str(r#"[{"command": "model.units", "params": {"system": "m-kg-s-N-Pa"}}, {"command": "mesh.node.add", "params": {"label": 1, "xyz": [0, 0, 0]}}, {"command": "mesh.node.add", "params": {"label": 2, "xyz": [1, 0, 0]}}, {"command": "mesh.element.add", "params": {"label": 1, "kind": "rod2", "nodes": [1, 2]}}, {"command": "model.set.create", "params": {"name": "root", "nodes": [1]}}, {"command": "model.set.create", "params": {"name": "tip", "nodes": [2]}}, {"command": "model.set.create", "params": {"name": "member", "elements": [1]}}, {"command": "model.material.create", "params": {"name": "synthetic", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0.3}, "density": 1}}, {"command": "model.section.create", "params": {"name": "rod", "set": "member", "material": "synthetic", "kind": {"type": "rod", "area": 1}}}, {"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}}, {"command": "model.boundary_condition.create", "params": {"name": "guide", "type": "displacement", "target": {"set": "tip"}, "values": [null, 0, 0, null, null, null]}}, {"command": "model.load.create", "params": {"name": "unit", "type": "force", "target": {"set": "tip"}, "vector": [1, 0, 0]}}, {"command": "model.step.create", "params": {"name": "Sweep", "procedure": {"type": "harmonic", "frequencies": [1, 1.01, 1.02], "method": {"type": "direct"}}, "boundary_conditions": ["clamp", "guide"], "loads": ["unit"], "output": ["displacement", "stress"]}}, {"command": "solve.harmonic", "params": {"step": "Sweep", "threads": 1}}]"#).unwrap();
    session.run_script(&recipe).unwrap();
    let p=serde_json::from_str(r#"{"step": "Sweep", "nodes": [2], "input_psd": [[1, 4], [1.02, 4]], "stress_projection": {"weights": [1, 0, 0, 0, 0, 0], "input_model": "common_scalar_unit_pattern", "fatigue": {"curve": {"amplitude_ref": 2, "cycles_ref": 100, "exponent": 2, "endurance_limit": 0}, "duration_seconds": 10, "mean_stress_pa": 0, "correction": {"type": "none"}, "stationary_gaussian": true, "method": "narrow_band"}}}"#).unwrap();
    (session, p)
}
#[test]
fn spectral_command_accepts_table_and_rejects_tb() {
    let (mut s, mut p) = setup();
    let f = &mut p["stress_projection"]["fatigue"];
    f["curve"] = json!({"points":[{"cycles":100,"amplitude":10},{"cycles":10000,"amplitude":1},{"cycles":10000000,"amplitude":0.1}],"endurance_limit":0,"surface_factor":0.8,"size_factor":0.5,"extrapolation":"end_slopes"});
    for method in ["narrow_band", "dirlik"] {
        p["stress_projection"]["fatigue"]["method"] = json!(method);
        let r = s.execute("solve.random", &p).unwrap();
        assert!(r["projected_stress"][0]["fatigue"]["damage"].as_f64().unwrap() > 0.);
    }
    p["stress_projection"]["fatigue"]["method"] = json!("benasciutti_tovo");
    assert!(s.execute("solve.random", &p).unwrap_err().to_string().contains("benasciutti_tovo"));
}
