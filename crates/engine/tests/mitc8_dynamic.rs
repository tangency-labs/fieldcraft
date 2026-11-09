//! Actual command assembly, analytical oracle recipes, and native result persistence.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

#[test]
fn mitc8_dynamic_and_buckling_results_round_trip() {
    let cases = [
        (
            include_str!("../../../tests/oracle/mitc8-navier-frequency/recipe.json"),
            vec![json!({"step":"consistent","mode":1,"nodes":[33]}), json!({"step":"lumped","mode":1,"nodes":[33]})],
        ),
        (
            include_str!("../../../tests/oracle/mitc8-navier-harmonic/recipe.json"),
            vec![json!({"step":"consistent_direct","nodes":[33]}), json!({"step":"lumped_modal","nodes":[33]})],
        ),
        (include_str!("../../../tests/oracle/mitc8-bryan-buckle/recipe.json"), vec![json!({"step":"Buckle","mode":1,"nodes":[33]})]),
    ];
    for (i, (text, probes)) in cases.into_iter().enumerate() {
        let mut s = Session::default();
        s.run_script(&serde_json::from_str::<Value>(text).unwrap()).unwrap();
        let before: Vec<_> = probes.iter().map(|p| s.execute("solve.probe", p).unwrap()).collect();
        let path = std::env::temp_dir().join(format!("fieldcraft-mitc8-dynamic-{}-{i}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":path})).unwrap();
        for (p, old) in probes.iter().zip(before) {
            assert_eq!(reopened.execute("solve.probe", p).unwrap(), old);
        }
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn curved_free_mitc8_has_six_rigid_modes_through_commands() {
    let mut s = Session::default();
    let original: Value = serde_json::from_str(include_str!("../../../tests/oracle/mitc8-pure-bending/recipe.json")).unwrap();
    let mut commands = Vec::new();
    for c in original.as_array().unwrap() {
        let mut c = c.clone();
        match c["command"].as_str().unwrap() {
            "model.step.create" | "solve.static" | "solve.probe" => continue,
            "mesh.node.add" => {
                let x = c["params"]["xyz"][0].as_f64().unwrap();
                c["params"]["xyz"][2] = json!(0.1 * x * x);
            }
            "model.material.create" => {
                c["params"]["density"] = json!(3.0);
            }
            _ => {}
        }
        commands.push(c);
    }
    s.run_script(&json!(commands)).unwrap();
    s.execute("model.step.create", &json!({"name":"Free","procedure":{"type":"frequency","modes":7}})).unwrap();
    let mut total: Option<f64> = None;
    for mass in ["consistent", "lumped"] {
        let r = s.execute("solve.frequency", &json!({"step":"Free","mass":mass})).unwrap();
        let modes = r["modes"].as_array().unwrap();
        assert_eq!(modes.len(), 7);
        let elastic = modes[6]["frequency"].as_f64().unwrap();
        assert!(elastic > 0.0);
        for m in &modes[..6] {
            assert!(m["frequency"].as_f64().unwrap() < 1e-5 * elastic, "{r}");
        }
        let mass = r["total_mass"].as_f64().unwrap();
        assert!(mass > 0.6);
        if let Some(old) = total {
            assert!((mass - old).abs() < 1e-12);
        }
        total = Some(mass);
    }
}
