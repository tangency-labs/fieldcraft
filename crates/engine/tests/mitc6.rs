//! Actual command assembly, analytical oracle recipes, and native result persistence.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn mitc6a_dynamic_and_buckling_results_round_trip() {
    let cases = [
        (include_str!("../../../tests/oracle/mitc6a-pure-bending/recipe.json"), vec![json!({"step":"Load","nodes":[1,8],"elements":[1]})]),
        (include_str!("../../../tests/oracle/mitc6a-coupled-strip/recipe.json"), vec![json!({"step":"Load","nodes":[1,8],"elements":[1]})]),
        (
            include_str!("../../../tests/oracle/mitc6a-navier-frequency/recipe.json"),
            vec![json!({"step":"consistent","mode":1,"nodes":[33]}), json!({"step":"lumped","mode":1,"nodes":[33]})],
        ),
        (
            include_str!("../../../tests/oracle/mitc6a-navier-harmonic-n12/recipe.json"),
            vec![json!({"step":"consistent_direct","nodes":[313]}), json!({"step":"lumped_modal","nodes":[313]})],
        ),
        (include_str!("../../../tests/oracle/mitc6a-bryan-buckle-n8/recipe.json"), vec![json!({"step":"Buckle","mode":1,"nodes":[145]})]),
    ];
    for (i, (text, probes)) in cases.into_iter().enumerate() {
        let mut s = Session::default();
        let replies = s.run_script(&serde_json::from_str::<Value>(text).unwrap()).unwrap();
        let before: Vec<_> = probes.iter().map(|p| s.execute("solve.probe", p).unwrap()).collect();
        let path = std::env::temp_dir().join(format!("fieldcraft-mitc6a-dynamic-{}-{i}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":path})).unwrap();
        for (p, old) in probes.iter().zip(before) {
            assert_eq!(reopened.execute("solve.probe", p).unwrap(), old);
        }
        if let Ok(folder) = std::env::var("FIELDCRAFT_TEST_ARTIFACTS") {
            let folder = std::path::Path::new(&folder);
            std::fs::create_dir_all(folder).unwrap();
            std::fs::copy(&path, folder.join(format!("mitc6a-case-{i}.fieldcraft"))).unwrap();
            std::fs::write(folder.join(format!("mitc6a-case-{i}-replies.json")), serde_json::to_vec_pretty(&replies).unwrap()).unwrap();
            let script: Value = serde_json::from_str(text).unwrap();
            let nodes: Vec<_> =
                script.as_array().unwrap().iter().filter(|c| c["command"] == "mesh.node.add").map(|c| c["params"]["label"].clone()).collect();
            let elements: Vec<_> =
                script.as_array().unwrap().iter().filter(|c| c["command"] == "mesh.element.add").map(|c| c["params"]["label"].clone()).collect();
            if i < 2 {
                let full = s.execute("solve.probe", &json!({"step":"Load","nodes":nodes,"elements":elements})).unwrap();
                std::fs::write(folder.join(format!("mitc6a-case-{i}-full-probe.json")), serde_json::to_vec_pretty(&full).unwrap()).unwrap();
            }
        }
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn curved_free_mitc6a_has_six_rigid_modes_through_commands() {
    let mut s = Session::default();
    let original: Value = serde_json::from_str(include_str!("../../../tests/oracle/mitc6a-pure-bending/recipe.json")).unwrap();
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

#[test]
fn mitc6a_static_round_trip() {
    let mut s = Session::default();
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/mitc6a-coupled-strip/recipe.json")).unwrap();
    s.run_script(&recipe).unwrap();
    let p = json!({"step":"Load","nodes":[1,8],"elements":[1]});
    let before = s.execute("solve.probe", &p).unwrap();
    let path = std::env::temp_dir().join(format!("fieldcraft-mitc6a-static-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut reopened = Session::default();
    reopened.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(reopened.execute("solve.probe", &p).unwrap(), before);
    std::fs::remove_file(path).unwrap();
}
