//! Physical swept layers through native command assembly and result persistence.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
#[test]
fn swept_layers_static_recovery_and_reopen() {
    let cases = [
        (
            "solid-shell-pure-bending",
            include_str!("../../../tests/oracle/solid-shell-pure-bending/recipe.json"),
            include_str!("../../../tests/oracle/solid-shell-pure-bending/expected.json"),
        ),
        (
            "solid-shell-coupled-layers",
            include_str!("../../../tests/oracle/solid-shell-coupled-layers/recipe.json"),
            include_str!("../../../tests/oracle/solid-shell-coupled-layers/expected.json"),
        ),
        (
            "solid-shell-thick-bending",
            include_str!("../../../tests/oracle/solid-shell-thick-bending/recipe.json"),
            include_str!("../../../tests/oracle/solid-shell-thick-bending/expected.json"),
        ),
    ];
    for (name, text, reference) in cases {
        let recipe: Value = serde_json::from_str(text).unwrap();
        let expected: Value = serde_json::from_str(reference).unwrap();
        let mut s = Session::default();
        s.run_script(&recipe).unwrap();
        let probe = recipe.as_array().unwrap().last().unwrap()["params"].clone();
        let before = s.execute("solve.probe", &probe).unwrap();
        let w = before["nodes"][0]["displacement"][2].as_f64().unwrap();
        let exact = expected["quantities"][0]["value"].as_f64().unwrap();
        assert!((w / exact - 1.0).abs() < 1e-7, "{name}: {w} vs {exact}");
        println!("{name}: w={w:.17e}, exact={exact:.17e}, relative={:.8e}", w / exact - 1.0);
        let root = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
        let path = root.join(format!("{name}-{}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut opened = Session::default();
        opened.execute("file.open", &json!({"path":path})).unwrap();
        assert_eq!(before, opened.execute("solve.probe", &probe).unwrap());
        if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_none() {
            std::fs::remove_file(path).unwrap();
        }
        let deck = fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MKgS);
        assert!(deck.is_err());
    }
}

#[test]
fn swept_mass_modes_and_low_frequency_command_limit() {
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/solid-shell-pure-bending/recipe.json")).unwrap();
    let mut s = Session::default();
    for c in recipe.as_array().unwrap() {
        let command = c["command"].as_str().unwrap();
        if !command.starts_with("solve.") {
            s.execute(command, &c["params"]).unwrap();
        }
    }
    let expected = 1.875104068711961_f64.powi(2) / (2.0 * std::f64::consts::PI) * (70e9_f64 * 0.02_f64.powi(2) / (12.0 * 2700.0)).sqrt();
    let step = recipe.as_array().unwrap().iter().find(|c| c["command"] == "model.step.create").unwrap();
    let loads = step["params"]["loads"].clone();
    s.execute("model.step.create", &json!({"name":"Modes","procedure":{"type":"frequency","modes":3},"boundary_conditions":["clamp"]})).unwrap();
    s.execute("model.step.create",&json!({"name":"Low","procedure":{"type":"harmonic","frequencies":[0.01],"method":{"type":"direct"},"damping":{"alpha":2.0}},"boundary_conditions":["clamp"],"loads":loads,"output":["displacement"]})).unwrap();
    s.execute("model.step.create",&json!({"name":"LowModal","procedure":{"type":"harmonic","frequencies":[0.01],"method":{"type":"modal","modes":48},"damping":{"alpha":2.0}},"boundary_conditions":["clamp"],"loads":loads,"output":["displacement"]})).unwrap();
    for mass in ["consistent", "lumped"] {
        let f = s.execute("solve.frequency", &json!({"step":"Modes","mass":mass})).unwrap();
        let hz = f["modes"][0]["frequency"].as_f64().unwrap();
        let total = f["total_mass"].as_f64().unwrap();
        assert!((total / 10.8 - 1.0).abs() < 1e-12);
        assert!((hz / expected - 1.0).abs() < 0.03, "{mass}: {hz} vs {expected}");
        s.execute("solve.harmonic", &json!({"step":"Low","mass":mass})).unwrap();
        let p = s.execute("solve.probe", &json!({"step":"Low","nodes":[15]})).unwrap();
        let w = p["nodes"][0]["samples"][0]["displacement_real"][2].as_f64().unwrap();
        assert!((w / (-5.3571428571428555e-5) - 1.0).abs() < 1e-4, "{p}");
        s.execute("solve.harmonic", &json!({"step":"LowModal","mass":mass})).unwrap();
        let mp = s.execute("solve.probe", &json!({"step":"LowModal","nodes":[15]})).unwrap();
        for field in ["displacement_real", "displacement_imaginary"] {
            let a = p["nodes"][0]["samples"][0][field][2].as_f64().unwrap();
            let b = mp["nodes"][0]["samples"][0][field][2].as_f64().unwrap();
            assert!((a - b).abs() < 1e-7 * a.abs(), "{mass} {field}: {a} vs {b}");
        }
        println!("solid-shell {mass}: totalmass={total}, first_hz={hz}, Euler={expected}, low_w={w}");
        let path = std::env::temp_dir().join(format!("solid-shell-dynamic-{mass}-{}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":path})).unwrap();
        assert_eq!(p, reopened.execute("solve.probe", &json!({"step":"Low","nodes":[15]})).unwrap());
        std::fs::remove_file(path).unwrap();
    }
    s.execute("model.step.create", &json!({"name":"Buckle","procedure":{"type":"buckle","modes":1},"boundary_conditions":["clamp"],"loads":loads}))
        .unwrap();
    let e = s.execute("solve.buckle", &json!({"step":"Buckle"})).unwrap_err().to_string();
    assert!(e.contains("solid-shell buckling"), "{e}");
}
