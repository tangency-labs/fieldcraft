use fieldcraft_engine::Session;
use serde_json::{Value, json};
#[test]
fn axisymmetric_ring_static_modes_harmonic_and_native_reopen() {
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/axishell-cylinder/recipe.json")).unwrap();
    let mut s = Session::default();
    s.run_script(&recipe).unwrap();
    let probe = json!({"step":"Static","nodes":[5]});
    let stat = s.execute("solve.probe", &probe).unwrap();
    let ur = stat["nodes"][0]["displacement"][0].as_f64().unwrap();
    assert!((ur / 0.0001 - 1.0).abs() < 1e-8, "{stat}");
    let loads: Vec<_> = (0..9).map(|i| format!("F{i}")).collect();
    s.execute("model.step.create", &json!({"name":"Modes","procedure":{"type":"frequency","modes":1},"boundary_conditions":["hold"]})).unwrap();
    for (name, method) in [("Direct", json!({"type":"direct"})), ("Modal", json!({"type":"modal","modes":1}))] {
        s.execute("model.step.create",&json!({"name":name,"procedure":{"type":"harmonic","frequencies":[30],"method":method,"damping":{"alpha":2}},"boundary_conditions":["hold"],"loads":loads,"output":["displacement","stress"]})).unwrap();
    }
    for mass in ["consistent", "lumped"] {
        let freq = s.execute("solve.frequency", &json!({"step":"Modes","mass":mass})).unwrap();
        let hz = freq["modes"][0]["frequency"].as_f64().unwrap();
        let expected = (1e8_f64 / 1000.0).sqrt() / (2.0 * std::f64::consts::PI);
        let total = freq["total_mass"].as_f64().unwrap();
        assert!((total / (20.0 * std::f64::consts::PI) - 1.0).abs() < 1e-10, "{freq}");
        assert!((hz / expected - 1.0).abs() < 1e-7, "{freq}");
        for name in ["Direct", "Modal"] {
            s.execute("solve.harmonic", &json!({"step":name,"mass":mass})).unwrap();
            let p = json!({"step":name,"nodes":[5]});
            let h = s.execute("solve.probe", &p).unwrap();
            let omega = 60.0 * std::f64::consts::PI;
            let a = 1e5 - omega * omega;
            let b = 2.0 * omega;
            let er = 10.0 * a / (a * a + b * b);
            let ei = -10.0 * b / (a * a + b * b);
            let sample = &h["nodes"][0]["samples"][0];
            for (field, expected) in [("displacement_real", er), ("displacement_imaginary", ei)] {
                let got = sample[field][0].as_f64().unwrap();
                assert!((got / expected - 1.0).abs() < 1e-7, "{field}: {got} vs {expected}");
            }
            println!("axishell {mass} {name}: static={ur:.17e}, hz={hz:.17e}, mass={total:.17e}, response={sample}");
            let root = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
            let path = root.join(format!("axishell-{mass}-{name}-{}.fieldcraft", std::process::id()));
            s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
            let mut reopened = Session::default();
            reopened.execute("file.open", &json!({"path":path})).unwrap();
            assert_eq!(h, reopened.execute("solve.probe", &p).unwrap());
            assert_eq!(stat, reopened.execute("solve.probe", &probe).unwrap());
            if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_none() {
                std::fs::remove_file(path).unwrap();
            }
        }
    }
    s.execute("model.step.create", &json!({"name":"Buckle","procedure":{"type":"buckle","modes":1},"boundary_conditions":["hold"],"loads":loads}))
        .unwrap();
    let err = s.execute("solve.buckle", &json!({"step":"Buckle"})).unwrap_err().to_string();
    assert!(err.contains("axisymmetric"), "{err}");
    let deck = fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MKgS).unwrap_err().to_string();
    assert!(deck.contains("SAX1"), "{deck}");
}

#[test]
fn annular_laminate_bending_patch_through_commands() {
    let mut s = Session::default();
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/axishell-annular-laminate/recipe.json")).unwrap();
    s.run_script(&recipe).unwrap();
    let probe = recipe.as_array().unwrap().last().unwrap()["params"].clone();
    let p = s.execute("solve.probe", &probe).unwrap();
    let d = &p["nodes"][0]["displacement"];
    assert!((d[1].as_f64().unwrap() / 0.003 - 1.0).abs() < 1e-7, "{p}");
    assert!((d[0].as_f64().unwrap() / (1.0 / 300000.0) - 1.0).abs() < 1e-7, "{p}");
    println!("axishell unsymmetric annulus: {p}");
    let root = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
    let path = root.join(format!("axishell-ply-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut open = Session::default();
    open.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(p, open.execute("solve.probe", &probe).unwrap());
    if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_none() {
        std::fs::remove_file(path).unwrap();
    }
}

#[test]
fn annular_flat_preload_buckling_matches_bessel_root_and_reopens() {
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/axishell-annular-buckling/recipe.json")).unwrap();
    let mut s = Session::default();
    for c in recipe.as_array().unwrap() {
        if c["command"] != "solve.buckle" {
            s.execute(c["command"].as_str().unwrap(), &c["params"]).unwrap();
        }
    }
    let out = s.execute("solve.buckle", &json!({"step":"Buckle"})).unwrap();
    println!("axisymmetric annular buckling: {out}");
    let factor = out["modes"][0]["load_factor"].as_f64().unwrap();
    assert!((factor / 93.80004147396585 - 1.0).abs() < 0.01, "{out}");
    let p = json!({"step":"Buckle","nodes":[17]});
    let before = s.execute("solve.probe", &p).unwrap();
    let root = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
    let path = root.join(format!("axishell-buckle-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut open = Session::default();
    open.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(before, open.execute("solve.probe", &p).unwrap());
    if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_none() {
        std::fs::remove_file(path).unwrap();
    }
}
