//! Original beam-end torsional oscillator with an explicit global rotary tensor.
use fieldcraft_engine::Session;
use serde_json::json;
const INPUT: &str = "*NODE\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=B31,ELSET=Beam\n1,1,2\n*ELEMENT,TYPE=ROTARYI,ELSET=Rotor\n10,2\n*MATERIAL,NAME=Elastic\n*ELASTIC\n1000000,.25\n*DENSITY\n1000\n*BEAM SECTION,ELSET=Beam,MATERIAL=Elastic,SECTION=CIRC\n.1\n0,0,-1\n*ROTARY INERTIA,ELSET=Rotor\n.1,.2,.25,.01,-.02,.03\n*BOUNDARY\n1,1,6\n2,1,3\n2,5,6\n*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY\n1\n*END STEP\n*STEP,NAME=Response,PERTURBATION\n*STEADY STATE DYNAMICS,DIRECT\n1,0\n*CLOAD\n2,4,.01\n*END STEP\n*STEP,NAME=Static\n*STATIC\n1,1\n*CLOAD,OP=NEW\n2,4,.01\n*END STEP\n";
#[test]
fn rotary_inertia_import_solves_and_roundtrips() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-rotary-import-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.inp");
    std::fs::write(&input, INPUT).unwrap();
    let mut s = Session::default();
    assert_eq!(s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap()["analysis_complete"], true);
    let j = std::f64::consts::PI * 0.1_f64.powi(4) / 2.0;
    let expected = ((400000.0 * j) / (0.1 + 1000.0 * j / 3.0)).sqrt() / (2.0 * std::f64::consts::PI);
    let modes = s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent","threads":1})).unwrap();
    assert!((modes["modes"][0]["frequency"].as_f64().unwrap() - expected).abs() < 1e-11);
    let probe = s.execute("solve.probe", &json!({"step":"Modes","nodes":[1,2]})).unwrap();
    let saved = dir.join("saved.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    assert_eq!(s.execute("solve.probe", &json!({"step":"Modes","nodes":[1,2]})).unwrap(), probe);
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), INPUT);
    let normalized = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap()["analysis_complete"], true);
    s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent","threads":1})).unwrap();
    assert_eq!(s.execute("solve.probe", &json!({"step":"Modes","nodes":[1,2]})).unwrap(), probe);
    assert_eq!(s.model.masses().values().next().unwrap().inertia, Some([0.1, 0.2, 0.25, 0.01, -0.02, 0.03]));
    std::fs::remove_dir_all(dir).unwrap();
}
