//! Original two-node axial oscillator: K=4, consistent free mass=1.
//! A linear equilibrium remains a dependency without introducing geometric stiffness.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
const INPUT: &str = r#"*NODE
1,0,0,0
2,1,0,0
*ELEMENT,TYPE=T3D2,ELSET=Bar
1,1,2
*MATERIAL,NAME=Elastic
*ELASTIC
4,0.25
*DENSITY
3
*SOLID SECTION,ELSET=Bar,MATERIAL=Elastic
1
*BOUNDARY
1,1,3
2,2,3
*STEP,NAME=Static
*STATIC
*CLOAD
2,1,1
*END STEP
*STEP,NAME=Modes,PERTURBATION
*FREQUENCY
1
*END STEP
*STEP,NAME=Response
*STEADY STATE DYNAMICS,DIRECT
0.1,0
*CLOAD
2,1,1
*END STEP
"#;
fn solve(s: &mut Session) -> Vec<Value> {
    let equilibrium = s.execute("solve.static", &json!({"step":"Static"})).unwrap();
    assert!((equilibrium["max_displacement"]["value"].as_f64().unwrap() - 0.25).abs() < 1e-12);
    let modes = s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"})).unwrap();
    assert!((modes["modes"][0]["frequency"].as_f64().unwrap() * std::f64::consts::PI - 1.).abs() < 1e-12);
    s.execute("solve.harmonic", &json!({"step":"Response","mass":"consistent"})).unwrap();
    let fieldcraft_engine::StepResult::Implicit(fieldcraft_solver::StepResult::Harmonic(h)) = s.result("Response").unwrap().as_ref() else {
        panic!("harmonic result")
    };
    let tip = h.nodes.iter().position(|n| *n == 2).unwrap();
    let omega = 2. * std::f64::consts::PI * 0.1;
    assert!((h.frames[0].real[tip][0] - 1. / (4. - omega * omega)).abs() < 1e-12);
    assert!(h.frames[0].relative_residual < 1e-12);
    ["Static", "Modes", "Response"].map(|step| s.execute("solve.probe", &json!({"step":step,"nodes":[1,2]})).unwrap()).to_vec()
}
#[test]
fn static_history_harmonic_keeps_linear_base_and_independent_loads_through_roundtrip() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-linear-history-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // Extra load at a constrained node acts only in the static equilibrium.
    // A stale inherited load would change the harmonic reaction.
    let source = INPUT.replacen("*CLOAD\n2,1,1", "*CLOAD\n1,1,7\n2,1,1", 1);
    let input = dir.join("input.inp");
    std::fs::write(&input, &source).unwrap();
    let mut s = Session::default();
    let r = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(r["analysis_complete"], true, "{r}");
    assert_eq!(s.model.steps().iter().map(|s| s.name.as_str()).collect::<Vec<_>>(), vec!["Static", "Modes", "Response"]);
    for name in ["Modes", "Response"] {
        let step = s.model.step(name).unwrap();
        assert!(step.linear_base);
        assert_eq!(step.base_state.as_deref(), Some("Static"));
    }
    assert_eq!(s.model.step("Static").unwrap().loads.len(), 2);
    assert_eq!(s.model.step("Response").unwrap().loads.len(), 1);
    let expected = solve(&mut s);
    let saved = dir.join("model.fieldcraft");
    s.execute("file.save_as", &json!({"path":saved,"results":true})).unwrap();
    s.execute("file.new", &json!({})).unwrap();
    s.execute("file.open", &json!({"path":saved})).unwrap();
    for (step, probe) in ["Static", "Modes", "Response"].into_iter().zip(&expected) {
        assert_eq!(&s.execute("solve.probe", &json!({"step":step,"nodes":[1,2]})).unwrap(), probe);
    }
    let replay = dir.join("replay.inp");
    s.execute("deck.source.write", &json!({"path":replay})).unwrap();
    assert_eq!(std::fs::read_to_string(replay).unwrap(), source);
    let normalized = dir.join("normalized.inp");
    s.execute("deck.write", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    let text = std::fs::read_to_string(&normalized).unwrap();
    assert!(text.contains("NAME=Static, NLGEOM=NO"));
    assert!(!text.contains("NLGEOM=YES"));
    let r = s.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(r["analysis_complete"], true, "{r}");
    assert_eq!(solve(&mut s), expected);
    let before = s.model.clone();
    let error = s.execute("model.step.edit", &json!({"name":"Response","base_state":null})).unwrap_err().to_string();
    assert!(error.contains("linear_base requires base_state"), "{error}");
    assert_eq!(s.model, before);
    s.execute("model.step.edit", &json!({"name":"Response","linear_base":false})).unwrap();
    let error = s.execute("deck.write", &json!({"format":"abaqus","path":dir.join("mixed.inp")})).unwrap_err().to_string();
    assert!(error.contains("cannot share"), "{error}");
    std::fs::remove_dir_all(dir).unwrap();
}
