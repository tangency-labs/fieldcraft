//! Original axial oscillator crosses deck import, three solvers, export and native persistence.
//! K=EA/L=4 and consistent free mass=rho*A*L/3=1: u_static=1/4,
//! omega=2 and H(omega)=1/(4-omega^2). No external solver data.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

const DECK: &str = "*NODE\n1,0,0,0\n2,1,0,0\n*ELEMENT,TYPE=T3D2,ELSET=Bar\n1,1,2\n*MATERIAL,NAME=Elastic\n*ELASTIC\n4,0.25\n*DENSITY\n3\n*SOLID SECTION,ELSET=Bar,MATERIAL=Elastic\n1\n*BOUNDARY\n1,1,3\n2,2,3\n*STEP,NAME=Modes,PERTURBATION\n*FREQUENCY\n1\n*END STEP\n*STEP,NAME=Response\n*STEADY STATE DYNAMICS,DIRECT\n0.1,0\n*CLOAD\n2,1,1\n*END STEP\n*STEP,NAME=Static\n*STATIC\n*CLOAD\n2,1,1\n*END STEP\n";

fn solve(s: &mut Session) -> Vec<Value> {
    let static_result = s.execute("solve.static", &json!({"step":"Static","threads":1})).unwrap();
    assert!((static_result["max_displacement"]["value"].as_f64().unwrap() - 0.25).abs() < 1e-12);
    let modes = s.execute("solve.frequency", &json!({"step":"Modes","threads":1,"mass":"consistent"})).unwrap();
    let frequency = modes["modes"][0]["frequency"].as_f64().unwrap();
    assert!((frequency * std::f64::consts::PI - 1.).abs() < 1e-12);
    s.execute("solve.harmonic", &json!({"step":"Response","threads":1,"mass":"consistent"})).unwrap();
    let fieldcraft_engine::StepResult::Implicit(fieldcraft_solver::StepResult::Harmonic(h)) = s.result("Response").unwrap().as_ref() else {
        panic!("harmonic result");
    };
    let tip = h.nodes.iter().position(|n| *n == 2).unwrap();
    let omega = 2. * std::f64::consts::PI * 0.1;
    assert!((h.frames[0].real[tip][0] - 1. / (4. - omega * omega)).abs() < 1e-12);
    assert!(h.frames[0].relative_residual < 1e-12);
    ["Static", "Modes", "Response"].map(|step| s.execute("solve.probe", &json!({"step":step,"nodes":[2]})).unwrap()).to_vec()
}

#[test]
fn imported_three_step_model_survives_export_and_native_result_roundtrip() {
    let dir = std::env::temp_dir().join(format!("fieldcraft-integrated-deck-results-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("original.inp");
    std::fs::write(&input, DECK).unwrap();
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let imported = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(imported["analysis_complete"], true, "{imported}");
    assert_eq!(s.model.steps().len(), 3);
    let expected = solve(&mut s);
    // Each solved kind uses the merged version-07 envelope, with exact decoder roundtrip.
    for step in ["Static", "Modes", "Response"] {
        let fieldcraft_engine::StepResult::Implicit(result) = s.result(step).unwrap().as_ref() else { panic!("implicit result") };
        use fieldcraft_io::results::{Decoded, Solved, Stored, StoredRef, decode, encode};
        let stored = Solved::Implicit(result.clone());
        let bytes = encode(&[(step, StoredRef { definition: [7; 32], result: stored.as_ref() })]).unwrap();
        assert_eq!(&bytes[..8], b"MCRSLT07");
        assert_eq!(decode(&bytes).unwrap(), Decoded::Results(vec![Stored { definition: [7; 32], result: stored }]));
    }
    let archive = dir.join("results.fieldcraft");
    s.execute("file.save_as", &json!({"path":archive,"results":true})).unwrap();
    let mut reopened = Session::default();
    reopened.execute("file.open", &json!({"path":archive})).unwrap();
    for (step, reference) in ["Static", "Modes", "Response"].into_iter().zip(&expected) {
        assert_eq!(&reopened.execute("solve.probe", &json!({"step":step,"nodes":[2]})).unwrap(), reference);
    }
    let exported = dir.join("roundtrip.inp");
    reopened.execute("deck.write", &json!({"format":"abaqus","path":exported,"units":"m-kg-s-N-Pa"})).unwrap();
    let mut reread = Session::default();
    reread.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let imported = reread.execute("deck.read", &json!({"format":"abaqus","path":exported,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(imported["analysis_complete"], true, "{imported}");
    assert_eq!(solve(&mut reread), expected);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn imported_harmonic_after_static_reports_incomplete_history_without_running_it() {
    let static_step = "*STEP,NAME=Static\n*STATIC\n*CLOAD\n2,1,1\n*END STEP\n";
    let history = DECK.replace(static_step, "").replace("*STEP,NAME=Modes", &format!("{static_step}*STEP,NAME=Modes"));
    let path = std::env::temp_dir().join(format!("fieldcraft-imported-base-state-{}.inp", std::process::id()));
    std::fs::write(&path, history).unwrap();
    let mut s = Session::default();
    let imported = s.execute("deck.read", &json!({"format":"abaqus","path":path,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(imported["analysis_complete"], false);
    assert!(
        imported["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|w| w.as_str().unwrap().contains("base-state perturbation requires the prestress interface"))
    );
    assert!(s.model.step("Response").is_none());
    assert!(s.execute("solve.harmonic", &json!({"step":"Response"})).is_err());
    std::fs::remove_file(path).unwrap();
}

#[test]
fn transient_base_state_is_explicitly_refused_without_changing_the_model() {
    let mut s = Session::default();
    s.run_script(&serde_json::from_str::<Value>(include_str!("data/transient-rod.json")).unwrap()).unwrap();
    s.execute(
        "model.step.create",
        &json!({"name":"Preload","procedure":{"type":"static"},"boundary_conditions":["clamp","guide"],"loads":["force"]}),
    )
    .unwrap();
    let before = s.model.clone();
    let error = s.execute("model.step.create", &json!({"name":"BasedMotion","base_state":"Preload","procedure":{"type":"transient","time_step":0.01,"increments":10},"boundary_conditions":["clamp","guide"],"loads":["force"]})).unwrap_err().to_string();
    assert!(error.contains("base_state requires a frequency or harmonic step"), "{error}");
    assert_eq!(s.model, before);
}
