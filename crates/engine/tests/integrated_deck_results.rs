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

fn history_outputs(s: &mut Session, displacement: f64, reaction: f64) -> Vec<Value> {
    s.execute("solve.static", &json!({"step":"Static","threads":1})).unwrap();
    let fieldcraft_engine::StepResult::Implicit(fieldcraft_solver::StepResult::Static(r)) = s.result("Static").unwrap().as_ref() else {
        panic!("static")
    };
    assert!((r.displacement[1][0] - displacement).abs() < 1e-12);
    assert!((r.reaction[0][0] - reaction).abs() < 1e-12);
    assert!((r.stress[1][0] - 4. * displacement).abs() < 1e-12);
    let modes = s.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent","threads":1})).unwrap();
    assert!((modes["modes"][0]["frequency"].as_f64().unwrap() - 1. / std::f64::consts::PI).abs() < 1e-12);
    s.execute("solve.harmonic", &json!({"step":"Response","mass":"consistent","threads":1})).unwrap();
    let fieldcraft_engine::StepResult::Implicit(fieldcraft_solver::StepResult::Harmonic(h)) = s.result("Response").unwrap().as_ref() else {
        panic!("harmonic")
    };
    let omega = 2. * std::f64::consts::PI * 0.1;
    assert!((h.frames[0].real[1][0] - 1. / (4. - omega * omega)).abs() < 1e-12);
    assert_eq!(h.mass, fieldcraft_solver::MassKind::Consistent);
    assert!(h.frames[0].relative_residual < 1e-12);
    ["Static", "Modes", "Response"].map(|step| s.execute("solve.probe", &json!({"step":step,"nodes":[1,2]})).unwrap()).to_vec()
}

#[test]
fn imported_static_history_harmonic_preserves_incremental_loads_and_linear_base() {
    let static_step = "*STEP,NAME=Static\n*STATIC\n*CLOAD\n2,1,1\n*END STEP\n";
    // Preserve the original three-step source and order, formerly explicitly refused.
    let original = DECK.replace(static_step, "").replace("*STEP,NAME=Modes", &format!("{static_step}*STEP,NAME=Modes"));
    let gravity = original.replacen("*CLOAD\n2,1,1", "*DLOAD\nBar,GRAV,1,1,0,0\n*CLOAD\n2,1,1", 1);
    let amplitude = original.replacen("*STEP,NAME=Static", "*AMPLITUDE,NAME=Ramp,TIME=TOTAL TIME\n0,0,1,2,2,3,4,9\n*STEP,NAME=Static", 1).replacen(
        "*CLOAD\n2,1,1",
        "*CLOAD,AMPLITUDE=Ramp\n2,1,1",
        1,
    ) + "*STEP,NAME=Resume\n*STATIC\n*END STEP\n";
    let root = std::env::temp_dir().join(format!("fieldcraft-integrated-static-history-{}", std::process::id()));
    for (case, source, u, rf) in [("original", original.clone(), 0.25, -1.), ("gravity", gravity, 0.625, -4.), ("amplitude", amplitude, 0.5, -2.)] {
        let dir = root.join(case);
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("input.inp");
        std::fs::write(&input, &source).unwrap();
        let mut s = Session::default();
        let imported = s.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(imported["analysis_complete"], true, "{imported}");
        for name in ["Modes", "Response"] {
            let step = s.model.step(name).unwrap();
            assert_eq!(step.base_state.as_deref(), Some("Static"));
            assert!(step.linear_base);
        }
        assert!(s.model.step("Modes").unwrap().loads.is_empty());
        assert_eq!(s.model.step("Response").unwrap().loads.len(), 1);
        let expected = history_outputs(&mut s, u, rf);
        // Perturbation steps do not overwrite general-step load/amplitude state or time.
        if case == "amplitude" {
            let r = s.execute("solve.static", &json!({"step":"Resume","threads":1})).unwrap();
            assert!((r["max_displacement"]["value"].as_f64().unwrap() - 0.75).abs() < 1e-12);
        }
        let definition = s.model.clone();
        let archive = dir.join("results.fieldcraft");
        s.execute("file.save_as", &json!({"path":archive,"results":true})).unwrap();
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":archive})).unwrap();
        assert_eq!(reopened.model, definition);
        for (step, old) in ["Static", "Modes", "Response"].into_iter().zip(&expected) {
            assert_eq!(&reopened.execute("solve.probe", &json!({"step":step,"nodes":[1,2]})).unwrap(), old);
            assert_eq!(reopened.result(step).unwrap().as_ref(), s.result(step).unwrap().as_ref());
            use fieldcraft_io::results::{Decoded, Solved, Stored, StoredRef, decode, encode};
            let fieldcraft_engine::StepResult::Implicit(r) = s.result(step).unwrap().as_ref() else { panic!("implicit") };
            let stored = Solved::Implicit(r.clone());
            let bytes = encode(&[(step, StoredRef { definition: [7; 32], result: stored.as_ref() })]).unwrap();
            assert_eq!(&bytes[..8], b"MCRSLT07");
            assert_eq!(decode(&bytes).unwrap(), Decoded::Results(vec![Stored { definition: [7; 32], result: stored }]));
        }
        let replay = dir.join("replay.inp");
        reopened.execute("deck.source.write", &json!({"path":replay})).unwrap();
        assert_eq!(std::fs::read_to_string(replay).unwrap(), source);
        let normalized = dir.join("normalized.inp");
        reopened.execute("deck.write", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
        let text = std::fs::read_to_string(&normalized).unwrap();
        assert!(text.contains("NAME=Static, NLGEOM=NO"));
        assert!(!text.contains("NLGEOM=YES"));
        let r = reopened.execute("deck.read", &json!({"format":"abaqus","path":normalized,"units":"m-kg-s-N-Pa"})).unwrap();
        assert_eq!(r["analysis_complete"], true, "{r}");
        assert_eq!(history_outputs(&mut reopened, u, rf), expected);
        if case == "amplitude" {
            let r = reopened.execute("solve.static", &json!({"step":"Resume","threads":1})).unwrap();
            assert!((r["max_displacement"]["value"].as_f64().unwrap() - 0.75).abs() < 1e-12);
        }
        if let Ok(folder) = std::env::var("FIELDCRAFT_TEST_ARTIFACTS") {
            let folder = std::path::Path::new(&folder).join(case);
            std::fs::create_dir_all(&folder).unwrap();
            for name in ["input.inp", "normalized.inp", "results.fieldcraft"] {
                std::fs::copy(dir.join(name), folder.join(name)).unwrap();
            }
            std::fs::write(folder.join("probes.json"), serde_json::to_vec_pretty(&expected).unwrap()).unwrap();
        }
    }
    // Explicit current-step amplitude is still unsupported; do not silently discard it.
    let invalid = original
        .replacen("*STEP,NAME=Static", "*AMPLITUDE,NAME=Bad\n0,1,1,1\n*STEP,NAME=Static", 1)
        .replace("*STEADY STATE DYNAMICS,DIRECT\n0.1,0\n*CLOAD", "*STEADY STATE DYNAMICS,DIRECT\n0.1,0\n*CLOAD,AMPLITUDE=Bad");
    let input = root.join("explicit-amplitude.inp");
    std::fs::write(&input, invalid).unwrap();
    let mut refused = Session::default();
    let reply = refused.execute("deck.read", &json!({"format":"abaqus","path":input,"units":"m-kg-s-N-Pa"})).unwrap();
    assert_eq!(reply["analysis_complete"], false);
    assert!(reply["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("amplitude on non-static step")));
    assert!(refused.model.step("Response").is_none());
    std::fs::remove_dir_all(root).unwrap();
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
