//! Original mass-spring model through commands, assembly, adaptive integration and storage.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn session(tolerance: f64, estimator: &str) -> Session {
    let mut s = Session::default();
    s.run_script(&serde_json::from_str(include_str!("data/transient-sdof.json")).unwrap()).unwrap();
    let mut p = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    p["time_step"] = json!(0.2);
    p["increments"] = json!(5);
    p["output_every"] = json!(1);
    p["adaptive"] = json!({"min_step":1e-5,"max_step":0.2,"absolute_tolerance":tolerance,"relative_tolerance":0,
        "max_attempts":10000,"estimator":estimator});
    s.execute("model.step.edit", &json!({"name":"STEP_FORCE","procedure":p})).unwrap();
    s
}
fn solve(s: &mut Session) -> (Value, f64) {
    let result = s.execute("solve.transient", &json!({"threads":1})).unwrap();
    let p = s.execute("solve.probe", &json!({"step":"STEP_FORCE","nodes":[1]})).unwrap();
    let error = p["nodes"][0]["samples"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            let t = v["time"].as_f64().unwrap();
            (v["displacement"][0].as_f64().unwrap() - (1.0 - (2.0 * t).cos()) / 4.0).abs()
        })
        .fold(0.0_f64, f64::max);
    println!("ADAPTIVE_COMMAND {result}, max_u_error_m={error:.16e}");
    (p, error)
}
#[test]
fn error_control_converges_and_persists_on_the_requested_grid() {
    for estimator in ["zienkiewicz_xie", "quadratic_acceleration"] {
        let (_, coarse) = solve(&mut session(1e-6, estimator));
        let mut s = session(1e-8, estimator);
        let (p, fine) = solve(&mut s);
        assert!(fine < 1e-5 && fine < 0.3 * coarse, "{estimator}: {fine} vs {coarse}");
        assert_eq!(p["nodes"][0]["samples"].as_array().unwrap().len(), 6);
        assert_eq!(p["nodes"][0]["samples"][0]["acceleration"][0], 1.0);
        for frame in p["frames"].as_array().unwrap() {
            assert!(frame["energy_balance"].as_f64().unwrap().abs() < 1e-12);
        }
        let path = std::env::temp_dir().join(format!("adaptive-{estimator}-{}.fieldcraft", std::process::id()));
        s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
        let mut reopened = Session::default();
        reopened.execute("file.open", &json!({"path":path})).unwrap();
        assert_eq!(reopened.execute("solve.probe", &json!({"step":"STEP_FORCE","nodes":[1]})).unwrap(), p);
        let _ = std::fs::remove_file(path);
    }
}
#[test]
fn load_corner_and_dimensioned_tolerances_reach_the_command_path() {
    let mut s = session(1e-8, "quadratic_acceleration");
    s.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let p = json!({"type":"transient","time_step":"0.2 s","increments":5,
        "adaptive":{"min_step":"0.00001 s","max_step":"0.2 s","absolute_tolerance":"0.00001 mm",
            "max_attempts":10000,"estimator":"quadratic_acceleration"},
        "load_amplitudes":{"F":{"type":"tabular","points":[["0 s",0],["0.37 s",1],["1 s",0]]}}});
    s.execute("model.step.edit", &json!({"name":"STEP_FORCE","procedure":p})).unwrap();
    s.execute("solve.transient", &json!({"threads":1})).unwrap();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let p = s.execute("solve.probe", &json!({"step":"STEP_FORCE","nodes":[1]})).unwrap();
    let g = |t: f64| if t > 0.0 { t / 4.0 - (2.0 * t).sin() / 8.0 } else { 0.0 };
    for v in p["nodes"][0]["samples"].as_array().unwrap() {
        let t = v["time"].as_f64().unwrap();
        let u = g(t) / 0.37 - (1.0 / 0.63 + 1.0 / 0.37) * g(t - 0.37) + g(t - 1.0) / 0.63;
        assert!((v["displacement"][0].as_f64().unwrap() - u).abs() < 1e-5);
    }
}
#[test]
fn exhausted_budgets_and_unrepresentable_options_are_explicit_errors() {
    for (key, value, expected) in [("min_step", json!(0.2), "min_step"), ("max_attempts", json!(1), "attempts")] {
        let mut s = session(1e-12, "zienkiewicz_xie");
        let mut p = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
        p["adaptive"][key] = value;
        s.execute("model.step.edit", &json!({"name":"STEP_FORCE","procedure":p})).unwrap();
        let e = s.execute("solve.transient", &json!({"threads":1})).unwrap_err().to_string();
        assert!(e.contains(expected), "{e}");
        assert!(s.result("STEP_FORCE").is_none());
    }
    let mut s = session(1e-8, "zienkiewicz_xie");
    let mut p = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    p["integration"] = json!({"type":"hht","alpha":-0.1});
    assert!(s.execute("model.step.edit", &json!({"name":"STEP_FORCE","procedure":p})).unwrap_err().to_string().contains("adaptive"));
    let path = std::env::temp_dir().join("adaptive-refused.inp");
    assert!(s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap_err().to_string().contains("adaptive"));
}
