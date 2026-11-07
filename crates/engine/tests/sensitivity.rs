use fieldcraft_engine::Session;
use serde_json::json;
fn oscillator() -> Session {
    let mut s = Session::default();
    s.run_script(&json!([
        {"command":"model.units","params":{"system":"m-kg-s-N-Pa"}},
        {"command":"mesh.node.add","params":{"label":1,"xyz":[0,0,0]}},
        {"command":"model.set.create","params":{"name":"tip","nodes":[1]}},
        {"command":"model.mass.create","params":{"name":"m","target":{"set":"tip"},"mass":1}},
        {"command":"model.spring.create","params":{"name":"k","type":"spring","first":{"set":"tip"},"direction":{"dof":1},"stiffness":4}},
        {"command":"model.boundary_condition.create","params":{"name":"guide","type":"displacement","target":{"set":"tip"},"values":[null,0,0,null,null,null]}},
        {"command":"model.load.create","params":{"name":"force","type":"force","target":{"set":"tip"},"vector":[1,0,0]}},
        {"command":"model.step.create","params":{"name":"Load","procedure":{"type":"static"},"boundary_conditions":["guide"],"loads":["force"]}},
        {"command":"model.step.create","params":{"name":"Modes","procedure":{"type":"frequency","modes":1},"boundary_conditions":["guide"]}}
    ])).unwrap();
    s
}
#[test]
fn sensitivity_commands_convert_quotient_units_and_leave_model_unchanged() {
    let mut s = oscillator();
    let before = s.model.clone();
    let params = json!({"step":"Load","parameter":{"type":"spring_stiffness","name":"k"},"node":1,"threads":1});
    let r = s.execute("solve.sensitivity.static", &params).unwrap();
    assert!((r["direct"].as_f64().unwrap() + 0.0625).abs() < 1e-10);
    assert_eq!(s.model, before);
    let m = s.execute("solve.sensitivity.frequency", &json!({"step":"Modes","parameter":{"type":"spring_stiffness","name":"k"}})).unwrap();
    let f = 1.0 / std::f64::consts::PI;
    assert!((m["modes"][0]["frequency_derivative"].as_f64().unwrap() - f / 8.0).abs() < 1e-9);
    s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
    let r = s.execute("solve.sensitivity.static", &params).unwrap();
    assert!((r["displacement"].as_f64().unwrap() - 250.0).abs() < 1e-9);
    assert!((r["direct"].as_f64().unwrap() + 62500.0).abs() < 1e-4);
    assert_eq!(r["units"]["direct"], "mm/(N/mm)");
    let before = s.model.clone();
    assert!(
        s.execute(
            "solve.sensitivity.static",
            &json!({"step":"Load","parameter":{"type":"spring_stiffness","name":"k"},"node":1,"relative_increment":0})
        )
        .is_err()
    );
    assert_eq!(s.model, before);
}

#[test]
fn bounded_frequency_fit_recovers_known_spring_and_reports_interval() {
    let mut s = oscillator();
    let before = s.model.clone();
    let target = 9.0_f64.sqrt() / (2.0 * std::f64::consts::PI);
    let p = json!({"step":"Modes","parameter":{"type":"spring_stiffness","name":"k"},"lower":"1 N/m","upper":"16 N/m","observed":[target],"standard_deviation":[0.001]});
    let result = s.execute("solve.inverse.frequency", &p).unwrap();
    assert!((result["estimate"].as_f64().unwrap() / 9.0 - 1.0).abs() < 1e-6, "{result}");
    let ci = result["confidence_interval"].as_array().unwrap();
    assert!(ci[0].as_f64().unwrap() < 9.0 && ci[1].as_f64().unwrap() > 9.0);
    assert_eq!(s.model, before);
    let mut noisy = p.clone();
    noisy["observed"] = json!([target + 0.0002]);
    let result = s.execute("solve.inverse.frequency", &noisy).unwrap();
    let ci = result["confidence_interval"].as_array().unwrap();
    assert!(ci[0].as_f64().unwrap() < 9.0 && ci[1].as_f64().unwrap() > 9.0);
    let mut bound = p.clone();
    bound["upper"] = json!(8.0);
    let result = s.execute("solve.inverse.frequency", &bound).unwrap();
    assert_eq!(result["estimate"], 8.0);
    assert_eq!(result["at_bound"], true);
    assert!(result["confidence_interval"].is_null());
}

#[test]
fn inverse_bounds_and_estimate_follow_parameter_display_units() {
    let mut s = oscillator();
    s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
    let before = s.model.clone();
    let target = 3.0 / (2.0 * std::f64::consts::PI);
    let p = json!({"step":"Modes","parameter":{"type":"spring_stiffness","name":"k"},
        "lower":0.001,"upper":0.016,"observed":[target],"standard_deviation":[0.001]});
    let result = s.execute("solve.inverse.frequency", &p).unwrap();
    assert!((result["estimate"].as_f64().unwrap() / 0.009 - 1.0).abs() < 1e-6, "{result}");
    assert_eq!(result["units"]["estimate"], "N/mm");
    let ci = result["confidence_interval"].as_array().unwrap();
    assert!(ci[0].as_f64().unwrap() < 0.009 && ci[1].as_f64().unwrap() > 0.009);
    assert_eq!(s.model, before);
}

#[test]
fn optimizer_iteration_cap_does_not_limit_the_modal_basis() {
    let mut s = oscillator();
    for node in 2..=12 {
        let name = format!("p{node}");
        s.execute("mesh.node.add", &json!({"label":node,"xyz":[node,0,0]})).unwrap();
        s.execute("model.set.create", &json!({"name":name,"nodes":[node]})).unwrap();
        s.execute("model.mass.create", &json!({"name":name,"target":{"set":name},"mass":1})).unwrap();
        s.execute(
            "model.spring.create",
            &json!({"name":name,"type":"spring","first":{"set":name},"direction":{"dof":1},"stiffness":(node+1)*(node+1)}),
        )
        .unwrap();
    }
    s.execute("model.set.create", &json!({"name":"all","nodes":(1..=12).collect::<Vec<_>>()})).unwrap();
    s.execute(
        "model.boundary_condition.create",
        &json!({"name":"allguide","type":"displacement","target":{"set":"all"},"values":[null,0,0,null,null,null]}),
    )
    .unwrap();
    s.execute("model.step.create", &json!({"name":"Many","procedure":{"type":"frequency","modes":1},"boundary_conditions":["allguide"]})).unwrap();
    let p = json!({"step":"Many","parameter":{"type":"spring_stiffness","name":"k"},"lower":1,"upper":16,"observed":[1.0/std::f64::consts::PI],"max_iterations":1});
    let result = s.execute("solve.inverse.frequency", &p).unwrap();
    assert!((result["estimate"].as_f64().unwrap() - 4.0).abs() < 1e-8, "{result}");
}
