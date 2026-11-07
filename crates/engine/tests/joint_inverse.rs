use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn springs() -> Session {
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    for i in 1..=3 {
        s.run_script(&json!([
 {"command":"mesh.node.add","params":{"label":i,"xyz":[i,0,0]}},
 {"command":"model.set.create","params":{"name":format!("n{i}"),"nodes":[i]}},
 {"command":"model.mass.create","params":{"name":format!("m{i}"),"target":{"set":format!("n{i}")},"mass":1}},
 {"command":"model.spring.create","params":{"name":format!("k{i}"),"type":"spring","first":{"set":format!("n{i}")},"direction":{"dof":1},"stiffness":4*i}},
 {"command":"model.boundary_condition.create","params":{"name":format!("g{i}"),"type":"displacement","target":{"set":format!("n{i}")},"values":[null,0,0,null,null,null]}},
 {"command":"model.load.create","params":{"name":format!("f{i}"),"type":"force","target":{"set":format!("n{i}")},"vector":[1,0,0]}}
 ])).unwrap();
    }
    s.execute(
        "model.step.create",
        &json!({"name":"S","procedure":{"type":"static"},"boundary_conditions":["g1","g2","g3"],"loads":["f1","f2","f3"]}),
    )
    .unwrap();
    s
}
fn request(n: usize, noisy: bool) -> Value {
    let noise = [0.2, -0.3, 0.1];
    json!({"parameters":(0..n).map(|i|json!({"parameter":{"type":"spring_stiffness","name":format!("k{}",i+1)},"lower_si":1,"upper_si":30,"scale_si":4*(i+1)})).collect::<Vec<_>>(),
 "observations":(0..n).map(|i|json!({"response":{"kind":"static","step":"S","node":i+1,"dof":1},"value_si":1.0/(5.0*(i+1) as f64)+if noisy {noise[i]*0.0001}else{0.0},"standard_deviation_si":0.0001})).collect::<Vec<_>>()})
}
#[test]
fn joint_two_and_three_springs_recover_truth_and_joint_noise_ellipsoid() {
    for n in [2, 3] {
        for noisy in [false, true] {
            let mut s = springs();
            let before = s.model.clone();
            let r = s.execute("solve.inverse.joint", &request(n, noisy)).unwrap();
            assert_eq!(s.model, before);
            let mut distance = 0.0;
            for i in 0..n {
                let x = r["estimate_si"][i].as_f64().unwrap();
                let truth = 5.0 * (i + 1) as f64;
                if !noisy {
                    assert!((x / truth - 1.0).abs() < 1e-6);
                }
                distance += (x - truth) * (x - truth) / r["covariance_si"][i][i].as_f64().unwrap();
            }
            assert!(distance < r["joint_chi_square"].as_f64().unwrap());
            assert_eq!(r["rank"], n);
        }
    }
}
#[test]
fn regularization_cannot_hide_unobserved_parameter() {
    let mut s = springs();
    let mut p = request(2, false);
    p["observations"][1] = p["observations"][0].clone();
    p["regularization"] = json!(10);
    let e = s.execute("solve.inverse.joint", &p).unwrap_err().to_string();
    assert!(e.contains("unidentifiable"), "{e}");
}

#[test]
fn joint_result_is_bit_identical_for_one_and_eight_threads() {
    let mut s = springs();
    let mut p = request(3, true);
    p["threads"] = json!(1);
    let one = s.execute("solve.inverse.joint", &p).unwrap();
    p["threads"] = json!(8);
    let eight = s.execute("solve.inverse.joint", &p).unwrap();
    assert_eq!(one, eight);
    p["threads"] = json!(48);
    let forty_eight = s.execute("solve.inverse.joint", &p).unwrap();
    assert_eq!(one, forty_eight);
}

#[test]
fn identified_regularized_fit_does_not_claim_sampling_covariance() {
    let mut s = springs();
    let mut p = request(2, false);
    p["regularization"] = json!(1.0);
    let r = s.execute("solve.inverse.joint", &p).unwrap();
    assert_eq!(r["rank"], 2);
    assert!(r["covariance_si"].is_null());
    assert!(r["joint_projection_intervals_si"].is_null());
}

#[test]
fn explicit_si_joint_contract_is_independent_of_display_units() {
    let mut s = springs();
    let p = request(2, false);
    let si = s.execute("solve.inverse.joint", &p).unwrap();
    s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
    let display = s.execute("solve.inverse.joint", &p).unwrap();
    assert_eq!(si, display);
}

#[test]
fn joint_iteration_cap_does_not_limit_modal_basis() {
    let mut s = springs();
    s.execute("model.step.create", &json!({"name":"F","procedure":{"type":"frequency","modes":2},"boundary_conditions":["g1","g2","g3"]})).unwrap();
    let mut p = request(2, false);
    p["max_iterations"] = json!(1);
    p["threads"] = json!(48);
    for i in 0..2 {
        p["observations"][i]["value_si"] = json!(1.0 / (4.0 * (i + 1) as f64));
    }
    let frequencies = s.execute("solve.frequency", &json!({"step":"F","threads":48})).unwrap();
    p["observations"].as_array_mut().unwrap().push(
        json!({"response":{"kind":"frequency","step":"F","mode":0},"value_si":frequencies["modes"][0]["frequency"],"standard_deviation_si":0.001}),
    );
    let r = s.execute("solve.inverse.joint", &p).unwrap();
    assert_eq!(r["rank"], 2);
    assert_eq!(r["iterations"], 0);
}
