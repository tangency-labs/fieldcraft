//! Closed-form two-coordinate spring matrix, implemented independently of solver output.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn system(k1: f64, k2: f64, c: f64) -> Session {
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    for (i, k) in [(1, k1), (2, k2)] {
        s.run_script(&json!([
   {"command":"mesh.node.add","params":{"label":i,"xyz":[i,0,0]}},
   {"command":"model.set.create","params":{"name":format!("n{i}"),"nodes":[i]}},
   {"command":"model.mass.create","params":{"name":format!("m{i}"),"target":{"set":format!("n{i}")},"mass":1}},
   {"command":"model.spring.create","params":{"name":format!("k{i}"),"type":"spring","first":{"set":format!("n{i}")},"direction":{"dof":1},"stiffness":k}},
   {"command":"model.boundary_condition.create","params":{"name":format!("g{i}"),"type":"displacement","target":{"set":format!("n{i}")},"values":[null,0,0,null,null,null]}}
  ])).unwrap();
    }
    if c > 0.0 {
        s.execute(
            "model.spring.create",
            &json!({"name":"coupling","type":"spring","first":{"set":"n1"},"second":{"set":"n2"},"direction":{"dof":1},"stiffness":c}),
        )
        .unwrap();
    }
    s.execute("model.step.create", &json!({"name":"F","procedure":{"type":"frequency","modes":2},"boundary_conditions":["g1","g2"]})).unwrap();
    s
}
fn parameters() -> Value {
    json!([
 {"parameter":{"type":"spring_stiffness","name":"k1"},"lower_si":1,"upper_si":9,"scale_si":4},
 {"parameter":{"type":"spring_stiffness","name":"k2"},"lower_si":1,"upper_si":9,"scale_si":4}])
}
fn request(k1: f64, k2: f64, c: f64) -> Value {
    let d = ((k1 - k2).powi(2) + 4.0 * c * c).sqrt();
    let lambda = [(k1 + k2 + 2.0 * c - d) / 2.0, (k1 + k2 + 2.0 * c + d) / 2.0];
    let ratio = lambda.map(|l| (k1 + c - l) / c);
    let mut observations = Vec::new();
    for (i, lambda) in lambda.iter().enumerate() {
        observations.push(json!({"response":{"kind":"matched_frequency","group":"g","measured_mode":i},"value_si":lambda.sqrt()/std::f64::consts::TAU,"standard_deviation_si":0.001}));
    }
    observations
        .push(json!({"response":{"kind":"mode_shape","group":"g","measured_mode":0,"sensor":1},"value_si":ratio[0],"standard_deviation_si":0.001}));
    json!({"parameters":parameters(),"observations":observations,"modal_observations":[{"id":"g","step":"F","sensors":[{"node":1,"dof":1},{"node":2,"dof":1}],"measured_modes":[[1,ratio[0]],[1,ratio[1]]],"min_mac":0.5,"min_assignment_margin":0.01,"normalization":{"kind":"reference","reference_sensor":0}}]})
}
#[test]
fn normalized_nelson_components_drive_joint_command_and_match_differences() {
    let mut s = system(4.0, 5.0, 0.3);
    let original = s.model.clone();
    let p = request(4.2, 5.3, 0.3);
    let r = s.execute("solve.inverse.joint", &p).unwrap();
    assert_eq!(s.model, original);
    for (i, t) in [4.2, 5.3].iter().enumerate() {
        assert!((r["estimate_si"][i].as_f64().unwrap() / t - 1.0).abs() < 1e-7, "{r}");
    }
    assert_eq!(r["modal_pairings"][0]["pairing"]["pairs"][0]["computed_mode"], 0);
    for col in 0..2 {
        let mut a = [4.2, 5.3];
        let h = 1e-5;
        a[col] += h;
        let plus = request(a[0], a[1], 0.3);
        a[col] -= 2.0 * h;
        let minus = request(a[0], a[1], 0.3);
        let fd = (plus["observations"][2]["value_si"].as_f64().unwrap() - minus["observations"][2]["value_si"].as_f64().unwrap()) / (2.0 * h);
        let actual = r["jacobian_si"][2][col].as_f64().unwrap();
        assert!((actual - fd).abs() < 1e-7, "{actual} vs {fd}");
    }
    let mut zero_ref = p.clone();
    zero_ref["modal_observations"][0]["measured_modes"][0] = json!([0, 1]);
    assert!(s.execute("solve.inverse.joint", &zero_ref).is_err());
    let mut ambiguous = p.clone();
    ambiguous["modal_observations"][0]["measured_modes"][1] = ambiguous["modal_observations"][0]["measured_modes"][0].clone();
    assert!(s.execute("solve.inverse.joint", &ambiguous).is_err());
}
#[test]
fn crossing_frequency_observations_follow_shapes_not_order() {
    // Diagonal two-DOF model: ordered frequencies are ambiguous about which spring changed.
    let mut s = system(4.01, 3.99, 0.0);
    let p = json!({"parameters":parameters(),"observations":[
 {"response":{"kind":"matched_frequency","group":"x","measured_mode":0},"value_si":3.98_f64.sqrt()/std::f64::consts::TAU,"standard_deviation_si":0.001},
 {"response":{"kind":"matched_frequency","group":"y","measured_mode":0},"value_si":4.02_f64.sqrt()/std::f64::consts::TAU,"standard_deviation_si":0.001}],
 "modal_observations":[
 {"id":"x","step":"F","sensors":[{"node":1,"dof":1},{"node":2,"dof":1}],"measured_modes":[[1,0]],"min_mac":0.99,"min_assignment_margin":0.01,"normalization":{"kind":"reference","reference_sensor":0}},
 {"id":"y","step":"F","sensors":[{"node":1,"dof":1},{"node":2,"dof":1}],"measured_modes":[[0,-2]],"min_mac":0.99,"min_assignment_margin":0.01,"normalization":{"kind":"reference","reference_sensor":1}}]});
    let r = s.execute("solve.inverse.joint", &p).unwrap();
    assert!((r["estimate_si"][0].as_f64().unwrap() - 3.98).abs() < 1e-7);
    assert!((r["estimate_si"][1].as_f64().unwrap() - 4.02).abs() < 1e-7);
    assert_eq!(r["modal_pairings"][0]["pairing"]["pairs"][0]["computed_mode"], 0);
    // Initial x was index1, final x index0; order-only fitting would retain the swapped identity.
}

#[test]
fn harmonic_batch_reuses_only_matching_invocation_local_bases() {
    let mut s = system(4.0, 5.0, 0.3);
    s.execute("model.load.create", &json!({"name":"P","type":"force","target":{"set":"n1"},"vector":[1,0,0]})).unwrap();
    for (name, frequencies, loss, amplitude) in
        [("Zero", vec![0.0], 0.0, 1.0), ("Positive", vec![0.1, 0.2], 0.03, 1.0), ("OtherLoad", vec![0.1], 0.03, 2.0)]
    {
        s.execute("model.step.create", &json!({"name":name,"procedure":{"type":"harmonic","frequencies":frequencies,"method":{"type":"modal","modes":2,"residual_vectors":1},"damping":{"alpha":0.2,"beta":0.01,"loss_factor":loss},"load_factors":{"P":{"amplitude":amplitude,"phase":0}}},"boundary_conditions":["g1","g2"],"loads":["P"]})).unwrap();
    }
    let mut requests = Vec::new();
    for k in ["k1", "k2"] {
        for step in ["Zero", "Positive"] {
            requests.push(json!({"step":step,"parameter":{"type":"spring_stiffness","name":k},"node":1,"dof":1}));
        }
    }
    let original = s.model.clone();
    let batch = s.execute("solve.sensitivity.harmonic", &json!({"requests":requests})).unwrap();
    assert_eq!(s.model, original);
    assert_eq!(batch["modal_basis_builds"], 1);
    for (i, r) in requests.iter().enumerate() {
        let single = s.execute("solve.sensitivity.harmonic", r).unwrap();
        assert_eq!(batch["results"][i]["frames"], single["frames"]);
    }
    requests.push(json!({"step":"OtherLoad","parameter":{"type":"spring_stiffness","name":"k1"},"node":1}));
    let other = s.execute("solve.sensitivity.harmonic", &json!({"requests":requests})).unwrap();
    assert_eq!(other["modal_basis_builds"], 2);
    s.execute("model.units", &json!({"system":"mm-t-s-N-MPa"})).unwrap();
    let displayed = s.execute("solve.sensitivity.harmonic", &json!({"requests":requests})).unwrap();
    for (i, r) in requests.iter().enumerate() {
        let single = s.execute("solve.sensitivity.harmonic", r).unwrap();
        assert_eq!(displayed["results"][i]["frames"], single["frames"]);
        assert_eq!(displayed["results"][i]["units"], single["units"]);
    }

    assert!(s.execute("solve.sensitivity.harmonic", &json!({"requests":[]})).is_err());
}
