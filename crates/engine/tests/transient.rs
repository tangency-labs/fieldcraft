//! End-to-end original rod model: consistent free mass rho*A*L/3=1, K=EA/L=4.
//! This checks command -> FE assembly -> integration -> stress -> file -> probe.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn rod() -> Session {
    let mut s = Session::default();
    let recipe: Value = serde_json::from_str(include_str!("data/transient-rod.json")).unwrap();
    s.run_script(&recipe).unwrap();
    s
}
fn probe(s: &mut Session) -> Value {
    s.execute("solve.probe", &json!({"step":"Pulse","nodes":[2]})).unwrap()
}
#[test]
fn finite_element_transient_round_trip_and_display_units() {
    let mut s = rod();
    s.execute("solve.transient", &json!({"step":"Pulse","threads":1})).unwrap();
    let before = probe(&mut s);
    for sample in before["nodes"][0]["samples"].as_array().unwrap() {
        let t = sample["time"].as_f64().unwrap();
        let u = (1.0 - (2.0 * t).cos()) / 4.0;
        assert!((sample["displacement"][0].as_f64().unwrap() - u).abs() < 2e-7);
        assert!((sample["stress"][0].as_f64().unwrap() - 4.0 * u).abs() < 8e-7);
    }
    assert!((before["nodes"][0]["samples"][0]["acceleration"][0].as_f64().unwrap() - 1.0).abs() < 1e-12);
    let path = std::env::temp_dir().join(format!("fieldcraft-transient-roundtrip-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut restored = Session::default();
    restored.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(probe(&mut restored), before);
    restored.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let after = probe(&mut restored);
    let a = &after["nodes"][0]["samples"][10];
    let b = &before["nodes"][0]["samples"][10];
    assert_eq!(a["time"], 1000.0);
    assert!((a["displacement"][0].as_f64().unwrap() - 1000.0 * b["displacement"][0].as_f64().unwrap()).abs() < 1e-10);
    assert!((a["acceleration"][0].as_f64().unwrap() - 0.001 * b["acceleration"][0].as_f64().unwrap()).abs() < 1e-12);
    let _ = std::fs::remove_file(path);
}
#[test]
fn modal_transient_matches_direct_on_same_finite_element_model() {
    let mut s = rod();
    s.execute("solve.transient", &json!({"threads":1})).unwrap();
    let direct = probe(&mut s);
    s.execute("model.step.edit",&json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.001,"increments":1000,"output_every":100,"method":{"type":"modal","modes":1,"damping_ratios":[0]}}})).unwrap();
    s.execute("solve.run", &json!({})).unwrap();
    let modal = probe(&mut s);
    for (d, m) in direct["nodes"][0]["samples"].as_array().unwrap().iter().zip(modal["nodes"][0]["samples"].as_array().unwrap()) {
        assert!((d["displacement"][0].as_f64().unwrap() - m["displacement"][0].as_f64().unwrap()).abs() < 1e-11);
        assert!((d["stress"][0].as_f64().unwrap() - m["stress"][0].as_f64().unwrap()).abs() < 1e-11);
    }
}
#[test]
fn finite_element_spectrum_participation_and_persistence() {
    let mut s = rod();
    // Consistent full mass gives free row M*r=1.5, not free diagonal Mff=1.
    s.execute("model.step.edit",&json!({"name":"Pulse","loads":[],"procedure":{"type":"response_spectrum","modes":1,"damping_ratio":0.05,"combination":"cqc","direction":[1,0,0],"pseudo_acceleration":[[0.1,2.0],[1.0,2.0]]}})).unwrap();
    s.execute("solve.spectrum", &json!({"threads":1})).unwrap();
    let p = probe(&mut s);
    let sample = &p["nodes"][0]["samples"][0];
    assert!(sample["time"].is_null());
    assert!((sample["displacement"][0].as_f64().unwrap() - 0.75).abs() < 1e-10);
    assert!((sample["stress"][0].as_f64().unwrap() - 3.0).abs() < 1e-10);
    let r = s.result("Pulse").unwrap();
    let fieldcraft_engine::StepResult::Implicit(r) = r.as_ref() else { panic!() };
    let stored = fieldcraft_io::results::Solved::Implicit(r.clone());
    use fieldcraft_io::results::{Decoded, Stored, StoredRef, decode, encode};
    let data = encode(&[("Pulse", StoredRef { definition: [7; 32], result: stored.as_ref() })]).unwrap();
    assert_eq!(decode(&data).unwrap(), Decoded::Results(vec![Stored { definition: [7; 32], result: stored }]));
    let mut bad = r.clone();
    let fieldcraft_solver::StepResult::Dynamics(ref mut d) = bad else { panic!() };
    d.frames[0].time = Some(0.0);
    assert!(fieldcraft_io::results::check(&fieldcraft_io::results::Solved::Implicit(bad)).is_err());
}
#[test]
fn initial_constraints_and_inactive_amplitudes_are_rejected() {
    let mut s = rod();
    assert!(s.execute("model.step.edit",&json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.01,"increments":10,"load_amplitudes":{"missing":{"type":"tabular","points":[[0,1],[1,1]]}}}})).is_err());
    s.execute(
        "model.step.edit",
        &json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.01,"increments":10,"initial_displacement":{"1":[0.1,0,0]}}}),
    )
    .unwrap();
    assert!(s.execute("solve.transient", &json!({})).is_err());
    assert!(s.result("Pulse").is_none());
}

#[test]
fn damping_interpolation_period_units_and_direction_envelopes() {
    let mut s = rod();
    s.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let mut procedure = json!({"type":"response_spectrum","modes":1,"damping_ratio":0.05,
        "combination":"cqc","directions":[{"direction":[0.6,0.8,0],"scale":1},{"direction":[-0.8,0.6,0],"scale":1}],
        "direction_combination":"srss","period_tables":[
            {"damping_ratio":0.02,"points":[["1000 ms","2 m/s^2"],["5000 ms","2 m/s^2"]]},
            {"damping_ratio":0.08,"points":[["1000 ms","4 m/s^2"],["5000 ms","4 m/s^2"]]}]});
    s.execute("model.step.edit", &json!({"name":"Pulse","loads":[],"procedure":procedure})).unwrap();
    s.execute("solve.spectrum", &json!({})).unwrap();
    let p = probe(&mut s);
    assert!((p["nodes"][0]["samples"][0]["displacement"][0].as_f64().unwrap() - 1125.0).abs() < 1e-7);
    procedure["direction_combination"] = json!("percent1004040");
    s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).unwrap();
    s.execute("solve.spectrum", &json!({})).unwrap();
    let p = probe(&mut s);
    assert!((p["nodes"][0]["samples"][0]["displacement"][0].as_f64().unwrap() - 1170.0).abs() < 1e-7);
    // Extrapolation and non-orthogonal direction sets are rejected before solving.
    procedure["damping_ratio"] = json!(0.1);
    assert!(s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).is_err());
    procedure["damping_ratio"] = json!(0.05);
    procedure["directions"][1]["direction"] = json!([0.6, 0.8, 0]);
    assert!(s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).is_err());
}

#[test]
fn displacement_base_samples_match_acceleration_base_with_units() {
    let mut s = rod();
    s.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let mut p = json!({"type":"transient","time_step":"0.001 s","increments":1000,"output_every":100,
        "base_acceleration":{"direction":[1,0,0],"magnitude":"1 m/s^2","amplitude":{"type":"tabular","points":[["0 s",1],["1 s",1]]}}});
    s.execute("model.step.edit", &json!({"name":"Pulse","loads":[],"procedure":p})).unwrap();
    s.execute("solve.transient", &json!({"threads":1})).unwrap();
    let reference = probe(&mut s);
    p.as_object_mut().unwrap().remove("base_acceleration");
    p["base_displacement"] = json!({"direction":[1,0,0],"magnitude":"1 m","initial_velocity":"0 m/s","initial_acceleration":"1 m/s^2",
        "amplitude":{"type":"tabular","points":(0..=1000).map(|i|{let t=i as f64*0.001;json!([format!("{t} s"),0.5*t*t])}).collect::<Vec<_>>()}});
    for method in [json!({"type":"direct"}), json!({"type":"modal","modes":1})] {
        p["method"] = method;
        s.execute("model.step.edit", &json!({"name":"Pulse","procedure":p})).unwrap();
        s.execute("solve.transient", &json!({"threads":1})).unwrap();
        let result = probe(&mut s);
        for (a, b) in result["nodes"][0]["samples"].as_array().unwrap().iter().zip(reference["nodes"][0]["samples"].as_array().unwrap()) {
            let t = a["time"].as_f64().unwrap() * 0.001;
            assert!((a["displacement"][0].as_f64().unwrap() + 375.0 * (1.0 - (2.0 * t).cos())).abs() < 0.0003);
            for key in ["displacement", "velocity", "acceleration"] {
                assert!((a[key][0].as_f64().unwrap() - b[key][0].as_f64().unwrap()).abs() < 1e-6, "{key}: {a} != {b}");
            }
        }
    }
    p["base_acceleration"] = json!({"direction":[1,0,0],"magnitude":1,"amplitude":{"type":"tabular","points":[[0,1],[1,1]]}});
    assert!(s.execute("model.step.edit", &json!({"name":"Pulse","procedure":p})).is_err());
}

#[test]
fn finite_element_duhamel_pulse_and_fourier_load() {
    // Duhamel integral of a triangular pulse, M=1, K=4, peak F=1 N.
    let g = |t: f64| if t > 0.0 { t / 4.0 - (2.0 * t).sin() / 8.0 } else { 0.0 };
    for modal in [false, true] {
        for periodic in [false, true] {
            let mut s = rod();
            let amplitude = if periodic {
                json!({"type":"periodic","angular_frequency":"1 rad/s","origin":0,"offset":0,"harmonics":[[0,1]]})
            } else {
                json!({"type":"tabular","points":[[0,0],[0.2,1],[0.4,0],[1,0]]})
            };
            let method = if modal { json!({"type":"modal","modes":1}) } else { json!({"type":"direct"}) };
            s.execute("model.step.edit",&json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.001,"increments":1000,"output_every":10,"method":method,"load_amplitudes":{"force":amplitude}}})).unwrap();
            s.execute("solve.run", &json!({})).unwrap();
            let result = probe(&mut s);
            let mut error: f64 = 0.0;
            for a in result["nodes"][0]["samples"].as_array().unwrap() {
                let t = a["time"].as_f64().unwrap();
                let exact = if periodic { (t.sin() - 0.5 * (2.0 * t).sin()) / 3.0 } else { (g(t) - 2.0 * g(t - 0.2) + g(t - 0.4)) / 0.2 };
                error = error.max((a["displacement"][0].as_f64().unwrap() - exact).abs());
            }
            eprintln!("rod modal={modal} periodic={periodic} max_displacement_error={error:e}");
            assert!(error < 3e-7);
        }
    }
}

#[test]
fn cantilever_tip_pulse_matches_four_mode_euler_bernoulli_series() {
    // Original slender beam: L=1 m, A=1e-4 m2, I=1e-4/12*1e-4 m4,
    // E=12e4 Pa, rho=1e4 kg/m3 => rho*A=1 kg/m, EI=1e-4 Nm2.
    // Uniform cantilever eigenfunctions have phi_tip^2 / modal_mass = 4/(rho*A*L).
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    let n = 40;
    for i in 0..=n {
        s.execute("mesh.node.add", &json!({"label":i+1,"xyz":[i as f64/n as f64,0,0]})).unwrap();
        if i > 0 {
            s.execute("mesh.element.add", &json!({"label":i,"kind":"beam2","nodes":[i,i+1]})).unwrap();
        }
    }
    for (name, nodes, elements) in
        [("all", (1..=n + 1).collect::<Vec<_>>(), (1..=n).collect::<Vec<_>>()), ("fixed", vec![1], vec![]), ("tip", vec![n + 1], vec![])]
    {
        s.execute("model.set.create", &json!({"name":name,"nodes":nodes,"elements":elements})).unwrap();
    }
    s.execute("model.material.create", &json!({"name":"Elastic","elastic":{"type":"isotropic","young":120000,"poisson":0.3},"density":10000}))
        .unwrap();
    s.execute("model.section.create",&json!({"name":"Beam","set":"all","material":"Elastic","kind":{"type":"beam","profile":{"type":"rectangle","width":0.01,"height":0.01},"orientation":[0,1,0]}})).unwrap();
    s.execute("model.boundary_condition.create", &json!({"name":"clamp","type":"displacement","target":{"set":"fixed"},"values":[0,0,0,0,0,0]}))
        .unwrap();
    s.execute("model.boundary_condition.create", &json!({"name":"plane","type":"displacement","target":{"set":"all"},"values":[0,null,0,0,0,null]}))
        .unwrap();
    s.execute("model.load.create", &json!({"name":"force","type":"force","target":{"set":"tip"},"vector":[0,0.0001,0]})).unwrap();
    s.execute("model.step.create",&json!({"name":"Impact","procedure":{"type":"transient","time_step":0.01,"increments":1000,"output_every":10,"method":{"type":"modal","modes":4},"load_amplitudes":{"force":{"type":"tabular","points":[[0,0],[0.1,1],[0.2,0],[10,0]]}}},"boundary_conditions":["clamp","plane"],"loads":["force"],"output":["displacement"]})).unwrap();
    s.execute("solve.transient", &json!({"threads":1})).unwrap();
    let r = s.execute("solve.probe", &json!({"step":"Impact","nodes":[n+1]})).unwrap();
    let roots: [f64; 4] = [1.875104068711961, 4.694091132974175, 7.854757438237612, 10.995540734875467];
    let mut error: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for a in r["nodes"][0]["samples"].as_array().unwrap() {
        let t = a["time"].as_f64().unwrap();
        let exact: f64 = roots
            .iter()
            .map(|b| {
                let w = b * b * 0.01;
                let g = |x: f64| if x > 0.0 { x / (w * w) - (w * x).sin() / (w * w * w) } else { 0.0 };
                0.0004 * (g(t) - 2.0 * g(t - 0.1) + g(t - 0.2)) / 0.1
            })
            .sum();
        scale = scale.max(exact.abs());
        error = error.max((a["displacement"][1].as_f64().unwrap() - exact).abs());
    }
    eprintln!("cantilever four-mode tip impact max_error={error:e}, peak={scale:e}, relative={:e}", error / scale);
    // 1% allows finite mesh and Timoshenko/rotary-inertia departure from EB,
    // plus fixed-step dispersion. This is not a modal truncation test.
    assert!(error / scale < 0.01);
}

#[test]
fn rod_consistent_mass_and_frequency_are_explicit() {
    let truss = fieldcraft_fem::truss::Truss::new(fieldcraft_model::ElementKind::Rod2, &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]).unwrap();
    let mass = truss.mass(3.0, 1.0).unwrap();
    let lumped = truss.lumped_mass(3.0, 1.0).unwrap();
    let stiffness = truss.stiffness(4.0, 1.0).unwrap();
    for (i, j, want) in [(0, 0, 1.0), (0, 1, 0.5), (1, 0, 0.5), (1, 1, 1.0)] {
        assert!((mass.get(i, j) - want).abs() < 1e-14);
    }
    assert!((lumped[0] - 1.5).abs() < 1e-14 && (lumped[1] - 1.5).abs() < 1e-14);
    let mut s = rod();
    s.execute("model.step.edit", &json!({"name":"Pulse","loads":[],"procedure":{"type":"frequency","modes":1}})).unwrap();
    let modes = fieldcraft_solver::solve_frequency(&s.model, Some("Pulse"), &Default::default(), &Default::default()).unwrap();
    let f = modes.modes[0].frequency;
    assert!((f - 1.0 / std::f64::consts::PI).abs() < 1e-12);
    eprintln!(
        "ROD_MASS_AUDIT {}",
        json!({"total_mass_kg":3,"matrix_per_translation_kg":[[mass.get(0,0),mass.get(0,1)],[mass.get(1,0),mass.get(1,1)]],"free_dof":"node2 U1","free_mass_kg":mass.get(1,1),"free_stiffness_N_per_m":stiffness.get(3,3),"actual_frequency_hz":f,"alternative_lumped_mass_kg":lumped,"actual_transient_mass_kind":"consistent","dt_s":0.001,"increments":1000,"beta":0.25,"gamma":0.5,"HHT_alpha":0,"F0plus_N":1,"u0_m":0,"v0_m_per_s":0,"a0_m_per_s2":1})
    );
}

#[test]
fn two_close_modes_keep_signed_components_through_spectrum_commands() {
    // Original two-DOF rotated spring model; frequencies are those of Wilson's
    // published close pair (Table 15.1), but no building geometry/deck is copied.
    let mut s = Session::default();
    s.execute("model.units", &json!({"system":"m-kg-s-N-Pa"})).unwrap();
    s.execute("mesh.node.add", &json!({"label":1,"xyz":[0,0,0]})).unwrap();
    s.execute("model.set.create", &json!({"name":"all","nodes":[1]})).unwrap();
    s.execute("model.mass.create", &json!({"name":"M","target":{"set":"all"},"mass":1})).unwrap();
    s.execute(
        "model.boundary_condition.create",
        &json!({"name":"plane","type":"displacement","target":{"set":"all"},"values":[null,null,0,null,null,null]}),
    )
    .unwrap();
    let frequencies: [f64; 2] = [13.87, 13.93];
    for (i, w) in frequencies.iter().enumerate() {
        s.execute(
            "model.spring.create",
            &json!({"name":format!("K{i}"),"type":"spring","first":{"set":"all"},"direction":{"vector":[if i==0 {1} else {-1},1,0]},"stiffness":w*w}),
        )
        .unwrap();
    }
    let r = frequencies[0] / frequencies[1];
    let z: f64 = 0.05;
    let rho = 8.0 * z * z * (1.0 + r) * r.powf(1.5) / ((1.0 - r * r).powi(2) + 4.0 * z * z * r * (1.0 + r).powi(2));
    let a = 0.5 / frequencies[0].powi(2);
    let b = 0.5 / frequencies[1].powi(2);
    let mut p =
        json!({"type":"response_spectrum","modes":2,"damping_ratio":z,"combination":"cqc","direction":[1,0,0],"pseudo_acceleration":[[1,1],[5,1]]});
    s.execute("model.step.create", &json!({"name":"Spectrum","procedure":p,"boundary_conditions":["plane"],"output":["displacement"]})).unwrap();
    for cqc in [true, false] {
        p["combination"] = json!(if cqc { "cqc" } else { "srss" });
        s.execute("model.step.edit", &json!({"name":"Spectrum","procedure":p})).unwrap();
        s.execute("solve.spectrum", &json!({"threads":1})).unwrap();
        let response = s.execute("solve.probe", &json!({"step":"Spectrum","nodes":[1]})).unwrap();
        let u = &response["nodes"][0]["samples"][0]["displacement"];
        let cross = if cqc { 2.0 * rho * a * b } else { 0.0 };
        assert!((u[0].as_f64().unwrap() - (a * a + b * b + cross).sqrt()).abs() < 1e-10);
        assert!((u[1].as_f64().unwrap() - (a * a + b * b - cross).sqrt()).abs() < 1e-10);
    }
}

#[test]
fn nonzero_initial_state_round_trip_and_missing_node_refusal() {
    let mut s = rod();
    let procedure = json!({"type":"transient","time_step":0.01,"increments":10,"output_every":1,
        "integration":{"type":"hht","alpha":-0.1},
        "initial_displacement":{"2":[0.001,0,0]},"initial_velocity":{"2":[0.002,0,0]}});
    s.execute("model.step.edit", &json!({"name":"Pulse","loads":[],"procedure":procedure})).unwrap();
    s.execute("solve.transient", &json!({"step":"Pulse","threads":1})).unwrap();
    let before = probe(&mut s);
    let initial = &before["nodes"][0]["samples"][0];
    assert_eq!(initial["displacement"][0], 0.001);
    assert_eq!(initial["velocity"][0], 0.002);
    assert!((initial["acceleration"][0].as_f64().unwrap() + 0.004).abs() < 1e-12);
    let path = std::env::temp_dir().join(format!("fieldcraft-initial-roundtrip-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut restored = Session::default();
    restored.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(probe(&mut restored), before);
    restored.execute("solve.transient", &json!({"step":"Pulse","threads":1})).unwrap();
    assert_eq!(probe(&mut restored), before);
    let mut missing = procedure;
    missing["initial_displacement"] = json!({"999":[0.001,0,0]});
    restored.execute("model.step.edit", &json!({"name":"Pulse","procedure":missing})).unwrap();
    let error = restored.execute("solve.transient", &json!({"step":"Pulse","threads":1})).unwrap_err();
    assert!(error.to_string().contains("initial node 999 not present"), "{error}");
    let _ = std::fs::remove_file(path);
}

#[test]
fn direct_energy_ledger_closes_and_survives_save_open() {
    let mut s = rod();
    let mut procedure = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    procedure["damping"] = json!({"alpha":0.2,"beta":0.0,"loss_factor":0.0});
    s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).unwrap();
    s.execute("solve.transient", &json!({"threads":1})).unwrap();
    let before = probe(&mut s);
    let mut previous_d = 0.0;
    let mut theory_error = 0.0_f64;
    for (frame, sample) in before["frames"].as_array().unwrap().iter().zip(before["nodes"][0]["samples"].as_array().unwrap()) {
        let work = frame["external_work"].as_f64().unwrap();
        let loss = frame["damping_dissipation"].as_f64().unwrap();
        let u = sample["displacement"][0].as_f64().unwrap();
        assert!((work - u).abs() < 1e-12); // constant 1 N, zero u0
        assert!(loss >= previous_d && frame["energy_balance"].as_f64().unwrap().abs() < 1e-12);
        previous_d = loss;
        let time = sample["time"].as_f64().unwrap();
        let omega = 3.99_f64.sqrt();
        let decay = (-0.1 * time).exp();
        let u_exact = 0.25 * (1.0 - decay * ((omega * time).cos() + 0.1 / omega * (omega * time).sin()));
        let v_exact = decay / omega * (omega * time).sin();
        let exact = [0.5 * v_exact * v_exact, 2.0 * u_exact * u_exact, u_exact - 0.5 * v_exact * v_exact - 2.0 * u_exact * u_exact, u_exact];
        for (key, expected) in ["kinetic_energy", "strain_energy", "damping_dissipation", "external_work"].iter().zip(exact) {
            theory_error = theory_error.max((frame[*key].as_f64().unwrap() - expected).abs());
        }
    }
    assert!(previous_d > 0.01);
    assert!(theory_error < 1e-6);
    println!("ENERGY_THEORY maximum_error_J={theory_error:.16e}, tolerance_J=1e-6");
    let maximum = before["frames"].as_array().unwrap().iter().map(|f| f["energy_balance"].as_f64().unwrap().abs()).fold(0.0_f64, f64::max);
    println!("ENERGY_AUDIT max_balance_J={maximum:.16e}, tolerance_J=1e-12, final_frame={}", before["frames"].as_array().unwrap().last().unwrap());
    let path = std::env::temp_dir().join(format!("fieldcraft-energy-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut restored = Session::default();
    restored.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(probe(&mut restored), before);
    restored.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    let converted = probe(&mut restored);
    // kN mm = J: the ledger uses the same energy conversion as kinetic energy.
    assert_eq!(
        converted["frames"].as_array().unwrap().last().unwrap()["external_work"],
        before["frames"].as_array().unwrap().last().unwrap()["external_work"]
    );
    let _ = std::fs::remove_file(path);
}
