//! Original finite-element models, compared with independently integrated Cartesian ODEs.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn run(s: &mut Session, c: &str, p: Value) -> Value {
    s.execute(c, &p).unwrap()
}
fn model(string: bool, dt: f64, alpha: Option<f64>, mass: &str, damping: f64) -> (Session, Value) {
    let mut s = Session::default();
    let mut commands = Vec::new();
    let mut add = |command: &str, params: Value| commands.push(json!({"command":command,"params":params}));
    add("model.units", json!({"system":"m-kg-s-N-Pa"}));
    let xyz = if string { vec![[-1., 0., 0.], [0., 0., 0.], [1., 0., 0.]] } else { vec![[0., 0., 0.], [0., -1., 0.]] };
    for (i, x) in xyz.iter().enumerate() {
        add("mesh.node.add", json!({"label":i+1,"xyz":x}));
    }
    add("mesh.element.add", json!({"label":1,"kind":"rod2","nodes":[1,2]}));
    if string {
        add("mesh.element.add", json!({"label":2,"kind":"rod2","nodes":[2,3]}));
    }
    add("model.set.create", json!({"name":"members","elements":if string{vec![1,2]}else{vec![1]}}));
    add("model.set.create", json!({"name":"fixed","nodes":if string{vec![1,3]}else{vec![1]}}));
    add("model.set.create", json!({"name":"tip","nodes":[2]}));
    add("model.material.create", json!({"name":"elastic","elastic":{"type":"isotropic","young":1000.,"poisson":0.25},"density":3.}));
    add("model.section.create", json!({"name":"rod","set":"members","material":"elastic","kind":{"type":"rod","area":1.}}));
    add("model.boundary_condition.create", json!({"name":"fixed","type":"displacement","target":{"set":"fixed"},"values":[0,0,0,null,null,null]}));
    add(
        "model.boundary_condition.create",
        json!({"name":"plane","type":"displacement","target":{"set":"tip"},"values":[if string{json!(0)}else{Value::Null},null,0,null,null,null]}),
    );
    add("model.load.create", json!({"name":"force","type":"force","target":{"set":"tip"},"vector":[0.,if string{8.}else{-9.81},0.]}));
    add(
        "model.step.create",
        json!({"name":"Motion","procedure":{"type":"transient","time_step":dt,"increments":(1./dt).round() as u32,"output_every":1,"integration":alpha.map_or(json!({"type":"average_acceleration"}),|a|json!({"type":"hht","alpha":a})),"initial_displacement":{"2":if string{vec![0.,0.25,0.]}else{vec![(3f64).sqrt()/2.,0.5,0.]}},"damping":{"alpha":damping}},"loads":["force"],"boundary_conditions":["fixed","plane"],"output":["displacement"]}),
    );
    let input = json!(commands);
    s.run_script(&input).unwrap();
    let _ = mass;
    (s, input)
}
fn derivative(y: [f64; 4], string: bool, m: f64, c: f64) -> [f64; 4] {
    let [x, z, vx, vz] = y;
    let (fx, fz) = if string {
        (0., 1000. * z.powi(3))
    } else {
        let strain = (x * x + z * z - 1.) / 2.;
        (1000. * strain * x, 1000. * strain * z)
    };
    [vx, vz, (-fx - c * m * vx) / m, ((if string { 8. } else { -9.81 }) - fz - c * m * vz) / m]
}
fn reference(string: bool, m: f64, c: f64, dt: f64, steps: usize) -> Vec<[f64; 4]> {
    let mut y = if string { [0., 0.25, 0., 0.] } else { [(3f64).sqrt() / 2., -0.5, 0., 0.] };
    let mut out = vec![y];
    let h = dt / 100.;
    for _ in 0..steps {
        for _ in 0..100 {
            let a = derivative(y, string, m, c);
            let b = derivative(std::array::from_fn(|i| y[i] + h * a[i] / 2.), string, m, c);
            let d = derivative(std::array::from_fn(|i| y[i] + h * b[i] / 2.), string, m, c);
            let e = derivative(std::array::from_fn(|i| y[i] + h * d[i]), string, m, c);
            for i in 0..4 {
                y[i] += h * (a[i] + 2. * b[i] + 2. * d[i] + e[i]) / 6.;
            }
        }
        out.push(y);
    }
    out
}
#[test]
#[ignore = "heavy: 48 previously verified nonlinear convergence cases"]
fn large_rotation_pendulum_and_tensioned_string_converge_with_energy_accounts() {
    let mut reports = Vec::new();
    for string in [false, true] {
        for mass in ["consistent", "lumped"] {
            for alpha in [None, Some(-0.1)] {
                for damping in [0., 0.2] {
                    let mut errors = Vec::new();
                    for dt in [0.01, 0.005, 0.0025] {
                        let (mut s, input) = model(string, dt, alpha, mass, damping);
                        let out = run(&mut s, "solve.transient.nonlinear", json!({"step":"Motion","mass":mass}));
                        let m = (if mass == "consistent" { 1. } else { 1.5 }) * (if string { 2. } else { 1. });
                        let reference = reference(string, m, damping, dt, (1. / dt).round() as usize);
                        let frames = out["frames"].as_array().unwrap();
                        assert_eq!(frames.len(), reference.len());
                        let mut error = [0f64; 3];
                        let mut reference_peak = [0f64; 3];
                        let mut balance = 0f64;
                        let mut residual = 0f64;
                        for (f, y) in frames.iter().zip(&reference) {
                            let acc = derivative(*y, string, m, damping);
                            for k in 0..2 {
                                for (q, (key, expected)) in [
                                    ("displacement", y[k] + if !string && k == 1 { 1. } else { 0. }),
                                    ("velocity", y[k + 2]),
                                    ("acceleration", acc[k + 2]),
                                ]
                                .iter()
                                .enumerate()
                                {
                                    reference_peak[q] = reference_peak[q].max(expected.abs());
                                    error[q] = error[q].max((f[key][1][k].as_f64().unwrap() - expected).abs());
                                }
                            }
                            balance = balance.max(f["energy_balance"].as_f64().unwrap().abs());
                            residual = residual.max(f["relative_residual"].as_f64().unwrap());
                        }
                        if dt == 0.0025 {
                            for q in 0..3 {
                                assert!(error[q] < 0.01 * reference_peak[q], "component {q}: {error:?}, peaks {reference_peak:?}");
                            }
                            let escale = frames
                                .iter()
                                .map(|f| {
                                    f["external_work"].as_f64().unwrap().abs()
                                        + f["kinetic_energy"].as_f64().unwrap()
                                        + f["strain_energy"].as_f64().unwrap()
                                })
                                .fold(0f64, f64::max);
                            assert!(balance < 1e-3 * escale, "energy {balance} / {escale}");
                        }
                        let initial = &frames[0];
                        assert!(initial["relative_residual"].as_f64().unwrap() < 1e-10);
                        assert!(residual < 1e-8);
                        assert_eq!(initial["physical_dissipation"], 0.);
                        if damping == 0. {
                            assert_eq!(frames.last().unwrap()["physical_dissipation"], 0.);
                        } else {
                            assert!(frames.last().unwrap()["physical_dissipation"].as_f64().unwrap() > 0.);
                        }
                        let tag = format!(
                            "{}-{mass}-{}-c{damping}-dt{dt}",
                            if string { "string" } else { "pendulum" },
                            if alpha.is_some() { "hht" } else { "newmark" }
                        );
                        if let Ok(dir) = std::env::var("FIELDCRAFT_NONLINEAR_EVIDENCE") {
                            let p = std::path::Path::new(&dir).join(&tag);
                            std::fs::create_dir_all(&p).unwrap();
                            std::fs::write(p.join("input.json"), serde_json::to_vec_pretty(&input).unwrap()).unwrap();
                            std::fs::write(p.join("actual.json"), serde_json::to_vec_pretty(&out).unwrap()).unwrap();
                            std::fs::write(p.join("reference.json"), serde_json::to_vec_pretty(&reference).unwrap()).unwrap();
                        }
                        reports.push(
                            json!({"case":tag,"mass_free":m,"max_uva_error":error,"max_energy_balance":balance,"max_relative_residual":residual,"reference_peak_uva":reference_peak}),
                        );
                        errors.push(error[0]);
                    }
                    assert!(errors[1] < errors[0] * 0.4, "{errors:?}");
                    assert!(errors[2] < errors[1] * 0.4, "{errors:?}");
                    assert!(errors[2] < 0.003, "{errors:?}");
                }
            }
        }
    }
    if let Ok(dir) = std::env::var("FIELDCRAFT_NONLINEAR_EVIDENCE") {
        std::fs::write(std::path::Path::new(&dir).join("report.json"), serde_json::to_vec_pretty(&reports).unwrap()).unwrap();
    }
}
#[test]
fn failed_nonlinear_solve_preserves_previous_result_and_success_is_stored() {
    let (mut s, _) = model(false, 0.01, None, "consistent", 0.);
    run(&mut s, "solve.transient", json!({"step":"Motion"}));
    let before = run(&mut s, "solve.probe", json!({"step":"Motion","nodes":[2]}));
    assert!(s.execute("solve.transient.nonlinear", &json!({"step":"Motion","max_iterations":1})).is_err());
    assert_eq!(before, run(&mut s, "solve.probe", json!({"step":"Motion","nodes":[2]})));
    run(&mut s, "solve.transient.nonlinear", json!({"step":"Motion"}));
    assert!(matches!(s.result("Motion").unwrap().as_ref(), fieldcraft_engine::StepResult::Nonlinear(_)));
}

#[test]
fn nonfinite_energy_is_refused_before_returning_diagnostics() {
    let (_, mut input) = model(false, 0.01, None, "consistent", 0.);
    input.as_array_mut().unwrap().retain(|c| c["command"] != "model.load.create");
    for command in input.as_array_mut().unwrap() {
        if command["command"] == "model.step.create" {
            let p = &mut command["params"];
            p["boundary_conditions"] = json!([]);
            p["loads"] = json!([]);
            p["procedure"]["time_step"] = json!(1e-140);
            p["procedure"]["increments"] = json!(1);
            p["procedure"]["initial_displacement"] = json!({});
            p["procedure"]["initial_velocity"] = json!({"1":[1e155,0.,0.],"2":[1e155,0.,0.]});
        }
    }
    let mut session = Session::default();
    session.run_script(&input).unwrap();
    let result = session.execute("solve.transient.nonlinear", &json!({"step":"Motion"}));
    let error = result.expect_err("overflowing kinetic energy must not serialize as null");
    assert!(error.to_string().contains("energy overflow"), "{error}");
}

#[test]
fn native_kind11_roundtrip_stress_and_version_boundaries() {
    use fieldcraft_io::results::{self, Solved, StoredRef};
    for (mass, alpha) in [("consistent", None), ("lumped", Some(-0.1))] {
        let (_, mut input) = model(true, 0.0025, alpha, mass, 0.2);
        for command in input.as_array_mut().unwrap() {
            if command["command"] == "model.step.create" {
                command["params"]["procedure"]["increments"] = json!(4);
                command["params"]["procedure"]["output_every"] = json!(2);
            }
        }
        let mut session = Session::default();
        session.run_script(&input).unwrap();
        let reply = run(&mut session, "solve.transient.nonlinear", json!({"step":"Motion","mass":mass}));
        assert_eq!(reply["persisted"], true);
        let result = session.result("Motion").unwrap().clone();
        let fieldcraft_engine::StepResult::Nonlinear(r) = result.as_ref() else { panic!("wrong kind") };
        r.validate().unwrap();
        assert_eq!(r.frames.len(), 3);
        for f in &r.frames {
            let y = f.displacement[1][1];
            // Independent closed expression for either unit-length rod: E_G=y²/2.
            for stress in &f.rods {
                assert!((stress.green_strain - y * y / 2.).abs() < 1e-15);
                assert!((stress.second_piola_stress - 500. * y * y).abs() < 1e-12);
                assert!((stress.nominal_stress - 500. * y * y * (1. + y * y).sqrt()).abs() < 1e-12);
            }
            assert!((f.internal_force[1][1] - 1000. * y.powi(3)).abs() < 1e-11);
            assert_eq!(f.external_force[1], [0., 8., 0.]);
        }
        let bytes = results::encode(&[("Motion", StoredRef { definition: [0; 32], result: result.as_stored() })]).unwrap();
        assert_eq!(&bytes[..8], b"MCRSLT07");
        let kind_offset = 8 + 4 + 8 + "Motion".len() + 32;
        assert_eq!(bytes[kind_offset], 11);
        let results::Decoded::Results(decoded) = results::decode(&bytes).unwrap() else { panic!("missing result") };
        assert_eq!(decoded[0].result, Solved::Nonlinear(r.clone()));
        for magic in [b"MCRSLT05", b"MCRSLT06"] {
            for kind in [9, 10, 11] {
                let mut bad = bytes.clone();
                bad[..8].copy_from_slice(magic);
                bad[kind_offset] = kind;
                assert!(results::decode(&bad).unwrap_err().to_string().contains("ambiguous development"));
            }
        }
        for kind in [9, 10, 12] {
            let mut bad = bytes.clone();
            bad[kind_offset] = kind;
            assert!(results::decode(&bad).unwrap_err().to_string().contains("unknown kind"));
        }
        assert!(results::decode(&bytes[..bytes.len() - 1]).is_err());
        let mut bad = r.clone();
        bad.frames[1].rods[0].second_piola_stress += 1.;
        assert!(bad.validate().is_err());
        bad = r.clone();
        bad.frames[1].relative_residual = 1.;
        assert!(bad.validate().is_err());
        let root = std::env::var_os("FIELDCRAFT_KIND11_EVIDENCE")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| std::env::temp_dir().join(format!("fieldcraft-kind11-{}", std::process::id())));
        std::fs::create_dir_all(&root).unwrap();
        let file = root.join(format!("{mass}.fieldcraft"));
        run(&mut session, "file.save_as", json!({"path":file}));
        let archive = fieldcraft_io::read(&std::fs::read(&file).unwrap()).unwrap();
        assert_eq!(archive.results[0].result, Solved::Nonlinear(r.clone()));
        let mut opened = Session::default();
        run(&mut opened, "file.open", json!({"path":file}));
        assert_eq!(opened.result("Motion").unwrap().as_ref(), result.as_ref());
        let probe = run(&mut opened, "solve.probe", json!({"step":"Motion","nodes":[2]}));
        assert_eq!(probe["frames"], reply["frames"]);
        run(&mut opened, "model.units", json!({"system":"mm-kg-ms-kN-GPa"}));
        let converted = run(&mut opened, "solve.probe", json!({"step":"Motion","nodes":[2]}));
        assert_eq!(converted["frames"], reply["frames"]);
        assert_eq!(converted["units"], "SI");
        assert!(opened.execute("solve.probe", &json!({"step":"Motion","nodes":[2],"frame":0})).is_err());
        assert!(opened.execute("solve.probe", &json!({"step":"Motion","nodes":[999]})).is_err());
        if std::env::var_os("FIELDCRAFT_KIND11_EVIDENCE").is_some() {
            std::fs::write(root.join(format!("{mass}-input.json")), serde_json::to_vec_pretty(&input).unwrap()).unwrap();
            std::fs::write(root.join(format!("{mass}-actual.json")), serde_json::to_vec_pretty(&reply).unwrap()).unwrap();
            std::fs::write(root.join(format!("{mass}-results.bin")), bytes).unwrap();
        } else {
            std::fs::remove_file(file).unwrap();
        }
    }
}
