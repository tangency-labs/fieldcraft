//! Prestressed tangent dynamics: physical discrete references and frozen interfaces.
use fieldcraft_engine::{Session, StepResult};
use fieldcraft_solver::StepResult as Implicit;
use serde_json::{Value, json};
fn run(s: &mut Session, c: &str, p: Value) -> Value {
    s.execute(c, &p).unwrap()
}
fn response(s: &Session, n: &str) -> fieldcraft_solver::transient::Response {
    let StepResult::Implicit(Implicit::Dynamics(r)) = s.result(n).unwrap().as_ref() else { panic!("dynamics") };
    r.clone()
}
fn roundtrip(s: &mut Session, tag: &str, script: &Value) {
    let root = std::env::var_os("FIELDCRAFT_BASE_STATE_EVIDENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join(format!("based-{}", std::process::id())));
    std::fs::create_dir_all(&root).unwrap();
    let path = root.join(format!("{tag}.fieldcraft"));
    run(s, "file.save_as", json!({"path":path}));
    let mut read = Session::default();
    run(&mut read, "file.open", json!({"path":path}));
    for name in ["Direct", "Modal", "Zero"] {
        assert_eq!(s.result(name), read.result(name));
    }
    if std::env::var_os("FIELDCRAFT_BASE_STATE_EVIDENCE").is_some() {
        std::fs::write(root.join(format!("{tag}-input.json")), serde_json::to_vec_pretty(script).unwrap()).unwrap();
        for name in ["Direct", "Modal", "Zero"] {
            std::fs::write(root.join(format!("{tag}-{name}.json")), serde_json::to_vec_pretty(&response(s, name)).unwrap()).unwrap();
        }
        for name in ["Preload", "Locked", "Modes"] {
            if let Some(StepResult::Implicit(result)) = s.result(name).map(|r| r.as_ref()) {
                let data = match result {
                    Implicit::Static(r) => serde_json::to_vec_pretty(r).unwrap(),
                    Implicit::Frequency(r) => serde_json::to_vec_pretty(r).unwrap(),
                    _ => continue,
                };
                assert_eq!(s.result(name), read.result(name));
                std::fs::write(root.join(format!("{tag}-{name}.json")), data).unwrap();
            }
        }
    } else {
        std::fs::remove_file(path).unwrap();
    }
}
#[test]
fn tensioned_two_rods_direct_modal_and_zero_perturbation() {
    for mass in ["consistent", "lumped"] {
        let mut commands = vec![json!({"command":"model.units","params":{"system":"m-kg-s-N-Pa"}})];
        let mut add = |command: &str, params: Value| commands.push(json!({"command":command,"params":params}));
        for (label, x) in [(1, 0.), (2, 1.), (3, 2.)] {
            add("mesh.node.add", json!({"label":label,"xyz":[x,0,0]}));
        }
        for (label, nodes) in [(1, [1, 2]), (2, [2, 3])] {
            add("mesh.element.add", json!({"label":label,"kind":"rod2","nodes":nodes}));
        }
        for (name, nodes) in [("left", vec![1]), ("mid", vec![2]), ("right", vec![3])] {
            add("model.set.create", json!({"name":name,"nodes":nodes}));
        }
        add("model.set.create", json!({"name":"rods","elements":[1,2]}));
        add("model.material.create", json!({"name":"E","elastic":{"type":"isotropic","young":1000.,"poisson":0.},"density":3.}));
        add("model.section.create", json!({"name":"A","set":"rods","material":"E","kind":{"type":"rod","area":1.}}));
        add("model.boundary_condition.create", json!({"name":"left","type":"displacement","target":{"set":"left"},"values":[0,0,0,null,null,null]}));
        add(
            "model.boundary_condition.create",
            json!({"name":"right","type":"displacement","target":{"set":"right"},"values":[0.02,0,0,null,null,null]}),
        );
        add(
            "model.boundary_condition.create",
            json!({"name":"plane","type":"displacement","target":{"set":"mid"},"values":[null,null,0,null,null,null]}),
        );
        // Explicit physical weak transverse support makes the unloaded equilibrium nonsingular.
        // Its 1 N/m is retained in both FE and independent reference, never hidden regularization.
        add("model.spring.create", json!({"name":"lateral","type":"spring","stiffness":1.,"first":{"set":"mid"},"direction":{"dof":2}}));
        add("model.step.create", json!({"name":"Preload","procedure":{"type":"static"},"boundary_conditions":["left","right","plane"]}));
        add("solve.static", json!({"step":"Preload"}));
        for (name, modal, zero) in [("Direct", false, false), ("Modal", true, false), ("Zero", false, true)] {
            add(
                "model.step.create",
                json!({"name":name,"base_state":"Preload","procedure":{"type":"transient","time_step":0.001,"increments":100,"integration":{"type":"average_acceleration"},"method":if modal{json!({"type":"modal","modes":2})}else{json!({"type":"direct"})},"initial_displacement":if zero{json!({})}else{json!({"2":[0,0.001,0]})}},"boundary_conditions":["left","right","plane"],"output":["displacement"]}),
            );
            add("solve.transient", json!({"step":name,"mass":mass}));
        }
        let script = json!(commands);
        let mut s = Session::default();
        s.run_script(&script).unwrap();
        let d = response(&s, "Direct");
        let m = response(&s, "Modal");
        let z = response(&s, "Zero");
        let mass_value: f64 = if mass == "consistent" { 2. } else { 3. };
        let w = (21. / mass_value).sqrt();
        assert!((m.natural_frequencies[0] - w / (2. * std::f64::consts::PI)).abs() < 1e-10);
        for ((a, b), zero) in d.frames.iter().zip(&m.frames).zip(&z.frames) {
            let t = a.time.unwrap();
            assert!((a.displacement[1][1] - 0.001 * (w * t).cos()).abs() < 1e-9);
            for (u, v) in a.displacement.iter().flatten().zip(b.displacement.iter().flatten()) {
                assert!((u - v).abs() < 1e-12);
            }
            assert!(zero.displacement.iter().flatten().all(|x| *x == 0.));
            assert!(zero.velocity.iter().flatten().all(|x| *x == 0.));
        }
        let energy = d.energy.as_ref().unwrap();
        assert!(energy.iter().all(|e| e.energy_balance.abs() < 1e-14));
        for f in &d.frames {
            assert!((f.acceleration[1][1] + w * w * f.displacement[1][1]).abs() < 1e-11);
        }
        roundtrip(&mut s, &format!("string-{mass}"), &script);
        // Saved dependency must stale the dynamic result when its preload changes.
        let definition = fieldcraft_io::results::definition_hash(&s.model, "Direct", &Default::default()).unwrap();
        run(
            &mut s,
            "model.boundary_condition.edit",
            json!({"name":"right","type":"displacement","target":{"set":"right"},"values":[0.04,0,0,null,null,null]}),
        );
        let changed = fieldcraft_io::results::definition_hash(&s.model, "Direct", &Default::default()).unwrap();
        assert_ne!(definition, changed);
        assert!(s.result("Direct").is_none());
    }
}
#[test]
fn clamped_bolt_frozen_contact_direct_modal_and_zero() {
    let mut commands: Vec<Value> = serde_json::from_str(include_str!("fixtures/pretension.json")).unwrap();
    commands.retain(|c| c["command"] != "solve.harmonic");
    let mut s = Session::default();
    s.run_script(&json!(commands)).unwrap();
    run(&mut s, "model.step.edit", json!({"name":"Modes","procedure":{"type":"frequency","modes":100}}));
    run(&mut s, "solve.frequency", json!({"step":"Modes"}));
    commands.push(json!({"command":"model.step.edit","params":{"name":"Modes","procedure":{"type":"frequency","modes":100}}}));
    commands.push(json!({"command":"solve.frequency","params":{"step":"Modes"}}));
    let StepResult::Implicit(Implicit::Frequency(f)) = s.result("Modes").unwrap().as_ref() else { panic!("frequency") };
    let f = f.clone();
    for (name, modal, zero) in [("Direct", false, false), ("Modal", true, false), ("Zero", false, true)] {
        let p = json!({"name":name,"base_state":"Locked","procedure":{"type":"transient","time_step":1e-5,"increments":20,"method":if modal{json!({"type":"modal","modes":100})}else{json!({"type":"direct"})}},"boundary_conditions":["Bottom","CapGuide"],"loads":if zero{json!([])}else{json!(["Shake"])},"output":["displacement"]});
        run(&mut s, "model.step.create", p.clone());
        commands.push(json!({"command":"model.step.create","params":p}));
        run(&mut s, "solve.transient", json!({"step":name}));
        commands.push(json!({"command":"solve.transient","params":{"step":name}}));
    }
    let d = response(&s, "Direct");
    let m = response(&s, "Modal");
    let z = response(&s, "Zero");
    let cap = f.reference_points.iter().position(|(n, _)| n == "Cap").unwrap();
    for (i, ((a, b), zero)) in d.frames.iter().zip(&m.frames).zip(&z.frames).enumerate() {
        let mut scale = 0f64;
        let mut err = 0f64;
        for n in 0..d.nodes.len() {
            for c in 0..3 {
                let expected: f64 = f
                    .modes
                    .iter()
                    .map(|mode| {
                        let theta = 2. * (mode.eigenvalue.sqrt() * 1e-5 / 2.).atan();
                        mode.shape[n][c] * mode.point_shape[cap][2] / mode.eigenvalue * (1. - (i as f64 * theta).cos())
                    })
                    .sum();
                scale = scale.max(expected.abs());
                err = err.max((a.displacement[n][c] - expected).abs());
                assert!((a.displacement[n][c] - b.displacement[n][c]).abs() < 1e-10);
            }
        }
        assert!(err < 1e-12 + scale * 1e-5, "{err}/{scale}");
        assert!(zero.displacement.iter().flatten().all(|x| *x == 0.));
        for (a, b) in [(5, 9), (6, 10), (7, 11), (8, 12), (21, 25), (22, 26), (23, 27), (24, 28)] {
            let ia = d.nodes.iter().position(|n| *n == a).unwrap();
            let ib = d.nodes.iter().position(|n| *n == b).unwrap();
            for c in 0..3 {
                assert!((d.frames[i].displacement[ia][c] - d.frames[i].displacement[ib][c]).abs() < 1e-12);
            }
        }
    }
    roundtrip(&mut s, "bolt", &json!(commands));
}

#[test]
fn axial_rod_base_nullspace_is_not_a_transverse_support() {
    for force in [1000_f64, 10000.] {
        let mut rows = Vec::new();
        let mut add = |command: &str, params: Value| rows.push(json!({"command":command,"params":params}));
        add("model.units", json!({"system":"m-kg-s-N-Pa"}));
        for (label, x) in [(1, 0.), (2, 1.)] {
            add("mesh.node.add", json!({"label":label,"xyz":[x,0,0]}));
        }
        add("mesh.element.add", json!({"label":1,"kind":"rod2","nodes":[1,2]}));
        for (name, nodes) in [("fix", vec![1]), ("tip", vec![2])] {
            add("model.set.create", json!({"name":name,"nodes":nodes}));
        }
        add("model.set.create", json!({"name":"rod","elements":[1]}));
        add("model.material.create", json!({"name":"E","elastic":{"type":"isotropic","young":2.1e11,"poisson":0.3},"density":7850.}));
        add("model.section.create", json!({"name":"A","set":"rod","material":"E","kind":{"type":"rod","area":1e-4}}));
        add("model.boundary_condition.create", json!({"name":"fixed","type":"fixed","target":{"set":"fix"}}));
        add("model.load.create", json!({"name":"pull","type":"force","target":{"set":"tip"},"vector":[force,0,0]}));
        add("model.step.create", json!({"name":"Preload","procedure":{"type":"static"},"boundary_conditions":["fixed"],"loads":["pull"]}));
        add(
            "model.step.create",
            json!({"name":"Modes","base_state":"Preload","procedure":{"type":"frequency","modes":3},"boundary_conditions":["fixed"]}),
        );
        add("solve.static", json!({"step":"Preload"}));
        add("solve.frequency", json!({"step":"Modes","mass":"lumped"}));
        for (name, modal, zero) in [("Direct", false, false), ("Modal", true, false), ("Zero", false, true)] {
            add(
                "model.step.create",
                json!({"name":name,"base_state":"Preload","procedure":{"type":"transient","time_step":1e-5,"increments":100,"method":if modal{json!({"type":"modal","modes":3})}else{json!({"type":"direct"})},"initial_displacement":if zero{json!({})}else{json!({"2":[0,1e-5,0]})}},"boundary_conditions":["fixed"],"output":["displacement"]}),
            );
            add("solve.transient", json!({"step":name,"mass":"lumped"}));
        }
        let script = json!(rows);
        let mut s = Session::default();
        s.run_script(&script).unwrap();
        let StepResult::Implicit(Implicit::Frequency(f)) = s.result("Modes").unwrap().as_ref() else { panic!() };
        let w = (force / 0.3925).sqrt();
        assert!((f.modes[0].frequency - w / (2. * std::f64::consts::PI)).abs() < 1e-8);
        let d = response(&s, "Direct");
        let m = response(&s, "Modal");
        let z = response(&s, "Zero");
        assert!((m.natural_frequencies[0] - f.modes[0].frequency).abs() < 1e-8);
        for ((a, b), zero) in d.frames.iter().zip(&m.frames).zip(&z.frames) {
            assert!((a.displacement[1][1] - 1e-5 * (w * a.time.unwrap()).cos()).abs() < 1e-10);
            assert!((a.displacement[1][1] - b.displacement[1][1]).abs() < 1e-12);
            assert!(zero.displacement.iter().flatten().all(|v| *v == 0.));
        }
        roundtrip(&mut s, &format!("axial-{force}"), &script);
        // A force in the missing transverse stiffness is an actual unsolved mechanism.
        run(&mut s, "model.load.edit", json!({"name":"pull","type":"force","target":{"set":"tip"},"vector":[force,1,0]}));
        assert!(
            s.execute("solve.transient", &json!({"step":"Direct","mass":"lumped"})).unwrap_err().to_string().contains("zero-stiffness direction")
        );
    }
}

/// Same total-Lagrangian force law as geometric dynamics, independently differentiated
/// here; the equilibrium is obtained through commands, never supplied analytically.
#[test]
fn finite_rod_equilibrium_tangent_and_dynamic_chain() {
    for force in [0_f64, 1000., 10000.] {
        for mass in ["consistent", "lumped"] {
            let mut rows = Vec::new();
            let mut add = |command: &str, params: Value| rows.push(json!({"command":command,"params":params}));
            add("model.units", json!({"system":"m-kg-s-N-Pa"}));
            for (label, x) in [(1, 0.), (2, 1.)] {
                add("mesh.node.add", json!({"label":label,"xyz":[x,0,0]}));
            }
            add("mesh.element.add", json!({"label":1,"kind":"rod2","nodes":[1,2]}));
            for (name, nodes) in [("fix", vec![1]), ("tip", vec![2])] {
                add("model.set.create", json!({"name":name,"nodes":nodes}));
            }
            add("model.set.create", json!({"name":"rod","elements":[1]}));
            add("model.material.create", json!({"name":"E","elastic":{"type":"isotropic","young":2.1e11,"poisson":0.3},"density":7850.}));
            add("model.section.create", json!({"name":"A","set":"rod","material":"E","kind":{"type":"rod","area":1e-4}}));
            add("model.boundary_condition.create", json!({"name":"fixed","type":"fixed","target":{"set":"fix"}}));
            if force != 0. {
                add("model.load.create", json!({"name":"pull","type":"force","target":{"set":"tip"},"vector":[force,0,0]}));
            }
            add(
                "model.step.create",
                json!({"name":"Preload","procedure":{"type":"nonlinear_static","nonlinear_geometry":true},"boundary_conditions":["fixed"],"loads":if force==0. {json!([])} else {json!(["pull"])}}),
            );
            add("solve.nonlinear", json!({"step":"Preload"}));
            add(
                "model.step.create",
                json!({"name":"Modes","base_state":"Preload","procedure":{"type":"frequency","modes":3},"boundary_conditions":["fixed"]}),
            );
            add("solve.frequency", json!({"step":"Modes","mass":mass}));
            for (name, modal, zero) in [("Direct", false, false), ("Modal", true, false), ("Zero", false, true)] {
                add(
                    "model.step.create",
                    json!({"name":name,"base_state":"Preload","procedure":{"type":"transient","time_step":1e-4,"increments":100,"method":if modal{json!({"type":"modal","modes":3})}else{json!({"type":"direct"})},"initial_displacement":if zero{json!({})}else{json!({"2":[0,1e-5,0]})}},"boundary_conditions":["fixed"],"output":["displacement"]}),
                );
                add("solve.transient", json!({"step":name,"mass":mass}));
            }
            let script = json!(rows);
            let mut s = Session::default();
            s.run_script(&script).unwrap();
            let StepResult::Implicit(Implicit::Static(base)) = s.result("Preload").unwrap().as_ref() else { panic!() };
            let ea = 21e6;
            let stretch = 1. + base.displacement[1][0];
            let axial = 0.5 * ea * stretch * (stretch * stretch - 1.);
            let ue = [0., 0., 0., base.displacement[1][0], base.displacement[1][1], base.displacement[1][2]];
            let recovered = fieldcraft_solver::transient::nonlinear::recover_rod(&[[0., 0., 0.], [1., 0., 0.]], &ue, 2.1e11, 1e-4).unwrap();
            assert!((recovered.nodal_internal_force[3] - force).abs() < 1e-6);
            if let Some(root) = std::env::var_os("FIELDCRAFT_BASE_STATE_EVIDENCE") {
                let root = std::path::PathBuf::from(root);
                std::fs::create_dir_all(&root).unwrap();
                let value = json!({"source":"production recover_rod evaluated at actual solve.nonlinear saved state", "recovery":recovered,"external_force":[0.,0.,0.,force,0.,0.],"free_residual":[recovered.nodal_internal_force[3]-force,recovered.nodal_internal_force[4],recovered.nodal_internal_force[5]],"reaction":base.reaction,"displacement":base.displacement});
                std::fs::write(
                    root.join(format!("finite-axial-{force}-{mass}-equilibrium-diagnostic.json")),
                    serde_json::to_vec_pretty(&value).unwrap(),
                )
                .unwrap();
            }

            assert!((axial - force).abs() < 1e-6, "force {force}, residual {}", axial - force);
            assert!((base.reaction[0][0] + force).abs() < 1e-6);
            assert!((base.elements[0].points[0].stress[0] - 2.1e11 * (stretch * stretch - 1.) / 2.).abs() < 1e-3);
            assert!((base.elements[0].points[0].strain[0] - (stretch * stretch - 1.) / 2.).abs() < 1e-14);
            assert!(base.nonlinear.as_ref().unwrap().increments.iter().all(|i| i.converged));
            // Scalar independent equilibrium root and exact material + stress tangent.
            let mut root = 1.;
            for _ in 0..12 {
                root -= (0.5 * ea * root * (root * root - 1.) - force) / (0.5 * ea * (3. * root * root - 1.));
            }
            assert!((stretch - root).abs() < 1e-12);
            let expected = [0.5 * ea * (3. * stretch * stretch - 1.), 0.5 * ea * (stretch * stretch - 1.), 0.5 * ea * (stretch * stretch - 1.)];
            let modal_mass = 0.785 / if mass == "lumped" { 2. } else { 3. };
            let StepResult::Implicit(Implicit::Frequency(freq)) = s.result("Modes").unwrap().as_ref() else { panic!() };
            assert_eq!(freq.modes.len(), 3);
            let mut actual = [[0.; 3]; 3];
            for (i, row) in actual.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    *value = modal_mass * modal_mass * freq.modes.iter().map(|m| m.eigenvalue * m.shape[1][i] * m.shape[1][j]).sum::<f64>();
                }
            }
            let law = |d: [f64; 3]| {
                let g = 0.5 * (d.iter().map(|x| x * x).sum::<f64>() - 1.);
                d.map(|x| ea * g * x)
            };
            for h in [1e-5, 1e-6, 1e-7] {
                for j in 0..3 {
                    let mut p = [stretch, 0., 0.];
                    let mut m = p;
                    p[j] += h;
                    m[j] -= h;
                    let fp = law(p);
                    let fm = law(m);
                    for i in 0..3 {
                        let fd = (fp[i] - fm[i]) / (2. * h);
                        assert!((actual[i][j] - fd).abs() < 0.1, "P={force}, ({i},{j}) {} vs {fd}", actual[i][j]);
                    }
                }
            }
            for i in 0..3 {
                assert!(
                    (actual[i][i] - expected[i]).abs() < 1e-9 * expected[0],
                    "force={force} mass={mass} actual={actual:?} expected={expected:?} frequencies={:?}",
                    freq.modes.iter().map(|m| m.frequency).collect::<Vec<_>>()
                );
            }
            let w = (expected[1] / modal_mass).sqrt();
            assert!((freq.modes[0].frequency - w / (2. * std::f64::consts::PI)).abs() < 1e-6);
            let d = response(&s, "Direct");
            let m = response(&s, "Modal");
            let z = response(&s, "Zero");
            let theta = 2. * (w * 1e-4 / 2.).atan();
            for (index, ((a, b), zero)) in d.frames.iter().zip(&m.frames).zip(&z.frames).enumerate() {
                let expected_u = 1e-5 * (index as f64 * theta).cos();
                assert!((a.displacement[1][1] - expected_u).abs() < 1e-11);
                for (c, stiffness) in expected.iter().enumerate() {
                    assert!((a.displacement[1][c] - b.displacement[1][c]).abs() < 1e-12);
                    assert!((a.acceleration[1][c] - b.acceleration[1][c]).abs() < 1e-8);
                    assert!((modal_mass * a.acceleration[1][c] + stiffness * a.displacement[1][c]).abs() < 1e-7);
                }
                assert!(zero.displacement.iter().chain(&zero.velocity).chain(&zero.acceleration).flatten().all(|v| *v == 0.));
            }
            roundtrip(&mut s, &format!("finite-axial-{force}-{mass}"), &script);
            if force == 0. {
                run(&mut s, "model.load.create", json!({"name":"pull","type":"force","target":{"set":"tip"},"vector":[0,1,0]}));
                run(&mut s, "model.step.edit", json!({"name":"Preload","loads":["pull"]}));
            }
            run(&mut s, "model.load.edit", json!({"name":"pull","type":"force","target":{"set":"tip"},"vector":[force,1,0]}));
            assert!(s.execute("solve.nonlinear", &json!({"step":"Preload"})).unwrap_err().to_string().contains("zero-stiffness direction"));
        }
    }
}
