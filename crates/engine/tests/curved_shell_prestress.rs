//! Small-preload ring limit of a cylindrical shell: public Rayleigh ring frequency
//! (NASA MCR-68-87 vol. II, p. II-1), plus independently integrated initial-stress
//! and pressure-potential energies; see docs/solvers/curved-shell-modes.md.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
use std::f64::consts::PI;
#[path = "common/cylinder.rs"]
mod cylinder;

fn definition(tri: bool, divisions: usize) -> Value {
    let (mut recipe, _) = cylinder::model(tri, true, true, divisions);
    let commands = recipe.as_array_mut().unwrap();
    let buckle = commands.pop().unwrap();
    for command in commands.iter_mut() {
        if command["command"] == "model.material.create" {
            command["params"]["density"] = json!(1000.0);
        }
        if command["command"] == "model.load.create" && command["params"]["name"] == "pressure" {
            command["params"]["magnitude"] = json!(0.0);
        }
    }
    let supports = buckle["params"]["boundary_conditions"].clone();
    commands.extend([
        json!({"command":"model.step.create","params":{"name":"Base","procedure":{"type":"static"},"boundary_conditions":supports,"loads":["pressure"]}}),
        json!({"command":"model.step.create","params":{"name":"Plain","procedure":{"type":"frequency","modes":2},"boundary_conditions":supports}}),
        json!({"command":"model.step.create","params":{"name":"Modes","procedure":{"type":"frequency","modes":2},"base_state":"Base","perturbation_follower_pressure":true,"boundary_conditions":supports}})
    ]);
    recipe
}
fn run(tri: bool, divisions: usize) {
    let recipe = definition(tri, divisions);
    let mut session = Session::default();
    session.run_script(&recipe).unwrap();
    let plain = session.execute("solve.frequency", &json!({"step":"Plain","mass":"consistent"})).unwrap();
    let zero = plain["modes"][0]["frequency"].as_f64().unwrap();
    let mut previous = 0.0;
    for (pressure, follower) in [(0.0, true), (5.0, true), (10.0, true), (10.0, false)] {
        session.execute("model.load.edit", &json!({"name":"pressure","type":"pressure","set":"pressure_faces","magnitude":-pressure})).unwrap();
        session.execute("model.step.edit", &json!({"name":"Modes","perturbation_follower_pressure":follower})).unwrap();
        let preload = session.execute("solve.static", &json!({"step":"Base"})).unwrap();
        let result = session.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"})).unwrap();
        let f = result["modes"][0]["frequency"].as_f64().unwrap();
        // R=1, h=.02, E=1e7, nu=0, rho=1000, n=2; no closed-end axial load.
        let d = 1e7_f64 * 0.02_f64.powi(3) / 12.0;
        let omega2 = ((36.0 / 5.0) * d + if follower { 12.0 / 5.0 } else { 9.0 / 5.0 } * pressure) / 20.0;
        let reference = omega2.sqrt() / (2.0 * PI);
        // Fixed a priori: thin-ring approximation; refine the mesh, never loosen the gate.
        // includes the small rotary-inertia/shear effect omitted by the analytic model.
        assert!((f / reference - 1.0).abs() < 0.03, "tri={tri}, p={pressure}, follower={follower}, f={f}, ref={reference}");
        if pressure == 0.0 {
            assert!((f / zero - 1.0).abs() < 1e-9);
        }
        if follower {
            assert!(f > previous);
            previous = f;
        }
        for residual in result["convergence"]["modes"].as_array().unwrap() {
            assert!(residual["relative_residual"].as_f64().unwrap() < 1e-6, "{residual}");
        }
        let total = result["total_mass"].as_f64().unwrap();
        assert!((total / (1000.0 * 0.02 * 0.2 * PI / 2.0) - 1.0).abs() < 1e-4);
        let all: Vec<_> = session.model.nodes().keys().copied().collect();
        let elements: Vec<_> = session.model.elements().keys().copied().collect();
        let static_fields = session.execute("solve.probe", &json!({"step":"Base","nodes":all,"elements":elements})).unwrap();
        let modal_fields = session.execute("solve.probe", &json!({"step":"Modes","nodes":all,"mode":1})).unwrap();
        if let Some(dir) = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS") {
            let path = std::path::PathBuf::from(dir).join(format!("raw-tri{tri}-n{divisions}-p{pressure}-follower{follower}.json"));
            std::fs::write(path, serde_json::to_vec_pretty(&json!({"frequency":result,"static_fields":static_fields,"mode":modal_fields})).unwrap())
                .unwrap();
        }
        if !follower {
            assert!(f < previous);
        }
        let mut radial = Vec::new();
        let mut basis = Vec::new();
        for row in modal_fields["nodes"].as_array().unwrap() {
            let label = row["node"].as_u64().unwrap() as u32;
            let x = session.model.nodes()[&label].xyz;
            radial.push(x[0] * row["shape"][0].as_f64().unwrap() + x[1] * row["shape"][1].as_f64().unwrap());
            basis.push(x[0] * x[0] - x[1] * x[1]);
        }
        let fit = radial.iter().zip(&basis).map(|(w, b)| w * b).sum::<f64>() / basis.iter().map(|b| b * b).sum::<f64>();
        let shape_error =
            (radial.iter().zip(&basis).map(|(w, b)| (w - fit * b).powi(2)).sum::<f64>() / radial.iter().map(|w| w * w).sum::<f64>()).sqrt();
        assert!(shape_error < 0.01, "n=2 shape residual {shape_error}");
        let mut hoop = Vec::new();
        let mut axial = Vec::new();
        for element in static_fields["elements"].as_array().unwrap() {
            let label = element["element"].as_u64().unwrap() as u32;
            let conn = &session.model.elements()[&label].nodes;
            for (node, stress) in conn.iter().zip(element["middle"].as_array().unwrap()) {
                let x = session.model.nodes()[node].xyz;
                hoop.push(
                    x[1] * x[1] * stress[0].as_f64().unwrap() + x[0] * x[0] * stress[1].as_f64().unwrap()
                        - 2.0 * x[0] * x[1] * stress[3].as_f64().unwrap(),
                );
                axial.push(stress[2].as_f64().unwrap());
            }
        }
        let mean_hoop = hoop.iter().sum::<f64>() / hoop.len() as f64;
        let max_axial = axial.iter().map(|s| s.abs()).fold(0.0, f64::max);
        if pressure > 0.0 {
            assert!((mean_hoop / (pressure / 0.02) - 1.0).abs() < 0.01, "hoop {mean_hoop}");
            assert!(max_axial / (pressure / 0.02) < 0.01, "axial {max_axial}");
            for row in static_fields["nodes"].as_array().unwrap() {
                let x = session.model.nodes()[&(row["node"].as_u64().unwrap() as u32)].xyz;
                let ur = x[0] * row["displacement"][0].as_f64().unwrap() + x[1] * row["displacement"][1].as_f64().unwrap();
                assert!((ur / (pressure / 2e5) - 1.0).abs() < 0.01, "radial expansion {ur}");
            }
        }
        let dir = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
        let stem = format!("pressurized-tri{tri}-p{pressure}-follower{follower}");
        let archive = dir.join(format!("{stem}-{}.fieldcraft", std::process::id()));
        session.execute("file.save_as", &json!({"path":archive,"results":true})).unwrap();
        let mut opened = Session::default();
        opened.execute("file.open", &json!({"path":archive})).unwrap();
        assert_eq!(static_fields, opened.execute("solve.probe", &json!({"step":"Base","nodes":all,"elements":elements})).unwrap());
        assert_eq!(modal_fields, opened.execute("solve.probe", &json!({"step":"Modes","nodes":all,"mode":1})).unwrap());
        let row = json!({"tri":tri,"pressure":pressure,"follower":follower,"reference_hz":reference,"actual_hz":f,"relative_error":f/reference-1.0,"n2_shape_relative_l2":shape_error,"mean_hoop_stress":mean_hoop,"max_axial_stress":max_axial,"preload":preload,"frequency":result,"static_fields":static_fields,"mode1_fields":modal_fields});
        println!("pressurized shell tri={tri} p={pressure} follower={follower}: f={f}, reference={reference}");
        if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_some() {
            std::fs::write(dir.join(format!("{stem}.json")), serde_json::to_vec_pretty(&row).unwrap()).unwrap();
        } else {
            std::fs::remove_file(archive).unwrap();
        }
    }
    if let Some(dir) = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS") {
        let path = std::path::PathBuf::from(dir);
        std::fs::write(path.join(format!("pressurized-tri{tri}-recipe.json")), serde_json::to_vec_pretty(&recipe).unwrap()).unwrap();
    }
    // Named failures remain errors, never an unsymmetric eigensolve or a dropped flag.
    let error = session.execute("model.step.edit", &json!({"name":"Plain","perturbation_follower_pressure":true})).unwrap_err().to_string();
    assert!(error.contains("base_state"));
    session.execute("model.step.edit", &json!({"name":"Modes","boundary_conditions":["plane_strain"]})).unwrap();
    assert!(
        session
            .execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"}))
            .unwrap_err()
            .to_string()
            .contains("unchanged base-state supports")
    );
}
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn pressurized_mitc8_modes_and_native_reopen() {
    run(false, 4);
}
#[test]
fn curved_mitc6_modal_preload_is_refused_until_field_and_residual_checks_pass() {
    let mut session = Session::default();
    session.run_script(&definition(true, 1)).unwrap();
    let error = session.execute("solve.frequency", &json!({"step":"Modes","mass":"consistent"})).unwrap_err().to_string();
    assert!(error.contains("curved MITC6-a prestressed frequency/harmonic"));
}
#[test]
fn curved_modal_options_refuse_changed_supports_before_preload() {
    let mut s = Session::default();
    s.run_script(&definition(false, 1)).unwrap();
    s.execute("model.step.edit", &json!({"name":"Modes","boundary_conditions":["plane_strain"]})).unwrap();
    assert!(s.execute("solve.frequency", &json!({"step":"Modes"})).unwrap_err().to_string().contains("unchanged base-state supports"));
}

#[test]
fn shell_modal_preload_does_not_silently_drop_thermal_or_initial_fields() {
    let mut session = Session::default();
    session.run_script(&definition(false, 1)).unwrap();
    session.execute("model.load.create", &json!({"name":"Warm","type":"temperature","target":{"set":"all"},"value":300})).unwrap();
    session.execute("model.step.edit", &json!({"name":"Base","loads":["Warm"]})).unwrap();
    assert!(session.execute("solve.frequency", &json!({"step":"Modes"})).unwrap_err().to_string().contains("thermal preload"));
    session.execute("model.step.edit", &json!({"name":"Base","loads":["pressure"]})).unwrap();
    session
        .model
        .set_initial_field(
            "stress",
            fieldcraft_model::InitialField::Stress { values: [(1, [1.0; 6])].into(), source: "independent rejection check".into() },
        )
        .unwrap();
    assert!(session.execute("solve.frequency", &json!({"step":"Modes"})).unwrap_err().to_string().contains("initial stress/temperature fields"));
}
