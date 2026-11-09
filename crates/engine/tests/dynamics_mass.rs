//! Independent one-element discrete rod: rho*A*L=3 kg, Kff=4 N/m.
//! Consistent Mff=1 kg; row-sum lumped Mff=1.5 kg. Full M*r free row=1.5 kg.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn rod() -> Session {
    let mut s = Session::default();
    s.run_script(&serde_json::from_str::<Value>(include_str!("data/transient-rod.json")).unwrap()).unwrap();
    s
}
fn probe(s: &mut Session) -> Value {
    s.execute("solve.probe", &json!({"step":"Pulse","nodes":[2]})).unwrap()
}
fn reopen(s: &mut Session, tag: &str) {
    let before = probe(s);
    let path = std::env::temp_dir().join(format!("fieldcraft-mass-{tag}-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut restored = Session::default();
    restored.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(probe(&mut restored), before);
    let _ = std::fs::remove_file(path);
}
#[test]
fn direct_and_modal_use_selected_mass_for_initial_balance_and_energy() {
    for (mass, m) in [("consistent", 1.0_f64), ("lumped", 1.5_f64)] {
        for modal in [false, true] {
            let mut s = rod();
            // Nonzero initial displacement and velocity exercise mass in both initial energy and acceleration.
            let mut procedure = json!({"type":"transient","time_step":0.001,"increments":1000,"output_every":100,
                "initial_displacement":{"2":[0.1,0,0]},"initial_velocity":{"2":[0.2,0,0]}});
            if modal {
                procedure["method"] = json!({"type":"modal","modes":1,"damping_ratios":[0]});
            }
            s.execute("model.step.edit", &json!({"name":"Pulse","procedure":procedure})).unwrap();
            let command = if modal { "solve.run" } else { "solve.transient" };
            let summary = s.execute(command, &json!({"mass":mass})).unwrap();
            let p = probe(&mut s);
            assert_eq!(p["mass"], mass);
            let omega = (4.0 / m).sqrt();
            if modal {
                assert!((summary["steps"][0]["summary"]["natural_frequencies"][0].as_f64().unwrap() - omega / std::f64::consts::TAU).abs() < 1e-12);
            }
            for (sample, frame) in p["nodes"][0]["samples"].as_array().unwrap().iter().zip(p["frames"].as_array().unwrap()) {
                let t = sample["time"].as_f64().unwrap();
                let u = 0.25 - 0.15 * (omega * t).cos() + 0.2 / omega * (omega * t).sin();
                let v = 0.15 * omega * (omega * t).sin() + 0.2 * (omega * t).cos();
                let a = (1.0 - 4.0 * u) / m;
                for (key, value, tol) in [("displacement", u, 3e-7), ("velocity", v, 5e-7), ("acceleration", a, 1.2e-6)] {
                    assert!((sample[key][0].as_f64().unwrap() - value).abs() < tol, "{mass} modal={modal} {key} t={t}: {sample}");
                }
                let un = sample["displacement"][0].as_f64().unwrap();
                let vn = sample["velocity"][0].as_f64().unwrap();
                assert!((frame["kinetic_energy"].as_f64().unwrap() - 0.5 * m * vn * vn).abs() < 1e-12);
                assert!((frame["strain_energy"].as_f64().unwrap() - 2.0 * un * un).abs() < 1e-12);
                assert!(frame["energy_balance"].as_f64().unwrap().abs() < 1e-11);
            }
            assert!((p["nodes"][0]["samples"][0]["acceleration"][0].as_f64().unwrap() - 0.6 / m).abs() < 1e-12);
            reopen(&mut s, &format!("{mass}-{modal}"));
        }
    }
}
#[test]
fn spectrum_mass_changes_frequency_and_nonflat_ordinate_and_round_trips() {
    for (mass, m) in [("consistent", 1.0_f64), ("lumped", 1.5_f64)] {
        let mut s = rod();
        s.execute(
            "model.step.edit",
            &json!({"name":"Pulse","loads":[],"procedure":{"type":"response_spectrum","modes":1,
            "damping_ratio":0.05,"combination":"cqc","direction":[1,0,0],"pseudo_acceleration":[[0.1,1.1],[1,2]]}}),
        )
        .unwrap();
        let out = s.execute("solve.spectrum", &json!({"mass":mass,"threads":1})).unwrap();
        let f = (4.0 / m).sqrt() / std::f64::consts::TAU;
        assert!((out["natural_frequencies"][0].as_f64().unwrap() - f).abs() < 1e-12);
        let p = probe(&mut s);
        assert_eq!(p["mass"], mass);
        let sample = &p["nodes"][0]["samples"][0];
        assert!((sample["displacement"][0].as_f64().unwrap() - 1.5 * (1.0 + f) / 4.0).abs() < 1e-12);
        assert!((sample["acceleration"][0].as_f64().unwrap() - 1.5 * (1.0 + f) / m).abs() < 1e-12);
        reopen(&mut s, &format!("spectrum-{mass}"));
    }
}
#[test]
fn per_call_default_compatibility_invalid_mass_and_native_export_guard() {
    let mut s = rod();
    s.execute("solve.transient", &json!({"max_memory_gib":1})).unwrap();
    let default = probe(&mut s);
    s.execute("solve.transient", &json!({"mass":"consistent"})).unwrap();
    assert_eq!(probe(&mut s), default);
    assert!(s.execute("solve.transient", &json!({"mass":"row_sum_typo"})).is_err());
    let path = std::env::temp_dir().join(format!("fieldcraft-mass-{}.inp", std::process::id()));
    let err = s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap_err().to_string();
    assert!(err.contains("rod2") && err.contains("mass=Consistent"), "{err}");
    s.execute("solve.transient", &json!({"mass":"lumped"})).unwrap();
    assert!(s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap_err().to_string().contains("stored mass"));
    s.execute("deck.write", &json!({"format":"abaqus","path":path,"mass":"lumped"})).unwrap();
    assert!(std::fs::read_to_string(&path).unwrap().contains("TYPE=T3D2"));
    s.execute("solve.transient", &json!({})).unwrap();
    assert_eq!(probe(&mut s), default);
    let _ = std::fs::remove_file(path);
}
#[test]
fn material_and_step_mass_damping_match_discrete_damped_closed_form() {
    for (mass, m) in [("consistent", 1.0_f64), ("lumped", 1.5_f64)] {
        for material in [false, true] {
            let mut s = rod();
            if material {
                s.execute(
                    "model.material.edit",
                    &json!({"name":"Elastic","elastic":{"type":"isotropic","young":4,"poisson":0.25},"density":3,"damping":{"alpha":0.2}}),
                )
                .unwrap();
            } else {
                s.execute("model.step.edit",&json!({"name":"Pulse","procedure":{"type":"transient","time_step":0.001,"increments":1000,"output_every":100,"damping":{"alpha":0.2}}})).unwrap();
            }
            s.execute("solve.transient", &json!({"mass":mass})).unwrap();
            let p = probe(&mut s);
            let w = (4.0 / m - 0.01).sqrt();
            for sample in p["nodes"][0]["samples"].as_array().unwrap() {
                let t = sample["time"].as_f64().unwrap();
                let u = 0.25 * (1.0 - (-0.1 * t).exp() * ((w * t).cos() + 0.1 / w * (w * t).sin()));
                assert!((sample["displacement"][0].as_f64().unwrap() - u).abs() < 3e-7, "{mass} material={material}");
            }
            assert!(p["frames"].as_array().unwrap().last().unwrap()["damping_dissipation"].as_f64().unwrap() > 0.0);
        }
    }
}
