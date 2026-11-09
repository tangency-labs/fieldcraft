//! Original complete models reach the native writer through deck.write. No commercial execution.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
fn session(recipe: &str) -> Session {
    let mut s = Session::default();
    s.run_script(&serde_json::from_str::<Value>(recipe).unwrap()).unwrap();
    s
}
fn write(s: &mut Session, tag: &str, units: &str) -> String {
    let path = std::env::temp_dir().join(format!("fieldcraft-dynamic-export-{tag}-{}.inp", std::process::id()));
    s.execute("deck.write", &json!({"format":"abaqus","path":path,"units":units})).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let _ = std::fs::remove_file(path);
    text
}
#[test]
fn direct_writer_keeps_fixed_step_force_and_integrator() {
    let mut s = session(include_str!("data/transient-sdof.json"));
    let text = write(&mut s, "direct", "m-kg-s-N-Pa");
    assert!(text.contains("*DYNAMIC, DIRECT, APPLICATION=TRANSIENT FIDELITY, ALPHA=0.0, BETA=0.25, GAMMA=0.5"), "{text}");
    assert!(text.contains("0.001, 1.0\n"));
    assert!(text.contains("*CLOAD, AMPLITUDE=fieldcraft-dynamic-load"));
    assert!(text.contains("1, 1, 1.0\n"));
    assert!(text.contains("*NODE OUTPUT, NSET=fieldcraft-dynamic-nodes"));
    assert!(text.contains("U, V, A"));
    let mut proc = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    proc["integration"] = json!({"type":"hht","alpha":-0.1});
    let name = s.model.steps()[0].name.clone();
    s.execute("model.step.edit", &json!({"name":name,"procedure":proc})).unwrap();
    let text = write(&mut s, "hht", "m-kg-s-N-Pa");
    assert!(text.contains("ALPHA=-0.1, BETA=0.30250000000000005, GAMMA=0.6"));
}
#[test]
fn modal_and_spectrum_writer_convert_units_and_preserve_combinations() {
    let mut s = session(include_str!("../../../docs/manual/structural/examples/modal-close-modes.json"));
    let text = write(&mut s, "modal", "mm-kg-ms-kN-GPa");
    assert!(text.contains("*MODAL DYNAMIC, CONTINUE=NO\n0.5, 1000.0\n"), "{text}");
    assert!(text.contains("1, 1, 0.02\n2, 2, 0.05\n"));
    assert!(text.contains("1, 1, 0.001\n"));
    assert!(text.contains("FREQUENCY=10"));
    let mut s = session(include_str!("../../../docs/manual/structural/examples/spectrum-close-modes.json"));
    let text = write(&mut s, "spectrum", "mm-kg-ms-kN-GPa");
    assert!(text.contains("*RESPONSE SPECTRUM, SUM=CQC, COMP=R40"), "{text}");
    assert!(text.contains("0.001, 0.001, 0.05\n0.001, 0.005, 0.05\n"));
    assert!(text.contains("fieldcraft-dynamic-spectrum, 0.0, 1.0, 0.0, 0.6"));
    assert!(!text.contains("U, V, A"));
}
#[test]
fn nonexpressible_dynamics_options_are_named_not_dropped() {
    let mut s = session(include_str!("../../../docs/manual/structural/examples/spectrum-close-modes.json"));
    let mut proc = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    proc["pseudo_acceleration"] = json!([[1, 1], [5, 2]]);
    s.execute("model.step.edit", &json!({"name":"Spectrum","procedure":proc})).unwrap();
    let path = std::env::temp_dir().join("fieldcraft-dynamics-refusal.inp");
    let e = s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap_err().to_string();
    assert!(e.contains("Spectrum") && e.contains("log-log"), "{e}");
    let mut s = session(include_str!("data/transient-sdof.json"));
    let name = s.model.steps()[0].name.clone();
    let mut proc = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    proc["initial_velocity"] = json!({"1":[1,0,0]});
    s.execute("model.step.edit", &json!({"name":name,"procedure":proc})).unwrap();
    let e = s.execute("deck.write", &json!({"format":"abaqus","path":path})).unwrap_err().to_string();
    assert!(e.contains("initial_velocity at node 1"), "{e}");
}

#[test]
fn native_amplitude_definitions_keep_time_and_fourier_conventions() {
    let mut s = session(include_str!("data/transient-sdof.json"));
    let name = s.model.steps()[0].name.clone();
    let load = s.model.steps()[0].loads.first().unwrap().clone();
    let mut proc = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
    proc["load_amplitudes"] = json!({load.clone():{"type":"periodic","origin":0,"offset":0,
        "angular_frequency":2,"harmonics":[[0,1]]}});
    s.execute("model.step.edit", &json!({"name":name,"procedure":proc})).unwrap();
    let text = write(&mut s, "periodic", "mm-kg-ms-kN-GPa");
    assert!(text.contains("DEFINITION=PERIODIC, TIME=STEP TIME\n1, 0.002, 0.0, 0.0\n0.0, 1.0\n"), "{text}");
    proc["load_amplitudes"] = json!({load:{"type":"smooth_step","points":[[0,0],[1,1]]}});
    s.execute("model.step.edit", &json!({"name":name,"procedure":proc})).unwrap();
    let text = write(&mut s, "smooth", "mm-kg-ms-kN-GPa");
    assert!(text.contains("DEFINITION=SMOOTH STEP, TIME=STEP TIME\n0.0, 0.0\n1000.0, 1.0\n"), "{text}");
}
