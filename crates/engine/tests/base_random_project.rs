//! Project-only bounded BASE replay; opt-in paths and no production commands.
use fieldcraft_engine::{Session, StepResult};
use fieldcraft_solver::StepResult as R;
use serde_json::{Value, json};
use std::{fs, time::Instant};
#[test]
#[ignore = "project BASE replay; explicit BASE_RANDOM_DIR"]
fn base_random_project() {
    let Ok(d) = std::env::var("BASE_RANDOM_DIR") else {
        eprintln!("project replay skipped: BASE_RANDOM_DIR was not supplied (not a CI oracle)");
        return;
    };
    let save = |name: &str, value: &Value| {
        fs::write(format!("{d}/{name}.json"), serde_json::to_vec_pretty(value).unwrap()).unwrap();
    };
    let start = Instant::now();
    let mut s = Session::default();
    let script: Value = serde_json::from_slice(&fs::read(format!("{d}/recipe.json")).unwrap()).unwrap();
    for (i, row) in script.as_array().unwrap().iter().enumerate() {
        let command = row["command"].as_str().unwrap();
        match s.execute(command, &row["params"]) {
            Ok(v) => save(&format!("setup-{i}"), &v),
            Err(e) => {
                save("failure", &json!({"command":command,"error":e.to_string(),"elapsed_s":start.elapsed().as_secs_f64()}));
                panic!("{e}");
            }
        }
    }
    save("check", &s.execute("model.check", &json!({})).unwrap());
    println!("BASE prepared; starting modal harmonic");
    let response = match s.execute("solve.harmonic", &json!({"step":"BASE_RANDOM","threads":16,"max_factor_entries":2000000000u64})) {
        Ok(v) => v,
        Err(e) => {
            save("failure", &json!({"command":"solve.harmonic","error":e.to_string(),"elapsed_s":start.elapsed().as_secs_f64()}));
            panic!("{e}");
        }
    };
    save("harmonic-summary", &response);
    let probe = s.execute("solve.probe", &json!({"step":"BASE_RANDOM","at":[[-0.163,0.002,0.0219],[-0.175,0.0518,0.0189]]})).unwrap();
    save("probes", &probe);
    let StepResult::Implicit(R::Harmonic(h)) = s.result("BASE_RANDOM").unwrap().as_ref() else { panic!("not harmonic") };
    save(
        "frame-audit",
        &json!({"natural_frequencies":h.natural_frequencies,"warnings":h.warnings,"frames":h.frames.iter().map(|f|json!({"f":f.frequency,"relative_residual":f.relative_residual,"projected_residual":f.projected_residual})).collect::<Vec<_>>(),"elapsed_s":start.elapsed().as_secs_f64()}),
    );
}
