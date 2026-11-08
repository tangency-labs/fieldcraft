//! Opt-in project convergence replay; no benchmark pass is inferred by this driver.
use fieldcraft_engine::{Session, StepResult};
use fieldcraft_solver::StepResult as R;
use serde_json::{Value, json};
use std::{fs, time::Instant};
#[test]
#[ignore = "project Hex20 5R convergence; explicit QUADRATIC_5R_DIR"]
fn quadratic_5r_project() {
    let Ok(dir) = std::env::var("QUADRATIC_5R_DIR") else {
        eprintln!("project replay skipped: QUADRATIC_5R_DIR was not supplied (not a CI oracle)");
        return;
    };
    for case in ["q1", "q2", "q3"] {
        let start = Instant::now();
        let save = |suffix: &str, v: &Value| {
            fs::write(format!("{dir}/{case}-{suffix}.json"), serde_json::to_vec_pretty(v).unwrap()).unwrap();
        };
        let mut s = Session::default();
        let input: Value = serde_json::from_slice(&fs::read(format!("{dir}/{case}.json")).unwrap()).unwrap();
        s.run_script(&input).unwrap();
        save("check", &s.execute("model.check", &json!({})).unwrap());
        save("static", &s.execute("solve.static", &json!({"step":"S","threads":16})).unwrap());
        save("modes", &s.execute("solve.frequency", &json!({"step":"M","threads":16})).unwrap());
        let StepResult::Implicit(R::Frequency(m)) = s.result("M").unwrap().as_ref() else { panic!("frequency missing") };
        let mut frequencies: Vec<f64> = (0..=200).map(|k| k as f64 * 5.).collect();
        frequencies.push(42.65);
        for mode in &m.modes {
            for k in -60..=60 {
                let f = mode.frequency * (1. + k as f64 * 0.001);
                if f > 0. && f < 1000. {
                    frequencies.push(f);
                }
            }
        }
        frequencies.sort_by(f64::total_cmp);
        frequencies.dedup_by(|a, b| (*a - *b).abs() < 1e-7);
        if let Ok(path) = std::env::var("COMMON_5R_GRID") {
            frequencies = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
        }
        s.execute("model.step.edit", &json!({"name":"H5R","procedure":{"type":"harmonic","frequencies":frequencies,"method":{"type":"modal","modes":16,"damping_ratios":[0.02]}}})).unwrap();
        save("frequency-grid", &json!(frequencies));
        save("harmonic", &s.execute("solve.harmonic", &json!({"step":"H5R","threads":16})).unwrap());
        let meta: Value = serde_json::from_slice(&fs::read(format!("{dir}/{case}-input.json")).unwrap()).unwrap();
        save("probes", &s.execute("solve.probe", &json!({"step":"H5R","nodes":meta["probe_nodes"]})).unwrap());
        save("execution", &json!({"elapsed_s":start.elapsed().as_secs_f64(),"completed":true}));
        println!("{case} completed in {:.3}s", start.elapsed().as_secs_f64());
    }
}
