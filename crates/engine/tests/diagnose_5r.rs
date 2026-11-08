use fieldcraft_engine::{Session, StepResult};
use fieldcraft_solver::StepResult as R;
use serde_json::{Value, json};
#[test]
#[ignore = "project diagnostic; explicit input/output paths"]
fn diagnose_5r() {
    let Ok(dir) = std::env::var("FIVE_R_DIR") else {
        eprintln!("project replay skipped: FIVE_R_DIR was not supplied (not a CI oracle)");
        return;
    };
    for case in ["corner", "center", "line"] {
        let input = std::fs::read_to_string(format!("{dir}/{case}.json")).unwrap();
        let mut s = Session::default();
        s.run_script(&serde_json::from_str::<Value>(&input).unwrap()).unwrap();
        let check = s.execute("model.check", &json!({})).unwrap();
        let stat = s.execute("solve.static", &json!({"step":"S","threads":1})).unwrap();
        s.execute("solve.frequency", &json!({"step":"M","threads":1})).unwrap();
        let StepResult::Implicit(R::Frequency(m)) = s.result("M").unwrap().as_ref() else { panic!() };
        let modes = serde_json::to_value(m).unwrap();
        let mut freqs: Vec<f64> = (0..=500).map(|i| i as f64 * 2.).collect();
        freqs.push(42.65);
        for mode in &m.modes {
            for i in -80..=80 {
                let f = mode.frequency * (1. + i as f64 * 0.001);
                if f > 0. && f < 1000. {
                    freqs.push(f);
                }
            }
        }
        freqs.sort_by(f64::total_cmp);
        freqs.dedup();
        s.execute(
            "model.step.edit",
            &json!({"name":"H5R","procedure":{"type":"harmonic","frequencies":freqs,"method":{"type":"modal","modes":16,"damping_ratios":[0.02]}}}),
        )
        .unwrap();
        let hs = s.execute("solve.harmonic", &json!({"step":"H5R","threads":1})).unwrap();
        let StepResult::Implicit(R::Harmonic(h)) = s.result("H5R").unwrap().as_ref() else { panic!() };
        std::fs::write(format!("{dir}/{case}-harmonic.json"), serde_json::to_vec(h).unwrap()).unwrap();
        let random = s.execute("solve.random", &json!({"step":"H5R","input_psd":[[0.,1e12],[1000.,1e12]]})).unwrap();
        std::fs::write(format!("{dir}/{case}-random.json"), serde_json::to_vec(&random).unwrap()).unwrap();
        std::fs::write(
            format!("{dir}/{case}-summary.json"),
            serde_json::to_vec_pretty(&json!({"check":check,"static":stat,"modes":modes,"harmonic":hs})).unwrap(),
        )
        .unwrap();
        println!("{case}: complete");
    }
}
