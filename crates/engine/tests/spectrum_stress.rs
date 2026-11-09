//! Original orthogonal two-rod assembly. Mff=I, full base participation=1.5 I.
//! At TIP, each modal tensor is averaged over the two adjacent rods before CQC.
use fieldcraft_engine::Session;
use serde_json::{Value, json};
#[test]
fn signed_modal_stress_recovery_precedes_modal_and_direction_combinations() {
    let recipe: Value = serde_json::from_str(include_str!("../../../docs/manual/structural/examples/spectrum-stress-rods.json")).unwrap();
    let ratio = (194.0449_f64 / 192.3769).sqrt();
    let z = 0.05_f64;
    let rho = 8.0 * z * z * (1.0 + ratio) * ratio.powf(1.5) / ((1.0 - ratio * ratio).powi(2) + 4.0 * z * z * ratio * (1.0 + ratio).powi(2));
    let q = 0.5_f64.sqrt();
    for combination in ["srss", "cqc"] {
        for direction in ["single", "srss", "percent1004040"] {
            let mut s = Session::default();
            s.run_script(&recipe).unwrap();
            let mut proc = serde_json::to_value(&s.model.steps()[0].procedure).unwrap();
            proc["combination"] = json!(combination);
            if direction != "single" {
                proc["directions"] = json!([{"direction":[1,0,0],"scale":1},{"direction":[0,1,0],"scale":0.6}]);
                proc["direction_combination"] = json!(direction);
            }
            s.execute("model.step.edit", &json!({"name":"Spectrum","procedure":proc})).unwrap();
            s.execute("solve.spectrum", &json!({"threads":1})).unwrap();
            let p = s.execute("solve.probe", &json!({"step":"Spectrum","nodes":[1]})).unwrap();
            let stress = p["nodes"][0]["samples"][0]["stress"].as_array().unwrap();
            let mut maximum = 0.0_f64;
            for (component, actual) in stress.iter().enumerate() {
                let modal = |axis: usize| {
                    let mut v = [0.0; 2];
                    for (i, n) in [[q, q, 0.0], [-q, q, 0.0]].iter().enumerate() {
                        let tensor = [n[0] * n[0], n[1] * n[1], 0.0, n[0] * n[1], 0.0, 0.0];
                        v[i] = 0.75 * n[axis] * tensor[component];
                    }
                    (v[0] * v[0] + v[1] * v[1] + if combination == "cqc" { 2.0 * rho * v[0] * v[1] } else { 0.0 }).max(0.0).sqrt()
                };
                let x = modal(0);
                let y = 0.6 * modal(1);
                let expected = match direction {
                    "single" => x,
                    "srss" => (x * x + y * y).sqrt(),
                    _ => (x + 0.4 * y).max(y + 0.4 * x),
                };
                maximum = maximum.max((actual.as_f64().unwrap() - expected).abs());
            }
            println!("SPECTRUM_STRESS {combination}/{direction} max_error_Pa={maximum:.16e}, tolerance_Pa=1e-10, actual={stress:?}");
            assert!(maximum < 1e-10);
            assert!(p["nodes"][0]["samples"][0]["time"].is_null());
        }
    }
}
