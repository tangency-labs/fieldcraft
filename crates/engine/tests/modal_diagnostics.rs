use fieldcraft_engine::Session;
use serde_json::{Value, json};

#[test]
fn frequency_reports_actual_equation_residuals_and_iteration_count() {
    let mut s = Session::default();
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/cantilever-frequencies/recipe.json")).unwrap();
    s.run_script(&recipe).unwrap();
    let mut reference = None;
    for threads in [1, 2, 7] {
        let result = s.execute("solve.frequency", &json!({"step":"Modes", "threads":threads})).unwrap();
        let c = &result["convergence"];
        let iterations = c["iterations"].as_u64().unwrap();
        let vectors = result["lanczos_vectors"].as_u64().unwrap();
        assert!(iterations > 0 && iterations < vectors);
        assert_eq!(iterations, vectors.div_ceil(6));
        for mode in result["modes"].as_array().unwrap() {
            let r = &mode["equation_residual"];
            let ratio = r["relative_residual"].as_f64().unwrap();
            let den = r["stiffness_action_norm"].as_f64().unwrap();
            assert!(den > 0.0 && ratio < 1e-6, "{r}");
            assert_eq!(ratio, r["residual_norm"].as_f64().unwrap() / den);
        }
        if let Some(r) = &reference {
            assert_eq!(&result, r);
        } else {
            reference = Some(result);
        }
    }
    // A nonzero spectral shift must not leak into the reported physical K residual.
    let shifted = s.execute("solve.frequency", &json!({"step":"Modes", "shift":-0.1, "threads":1})).unwrap();
    for mode in shifted["modes"].as_array().unwrap() {
        assert!(mode["equation_residual"]["relative_residual"].as_f64().unwrap() < 1e-6, "{mode}");
    }
    let error = s.execute("solve.frequency", &json!({"step":"Modes", "max_iterations":6})).unwrap_err().to_string();
    assert!(error.contains("basis_limit_reached=true") && error.contains("1 Lanczos block iterations"), "{error}");
}
