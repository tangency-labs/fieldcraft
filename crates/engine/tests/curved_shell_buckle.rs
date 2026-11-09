//! Classical perfect-cylinder bifurcation, independently generated quarter meshes.
//! Axial: NASA SP-8007 classical elastic result (knockdown factor = 1).
//! Dead radial traction: inextensional n=2 ring, p=4D/R^3; NOT follower pressure
//! (whose ring limit is 3D/R^3). The shell solver's load vector is fixed.
use fieldcraft_engine::Session;
use serde_json::json;
#[path = "common/cylinder.rs"]
mod cylinder;
use cylinder::model;

fn run_case(tri: bool, pressure: bool, follower: bool) {
    let (recipe, reference) = model(tri, pressure, follower, 4);
    let mut s = Session::default();
    s.run_script(&recipe).unwrap();
    let result = s.execute("solve.buckle", &json!({"step":"Buckle"})).unwrap();
    let factor = result["modes"][0]["load_factor"].as_f64().unwrap();
    println!(
        "curved cylinder tri={tri} pressure={pressure} follower={follower} reference={reference} actual={factor} relative_error={} result={result}",
        factor / reference - 1.0
    );
    let dir = std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").map(std::path::PathBuf::from).unwrap_or_else(std::env::temp_dir);
    let stem = format!("cylinder-tri{tri}-pressure{pressure}-follower{follower}");
    let archive = dir.join(format!("{stem}-{}.fieldcraft", std::process::id()));
    let p = json!({"step":"Buckle","nodes":[1,2,3],"mode":1});
    let before = s.execute("solve.probe", &p).unwrap();
    s.execute("file.save_as", &json!({"path":archive,"results":true})).unwrap();
    let mut open = Session::default();
    open.execute("file.open", &json!({"path":archive})).unwrap();
    assert_eq!(before, open.execute("solve.probe", &p).unwrap());
    if std::env::var_os("FIELDCRAFT_TEST_ARTIFACTS").is_some() {
        std::fs::write(dir.join(format!("{stem}-recipe.json")), serde_json::to_vec_pretty(&recipe).unwrap()).unwrap();
        std::fs::write(dir.join(format!("{stem}-result.json")), serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    } else {
        std::fs::remove_file(archive).unwrap();
    }
    // Fixed before runs: coarse quadratic mesh, thin-shell asymptotic reference.
    assert!((factor / reference - 1.0).abs() < 0.05, "{result}");
}
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn cylinders_through_buckling_commands_and_reopen() {
    for tri in [false, true] {
        for pressure in [false, true] {
            run_case(tri, pressure, false);
        }
    }
}
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn follower_cylinder_pressure_and_reopen() {
    for tri in [false, true] {
        run_case(tri, true, true);
    }
}
