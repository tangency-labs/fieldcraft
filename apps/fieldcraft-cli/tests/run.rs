//! `fieldcraft-cli run` on a saved model: a project with a saved render writes its clip with no
//! script of its own beyond opening it and naming the render (#646).

use std::process::Command;

use fieldcraft_engine::Session;
use serde_json::json;

#[test]
fn a_saved_render_writes_its_clip_from_the_cli() {
    let d = std::env::temp_dir().join(format!("fieldcraft-cli-run-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    // A small sheared gap, solved, shown on a slice and kept with a view and a render.
    let periodic = json!({"type": "periodic"});
    let mut s = Session::default();
    s.run_script(&json!([
        {"command": "model.material.create", "params": {"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}}},
        {"command": "model.fluid_domain.create", "params": {"name": "Gap", "material": "Oil", "min": [0, 0, 0], "max": [0.125, 1, 0.125], "spacing": 0.125,
            "faces": {"x_min": periodic, "x_max": periodic, "z_min": periodic, "z_max": periodic}}},
        {"command": "model.fluid_domain.boundary", "params": {"domain": "Gap", "at": "y_max", "condition": {"type": "moving_wall", "velocity": [0.1, 0, 0]}}},
        {"command": "model.step.create", "params": {"name": "Shear", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "frames": 4}}},
        {"command": "solve.fluid", "params": {}},
        {"command": "results.slice", "params": {"axis": "z"}},
        {"command": "view.save", "params": {"name": "Gap", "camera": {"yaw": -30, "pitch": 20, "fit": true}}},
        {"command": "results.render.save", "params": {"name": "Shear clip", "view": "Gap", "size": [80, 48], "caption": "Couette flow"}},
    ]))
    .unwrap();
    let model = d.join("gap.fieldcraft");
    s.execute("file.save_as", &json!({"path": model})).unwrap();

    let script = d.join("render.json");
    let clip = d.join("clip.mp4");
    std::fs::write(
        &script,
        json!([
            {"command": "file.open", "params": {"path": model}},
            {"command": "results.export.animation", "params": {"render": "Shear clip", "path": clip}},
        ])
        .to_string(),
    )
    .unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_fieldcraft-cli")).arg("run").arg(&script).env("FIELDCRAFT_GPU", "cpu").output().unwrap();
    let text = format!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    assert!(out.status.success(), "{text}");
    for n in 0..5 {
        assert!(d.join(format!("clip_{n:05}.png")).is_file(), "picture {n}: {text}");
    }
    assert!(!d.join("clip_00005.png").exists());
    let _ = std::fs::remove_dir_all(&d);
}
