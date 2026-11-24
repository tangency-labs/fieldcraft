use serde_json::json;

use crate::Session;

/// The sample plate in a box of water flowing along x, written as an OpenFOAM case: a patch
/// per named face of the plate, lengths in metres, nothing of the model changed.
#[test]
fn writes_a_fluid_step_as_an_openfoam_case() {
    let (_guard, design) = crate::tests::sample_design_file(&fieldcraft_geometry::sample::Plate::default());
    let dir = std::path::Path::new(&design).parent().unwrap().to_path_buf();
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": design})).unwrap();
    let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
    let [lo, hi] = [0, 1].map(|k| bodies[0]["bbox"][k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>());
    let min: Vec<f64> = lo.iter().map(|v| (v / 5.0).floor() * 5.0 - 20.0).collect();
    let max: Vec<f64> = hi.iter().map(|v| (v / 5.0).ceil() * 5.0 + 20.0).collect();
    let script = json!([
        {"command": "model.material.create", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-9}}},
        {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"],
            "faces": {"x_min": {"type": "velocity", "velocity": [1000, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0},
                      "z_min": {"type": "slip"}}}},
        {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10}}},
    ]);
    s.run_script(&script).unwrap();
    let journal = s.journal.len();
    let case = dir.join("case");
    // The first lattice cell's centre, in the inlet's corner, away from the plate.
    let location: Vec<f64> = min.iter().map(|v| v + 2.5).collect();
    let out = s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "processes": 4, "location": location})).unwrap();
    // The reply, in the input units, is the point written: the centre of the finest cell holding
    // the given one, which here is that point itself (the first lattice cell's centre).
    let written: Vec<f64> = out["location"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    for (w, l) in written.iter().zip(&location) {
        assert!((w - l).abs() <= 1e-9 * l.abs().max(1.0), "{out}");
    }
    assert_eq!(s.journal.len(), journal, "writing a case does not change the model");
    assert_eq!(out["solver"], "pimpleFoam");
    let read = |rel: &str| std::fs::read_to_string(case.join(rel)).unwrap();
    // Lengths in metres: the 5 mm spacing is the background cell.
    let extent: Vec<f64> = (0..3).map(|i| (max[i] - min[i]) / 5.0).collect();
    assert_eq!(out["background_cells"], json!(extent.iter().map(|n| n.round() as u64).collect::<Vec<_>>()), "{out}");
    assert!(read("system/blockMeshDict").contains(&format!("({:e} ", min[0] / 1000.0)), "{}", read("system/blockMeshDict"));
    // A solid and a patch per named face of the plate; m/s and m²/s.
    let faces = bodies[0]["faces"].as_array().unwrap();
    let stl = read("constant/triSurface/Body1.stl");
    assert_eq!(stl.matches("endsolid").count(), faces.len(), "{out}");
    assert_eq!(out["patches"].as_object().unwrap().len(), 6 + faces.len());
    assert!(read("0/U").contains("x_min { type fixedValue; value uniform (1e0 0 0); }"), "{}", read("0/U"));
    assert!(read("constant/transportProperties").contains("nu 1e-6;"));
    assert!(read("system/decomposeParDict").contains("numberOfSubdomains 4;"));
    assert!(read("Allrun").contains("snappyHexMesh -overwrite"));
    // A second export into the same directory is refused and leaves the first case as it was.
    let before = read("Allrun");
    let e = s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap(), "location": location})).unwrap_err().to_string();
    assert!(e.contains("the case directory is not empty"), "{e}");
    assert_eq!(read("Allrun"), before);
    // A directory holding only a folder (an old case's time or system directory) is not empty either.
    let folders = dir.join("folders");
    std::fs::create_dir_all(folders.join("0.5")).unwrap();
    let e = s.execute("deck.openfoam", &json!({"path": folders.to_str().unwrap(), "location": location})).unwrap_err().to_string();
    assert!(e.contains("the case directory is not empty"), "{e}");
    // What OpenFOAM cannot be given yet is named, and nothing is written.
    s.execute(
        "model.fluid_domain.edit",
        &json!({"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"], "acceleration": [0, 0, -9810]}),
    )
    .unwrap();
    let other = dir.join("other");
    let e = s.execute("deck.openfoam", &json!({"path": other.to_str().unwrap(), "location": location})).unwrap_err().to_string();
    assert!(e.contains("acceleration"), "{e}");
    assert!(!other.exists());
}

#[test]
fn a_path_is_one_posix_shell_word() {
    assert_eq!(super::shell_quote("/tmp/my case's dir"), r"'/tmp/my case'\''s dir'");
    assert_eq!(super::shell_quote("plain"), "'plain'");
}
