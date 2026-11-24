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

/// A "run" of the plate's case made by hand in the documented ASCII format: OpenFOAM cells at the
/// centres of a 5 mm grid over the box (every fluid lattice cell near one), and the flow at two
/// times. It reads back as the step's fluid results, reconstructed and decomposed alike; a case
/// of another domain or one with a gap is refused, never read as zeros.
#[test]
fn reads_a_cases_results_back_onto_the_lattice() {
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
            "faces": {"x_min": {"type": "velocity", "velocity": [1000, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}}},
        {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10}}},
    ]);
    s.run_script(&script).unwrap();
    let case = dir.join("case");
    s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap()})).unwrap();
    let journal = s.journal.len();
    assert!(s.execute("results.openfoam", &json!({"path": case.to_str().unwrap()})).is_err(), "not run yet: nothing to read");

    // OpenFOAM cells at the 5 mm grid's centres (metres); x below the middle flows at u, above at 2u.
    let (lo_m, hi_m): (Vec<f64>, Vec<f64>) = (min.iter().map(|v| v / 1000.0).collect(), max.iter().map(|v| v / 1000.0).collect());
    let n: Vec<usize> = (0..3).map(|a| ((max[a] - min[a]) / 5.0).round() as usize).collect();
    let mut centres = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                centres.push([i, j, k].iter().enumerate().map(|(a, &c)| lo_m[a] + (c as f64 + 0.5) * 0.005).collect::<Vec<f64>>());
            }
        }
    }
    let middle = (lo_m[0] + hi_m[0]) / 2.0;
    let vector_file = |items: &[String]| {
        format!(
            "FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class volVectorField;\n    object U;\n}}\ninternalField nonuniform List<vector> {}\n(\n{}\n)\n;\nboundaryField\n{{\n}}\n",
            items.len(),
            items.join("\n")
        )
    };
    let point = |c: &Vec<f64>| format!("({} {} {})", c[0], c[1], c[2]);
    std::fs::write(case.join("0/C"), vector_file(&centres.iter().map(point).collect::<Vec<_>>())).unwrap();
    for (t, u) in [("0.25", 0.5), ("0.5", 1.0)] {
        std::fs::create_dir_all(case.join(t)).unwrap();
        let v: Vec<String> = centres.iter().map(|c| format!("({} 0 0)", if c[0] > middle { 2.0 * u } else { u })).collect();
        std::fs::write(case.join(t).join("U"), vector_file(&v)).unwrap();
        std::fs::write(case.join(t).join("p"), "internalField uniform 2;\n").unwrap();
    }
    let r = s.execute("results.openfoam", &json!({"path": case.to_str().unwrap()})).unwrap();
    assert_eq!(
        (r["frames"].clone(), r["source"].clone(), r["openfoam_cells"].clone()),
        (json!(2), json!("reconstructed"), json!(centres.len())),
        "{r}"
    );
    assert_eq!(s.journal.len(), journal, "reading results does not change the model");
    assert!(r["fluid_cells"].as_u64().unwrap() < centres.len() as u64, "the plate's cells are solid: {r}");
    // Probed through the ordinary fluid commands: the nearest cell's flow, and p × ρ as the gauge pressure.
    let probe = |s: &mut Session, x: f64, frame: u64| {
        let pt = [x, min[1] + 2.5, min[2] + 2.5];
        s.execute("solve.probe", &json!({"step": "Flow", "points": [pt], "frame": frame})).unwrap()
    };
    let high = probe(&mut s, max[0] - 2.5, 1);
    assert_eq!(high["points"][0]["velocity"][0].as_f64(), Some(2000.0), "2 m/s at the outlet side at 0.5 s, in mm/s: {high}");
    let low = probe(&mut s, min[0] + 2.5, 0);
    assert_eq!(low["points"][0]["velocity"][0].as_f64(), Some(500.0), "0.5 m/s at the inlet side at 0.25 s, in mm/s: {low}");

    // Fewer frames: a range and every n-th, always with the last.
    let read_with = |s: &mut Session, path: &std::path::Path, extra: serde_json::Value| {
        let mut params = json!({"path": path.to_str().unwrap()});
        if let (Some(o), Some(e)) = (params.as_object_mut(), extra.as_object()) {
            o.extend(e.clone());
        }
        s.execute("results.openfoam", &params)
    };
    assert_eq!(read_with(&mut s, &case, json!({"to": 0.3})).unwrap()["frames"], 1);
    let e = read_with(&mut s, &case, json!({"from": 0.3, "to": 0.4})).unwrap_err().to_string();
    assert!(e.contains("no written time") && e.contains("0.25 s to 0.5 s"), "times in the display units: {e}");
    // A bound matches a time directory's name within half a unit of its last digit: the case is
    // written with 9 significant digits (10 frames), so 0.5·10⁻⁸ relative, and no further.
    for (from, to) in [(0.25 * (1.0 + 4e-9), 0.3), (0.0, 0.25 * (1.0 - 4e-9))] {
        assert_eq!(read_with(&mut s, &case, json!({"from": from, "to": to})).unwrap()["frames"], 1, "{from} {to}");
    }
    for (from, to) in [(0.25 * (1.0 + 6e-9), 0.3), (0.0, 0.25 * (1.0 - 6e-9))] {
        assert!(read_with(&mut s, &case, json!({"from": from, "to": to})).is_err(), "{from} {to} is past the last digit");
    }
    assert_eq!(read_with(&mut s, &case, json!({"from": 0.26, "to": 0.5})).unwrap()["frames"], 1);
    assert!(read_with(&mut s, &case, json!({"every": 0})).is_err());

    // Decomposed (no time directory beside system): each processor's own centres and fields,
    // here in a permuted order, give the same flow as the reconstructed case.
    read_with(&mut s, &case, json!({})).unwrap();
    let expected = probe(&mut s, max[0] - 2.5, 1)["points"].clone();
    let split = dir.join("split");
    std::fs::create_dir_all(&split).unwrap();
    std::fs::copy(case.join("fieldcraft.json"), split.join("fieldcraft.json")).unwrap();
    let ranks: [Vec<usize>; 2] = [(0..centres.len()).filter(|i| i % 2 == 1).collect(), (0..centres.len()).filter(|i| i % 2 == 0).collect()];
    for (rank, ids) in ranks.iter().enumerate() {
        let part = split.join(format!("processor{rank}"));
        std::fs::create_dir_all(part.join("0")).unwrap();
        std::fs::write(part.join("0/C"), vector_file(&ids.iter().map(|&i| point(&centres[i])).collect::<Vec<_>>())).unwrap();
        for (t, u) in [("0.25", 0.5), ("0.5", 1.0)] {
            std::fs::create_dir_all(part.join(t)).unwrap();
            let v: Vec<String> = ids.iter().map(|&i| format!("({} 0 0)", if centres[i][0] > middle { 2.0 * u } else { u })).collect();
            std::fs::write(part.join(t).join("U"), vector_file(&v)).unwrap();
            std::fs::write(part.join(t).join("p"), "internalField uniform 2;\n").unwrap();
        }
    }
    let r = read_with(&mut s, &split, json!({})).unwrap();
    assert_eq!(r["source"], "decomposed (2 processors)", "{r}");
    assert_eq!(probe(&mut s, max[0] - 2.5, 1)["points"], expected, "the processors give the reconstructed flow");
    // A rank's field of another mesh than its centres is refused.
    std::fs::write(split.join("processor1/0.5/U"), vector_file(&["(1 0 0)".to_string()])).unwrap();
    assert!(read_with(&mut s, &split, json!({})).is_err());

    // A case whose cells lie elsewhere (another model's) is refused before any search.
    let moved: Vec<String> = centres.iter().map(|c| format!("({} {} {})", c[0] + 10.0, c[1], c[2])).collect();
    std::fs::write(case.join("0/C"), vector_file(&moved)).unwrap();
    let e = read_with(&mut s, &case, json!({})).unwrap_err().to_string();
    assert!(e.contains("no OpenFOAM cell centre lies in the step's domain"), "{e}");
    // A case that leaves part of the domain without cells is an error, never zero flow there.
    let half: Vec<String> = centres.iter().filter(|c| c[0] < middle).map(point).collect();
    std::fs::write(case.join("0/C"), vector_file(&half)).unwrap();
    for t in ["0.25", "0.5"] {
        std::fs::write(case.join(t).join("U"), vector_file(&vec!["(1 0 0)".to_string(); half.len()])).unwrap();
    }
    let e = read_with(&mut s, &case, json!({})).unwrap_err().to_string();
    assert!(e.contains("does not cover the step's domain") && e.contains(" mm"), "the position and reach in display units: {e}");
    // A marker missing a field is refused, not filled in with a default.
    std::fs::write(case.join("fieldcraft.json"), r#"{"step": "Flow", "solver": "pimpleFoam"}"#).unwrap();
    assert!(read_with(&mut s, &case, json!({})).unwrap_err().to_string().contains("fieldcraft.json"));
}

/// The plate's case written by `deck.openfoam`, with its session, domain box and a time directory
/// per name in `times` (empty folders: what these cases refuse comes before any file is read).
fn written_case(times: &[&str]) -> (crate::tests::TempDir, Session, std::path::PathBuf, Vec<f64>, Vec<f64>) {
    let (guard, design) = crate::tests::sample_design_file(&fieldcraft_geometry::sample::Plate::default());
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
            "faces": {"x_min": {"type": "velocity", "velocity": [1000, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}}},
        {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "end_time": 0.5, "frames": 10}}},
    ]);
    s.run_script(&script).unwrap();
    let case = dir.join("case");
    s.execute("deck.openfoam", &json!({"path": case.to_str().unwrap()})).unwrap();
    for t in times {
        std::fs::create_dir_all(case.join(t)).unwrap();
    }
    (guard, s, case, min, max)
}

fn refusal(s: &mut Session, case: &std::path::Path, extra: serde_json::Value) -> String {
    let mut params = json!({"path": case.to_str().unwrap()});
    if let (Some(o), Some(e)) = (params.as_object_mut(), extra.as_object()) {
        o.extend(e.clone());
    }
    s.execute("results.openfoam", &params).unwrap_err().to_string()
}

/// `max_cells` is raised only up to what the session's memory bound holds, as `solve.fluid`'s is:
/// a huge one cannot overflow the lattice's index.
#[test]
fn results_max_cells_is_capped_by_the_memory_bound() {
    let (_guard, mut s, case, ..) = written_case(&["0.25"]);
    s.limits.memory = 4_000_000_000;
    let e = refusal(&mut s, &case, json!({"max_cells": 1e15}));
    assert!(e.contains("max_cells") && e.contains("this machine allows") && e.contains("4000 MB"), "{e}");
    s.limits.memory = 0;
    let e = refusal(&mut s, &case, json!({"max_cells": 1e15}));
    assert!(e.contains("max_cells") && e.contains("cannot be raised"), "{e}");
}

/// A spacing the box is not a whole number of is refused, as the model's own lattice does.
#[test]
fn results_spacing_must_span_the_box_in_whole_cells() {
    let (_guard, mut s, case, ..) = written_case(&["0.25"]);
    let e = refusal(&mut s, &case, json!({"spacing": 30}));
    assert!(e.contains("whole number"), "{e}");
}

/// The lattice, the frames and the case's centres are checked against the session's memory bound
/// before they are allocated.
#[test]
fn results_peak_memory_is_checked_against_the_bound() {
    let (_guard, mut s, case, ..) = written_case(&["0.25"]);
    s.limits.memory = 10_000;
    let e = refusal(&mut s, &case, json!({}));
    assert!(e.contains("memory") && e.contains("peak"), "{e}");
}

/// Many frames of a fine lattice pass the cell cap and still fail the frame memory.
#[test]
fn results_frames_are_bounded_by_the_frame_memory() {
    let names: Vec<String> = (1..=40).map(|i| i.to_string()).collect();
    let (_guard, mut s, case, min, max) = written_case(&names.iter().map(String::as_str).collect::<Vec<_>>());
    s.limits.memory = 0;
    // The finest spacing (mm) with more than 3.2 million cells, 16 bytes a cell: 40 frames pass 2 GB.
    let base: f64 = (0..3).map(|a| (max[a] - min[a]) / 5.0).product();
    let k = (1..).find(|k| base * (*k as f64).powi(3) > 3.2e6).unwrap();
    assert!(base * (k as f64).powi(3) <= 16e6, "the cell cap holds");
    let e = refusal(&mut s, &case, json!({"spacing": 5.0 / k as f64}));
    assert!(e.contains("40 frames") && e.contains("take more than 2000000000 bytes"), "{e}");
}

/// Two directories that name the same time leave the frame order to the file system: refused.
#[test]
fn results_two_directories_of_one_time_are_refused() {
    let (_guard, mut s, case, ..) = written_case(&["0.5", "5e-01"]);
    let e = refusal(&mut s, &case, json!({}));
    assert!(e.contains("two time directories for the same time"), "{e}");
}

/// Processor directories must be numbered from 0 without a gap.
#[test]
fn results_processor_directories_must_be_consecutive() {
    let (_guard, mut s, case, ..) = written_case(&[]);
    for rank in [0, 2] {
        std::fs::create_dir_all(case.join(format!("processor{rank}/0.5"))).unwrap();
    }
    let e = refusal(&mut s, &case, json!({}));
    assert!(e.contains("numbered 0 to N"), "{e}");
}
