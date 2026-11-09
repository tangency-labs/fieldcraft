//! Native kind11 rod history -> reopen -> actual nominal stress -> rainflow/Miner.
//! The finite-strain constitutive identity is P=E*(lambda^2-1)*lambda/2.
use fieldcraft_engine::Session;
use serde_json::{Value, json};

fn probe(s: &mut Session) -> Value {
    s.execute("solve.probe", &json!({"step":"Pulse","nodes":[2],"elements":[1]})).unwrap()
}

#[test]
fn reopened_nonlinear_rod_nominal_history_drives_fatigue() {
    let mut s = Session::default();
    s.run_script(&serde_json::from_str::<Value>(include_str!("data/transient-rod.json")).unwrap()).unwrap();
    s.execute("model.step.edit", &json!({"name":"Pulse","loads":[],"procedure":{"type":"transient","time_step":0.01,"increments":10,"output_every":1,"initial_displacement":{"2":[0.25,0,0]}}})).unwrap();
    s.execute("solve.transient.nonlinear", &json!({"step":"Pulse","mass":"consistent","threads":1})).unwrap();
    let original = probe(&mut s);
    let saved = s.result("Pulse").unwrap().clone();
    let bytes =
        fieldcraft_io::results::encode(&[("Pulse", fieldcraft_io::results::StoredRef { definition: [7; 32], result: saved.as_stored() })]).unwrap();
    assert_eq!(&bytes[..8], b"MCRSLT07");
    assert_eq!(bytes[12 + 8 + "Pulse".len() + 32], 11);
    let path = std::env::temp_dir().join(format!("fieldcraft-kind11-fatigue-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut reopened = Session::default();
    reopened.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(reopened.model, s.model);
    assert_eq!(reopened.result("Pulse").unwrap().as_ref(), saved.as_ref());
    let history = probe(&mut reopened);
    assert_eq!(history, original);
    assert_eq!(history["units"], "SI");
    assert_eq!(history["kind"], 11);
    let element = &history["elements"][0];
    assert_eq!(element["label"], 1);
    assert_eq!(element["reference_coordinates"], json!([[0., 0., 0.], [1., 0., 0.]]));
    assert_eq!(element["reference_area"], 1.);
    let frames = history["frames"].as_array().unwrap();
    assert_eq!(frames.len(), 11);
    let mut samples = Vec::new();
    for (i, f) in frames.iter().enumerate() {
        let time = f["time"].as_f64().unwrap();
        assert!((time - i as f64 * 0.01).abs() < 1e-14);
        assert!(f["relative_residual"].as_f64().unwrap() < 1e-8);
        let stretch = 1. + f["displacement"][1][0].as_f64().unwrap();
        let stress = &f["rods"][0];
        let pk2 = 2. * (stretch * stretch - 1.);
        assert!((stress["second_piola_stress"].as_f64().unwrap() - pk2).abs() < 1e-12);
        assert!((stress["nominal_stress"].as_f64().unwrap() - pk2 * stretch).abs() < 1e-12);
        samples.push(json!({"time":f["time"],"stress":stress["nominal_stress"]}));
    }
    assert!((samples[0]["stress"].as_f64().unwrap() - 1.40625).abs() < 1e-12);
    assert!(samples.windows(2).all(|w| w[1]["stress"].as_f64().unwrap() < w[0]["stress"].as_f64().unwrap()));
    let p = json!({"source_id":"original kind11 rod, reopened native history, nominal stress","step":"Pulse","max_time_step_seconds":0.01,
        "axial_history":{"step":history["step"],"element":element["label"],"location":"constant_rod2_midpoint","stress_measure":"nominal_axial_reference_area","axis":"material_following_axial","reference_material_axis":[1,0,0],"reference_area_m2":element["reference_area"],"stress_unit":"Pa","time_unit":"s","samples":samples},
        "curve_stress_measure":"nominal_axial_reference_area","calibration_source":"synthetic nominal S-N verification, not measured material",
        "curve":{"amplitude_ref":1,"cycles_ref":100,"exponent":2,"endurance_limit":0},"correction":{"type":"none"}});
    let fatigue = reopened.execute("solve.fatigue.transient", &p).unwrap();
    let amplitude = (samples[0]["stress"].as_f64().unwrap() - samples.last().unwrap()["stress"].as_f64().unwrap()) / 2.;
    assert_eq!(fatigue["cycle_count"], 0.5);
    assert_eq!(fatigue["cycles"].as_array().unwrap().len(), 1);
    assert!((fatigue["damage"].as_f64().unwrap() - 0.5 * amplitude * amplitude / 100.).abs() < 1e-14);
    assert!((fatigue["duration_seconds"].as_f64().unwrap() - 0.1).abs() < 1e-14);
    assert_eq!(fatigue["persisted"], false);
    for (pointer, value) in
        [("/axial_history/stress_measure", json!("PK2")), ("/axial_history/location", json!("global nodal average; shell top surface"))]
    {
        let mut bad = p.clone();
        *bad.pointer_mut(pointer).unwrap() = value;
        assert!(reopened.execute("solve.fatigue.transient", &bad).is_err());
    }
    reopened.execute("model.units", &json!({"system":"mm-kg-ms-kN-GPa"})).unwrap();
    assert_eq!(probe(&mut reopened), history);
    assert_eq!(reopened.execute("solve.fatigue.transient", &p).unwrap(), fatigue);
    if let Ok(folder) = std::env::var("FIELDCRAFT_TEST_ARTIFACTS") {
        let folder = std::path::Path::new(&folder);
        std::fs::create_dir_all(folder).unwrap();
        std::fs::copy(&path, folder.join("kind11.fieldcraft")).unwrap();
        for (name, value) in [("probe", history), ("fatigue-input", p), ("fatigue", fatigue)] {
            std::fs::write(folder.join(format!("{name}.json")), serde_json::to_vec_pretty(&value).unwrap()).unwrap();
        }
    }
    std::fs::remove_file(path).unwrap();
}
