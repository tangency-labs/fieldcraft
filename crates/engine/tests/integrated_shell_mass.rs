//! Project-level persistence/re-execution using original Navier and tapered-solid models.
//! Physics references remain in the existing oracle and section-mass integration tests.
use fieldcraft_engine::{Session, commands, schema};
use serde_json::{Value, json};

fn project_chain(name: &str, definitions: Vec<Value>, calls: Vec<Value>) {
    let mut s = Session::default();
    s.run_script(&json!(definitions)).unwrap();
    let mut replies = serde_json::Map::new();
    for call in &calls {
        let reply = s.execute(call["command"].as_str().unwrap(), &call["params"]).unwrap();
        replies.insert(call["params"]["step"].as_str().unwrap().to_owned(), reply);
    }
    if name == "mixed-mass" {
        let c = replies["Modes-consistent"]["modes"][0]["frequency"].as_f64().unwrap();
        let l = replies["Modes-lumped"]["modes"][0]["frequency"].as_f64().unwrap();
        assert!((c - l).abs() > 1e-6, "policies must affect the actual operator");
        assert!((replies["Modes-consistent"]["total_mass"].as_f64().unwrap() - 3.0).abs() < 1e-12);
        assert!((replies["Modes-lumped"]["total_mass"].as_f64().unwrap() - 3.0).abs() < 1e-12);
    }
    let model = s.model.clone();
    let results: Vec<_> = calls
        .iter()
        .map(|c| (c["params"]["step"].as_str().unwrap().to_owned(), s.result(c["params"]["step"].as_str().unwrap()).unwrap().as_ref().clone()))
        .collect();
    let dir = std::env::temp_dir().join(format!("fieldcraft-project-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("project.fieldcraft");
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut open = Session::default();
    open.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(open.model, model, "all sections, orientations, sets, loads and steps");
    for (step, result) in &results {
        assert_eq!(open.result(step).unwrap().as_ref(), result, "stored {step}");
    }
    // Reverse execution detects accidental dependence on a previously active result or mass choice.
    for call in calls.iter().rev() {
        let step = call["params"]["step"].as_str().unwrap();
        let reply = open.execute(call["command"].as_str().unwrap(), &call["params"]).unwrap();
        println!("{name} {step}: {reply}");
        assert_eq!(open.model, model);
        for (old_step, result) in &results {
            assert_eq!(open.result(old_step).unwrap().as_ref(), result, "after {step}: {old_step}");
        }
    }
    if let Ok(out) = std::env::var("FIELDCRAFT_TEST_ARTIFACTS") {
        let out = std::path::Path::new(&out).join(name);
        std::fs::create_dir_all(&out).unwrap();
        std::fs::copy(&path, out.join("project.fieldcraft")).unwrap();
        std::fs::write(out.join("definitions.json"), serde_json::to_vec_pretty(&definitions).unwrap()).unwrap();
        std::fs::write(out.join("replies.json"), serde_json::to_vec_pretty(&replies).unwrap()).unwrap();
        std::fs::write(out.join("solve-calls.json"), serde_json::to_vec_pretty(&calls).unwrap()).unwrap();
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn composite_mitc8_multiple_analysis_results_reopen_and_reexecute() {
    let frequency: Vec<Value> = serde_json::from_str(include_str!("../../../tests/oracle/mitc8-navier-frequency/recipe.json")).unwrap();
    let harmonic: Vec<Value> = serde_json::from_str(include_str!("../../../tests/oracle/mitc8-navier-harmonic/recipe.json")).unwrap();
    let common = |r: &[Value]| {
        r.iter().filter(|c| c["command"] != "model.step.create" && !c["command"].as_str().unwrap().starts_with("solve.")).cloned().collect::<Vec<_>>()
    };
    assert_eq!(common(&frequency), common(&harmonic));
    let mut definitions = common(&frequency);
    let mut calls = Vec::new();
    for c in frequency.iter().chain(&harmonic) {
        match c["command"].as_str().unwrap() {
            "model.step.create" => definitions.push(c.clone()),
            "solve.frequency" | "solve.harmonic" => calls.push(c.clone()),
            _ => {}
        }
    }
    let mut stat = harmonic.iter().find(|c| c["command"] == "model.step.create").unwrap()["params"].clone();
    stat["name"] = json!("Static");
    stat["procedure"] = json!({"type":"static"});
    definitions.push(json!({"command":"model.step.create","params":stat}));
    calls.push(json!({"command":"solve.static","params":{"step":"Static","threads":1}}));
    project_chain("mitc8", definitions, calls);
}

#[test]
fn mixed_section_mass_policies_survive_multi_analysis_reexecution() {
    // Two independent tapered hexes: determinant 1+z; first row-sum, second consistent.
    // This is a project identity test; exact integration physics is in section_mass_tests.
    let mut d = vec![
        json!({"command":"model.units","params":{"system":"m-kg-s-N-Pa"}}),
        json!({"command":"model.material.create","params":{"name":"M","elastic":{"type":"isotropic","young":1000,"poisson":0},"density":1}}),
    ];
    for (part, policy) in ["row_sum", "consistent"].into_iter().enumerate() {
        let shift = 8 * part;
        for (i, xyz) in
            [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.], [0., 0., 1.], [2., 0., 1.], [2., 1., 1.], [0., 1., 1.]].into_iter().enumerate()
        {
            d.push(json!({"command":"mesh.node.add","params":{"label":shift+i+1,"xyz":[xyz[0]+3.0*part as f64,xyz[1],xyz[2]]}}));
        }
        d.extend([
            json!({"command":"mesh.element.add","params":{"label":part+1,"kind":"hex8","nodes":(shift+1..=shift+8).collect::<Vec<_>>()}}),
            json!({"command":"model.set.create","params":{"name":policy,"elements":[part+1]}}),
            json!({"command":"model.section.create","params":{"name":policy,"set":policy,"material":"M","kind":{"type":"solid"},"formulation":{"lumped_mass":policy}}}),
        ]);
    }
    d.extend([
        json!({"command":"model.set.create","params":{"name":"base","nodes":[1,2,3,4,9,10,11,12]}}),
        json!({"command":"model.set.create","params":{"name":"top","nodes":[5,6,7,8,13,14,15,16]}}),
        json!({"command":"model.boundary_condition.create","params":{"name":"fix","type":"displacement","target":{"set":"base"},"values":[0,0,0,null,null,null]}}),
        json!({"command":"model.load.create","params":{"name":"force","type":"force","target":{"set":"top"},"vector":[1,0,0],"distribution":"per_node"}}),
        json!({"command":"model.step.create","params":{"name":"Static","procedure":{"type":"static"},"boundary_conditions":["fix"],"loads":["force"]}}),
    ]);
    let mut calls = vec![json!({"command":"solve.static","params":{"step":"Static","threads":1}})];
    for mass in ["consistent", "lumped"] {
        for (name, command, procedure) in [
            (format!("Modes-{mass}"), "solve.frequency", json!({"type":"frequency","modes":3})),
            (
                format!("Response-{mass}"),
                "solve.harmonic",
                json!({"type":"harmonic","frequencies":[0.1,0.2],"method":{"type":"direct"},"damping":{"alpha":0.2}}),
            ),
        ] {
            d.push(json!({"command":"model.step.create","params":{"name":name,"procedure":procedure,"boundary_conditions":["fix"],"loads":if command=="solve.harmonic" {vec!["force"]}else{vec![]}}}));
            calls.push(json!({"command":command,"params":{"step":name,"mass":mass,"threads":1}}));
        }
    }
    project_chain("mixed-mass", d, calls);
}

#[test]
fn command_forms_expose_supported_shell_and_mass_options() {
    for (id, choices) in [
        ("mesh.element.add", vec!["quad8", "tri6", "beam2", "hex8"]),
        ("model.section.create", vec!["composite_shell", "solid_shell", "axisymmetric_shell", "lumped_mass", "row_sum", "consistent"]),
        ("model.step.create", vec!["static", "frequency", "harmonic", "buckle"]),
        ("solve.frequency", vec!["consistent", "lumped"]),
        ("solve.harmonic", vec!["consistent", "lumped"]),
    ] {
        let spec = commands().iter().find(|c| c.id == id).unwrap();
        let form = schema::json_schema(spec.form.unwrap()).to_string();
        for choice in choices {
            assert!(form.contains(&format!("\"{choice}\"")), "{id} missing {choice}");
        }
    }
}
