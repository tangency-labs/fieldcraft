//! Existing original cross-ply strip: one model, three analyses, version-07 persistence.
//! Independent pure-bending values and bounds come from composite-strip-pure-bending.
use fieldcraft_engine::Session;
use fieldcraft_solver::{MassKind, StepResult};
use serde_json::{Value, json};

fn result(s: &Session, step: &str) -> StepResult {
    let fieldcraft_engine::StepResult::Implicit(r) = s.result(step).unwrap().as_ref() else { panic!("implicit result") };
    r.clone()
}

#[test]
fn cross_ply_static_modes_and_harmonic_keep_model_and_result_semantics() {
    let mut s = Session::default();
    let recipe: Value = serde_json::from_str(include_str!("../../../tests/oracle/composite-strip-pure-bending/recipe.json")).unwrap();
    s.run_script(&recipe).unwrap();
    let static_probe = s.execute("solve.probe", &json!({"step":"Bend","nodes":[18],"elements":[8]})).unwrap();
    let expected: Value = serde_json::from_str(include_str!("../../../tests/oracle/composite-strip-pure-bending/expected.json")).unwrap();
    // The material-axis ply stresses make losing the 90-degree ply direction observable.
    for q in expected["quantities"].as_array().unwrap().iter().take(5) {
        let actual = static_probe.pointer(q["measure"]["pointer"].as_str().unwrap()).unwrap().as_f64().unwrap();
        assert!((actual - q["value"].as_f64().unwrap()).abs() <= q["tolerance"]["absolute"].as_f64().unwrap());
    }
    let StepResult::Static(st) = result(&s, "Bend") else { panic!("static") };
    let mut steps = vec![("Bend".to_owned(), 0)];
    let mut masses = Vec::new();
    for (name, mass) in [("consistent", MassKind::Consistent), ("lumped", MassKind::Lumped)] {
        let modes = format!("Modes_{name}");
        let response = format!("Response_{name}");
        s.execute("model.step.create", &json!({"name":modes,"procedure":{"type":"frequency","modes":2},"boundary_conditions":["clamp"]})).unwrap();
        let reply = s.execute("solve.frequency", &json!({"step":modes,"mass":mass,"threads":1})).unwrap();
        masses.push(reply["total_mass"].as_f64().unwrap());
        assert!(reply["modes"][0]["frequency"].as_f64().unwrap() > 0.);
        s.execute("model.step.create", &json!({"name":response,"procedure":{"type":"harmonic","frequencies":[0.,1.],"method":{"type":"direct"}},"boundary_conditions":["clamp"],"loads":["moment0","moment1","moment2"]})).unwrap();
        s.execute("solve.harmonic", &json!({"step":response,"mass":mass,"threads":1})).unwrap();
        let StepResult::Harmonic(h) = result(&s, &response) else { panic!("harmonic") };
        assert_eq!(h.mass, mass);
        assert_eq!(h.nodes, st.nodes);
        assert_eq!(h.frames.len(), 2);
        for (u, v) in h.frames[0].real.iter().flatten().zip(st.displacement.iter().flatten()) {
            assert!((u - v).abs() < 4e-14);
        }
        assert!(h.frames[0].imaginary.iter().flatten().all(|v| v.abs() < 1e-20));
        for f in &h.frames {
            assert!(f.relative_residual < 1e-8);
        }
        let StepResult::Frequency(f) = result(&s, &modes) else { panic!("frequency") };
        assert_eq!(f.mass_kind, Some(mass));
        assert!(f.section_mass_policy.as_deref().unwrap().contains("declared section"));
        let probe = s.execute("solve.probe", &json!({"step":modes,"nodes":[18]})).unwrap();
        assert_eq!(probe["mass_kind"], json!(mass));
        steps.extend([(modes, 12), (response, 5)]);
    }
    assert!((masses[0] - masses[1]).abs() < 1e-12);
    let definition = s.model.clone();
    let mut saved = Vec::new();
    let mut probes = Vec::new();
    for (step, kind) in &steps {
        let r = result(&s, step);
        use fieldcraft_io::results::{Decoded, Solved, Stored, StoredRef, decode, encode};
        let stored = Solved::Implicit(r.clone());
        let bytes = encode(&[(step, StoredRef { definition: [7; 32], result: stored.as_ref() })]).unwrap();
        assert_eq!(&bytes[..8], b"MCRSLT07");
        // One entry: magic, count, name length, name, definition hash, kind.
        assert_eq!(bytes[12 + 8 + step.len() + 32], *kind);
        assert_eq!(decode(&bytes).unwrap(), Decoded::Results(vec![Stored { definition: [7; 32], result: stored }]));
        let p = if *kind == 0 { json!({"step":step,"nodes":[1,18],"elements":[8]}) } else { json!({"step":step,"nodes":[1,18]}) };
        probes.push((p.clone(), s.execute("solve.probe", &p).unwrap()));
        saved.push(r);
    }
    let path = std::env::temp_dir().join(format!("fieldcraft-integrated-composite-{}.fieldcraft", std::process::id()));
    s.execute("file.save_as", &json!({"path":path,"results":true})).unwrap();
    let mut reopened = Session::default();
    reopened.execute("file.open", &json!({"path":path})).unwrap();
    assert_eq!(reopened.model, definition);
    for ((step, _), old) in steps.iter().zip(saved) {
        assert_eq!(result(&reopened, step), old);
    }
    for (p, old) in probes {
        assert_eq!(reopened.execute("solve.probe", &p).unwrap(), old);
    }
    // The actual solve mass and section policy survive kind12 and native reopening.
    std::fs::remove_file(path).unwrap();
}
