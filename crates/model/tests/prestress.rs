//! Equilibrium dependencies must survive persistence and participate in result invalidation.
use fieldcraft_model::{Model, Procedure, Step};

fn modes(base: &str) -> Step {
    let mut step = Step::new("Modes", Procedure::Frequency { modes: 2, min_hz: None, max_hz: None });
    step.base_state = Some(base.into());
    step
}

#[test]
fn equilibrium_dependency_is_ordered_and_round_trips() {
    let mut model = Model::default();
    assert!(model.add_step(modes("Load"), None).is_err());
    model.add_step(Step::new("Load", Procedure::Static {}), None).unwrap();
    assert!(model.add_step(modes("Load"), Some("Load")).is_err());
    model.add_step(modes("Load"), None).unwrap();
    assert!(model.remove_step("Load").is_err());
    assert!(model.replace_step(Step::new("Load", Procedure::Frequency { modes: 1, min_hz: None, max_hz: None })).is_err());
    let stored = serde_json::to_vec(&model).unwrap();
    let restored: Model = serde_json::from_slice(&stored).unwrap();
    assert_eq!(restored.perturbation_base(restored.step("Modes").unwrap()).unwrap().unwrap().name, "Load");
    let before = serde_json::to_value(restored.step_inputs("Modes").unwrap()).unwrap();
    let mut base = model.step("Load").unwrap().clone();
    base.output.clear();
    model.replace_step(base).unwrap();
    let after = serde_json::to_value(model.step_inputs("Modes").unwrap()).unwrap();
    assert_ne!(before, after, "the equilibrium definition is part of the perturbation fingerprint");
    model.remove_step("Modes").unwrap();
    model.remove_step("Load").unwrap();
}
