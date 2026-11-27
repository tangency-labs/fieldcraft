//! No deck format describes electromagnetics: a material's permittivity and permeability are
//! refused by every writer rather than left out silently, even in a model with no
//! electromagnetic step.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::Error;
use fieldcraft_model::{
    BoundaryCondition, Elastic, Electromagnetic, ElementKind, Material, Model, Procedure, Section, SectionKind, Set, Step, Target,
};

fn model(electromagnetic: Option<Electromagnetic>) -> Model {
    let mut m = Model::default();
    for (l, x) in [(1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]), (3, [0.0, 1.0, 0.0]), (4, [0.0, 0.0, 1.0])] {
        m.set_node(l, x).unwrap();
    }
    m.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 4]).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_set("base", Set { nodes: vec![1, 2, 3], ..Set::default() }).unwrap();
    let steel =
        Material { elastic: Some(Elastic::Isotropic { young: 210e9, poisson: 0.3 }), density: Some(7850.0), electromagnetic, ..Material::default() };
    m.add_material("Steel", steel).unwrap();
    m.add_section("all", Section { set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Default::default() })
        .unwrap();
    m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    let mut step = Step::new("Load", Procedure::Static {});
    step.boundary_conditions.insert("fix".into());
    m.add_step(step, None).unwrap();
    m
}

#[test]
fn every_format_refuses_a_material_with_electromagnetic_properties() {
    let m = model(Some(Electromagnetic { permittivity: 1.0, permeability: 100.0 }));
    for (format, result) in
        [("abaqus", abaqus::write(&m).map(drop)), ("nastran", nastran::write(&m).map(drop)), ("lsdyna", lsdyna::write(&m).map(drop))]
    {
        match result {
            Err(Error::Unsupported { what, .. }) => {
                assert!(what.contains("material `Steel`: electromagnetic properties (no deck format describes them)"), "{format}: {what}")
            }
            other => panic!("{format}: expected unsupported, got {other:?}"),
        }
    }
    // Without them the model is written.
    let m = model(None);
    for (format, result) in
        [("abaqus", abaqus::write(&m).map(drop)), ("nastran", nastran::write(&m).map(drop)), ("lsdyna", lsdyna::write(&m).map(drop))]
    {
        assert!(result.is_ok(), "{format}: {result:?}");
    }
}
