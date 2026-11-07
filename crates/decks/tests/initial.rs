//! Initial temperature and stress fields (#160) are not written to any deck yet: every format
//! reports them as unsupported rather than leaving them out silently.

mod common;

use common::{abaqus, lsdyna, nastran};
use fieldcraft_decks::Error;
use fieldcraft_model::{BoundaryCondition, Elastic, ElementKind, InitialField, Material, Model, Procedure, Section, SectionKind, Set, Step, Target};

fn model() -> Model {
    let mut m = Model::default();
    for (l, x) in [(1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0]), (3, [0.0, 1.0, 0.0]), (4, [0.0, 0.0, 1.0])] {
        m.set_node(l, x).unwrap();
    }
    m.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 4]).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_set("base", Set { nodes: vec![1, 2, 3], ..Set::default() }).unwrap();
    m.add_material("Steel", Material { elastic: Some(Elastic::Isotropic { young: 210_000.0, poisson: 0.3 }), ..Material::default() }).unwrap();
    m.add_section(
        "all",
        Section { orientation: None, set: "all".into(), material: "Steel".into(), kind: SectionKind::Solid {}, formulation: Default::default() },
    )
    .unwrap();
    m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
    let mut step = Step::new("Release", Procedure::Static {});
    step.boundary_conditions.insert("fix".into());
    m.add_step(step, None).unwrap();
    m.set_initial_field("Residual", InitialField::Stress { values: [(4, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0])].into(), source: String::new() }).unwrap();
    m
}

fn unsupported(e: std::result::Result<impl std::fmt::Debug, Error>) -> String {
    match e {
        Err(Error::Unsupported { what, .. }) => what,
        other => panic!("expected unsupported, got {other:?}"),
    }
}

#[test]
fn every_format_reports_initial_fields_as_unsupported() {
    let m = model();
    for what in [unsupported(abaqus::write(&m)), unsupported(nastran::write(&m)), unsupported(lsdyna::write(&m))] {
        assert!(what.contains("initial field `Residual`: initial stress fields (not written to decks yet)"), "{what}");
    }
}
