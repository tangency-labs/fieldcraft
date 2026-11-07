//! Until a format carries composite data, it must refuse that data by definition name.
use fieldcraft_decks::{Error, FORMATS};
use fieldcraft_model::units::Preset;
use fieldcraft_model::*;

fn plate() -> Model {
    let mut m = Model::default();
    for (i, p) in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]].into_iter().enumerate() {
        m.set_node(i as u32 + 1, p).unwrap();
    }
    m.set_element(1, ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
    m.add_set("all", Set { elements: vec![1], ..Set::default() }).unwrap();
    m.add_material(
        "Carbon",
        Material { elastic: Some(Elastic::Isotropic { young: 70e9, poisson: 0.25 }), density: Some(1600.0), ..Material::default() },
    )
    .unwrap();
    m.add_section(
        "plate",
        Section {
            set: "all".into(),
            material: "Carbon".into(),
            kind: SectionKind::Shell { thickness: 0.001, offset: 0.0, points: 3 },
            orientation: None,
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    m.add_boundary_condition("hold", BoundaryCondition::Fixed { target: Target::Set("all".into()) }).unwrap();
    let mut step = Step::new("Load", Procedure::Static {});
    step.boundary_conditions.insert("hold".into());
    m.add_step(step, None).unwrap();
    m
}

#[test]
fn each_writer_refuses_the_named_data_instead_of_silently_dropping_it() {
    let base = plate();
    let mut oriented = base.clone();
    let mut section = oriented.section("plate").unwrap().clone();
    section.orientation = Some(SectionOrientation {
        frame: MaterialOrientation::Rectangular { axis1: [1.0, 1.0, 0.0], axis2: [-1.0, 1.0, 0.0] },
        shell_reference_axis: 1,
    });
    oriented.replace_section("plate", section).unwrap();
    let mut composite = base.clone();
    let mut section = composite.section("plate").unwrap().clone();
    section.material.clear();
    section.kind = SectionKind::CompositeShell {
        layup: Layup { plies: vec![Ply { material: "Carbon".into(), thickness: 0.001, angle: 0.0, points: 3 }], offset: 0.0 },
    };
    composite.replace_section("plate", section).unwrap();
    let mut strengths = base;
    let mut material = strengths.material("Carbon").unwrap().clone();
    material.ply_strength =
        Some(PlyStrength { tension: [1e9, 50e6], compression: [800e6, 100e6], shear12: 60e6, shear23: 40e6, tsai_wu_interaction: -0.5 });
    strengths.replace_material("Carbon", material).unwrap();
    for (model, named) in [
        (oriented, "material orientation of section `plate`"),
        (composite, "composite shell section `plate`"),
        (strengths, "ply strengths of material `Carbon`"),
    ] {
        for format in FORMATS {
            if format.name == "abaqus" && !named.starts_with("ply strengths") {
                continue;
            }
            if format.name == "lsdyna" && !named.starts_with("ply strengths") {
                continue;
            }
            if format.name == "nastran" && named.starts_with("material orientation") {
                continue;
            }
            match format.write(&model, Preset::MKgS) {
                Err(Error::Unsupported { what, .. }) => assert!(what.contains(named), "{}: {what}", format.name),
                other => panic!("{} must name unsupported data: {other:?}", format.name),
            }
        }
    }
}
