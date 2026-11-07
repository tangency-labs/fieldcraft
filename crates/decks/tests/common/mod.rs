//! The test models are built with millimetre–tonne–second numbers, as decks show them; the
//! writers take models stored in SI. These wrappers store a test model in SI and write it in
//! mm-t-s-N-MPa, so every number in a deck is the one the test typed.
#![allow(dead_code, unused_imports)]

use fieldcraft_model::Model;
use fieldcraft_model::units::{DisplaySystem, Preset};

/// A model built in mm-t-s-N-MPa numbers, stored in SI.
pub fn si(m: &Model) -> Model {
    let mm = DisplaySystem::new(Preset::MmTonneS);
    m.converted(&|q, v| mm.to_si(q, v).ok()).unwrap()
}

pub mod abaqus {
    pub use fieldcraft_decks::abaqus::*;

    pub fn write(m: &fieldcraft_model::Model) -> fieldcraft_decks::Result<String> {
        fieldcraft_decks::abaqus::write(&super::si(m), super::Preset::MmTonneS)
    }
}

pub mod nastran {
    pub use fieldcraft_decks::nastran::*;

    pub fn write(m: &fieldcraft_model::Model) -> fieldcraft_decks::Result<String> {
        fieldcraft_decks::nastran::write(&super::si(m), super::Preset::MmTonneS)
    }
}

pub mod lsdyna {
    pub use fieldcraft_decks::lsdyna::*;

    pub fn write(m: &fieldcraft_model::Model) -> fieldcraft_decks::Result<fieldcraft_decks::lsdyna::Written> {
        fieldcraft_decks::lsdyna::write(&super::si(m), super::Preset::MmTonneS)
    }

    /// A model built in SI numbers, written in SI.
    pub fn write_si(m: &fieldcraft_model::Model) -> fieldcraft_decks::Result<fieldcraft_decks::lsdyna::Written> {
        fieldcraft_decks::lsdyna::write(m, super::Preset::MKgS)
    }
}

pub fn laminate_model() -> Model {
    use fieldcraft_model::{
        Elastic, ElementKind, Formulation, Layup, Material, MaterialOrientation, Model, Ply, Procedure, Section, SectionKind, SectionOrientation,
        Set, Step,
    };
    let mut model = Model::default();
    for (i, p) in [[0., 0., 0.], [2., 1., 0.], [3., 3., 0.], [0., 2., 0.]].into_iter().enumerate() {
        model.set_node(i as u32 + 1, p).unwrap();
    }
    model.set_element(1, ElementKind::Quad4, vec![1, 2, 3, 4]).unwrap();
    model.add_set("plate", Set { elements: vec![1], nodes: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
    model
        .add_material(
            "ply",
            Material {
                elastic: Some(Elastic::EngineeringConstants { young: [120e9, 9e9, 7e9], poisson: [0.25, 0.2, 0.3], shear: [5e9, 4e9, 3e9] }),
                density: Some(1600.),
                ..Material::default()
            },
        )
        .unwrap();
    model
        .add_section(
            "laminate",
            Section {
                set: "plate".into(),
                material: String::new(),
                formulation: Formulation::default(),
                kind: SectionKind::CompositeShell {
                    layup: Layup {
                        offset: 0.25,
                        plies: vec![
                            Ply { material: "ply".into(), thickness: 0.001, angle: std::f64::consts::FRAC_PI_6, points: 1 },
                            Ply { material: "ply".into(), thickness: 0.002, angle: -std::f64::consts::FRAC_PI_4, points: 1 },
                        ],
                    },
                },
                orientation: Some(SectionOrientation {
                    frame: MaterialOrientation::Rectangular { axis1: [0., 1., 0.], axis2: [-1., 0., 0.] },
                    shell_reference_axis: 1,
                }),
            },
        )
        .unwrap();
    model.add_step(Step::new("load", Procedure::Static {}), None).unwrap();
    model
}
