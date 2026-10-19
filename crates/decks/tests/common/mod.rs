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
