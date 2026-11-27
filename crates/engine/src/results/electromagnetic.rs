//! The electromagnetic results a step draws: a cavity mode's electric field
//! ([`fieldcraft_solver::CavityResult`]), its magnitude or one component, at the nodes.
//!
//! The field of a mode has no scale of its own (any multiple of it is the same mode): the solver
//! gives it with its largest |E| equal to 1, so these values have no unit.

use fieldcraft_solver::CavityResult;

use super::{Pick, magnitude};
use crate::{EngineError, Result};

/// An electric field quantity at the nodes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ElectricField {
    Magnitude,
    X,
    Y,
    Z,
}

impl ElectricField {
    pub const ALL: [ElectricField; 4] = [ElectricField::Magnitude, ElectricField::X, ElectricField::Y, ElectricField::Z];

    /// The id commands take.
    pub fn id(self) -> &'static str {
        match self {
            ElectricField::Magnitude => "electric_field",
            ElectricField::X => "electric_field_x",
            ElectricField::Y => "electric_field_y",
            ElectricField::Z => "electric_field_z",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            ElectricField::Magnitude => "Electric field magnitude",
            ElectricField::X => "Electric field X",
            ElectricField::Y => "Electric field Y",
            ElectricField::Z => "Electric field Z",
        }
    }

    /// Its value at every node of the cavity result, in the result's node order, for the mode
    /// `pick` shows.
    pub fn values(self, c: &CavityResult, pick: Pick) -> Result<Vec<f64>> {
        let m = pick
            .mode
            .checked_sub(1)
            .and_then(|k| c.modes.get(k))
            .ok_or_else(|| EngineError::Params(format!("step `{}` has modes 1 to {}; there is no mode {}", c.step, c.modes.len(), pick.mode)))?;
        let of = |e: &[f64; 3]| match self {
            ElectricField::Magnitude => magnitude(*e),
            ElectricField::X => e[0],
            ElectricField::Y => e[1],
            ElectricField::Z => e[2],
        };
        Ok(m.field.iter().map(of).collect())
    }
}
