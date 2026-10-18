//! Display systems: which unit each quantity is shown and typed in. Choosing one never changes a
//! stored value.

use serde::{Deserialize, Serialize};

use super::{Quantity, Unit, UnitError};

/// A consistent set of display units.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Preset {
    /// Strict SI: m, kg, s, N, Pa, K.
    #[serde(rename = "SI")]
    Si,
    /// SI with temperatures in °C and angles in degrees.
    #[serde(rename = "m-kg-s")]
    MKgS,
    /// Millimetre, tonne, second, newton, megapascal, °C: what most structural CAE users expect.
    #[default]
    #[serde(rename = "mm-t-s-N-MPa")]
    MmTonneS,
    /// Inch, lbf·s²/in, second, pound-force, psi, °F.
    #[serde(rename = "in-lbf-s-psi")]
    InLbfS,
}

impl Preset {
    pub const ALL: [Preset; 4] = [Preset::Si, Preset::MKgS, Preset::MmTonneS, Preset::InLbfS];

    pub fn label(self) -> &'static str {
        match self {
            Preset::Si => "SI",
            Preset::MKgS => "m-kg-s",
            Preset::MmTonneS => "mm-t-s-N-MPa",
            Preset::InLbfS => "in-lbf-s-psi",
        }
    }

    pub fn from_label(label: &str) -> Option<Preset> {
        Preset::ALL.into_iter().find(|p| p.label() == label)
    }

    /// The unit this preset shows `q` in.
    pub fn symbol(self, q: Quantity) -> &'static str {
        use Quantity as Q;
        match self {
            Preset::Si => q.si_symbol(),
            Preset::MKgS => match q {
                Q::Temperature => "°C",
                Q::Angle => "deg",
                _ => q.si_symbol(),
            },
            Preset::MmTonneS => match q {
                Q::Dimensionless => "1",
                Q::Angle => "deg",
                Q::Length => "mm",
                Q::Area => "mm²",
                Q::Volume => "mm³",
                Q::SecondMomentOfArea => "mm⁴",
                Q::Mass => "t",
                Q::Density => "t/mm³",
                Q::Time => "s",
                Q::Frequency => "Hz",
                Q::Velocity => "mm/s",
                Q::Acceleration => "mm/s²",
                Q::Force => "N",
                Q::Moment => "N·mm",
                Q::Pressure => "MPa",
                Q::Stiffness => "N/mm",
                Q::Damping => "N·s/mm",
                Q::Energy => "mJ",
                Q::Power => "mW",
                Q::DynamicViscosity => "MPa·s",
                Q::KinematicViscosity => "mm²/s",
                Q::Temperature => "°C",
                Q::TemperatureDifference => "K",
                Q::ThermalExpansion => "K⁻¹",
                Q::ThermalConductivity => "mW/(mm·K)",
                Q::SpecificHeat => "mJ/(t·K)",
                Q::HeatTransferCoefficient => "mW/(mm²·K)",
                Q::HeatFlux => "mW/mm²",
                Q::HeatGeneration => "mW/mm³",
            },
            Preset::InLbfS => match q {
                Q::Dimensionless => "1",
                Q::Angle => "deg",
                Q::Length => "in",
                Q::Area => "in²",
                Q::Volume => "in³",
                Q::SecondMomentOfArea => "in⁴",
                Q::Mass => "lbf·s²/in",
                Q::Density => "lbf·s²/in⁴",
                Q::Time => "s",
                Q::Frequency => "Hz",
                Q::Velocity => "in/s",
                Q::Acceleration => "in/s²",
                Q::Force => "lbf",
                Q::Moment => "lbf·in",
                Q::Pressure => "psi",
                Q::Stiffness => "lbf/in",
                Q::Damping => "lbf·s/in",
                Q::Energy => "in·lbf",
                Q::Power => "in·lbf/s",
                Q::DynamicViscosity => "psi·s",
                Q::KinematicViscosity => "in²/s",
                Q::Temperature => "°F",
                Q::TemperatureDifference => "°F",
                Q::ThermalExpansion => "°F⁻¹",
                Q::ThermalConductivity => "lbf/(s·°F)",
                Q::SpecificHeat => "in²/(s²·°F)",
                Q::HeatTransferCoefficient => "lbf/(in·s·°F)",
                Q::HeatFlux => "lbf/(in·s)",
                Q::HeatGeneration => "lbf/(in²·s)",
            },
        }
    }
}

/// The units values are shown and typed in: a preset with optional per-quantity overrides (say
/// mm-t-s-N-MPa with energies in J).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DisplaySystem {
    preset: Preset,
    overrides: Vec<(Quantity, Unit)>,
}

impl DisplaySystem {
    pub fn new(preset: Preset) -> DisplaySystem {
        DisplaySystem { preset, overrides: Vec::new() }
    }

    pub fn preset(&self) -> Preset {
        self.preset
    }

    /// Shows `q` in `unit` instead of the preset's unit; refused when the dimension differs.
    pub fn with_override(mut self, q: Quantity, unit: &str) -> Result<DisplaySystem, UnitError> {
        let unit = Unit::parse(unit)?;
        q.express(0.0, &unit)?;
        self.overrides.retain(|(o, _)| *o != q);
        self.overrides.push((q, unit));
        Ok(self)
    }

    /// The unit `q` is shown and typed in.
    pub fn unit(&self, q: Quantity) -> Result<Unit, UnitError> {
        match self.overrides.iter().find(|(o, _)| *o == q) {
            Some((_, u)) => Ok(u.clone()),
            None => Unit::parse(self.preset.symbol(q)),
        }
    }
}
