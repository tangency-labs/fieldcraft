//! Units: the unit system label a model carries today, and the SI unit table every conversion
//! goes through.
//!
//! **The table** ([`Dimension`], [`Unit`], [`Quantity`], [`DisplaySystem`]) is the one place
//! FieldCraft writes down a conversion factor. Each [`Quantity`] has a dimension and an SI unit;
//! a [`Unit`] parses expressions such as `MPa`, `t/mm³` or `W/(m·K)` to a dimension and a factor
//! to SI; [`Quantity::parse`] reads `"210 GPa"` or `"20 °C"` with a dimension check, and
//! [`Quantity::format`] shows an SI value in a [`DisplaySystem`]. Absolute temperatures
//! ([`Quantity::Temperature`]) take the °C/°F offset; temperature differences and every compound
//! unit never do. Issue #579 moves model storage to SI on top of this table; until then the
//! model still stores numbers as entered.
//!
//! **The label** ([`UnitSystem`]): the model never converts yet: numbers are stored as entered.
//! The label records which consistent system they are in, so the check can warn when a magnitude
//! is implausible for it (a steel modulus of 210 000 in a metre–pascal model, a density of 7850
//! in a millimetre–tonne one), and so writers can say which system a deck is in.

mod dimension;
mod display;
mod quantity;
mod unit;

pub use dimension::Dimension;
pub use display::{DisplaySystem, Preset};
pub use quantity::{Quantity, ValueInput, WithUnit, format_number};
pub use unit::Unit;

/// Why a unit or a value with a unit was refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum UnitError {
    #[error("no value or unit given")]
    Empty,
    #[error("unit text is too long")]
    TooLong,
    #[error("`{0}` is not a finite number")]
    NotFinite(String),
    #[error("unknown unit `{0}`")]
    UnknownUnit(String),
    #[error("cannot read `{0}` as a value with a unit (for example `210 GPa`, `t/mm³`, `W/(m·K)`)")]
    Syntax(String),
    #[error("`{0}` is out of range")]
    OutOfRange(String),
    #[error("`{0}` is below absolute zero")]
    BelowAbsoluteZero(String),
    #[error("`{unit}` is not a unit of {quantity} (angles take rad, deg or °; ratios take 1 or %)")]
    WrongKind { unit: String, quantity: &'static str },
    #[error("`{unit}` is {found}, but a {quantity} is {expected}")]
    WrongDimension { unit: String, found: Dimension, quantity: &'static str, expected: Dimension },
}

use serde::{Deserialize, Serialize};

/// A consistent unit system: length, mass, time, force, stress.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum UnitSystem {
    /// Not chosen: magnitudes are not checked.
    #[default]
    #[serde(rename = "unspecified")]
    Unspecified,
    /// Metre, kilogram, second, newton, pascal.
    #[serde(rename = "m-kg-s-N-Pa")]
    MKgS,
    /// Millimetre, tonne, second, newton, megapascal: the usual choice for parts in millimetres.
    #[serde(rename = "mm-t-s-N-MPa")]
    MmTonneS,
    /// Millimetre, kilogram, millisecond, kilonewton, gigapascal.
    #[serde(rename = "mm-kg-ms-kN-GPa")]
    MmKgMs,
    /// Millimetre, gram, millisecond, newton, megapascal.
    #[serde(rename = "mm-g-ms-N-MPa")]
    MmGMs,
    /// Inch, pound-force second squared per inch, second, pound-force, psi.
    #[serde(rename = "in-lbf-s-psi")]
    InLbfS,
}

impl UnitSystem {
    pub const ALL: [UnitSystem; 6] =
        [UnitSystem::Unspecified, UnitSystem::MKgS, UnitSystem::MmTonneS, UnitSystem::MmKgMs, UnitSystem::MmGMs, UnitSystem::InLbfS];

    /// The label used in commands and files (`mm-t-s-N-MPa`, …).
    pub fn label(self) -> &'static str {
        match self {
            UnitSystem::Unspecified => "unspecified",
            UnitSystem::MKgS => "m-kg-s-N-Pa",
            UnitSystem::MmTonneS => "mm-t-s-N-MPa",
            UnitSystem::MmKgMs => "mm-kg-ms-kN-GPa",
            UnitSystem::MmGMs => "mm-g-ms-N-MPa",
            UnitSystem::InLbfS => "in-lbf-s-psi",
        }
    }

    pub fn from_label(label: &str) -> Option<UnitSystem> {
        UnitSystem::ALL.into_iter().find(|u| u.label() == label)
    }

    /// The length unit (`m`, `mm`, `in`), or `None` when unspecified.
    pub fn length(self) -> Option<&'static str> {
        match self {
            UnitSystem::Unspecified => None,
            UnitSystem::MKgS => Some("m"),
            UnitSystem::MmTonneS | UnitSystem::MmKgMs | UnitSystem::MmGMs => Some("mm"),
            UnitSystem::InLbfS => Some("in"),
        }
    }

    /// Size of one stress unit in pascals, one density unit in kg/m³ and one acceleration unit in
    /// m/s², or `None` when unspecified.
    fn scales(self) -> Option<Scales> {
        const IN: f64 = 0.0254;
        const LBF: f64 = 4.448_221_615_260_5;
        Some(match self {
            UnitSystem::Unspecified => return None,
            UnitSystem::MKgS => Scales { stress: 1.0, density: 1.0, acceleration: 1.0 },
            UnitSystem::MmTonneS => Scales { stress: 1e6, density: 1e12, acceleration: 1e-3 },
            UnitSystem::MmKgMs => Scales { stress: 1e9, density: 1e9, acceleration: 1e3 },
            UnitSystem::MmGMs => Scales { stress: 1e6, density: 1e6, acceleration: 1e3 },
            // Mass unit lbf·s²/in = LBF / IN kg; density unit = that per cubic inch.
            UnitSystem::InLbfS => Scales { stress: LBF / (IN * IN), density: LBF / IN / (IN * IN * IN), acceleration: IN },
        })
    }

    /// A warning when a Young's or shear modulus is implausible for any engineering material in
    /// this system (outside 10 kPa – 2 TPa).
    pub fn check_modulus(self, what: &str, value: f64) -> Option<String> {
        self.check(what, value, |s| s.stress, (1e4, 2e12), "Pa")
    }

    /// A warning when a mass density is implausible (outside 1 – 25 000 kg/m³).
    pub fn check_density(self, value: f64) -> Option<String> {
        self.check("density", value, |s| s.density, (1.0, 25_000.0), "kg/m³")
    }

    /// A warning when a gravitational acceleration is implausible (outside 0.1 – 100 m/s²).
    pub fn check_gravity(self, magnitude: f64) -> Option<String> {
        self.check("gravity", magnitude, |s| s.acceleration, (0.1, 100.0), "m/s²")
    }

    fn check(self, what: &str, value: f64, scale: fn(&Scales) -> f64, (lo, hi): (f64, f64), si: &str) -> Option<String> {
        let s = self.scales()?;
        let si_value = value.abs() * scale(&s);
        if (lo..=hi).contains(&si_value) {
            return None;
        }
        let k = scale(&s);
        Some(format!(
            "{what} {value} is {si_value:.3e} {si}, outside the plausible {lo:.0e}–{hi:.0e} {si} ({:.3e}–{:.3e} in {}); check the unit system",
            lo / k,
            hi / k,
            self.label()
        ))
    }
}

struct Scales {
    stress: f64,
    density: f64,
    acceleration: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steel_is_plausible_in_every_system_and_wrong_systems_warn() {
        // Steel: E = 210 GPa, ρ = 7850 kg/m³, g = 9.81 m/s², in each system's own units.
        let steel = [
            (UnitSystem::MKgS, 2.1e11, 7850.0, 9.81),
            (UnitSystem::MmTonneS, 2.1e5, 7.85e-9, 9810.0),
            (UnitSystem::MmKgMs, 210.0, 7.85e-6, 9.81e-3),
            (UnitSystem::MmGMs, 2.1e5, 7.85e-3, 9.81e-3),
            (UnitSystem::InLbfS, 3.05e7, 7.35e-4, 386.1),
        ];
        for (u, e, rho, g) in steel {
            assert_eq!(u.check_modulus("E", e), None, "{u:?}");
            assert_eq!(u.check_density(rho), None, "{u:?}");
            assert_eq!(u.check_gravity(g), None, "{u:?}");
            assert_eq!(UnitSystem::from_label(u.label()), Some(u));
            assert_eq!(serde_json::to_value(u).unwrap(), serde_json::json!(u.label()));
        }
        // SI numbers in a millimetre–tonne model.
        assert!(UnitSystem::MmTonneS.check_modulus("E", 2.1e11).is_some());
        assert!(UnitSystem::MmTonneS.check_density(7850.0).is_some());
        assert!(UnitSystem::MKgS.check_gravity(9810.0).is_some());
        assert_eq!(UnitSystem::Unspecified.check_density(1e30), None);
        assert!(UnitSystem::MKgS.check_density(f64::NAN).is_some());
    }
}
