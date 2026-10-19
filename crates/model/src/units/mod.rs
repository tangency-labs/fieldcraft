//! Units: the SI unit table every conversion goes through.
//!
//! The model stores every quantity in SI (m, kg, s, N, Pa, K, rad). This table is the one place
//! FieldCraft writes down a conversion factor: each [`Quantity`] has a dimension and an SI unit;
//! a [`Unit`] parses expressions such as `MPa`, `t/mm³` or `W/(m·K)` to a dimension and a factor
//! to SI; [`Quantity::parse`] reads `"210 GPa"` or `"20 °C"` with a dimension check, and
//! [`Quantity::format`] shows an SI value in a [`DisplaySystem`] (a consistent [`Preset`] with
//! optional per-quantity overrides). Absolute temperatures ([`Quantity::Temperature`]) take the
//! °C/°F offset; temperature differences and every compound unit never do. Conversion happens
//! only at the boundaries (commands, forms, replies, files of other formats); docs/units.md.

mod decimal;
mod dimension;
mod display;
mod quantity;
mod unit;

pub use dimension::Dimension;
pub use display::{DisplaySystem, Preset};
pub use quantity::{Quantity, ValueInput, WithUnit, format_number};
pub use unit::{CELSIUS_ZERO, Unit};

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
    #[error(
        "`{unit}` is not a unit of {quantity} (angles, angular velocities and torsional springs take rad, deg, ° or rpm; no other quantity does)"
    )]
    WrongKind { unit: String, quantity: &'static str },
    #[error("`{unit}` is {found}, but a {quantity} is {expected}")]
    WrongDimension { unit: String, found: Dimension, quantity: &'static str, expected: Dimension },
}

/// `kelvin` less `zero`, in decimal when both are short decimals (so 293.15 K less 273.15 is 20).
pub fn kelvin_above(kelvin: f64, zero: f64) -> f64 {
    decimal::add(kelvin, -zero)
}

/// `above` plus `zero` (the inverse of [`kelvin_above`]).
pub fn kelvin_from(above: f64, zero: f64) -> f64 {
    decimal::add(above, zero)
}

/// An absolute temperature in kelvin as degrees Celsius (in decimal: a typed value converts
/// exactly).
pub fn kelvin_to_celsius(kelvin: f64) -> f64 {
    decimal::add(kelvin, -unit::CELSIUS_ZERO)
}

/// An absolute temperature in degrees Celsius as kelvin.
pub fn celsius_to_kelvin(celsius: f64) -> f64 {
    decimal::add(celsius, unit::CELSIUS_ZERO)
}

/// A warning when `si`, a value of quantity `q` called `what`, lies outside `lo..=hi` (SI): a
/// magnitude no engineering value has, which usually means a number typed in the wrong unit.
/// The values are shown in `display`.
pub fn check_magnitude(what: &str, q: Quantity, si: f64, (lo, hi): (f64, f64), display: &DisplaySystem) -> Option<String> {
    if (lo..=hi).contains(&si.abs()) {
        return None;
    }
    let show = |v: f64| q.format(v, display).unwrap_or_else(|_| format!("{v:e} {}", q.si_symbol()));
    Some(format!("{what} {} is outside the plausible {} to {}; check the unit it was given in", show(si), show(lo), show(hi)))
}

/// Plausible moduli: 10 kPa (soft foams) to 2 TPa (diamond, graphene), in Pa.
pub const MODULUS_RANGE: (f64, f64) = (1e4, 2e12);
/// Plausible densities, kg/m³.
pub const DENSITY_RANGE: (f64, f64) = (1.0, 25_000.0);
/// Plausible gravitational accelerations, m/s².
pub const GRAVITY_RANGE: (f64, f64) = (0.1, 100.0);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn magnitudes_are_checked_in_si_and_shown_in_the_display_system() {
        let mm = DisplaySystem::default();
        assert_eq!(check_magnitude("E", Quantity::Pressure, 2.1e11, MODULUS_RANGE, &mm), None);
        // 210000 typed as Pa where MPa was meant.
        let w = check_magnitude("E", Quantity::Pressure, 2.1e5 * 1e6 * 1e6, MODULUS_RANGE, &mm).unwrap();
        assert_eq!(w, "E 2.1e11 MPa is outside the plausible 0.01 MPa to 2000000 MPa; check the unit it was given in");
        assert!(check_magnitude("density", Quantity::Density, 7.85e-9, DENSITY_RANGE, &mm).is_some());
        assert!(check_magnitude("g", Quantity::Acceleration, f64::NAN, GRAVITY_RANGE, &mm).is_some());
    }
}
