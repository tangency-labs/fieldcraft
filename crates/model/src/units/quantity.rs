//! The kinds of physical quantity FieldCraft stores, and parsing and formatting of values with
//! units ("210 GPa", `{"value": 20, "unit": "°C"}`) against them.

use serde::{Deserialize, Serialize};

use super::unit::MAX_TEXT;
use super::{Dimension, DisplaySystem, Unit, UnitError};

macro_rules! quantities {
    ($($variant:ident = $name:literal, $dim:expr, $si:literal;)*) => {
        /// A kind of physical quantity. Each has one dimension and one SI unit, which is how it is
        /// stored.
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
        #[serde(rename_all = "snake_case")]
        pub enum Quantity {
            $($variant,)*
        }

        impl Quantity {
            pub const ALL: &[Quantity] = &[$(Quantity::$variant,)*];

            /// The name used in commands and files (`thermal_conductivity`, …).
            pub fn name(self) -> &'static str {
                match self { $(Quantity::$variant => $name,)* }
            }

            pub fn dimension(self) -> Dimension {
                match self { $(Quantity::$variant => $dim,)* }
            }

            /// The SI unit the quantity is stored in.
            pub fn si_symbol(self) -> &'static str {
                match self { $(Quantity::$variant => $si,)* }
            }
        }
    };
}

const fn d(length: i8, mass: i8, time: i8, temperature: i8) -> Dimension {
    Dimension::new(length, mass, time, temperature)
}

quantities! {
    Dimensionless = "dimensionless", d(0, 0, 0, 0), "1";
    Angle = "angle", d(0, 0, 0, 0), "rad";
    Length = "length", d(1, 0, 0, 0), "m";
    Area = "area", d(2, 0, 0, 0), "m²";
    Volume = "volume", d(3, 0, 0, 0), "m³";
    SecondMomentOfArea = "second_moment_of_area", d(4, 0, 0, 0), "m⁴";
    InverseLength = "inverse_length", d(-1, 0, 0, 0), "1/m";
    Mass = "mass", d(0, 1, 0, 0), "kg";
    RotaryInertia = "rotary_inertia", d(2, 1, 0, 0), "kg·m²";
    Density = "density", d(-3, 1, 0, 0), "kg/m³";
    Time = "time", d(0, 0, 1, 0), "s";
    Frequency = "frequency", d(0, 0, -1, 0), "Hz";
    AngularVelocity = "angular_velocity", d(0, 0, -1, 0), "rad/s";
    StrainRate = "strain_rate", d(0, 0, -1, 0), "1/s";
    Eigenvalue = "eigenvalue", d(0, 0, -2, 0), "1/s²";
    Velocity = "velocity", d(1, 0, -1, 0), "m/s";
    InverseVelocity = "inverse_velocity", d(-1, 0, 1, 0), "s/m";
    Acceleration = "acceleration", d(1, 0, -2, 0), "m/s²";
    VolumeFlow = "volume_flow", d(3, 0, -1, 0), "m³/s";
    MassFlow = "mass_flow", d(0, 1, -1, 0), "kg/s";
    Momentum = "momentum", d(1, 1, -1, 0), "N·s";
    Force = "force", d(1, 1, -2, 0), "N";
    Moment = "moment", d(2, 1, -2, 0), "N·m";
    Pressure = "pressure", d(-1, 1, -2, 0), "Pa";
    ForcePerVolume = "force_per_volume", d(-2, 1, -2, 0), "N/m³";
    Stiffness = "stiffness", d(0, 1, -2, 0), "N/m";
    ForcePerLength = "force_per_length", d(0, 1, -2, 0), "N/m";
    RotationalStiffness = "rotational_stiffness", d(2, 1, -2, 0), "N·m/rad";
    Damping = "damping", d(0, 1, -1, 0), "N·s/m";
    RotationalDamping = "rotational_damping", d(2, 1, -1, 0), "N·m·s/rad";
    Energy = "energy", d(2, 1, -2, 0), "J";
    Power = "power", d(2, 1, -3, 0), "W";
    DynamicViscosity = "dynamic_viscosity", d(-1, 1, -1, 0), "Pa·s";
    KinematicViscosity = "kinematic_viscosity", d(2, 0, -1, 0), "m²/s";
    Temperature = "temperature", d(0, 0, 0, 1), "K";
    TemperatureDifference = "temperature_difference", d(0, 0, 0, 1), "K";
    ThermalExpansion = "thermal_expansion", d(0, 0, 0, -1), "K⁻¹";
    ThermalConductivity = "thermal_conductivity", d(1, 1, -3, -1), "W/(m·K)";
    SpecificHeat = "specific_heat", d(2, 0, -2, -1), "J/(kg·K)";
    HeatTransferCoefficient = "heat_transfer_coefficient", d(0, 1, -3, -1), "W/(m²·K)";
    HeatFlux = "heat_flux", d(0, 1, -3, 0), "W/m²";
    HeatGeneration = "heat_generation", d(-1, 1, -3, 0), "W/m³";
}

impl Quantity {
    pub fn from_name(name: &str) -> Option<Quantity> {
        Quantity::ALL.iter().copied().find(|q| q.name() == name)
    }

    /// Whether the unit's zero offset applies (°C and °F as absolute temperatures).
    fn absolute(self) -> bool {
        self == Quantity::Temperature
    }

    /// Converts `value` in `unit` to this quantity's SI value, checking the dimension and the
    /// result.
    pub fn to_si(self, value: f64, unit: &Unit) -> Result<f64, UnitError> {
        if !value.is_finite() {
            return Err(UnitError::NotFinite(value.to_string()));
        }
        self.check(unit)?;
        let si = unit.to_si(value, self.absolute());
        if !si.is_finite() {
            return Err(UnitError::OutOfRange(format!("{value} {}", unit.symbol())));
        }
        if self.absolute() && si < 0.0 {
            return Err(UnitError::BelowAbsoluteZero(format!("{value} {}", unit.symbol())));
        }
        Ok(si)
    }

    /// Converts an SI value of this quantity to `unit`.
    pub fn express(self, si: f64, unit: &Unit) -> Result<f64, UnitError> {
        self.check(unit)?;
        let v = unit.express(si, self.absolute());
        if v.is_finite() { Ok(v) } else { Err(UnitError::NotFinite(v.to_string())) }
    }

    /// Whether the quantity is per radian or in radians (an angle, an angular velocity, a
    /// torsional spring): its units name the angle (`rad`, `deg`, `°`, `rpm`), and other
    /// quantities' units never do.
    pub fn is_angular(self) -> bool {
        matches!(self, Quantity::Angle | Quantity::AngularVelocity | Quantity::RotationalStiffness | Quantity::RotationalDamping)
    }

    fn check(self, unit: &Unit) -> Result<(), UnitError> {
        if unit.dimension() == self.dimension() {
            // Angles and ratios share the empty dimension (and an angular velocity a frequency's):
            // only angular quantities take rad, deg or °, and they must.
            if unit.is_angle() == self.is_angular() {
                Ok(())
            } else {
                Err(UnitError::WrongKind { unit: unit.symbol().to_string(), quantity: self.name() })
            }
        } else {
            Err(UnitError::WrongDimension {
                unit: unit.symbol().to_string(),
                found: unit.dimension(),
                quantity: self.name(),
                expected: self.dimension(),
            })
        }
    }

    /// Parses `"210 GPa"`, `"25 mm/s"`, `"-5 °C"` (or a bare `"210"`, meaning `display`'s unit)
    /// to the SI value.
    pub fn parse(self, text: &str, display: &DisplaySystem) -> Result<f64, UnitError> {
        let text = text.trim();
        if text.is_empty() {
            return Err(UnitError::Empty);
        }
        if text.len() > MAX_TEXT {
            return Err(UnitError::TooLong);
        }
        let (number, unit) = split_number(text).ok_or_else(|| UnitError::Syntax(text.to_string()))?;
        if unit.trim().is_empty() { self.to_si(number, &display.unit(self)?) } else { self.to_si(number, &Unit::parse(unit)?) }
    }

    /// Formats an SI value in `display`'s unit for this quantity, with the unit: `210000 MPa`.
    pub fn format(self, si: f64, display: &DisplaySystem) -> Result<String, UnitError> {
        let unit = display.unit(self)?;
        let v = self.express(si, &unit)?;
        Ok(if unit.symbol() == "1" { format_number(v) } else { format!("{} {}", format_number(v), unit.symbol()) })
    }
}

/// The leading number of `text` and the rest. `inf` and `nan` are not numbers here.
fn split_number(text: &str) -> Option<(f64, &str)> {
    let end = text.find(|c: char| !(c.is_ascii_digit() || matches!(c, '+' | '-' | '.' | 'e' | 'E'))).unwrap_or(text.len());
    // Shrink until the prefix parses: "2e" in "2em" is not an exponent.
    (1..=end).rev().find_map(|i| {
        let head = text.get(..i)?;
        let v: f64 = head.parse().ok()?;
        Some((v, text.get(i..)?))
    })
}

/// A value with six significant digits and no trailing zeros: `210000`, `7.85e-9`, `0.3`.
pub fn format_number(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return if v == 0.0 { "0".to_string() } else { v.to_string() };
    }
    let exp = v.abs().log10().floor();
    if (-4.0..7.0).contains(&exp) {
        let decimals = (5.0 - exp).max(0.0) as usize;
        let s = format!("{v:.decimals$}");
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s }
    } else {
        let s = format!("{v:.5e}");
        match s.split_once('e') {
            Some((mantissa, e)) => {
                let m = if mantissa.contains('.') { mantissa.trim_end_matches('0').trim_end_matches('.') } else { mantissa };
                format!("{m}e{e}")
            }
            None => s,
        }
    }
}

/// A value as commands, scripts and the MCP give it: a bare number (in the display system), a
/// string with its unit (`"210 GPa"`), or `{"value": 210, "unit": "GPa"}`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ValueInput {
    Number(f64),
    Text(String),
    WithUnit(WithUnit),
}

/// `{"value": 210, "unit": "GPa"}`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WithUnit {
    pub value: f64,
    pub unit: String,
}

impl ValueInput {
    /// The SI value of this input as quantity `q`; a bare number is in `display`'s unit.
    pub fn to_si(&self, q: Quantity, display: &DisplaySystem) -> Result<f64, UnitError> {
        match self {
            ValueInput::Number(v) => q.to_si(*v, &display.unit(q)?),
            ValueInput::Text(t) => q.parse(t, display),
            ValueInput::WithUnit(w) => q.to_si(w.value, &Unit::parse(&w.unit)?),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_si_symbol_parses_to_its_dimension_with_factor_one() {
        for &q in Quantity::ALL {
            let u = Unit::parse(q.si_symbol()).unwrap_or_else(|e| panic!("{q:?}: {e}"));
            assert_eq!(u.dimension(), q.dimension(), "{q:?}");
            assert!((u.factor() - 1.0).abs() < 1e-15, "{q:?}");
            assert_eq!(Quantity::from_name(q.name()), Some(q));
        }
    }

    #[test]
    fn numbers_split_from_units() {
        assert_eq!(split_number("210 GPa"), Some((210.0, " GPa")));
        assert_eq!(split_number("2.1e11Pa"), Some((2.1e11, "Pa")));
        assert_eq!(split_number("-5 °C"), Some((-5.0, " °C")));
        assert_eq!(split_number("3e"), Some((3.0, "e")));
        assert_eq!(split_number("inf Pa"), None);
        assert_eq!(split_number("NaN"), None);
        assert_eq!(split_number("GPa"), None);
    }

    #[test]
    fn numbers_format_with_six_significant_digits() {
        assert_eq!(format_number(210000.0), "210000");
        assert_eq!(format_number(7.85e-9), "7.85e-9");
        assert_eq!(format_number(0.3), "0.3");
        assert_eq!(format_number(2.1e11), "2.1e11");
        assert_eq!(format_number(-20.0), "-20");
        assert_eq!(format_number(1.0 / 3.0), "0.333333");
        assert_eq!(format_number(0.0), "0");
    }
}
