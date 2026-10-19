//! Display systems: which unit each quantity is shown and typed in. Choosing one never changes a
//! stored value.

use serde::{Deserialize, Serialize};

use super::{Quantity, Unit, UnitError};

/// A consistent unit system to show and type values in: every unit in it is a product of its
/// base units (length, mass, time, temperature step), so force = mass × acceleration and
/// stress = force / area hold with no factor, and a deck written in it needs no conversion
/// constants. These are the systems FieldCraft models carried as a label before storage moved
/// to SI (a model of format 2 or older is shown in the system it declared). Angles are in
/// radians in every one, as rotations were before.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Preset {
    /// Metre, kilogram, second, newton, pascal, kelvin: strict SI.
    #[serde(rename = "m-kg-s-N-Pa")]
    MKgS,
    /// Millimetre, tonne, second, newton, megapascal, °C: what most structural CAE users expect,
    /// and the default.
    #[default]
    #[serde(rename = "mm-t-s-N-MPa")]
    MmTonneS,
    /// Millimetre, kilogram, millisecond, kilonewton, gigapascal, °C.
    #[serde(rename = "mm-kg-ms-kN-GPa")]
    MmKgMs,
    /// Millimetre, gram, millisecond, newton, megapascal, °C.
    #[serde(rename = "mm-g-ms-N-MPa")]
    MmGMs,
    /// Inch, lbf·s²/in, second, pound-force, psi, °F.
    #[serde(rename = "in-lbf-s-psi")]
    InLbfS,
}

/// One quantity's unit in each preset, in [`Preset::ALL`] order.
type Row = (Quantity, [&'static str; 5]);

/// The presets' units. `presets_are_consistent` checks every entry against the preset's base
/// units, so a typo here fails the tests rather than a deck.
const TABLE: &[Row] = {
    use Quantity as Q;
    &[
        (Q::Dimensionless, ["1", "1", "1", "1", "1"]),
        (Q::Angle, ["rad", "rad", "rad", "rad", "rad"]),
        (Q::Length, ["m", "mm", "mm", "mm", "in"]),
        (Q::Area, ["m²", "mm²", "mm²", "mm²", "in²"]),
        (Q::Volume, ["m³", "mm³", "mm³", "mm³", "in³"]),
        (Q::SecondMomentOfArea, ["m⁴", "mm⁴", "mm⁴", "mm⁴", "in⁴"]),
        (Q::InverseLength, ["1/m", "1/mm", "1/mm", "1/mm", "1/in"]),
        (Q::Mass, ["kg", "t", "kg", "g", "lbf·s²/in"]),
        (Q::RotaryInertia, ["kg·m²", "t·mm²", "kg·mm²", "g·mm²", "lbf·s²·in"]),
        (Q::Density, ["kg/m³", "t/mm³", "kg/mm³", "g/mm³", "lbf·s²/in⁴"]),
        (Q::Time, ["s", "s", "ms", "ms", "s"]),
        (Q::Frequency, ["Hz", "Hz", "kHz", "kHz", "Hz"]),
        (Q::AngularVelocity, ["rad/s", "rad/s", "rad/ms", "rad/ms", "rad/s"]),
        (Q::StrainRate, ["1/s", "1/s", "1/ms", "1/ms", "1/s"]),
        (Q::Eigenvalue, ["1/s²", "1/s²", "1/ms²", "1/ms²", "1/s²"]),
        (Q::Velocity, ["m/s", "mm/s", "mm/ms", "mm/ms", "in/s"]),
        (Q::InverseVelocity, ["s/m", "s/mm", "ms/mm", "ms/mm", "s/in"]),
        (Q::Acceleration, ["m/s²", "mm/s²", "mm/ms²", "mm/ms²", "in/s²"]),
        (Q::VolumeFlow, ["m³/s", "mm³/s", "mm³/ms", "mm³/ms", "in³/s"]),
        (Q::MassFlow, ["kg/s", "t/s", "kg/ms", "g/ms", "lbf·s/in"]),
        (Q::Momentum, ["N·s", "N·s", "kN·ms", "N·ms", "lbf·s"]),
        (Q::Force, ["N", "N", "kN", "N", "lbf"]),
        (Q::Moment, ["N·m", "N·mm", "kN·mm", "N·mm", "lbf·in"]),
        (Q::Pressure, ["Pa", "MPa", "GPa", "MPa", "psi"]),
        (Q::ForcePerVolume, ["N/m³", "N/mm³", "kN/mm³", "N/mm³", "lbf/in³"]),
        (Q::Stiffness, ["N/m", "N/mm", "kN/mm", "N/mm", "lbf/in"]),
        (Q::ForcePerLength, ["N/m", "N/mm", "kN/mm", "N/mm", "lbf/in"]),
        (Q::RotationalStiffness, ["N·m/rad", "N·mm/rad", "kN·mm/rad", "N·mm/rad", "lbf·in/rad"]),
        (Q::Damping, ["N·s/m", "N·s/mm", "kN·ms/mm", "N·ms/mm", "lbf·s/in"]),
        (Q::RotationalDamping, ["N·m·s/rad", "N·mm·s/rad", "kN·mm·ms/rad", "N·mm·ms/rad", "lbf·in·s/rad"]),
        (Q::Energy, ["J", "mJ", "J", "mJ", "in·lbf"]),
        (Q::Power, ["W", "mW", "kW", "W", "in·lbf/s"]),
        (Q::DynamicViscosity, ["Pa·s", "MPa·s", "GPa·ms", "MPa·ms", "psi·s"]),
        (Q::KinematicViscosity, ["m²/s", "mm²/s", "mm²/ms", "mm²/ms", "in²/s"]),
        (Q::Temperature, ["K", "°C", "°C", "°C", "°F"]),
        (Q::TemperatureDifference, ["K", "K", "K", "K", "°F"]),
        (Q::ThermalExpansion, ["K⁻¹", "K⁻¹", "K⁻¹", "K⁻¹", "°F⁻¹"]),
        (Q::ThermalConductivity, ["W/(m·K)", "mW/(mm·K)", "kW/(mm·K)", "W/(mm·K)", "lbf/(s·°F)"]),
        (Q::SpecificHeat, ["J/(kg·K)", "mJ/(t·K)", "J/(kg·K)", "mJ/(g·K)", "in²/(s²·°F)"]),
        (Q::HeatTransferCoefficient, ["W/(m²·K)", "mW/(mm²·K)", "kW/(mm²·K)", "W/(mm²·K)", "lbf/(in·s·°F)"]),
        (Q::HeatFlux, ["W/m²", "mW/mm²", "kW/mm²", "W/mm²", "lbf/(in·s)"]),
        (Q::HeatGeneration, ["W/m³", "mW/mm³", "kW/mm³", "W/mm³", "lbf/(in²·s)"]),
    ]
};

impl Preset {
    pub const ALL: [Preset; 5] = [Preset::MKgS, Preset::MmTonneS, Preset::MmKgMs, Preset::MmGMs, Preset::InLbfS];

    /// The label used in commands and files.
    pub fn label(self) -> &'static str {
        match self {
            Preset::MKgS => "m-kg-s-N-Pa",
            Preset::MmTonneS => "mm-t-s-N-MPa",
            Preset::MmKgMs => "mm-kg-ms-kN-GPa",
            Preset::MmGMs => "mm-g-ms-N-MPa",
            Preset::InLbfS => "in-lbf-s-psi",
        }
    }

    /// The preset of a label. `SI` is the metre–kilogram–second one; `unspecified`, the label of
    /// models that never chose a system, reads as millimetre–tonne–second (their geometry is in
    /// millimetres).
    pub fn from_label(label: &str) -> Option<Preset> {
        match label {
            "SI" => Some(Preset::MKgS),
            "unspecified" => Some(Preset::MmTonneS),
            _ => Preset::ALL.into_iter().find(|p| p.label() == label),
        }
    }

    fn column(self) -> usize {
        Preset::ALL.iter().position(|p| *p == self).unwrap_or(0)
    }

    /// The unit this preset shows `q` in.
    pub fn symbol(self, q: Quantity) -> &'static str {
        TABLE.iter().find(|(r, _)| *r == q).and_then(|(_, units)| units.get(self.column()).copied()).unwrap_or_else(|| q.si_symbol())
    }
}

/// Most per-quantity overrides a display system holds (one per quantity).
const MAX_OVERRIDES: usize = 64;

/// The units values are shown and typed in: a preset with optional per-quantity overrides (say
/// mm-t-s-N-MPa with energies in J and angles in degrees).
#[derive(Debug, Clone, PartialEq)]
pub struct DisplaySystem {
    preset: Preset,
    /// Sorted by quantity, one per quantity.
    overrides: Vec<(Quantity, Unit)>,
    /// Every quantity's unit, in [`Quantity::ALL`] order, parsed once.
    units: Vec<Unit>,
}

impl Default for DisplaySystem {
    fn default() -> Self {
        DisplaySystem::new(Preset::default())
    }
}

/// The quantities whose unit holds the temperature's step: they follow a temperature override.
const PER_DEGREE: [Quantity; 5] = [
    Quantity::TemperatureDifference,
    Quantity::ThermalExpansion,
    Quantity::ThermalConductivity,
    Quantity::SpecificHeat,
    Quantity::HeatTransferCoefficient,
];

/// The degree of a temperature unit as the per-degree units write it: `K`, `°C` or `°F`
/// (`None` for another scale).
fn degree_symbol(u: &Unit) -> Option<&'static str> {
    let t = Quantity::Temperature;
    let (zero, step) = (t.to_si(0.0, u).ok()?, t.to_si(1.0, u).ok()? - t.to_si(0.0, u).ok()?);
    [("K", 0.0, 1.0), ("°C", 273.15, 1.0), ("°F", 255.3722222222222, 5.0 / 9.0)]
        .into_iter()
        .find(|(_, z, s)| (zero - z).abs() < 1e-9 && (step - s).abs() < 1e-12)
        .map(|(d, ..)| d)
}

impl DisplaySystem {
    pub fn new(preset: Preset) -> DisplaySystem {
        let mut d = DisplaySystem { preset, overrides: Vec::new(), units: Vec::new() };
        d.build();
        d
    }

    /// The display system of a preset label (see [`Preset::from_label`]).
    pub fn from_label(label: &str) -> Option<DisplaySystem> {
        Preset::from_label(label).map(DisplaySystem::new)
    }

    fn build(&mut self) {
        // A temperature shown in another scale carries to the temperature difference and the
        // per-degree quantities, unless they have overrides of their own: with mm-t-s and °F,
        // ΔT is in °F and an expansion coefficient per °F.
        let degree = self.overrides.iter().find(|(o, _)| *o == Quantity::Temperature).and_then(|(_, u)| degree_symbol(u));
        let step = self.preset.symbol(Quantity::TemperatureDifference);
        self.units = Quantity::ALL
            .iter()
            .map(|&q| match self.overrides.iter().find(|(o, _)| *o == q) {
                Some((_, u)) => u.clone(),
                None => {
                    let symbol = match degree {
                        Some(d) if PER_DEGREE.contains(&q) && d != step => self.preset.symbol(q).replace(step, d),
                        _ => self.preset.symbol(q).to_string(),
                    };
                    // Every table entry parses (`every_preset_unit_has_its_quantitys_dimension`).
                    Unit::parse(&symbol)
                        .or_else(|_| Unit::parse(self.preset.symbol(q)))
                        .or_else(|_| Unit::parse(q.si_symbol()))
                        .unwrap_or_else(|_| Unit::one())
                }
            })
            .collect();
    }

    pub fn preset(&self) -> Preset {
        self.preset
    }

    /// The overrides, by quantity.
    pub fn overrides(&self) -> impl Iterator<Item = (Quantity, &str)> {
        self.overrides.iter().map(|(q, u)| (*q, u.symbol()))
    }

    /// Shows `q` in `unit` instead of the preset's unit; refused when the dimension differs.
    pub fn with_override(mut self, q: Quantity, unit: &str) -> Result<DisplaySystem, UnitError> {
        let unit = Unit::parse(unit)?;
        q.express(0.0, &unit)?;
        self.overrides.retain(|(o, _)| *o != q);
        if unit.symbol() != self.preset.symbol(q) {
            self.overrides.push((q, unit));
            self.overrides.sort_by_key(|(o, _)| *o);
        }
        self.build();
        Ok(self)
    }

    /// The unit `q` is shown and typed in.
    pub fn unit(&self, q: Quantity) -> Result<Unit, UnitError> {
        Ok(self.unit_ref(q).clone())
    }

    /// The unit `q` is shown and typed in, borrowed.
    pub fn unit_ref(&self, q: Quantity) -> &Unit {
        static ONE: std::sync::OnceLock<Unit> = std::sync::OnceLock::new();
        Quantity::ALL.iter().position(|o| *o == q).and_then(|i| self.units.get(i)).unwrap_or_else(|| ONE.get_or_init(Unit::one))
    }

    /// The symbol of the unit `q` is shown in (`MPa`).
    pub fn symbol(&self, q: Quantity) -> &str {
        self.unit_ref(q).symbol()
    }

    /// The SI value of `value` given in this system's unit for `q`.
    pub fn to_si(&self, q: Quantity, value: f64) -> Result<f64, UnitError> {
        q.to_si(value, self.unit_ref(q))
    }

    /// An SI value of `q` in this system's unit.
    pub fn express(&self, q: Quantity, si: f64) -> Result<f64, UnitError> {
        q.express(si, self.unit_ref(q))
    }

    /// The preset's label, with the overrides: `mm-t-s-N-MPa (energy: J)`.
    pub fn label(&self) -> String {
        if self.overrides.is_empty() {
            return self.preset.label().to_string();
        }
        let list: Vec<String> = self.overrides.iter().map(|(q, u)| format!("{}: {}", q.name(), u.symbol())).collect();
        format!("{} ({})", self.preset.label(), list.join(", "))
    }
}

/// A display system as stored: `"mm-t-s-N-MPa"`, or `{"preset": …, "overrides": {"energy": "J"}}`.
#[derive(Serialize, Deserialize)]
#[serde(untagged)]
enum Repr {
    Preset(Preset),
    Full(Full),
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Full {
    preset: Preset,
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    overrides: std::collections::BTreeMap<String, String>,
}

impl Serialize for DisplaySystem {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        if self.overrides.is_empty() {
            return Repr::Preset(self.preset).serialize(s);
        }
        let overrides = self.overrides.iter().map(|(q, u)| (q.name().to_string(), u.symbol().to_string())).collect();
        Repr::Full(Full { preset: self.preset, overrides }).serialize(s)
    }
}

impl<'de> Deserialize<'de> for DisplaySystem {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        use serde::de::Error;
        match Repr::deserialize(d)? {
            Repr::Preset(p) => Ok(DisplaySystem::new(p)),
            Repr::Full(Full { preset, overrides }) => {
                if overrides.len() > MAX_OVERRIDES {
                    return Err(D::Error::custom(format!("a display system overrides at most {MAX_OVERRIDES} quantities")));
                }
                let mut out = DisplaySystem::new(preset);
                for (name, unit) in overrides {
                    let q = Quantity::from_name(&name)
                        .ok_or_else(|| D::Error::custom(format!("unknown quantity `{}`", name.chars().take(64).collect::<String>())))?;
                    out = out.with_override(q, &unit).map_err(|e| D::Error::custom(format!("{name}: {e}")))?;
                }
                Ok(out)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_quantity_is_in_the_table_once() {
        for &q in Quantity::ALL {
            assert_eq!(TABLE.iter().filter(|(r, _)| *r == q).count(), 1, "{q:?}");
        }
    }

    /// Every preset unit is the product of the preset's base units raised to the quantity's
    /// dimension: the system is consistent.
    #[test]
    fn presets_are_consistent() {
        for preset in Preset::ALL {
            let base = |q: Quantity| Unit::parse(preset.symbol(q)).unwrap().factor();
            let (l, m, t, th) = (base(Quantity::Length), base(Quantity::Mass), base(Quantity::Time), base(Quantity::TemperatureDifference));
            for &q in Quantity::ALL {
                let d = q.dimension();
                let want = l.powi(d.length.into()) * m.powi(d.mass.into()) * t.powi(d.time.into()) * th.powi(d.temperature.into());
                let got = Unit::parse(preset.symbol(q)).unwrap().factor();
                assert!((got - want).abs() <= 1e-12 * want, "{preset:?} {q:?}: {} is {got}, the base units give {want}", preset.symbol(q));
            }
        }
    }

    #[test]
    fn stored_forms() {
        let d = DisplaySystem::new(Preset::MmTonneS);
        assert_eq!(serde_json::to_value(&d).unwrap(), serde_json::json!("mm-t-s-N-MPa"));
        let o = d.clone().with_override(Quantity::Energy, "J").unwrap();
        let v = serde_json::to_value(&o).unwrap();
        assert_eq!(v, serde_json::json!({"preset": "mm-t-s-N-MPa", "overrides": {"energy": "J"}}));
        assert_eq!(serde_json::from_value::<DisplaySystem>(v).unwrap(), o);
        // An override back to the preset's own unit is no override.
        assert_eq!(o.with_override(Quantity::Energy, "mJ").unwrap(), d);
        for bad in [
            serde_json::json!("mm-t-s"),
            serde_json::json!({"preset": "SI"}),
            serde_json::json!({"preset": "mm-t-s-N-MPa", "overrides": {"energy": "MPa"}}),
            serde_json::json!({"preset": "mm-t-s-N-MPa", "overrides": {"warp": "m"}}),
            serde_json::json!({"preset": "mm-t-s-N-MPa", "extra": 1}),
        ] {
            assert!(serde_json::from_value::<DisplaySystem>(bad.clone()).is_err(), "{bad}");
        }
        assert_eq!(DisplaySystem::from_label("unspecified"), Some(DisplaySystem::new(Preset::MmTonneS)));
        assert_eq!(DisplaySystem::from_label("SI"), Some(DisplaySystem::new(Preset::MKgS)));
    }
}
