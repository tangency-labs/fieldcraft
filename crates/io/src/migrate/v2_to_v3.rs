//! Format 2 → 3: every stored quantity in SI. **Released with format 3: never edit.** A mistake
//! is fixed by a new step.
//!
//! Formats 1 and 2 stored numbers as they were typed, in the unit system the model declared in
//! `model.json`'s `units` (`mm-t-s-N-MPa`, `m-kg-s-N-Pa`, `mm-kg-ms-kN-GPa`, `mm-g-ms-N-MPa`,
//! `in-lbf-s-psi`, or `unspecified` when none was chosen). Format 3 stores SI (m, kg, s, N, Pa,
//! K, rad) and shows the model in a display system. This step:
//!
//! - **`model.json`**: scales every dimensional number to SI by the declared system, through its
//!   own field map ([`Walk::model`]) and its own factors ([`System`]), which stand for format 2
//!   as released and never follow later changes to the model's types or the unit table. A
//!   number the map does not know refuses the file rather than pass unscaled. `units` becomes
//!   `display_units`, the declared system, so every number a user sees stays the same.
//!   - `unspecified` is read as mm-t-s-N-MPa: a model's geometry was always in millimetres, and
//!     every default and example was in that system.
//!   - Absolute temperatures were in °C in the millimetre systems, K in m-kg-s-N-Pa and °F in
//!     in-lbf-s-psi (the systems' own temperature units; every solver result depended only on
//!     temperature differences, so nothing computed changes). They go to kelvin with the offset;
//!     temperature differences, expansion coefficients and every other per-degree quantity
//!     scale without it.
//!   - Lengths that the commands always took in millimetres, whatever the system (datum points,
//!     mesh plans: geometry from SolidCraft is in millimetres), scale from millimetres.
//!   - Defaults the model read for absent dimensional fields (a reference temperature of 0, a
//!     plane strain thickness and a Johnson–Cook reference strain rate of 1) are written out,
//!     converted, since a default means the same number in SI no longer.
//! - **`mesh.bin`**: `MCMESH01` → `MCMESH02`, the same layout with node coordinates scaled from
//!   the declared length unit to metres ([`mesh`]; its own reader of the `MCMESH01` layout).
//! - **`journal.json`**: its numbers are not touched. The list of commands gets a header naming
//!   the system its bare numbers were typed in (`units`, the declared system) and the system a
//!   model started in before the first command (`model_units`: `unspecified`, read as
//!   mm-t-s-N-MPa), so replaying it gives the model this step makes:
//!   `{"units": "mm-t-s-N-MPa", "model_units": "mm-t-s-N-MPa", "commands": [...]}`.
//!
//! Scaling by a power of ten works on the decimal a number was written as ([`decimal`]), as the
//! unit table of the release that made format 3 did, so a typed 3.97 mm is the double nearest
//! 0.00397 m.

use std::collections::BTreeSet;

use serde_json::{Map, Value};

use super::{Parts, set_format};
use crate::{Error, Result, invalid};

/// What a format 2 number is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Q {
    Length,
    /// A length the commands took in millimetres whatever the system.
    Millimetre,
    Area,
    InverseLength,
    Mass,
    RotaryInertia,
    Density,
    Time,
    Frequency,
    AngularVelocity,
    StrainRate,
    Velocity,
    InverseVelocity,
    Acceleration,
    VolumeFlow,
    Force,
    Moment,
    Pressure,
    ForcePerVolume,
    Stiffness,
    RotationalStiffness,
    Damping,
    RotationalDamping,
    DynamicViscosity,
    Temperature,
    ThermalExpansion,
    ThermalConductivity,
    SpecificHeat,
    HeatTransferCoefficient,
    HeatFlux,
    HeatGeneration,
}

/// The exponents of length, mass, time and temperature of a quantity.
fn dimension(q: Q) -> [i32; 4] {
    match q {
        Q::Length | Q::Millimetre => [1, 0, 0, 0],
        Q::Area => [2, 0, 0, 0],
        Q::InverseLength => [-1, 0, 0, 0],
        Q::Mass => [0, 1, 0, 0],
        Q::RotaryInertia => [2, 1, 0, 0],
        Q::Density => [-3, 1, 0, 0],
        Q::Time => [0, 0, 1, 0],
        Q::Frequency | Q::AngularVelocity | Q::StrainRate => [0, 0, -1, 0],
        Q::Velocity => [1, 0, -1, 0],
        Q::InverseVelocity => [-1, 0, 1, 0],
        Q::Acceleration => [1, 0, -2, 0],
        Q::VolumeFlow => [3, 0, -1, 0],
        Q::Force => [1, 1, -2, 0],
        Q::Moment | Q::RotationalStiffness => [2, 1, -2, 0],
        Q::Pressure => [-1, 1, -2, 0],
        Q::ForcePerVolume => [-2, 1, -2, 0],
        Q::Stiffness => [0, 1, -2, 0],
        Q::Damping => [0, 1, -1, 0],
        Q::RotationalDamping => [2, 1, -1, 0],
        Q::DynamicViscosity => [-1, 1, -1, 0],
        Q::Temperature => [0, 0, 0, 1],
        Q::ThermalExpansion => [0, 0, 0, -1],
        Q::ThermalConductivity => [1, 1, -3, -1],
        Q::SpecificHeat => [2, 0, -2, -1],
        Q::HeatTransferCoefficient => [0, 1, -3, -1],
        Q::HeatFlux => [0, 1, -3, 0],
        Q::HeatGeneration => [-1, 1, -3, 0],
    }
}

/// A unit system of format 2, as a conversion to SI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum System {
    MKgS,
    MmTonneS,
    MmKgMs,
    MmGMs,
    InLbfS,
}

/// A factor to SI, `mantissa × 10^exp10`, and the SI value of the unit's zero.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Factor {
    mantissa: f64,
    exp10: i32,
    offset: f64,
}

/// One degree Fahrenheit in kelvin, and Fahrenheit zero.
const FAHRENHEIT_STEP: f64 = 5.0 / 9.0;
const FAHRENHEIT_ZERO: f64 = 459.67 * FAHRENHEIT_STEP;
const CELSIUS_ZERO: f64 = 273.15;

impl System {
    fn of(label: &str) -> Option<System> {
        Some(match label {
            "m-kg-s-N-Pa" => System::MKgS,
            "mm-t-s-N-MPa" | "unspecified" => System::MmTonneS,
            "mm-kg-ms-kN-GPa" => System::MmKgMs,
            "mm-g-ms-N-MPa" => System::MmGMs,
            "in-lbf-s-psi" => System::InLbfS,
            _ => return None,
        })
    }

    /// The label format 3 shows the model in (`unspecified` becomes mm-t-s-N-MPa).
    fn label(self) -> &'static str {
        match self {
            System::MKgS => "m-kg-s-N-Pa",
            System::MmTonneS => "mm-t-s-N-MPa",
            System::MmKgMs => "mm-kg-ms-kN-GPa",
            System::MmGMs => "mm-g-ms-N-MPa",
            System::InLbfS => "in-lbf-s-psi",
        }
    }

    /// The factor of quantity `q` in this system.
    fn factor(self, q: Q) -> Factor {
        // The metric systems' base units are powers of ten of SI's: (length, mass, time).
        let base = match self {
            System::MKgS => [0, 0, 0],
            System::MmTonneS => [-3, 3, 0],
            System::MmKgMs => [-3, 0, -3],
            System::MmGMs => [-3, -3, -3],
            System::InLbfS => return inch_pound(q),
        };
        let [l, m, t, _] = dimension(q);
        let exp10 = if q == Q::Millimetre { -3 } else { l * base[0] + m * base[1] + t * base[2] };
        // Absolute temperatures: kelvin in m-kg-s-N-Pa, °C in the millimetre systems.
        let offset = if q == Q::Temperature && self != System::MKgS { CELSIUS_ZERO } else { 0.0 };
        Factor { mantissa: 1.0, exp10, offset }
    }
}

/// The in-lbf-s-psi factors (inch, lbf·s²/in, second, °F), as the unit table of the release
/// that made format 3 computed them.
fn inch_pound(q: Q) -> Factor {
    let mantissa = match q {
        Q::Millimetre => return Factor { mantissa: 1.0, exp10: -3, offset: 0.0 },
        Q::Length | Q::Velocity | Q::Acceleration => 0.0254,
        Q::Area => 0.00064516,
        Q::InverseLength | Q::InverseVelocity => 39.37007874015748,
        Q::Mass | Q::Stiffness | Q::Damping | Q::HeatFlux => 175.12683524647636,
        Q::RotaryInertia | Q::Moment | Q::RotationalStiffness | Q::RotationalDamping => 0.11298482902761668,
        Q::Density => 10686895.178201316,
        Q::Time | Q::Frequency | Q::AngularVelocity | Q::StrainRate => 1.0,
        Q::VolumeFlow => 1.6387064e-5,
        Q::Force => 4.4482216152605,
        Q::Pressure | Q::DynamicViscosity | Q::HeatGeneration => 6894.757293168361,
        Q::ForcePerVolume => 271447.1375263134,
        Q::Temperature => FAHRENHEIT_STEP,
        Q::ThermalExpansion => 1.7999999999999998,
        Q::ThermalConductivity => 8.006798907468898,
        Q::SpecificHeat => 0.001161288,
        Q::HeatTransferCoefficient => 315.22830344365747,
    };
    let offset = if q == Q::Temperature { FAHRENHEIT_ZERO } else { 0.0 };
    Factor { mantissa, exp10: 0, offset }
}

impl Factor {
    fn to_si(self, v: f64) -> f64 {
        let scaled = if self.mantissa == 1.0 { v } else { v * self.mantissa };
        let scaled = decimal::shift(scaled, self.exp10);
        if self.offset != 0.0 { decimal::add(scaled, self.offset) } else { scaled }
    }
}

/// Run the step.
pub(super) fn v2_to_v3(parts: &mut Parts) -> Result<()> {
    let Some(model) = parts.model.as_object_mut() else { return invalid("`model.json` is not a JSON object") };
    let label = match model.remove("units") {
        None => "unspecified".to_string(),
        Some(Value::String(s)) => s,
        Some(_) => return invalid("`model.json`: `units` is not a unit system label"),
    };
    let system = System::of(&label)
        .ok_or_else(|| Error::Invalid(format!("`model.json`: unknown unit system `{}`", label.chars().take(64).collect::<String>())))?;
    model.insert("display_units".into(), Value::String(system.label().into()));
    let mut walk = Walk { system, handled: BTreeSet::new() };
    walk.model(&mut parts.model)?;
    walk.check_all_handled(&parts.model)?;
    parts.mesh = mesh(&parts.mesh, system.factor(Q::Length))?;
    let commands = std::mem::take(&mut parts.journal);
    if !commands.is_array() {
        return invalid("`journal.json` is not a list of commands");
    }
    let mut journal = Map::new();
    journal.insert("units".into(), Value::String(system.label().into()));
    journal.insert("model_units".into(), Value::String(System::MmTonneS.label().into()));
    journal.insert("commands".into(), commands);
    parts.journal = Value::Object(journal);
    set_format(parts, 3)
}

/// Walks `model.json`, converting each number it knows and recording the JSON pointer of every
/// number it has dealt with (converted or dimensionless).
struct Walk {
    system: System,
    handled: BTreeSet<String>,
}

/// A JSON pointer token (RFC 6901).
fn token(key: &str) -> String {
    key.replace('~', "~0").replace('/', "~1")
}

impl Walk {
    fn convert(&self, v: f64, q: Q, at: &str) -> Result<Value> {
        let si = self.system.factor(q).to_si(v);
        serde_json::Number::from_f64(si)
            .map(Value::Number)
            .ok_or_else(|| Error::Invalid(format!("`model.json`: the number at {at} is out of range in SI")))
    }

    /// The value at `at`: a number, a list of them (nested) or null, of quantity `q` (`None`:
    /// dimensionless, a count or a label).
    fn value(&mut self, v: &mut Value, at: &str, q: Option<Q>) -> Result<()> {
        match v {
            Value::Number(n) => {
                if let Some(q) = q {
                    let x = n.as_f64().ok_or_else(|| Error::Invalid(format!("`model.json`: the number at {at} is not finite")))?;
                    *v = self.convert(x, q, at)?;
                }
                self.handled.insert(at.to_string());
            }
            Value::Array(items) => {
                for (i, item) in items.iter_mut().enumerate() {
                    self.value(item, &format!("{at}/{i}"), q)?;
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// Field `key` of `parent` (at `at`), when present.
    fn field(&mut self, parent: &mut Value, at: &str, key: &str, q: Option<Q>) -> Result<()> {
        match parent.get_mut(key) {
            Some(v) => self.value(v, &format!("{at}/{}", token(key)), q),
            None => Ok(()),
        }
    }

    /// Field `key` of quantity `q`, written out as `default` (in the declared system) when
    /// absent: a default that meant this number in the declared system.
    fn field_or(&mut self, parent: &mut Value, at: &str, key: &str, q: Q, default: f64) -> Result<()> {
        if let Some(obj) = parent.as_object_mut()
            && !obj.contains_key(key)
            && let Some(n) = serde_json::Number::from_f64(default)
        {
            obj.insert(key.into(), Value::Number(n));
        }
        self.field(parent, at, key, Some(q))
    }

    /// Every number under `parent`'s `key` is dimensionless (a count, a label, a ratio).
    fn plain(&mut self, parent: &mut Value, at: &str, key: &str) {
        if let Some(v) = parent.get_mut(key) {
            self.mark_all(v, &format!("{at}/{}", token(key)));
        }
    }

    fn mark_all(&mut self, v: &Value, at: &str) {
        match v {
            Value::Number(_) => {
                self.handled.insert(at.to_string());
            }
            Value::Array(items) => items.iter().enumerate().for_each(|(i, x)| self.mark_all(x, &format!("{at}/{i}"))),
            Value::Object(m) => m.iter().for_each(|(k, x)| self.mark_all(x, &format!("{at}/{}", token(k)))),
            _ => {}
        }
    }

    /// The members of a map of definitions by name, each walked by `f`.
    fn each_named(&mut self, model: &mut Value, key: &str, mut f: impl FnMut(&mut Self, &mut Value, &str) -> Result<()>) -> Result<()> {
        if let Some(Value::Object(items)) = model.get_mut(key) {
            for (name, item) in items.iter_mut() {
                f(self, item, &format!("/{}/{}", token(key), token(name)))?;
            }
        }
        Ok(())
    }

    /// The items of a list, each walked by `f`.
    fn each(&mut self, parent: &mut Value, at: &str, key: &str, mut f: impl FnMut(&mut Self, &mut Value, &str) -> Result<()>) -> Result<()> {
        if let Some(Value::Array(items)) = parent.get_mut(key) {
            for (i, item) in items.iter_mut().enumerate() {
                f(self, item, &format!("{at}/{}/{i}", token(key)))?;
            }
        }
        Ok(())
    }

    fn model(&mut self, m: &mut Value) -> Result<()> {
        self.each_named(m, "sets", |w, set, at| {
            for key in ["nodes", "elements", "faces", "edges"] {
                w.plain(set, at, key);
            }
            Ok(())
        })?;
        self.each_named(m, "materials", Walk::material)?;
        self.each_named(m, "sections", Walk::section)?;
        self.each_named(m, "reference_points", |w, rp, at| {
            w.field(rp, at, "xyz", Some(Q::Length))?;
            w.plain(rp, at, "label");
            Ok(())
        })?;
        self.each_named(m, "couplings", |w, c, at| {
            w.plain(c, at, "dofs");
            Ok(())
        })?;
        self.each_named(m, "contacts", |w, c, at| {
            w.friction(c, at)?;
            w.plain(c, at, "penalty_scale");
            Ok(())
        })?;
        self.each_named(m, "ties", |w, t, at| w.field(t, at, "position_tolerance", Some(Q::Length)))?;
        self.each_named(m, "masses", |w, p, at| {
            w.field(p, at, "mass", Some(Q::Mass))?;
            w.field(p, at, "inertia", Some(Q::RotaryInertia))
        })?;
        self.each_named(m, "springs", Walk::spring)?;
        self.each_named(m, "boundary_conditions", Walk::boundary_condition)?;
        self.each_named(m, "loads", Walk::load)?;
        self.each_named(m, "fluid_domains", Walk::fluid_domain)?;
        self.each(m, "", "steps", Walk::step)?;
        self.each_named(m, "datums", |w, d, at| {
            w.field(d, at, "point", Some(Q::Millimetre))?;
            w.plain(d, at, "normal");
            w.plain(d, at, "direction");
            Ok(())
        })?;
        self.each(m, "", "seeds", |w, s, at| {
            w.plain(s, at, "count");
            w.field(s, at, "size", Some(Q::Length))?;
            w.plain(s, at, "bias");
            Ok(())
        })?;
        self.each_named(m, "size_controls", |w, c, at| {
            w.field(c, at, "size", Some(Q::Length))?;
            if let Some(region) = c.get_mut("region") {
                let at = format!("{at}/region");
                for key in ["min", "max", "center", "radius"] {
                    w.field(region, &at, key, Some(Q::Length))?;
                }
            }
            Ok(())
        })?;
        self.plain(m, "", "size_growth");
        self.each_named(m, "relaxations", |w, r, at| {
            w.each(r, at, "cells", |w, cell, at| {
                w.field(cell, at, "at", Some(Q::Length))?;
                w.field(cell, at, "moved", Some(Q::Length))
            })
        })?;
        self.each_named(m, "mesh_plans", |w, p, at| {
            w.field(p, at, "size", Some(Q::Millimetre))?;
            w.field(p, at, "tolerance", Some(Q::Millimetre))?;
            for key in ["layers", "direction", "directions", "features", "unmapped"] {
                w.plain(p, at, key);
            }
            w.each(p, at, "cuts", |w, cut, at| {
                w.field(cut, at, "point", Some(Q::Millimetre))?;
                w.plain(cut, at, "normal");
                w.plain(cut, at, "feature");
                Ok(())
            })?;
            w.each(p, at, "rings", |w, ring, at| {
                for key in ["center", "radius", "frame"] {
                    w.field(ring, at, key, Some(Q::Millimetre))?;
                }
                w.plain(ring, at, "axis");
                w.plain(ring, at, "feature");
                Ok(())
            })
        })?;
        self.each(m, "", "stitches", |w, s, at| w.field(s, at, "tolerance", Some(Q::Length)))?;
        self.each(m, "", "defeatures", |w, d, at| {
            w.plain(d, at, "features");
            Ok(())
        })?;
        Ok(())
    }

    fn material(&mut self, m: &mut Value, at: &str) -> Result<()> {
        if let Some(e) = m.get_mut("elastic") {
            let at = format!("{at}/elastic");
            self.field(e, &at, "young", Some(Q::Pressure))?;
            self.field(e, &at, "shear", Some(Q::Pressure))?;
            self.plain(e, &at, "poisson");
        }
        self.field(m, at, "density", Some(Q::Density))?;
        if let Some(f) = m.get_mut("fluid") {
            self.field(f, &format!("{at}/fluid"), "viscosity", Some(Q::DynamicViscosity))?;
        }
        if let Some(h) = m.get_mut("hyperelastic") {
            let at = format!("{at}/hyperelastic");
            for key in ["c01", "c20", "c30"] {
                self.field(h, &at, key, Some(Q::Pressure))?;
            }
        }
        if let Some(p) = m.get_mut("plastic").filter(|p| p.is_object()) {
            let at = format!("{at}/plastic");
            let johnson_cook = p.get("type").and_then(Value::as_str) == Some("johnson_cook");
            for key in ["yield_stress", "hardening_modulus", "a", "b"] {
                self.field(p, &at, key, Some(Q::Pressure))?;
            }
            self.each(p, &at, "points", |w, point, at| {
                w.plain(point, at, "strain");
                w.field(point, at, "stress", Some(Q::Pressure))
            })?;
            for key in ["n", "c", "m", "taylor_quinney"] {
                self.plain(p, &at, key);
            }
            if johnson_cook {
                self.field_or(p, &at, "reference_strain_rate", Q::StrainRate, 1.0)?;
            }
            self.field(p, &at, "room_temperature", Some(Q::Temperature))?;
            self.field(p, &at, "melt_temperature", Some(Q::Temperature))?;
            self.field(p, &at, "specific_heat", Some(Q::SpecificHeat))?;
        }
        if let Some(d) = m.get_mut("damage") {
            self.mark_all(d, &format!("{at}/damage"));
        }
        if let Some(e) = m.get_mut("erosion") {
            self.mark_all(e, &format!("{at}/erosion"));
        }
        if let Some(t) = m.get_mut("thermal").filter(|t| t.is_object()) {
            let at = format!("{at}/thermal");
            self.field(t, &at, "conductivity", Some(Q::ThermalConductivity))?;
            self.field(t, &at, "specific_heat", Some(Q::SpecificHeat))?;
            self.field(t, &at, "expansion", Some(Q::ThermalExpansion))?;
            self.field_or(t, &at, "reference_temperature", Q::Temperature, 0.0)?;
        }
        Ok(())
    }

    fn section(&mut self, s: &mut Value, at: &str) -> Result<()> {
        let Some(kind) = s.get_mut("kind") else { return Ok(()) };
        let at = format!("{at}/kind");
        match kind.get("type").and_then(Value::as_str) {
            Some("plane_strain") => self.field_or(kind, &at, "thickness", Q::Length, 1.0)?,
            _ => self.field(kind, &at, "thickness", Some(Q::Length))?,
        }
        self.plain(kind, &at, "offset");
        self.plain(kind, &at, "points");
        self.plain(kind, &at, "orientation");
        self.field(kind, &at, "area", Some(Q::Area))?;
        if let Some(profile) = kind.get_mut("profile") {
            let at = format!("{at}/profile");
            for key in ["width", "height", "radius", "thickness", "web", "flange"] {
                self.field(profile, &at, key, Some(Q::Length))?;
            }
        }
        Ok(())
    }

    fn friction(&mut self, parent: &mut Value, at: &str) -> Result<()> {
        if let Some(f) = parent.get_mut("friction") {
            let at = format!("{at}/friction");
            self.plain(f, &at, "static");
            self.plain(f, &at, "kinetic");
            self.field(f, &at, "decay", Some(Q::InverseVelocity))?;
        }
        Ok(())
    }

    fn spring(&mut self, s: &mut Value, at: &str) -> Result<()> {
        // A spring along degree of freedom 4 to 6 is torsional.
        let torsional = s.get("direction").and_then(|d| d.get("dof")).and_then(Value::as_u64).is_some_and(|d| (4..=6).contains(&d));
        self.field(s, at, "stiffness", Some(if torsional { Q::RotationalStiffness } else { Q::Stiffness }))?;
        self.field(s, at, "coefficient", Some(if torsional { Q::RotationalDamping } else { Q::Damping }))?;
        if let Some(d) = s.get_mut("direction") {
            let at = format!("{at}/direction");
            self.plain(d, &at, "dof");
            self.plain(d, &at, "vector");
        }
        Ok(())
    }

    fn boundary_condition(&mut self, bc: &mut Value, at: &str) -> Result<()> {
        // Displacements: three lengths, then three rotations in radians.
        if let Some(Value::Array(values)) = bc.get_mut("values") {
            for (i, v) in values.iter_mut().enumerate() {
                self.value(v, &format!("{at}/values/{i}"), (i < 3).then_some(Q::Length))?;
            }
        }
        if let Some(plane) = bc.get_mut("plane") {
            self.mark_all(plane, &format!("{at}/plane"));
        }
        self.field(bc, at, "value", Some(Q::Temperature))
    }

    fn load(&mut self, load: &mut Value, at: &str) -> Result<()> {
        let kind = load.get("type").and_then(Value::as_str).unwrap_or_default().to_string();
        let vector = match kind.as_str() {
            "force" => Some(Q::Force),
            "moment" => Some(Q::Moment),
            "traction" => Some(Q::Pressure),
            "body_force" => Some(Q::ForcePerVolume),
            _ => None,
        };
        if let Some(q) = vector {
            self.field(load, at, "vector", Some(q))?;
        }
        self.field(load, at, "magnitude", Some(Q::Pressure))?;
        self.field(load, at, "acceleration", Some(Q::Acceleration))?;
        self.field(load, at, "angular_velocity", Some(Q::AngularVelocity))?;
        self.field(load, at, "point", Some(Q::Length))?;
        self.plain(load, at, "axis");
        self.field(load, at, "flux", Some(Q::HeatFlux))?;
        self.field(load, at, "film", Some(Q::HeatTransferCoefficient))?;
        self.field(load, at, "sink", Some(Q::Temperature))?;
        self.field(load, at, "power", Some(Q::HeatGeneration))?;
        self.field(load, at, "value", Some(Q::Temperature))
    }

    /// A fluid face or surface condition: velocities and a pressure.
    fn face_condition(&mut self, c: &mut Value, at: &str) -> Result<()> {
        self.field(c, at, "velocity", Some(Q::Velocity))?;
        self.field(c, at, "pressure", Some(Q::Pressure))
    }

    /// A thermal condition: a held temperature or a heat flux, both as `value`.
    fn thermal_condition(&mut self, c: &mut Value, at: &str) -> Result<()> {
        let q = match c.get("type").and_then(Value::as_str) {
            Some("flux") => Q::HeatFlux,
            _ => Q::Temperature,
        };
        self.field(c, at, "value", Some(q))
    }

    fn fluid_domain(&mut self, d: &mut Value, at: &str) -> Result<()> {
        for key in ["min", "max", "spacing", "tolerance"] {
            self.field(d, at, key, Some(Q::Length))?;
        }
        self.field(d, at, "initial_velocity", Some(Q::Velocity))?;
        self.field(d, at, "acceleration", Some(Q::Acceleration))?;
        if let Some(Value::Object(faces)) = d.get_mut("faces") {
            for (face, c) in faces.iter_mut() {
                self.face_condition(c, &format!("{at}/faces/{}", token(face)))?;
            }
        }
        self.each(d, at, "surfaces", |w, s, at| match s.get_mut("condition") {
            Some(c) => w.face_condition(c, &format!("{at}/condition")),
            None => Ok(()),
        })?;
        if let Some(t) = d.get_mut("thermal").filter(|t| t.is_object()) {
            let at = format!("{at}/thermal");
            self.field_or(t, &at, "initial_temperature", Q::Temperature, 0.0)?;
            self.field(t, &at, "gravity", Some(Q::Acceleration))?;
            self.field(t, &at, "expansion", Some(Q::ThermalExpansion))?;
            if let Some(Value::Object(faces)) = t.get_mut("faces") {
                for (face, c) in faces.iter_mut() {
                    self.thermal_condition(c, &format!("{at}/faces/{}", token(face)))?;
                }
            }
            self.each(t, &at, "surfaces", |w, s, at| match s.get_mut("condition") {
                Some(c) => w.thermal_condition(c, &format!("{at}/condition")),
                None => Ok(()),
            })?;
        }
        self.each(d, at, "porous", |w, z, at| {
            w.field(z, at, "min", Some(Q::Length))?;
            w.field(z, at, "max", Some(Q::Length))?;
            if let Some(r) = z.get_mut("resistance") {
                let at = format!("{at}/resistance");
                w.field(r, &at, "permeability", Some(Q::Area))?;
                w.field(r, &at, "inertial", Some(Q::InverseLength))?;
                w.plain(r, &at, "coefficient");
            }
            Ok(())
        })?;
        self.each(d, at, "fans", |w, f, at| {
            w.field(f, at, "min", Some(Q::Length))?;
            w.field(f, at, "max", Some(Q::Length))?;
            // `[Q, Δp]` points: a volume flow and a pressure rise.
            if let Some(Value::Array(points)) = f.get_mut("curve") {
                for (i, point) in points.iter_mut().enumerate() {
                    if let Value::Array(pair) = point {
                        for (j, v) in pair.iter_mut().enumerate() {
                            let q = if j == 0 { Q::VolumeFlow } else { Q::Pressure };
                            w.value(v, &format!("{at}/curve/{i}/{j}"), Some(q))?;
                        }
                    }
                }
            }
            w.plain(f, at, "relaxation");
            Ok(())
        })
    }

    fn step(&mut self, step: &mut Value, at: &str) -> Result<()> {
        let Some(p) = step.get_mut("procedure") else { return Ok(()) };
        let at = format!("{at}/procedure");
        let kind = p.get("type").and_then(Value::as_str).unwrap_or_default().to_string();
        match kind.as_str() {
            "frequency" => {
                self.plain(p, &at, "modes");
                self.field(p, &at, "min_hz", Some(Q::Frequency))?;
                self.field(p, &at, "max_hz", Some(Q::Frequency))?;
            }
            "buckle" => {
                self.plain(p, &at, "modes");
                self.plain(p, &at, "shift");
            }
            "explicit_dynamic" => {
                for key in ["end_time", "output_interval", "history_interval"] {
                    self.field(p, &at, key, Some(Q::Time))?;
                }
                for key in ["dt_scale", "max_increments", "hourglass", "bulk_viscosity"] {
                    self.plain(p, &at, key);
                }
                self.each(p, &at, "initial_velocities", |w, v, at| w.field(v, at, "velocity", Some(Q::Velocity)))?;
                self.each(p, &at, "rigid_walls", |w, wall, at| {
                    w.field(wall, at, "point", Some(Q::Length))?;
                    w.plain(wall, at, "normal");
                    w.friction(wall, at)
                })?;
                if let Some(ms) = p.get_mut("mass_scaling") {
                    self.field(ms, &format!("{at}/mass_scaling"), "target_increment", Some(Q::Time))?;
                }
            }
            "fluid" => {
                self.field(p, &at, "end_time", Some(Q::Time))?;
                self.field(p, &at, "time_step", Some(Q::Time))?;
                for key in ["steps", "frames", "magic"] {
                    self.plain(p, &at, key);
                }
                if let Some(s) = p.get_mut("steady") {
                    let at = format!("{at}/steady");
                    self.field(s, &at, "window", Some(Q::Time))?;
                    self.field(s, &at, "min_time", Some(Q::Time))?;
                    self.plain(s, &at, "tolerance");
                    self.plain(s, &at, "divergence");
                }
            }
            "heat_transfer" => {
                if let Some(t) = p.get_mut("transient").filter(|t| t.is_object()) {
                    let at = format!("{at}/transient");
                    self.field(t, &at, "end_time", Some(Q::Time))?;
                    self.field(t, &at, "increment", Some(Q::Time))?;
                    self.plain(t, &at, "theta");
                    self.field_or(t, &at, "initial_temperature", Q::Temperature, 0.0)?;
                }
            }
            "nonlinear_static" => {
                // Increments are fractions of the step, tolerances relative, the load path load
                // factors.
                for key in ["initial_increment", "min_increment", "max_increment", "max_increments", "max_iterations", "tolerances", "load_path"] {
                    self.plain(p, &at, key);
                }
                if let Some(imp) = p.get_mut("imperfection") {
                    let at = format!("{at}/imperfection");
                    self.plain(imp, &at, "mode");
                    // The buckling mode, scaled to a largest component of one, times `scale`.
                    self.field(imp, &at, "scale", Some(Q::Length))?;
                }
                if let Some(arc) = p.get_mut("arc_length") {
                    let at = format!("{at}/arc_length");
                    self.plain(arc, &at, "max_load_factor");
                    if let Some(limit) = arc.get_mut("displacement_limit") {
                        let at = format!("{at}/displacement_limit");
                        // Degrees of freedom 4 to 6 are rotations, in radians.
                        let rotation = limit.get("direction").and_then(Value::as_u64).is_some_and(|d| d >= 4);
                        self.plain(limit, &at, "node");
                        self.plain(limit, &at, "direction");
                        if rotation {
                            self.plain(limit, &at, "displacement");
                        } else {
                            self.field(limit, &at, "displacement", Some(Q::Length))?;
                        }
                    }
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// Refuse a number the map does not know: it would be read as SI unscaled.
    fn check_all_handled(&self, model: &Value) -> Result<()> {
        let mut numbers = Vec::new();
        collect_numbers(model, String::new(), &mut numbers);
        match numbers.into_iter().find(|p| !self.handled.contains(p)) {
            None => Ok(()),
            Some(p) => invalid(format!(
                "cannot bring `model.json` to format 3: no unit is known for the number at {}",
                p.chars().take(200).collect::<String>()
            )),
        }
    }
}

fn collect_numbers(v: &Value, at: String, out: &mut Vec<String>) {
    match v {
        Value::Number(_) => out.push(at),
        Value::Array(items) => items.iter().enumerate().for_each(|(i, x)| collect_numbers(x, format!("{at}/{i}"), out)),
        Value::Object(m) => m.iter().for_each(|(k, x)| collect_numbers(x, format!("{at}/{}", token(k)), out)),
        _ => {}
    }
}

/// `MCMESH01` → `MCMESH02`: node coordinates scaled by `length` to metres, every other byte as
/// it was. This reads only the `MCMESH01` header and node table (layout of format 2's
/// docs/formats/fieldcraft.md); the rest is checked when the current reader decodes the result.
fn mesh(bytes: &[u8], length: Factor) -> Result<Vec<u8>> {
    const HEADER: usize = 8 + 5 * 8;
    const NODE: usize = 4 + 3 * 8;
    if bytes.get(..8) != Some(b"MCMESH01".as_slice()) {
        return invalid("`mesh.bin` of format 2 must start with MCMESH01");
    }
    let header = bytes.get(..HEADER).ok_or_else(|| Error::Invalid("`mesh.bin` is truncated".into()))?;
    let count = header.get(8..16).and_then(|b| <[u8; 8]>::try_from(b).ok()).map(u64::from_le_bytes).unwrap_or(0);
    let nodes_end = usize::try_from(count)
        .ok()
        .and_then(|n| n.checked_mul(NODE))
        .and_then(|n| n.checked_add(HEADER))
        .filter(|end| *end <= bytes.len())
        .ok_or_else(|| Error::Invalid("`mesh.bin` lists more nodes than it holds".into()))?;
    let mut out = Vec::with_capacity(bytes.len());
    out.extend_from_slice(b"MCMESH02");
    out.extend_from_slice(header.get(8..).unwrap_or_default());
    for record in bytes.get(HEADER..nodes_end).unwrap_or_default().as_chunks::<NODE>().0 {
        let (label, coordinates) = record.split_at(4);
        out.extend_from_slice(label);
        for c in coordinates.as_chunks::<8>().0 {
            let x = f64::from_le_bytes(*c);
            let si = length.to_si(x);
            if !si.is_finite() {
                return invalid("`mesh.bin` holds a node coordinate that is not finite in metres");
            }
            out.extend_from_slice(&si.to_le_bytes());
        }
    }
    out.extend_from_slice(bytes.get(nodes_end..).unwrap_or_default());
    Ok(out)
}

/// Scaling by powers of ten and adding offsets on the decimal a value was written as: the
/// algorithm of the unit table of the release that made format 3 (fieldcraft-model
/// `units::decimal`), kept here as it was.
mod decimal {
    const EXACT_DIGITS: usize = 15;
    const POW10: [f64; 23] =
        [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22];

    fn short_decimal(v: f64) -> Option<(i128, i32)> {
        if v == 0.0 || !v.is_finite() {
            return None;
        }
        let text = format!("{v:e}");
        let (mantissa, exponent) = text.split_once('e')?;
        let exponent: i32 = exponent.parse().ok()?;
        let (negative, mantissa) = match mantissa.strip_prefix('-') {
            Some(m) => (true, m),
            None => (false, mantissa),
        };
        let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        if whole.len() + fraction.len() > EXACT_DIGITS {
            return None;
        }
        let mut m: i128 = 0;
        for c in whole.chars().chain(fraction.chars()) {
            m = m.checked_mul(10)?.checked_add(i128::from(c.to_digit(10)?))?;
        }
        let fraction_len = i32::try_from(fraction.len()).ok()?;
        Some((if negative { -m } else { m }, exponent.checked_sub(fraction_len)?))
    }

    fn from_decimal(mantissa: i128, exponent: i32) -> f64 {
        format!("{mantissa}e{exponent}").parse().unwrap_or(f64::NAN)
    }

    pub(super) fn shift(v: f64, k: i32) -> f64 {
        if k == 0 || v == 0.0 || !v.is_finite() {
            return v;
        }
        if let Some((m, e)) = short_decimal(v)
            && let Some(e) = e.checked_add(k)
        {
            return from_decimal(m, e);
        }
        let (mut v, mut k) = (v, k);
        let step = POW10.len() as i32 - 1;
        while k != 0 && v.is_finite() && v != 0.0 {
            let n = k.clamp(-step, step);
            let p = POW10.get(n.unsigned_abs() as usize).copied().unwrap_or(1.0);
            v = if n > 0 { v * p } else { v / p };
            k -= n;
        }
        v
    }

    pub(super) fn add(a: f64, b: f64) -> f64 {
        if a == 0.0 || b == 0.0 {
            return a + b;
        }
        let (Some((ma, ea)), Some((mb, eb))) = (short_decimal(a), short_decimal(b)) else { return a + b };
        let e = ea.min(eb);
        let scale = |m: i128, from: i32| -> Option<i128> {
            let n = u32::try_from(from - e).ok()?;
            m.checked_mul(10i128.checked_pow(n)?)
        };
        match (scale(ma, ea), scale(mb, eb)) {
            (Some(x), Some(y)) => match x.checked_add(y) {
                Some(0) => 0.0,
                Some(sum) => from_decimal(sum, e),
                None => a + b,
            },
            _ => a + b,
        }
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    fn parts(model: Value, mesh: Vec<u8>) -> Parts {
        Parts { manifest: json!({"format": "fieldcraft/2", "generator": "FieldCraft 0.7.0"}), model, mesh, journal: json!([]) }
    }

    /// An `MCMESH01` part with `nodes` and no elements.
    fn mesh01(nodes: &[(u32, [f64; 3])]) -> Vec<u8> {
        let mut b = b"MCMESH01".to_vec();
        for n in [nodes.len() as u64, 0, 0, 0, 0] {
            b.extend_from_slice(&n.to_le_bytes());
        }
        for (label, xyz) in nodes {
            b.extend_from_slice(&label.to_le_bytes());
            xyz.iter().for_each(|c| b.extend_from_slice(&c.to_le_bytes()));
        }
        b
    }

    #[test]
    fn scales_by_the_declared_system() {
        let mut p = parts(
            json!({"units": "mm-t-s-N-MPa", "materials": {"Steel": {"elastic": {"type": "isotropic", "young": 210000.0, "poisson": 0.3}, "density": 7.85e-9,
                "thermal": {"expansion": 1.2e-5, "conductivity": 0.05, "specific_heat": 4.6e8}}},
                "loads": {"g": {"type": "gravity", "acceleration": [0.0, 0.0, -9810.0]}, "hot": {"type": "temperature", "target": {"set": "all"}, "value": 20.0}},
                "datums": {"Mid": {"type": "plane", "point": [10.0, 0.0, 0.0], "normal": [1.0, 0.0, 0.0]}}}),
            mesh01(&[(1, [40.0, -0.0, 3.97])]),
        );
        v2_to_v3(&mut p).unwrap();
        let m = &p.model;
        assert_eq!(m["display_units"], "mm-t-s-N-MPa");
        assert!(m.get("units").is_none());
        let steel = &m["materials"]["Steel"];
        assert_eq!((steel["elastic"]["young"].as_f64(), steel["elastic"]["poisson"].as_f64()), (Some(2.1e11), Some(0.3)));
        assert_eq!(steel["density"].as_f64(), Some(7850.0));
        assert_eq!(steel["thermal"]["specific_heat"].as_f64(), Some(460.0));
        assert_eq!(steel["thermal"]["conductivity"].as_f64(), Some(0.05));
        assert_eq!(steel["thermal"]["expansion"].as_f64(), Some(1.2e-5));
        // A reference temperature left out meant 0 °C.
        assert_eq!(steel["thermal"]["reference_temperature"].as_f64(), Some(273.15));
        assert_eq!(m["loads"]["g"]["acceleration"], json!([0.0, 0.0, -9.81]));
        assert_eq!(m["loads"]["hot"]["value"].as_f64(), Some(293.15));
        assert_eq!(m["datums"]["Mid"]["point"], json!([0.01, 0.0, 0.0]));
        assert_eq!(&p.mesh[..8], b"MCMESH02");
        let coords: Vec<f64> = p.mesh[48 + 4..48 + 28].chunks(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
        assert_eq!(coords, [0.04, -0.0, 0.00397]);
        assert!(coords[1].is_sign_negative());
        assert_eq!(p.journal, json!({"units": "mm-t-s-N-MPa", "model_units": "mm-t-s-N-MPa", "commands": []}));
        assert_eq!(p.manifest["format"], "fieldcraft/3");
    }

    #[test]
    fn metres_stay_and_datums_were_millimetres() {
        let mut p = parts(
            json!({"units": "m-kg-s-N-Pa", "reference_points": {"RP": {"xyz": [0.1, 0.0, 2.0], "label": 7}},
                "boundary_conditions": {"b": {"type": "displacement", "target": {"set": "s"}, "values": [0.001, null, 0.0, 0.5, null, null]},
                    "t": {"type": "temperature", "target": {"set": "s"}, "value": 300.0}},
                "datums": {"P": {"type": "point", "point": [5.0, 0.0, 0.0]}}}),
            mesh01(&[(1, [0.1, 0.05, 0.005])]),
        );
        v2_to_v3(&mut p).unwrap();
        assert_eq!(p.model["reference_points"]["RP"]["xyz"], json!([0.1, 0.0, 2.0]));
        assert_eq!(p.model["reference_points"]["RP"]["label"], 7);
        // Rotations are radians; temperatures were kelvin.
        assert_eq!(p.model["boundary_conditions"]["b"]["values"], json!([0.001, null, 0.0, 0.5, null, null]));
        assert_eq!(p.model["boundary_conditions"]["t"]["value"].as_f64(), Some(300.0));
        assert_eq!(p.model["datums"]["P"]["point"], json!([0.005, 0.0, 0.0]));
        assert_eq!(p.mesh[..8], *b"MCMESH02");
        assert_eq!(p.mesh[8..], mesh01(&[(1, [0.1, 0.05, 0.005])])[8..]);
    }

    #[test]
    fn unspecified_reads_as_millimetre_tonne_second() {
        let mut p = parts(json!({"sections": {"s": {"set": "a", "material": "m", "kind": {"type": "plane_strain"}}}}), mesh01(&[]));
        v2_to_v3(&mut p).unwrap();
        assert_eq!(p.model["display_units"], "mm-t-s-N-MPa");
        // The default thickness of 1 was 1 mm.
        assert_eq!(p.model["sections"]["s"]["kind"]["thickness"].as_f64(), Some(0.001));
        assert_eq!(p.journal["units"], "mm-t-s-N-MPa");
    }

    #[test]
    fn inch_pound_and_millisecond_systems() {
        let mut p = parts(
            json!({"units": "in-lbf-s-psi", "loads": {"t": {"type": "temperature", "target": {"set": "s"}, "value": 68.0}, "p": {"type": "pressure", "set": "s", "magnitude": 1.0}}}),
            mesh01(&[(1, [1.0, 0.0, 0.0])]),
        );
        v2_to_v3(&mut p).unwrap();
        assert!((p.model["loads"]["t"]["value"].as_f64().unwrap() - 293.15).abs() < 1e-12);
        assert_eq!(p.model["loads"]["p"]["magnitude"].as_f64(), Some(6894.757293168361));
        let x = f64::from_le_bytes(p.mesh[52..60].try_into().unwrap());
        assert_eq!(x, 0.0254);
        let mut p = parts(
            json!({"units": "mm-kg-ms-kN-GPa", "steps": [{"name": "Modes", "procedure": {"type": "frequency", "modes": 3, "max_hz": 2.0}, "output": []}],
                "materials": {"Steel": {"density": 7.85e-6, "elastic": {"type": "isotropic", "young": 210.0, "poisson": 0.3}}}}),
            mesh01(&[]),
        );
        v2_to_v3(&mut p).unwrap();
        // Cycles per millisecond.
        assert_eq!(p.model["steps"][0]["procedure"]["max_hz"].as_f64(), Some(2000.0));
        assert_eq!(p.model["steps"][0]["procedure"]["modes"], 3);
        assert_eq!(p.model["materials"]["Steel"]["density"].as_f64(), Some(7850.0));
        assert_eq!(p.model["materials"]["Steel"]["elastic"]["young"].as_f64(), Some(2.1e11));
    }

    #[test]
    fn refuses_what_it_cannot_scale() {
        for (model, why) in [
            (json!({"units": "furlong-firkin-fortnight"}), "unknown unit system"),
            (json!({"units": 3}), "not a unit system"),
            (json!({"materials": {"m": {"elastic": {"type": "isotropic", "young": 1.0, "poisson": 0.3, "warp": 2.0}}}}), "/materials/m/elastic/warp"),
            (json!({"novel": {"x": 1.0}}), "/novel/x"),
            (json!({"loads": {"l": {"type": "force", "target": {"set": "s"}, "vector": [1e308, 0.0, 0.0]}}}), ""),
            (json!([]), "not a JSON object"),
        ] {
            let mut p = parts(model.clone(), mesh01(&[]));
            let before = p.clone();
            let e = v2_to_v3(&mut p);
            if why.is_empty() {
                // 1e308 N is finite in SI too.
                assert!(e.is_ok(), "{model}");
                continue;
            }
            let e = e.unwrap_err().to_string();
            assert!(e.contains(why), "{model}: {e}");
            let _ = before;
        }
        let mut p = parts(json!({"units": "mm-t-s-N-MPa", "loads": {"l": {"type": "pressure", "set": "s", "magnitude": 1e305}}}), mesh01(&[]));
        assert!(v2_to_v3(&mut p).unwrap_err().to_string().contains("out of range"));
        for bad in [b"MCMESH02".to_vec(), b"MCMESH01".to_vec(), {
            let mut m = mesh01(&[(1, [0.0; 3])]);
            m[8] = 9;
            m
        }] {
            let mut p = parts(json!({}), bad);
            assert!(v2_to_v3(&mut p).is_err());
        }
        let mut p = parts(json!({}), mesh01(&[]));
        p.journal = json!({"not": "a list"});
        assert!(v2_to_v3(&mut p).is_err());
    }

    /// Every released format 1 file (through 1 → 2) takes the step: the map knows every number
    /// in it, and the mesh and the journal come out in format 3's form.
    #[test]
    fn every_format_1_fixture_takes_the_step() {
        use std::io::Read;
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/formats/v1");
        let mut seen = 0;
        for writer in std::fs::read_dir(&root).unwrap().map(|e| e.unwrap().path()).filter(|p| p.is_dir()) {
            for path in std::fs::read_dir(&writer).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|x| x == "fieldcraft")) {
                let bytes = std::fs::read(&path).unwrap();
                let mut zip = zip::ZipArchive::new(std::io::Cursor::new(&bytes)).unwrap();
                let mut part = |name: &str| {
                    let mut b = Vec::new();
                    zip.by_name(name).unwrap().read_to_end(&mut b).unwrap();
                    b
                };
                let json = |b: Vec<u8>| -> Value { serde_json::from_slice(&b).unwrap() };
                let mut p = Parts {
                    manifest: json(part("fieldcraft.json")),
                    model: json(part("model.json")),
                    mesh: part("mesh.bin"),
                    journal: json(part("journal.json")),
                };
                super::super::step(&mut p, 1).unwrap();
                let journal = p.journal.clone();
                v2_to_v3(&mut p).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
                assert_eq!(p.manifest["format"], "fieldcraft/3");
                assert_eq!(p.journal["commands"], journal, "{}: the journal's commands are untouched", path.display());
                assert_eq!(&p.mesh[..8], b"MCMESH02");
                seen += 1;
            }
        }
        assert!(seen >= 30, "{seen}");
    }
}
