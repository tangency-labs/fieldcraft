//! A small, cited catalogue of metal and fluid properties. Values are stored in SI, like the
//! model; display units are only used when reporting them. Strengths are information, never an
//! automatically selected plasticity law or failure criterion. Applicability belongs to each
//! value: a room-temperature constant does not describe temperature-dependent behaviour.
//!
//! Temperature-dependent properties are curves: points (kelvin, SI value) taken from a published
//! table, or a published equation evaluated at the listed temperatures, cited like the constants.
//! A material built from an entry takes each curve's value at a temperature the user chooses as
//! its constant and carries the curve as a [`crate::curves::TemperatureTable`], saved with the
//! model.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use crate::curves::{TemperatureDependent, TemperatureTable};
use crate::units::{DisplaySystem, Quantity};
use crate::{Directional, Elastic, Error, Fluid, Material, Result, Thermal};

const MIN_YOUNG: f64 = 1e9; // unit-factors: SI plausibility lower bound in Pa, not a conversion
const MAX_YOUNG: f64 = 1e12; // unit-factors: SI plausibility upper bound in Pa, not a conversion
const MAX_EXPANSION: f64 = 1e-3; // unit-factors: SI plausibility upper bound in inverse kelvin, not a conversion
const MIN_STRENGTH: f64 = 1e6; // unit-factors: SI plausibility lower bound in Pa, not a conversion
const MIN_VISCOSITY: f64 = 1e-6; // unit-factors: SI plausibility lower bound in Pa s, not a conversion
const MAX_VISCOSITY: f64 = 1e3; // unit-factors: SI plausibility upper bound in Pa s, not a conversion

/// The temperature the constants of a material are taken at unless the user chooses one, and the
/// default stress-free temperature: 20 degrees Celsius, a modelling default.
pub const DEFAULT_TEMPERATURE: f64 = crate::units::CELSIUS_ZERO + 20.0;

/// Stable entry names and labels, also used by command forms.
pub const CHOICES: &[(&str, &str)] = &[
    ("EN10025-S235", "S235 structural steel"),
    ("EN10025-S355", "S355 structural steel"),
    ("EN10088-14301", "1.4301 / 304 stainless steel"),
    ("EN10088-14401", "1.4401 / 316 stainless steel"),
    ("AA6061-T6", "6061-T6 aluminium"),
    ("AA6082-T6", "6082-T6 aluminium"),
    ("AA7075-T6", "7075-T6 aluminium"),
    ("AA2024-T3", "2024-T3 aluminium"),
    ("Ti6Al4V-annealed", "Ti-6Al-4V annealed"),
    ("air", "Air, dry"),
];

pub const FAMILIES: &[(&str, &str)] =
    &[("steel", "Structural steel"), ("stainless", "Stainless steel"), ("aluminium", "Aluminium"), ("titanium", "Titanium"), ("fluid", "Fluid")];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Property {
    Young,
    Poisson,
    Density,
    Conductivity,
    SpecificHeat,
    Expansion,
    Yield,
    Tensile,
    Viscosity,
}

impl Property {
    pub fn quantity(self) -> Quantity {
        match self {
            Self::Young | Self::Yield | Self::Tensile => Quantity::Pressure,
            Self::Poisson => Quantity::Dimensionless,
            Self::Density => Quantity::Density,
            Self::Conductivity => Quantity::ThermalConductivity,
            Self::SpecificHeat => Quantity::SpecificHeat,
            Self::Expansion => Quantity::ThermalExpansion,
            Self::Viscosity => Quantity::DynamicViscosity,
        }
    }

    /// The property's name in messages.
    pub fn name(self) -> &'static str {
        match self {
            Self::Young => "Young's modulus",
            Self::Poisson => "Poisson's ratio",
            Self::Density => "density",
            Self::Conductivity => "conductivity",
            Self::SpecificHeat => "specific heat",
            Self::Expansion => "expansion",
            Self::Yield => "yield strength",
            Self::Tensile => "tensile strength",
            Self::Viscosity => "viscosity",
        }
    }

    /// Every property, in report order.
    pub const ALL: [Property; 9] = [
        Property::Young,
        Property::Poisson,
        Property::Density,
        Property::Conductivity,
        Property::SpecificHeat,
        Property::Expansion,
        Property::Viscosity,
        Property::Yield,
        Property::Tensile,
    ];

    pub fn applied(self) -> bool {
        !matches!(self, Self::Yield | Self::Tensile)
    }

    /// Whether a solver follows this property's temperature curve (a material's
    /// [`crate::TemperatureDependent`] tables): any other curve would be read as its constant
    /// alone.
    pub fn follows_temperature(self) -> bool {
        matches!(self, Self::Density | Self::Conductivity | Self::SpecificHeat | Self::Viscosity)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Basis {
    Design,
    Minimum,
    Typical,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub document: String,
    pub edition: String,
    pub location: String,
    /// The page, when it was checked against the document itself; without one the location
    /// (clause, table or subject) identifies the place.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub page: String,
    /// Optional publisher/catalogue link; document, edition, location and page identify the source.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Datum {
    pub value: f64,
    pub unit: String,
    pub basis: Basis,
    pub source: Source,
    pub applicability: String,
}

/// A property as a function of the temperature: `points` are (kelvin, SI value), the source's
/// table or its equation evaluated at the listed temperatures; linear between them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Curve {
    pub unit: String,
    pub points: TemperatureTable,
    pub basis: Basis,
    pub source: Source,
    pub applicability: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Designation {
    pub standard: String,
    pub grade: String,
    pub condition: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Entry {
    pub id: String,
    pub name: String,
    pub family: String,
    pub designation: Designation,
    pub properties: BTreeMap<Property, Datum>,
    /// Temperature-dependent properties; a property is a constant or a curve, never both.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub curves: BTreeMap<Property, Curve>,
}

/// A fixed, embedded catalogue; no files or network are consulted at runtime.
pub fn entries() -> Result<&'static [Entry]> {
    static ENTRIES: OnceLock<std::result::Result<Vec<Entry>, String>> = OnceLock::new();
    ENTRIES
        .get_or_init(|| {
            let values: Vec<Entry> = serde_json::from_str(include_str!("../data/materials.json")).map_err(|e| e.to_string())?;
            validate(&values).map_err(|e| e.to_string())?;
            Ok(values)
        })
        .as_ref()
        .map(Vec::as_slice)
        .map_err(|e| Error::Invalid(format!("the bundled material library is not valid: {e}")))
}

fn validate(values: &[Entry]) -> Result<()> {
    let mut ids = BTreeSet::new();
    if values.len() != CHOICES.len() {
        return Err(Error::Invalid("the material library and its command choices differ".into()));
    }
    for e in values {
        if !ids.insert(&e.id) || !CHOICES.iter().any(|(id, _)| *id == e.id) {
            return Err(Error::Invalid("material library ids must be unique and match the fixed catalogue".into()));
        }
        if e.name.trim().is_empty()
            || e.designation.standard.trim().is_empty()
            || e.designation.grade.trim().is_empty()
            || e.designation.condition.trim().is_empty()
            || !FAMILIES.iter().any(|(id, _)| *id == e.family)
        {
            return Err(Error::Invalid(format!("material library entry {} needs a designation and family", e.id)));
        }
        let cited = |unit: &str, q: Quantity, source: &Source, applicability: &str| {
            unit == q.si_symbol()
                && !applicability.trim().is_empty()
                && [&source.document, &source.edition, &source.location].iter().all(|s| !s.trim().is_empty())
                && (source.page.is_empty() || !source.page.trim().is_empty())
                && source.url.as_ref().is_none_or(|url| url.starts_with("https://"))
        };
        let uncited = || Error::Invalid(format!("material library entry {} has an uncited or non-SI property", e.id));
        for (&p, d) in &e.properties {
            if !(cited(&d.unit, p.quantity(), &d.source, &d.applicability) && d.value.is_finite()) {
                return Err(uncited());
            }
            plausible(e, p, d.value)?;
        }
        for (&p, c) in &e.curves {
            if !cited(&c.unit, p.quantity(), &c.source, &c.applicability) || e.properties.contains_key(&p) || !p.follows_temperature() {
                return Err(uncited());
            }
            c.points.validate(&format!("{} {}", e.id, p.name()))?;
            for &[_, v] in &c.points.0 {
                plausible(e, p, v)?;
            }
        }
        e.material(crate::units::CELSIUS_ZERO, DEFAULT_TEMPERATURE)?;
    }
    Ok(())
}

fn plausible(e: &Entry, p: Property, value: f64) -> Result<()> {
    let plausible = match p {
        Property::Young => (MIN_YOUNG..=MAX_YOUNG).contains(&value),
        Property::Poisson => (0.0..0.5).contains(&value),
        // A gas is light; a solid or a liquid is not.
        Property::Density if e.is_fluid() => (0.01..=30_000.0).contains(&value),
        Property::Density => (100.0..=30_000.0).contains(&value),
        Property::Conductivity => (0.01..=2000.0).contains(&value),
        Property::SpecificHeat => (1.0..=20_000.0).contains(&value),
        Property::Expansion => (0.0..=MAX_EXPANSION).contains(&value),
        Property::Viscosity => (MIN_VISCOSITY..=MAX_VISCOSITY).contains(&value),
        Property::Yield | Property::Tensile => (MIN_STRENGTH..=1e10).contains(&value),
    };
    if plausible { Ok(()) } else { Err(Error::Invalid(format!("material library entry {} has an implausible property", e.id))) }
}

/// Find an exact stable id. Errors include the small, bounded catalogue, never hostile input.
pub fn entry(id: &str) -> Result<&'static Entry> {
    if id.len() <= 64
        && let Some(e) = entries()?.iter().find(|e| e.id == id)
    {
        return Ok(e);
    }
    Err(Error::Invalid(format!(
        "unknown material library entry; choose an id from model.material.library.list: {}",
        CHOICES.iter().map(|(id, _)| *id).collect::<Vec<_>>().join(", ")
    )))
}

impl Entry {
    /// A constant, or a curve's value at `temperature` (kelvin; its end value outside it).
    fn value_at(&self, property: Property, temperature: f64) -> Option<f64> {
        self.value(property).or_else(|| self.curves.get(&property).map(|c| c.points.value(temperature)))
    }

    fn value(&self, property: Property) -> Option<f64> {
        self.properties.get(&property).map(|d| d.value)
    }

    /// The curves `temperature` (kelvin) lies outside of, where their end values are taken.
    pub fn outside(&self, temperature: f64) -> Vec<Property> {
        self.curves.iter().filter(|(_, c)| !c.points.covers(temperature)).map(|(p, _)| *p).collect()
    }

    /// A fluid: no elastic constants, a viscosity.
    pub fn is_fluid(&self) -> bool {
        self.family == "fluid"
    }

    /// Build an SI material: constants as given, each curve's value at `temperature` (kelvin) as
    /// its constant with the curve kept as the material's temperature table. The stress-free
    /// reference temperature is a modelling choice made by the caller, not an asserted material
    /// property from a handbook, and so is the temperature the constants are taken at.
    pub fn material(&self, reference_temperature: f64, temperature: f64) -> Result<Material> {
        if !reference_temperature.is_finite() || reference_temperature < 0.0 {
            return Err(Error::Invalid("the stress-free reference temperature must be at or above zero kelvin".into()));
        }
        if !temperature.is_finite() || temperature < 0.0 {
            return Err(Error::Invalid("the temperature of the constant properties must be at or above zero kelvin".into()));
        }
        let (young, poisson) = (self.value(Property::Young), self.value(Property::Poisson));
        let elastic = match (young, poisson) {
            (Some(young), Some(poisson)) => Some(Elastic::Isotropic { young, poisson }),
            (None, None) if self.is_fluid() => None,
            _ => return Err(Error::Invalid(format!("material library entry {} needs both isotropic elastic constants", self.id))),
        };
        let fluid = match self.value_at(Property::Viscosity, temperature) {
            Some(viscosity) => Some(Fluid::Newtonian { viscosity }),
            None if self.is_fluid() => return Err(Error::Invalid(format!("material library entry {} needs a viscosity", self.id))),
            None => None,
        };
        let conductivity = self.value_at(Property::Conductivity, temperature).map(Directional::Isotropic);
        let specific_heat = self.value_at(Property::SpecificHeat, temperature);
        let expansion = self.value(Property::Expansion).map(Directional::Isotropic);
        let thermal = (conductivity.is_some() || specific_heat.is_some() || expansion.is_some()).then_some(Thermal {
            conductivity,
            specific_heat,
            expansion,
            reference_temperature,
        });
        let table = |p| self.curves.get(&p).map(|c| c.points.clone());
        let tables = TemperatureDependent {
            conductivity: table(Property::Conductivity),
            specific_heat: table(Property::SpecificHeat),
            density: table(Property::Density),
            viscosity: table(Property::Viscosity),
        };
        let material = Material {
            elastic,
            density: self.value_at(Property::Density, temperature),
            fluid,
            thermal,
            temperature_dependent: (!tables.is_empty()).then_some(tables),
            ..Material::default()
        };
        material.validate()?;
        Ok(material)
    }

    /// Every supplied datum, expressed through the same public quantity boundary as commands.
    pub fn report(&self, display: &DisplaySystem) -> Result<Vec<ReportedValue>> {
        self.properties
            .iter()
            .map(|(&property, datum)| {
                let q = property.quantity();
                Ok(ReportedValue {
                    property,
                    value: display.express(q, datum.value).map_err(|e| Error::Invalid(e.to_string()))?,
                    unit: display.symbol(q).to_string(),
                    basis: datum.basis,
                    source: datum.source.clone(),
                    applicability: datum.applicability.clone(),
                })
            })
            .collect()
    }

    /// Every curve, its points expressed through the same public quantity boundary.
    pub fn report_curves(&self, display: &DisplaySystem) -> Result<Vec<ReportedCurve>> {
        let express = |q, v| display.express(q, v).map_err(|e| Error::Invalid(e.to_string()));
        self.curves
            .iter()
            .map(|(&property, c)| {
                let points = c
                    .points
                    .0
                    .iter()
                    .map(|&[t, v]| Ok([express(Quantity::Temperature, t)?, express(property.quantity(), v)?]))
                    .collect::<Result<_>>()?;
                Ok(ReportedCurve {
                    property,
                    points,
                    temperature_unit: display.symbol(Quantity::Temperature).to_string(),
                    unit: display.symbol(property.quantity()).to_string(),
                    basis: c.basis,
                    source: c.source.clone(),
                    applicability: c.applicability.clone(),
                })
            })
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ReportedCurve {
    pub property: Property,
    /// (temperature, value) in the display units.
    pub points: Vec<[f64; 2]>,
    pub temperature_unit: String,
    pub unit: String,
    pub basis: Basis,
    pub source: Source,
    pub applicability: String,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ReportedValue {
    pub property: Property,
    pub value: f64,
    pub unit: String,
    pub basis: Basis,
    pub source: Source,
    pub applicability: String,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::{CELSIUS_ZERO, Preset, Unit};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1e-20)
    }

    #[test]
    fn catalogue_is_cited_and_materials_remain_si_in_every_display_system() {
        let all = entries().unwrap();
        assert_eq!(all.len(), 10);
        validate(all).unwrap();
        for e in all {
            let material = e.material(293.15, 293.15).unwrap();
            assert!(material.plastic.is_none() && material.damage.is_none() && material.erosion.is_none());
            for preset in Preset::ALL {
                let display = DisplaySystem::new(preset);
                for r in e.report(&display).unwrap() {
                    let datum = e.properties.get(&r.property).unwrap();
                    assert!(close(display.to_si(r.property.quantity(), r.value).unwrap(), datum.value), "{}: {:?}", e.id, r.property);
                    assert_eq!((&r.source, r.basis, &r.applicability), (&datum.source, datum.basis, &datum.applicability));
                }
                assert_eq!(material, e.material(293.15, 293.15).unwrap());
                material.validate().unwrap();
            }
        }
    }

    #[test]
    fn a_missing_citation_invalid_property_or_duplicate_id_is_refused() {
        let original = entries().unwrap().to_vec();
        let mut bad = original.clone();
        bad.first_mut().unwrap().properties.get_mut(&Property::Young).unwrap().source.location.clear();
        assert!(validate(&bad).is_err());
        let mut bad = original.clone();
        bad.first_mut().unwrap().properties.get_mut(&Property::Young).unwrap().value = f64::NAN;
        assert!(validate(&bad).is_err());
        let mut bad = original.clone();
        bad.first_mut().unwrap().properties.remove(&Property::Poisson);
        assert!(validate(&bad).is_err());
        let mut bad = original;
        bad.get_mut(1).unwrap().id = bad.first().unwrap().id.clone();
        assert!(validate(&bad).is_err());
        assert!(entry(&"x".repeat(100_000)).unwrap_err().to_string().len() < 1000);
        assert!(entry("S355").is_err());
    }

    #[test]
    fn citations_without_a_download_link_survive_reports_and_serialization() {
        let steel = entry("EN10025-S235").unwrap();
        let source = &steel.properties[&Property::Young].source;
        assert_eq!(source.document, "EN 1993-1-1");
        assert!(source.url.is_none());
        let json = serde_json::to_value(source).unwrap();
        assert!(json.get("url").is_none());
        assert_eq!(&serde_json::from_value::<Source>(json).unwrap(), source);
        let report = steel.report(&DisplaySystem::new(Preset::MmTonneS)).unwrap();
        assert_eq!(&report.iter().find(|r| r.property == Property::Young).unwrap().source, source);

        let mut bad = entries().unwrap().to_vec();
        bad.first_mut().unwrap().properties.get_mut(&Property::Young).unwrap().source.url = Some("http://example.org".into());
        assert!(validate(&bad).is_err());
    }

    #[test]
    fn structural_and_thermal_reports_use_existing_quantity_conversions() {
        let mm = DisplaySystem::new(Preset::MmTonneS);
        let steel = entry("EN10025-S355").unwrap();
        let report = steel.report(&mm).unwrap();
        let value = |property| report.iter().find(|r| r.property == property).unwrap().value;
        assert!(close(value(Property::Young), 210_000.0));
        assert!(close(value(Property::Density), 7.85e-9));
        assert!(close(value(Property::Expansion), 12e-6));
        let stainless = entry("EN10088-14301").unwrap();
        let report = stainless.report(&mm).unwrap();
        assert!(close(report.iter().find(|r| r.property == Property::Conductivity).unwrap().value, 15.0));
        assert!(close(report.iter().find(|r| r.property == Property::SpecificHeat).unwrap().value, 5e8));
        let fahrenheit = mm.with_override(Quantity::Temperature, "degF").unwrap();
        for r in stainless.report(&fahrenheit).unwrap() {
            let datum = stainless.properties.get(&r.property).unwrap();
            assert!(close(fahrenheit.to_si(r.property.quantity(), r.value).unwrap(), datum.value));
        }
        // Original published non-SI values are converted by the shared parser, not local factors.
        let ksi = Unit::parse("ksi").unwrap();
        assert!(close(entry("AA2024-T3").unwrap().value(Property::Young).unwrap(), Quantity::Pressure.to_si(10_500.0, &ksi).unwrap()));
        assert!(close(entry("AA7075-T6").unwrap().value(Property::Young).unwrap(), Quantity::Pressure.to_si(10_300.0, &ksi).unwrap()));
        assert!(close(
            entry("AA2024-T3").unwrap().value(Property::Density).unwrap(),
            Quantity::Density.to_si(0.100, &Unit::parse("lb/in^3").unwrap()).unwrap()
        ));
        assert!(stainless.material(-1.0, 293.15).is_err());
        assert!(stainless.material(f64::INFINITY, 293.15).is_err());
        assert!(stainless.material(293.15, f64::NAN).is_err());
    }

    /// EN 1993-1-2's carbon-steel specific heat (J/(kg K)) and conductivity (W/(m K)) at `c` °C,
    /// each branch as the standard bounds it.
    fn steel_specific_heat(c: f64) -> f64 {
        match c {
            c if c < 600.0 => 425.0 + 7.73e-1 * c - 1.69e-3 * c * c + 2.22e-6 * c * c * c,
            c if c < 735.0 => 666.0 + 13002.0 / (738.0 - c),
            c if c < 900.0 => 545.0 + 17820.0 / (c - 731.0),
            _ => 650.0,
        }
    }

    fn steel_conductivity(c: f64) -> f64 {
        if c < 800.0 { 54.0 - 3.33e-2 * c } else { 27.3 }
    }

    fn curve(id: &str, p: Property) -> &'static Curve {
        entry(id).unwrap().curves.get(&p).unwrap()
    }

    /// Each point against its source: the published table, or the published equation.
    fn each_point(id: &str, p: Property, source: impl Fn(f64) -> f64) {
        let c = curve(id, p);
        for &[t, v] in &c.points.0 {
            assert!(close(v, source(t)), "{id} {p:?} at {t} K: {v} against {}", source(t));
        }
    }

    #[test]
    fn curves_match_their_sources_at_the_tabulated_points() {
        // U.S. Standard Atmosphere 1976: eq. (51) with β and S of §1.3.11, eq. (53), and the
        // perfect gas eq. (42) at P0 with R* and M0.
        let (beta, s, r, m0, p0) = (1.458e-6, 110.4, 8.31432e3, 28.9644, 1.013250e5);
        each_point("air", Property::Viscosity, |t| beta * t.powf(1.5) / (t + s));
        each_point("air", Property::Conductivity, |t| 2.64638e-3 * t.powf(1.5) / (t + 245.4 * 10f64.powf(-12.0 / t)));
        each_point("air", Property::Density, |t| p0 * m0 / (r * t));
        assert!(close(entry("air").unwrap().properties[&Property::SpecificHeat].value, 1.4 / 0.4 * r / m0));

        // EN 1993-1-2 3.4.1 (carbon steel), in °C.
        let celsius = |f: fn(f64) -> f64| move |t: f64| f(t - CELSIUS_ZERO);
        for id in ["EN10025-S235", "EN10025-S355"] {
            each_point(id, Property::SpecificHeat, celsius(steel_specific_heat));
            each_point(id, Property::Conductivity, celsius(steel_conductivity));
            // Between the points too: within 0.1 % of the equations, through the 735 °C peak.
            let c = curve(id, Property::SpecificHeat);
            for tenth in 200..12000 {
                let t = CELSIUS_ZERO + f64::from(tenth) / 10.0;
                let want = steel_specific_heat(t - CELSIUS_ZERO);
                assert!((c.points.value(t) - want).abs() <= 1e-3 * want + 1e-9, "{id} at {t} K: {} against {want}", c.points.value(t));
            }
            assert_eq!(c.points.value(CELSIUS_ZERO + 735.0), 5000.0);
        }
        for (id, p, [lo, hi]) in [
            ("EN10025-S235", Property::SpecificHeat, [20.0, 1200.0]),
            ("EN10025-S235", Property::Conductivity, [20.0, 1200.0]),
            ("air", Property::Density, [200.0 - CELSIUS_ZERO, 400.0 - CELSIUS_ZERO]),
        ] {
            let range = curve(id, p).points.range();
            assert!(close(range[0], lo + CELSIUS_ZERO) && close(range[1], hi + CELSIUS_ZERO), "{id} {p:?} {range:?}");
        }
    }

    #[test]
    fn materials_take_their_constants_at_the_chosen_temperature_and_keep_the_curves() {
        let air = entry("air").unwrap();
        // Inside the table: linear between 300 K and 310 K.
        let m = air.material(293.15, 305.0).unwrap();
        let mu = |t: f64| 1.458e-6 * t.powf(1.5) / (t + 110.4);
        assert!(close(m.fluid.unwrap().newtonian().unwrap(), 0.5 * (mu(300.0) + mu(310.0))));
        assert!(m.elastic.is_none() && m.plastic.is_none());
        let tables = m.temperature_dependent.as_ref().unwrap();
        assert_eq!(tables.viscosity.as_ref().unwrap(), &curve("air", Property::Viscosity).points);
        assert!(tables.specific_heat.is_none(), "the perfect-gas specific heat is a constant");
        assert!(air.outside(305.0).is_empty());
        // Outside: the end value, and the curves say so.
        let hot = air.material(293.15, 450.0).unwrap();
        assert!(close(hot.fluid.unwrap().newtonian().unwrap(), mu(400.0)));
        assert_eq!(air.outside(450.0), vec![Property::Density, Property::Conductivity, Property::Viscosity]);
        // Solids: a curve gives the constant, the other constants stay.
        let steel = entry("EN10025-S355").unwrap().material(293.15, CELSIUS_ZERO + 100.0).unwrap();
        let thermal = steel.thermal.unwrap();
        assert!(close(thermal.conductivity.unwrap().axes()[0], 54.0 - 3.33));
        assert_eq!(steel.density, Some(7850.0));
        assert!(steel.temperature_dependent.as_ref().unwrap().conducts());
        // Reports carry every curve in display units with its citation.
        let mm = DisplaySystem::new(Preset::MmTonneS).with_override(Quantity::Temperature, "degC").unwrap();
        let report = entry("air").unwrap().report_curves(&mm).unwrap();
        let mu = report.iter().find(|c| c.property == Property::Viscosity).unwrap();
        assert!(close(mu.points[0][0], 200.0 - CELSIUS_ZERO));
        assert!(close(mm.to_si(Quantity::DynamicViscosity, mu.points[0][1]).unwrap(), curve("air", Property::Viscosity).points.0[0][1]));
        assert_eq!(mu.source.page, "19");
    }

    #[test]
    fn a_curve_beside_its_constant_or_an_implausible_point_is_refused() {
        let original = entries().unwrap().to_vec();
        let air = original.iter().position(|e| e.id == "air").unwrap();
        let mut bad = original.clone();
        let datum = bad[0].properties[&Property::Density].clone();
        bad[air].properties.insert(Property::Density, datum);
        assert!(validate(&bad).is_err());
        let mut bad = original.clone();
        bad[air].curves.get_mut(&Property::Viscosity).unwrap().points.0[3][1] = 1e9;
        assert!(validate(&bad).is_err());
        let mut bad = original.clone();
        bad[air].curves.get_mut(&Property::Viscosity).unwrap().points.0.swap(2, 3);
        assert!(validate(&bad).is_err());
        let mut bad = original.clone();
        bad[air].curves.get_mut(&Property::Viscosity).unwrap().source.document.clear();
        assert!(validate(&bad).is_err());
        let mut bad = original.clone();
        bad[air].curves.remove(&Property::Viscosity);
        assert!(validate(&bad).is_err(), "a fluid needs its viscosity");
        // A curve no solver follows (expansion here) would be taken as its constant alone.
        let mut bad = original.clone();
        let steel = bad.iter().position(|e| e.id == "EN10025-S355").unwrap();
        let alpha = bad[steel].properties.remove(&Property::Expansion).unwrap();
        let mut curve = bad[steel].curves[&Property::Conductivity].clone();
        curve.unit = alpha.unit.clone();
        curve.points = TemperatureTable(vec![[293.15, alpha.value], [473.15, alpha.value]]);
        bad[steel].curves.insert(Property::Expansion, curve);
        assert!(validate(&bad).is_err(), "an expansion curve");
        // A gas's density is plausible for a fluid only, not for a metal.
        let mut bad = original;
        bad[0].properties.get_mut(&Property::Density).unwrap().value = 1.2;
        assert!(validate(&bad).is_err(), "a steel as light as air");
    }
}
