//! A small, cited catalogue of constant metal properties. Values are stored in SI, like the
//! model; display units are only used when reporting them. Strengths are information, never an
//! automatically selected plasticity law or failure criterion. Applicability belongs to each
//! value: a room-temperature constant does not describe temperature-dependent behaviour.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use crate::units::{DisplaySystem, Quantity};
use crate::{Directional, Elastic, Error, Material, Result, Thermal};

const MIN_YOUNG: f64 = 1e9; // unit-factors: SI plausibility lower bound in Pa, not a conversion
const MAX_YOUNG: f64 = 1e12; // unit-factors: SI plausibility upper bound in Pa, not a conversion
const MAX_EXPANSION: f64 = 1e-3; // unit-factors: SI plausibility upper bound in inverse kelvin, not a conversion
const MIN_STRENGTH: f64 = 1e6; // unit-factors: SI plausibility lower bound in Pa, not a conversion

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
];

pub const FAMILIES: &[(&str, &str)] =
    &[("steel", "Structural steel"), ("stainless", "Stainless steel"), ("aluminium", "Aluminium"), ("titanium", "Titanium")];

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
        }
    }

    pub fn applied(self) -> bool {
        !matches!(self, Self::Yield | Self::Tensile)
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
        for (p, d) in &e.properties {
            let q = p.quantity();
            if d.unit != q.si_symbol()
                || !d.value.is_finite()
                || d.applicability.trim().is_empty()
                || [&d.source.document, &d.source.edition, &d.source.location, &d.source.page].iter().any(|s| s.trim().is_empty())
                || d.source.url.as_ref().is_some_and(|url| !url.starts_with("https://"))
            {
                return Err(Error::Invalid(format!("material library entry {} has an uncited or non-SI property", e.id)));
            }
            let plausible = match p {
                Property::Young => (MIN_YOUNG..=MAX_YOUNG).contains(&d.value),
                Property::Poisson => (0.0..0.5).contains(&d.value),
                Property::Density => (100.0..=30_000.0).contains(&d.value),
                Property::Conductivity => (0.01..=2000.0).contains(&d.value),
                Property::SpecificHeat => (1.0..=10_000.0).contains(&d.value),
                Property::Expansion => (0.0..=MAX_EXPANSION).contains(&d.value),
                Property::Yield | Property::Tensile => (MIN_STRENGTH..=1e10).contains(&d.value),
            };
            if !plausible {
                return Err(Error::Invalid(format!("material library entry {} has an implausible property", e.id)));
            }
        }
        e.material(crate::units::CELSIUS_ZERO)?;
    }
    Ok(())
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
    fn value(&self, property: Property) -> Option<f64> {
        self.properties.get(&property).map(|d| d.value)
    }

    /// Build an SI material. The stress-free reference temperature is a modelling choice made
    /// by the caller, not an asserted material property from a handbook.
    pub fn material(&self, reference_temperature: f64) -> Result<Material> {
        if !reference_temperature.is_finite() || reference_temperature < 0.0 {
            return Err(Error::Invalid("the stress-free reference temperature must be at or above zero kelvin".into()));
        }
        let (young, poisson) = (self.value(Property::Young), self.value(Property::Poisson));
        let elastic = match (young, poisson) {
            (Some(young), Some(poisson)) => Some(Elastic::Isotropic { young, poisson }),
            _ => return Err(Error::Invalid(format!("material library entry {} needs both isotropic elastic constants", self.id))),
        };
        let conductivity = self.value(Property::Conductivity).map(Directional::Isotropic);
        let specific_heat = self.value(Property::SpecificHeat);
        let expansion = self.value(Property::Expansion).map(Directional::Isotropic);
        let thermal = (conductivity.is_some() || specific_heat.is_some() || expansion.is_some()).then_some(Thermal {
            conductivity,
            specific_heat,
            expansion,
            reference_temperature,
        });
        let material = Material { elastic, density: self.value(Property::Density), thermal, ..Material::default() };
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
    use crate::units::{Preset, Unit};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1e-20)
    }

    #[test]
    fn catalogue_is_cited_and_materials_remain_si_in_every_display_system() {
        let all = entries().unwrap();
        assert_eq!(all.len(), 9);
        validate(all).unwrap();
        for e in all {
            let material = e.material(293.15).unwrap();
            assert!(material.plastic.is_none() && material.damage.is_none() && material.erosion.is_none());
            for preset in Preset::ALL {
                let display = DisplaySystem::new(preset);
                for r in e.report(&display).unwrap() {
                    let datum = e.properties.get(&r.property).unwrap();
                    assert!(close(display.to_si(r.property.quantity(), r.value).unwrap(), datum.value), "{}: {:?}", e.id, r.property);
                    assert_eq!((&r.source, r.basis, &r.applicability), (&datum.source, datum.basis, &datum.applicability));
                }
                assert_eq!(material, e.material(293.15).unwrap());
                material.validate().unwrap();
            }
        }
    }

    #[test]
    fn a_missing_citation_invalid_property_or_duplicate_id_is_refused() {
        let original = entries().unwrap().to_vec();
        let mut bad = original.clone();
        bad.first_mut().unwrap().properties.get_mut(&Property::Young).unwrap().source.page.clear();
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
        assert!(stainless.material(-1.0).is_err());
        assert!(stainless.material(f64::INFINITY).is_err());
    }
}
