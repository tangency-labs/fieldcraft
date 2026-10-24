//! The fixed, independently curated material catalogue. No runtime import or network lookup.

use fieldcraft_model::library::{self, Entry, Property};
use fieldcraft_model::units::{CELSIUS_ZERO, Quantity};
use serde_json::{Value, json};

use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const LIST: CommandSpec = CommandSpec {
    id: "model.material.library.list",
    label: "Browse Material Library",
    params: r#"{"family"?: "steel|stainless|aluminium|titanium"} — list the fixed cited catalogue, with properties in display units, conditions and unavailable properties. Read-only. Strengths are informational, never selected plasticity or failure criteria."#,
    form: Some(super::forms::MATERIAL_LIBRARY_LIST),
    kind: CommandKind::Model,
    requires: &[],
    run: list,
};

pub const ADD: CommandSpec = CommandSpec {
    id: "model.material.library",
    label: "Add Library Material",
    params: r#"{"entry": "EN10025-S235|EN10025-S355|EN10088-14301|EN10088-14401|AA6061-T6|AA6082-T6|AA7075-T6|AA2024-T3|Ti6Al4V-annealed", "name"?: "Steel" (default: entry id), "reference_temperature"?: 20 (default 20 °C; bare numbers in display units, or a value with its unit)} — create an isotropic material from cited constant properties, refusing duplicate names. The stress-free temperature is your modelling choice. Missing data stays absent; strength data is reported but never applied as plasticity/damage/failure criteria."#,
    form: Some(super::forms::MATERIAL_LIBRARY),
    kind: CommandKind::Model,
    requires: &[],
    run: add,
};

fn parameters<'a>(p: &'a Value, allowed: &[&str]) -> Result<&'a serde_json::Map<String, Value>> {
    let obj = p.as_object().ok_or_else(|| EngineError::Params("material library parameters must be an object".into()))?;
    if obj.keys().any(|k| !allowed.contains(&k.as_str())) {
        return Err(EngineError::Params("unknown material library parameter; use the command's documented fields".into()));
    }
    Ok(obj)
}

fn text<'a>(p: &'a serde_json::Map<String, Value>, key: &str) -> Result<Option<&'a str>> {
    p.get(key)
        .map(|v| {
            v.as_str()
                .filter(|s| !s.is_empty() && s.len() <= 64)
                .ok_or_else(|| EngineError::invalid(key, "must be a nonempty string of at most 64 bytes"))
        })
        .transpose()
}

fn report(s: &Session, e: &Entry) -> Result<Value> {
    let (applied, information): (Vec<_>, Vec<_>) = e.report(s.units())?.into_iter().partition(|v| v.property.applied());
    let unavailable: Vec<_> = [
        Property::Young,
        Property::Poisson,
        Property::Density,
        Property::Conductivity,
        Property::SpecificHeat,
        Property::Expansion,
        Property::Yield,
        Property::Tensile,
    ]
    .into_iter()
    .filter(|p| !e.properties.contains_key(p))
    .collect();
    Ok(json!({
        "id": e.id, "name": e.name, "family": e.family, "designation": e.designation,
        "display_units": s.units().label(), "values": applied,
        "not_applied": information, "unavailable": unavailable,
        "limits": "Constant isotropic approximation within each source's stated conditions. Strengths are informational; no plasticity, damage, safety factor or failure criterion is selected. Missing properties must be supplied before analyses that need them."
    }))
}

fn list(s: &mut Session, p: &Value) -> Result<Value> {
    let p = parameters(p, &["family"])?;
    let family = text(p, "family")?;
    if family.is_some_and(|f| !library::FAMILIES.iter().any(|(id, _)| *id == f)) {
        return Err(EngineError::invalid("family", "choose steel, stainless, aluminium or titanium"));
    }
    let entries = library::entries()?.iter().filter(|e| family.is_none_or(|f| e.family == f)).map(|e| report(s, e)).collect::<Result<Vec<_>>>()?;
    Ok(json!({"count": entries.len(), "entries": entries}))
}

fn add(s: &mut Session, p: &Value) -> Result<Value> {
    let p = parameters(p, &["entry", "name", "reference_temperature"])?;
    let id = text(p, "entry")?.ok_or_else(|| EngineError::param("entry", crate::schema::ProblemKind::Required))?;
    let entry = library::entry(id)?;
    let name = if p.get("name").and_then(Value::as_str) == Some("") { &entry.id } else { text(p, "name")?.unwrap_or(&entry.id) };
    let reference = match p.get("reference_temperature") {
        Some(v) => {
            v.as_f64().filter(|v| v.is_finite()).ok_or_else(|| EngineError::invalid("reference_temperature", "must be a finite temperature"))?
        }
        None => CELSIUS_ZERO + 20.0, // A modelling default of 20 degrees Celsius, not a conversion factor.
    };
    let material = entry.material(reference)?;
    let mut reply = report(s, entry)?;
    reply["name"] = json!(name);
    reply["material_si"] = json!(material);
    reply["reference_temperature"] = json!({
        "value": s.units().express(Quantity::Temperature, reference).map_err(|e| EngineError::Params(e.to_string()))?,
        "unit": s.units().symbol(Quantity::Temperature),
        "basis": "user modelling choice: temperature of zero thermal strain; default 20 degrees Celsius",
        "applied": material.thermal.as_ref().is_some_and(|t| t.expansion.is_some())
    });
    s.model.add_material(name, material)?;
    Ok(reply)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::units::{DisplaySystem, Preset};

    #[test]
    fn list_is_read_only_and_filters_cited_conditions() {
        let mut s = Session::default();
        let before = s.revision();
        let all = s.execute(LIST.id, &json!({})).unwrap();
        assert_eq!(all["count"], 9);
        let aluminium = s.execute(LIST.id, &json!({"family": "aluminium"})).unwrap();
        assert_eq!(aluminium["count"], 4);
        assert_eq!(s.revision(), before);
        for e in aluminium["entries"].as_array().unwrap() {
            assert_eq!(e["family"], "aluminium");
            assert!(e["designation"]["condition"].as_str().unwrap().len() > 20);
            assert!(e["values"][0]["source"]["url"].as_str().unwrap().starts_with("https://"));
        }
    }

    #[test]
    fn every_entry_is_the_same_si_material_in_all_presets_and_temperature_scales() {
        for preset in Preset::ALL {
            for temperature in ["K", "degC", "degF"] {
                let mut s = Session {
                    default_units: DisplaySystem::new(preset).with_override(Quantity::Temperature, temperature).unwrap(),
                    ..Session::default()
                };
                for &(id, _) in library::CHOICES {
                    let reply = s.execute(ADD.id, &json!({"entry": id, "reference_temperature": "68 degF"})).unwrap();
                    let material = s.model.material(id).unwrap();
                    assert_eq!(*material, library::entry(id).unwrap().material(293.15).unwrap());
                    assert!(material.plastic.is_none() && material.damage.is_none() && material.erosion.is_none());
                    assert_eq!(reply["name"], id);
                    assert_eq!(reply["material_si"], json!(material));
                    for value in reply["values"].as_array().unwrap() {
                        let property: Property = serde_json::from_value(value["property"].clone()).unwrap();
                        let shown = value["value"].as_f64().unwrap();
                        let si = s.units().to_si(property.quantity(), shown).unwrap();
                        let expected = library::entry(id).unwrap().properties[&property].value;
                        assert!((si - expected).abs() < 1e-12 * expected.abs());
                    }
                }
            }
        }
    }

    #[test]
    fn default_reference_and_optional_data_are_honest_and_failures_atomic() {
        let mut s = Session::default();
        let r = s.execute(ADD.id, &json!({"entry": "AA2024-T3", "name": "Skin"})).unwrap();
        assert!(r["unavailable"].as_array().unwrap().contains(&json!("conductivity")));
        assert_eq!(s.model.material("Skin").unwrap().thermal, None);
        s.execute(ADD.id, &json!({"entry": "EN10025-S355"})).unwrap();
        assert_eq!(s.model.material("EN10025-S355").unwrap().thermal.unwrap().reference_temperature, 293.15);
        let before = s.model.clone();
        for bad in [
            json!({"entry": "AA2024-T3", "name": "Skin"}),
            json!({"entry": "no-such-entry"}),
            json!({"entry": 3}),
            json!({"entry": null}),
            json!({"entry": "x".repeat(100_000)}),
            json!({"entry": "AA2024-T3", "name": "invalid name"}),
            json!({"entry": "AA2024-T3", "name": "x".repeat(100_000)}),
            json!({"entry": "AA2024-T3", "reference_temperature": "-1 K"}),
            json!({"entry": "AA2024-T3", "density": 1}),
            json!([]),
        ] {
            assert!(s.execute(ADD.id, &bad).is_err());
            assert_eq!(s.model, before);
        }
        for bad in [json!({"family": []}), json!({"family": "x".repeat(100_000)}), json!({"family": "plastic"})] {
            assert!(s.execute(LIST.id, &bad).is_err());
            assert_eq!(s.model, before);
        }
    }

    #[test]
    fn library_steel_preserves_the_cantilever_solution() {
        let original = crate::sample::cantilever();
        let mut baseline = Session::default();
        baseline.run_script(&original).unwrap();
        let expected = baseline.execute("solve.static", &json!({})).unwrap();
        let mut replacement = original;
        let steps = replacement.as_array_mut().unwrap();
        let steel = steps.iter_mut().find(|s| s["command"] == "model.material.create").unwrap();
        *steel = json!({"command": ADD.id, "params": {"entry": "EN10025-S355", "name": "Steel"}});
        let mut s = Session::default();
        s.run_script(&replacement).unwrap();
        let actual = s.execute("solve.static", &json!({})).unwrap();
        assert_eq!(actual["max_displacement"], expected["max_displacement"]);
        assert!(s.model.check().iter().all(|d| d.severity != fieldcraft_model::Severity::Warning));
    }

    #[test]
    fn cited_stainless_conductivity_gives_the_linear_temperature_patch() {
        let recipe: Value = serde_json::from_str(include_str!("../../../../docs/materials/conduction.json")).unwrap();
        let mut s = Session::default();
        s.run_script(&recipe).unwrap();
        let probe = s.execute("solve.probe", &json!({"step": "Heat", "at": [[0.01, 0.01, 0.06]]})).unwrap();
        let node = &probe["nodes"][0];
        assert!((node["temperature"].as_f64().unwrap() - 298.15).abs() < 1e-10, "{probe}");
        assert!((node["heat_flux"][2].as_f64().unwrap() + 1250.0).abs() < 1e-8, "{probe}");
        let cold = s.execute("solve.probe", &json!({"step": "Heat", "set": "Cold"})).unwrap();
        let hot = s.execute("solve.probe", &json!({"step": "Heat", "set": "Hot"})).unwrap();
        let sum = |v: &Value| v["nodes"].as_array().unwrap().iter().map(|n| n["reaction"].as_f64().unwrap()).sum::<f64>();
        assert!((sum(&cold).abs() - 0.75).abs() < 1e-10, "{cold}");
        assert!((sum(&cold) + sum(&hot)).abs() < 1e-10);
    }
}
