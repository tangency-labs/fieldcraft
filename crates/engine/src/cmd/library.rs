//! The fixed, independently curated material catalogue. No runtime import or network lookup.

use fieldcraft_model::library::{self, Entry, Property};
use fieldcraft_model::units::Quantity;
use serde_json::{Value, json};

use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const LIST: CommandSpec = CommandSpec {
    id: "model.material.library.list",
    label: "Browse Material Library",
    params: r#"{"family"?: "steel|stainless|aluminium|titanium|fluid"} — list the fixed cited catalogue, with properties in display units, temperature curves (cited points, linear between them), conditions and unavailable properties. Read-only. Strengths are informational, never selected plasticity or failure criteria."#,
    form: Some(super::forms::MATERIAL_LIBRARY_LIST),
    kind: CommandKind::Model,
    requires: &[],
    run: list,
};

pub const ADD: CommandSpec = CommandSpec {
    id: "model.material.library.add",
    label: "Add Library Material",
    params: r#"{"entry": "EN10025-S235|EN10025-S355|EN10088-14301|EN10088-14401|AA6061-T6|AA6082-T6|AA7075-T6|AA2024-T3|Ti6Al4V-annealed|air", "name"?: "Steel" (default: entry id), "reference_temperature"?: 20 (default 20 °C; bare numbers in display units, or a value with its unit), "temperature"?: 20 (default 20 °C: where temperature-dependent properties give their constants)} — create an isotropic material (a Newtonian fluid for air) from cited properties, refusing duplicate names. Temperature curves are saved with the material; the solvers use the constants taken at `temperature`, with a warning when it lies outside a curve (the end value is taken). The stress-free temperature is your modelling choice. Missing data stays absent; strength data is reported but never applied as plasticity/damage/failure criteria."#,
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
    // A fluid has no elastic constants, expansion or strengths to miss; a solid no viscosity.
    let relevant = |p: &Property| {
        if e.is_fluid() {
            matches!(p, Property::Density | Property::Conductivity | Property::SpecificHeat | Property::Viscosity)
        } else {
            *p != Property::Viscosity
        }
    };
    let unavailable: Vec<_> =
        Property::ALL.into_iter().filter(relevant).filter(|p| !e.properties.contains_key(p) && !e.curves.contains_key(p)).collect();
    Ok(json!({
        "id": e.id, "name": e.name, "family": e.family, "designation": e.designation,
        "display_units": s.units().label(), "values": applied, "curves": e.report_curves(s.units())?,
        "not_applied": information, "unavailable": unavailable,
        "limits": "Isotropic properties within each source's stated conditions; curves are linear between their points and held at their end values outside them. Strengths are informational; no plasticity, damage, safety factor or failure criterion is selected. Missing properties must be supplied before analyses that need them."
    }))
}

fn list(s: &mut Session, p: &Value) -> Result<Value> {
    let p = parameters(p, &["family"])?;
    let family = text(p, "family")?;
    if family.is_some_and(|f| !library::FAMILIES.iter().any(|(id, _)| *id == f)) {
        return Err(EngineError::invalid("family", "choose steel, stainless, aluminium, titanium or fluid"));
    }
    let entries = library::entries()?.iter().filter(|e| family.is_none_or(|f| e.family == f)).map(|e| report(s, e)).collect::<Result<Vec<_>>>()?;
    Ok(json!({"count": entries.len(), "entries": entries}))
}

fn add(s: &mut Session, p: &Value) -> Result<Value> {
    let p = parameters(p, &["entry", "name", "reference_temperature", "temperature"])?;
    let id = text(p, "entry")?.ok_or_else(|| EngineError::param("entry", crate::schema::ProblemKind::Required))?;
    let entry = library::entry(id)?;
    let name = if p.get("name").and_then(Value::as_str) == Some("") { &entry.id } else { text(p, "name")?.unwrap_or(&entry.id) };
    let temperature = |key: &str| match p.get(key) {
        Some(v) => v.as_f64().filter(|v| v.is_finite()).ok_or_else(|| EngineError::invalid(key, "must be a finite temperature")),
        None => Ok(library::DEFAULT_TEMPERATURE),
    };
    let (reference, at) = (temperature("reference_temperature")?, temperature("temperature")?);
    let material = entry.material(reference, at)?;
    let shown = |t: f64| s.units().express(Quantity::Temperature, t).map_err(|e| EngineError::Params(e.to_string()));
    let mut reply = report(s, entry)?;
    reply["name"] = json!(name);
    reply["material_si"] = json!(material);
    reply["reference_temperature"] = json!({
        "value": shown(reference)?,
        "unit": s.units().symbol(Quantity::Temperature),
        "basis": "user modelling choice: temperature of zero thermal strain; default 20 degrees Celsius",
        "applied": material.thermal.as_ref().is_some_and(|t| t.expansion.is_some())
    });
    reply["temperature"] = json!({
        "value": shown(at)?,
        "unit": s.units().symbol(Quantity::Temperature),
        "basis": "user modelling choice: the temperature the constants of temperature-dependent properties are taken at; default 20 degrees Celsius",
        "applied": !entry.curves.is_empty()
    });
    let outside = entry.outside(at);
    if !outside.is_empty() {
        reply["warnings"] = json!([format!(
            "{} lies outside the {} curve(s) of {}: their end values were taken as the constants, an extrapolation the source does not cover",
            Quantity::Temperature.format(at, s.units()).map_err(|e| EngineError::Params(e.to_string()))?,
            outside.iter().map(|p| p.name()).collect::<Vec<_>>().join(", "),
            entry.id
        )]);
    }
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
        assert_eq!(all["count"], 10);
        let aluminium = s.execute(LIST.id, &json!({"family": "aluminium"})).unwrap();
        assert_eq!(aluminium["count"], 4);
        assert_eq!(s.revision(), before);
        for e in aluminium["entries"].as_array().unwrap() {
            assert_eq!(e["family"], "aluminium");
            assert!(e["designation"]["condition"].as_str().unwrap().len() > 20);
            for value in e["values"].as_array().unwrap() {
                let source = &value["source"];
                for key in ["document", "edition", "location"] {
                    assert!(!source[key].as_str().unwrap().is_empty());
                }
                if let Some(page) = source.get("page") {
                    assert!(!page.as_str().unwrap().is_empty());
                }
                if let Some(url) = source.get("url") {
                    assert!(url.as_str().unwrap().starts_with("https://"));
                }
            }
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
                    assert_eq!(*material, library::entry(id).unwrap().material(293.15, 293.15).unwrap());
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
    fn fluids_carry_their_cited_curves_and_warn_outside_them() {
        let mut s = Session::default();
        let r = s.execute(ADD.id, &json!({"entry": "air", "temperature": "30 degC"})).unwrap();
        assert!(r.get("warnings").is_none(), "{r}");
        assert_eq!(r["unavailable"], json!([]));
        let air = s.model.material("air").unwrap();
        let mu = air.fluid.unwrap().newtonian().unwrap();
        // 303.15 K, between the 300 K and 310 K points of the Standard's eq. (51).
        let sutherland = |t: f64| 1.458e-6 * t.powf(1.5) / (t + 110.4);
        let expected = sutherland(300.0) + (303.15 - 300.0) / 10.0 * (sutherland(310.0) - sutherland(300.0));
        assert!((mu - expected).abs() < 1e-15, "{mu}");
        assert_eq!(air.temperature_dependent.as_ref().unwrap().viscosity.as_ref().unwrap().0.len(), 21);
        let curves = r["curves"].as_array().unwrap();
        assert_eq!(curves.len(), 3);
        assert!(curves.iter().all(|c| c["source"]["document"].as_str().unwrap().contains("Standard Atmosphere")));
        // The default display system shows °C: the first point at 200 K is −73.15 °C.
        assert!((curves[0]["points"][0][0].as_f64().unwrap() + 73.15).abs() < 1e-9, "{}", curves[0]);

        let r = s.execute(ADD.id, &json!({"entry": "air", "name": "Hot_air", "temperature": "500 K"})).unwrap();
        let warning = r["warnings"][0].as_str().unwrap();
        assert!(warning.contains("outside") && warning.contains("viscosity") && warning.contains("air"), "{warning}");
        let air = s.model.material("Hot_air").unwrap();
        assert_eq!(air.density, Some(library::entry("air").unwrap().curves[&Property::Density].points.0.last().unwrap()[1]));
        assert!(air.elastic.is_none());
        let fluids = s.execute(LIST.id, &json!({"family": "fluid"})).unwrap();
        assert_eq!(fluids["count"], 1);
        let before = s.model.clone();
        for bad in [json!({"entry": "air", "name": "A2", "temperature": "-1 K"}), json!({"entry": "air", "name": "A3", "temperature": "hot"})] {
            assert!(s.execute(ADD.id, &bad).is_err());
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
        let recipe: Value = serde_json::from_str(include_str!("../../../../docs/recipes/material-library-conduction.json")).unwrap();
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
