use fieldcraft_model::Severity;
use fieldcraft_model::units::{DisplaySystem, Preset, Quantity};
use serde_json::{Value, json};

use super::{listed, name};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const INFO: CommandSpec =
    CommandSpec { id: "model.info", label: "Model Info", params: "{}", form: Some(&[]), kind: CommandKind::Model, requires: &[], run: info };

pub const UNITS: CommandSpec = CommandSpec {
    id: "model.units",
    label: "Unit System",
    params: r#"{"system": "mm-t-s-N-MPa|m-kg-s-N-Pa|mm-kg-ms-kN-GPa|mm-g-ms-N-MPa|in-lbf-s-psi|default", "temperature"?: "°C|K|°F" (temperature differences and per-degree units follow it), "overrides": {"energy": "J", "angle": "deg", …}} — the units this model is shown and typed in (bare numbers of later commands are in them); `default` follows the user's setting. Stored values never change: every quantity is stored in SI. `unspecified` (older models) reads as mm-t-s-N-MPa"#,
    form: Some(super::forms::UNITS),
    kind: CommandKind::Model,
    requires: &[],
    run: units,
};

pub const CHECK: CommandSpec =
    CommandSpec { id: "model.check", label: "Check Model", params: "{}", form: Some(&[]), kind: CommandKind::Model, requires: &[], run: check };

pub const SET_RESOLVE: CommandSpec = CommandSpec {
    id: "model.set.resolve",
    label: "Resolve Set",
    params: r#"{"name": "root"} — the set's nodes, elements and element faces on the current mesh (the first 100 of each); fails naming every member that is not on the mesh"#,
    form: Some(super::forms::SET_RESOLVE),
    kind: CommandKind::Model,
    requires: &[Requirement::Set],
    run: set_resolve,
};

fn info(s: &mut Session, _: &Value) -> Result<Value> {
    let counts: serde_json::Map<String, Value> = s.model.element_counts().into_iter().map(|(k, n)| (k.name().to_string(), json!(n))).collect();
    let definitions: serde_json::Map<String, Value> = s.model.definition_counts().into_iter().map(|(k, n)| (k.to_string(), json!(n))).collect();
    Ok(json!({
        "nodes": s.model.nodes().len(),
        "elements": counts,
        "display_units": s.units().label(),
        "own_units": s.model.display_units().is_some(),
        "definitions": definitions,
    }))
}

/// The units a display system shows the common quantities in, for replies.
pub(crate) fn unit_table(d: &DisplaySystem) -> Value {
    let q = [
        Quantity::Length,
        Quantity::Mass,
        Quantity::Time,
        Quantity::Force,
        Quantity::Pressure,
        Quantity::Density,
        Quantity::Energy,
        Quantity::Temperature,
        Quantity::TemperatureDifference,
        Quantity::ThermalExpansion,
        Quantity::Angle,
    ];
    Value::Object(q.iter().map(|q| (q.name().to_string(), json!(d.symbol(*q)))).collect())
}

fn units(s: &mut Session, p: &Value) -> Result<Value> {
    let label = p.get("system").and_then(Value::as_str).unwrap_or_default();
    let base = match label {
        "default" => None,
        _ => Some(Preset::from_label(label).ok_or_else(|| {
            let known: Vec<&str> = Preset::ALL.iter().map(|u| u.label()).collect();
            EngineError::Params(format!("`system` must be one of {} or default", known.join(", ")))
        })?),
    };
    let overrides = p.get("overrides").and_then(Value::as_object);
    if base.is_none() && overrides.is_some_and(|o| !o.is_empty()) {
        return Err(EngineError::invalid("overrides", "need a `system` they change"));
    }
    let mut system = base.map(DisplaySystem::new);
    if let (Some(d), Some(o)) = (system.as_mut(), overrides) {
        for (name, unit) in o {
            let q = Quantity::from_name(name).ok_or_else(|| EngineError::invalid(format!("overrides/{name}"), "is not a quantity"))?;
            let unit = unit.as_str().ok_or_else(|| EngineError::invalid(format!("overrides/{name}"), "must be a unit, such as \"J\""))?;
            *d = d.clone().with_override(q, unit).map_err(|e| EngineError::invalid(format!("overrides/{name}"), e.to_string()))?;
        }
    }
    // A temperature scale of its own (the form's choice): the temperature differences and the
    // per-degree quantities follow it (`DisplaySystem`), unless overridden too.
    if let Some(t) = p.get("temperature").filter(|t| !t.is_null()) {
        let t = t.as_str().filter(|t| ["°C", "K", "°F"].contains(t)).ok_or_else(|| EngineError::invalid("temperature", "must be °C, K or °F"))?;
        if overrides.is_some_and(|o| o.contains_key("temperature")) {
            return Err(EngineError::invalid("temperature", "is given twice (also in overrides)"));
        }
        let d = system.as_mut().ok_or_else(|| EngineError::invalid("temperature", "needs a `system` it changes"))?;
        *d = d.clone().with_override(Quantity::Temperature, t).map_err(|e| EngineError::invalid("temperature", e.to_string()))?;
    }
    s.model.set_display_units(system);
    Ok(json!({"system": s.units().label(), "own_units": s.model.display_units().is_some(), "units": unit_table(s.units())}))
}

fn check(s: &mut Session, _: &Value) -> Result<Value> {
    let mut diagnostics = s.model.check();
    if let Some(design) = s.design.as_deref() {
        for (body, n) in super::replace::stale(&s.model, design) {
            diagnostics.push(fieldcraft_model::Diagnostic {
                severity: Severity::Warning,
                subject: format!("body `{body}`"),
                message: format!(
                    "{n} elements mesh a body the open design no longer has (its partitions changed since it was meshed); mesh.clear removes them"
                ),
            });
        }
    }
    if let Some(why) = s.mesh_outdated() {
        diagnostics.push(fieldcraft_model::Diagnostic { severity: Severity::Warning, subject: "mesh".into(), message: why });
    }
    let errors = diagnostics.iter().filter(|d| d.severity == Severity::Error).count();
    let first_error = diagnostics.iter().find(|d| d.severity == Severity::Error).map(ToString::to_string);
    let progress = crate::progress::to_json(&s.progress_given(first_error));
    Ok(json!({
        "ready": errors == 0,
        "errors": errors,
        "warnings": diagnostics.len() - errors,
        "diagnostics": diagnostics,
        "progress": progress,
    }))
}

fn set_resolve(s: &mut Session, p: &Value) -> Result<Value> {
    let r = s.model.resolve_set(name(p)?)?;
    Ok(json!({"nodes": listed(&r.nodes), "elements": listed(&r.elements), "faces": listed(&r.faces)}))
}

#[cfg(test)]
mod temperature_tests {
    use serde_json::json;

    use crate::Session;

    /// A temperature scale of its own carries to ΔT and the per-degree units; numbers typed
    /// after follow it (a 1.8 °F rise is 1 K).
    #[test]
    fn a_temperature_scale_carries_to_differences_and_per_degree_units() {
        let mut s = Session::default();
        let r = s.execute("model.units", &json!({"system": "mm-t-s-N-MPa", "temperature": "°F"})).unwrap();
        assert_eq!(
            (&r["units"]["temperature"], &r["units"]["temperature_difference"], &r["units"]["thermal_expansion"], &r["units"]["pressure"]),
            (&json!("°F"), &json!("°F"), &json!("°F⁻¹"), &json!("MPa")),
            "{r}"
        );
        s.execute("model.material.create", &json!({"name": "Steel", "thermal": {"expansion": 1e-5, "reference_temperature": 68}})).unwrap();
        let t = s.model.material("Steel").unwrap().thermal.unwrap();
        assert_eq!(t.reference_temperature, 293.15);
        assert!((t.expansion.unwrap().axes()[0] - 1.8e-5).abs() < 1e-20, "{t:?}");
        for bad in [
            json!({"system": "mm-t-s-N-MPa", "temperature": "°R"}),
            json!({"system": "mm-t-s-N-MPa", "temperature": "°F", "overrides": {"temperature": "K"}}),
        ] {
            assert!(s.execute("model.units", &bad).is_err(), "{bad}");
        }
    }
}
