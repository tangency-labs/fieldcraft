//! Values as people read them: a definition's fields with the labels and units of its form
//! ([`rows`], the context panel's properties), and numbers rounded and grouped ([`quantity`]).
//! The command's form schema ([`fieldcraft_engine::schema`]) knows each field's label, unit and
//! choices, so a definition is described by the same words as the form that edits it, never as
//! its JSON.

use fieldcraft_engine::schema::{self, Field, Kind};
use fieldcraft_model::{DisplaySystem, Model};
use serde_json::{Map, Value};

use crate::browser::Def;

/// Lists longer than this are given by their length.
const LIST: usize = 6;
/// Nesting followed at most.
const MAX_DEPTH: usize = 6;

/// Whole numbers below this are written in full.
const IN_FULL_BELOW: f64 = 1e9; // unit-factors: a count of digits, not a unit

/// A number in a few significant figures: whole numbers below a billion in full, with thin spaces
/// between groups of three digits ("210 000"); others as [`crate::results::number`].
pub fn quantity(v: f64) -> String {
    if v.is_finite() && v.fract() == 0.0 && v.abs() < IN_FULL_BELOW {
        let digits = format!("{}", v.abs() as u64);
        let mut out = String::new();
        for (i, c) in digits.chars().enumerate() {
            if i > 0 && (digits.len() - i) % 3 == 0 {
                out.push('\u{2009}');
            }
            out.push(c);
        }
        return if v < 0.0 { format!("−{out}") } else { out };
    }
    // Trailing zeros say nothing here: "0.3", "7.85e−9".
    let s = crate::results::number(v);
    let (mantissa, exponent) = s.split_once('e').map_or((s.as_str(), None), |(m, e)| (m, Some(e)));
    let mantissa = if mantissa.contains('.') { mantissa.trim_end_matches('0').trim_end_matches('.') } else { mantissa };
    let out = exponent.map_or_else(|| mantissa.to_string(), |e| format!("{mantissa}e{e}"));
    out.replace('-', "−")
}

/// `v` with its unit, when it has one.
fn with_unit(v: f64, unit: Option<&str>) -> String {
    match unit {
        Some(u) => format!("{} {u}", quantity(v)),
        None => quantity(v),
    }
}

/// The definition `name` of `kind` as (label, text) rows, in its form's words and in `units`;
/// `None` when it is gone (or its numbers cannot be told apart).
pub fn rows(model: &Model, kind: Def, name: &str, units: &DisplaySystem) -> Option<Vec<(String, String)>> {
    let v = crate::browser::definition(model, kind, name)?;
    if kind == Def::Set {
        return Some(set_rows(model, name, v.as_object()?));
    }
    // A datum reads as its edit command's parameters (its geometry in millimetres).
    if kind == Def::Datum {
        let (command, params) = kind.edit(model, name, units)?;
        return Some(command_rows(&command, params.as_object()?, units));
    }
    let shown = fieldcraft_engine::units::shown_params(&kind.command("create"), &v, units).ok()?;
    let obj = shown.as_object()?;
    Some(command_rows(&kind.command("create"), obj, units))
}

/// Parameters of `command`, given in `units`, as (label, text) rows in its form's words and units;
/// fields its form does not show stay readable by their key.
pub fn command_rows(command: &str, obj: &Map<String, Value>, units: &DisplaySystem) -> Vec<(String, String)> {
    let form = fieldcraft_engine::find_command(command).and_then(|spec| spec.form);
    let mut out = Vec::new();
    let mut used = Vec::new();
    if let Some(fields) = form {
        object(fields, obj, units, &mut out, &mut used, 0);
    }
    for (k, field) in obj {
        if k != "name" && !used.iter().any(|u| u == k) {
            out.push((sentence(k), crate::context::value_text(field)));
        }
    }
    out
}

/// The seed on `edge` as (label, text) rows: the edge, then `mesh.seed.edge`'s fields in `units`.
pub fn seed_rows(seed: &fieldcraft_model::Seed, units: &DisplaySystem) -> Vec<(String, String)> {
    let mut out = vec![("Edge".to_string(), format!("{} · {}", seed.edge.body, seed.edge.entity.name().unwrap_or_default()))];
    let params = crate::menu::seed_params(seed, units).unwrap_or_default();
    if let (Some(fields), Some(obj)) = (fieldcraft_engine::find_command("mesh.seed.edge").and_then(|s| s.form), params.as_object()) {
        let mut obj = obj.clone();
        obj.remove("edges");
        object(fields, &obj, units, &mut out, &mut Vec::new(), 0);
    }
    out
}

/// "contact_pressure" → "Contact pressure".
fn sentence(key: &str) -> String {
    let s = key.replace('_', " ");
    let mut c = s.chars();
    c.next().map_or_else(String::new, |f| f.to_uppercase().chain(c).collect())
}

/// The rows of one object's fields; `used` collects the top-level keys described.
fn object(fields: &[Field], obj: &Map<String, Value>, units: &DisplaySystem, out: &mut Vec<(String, String)>, used: &mut Vec<String>, depth: usize) {
    if depth > MAX_DEPTH {
        return;
    }
    for f in fields {
        if f.key == "name" {
            continue;
        }
        if let Kind::Variants { tag, variants } = f.kind
            && f.key.is_empty()
        {
            let Some(t) = obj.get(tag).and_then(Value::as_str) else { continue };
            used.push(tag.to_string());
            match variants.iter().find(|v| v.tag == t) {
                Some(v) => {
                    out.push((f.label.to_string(), v.label.to_string()));
                    object(v.fields, obj, units, out, used, depth + 1);
                }
                None => out.push((f.label.to_string(), t.to_string())),
            }
            continue;
        }
        let Some(value) = obj.get(f.key).filter(|v| !v.is_null()) else { continue };
        used.push(f.key.to_string());
        match f.kind {
            Kind::Group(inner) => {
                if let Some(o) = value.as_object() {
                    object(inner, o, units, out, &mut Vec::new(), depth + 1);
                }
            }
            Kind::Variants { tag, variants } => {
                let Some(o) = value.as_object() else { continue };
                let t = o.get(tag).and_then(Value::as_str).unwrap_or_default();
                match variants.iter().find(|v| v.tag == t) {
                    Some(v) => {
                        out.push((f.label.to_string(), v.label.to_string()));
                        object(v.fields, o, units, out, &mut Vec::new(), depth + 1);
                    }
                    None => out.push((f.label.to_string(), crate::context::value_text(value))),
                }
            }
            _ => out.push((f.label.to_string(), scalar(f, value, units))),
        }
    }
}

/// One field's value as text, with its unit.
fn scalar(f: &Field, v: &Value, units: &DisplaySystem) -> String {
    let unit = f.unit_label(units);
    match f.kind {
        Kind::Number { .. } => v.as_f64().map_or_else(|| crate::context::value_text(v), |x| with_unit(x, unit.as_deref())),
        Kind::Integer { .. } => v.as_i64().map_or_else(|| crate::context::value_text(v), |x| quantity(x as f64)),
        Kind::Bool => (if v.as_bool() == Some(true) { "yes" } else { "no" }).into(),
        Kind::Choice(options) => {
            let s = v.as_str().unwrap_or_default();
            options.iter().find(|(o, _)| *o == s).map_or_else(|| s.to_string(), |(_, l)| (*l).to_string())
        }
        Kind::Vector { labels, .. } => vector(labels, v, unit.as_deref()),
        Kind::Tuple { labels, .. } => vector(labels, v, None),
        Kind::Ref(_) | Kind::Text | Kind::Path { .. } => v.as_str().map_or_else(|| crate::context::value_text(v), str::to_string),
        Kind::RefList(_) => list(v, |x| x.as_str().map(str::to_string)),
        // `Bracket.3: Tet, Bracket.7: Swept`.
        Kind::RefChoices(_, options) => v.as_object().map_or_else(
            || crate::context::value_text(v),
            |o| {
                o.iter()
                    .map(|(k, x)| {
                        let s = x.as_str().unwrap_or_default();
                        format!("{k}: {}", options.iter().find(|(o, _)| *o == s).map_or(s, |(_, l)| *l))
                    })
                    .collect::<Vec<_>>()
                    .join(", ")
            },
        ),
        Kind::NameOr { .. } => match v {
            Value::String(s) => s.clone(),
            _ => crate::context::value_text(v),
        },
        Kind::RefText(_) => v.as_object().map_or_else(
            || crate::context::value_text(v),
            |o| o.iter().map(|(k, x)| format!("{k}: {}", x.as_str().unwrap_or_default())).collect::<Vec<_>>().join(", "),
        ),
        Kind::Either(options) => match v {
            Value::String(s) => options.iter().find(|o| o.key == s).map_or_else(|| s.clone(), |o| o.label.to_string()),
            Value::Object(o) => o
                .iter()
                .next()
                .map(|(k, x)| match options.iter().find(|op| op.key == k) {
                    Some(op) => format!("{} {}", op.label, scalar(op, x, units)),
                    None => crate::context::value_text(v),
                })
                .unwrap_or_default(),
            _ => crate::context::value_text(v),
        },
        Kind::Rows { .. } => v.as_array().map_or_else(|| crate::context::value_text(v), |a| plural(a.len(), "row")),
        Kind::Sides(side) => match v
            .as_array()
            .and_then(|a| a.iter().map(|x| Some((x.get("element")?.as_u64()?, x.get(side.key())?.as_u64()?))).collect::<Option<Vec<_>>>())
        {
            Some(l) => crate::fields::format_sides(&l),
            None => crate::context::value_text(v),
        },
        Kind::Labels(_) => match v.as_array().and_then(|a| a.iter().map(Value::as_u64).collect::<Option<Vec<_>>>()) {
            Some(l) => crate::fields::format_labels(&l),
            None => crate::context::value_text(v),
        },
        Kind::Points { columns, .. } => {
            v.as_array().map_or_else(|| crate::context::value_text(v), |a| plural(a.len(), if columns == schema::XYZ { "point" } else { "row" }))
        }
        Kind::Multi(options) => {
            list(v, |x| options.iter().find(|(o, _)| fieldcraft_engine::schema::multi_value(o) == *x).map(|(_, l)| (*l).to_string()))
        }
        Kind::PerAxis { .. } => match v {
            Value::Array(a) => a
                .iter()
                .map(|x| x.as_f64().map_or_else(|| crate::context::value_text(x), |x| with_unit(x, unit.as_deref())))
                .collect::<Vec<_>>()
                .join(", "),
            _ => v.as_f64().map_or_else(|| crate::context::value_text(v), |x| with_unit(x, unit.as_deref())),
        },
        Kind::NumberOrAuto { .. } => match v {
            Value::String(s) => s.clone(),
            _ => v.as_f64().map_or_else(|| crate::context::value_text(v), |x| with_unit(x, unit.as_deref())),
        },
        Kind::Json | Kind::Unit | Kind::Group(_) | Kind::Variants { .. } => crate::context::value_text(v),
    }
}

/// A vector by its components: "x 0, y −1 000, z 0 N". Components left free (null) are listed
/// last as free; rotations (`ur…`) are in radians.
fn vector(labels: &[&str], v: &Value, unit: Option<&str>) -> String {
    let Some(a) = v.as_array() else { return crate::context::value_text(v) };
    let (mut set, mut free) = (Vec::new(), Vec::new());
    for (i, x) in a.iter().enumerate() {
        let label = labels.get(i).copied().unwrap_or("?");
        match x.as_f64() {
            Some(x) => {
                let u = if label.starts_with("ur") { Some("rad") } else { unit };
                set.push(format!("{} {}", label.to_uppercase(), with_unit(x, u)));
            }
            None => free.push(label.to_uppercase()),
        }
    }
    // A force or moment vector shares one unit: given once, at the end.
    if free.is_empty() && labels.iter().all(|l| !l.starts_with("ur")) {
        let parts: Vec<String> = a
            .iter()
            .enumerate()
            .map(|(i, x)| format!("{} {}", labels.get(i).copied().unwrap_or("?"), x.as_f64().map_or("—".into(), quantity)))
            .collect();
        return match unit {
            Some(u) => format!("{} {u}", parts.join(", ")),
            None => parts.join(", "),
        };
    }
    let mut out = set.join(", ");
    if !free.is_empty() {
        out.push_str(&format!("; {} free", free.join(", ")));
    }
    out
}

fn list(v: &Value, item: impl Fn(&Value) -> Option<String>) -> String {
    let Some(a) = v.as_array() else { return crate::context::value_text(v) };
    let names: Vec<String> = a.iter().filter_map(item).collect();
    if names.len() > LIST { format!("{}, … ({} in all)", names[..LIST].join(", "), names.len()) } else { names.join(", ") }
}

fn plural(n: usize, word: &str) -> String {
    if n == 1 { format!("1 {word}") } else { format!("{n} {word}s") }
}

/// A set: what it is defined by (design entities by kind and body, mesh labels) and what that
/// resolves to in the mesh.
fn set_rows(model: &Model, name: &str, obj: &Map<String, Value>) -> Vec<(String, String)> {
    let mut out = Vec::new();
    if let Some(geometry) = obj.get("geometry").and_then(Value::as_array).filter(|g| !g.is_empty()) {
        let mut kinds: std::collections::BTreeMap<&str, usize> = Default::default();
        let mut bodies = std::collections::BTreeSet::new();
        for g in geometry {
            let kind = ["face", "edge", "vertex"].into_iter().find(|k| g.get(*k).is_some()).unwrap_or("body");
            *kinds.entry(kind).or_default() += 1;
            if let Some(b) = g.get("body").and_then(Value::as_str) {
                bodies.insert(b);
            }
        }
        let what: Vec<String> = kinds
            .iter()
            .map(|(k, n)| {
                if *k == "vertex" && *n != 1 {
                    format!("{n} vertices")
                } else if *k == "body" && *n != 1 {
                    format!("{n} bodies")
                } else {
                    plural(*n, k)
                }
            })
            .collect();
        let on = if bodies.len() == 1 && !kinds.contains_key("body") {
            bodies.iter().next().map(|b| format!(" of {b}")).unwrap_or_default()
        } else if !kinds.contains_key("body") {
            format!(" on {} bodies", bodies.len())
        } else {
            String::new()
        };
        out.push(("Geometry".into(), format!("{}{on}", what.join(", "))));
    }
    for (key, label, word) in
        [("nodes", "Nodes", "node"), ("elements", "Elements", "element"), ("faces", "Element faces", "face"), ("edges", "Element edges", "edge")]
    {
        if let Some(a) = obj.get(key).and_then(Value::as_array).filter(|a| !a.is_empty()) {
            out.push((label.into(), plural(a.len(), word)));
        }
    }
    if let Some(r) = model.resolve_set_partial(name) {
        let mut parts = vec![plural(r.nodes.len(), "node")];
        if !r.elements.is_empty() {
            parts.push(plural(r.elements.len(), "element"));
        }
        if !r.faces.is_empty() {
            parts.push(plural(r.faces.len(), "element face"));
        }
        if !r.edges.is_empty() {
            parts.push(plural(r.edges.len(), "element edge"));
        }
        out.push(("In the mesh".into(), parts.join(", ")));
    }
    out
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;

    use super::*;

    fn text(rows: &[(String, String)]) -> String {
        rows.iter().map(|(k, v)| format!("{k}: {v}")).collect::<Vec<_>>().join("\n")
    }

    #[test]
    fn numbers_are_grouped_and_rounded() {
        assert_eq!(quantity(210000.0), "210\u{2009}000");
        assert_eq!(quantity(-1000.0), "−1\u{2009}000");
        assert_eq!(quantity(0.3), "0.3");
        assert_eq!(quantity(-2.5), "−2.5");
        assert_eq!(quantity(7.85e-9), "7.85e−9");
    }

    /// Definitions read in their forms' words and units, never as JSON (#335).
    #[test]
    fn definitions_read_with_labels_and_units() {
        let mut s = Session::default();
        s.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let m = &s.model;
        let mm = DisplaySystem::new(fieldcraft_model::units::Preset::MmTonneS);
        let steel = text(&rows(m, Def::Material, "Steel", &mm).unwrap());
        assert!(steel.contains("Elasticity: Isotropic") && steel.contains("210\u{2009}000 MPa") && steel.contains("7.85e−9 t/mm³"), "{steel}");
        // The same steel in SI: the stored numbers, in their units.
        let si = text(&rows(m, Def::Material, "Steel", &DisplaySystem::new(fieldcraft_model::units::Preset::MKgS)).unwrap());
        assert!(si.contains("2.1e11 Pa") && si.contains("7\u{2009}850 kg/m³"), "{si}");
        assert!(!steel.contains('{'), "{steel}");
        let load = Def::Load.names(m).into_iter().next().unwrap();
        let load = text(&rows(m, Def::Load, &load, &mm).unwrap());
        assert!(load.contains("Type: Force") && load.contains("Target: Set ") && load.contains(" N"), "{load}");
        assert!(!load.contains('{') && !load.contains('['), "{load}");
        let set = Def::Set.names(m).into_iter().next().unwrap();
        let set = text(&rows(m, Def::Set, &set, &mm).unwrap());
        assert!(set.contains("In the mesh: "), "{set}");

        let mut s = Session::default();
        s.run_script(&fieldcraft_engine::sample::rubber_seal()).unwrap();
        let press = text(&rows(&s.model, Def::BoundaryCondition, "Press", &mm).unwrap());
        assert!(press.contains("Type: Displacement") && press.contains("U3 −2.5 mm") && press.contains("UR1, UR2, UR3 free"), "{press}");
    }
}
