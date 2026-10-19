//! `results.xy`: one quantity of a solved step against its natural axis — a nodal field per
//! frame against the load factor (nonlinear static) or time (explicit, transient heat), a mode
//! shape's value or the frequencies against the mode number, and an explicit step's histories
//! (energies, momentum, constraint force, history nodes) against time. Read-only; with `path`
//! the series is also written as a two-column CSV.

use serde_json::{Value, json};

use super::opt;
use super::results::{Cell, Table, Want, step_for};
use crate::requirement::Requirement;
use crate::results::{Field, Kind, Pick, Quantity, field_values, with_unit};
use crate::schema::ProblemKind;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};

pub const XY: CommandSpec = CommandSpec {
    id: "results.xy",
    label: "XY Plot",
    params: r#"{"step": "Load", "y": {"field": "displacement_z", "node": 21}, "x": "auto", "every": 1, "path": "tip.csv"} — a quantity of a solved step against its natural axis: `x` "load_factor" (nonlinear static), "time" (explicit, transient heat), "mode" (frequency and buckling) or "auto" (the step's own, the default). `y` is {"field": <a contour field id>, "node": <label>} (a nodal field per frame or mode: on nonlinear static steps the displacements, whose frames start at load factor 0; on heat steps the temperature), {"energy": "kinetic|internal|hourglass|viscous|external|constraint|contact|friction|plastic|eroded|added_kinetic|balance"}, {"history": "constraint_force_x|_y|_z|momentum_x|_y|_z|increment"}, {"history_node": 21, "quantity": "displacement|displacement_x|…|velocity_z"} (explicit histories), {"reaction": "top", "component": "x|y|z|magnitude"} (nonlinear static: the reaction sum on a set per frame), "frequency" (frequency steps) or "load_factor" (buckling steps). On nonlinear static steps `x` may instead be {"displacement": "top", "component": "x|y|z|magnitude"} (the mean displacement of a set's nodes per frame) or {"displacement_node": 21, "component": …}: a reaction against it is a force–travel curve. An explicit step's `constraint_force` history is its reaction against time. At most 100000 points: `every` keeps one in so many (by default enough to fit). A stale result is refused. With `path` the series is written as CSV (a header line, then x,y rows). Replies {step, x: {label, values}, y: {label, values}, points, every, frames?} (`frames`: the field frame of each point, for results.frame)"#,
    form: Some(super::forms::RESULTS_XY),
    kind: CommandKind::Session,
    requires: &[Requirement::Result],
    run: xy,
};

/// Most points a series holds.
pub const MAX_POINTS: usize = 100_000;

/// What a series is plotted against.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Axis {
    LoadFactor,
    Time,
    Mode,
}

impl Axis {
    fn id(self) -> &'static str {
        match self {
            Axis::LoadFactor => "load_factor",
            Axis::Time => "time",
            Axis::Mode => "mode",
        }
    }

    fn parse(id: &str) -> Option<Axis> {
        [Axis::LoadFactor, Axis::Time, Axis::Mode].into_iter().find(|a| a.id() == id)
    }

    fn label(self, s: &Session) -> String {
        match self {
            Axis::LoadFactor => "Load factor".into(),
            Axis::Time => with_unit("Time", Quantity::Time.unit(s.model.units()).as_deref()),
            Axis::Mode => "Mode".into(),
        }
    }
}

/// The energies of an explicit step (`y: {"energy": …}`).
pub const ENERGIES: [&str; 12] =
    ["kinetic", "internal", "hourglass", "viscous", "external", "constraint", "contact", "friction", "plastic", "eroded", "added_kinetic", "balance"];
/// The other histories of an explicit step (`y: {"history": …}`).
pub const HISTORIES: [&str; 7] =
    ["constraint_force_x", "constraint_force_y", "constraint_force_z", "momentum_x", "momentum_y", "momentum_z", "increment"];
/// What a history node keeps (`y: {"history_node": n, "quantity": …}`).
pub const NODE_QUANTITIES: [&str; 8] =
    ["displacement", "displacement_x", "displacement_y", "displacement_z", "velocity", "velocity_x", "velocity_y", "velocity_z"];

/// A series: its labels, its x values and a function giving the point at an index.
struct Series<'a> {
    x_label: String,
    y_label: String,
    len: usize,
    point: Box<dyn Fn(usize) -> Result<(f64, f64)> + 'a>,
    /// Whether point `i` is field frame `i` (results.frame can show it).
    frames: bool,
}

fn y_error(why: impl Into<String>) -> EngineError {
    EngineError::invalid("y", why)
}

fn magnitude(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// A vector's component or magnitude by suffix: `""` the magnitude, `"_x"`… a component.
fn component(v: &[f64; 3], suffix: &str) -> Option<f64> {
    match suffix {
        "" => Some(magnitude(v)),
        "_x" => Some(v[0]),
        "_y" => Some(v[1]),
        "_z" => Some(v[2]),
        _ => None,
    }
}

/// A node label from `y.node` or `y.history_node`.
fn node_label(y: &Value, key: &str) -> Result<u32> {
    let v = y.get(key).ok_or_else(|| y_error(format!("needs `{key}`: a node label")))?;
    v.as_u64()
        .and_then(|n| u32::try_from(n).ok())
        .filter(|n| *n > 0)
        .ok_or_else(|| y_error(format!("`{key}` must be a node label from 1 to {}", u32::MAX)))
}

/// The fields a step plots per frame or mode at a node.
pub fn plottable(r: &StepResult) -> Vec<Field> {
    let shapes = [Field::Displacement, Field::DisplacementX, Field::DisplacementY, Field::DisplacementZ];
    match r.kind() {
        // Nonlinear frames hold displacements only; the other fields are the final state's.
        Kind::Static(_) | Kind::Frequency(_) => shapes.to_vec(),
        Kind::Explicit(_) => Field::available(r).collect(),
        Kind::Heat(_) => vec![Field::Temperature],
        Kind::Fluid(_) => Vec::new(),
    }
}

fn ids(list: &[&str]) -> String {
    list.join(", ")
}

/// The series `y` asks for on step `step`'s result `r`, against `axis`.
fn series<'a>(s: &'a Session, step: &str, r: &'a StepResult, axis: Axis, y: &Value) -> Result<Series<'a>> {
    let x_label = axis.label(s);
    let units = s.model.units();
    let not_explicit = || y_error(format!("step `{step}` is not an explicit dynamic step; it has no histories"));
    if let Some(name) = y.as_str() {
        let Kind::Frequency(f) = r.kind() else {
            return Err(y_error(format!("`{name}` needs a frequency or buckling step; `{step}` is not one")));
        };
        let buckling = f.buckling.is_some();
        let (y_label, field): (String, fn(&fieldcraft_solver::Mode) -> f64) = match (name, buckling) {
            ("frequency", false) => (with_unit("Frequency", Some("cycles per time unit")), |m| m.frequency),
            ("load_factor", true) => ("Buckling load factor".into(), |m| m.eigenvalue),
            ("frequency", true) => return Err(y_error(format!("step `{step}` is a buckling step: its modes are load factors (y: \"load_factor\")"))),
            ("load_factor", false) => return Err(y_error(format!("step `{step}` is a frequency step: y: \"frequency\""))),
            _ => {
                return Err(y_error(format!(
                    "`{name}` is not a series; use \"frequency\", \"load_factor\" or an object ({{\"field\", \"node\"}}, …)"
                )));
            }
        };
        return Ok(Series {
            x_label,
            y_label,
            len: f.modes.len(),
            point: Box::new(move |i| Ok((i as f64 + 1.0, f.modes.get(i).map_or(f64::NAN, field)))),
            frames: false,
        });
    }
    if !y.is_object() {
        return Err(y_error(
            "must be an object ({\"field\", \"node\"}, {\"energy\"}, {\"history\"} or {\"history_node\", \"quantity\"}) or \"frequency\"",
        ));
    }
    if let Some(id) = y.get("field") {
        let available = plottable(r);
        let list = || available.iter().map(|f| f.id()).collect::<Vec<_>>().join(", ");
        let field = id.as_str().and_then(Field::parse).filter(|f| available.contains(f)).ok_or_else(|| {
            if matches!(r.kind(), Kind::Fluid(_)) {
                y_error(format!("step `{step}` is a fluid step: its fields are per cell (solve.probe)"))
            } else {
                y_error(format!("`field` must be one step `{step}` has per {}: {}", if axis == Axis::Mode { "mode" } else { "frame" }, list()))
            }
        })?;
        let label = node_label(y, "node")?;
        let at = r.nodes().iter().position(|n| *n == label).ok_or_else(|| y_error(format!("node {label} has no result in step `{step}`")))?;
        let model = &s.model;
        let value = move |pick: Pick| -> Result<f64> { Ok(field_values(model, r, field, pick)?.get(at).copied().unwrap_or(f64::NAN)) };
        let y_label = format!("{} at node {label}", field.title(r, units));
        let (len, frames, x): (usize, bool, Box<dyn Fn(usize) -> f64 + 'a>) = match r.kind() {
            Kind::Static(st) => {
                let h = st.nonlinear.as_ref().ok_or_else(|| {
                    y_error(format!("step `{step}` is a linear static step: one state, nothing to plot against (a nonlinear static step keeps a frame per increment)"))
                })?;
                (h.frames.len(), true, Box::new(move |i| h.frames.get(i).map_or(f64::NAN, |f| f.load_factor)))
            }
            Kind::Explicit(e) => (e.frames.len(), true, Box::new(move |i| e.frames.get(i).map_or(f64::NAN, |f| f.time))),
            Kind::Heat(h) => (h.frames.len(), true, Box::new(move |i| h.times.get(i).copied().unwrap_or(f64::NAN))),
            Kind::Frequency(f) => (f.modes.len(), false, Box::new(|i| i as f64 + 1.0)),
            Kind::Fluid(_) => return Err(y_error(format!("step `{step}` is a fluid step: its fields are per cell (solve.probe)"))),
        };
        let point = move |i: usize| -> Result<(f64, f64)> {
            let pick = if frames { Pick { mode: 1, frame: Some(i) } } else { Pick { mode: i + 1, frame: None } };
            Ok((x(i), value(pick)?))
        };
        return Ok(Series { x_label, y_label, len, point: Box::new(point), frames });
    }
    if let Some(set) = y.get("reaction") {
        let frames = nonlinear_frames(step, r).map_err(y_error)?;
        let set = set.as_str().ok_or_else(|| y_error("`reaction` must be a set name"))?;
        let (suffix, title) = vector_component(y).map_err(y_error)?;
        let at = set_indices(s, r, set)?;
        let point = move |i: usize| -> Result<(f64, f64)> {
            let Some(f) = frames.get(i) else { return Ok((f64::NAN, f64::NAN)) };
            let sum = at.iter().map(|k| f.reaction_at(*k)).fold([0.0; 3], |a, v| std::array::from_fn(|c| a[c] + v[c]));
            Ok((f.load_factor, component(&sum, suffix).unwrap_or(f64::NAN)))
        };
        return Ok(Series {
            x_label,
            y_label: format!("{} on {set}", with_unit(&format!("Reaction {title}"), Quantity::Force.unit(units).as_deref())),
            len: frames.len(),
            point: Box::new(point),
            frames: true,
        });
    }
    let Kind::Explicit(e) = r.kind() else { return Err(not_explicit()) };
    let h = &e.history;
    let time = move |i: usize| h.time.get(i).copied().unwrap_or(f64::NAN);
    if let Some(name) = y.get("energy") {
        let name = name.as_str().filter(|n| ENERGIES.contains(n)).ok_or_else(|| y_error(format!("`energy` must be one of {}", ids(&ENERGIES))))?;
        let pick: fn(&fieldcraft_explicit::Energies) -> f64 = match name {
            "kinetic" => |e| e.kinetic,
            "internal" => |e| e.internal,
            "hourglass" => |e| e.hourglass,
            "viscous" => |e| e.viscous,
            "external" => |e| e.external,
            "constraint" => |e| e.constraint,
            "contact" => |e| e.contact,
            "friction" => |e| e.friction,
            "plastic" => |e| e.plastic,
            "eroded" => |e| e.eroded,
            "added_kinetic" => |e| e.added_kinetic,
            _ => |e| e.balance,
        };
        let label = format!("{} energy", name.replace('_', " "));
        let label = label.get(..1).map(str::to_uppercase).unwrap_or_default() + label.get(1..).unwrap_or_default();
        return Ok(Series {
            x_label,
            y_label: with_unit(&label, Quantity::Energy.unit(units).as_deref()),
            len: h.time.len(),
            point: Box::new(move |i| Ok((time(i), h.energies.get(i).map_or(f64::NAN, pick)))),
            frames: false,
        });
    }
    if let Some(name) = y.get("history") {
        let name = name.as_str().filter(|n| HISTORIES.contains(n)).ok_or_else(|| y_error(format!("`history` must be one of {}", ids(&HISTORIES))))?;
        let axis_of = |suffix: &str| match suffix {
            "x" => 0,
            "y" => 1,
            _ => 2,
        };
        let (y_label, value): (String, Box<dyn Fn(usize) -> f64 + 'a>) = if let Some(c) = name.strip_prefix("constraint_force_") {
            let k = axis_of(c);
            (
                with_unit(&format!("Constraint force {}", c.to_uppercase()), Quantity::Force.unit(units).as_deref()),
                Box::new(move |i| h.constraint_force.get(i).map_or(f64::NAN, |v| v[k])),
            )
        } else if let Some(c) = name.strip_prefix("momentum_") {
            let k = axis_of(c);
            (format!("Momentum {}", c.to_uppercase()), Box::new(move |i| h.momentum.get(i).map_or(f64::NAN, |v| v[k])))
        } else {
            (with_unit("Time increment", Quantity::Time.unit(units).as_deref()), Box::new(move |i| h.increment.get(i).copied().unwrap_or(f64::NAN)))
        };
        return Ok(Series { x_label, y_label, len: h.time.len(), point: Box::new(move |i| Ok((time(i), value(i)))), frames: false });
    }
    if y.get("history_node").is_some() {
        let label = node_label(y, "history_node")?;
        let k = h.nodes.iter().position(|n| *n == label).ok_or_else(|| {
            let listed: Vec<String> = h.nodes.iter().take(super::LISTED).map(u32::to_string).collect();
            y_error(format!(
                "node {label} has no history in step `{step}`; {}",
                if listed.is_empty() { "the step keeps no node histories".to_string() } else { format!("history nodes: {}", listed.join(", ")) }
            ))
        })?;
        let quantity = y
            .get("quantity")
            .and_then(Value::as_str)
            .filter(|q| NODE_QUANTITIES.contains(q))
            .ok_or_else(|| y_error(format!("`quantity` must be one of {}", ids(&NODE_QUANTITIES))))?;
        let (values, suffix, title, q) = match quantity.strip_prefix("displacement") {
            Some(suffix) => (&h.displacement, suffix, "Displacement", Quantity::Length),
            None => (&h.velocity, quantity.trim_start_matches("velocity"), "Velocity", Quantity::Velocity),
        };
        let name =
            if suffix.is_empty() { format!("{title} magnitude") } else { format!("{title} {}", suffix.trim_start_matches('_').to_uppercase()) };
        let suffix = suffix.to_string();
        return Ok(Series {
            x_label,
            y_label: format!("{} at node {label}", with_unit(&name, q.unit(units).as_deref())),
            len: h.time.len(),
            point: Box::new(move |i| Ok((time(i), values.get(i).and_then(|row| row.get(k)).and_then(|v| component(v, &suffix)).unwrap_or(f64::NAN)))),
            frames: false,
        });
    }
    Err(y_error("must have one of `field` (with `node`), `reaction`, `energy`, `history` or `history_node` (with `quantity`)"))
}

/// The frames of a nonlinear static step.
fn nonlinear_frames<'a>(step: &str, r: &'a StepResult) -> std::result::Result<&'a [fieldcraft_solver::nonlinear::Frame], String> {
    match r.kind() {
        Kind::Static(st) => st.nonlinear.as_ref().map(|h| h.frames.as_slice()).ok_or_else(|| {
            format!("step `{step}` is a linear static step: one state (a nonlinear static step keeps its reactions and displacements per increment)")
        }),
        _ => Err(format!(
            "step `{step}` is not a nonlinear static step: reactions and set displacements per frame are kept for nonlinear static steps"
        )),
    }
}

/// A vector's `component` (x, y, z or magnitude, the default): its [`component`] suffix and
/// its title.
fn vector_component(v: &Value) -> std::result::Result<(&'static str, &'static str), String> {
    match v.get("component").map(|c| c.as_str()) {
        None | Some(Some("magnitude")) => Ok(("", "magnitude")),
        Some(Some("x")) => Ok(("_x", "X")),
        Some(Some("y")) => Ok(("_y", "Y")),
        Some(Some("z")) => Ok(("_z", "Z")),
        _ => Err("`component` must be \"x\", \"y\", \"z\" or \"magnitude\"".into()),
    }
}

/// The result indices of a set's nodes (an error for an unknown set or one without results).
fn set_indices(s: &Session, r: &StepResult, set: &str) -> Result<Vec<usize>> {
    let nodes = r.nodes();
    let at: Vec<usize> = s.model.resolve_set(set)?.nodes.iter().filter_map(|l| nodes.binary_search(l).ok()).collect();
    if at.is_empty() {
        return Err(EngineError::Params(format!("set `{set}` has no nodes with results")));
    }
    Ok(at)
}

/// `x` as a travel on a nonlinear static step: {"displacement": set} (the mean over its nodes)
/// or {"displacement_node": label}, and a component. Its label and value per frame.
fn travel<'a>(s: &'a Session, step: &str, r: &'a StepResult, x: &Value) -> Result<(String, Box<dyn Fn(usize) -> f64 + 'a>)> {
    let x_error = |why: String| EngineError::invalid("x", why);
    let frames = nonlinear_frames(step, r).map_err(x_error)?;
    let (suffix, title) = vector_component(x).map_err(x_error)?;
    let (at, of) = match (x.get("displacement"), x.get("displacement_node")) {
        (Some(set), None) => {
            let set = set.as_str().ok_or_else(|| x_error("`displacement` must be a set name".into()))?;
            (set_indices(s, r, set)?, format!("of {set} (mean)"))
        }
        (None, Some(node)) => {
            let label =
                node.as_u64().and_then(|n| u32::try_from(n).ok()).ok_or_else(|| x_error("`displacement_node` must be a node label".into()))?;
            let i = r.nodes().binary_search(&label).map_err(|_| x_error(format!("node {label} has no result in step `{step}`")))?;
            (vec![i], format!("at node {label}"))
        }
        _ => return Err(x_error("an object must have one of `displacement` (a set) or `displacement_node` (a node label)".into())),
    };
    let label = format!("{} {of}", with_unit(&format!("Displacement {title}"), Quantity::Length.unit(s.model.units()).as_deref()));
    let n = at.len() as f64;
    let value = move |i: usize| -> f64 {
        let Some(f) = frames.get(i) else { return f64::NAN };
        let sum = at.iter().map(|k| f.displacement.get(*k).copied().unwrap_or_default()).fold([0.0; 3], |a, v| std::array::from_fn(|c| a[c] + v[c]));
        component(&sum.map(|c| c / n), suffix).unwrap_or(f64::NAN)
    };
    Ok((label, Box::new(value)))
}

/// The step's own axis.
fn natural(step: &str, r: &StepResult) -> Result<Axis> {
    match r.kind() {
        Kind::Static(st) if st.nonlinear.is_some() => Ok(Axis::LoadFactor),
        Kind::Static(_) => Err(EngineError::Params(format!(
            "step `{step}` is a linear static step: one state, nothing to plot against (a nonlinear static step keeps a frame per increment)"
        ))),
        Kind::Explicit(_) | Kind::Heat(_) => Ok(Axis::Time),
        Kind::Frequency(_) => Ok(Axis::Mode),
        Kind::Fluid(_) => Err(EngineError::Params(format!("step `{step}` is a fluid step: it has no series to plot yet"))),
    }
}

fn xy(s: &mut Session, p: &Value) -> Result<Value> {
    let step = step_for(s, p, Want::Any)?;
    let (r, current) = s.latest_result(&step).ok_or_else(|| EngineError::Params(format!("step `{step}` has no result")))?;
    let x_travel = p.get("x").filter(|x| x.is_object());
    let axis = match x_travel {
        Some(_) => None,
        None => match opt::<String>(p, "x")?.as_deref() {
            None | Some("auto") => None,
            Some(id) => Some(Axis::parse(id).ok_or_else(|| {
                EngineError::invalid(
                    "x",
                    "must be \"auto\", \"load_factor\", \"time\", \"mode\" or a travel ({\"displacement\": set} or {\"displacement_node\": label})",
                )
            })?),
        },
    };
    let y = p.get("y").filter(|v| !v.is_null()).ok_or_else(|| EngineError::param("y", ProblemKind::Required))?;
    if !current {
        return Err(EngineError::Params(format!("the result of step `{step}` is stale (the model changed since it was solved); solve it again")));
    }
    let own = natural(&step, r)?;
    if let Some(a) = axis.filter(|a| *a != own) {
        return Err(EngineError::invalid("x", format!("step `{step}` is plotted against {}, not {}", own.id(), a.id())));
    }
    let every = match opt::<u64>(p, "every")? {
        Some(0) => return Err(EngineError::param("every", ProblemKind::Below(1.0))),
        Some(k) => Some(usize::try_from(k).unwrap_or(usize::MAX)),
        None => None,
    };
    let path = match p.get("path") {
        Some(_) => Some(super::file_param(s, p, "path", fieldcraft_io::roots::Access::Write)?),
        None => None,
    };
    let mut series = series(s, &step, r, own, y)?;
    let mut x_axis = own.id();
    if let Some(x) = x_travel {
        let (label, value) = travel(s, &step, r, x)?;
        if !series.frames {
            return Err(EngineError::invalid("x", "a travel needs a series per frame of the nonlinear static step"));
        }
        let point = series.point;
        series = Series { x_label: label, point: Box::new(move |i| Ok((value(i), point(i)?.1))), ..series };
        x_axis = "displacement";
    }
    let every = every.unwrap_or_else(|| series.len.div_ceil(MAX_POINTS).max(1));
    if series.len.div_ceil(every) > MAX_POINTS {
        return Err(EngineError::invalid("every", format!("{} points in steps of {every} are more than {MAX_POINTS}; raise `every`", series.len)));
    }
    let picked: Vec<usize> = (0..series.len).step_by(every).collect();
    let (mut xs, mut ys) = (Vec::with_capacity(picked.len()), Vec::with_capacity(picked.len()));
    for &i in &picked {
        let (x, y) = (series.point)(i)?;
        xs.push(x);
        ys.push(y);
    }
    let mut out = json!({
        "step": step,
        "x": {"axis": x_axis, "label": series.x_label, "values": xs},
        "y": {"label": series.y_label, "values": ys},
        "points": picked.len(),
        "every": every,
    });
    if series.frames {
        out["frames"] = json!(picked);
    }
    if let Some(path) = path {
        let cell = |v: f64| if v.is_finite() { Cell::Number(v) } else { Cell::Empty };
        let table = Table { header: vec![series.x_label, series.y_label], rows: xs.iter().zip(&ys).map(|(x, y)| vec![cell(*x), cell(*y)]).collect() };
        fieldcraft_io::fs::write_output(std::path::Path::new(&path), table.csv().as_bytes())?;
        out["path"] = json!(path);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cmd::solve::tests::cantilever;

    fn values(r: &Value, axis: &str) -> Vec<f64> {
        r[axis]["values"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap_or(f64::NAN)).collect()
    }

    fn impact() -> Session {
        let mut s = Session::default();
        s.run_script(&crate::sample::impact()).unwrap();
        s.execute("solve.explicit", &json!({})).unwrap();
        s
    }

    /// Energies and constraint forces are solve.history's numbers, sample for sample.
    #[test]
    fn explicit_series_are_the_histories() {
        let mut s = impact();
        let history = s.execute("solve.history", &json!({"step": "Impact"})).unwrap();
        let r = s.execute("results.xy", &json!({"step": "Impact", "y": {"energy": "kinetic"}})).unwrap();
        let time: Vec<f64> = history["time"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        let kinetic: Vec<f64> = history["energies"].as_array().unwrap().iter().map(|e| e["kinetic"].as_f64().unwrap()).collect();
        assert_eq!(values(&r, "x"), time);
        assert_eq!(values(&r, "y"), kinetic);
        assert_eq!(r["x"]["axis"], "time");
        assert!(r.get("frames").is_none(), "history samples are not field frames");
        let r = s.execute("results.xy", &json!({"step": "Impact", "x": "time", "y": {"history": "constraint_force_x"}})).unwrap();
        let fx: Vec<f64> = history["constraint_force"].as_array().unwrap().iter().map(|f| f[0].as_f64().unwrap()).collect();
        assert_eq!(values(&r, "y"), fx);
        // A nodal field per field frame, with the frames for results.frame.
        let r = s.execute("results.xy", &json!({"step": "Impact", "y": {"field": "velocity", "node": 1}})).unwrap();
        assert_eq!(r["points"], 11);
        assert_eq!(r["frames"].as_array().unwrap().len(), 11);
        assert!((values(&r, "y")[0] - 10_000.0).abs() < 1e-6, "every node starts at 10 m/s: {r}");
        // `every` strides.
        let r = s.execute("results.xy", &json!({"step": "Impact", "y": {"energy": "balance"}, "every": 2})).unwrap();
        assert_eq!(r["points"].as_u64().unwrap() as usize, time.len().div_ceil(2));
    }

    /// The frequency series is the solve's frequencies against modes 1 to n.
    #[test]
    fn frequencies_against_the_mode_number() {
        let mut s = Session::default();
        cantilever(&mut s);
        let summary = s.execute("solve.frequency", &json!({})).unwrap();
        let r = s.execute("results.xy", &json!({"step": "Modes", "y": "frequency"})).unwrap();
        let expected: Vec<f64> = summary["modes"].as_array().unwrap().iter().map(|m| m["frequency"].as_f64().unwrap()).collect();
        assert_eq!(values(&r, "y"), expected);
        assert_eq!(values(&r, "x"), (1..=expected.len()).map(|m| m as f64).collect::<Vec<_>>());
        assert!(s.execute("results.xy", &json!({"step": "Modes", "y": "load_factor"})).is_err());
        assert!(s.execute("results.xy", &json!({"step": "Modes", "x": "time", "y": "frequency"})).is_err());
    }

    /// A small tip load on the nonlinear cantilever stays linear: the tip deflection grows in
    /// proportion to the load factor, from 0, to the linear static answer (P L³ / 3 E I plus
    /// shear), and its last point is solve.probe's at the last frame. The CSV holds the series.
    #[test]
    fn tip_deflection_against_the_load_factor() {
        let mut s = Session::default();
        cantilever(&mut s);
        s.run_script(&json!([
            {"command": "model.section.edit", "params": {"name": "bar", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {}}},
            {"command": "model.load.edit", "params": {"name": "push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -0.001]}},
            {"command": "model.step.create", "params": {"name": "Large", "procedure": {"type": "nonlinear_static", "initial_increment": 0.25},
                "boundary_conditions": ["clamp"], "loads": ["push"]}},
        ]))
        .unwrap();
        s.execute("solve.static", &json!({"step": "Bend"})).unwrap();
        let linear = -s.execute("solve.probe", &json!({"step": "Bend", "nodes": [11]})).unwrap()["nodes"][0]["displacement"][2].as_f64().unwrap();
        s.execute("solve.nonlinear", &json!({"step": "Large"})).unwrap();
        let path = std::env::temp_dir().join(format!("fieldcraft-xy-{}.csv", std::process::id())).to_string_lossy().into_owned();
        let r = s.execute("results.xy", &json!({"step": "Large", "y": {"field": "displacement_z", "node": 11}, "path": path})).unwrap();
        let (x, y) = (values(&r, "x"), values(&r, "y"));
        assert_eq!(r["x"]["axis"], "load_factor");
        assert_eq!((x[0], y[0]), (0.0, 0.0), "the undeformed state first");
        assert!((x.last().unwrap() - 1.0).abs() < 1e-12, "{r}");
        for (l, d) in x.iter().zip(&y).skip(1) {
            assert!((-d / l - linear).abs() / linear < 1e-3, "λ {l}: {d} against {linear} λ");
        }
        let frames = r["frames"].as_array().unwrap().len();
        let probe = s.execute("solve.probe", &json!({"step": "Large", "nodes": [11], "frame": frames - 1})).unwrap();
        assert_eq!(*y.last().unwrap(), probe["nodes"][0]["displacement"][2].as_f64().unwrap());
        let csv = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], "Load factor,Displacement Z at node 11");
        assert_eq!(lines.len(), x.len() + 1);
        assert_eq!(lines.last().unwrap().split(',').nth(1).unwrap().parse::<f64>().unwrap(), *y.last().unwrap());
        let _ = std::fs::remove_file(&path);
        // Stresses are the final state's only; a linear step has one state.
        assert!(s.execute("results.xy", &json!({"step": "Large", "y": {"field": "von_mises", "node": 11}})).is_err());
        assert!(s.execute("results.xy", &json!({"step": "Bend", "y": {"field": "displacement_z", "node": 11}})).is_err());
        assert!(s.execute("results.xy", &json!({"step": "Bend", "y": {"reaction": "root"}})).is_err());
        // The root's reaction balances the dead end load at every frame: proportional to λ.
        let r = s.execute("results.xy", &json!({"step": "Large", "y": {"reaction": "root", "component": "z"}})).unwrap();
        let (x, y) = (values(&r, "x"), values(&r, "y"));
        let end = *y.last().unwrap();
        assert!(end.abs() > 0.0 && r["frames"].as_array().unwrap().len() == x.len(), "{r}");
        for (l, f) in x.iter().zip(&y) {
            assert!((f - l * end).abs() < 1e-6 * end.abs(), "λ {l}: {f} against {}", l * end);
        }
        // Against the tip's travel: the same reactions.
        let t = s
            .execute(
                "results.xy",
                &json!({"step": "Large", "y": {"reaction": "root", "component": "z"}, "x": {"displacement_node": 11, "component": "z"}}),
            )
            .unwrap();
        assert_eq!(values(&t, "y"), y);
        assert_eq!(t["x"]["label"], "Displacement Z at node 11");
        // A stale result is refused.
        s.execute("mesh.node.add", &json!({"label": 99_999, "xyz": [0, 0, 0]})).unwrap();
        let e = s.execute("results.xy", &json!({"step": "Large", "y": {"field": "displacement_z", "node": 11}})).unwrap_err().to_string();
        assert!(e.contains("stale"), "{e}");
    }

    #[test]
    fn xy_hostile_params_never_panic() {
        let mut modes = Session::default();
        cantilever(&mut modes);
        modes.execute("solve.frequency", &json!({"step": "Modes"})).unwrap();
        let hostile = [
            json!(null),
            json!([]),
            json!({}),
            json!({"y": null}),
            json!({"y": 7}),
            json!({"y": []}),
            json!({"y": "nope"}),
            json!({"step": "Nope", "y": "frequency"}),
            json!({"step": "Impact", "y": {"field": "nope", "node": 1}}),
            json!({"step": "Impact", "y": {"field": "velocity", "node": 0}}),
            json!({"step": "Impact", "y": {"field": "velocity", "node": 4_294_967_296_u64}}),
            json!({"step": "Impact", "y": {"field": "velocity", "node": -1}}),
            json!({"step": "Impact", "y": {"field": "velocity", "node": f64::NAN}}),
            json!({"step": "Impact", "y": {"field": "velocity"}}),
            json!({"step": "Impact", "y": {"energy": 3}}),
            json!({"step": "Impact", "y": {"history": "nope"}}),
            json!({"step": "Impact", "y": {"history_node": 1e300, "quantity": "velocity"}}),
            json!({"step": "Impact", "y": {"history_node": 1, "quantity": "nope"}}),
            json!({"step": "Impact", "x": "load_factor", "y": {"energy": "kinetic"}}),
            json!({"step": "Impact", "x": 7, "y": {"energy": "kinetic"}}),
            json!({"step": "Impact", "y": {"energy": "kinetic"}, "every": 0}),
            json!({"step": "Impact", "y": {"energy": "kinetic"}, "every": 1e300}),
            json!({"step": "Impact", "y": {"energy": "kinetic"}, "path": ""}),
            json!({"step": "Impact", "y": {"energy": "kinetic"}, "path": "/no/such/dir/x.csv"}),
            json!({"step": "Modes", "y": {"energy": "kinetic"}}),
            json!({"step": "Modes", "y": {"field": "von_mises", "node": 1}}),
        ];
        for s in [impact(), modes] {
            for p in &hostile {
                let mut t = s.clone();
                assert!(t.execute("results.xy", p).is_err(), "{p}");
                assert_eq!(t.model, s.model);
            }
        }
        let mut empty = Session::default();
        let e = empty.execute("results.xy", &json!({"y": "frequency"})).unwrap_err().to_string();
        assert!(e.contains("solve"), "{e}");
    }
}
