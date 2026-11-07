//! Initial fields mapped from another analysis (#160): `model.initial_field.map` interpolates a
//! solved step's nodal temperature or stress (of another `.fieldcraft` model, on another mesh,
//! or of this one) onto this model's nodes, as the temperature or residual stress a static step
//! starts from ([`fieldcraft_model::initial`]); `model.initial_field.delete` removes one.
//!
//! Each target node is located in a solid element of the source mesh and takes the element's
//! shape-function interpolation of the source's nodal values there
//! ([`fieldcraft_fem::transfer`]): exact for uniform and linear fields on any two meshes. A node
//! just outside the source (a curved face meshed coarser there) takes the value at the closest
//! point of the closest element within `tolerance`; a node farther out gets no value and is
//! counted as unmapped.

use std::collections::BTreeMap;

use fieldcraft_fem::transfer::Locator;
use fieldcraft_io::roots::Access;
use fieldcraft_model::{InitialField, Model};
use fieldcraft_solver::StepResult as Implicit;
use serde_json::{Value, json};

use super::{file_param, name, opt, req};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};

pub const MAP: CommandSpec = CommandSpec {
    id: "model.initial_field.map",
    label: "Map Initial Field",
    params: r#"{"name": "Residual", "quantity": "temperature" | "stress", "step": "Cool", "file": "moulding.fieldcraft", "frame": 12, "set": "part", "tolerance": 0.5} — interpolates step `step`'s nodal temperature (a heat transfer step's frame, the last when `frame` is left out, or a coupled temperature–displacement step's) or nodal stress (a static step's) from the model in `file` (its saved results, solved from its current definition; this model's when `file` is left out) onto this model's nodes (those of `set`, or every node of a solid element), by the shape functions of the source element holding each node: the temperature static steps measure thermal strain from, or the residual stress they start from. Nodes outside the source mesh by at most `tolerance` (a length; a quarter of the source's mean element size by default) take the value at the closest point; nodes farther out get none. Replaces a field of that name. Replies with the counts mapped `inside` and `nearby`, `unmapped`, and the largest distance"#,
    form: Some(super::forms::INITIAL_FIELD_MAP),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: map,
};

pub const DELETE: CommandSpec = CommandSpec {
    id: "model.initial_field.delete",
    label: "Delete Initial Field",
    params: r#"{"name": "Residual"} — removes the initial field"#,
    form: Some(super::forms::INITIAL_FIELD_DELETE),
    kind: CommandKind::Model,
    requires: &[],
    run: |s: &mut Session, p: &Value| {
        s.model.remove_initial_field(name(p)?)?;
        Ok(Value::Null)
    },
};

/// The default tolerance, as a fraction of the source's mean element size.
const DEFAULT_TOLERANCE: f64 = 0.25;

/// Target nodes located between two checks for cancellation.
const CHECK_EVERY: usize = 4096;

fn map(s: &mut Session, p: &Value) -> Result<Value> {
    let field_name = name(p)?.to_string();
    let stress = match req::<String>(p, "quantity")?.as_str() {
        "temperature" => false,
        "stress" => true,
        _ => return Err(EngineError::invalid("quantity", "must be `temperature` or `stress`")),
    };
    let step: String = req(p, "step")?;
    let frame: Option<usize> = opt(p, "frame")?;
    let set: Option<String> = opt(p, "set")?;
    let tolerance: Option<f64> = opt(p, "tolerance")?;
    if let Some(t) = tolerance
        && !(t >= 0.0 && t.is_finite())
    {
        return Err(EngineError::invalid("tolerance", "must be a length of zero or more"));
    }
    // The source: another model's file, or this model.
    let (file_model, source, label) = match opt::<String>(p, "file")? {
        Some(_) => {
            let path = super::file::absolute(&file_param(s, p, "file", Access::Read)?);
            let bytes = super::file::read_file(&path, fieldcraft_io::MAX_FILE_BYTES)?;
            let file = fieldcraft_io::read(&bytes).map_err(|e| EngineError::Params(format!("{path}: {e}")))?;
            let stored = file
                .results
                .iter()
                .find(|r| r.result.step() == step)
                .ok_or_else(|| EngineError::invalid("step", format!("`{path}` keeps no results of a step `{}`", short(&step))))?;
            if fieldcraft_io::results::definition_hash(&file.model, &step, &file.geometry) != Some(stored.definition) {
                return Err(EngineError::invalid(
                    "step",
                    format!("the results of step `{step}` in `{path}` were solved from another definition than its model's; solve it again and save"),
                ));
            }
            let fieldcraft_io::results::Solved::Implicit(result) = &stored.result else {
                return Err(not_implicit(&step));
            };
            (Some(file.model), values(result, stress, frame)?, format!("{path}, step {step}"))
        }
        None => {
            let result = s.result(&step).ok_or_else(|| EngineError::invalid("step", format!("step `{}` has no current result", short(&step))))?;
            let StepResult::Implicit(result) = result.as_ref() else {
                return Err(not_implicit(&step));
            };
            (None, values(result, stress, frame)?, format!("step {step}"))
        }
    };
    let source_model = file_model.as_ref().unwrap_or(&s.model);
    let label = match frame {
        Some(f) => format!("{label}, frame {f}"),
        None => label,
    };
    // The source's solid elements whose nodes all have a value, and those values (a
    // temperature in the first component).
    let value_at = |n: &u32| match &source {
        InitialField::Temperature { values, .. } => values.get(n).map(|t| [*t, 0.0, 0.0, 0.0, 0.0, 0.0]),
        InitialField::Stress { values, .. } => values.get(n).copied(),
    };
    let mut elements = Vec::new();
    let mut element_values = Vec::new();
    for e in source_model.elements().values().filter(|e| e.kind.dimension() == 3) {
        let xyz: Option<Vec<[f64; 3]>> = e.nodes.iter().map(|n| source_model.node(*n).map(|x| x.xyz)).collect();
        let at: Option<Vec<[f64; 6]>> = e.nodes.iter().map(value_at).collect();
        if let (Some(xyz), Some(at)) = (xyz, at) {
            elements.push((e.kind, xyz));
            element_values.push(at);
        }
    }
    if elements.is_empty() {
        return Err(EngineError::invalid("step", "the source has no solid elements with values to map from"));
    }
    let locator = Locator::new(elements).map_err(|e| EngineError::Params(format!("the source mesh: {e}")))?;
    let tolerance = tolerance.unwrap_or(DEFAULT_TOLERANCE * locator.mean_size());
    // The target nodes.
    let targets = target_nodes(&s.model, set.as_deref())?;
    let (mut inside, mut nearby, mut unmapped, mut farthest) = (0usize, 0usize, 0usize, 0.0f64);
    let mut mapped = BTreeMap::new();
    for (i, label) in targets.iter().enumerate() {
        if i % CHECK_EVERY == 0 {
            fieldcraft_model::progress::check()?;
        }
        let Some(x) = s.model.node(*label).map(|n| n.xyz) else { continue };
        let Some(hit) = locator.locate(x, tolerance) else {
            unmapped += 1;
            continue;
        };
        if hit.distance > 0.0 {
            nearby += 1;
            farthest = farthest.max(hit.distance);
        } else {
            inside += 1;
        }
        let Some(at) = element_values.get(hit.element) else { continue };
        let mut v = [0.0; 6];
        for (a, w) in at.iter().zip(&hit.weights) {
            for k in 0..6 {
                v[k] += w * a[k];
            }
        }
        mapped.insert(*label, v);
    }
    if inside + nearby == 0 {
        return Err(EngineError::Params(format!(
            "none of the {} target nodes lies in the source mesh or within the tolerance of it; are the two in the same place and units?",
            targets.len()
        )));
    }
    let field = if stress {
        InitialField::Stress { values: mapped, source: label.clone() }
    } else {
        InitialField::Temperature { values: mapped.into_iter().map(|(n, v)| (n, v[0])).collect(), source: label.clone() }
    };
    let kind = field.kind();
    s.model.set_initial_field(&field_name, field)?;
    let mut warnings = Vec::new();
    if unmapped > 0 {
        warnings
            .push(format!("{unmapped} of the {} target nodes are farther than the tolerance from the source mesh and got no value", targets.len()));
    }
    Ok(json!({
        "name": field_name,
        "quantity": kind,
        "source": label,
        "inside": inside,
        "nearby": nearby,
        "unmapped": unmapped,
        "max_distance": farthest,
        "tolerance": tolerance,
        "warnings": warnings,
    }))
}

/// At most 64 characters of a name, for messages.
fn short(name: &str) -> String {
    name.chars().take(64).collect()
}

fn not_implicit(step: &str) -> EngineError {
    EngineError::Params(format!(
        "not supported yet: mapping from an explicit or fluid step (`{}`); map from a static, heat transfer or coupled step",
        short(step)
    ))
}

/// The nodal stresses (`stress`) or temperatures in a step's result, as a field with no source.
fn values(r: &Implicit, stress: bool, frame: Option<usize>) -> Result<InitialField> {
    let what = if stress { "stresses" } else { "temperatures" };
    let wrong = |is: &str| EngineError::invalid("step", format!("step `{}` is {is}; it has no nodal {what}", short(r.step())));
    let source = String::new();
    match (stress, r) {
        (false, Implicit::Heat(h)) => {
            let at = match frame {
                Some(f) => h.frames.get(f).ok_or_else(|| {
                    EngineError::invalid("frame", format!("step `{}` has frames 0 to {}", short(&h.step), h.frames.len().saturating_sub(1)))
                })?,
                None => h.frames.last().ok_or_else(|| wrong("empty"))?,
            };
            Ok(InitialField::Temperature { values: h.nodes.iter().copied().zip(at.iter().copied()).collect(), source })
        }
        (false, Implicit::Static(st)) => {
            let t = st.temperature.as_ref().ok_or_else(|| wrong("a static step without temperatures"))?;
            refuse_frame(frame)?;
            Ok(InitialField::Temperature { values: st.nodes.iter().copied().zip(t.iter().copied()).collect(), source })
        }
        (true, Implicit::Static(st)) => {
            refuse_frame(frame)?;
            Ok(InitialField::Stress { values: st.nodes.iter().copied().zip(st.stress.iter().copied()).collect(), source })
        }
        (_, Implicit::Heat(_)) => Err(wrong("a heat transfer step")),
        (_, Implicit::Frequency(_)) => Err(wrong("a frequency or buckling step")),
        (_, Implicit::Harmonic(_)) => Err(wrong("a harmonic phasor step; initial fields require real static values")),
    }
}

fn refuse_frame(frame: Option<usize>) -> Result<()> {
    match frame {
        Some(_) => Err(EngineError::invalid("frame", "only a heat transfer step has frames to choose from")),
        None => Ok(()),
    }
}

/// The nodes of `set`, or every node of a solid element, ascending.
fn target_nodes(model: &Model, set: Option<&str>) -> Result<Vec<u32>> {
    let nodes: Vec<u32> = match set {
        Some(name) => model.resolve_set(name).map_err(|e| EngineError::invalid("set", e.to_string()))?.nodes,
        None => {
            let mut all: Vec<u32> = model.elements().values().filter(|e| e.kind.dimension() == 3).flat_map(|e| e.nodes.iter().copied()).collect();
            all.sort_unstable();
            all.dedup();
            all
        }
    };
    if nodes.is_empty() {
        return Err(EngineError::invalid("set", "there are no target nodes (the model has no solid elements, or the set no nodes)"));
    }
    Ok(nodes)
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use crate::{Session, StepResult};

    const L: f64 = 0.4;
    const W: f64 = 0.1;
    const E: f64 = 2e11;
    const NU: f64 = 0.3;
    const ALPHA: f64 = 1.2e-5;
    const T_REF: f64 = 293.15;

    fn scratch(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("fieldcraft-initial-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// A bar [0, L] × [0, W] × [0, W] in SI meshed with `kind` in `divisions`, its material,
    /// a solid section and its six faces named f1…f6 (f6 at x = 0, f4 at x = L).
    fn bar(s: &mut Session, kind: &str, divisions: [u32; 3]) {
        bar_of(s, kind, divisions, L);
    }

    fn bar_of(s: &mut Session, kind: &str, divisions: [u32; 3], l: f64) {
        let corners = [[0.0, 0.0, 0.0], [l, 0.0, 0.0], [l, W, 0.0], [0.0, W, 0.0], [0.0, 0.0, W], [l, 0.0, W], [l, W, W], [0.0, W, W]];
        let faces = json!({"1": "f1", "2": "f2", "3": "f3", "4": "f4", "5": "f5", "6": "f6"});
        let face = |f: &str| json!({"name": f, "geometry": [{"body": "Bar", "face": f}]});
        s.run_script(&json!([
            {"command": "model.units", "params": {"system": "m-kg-s-N-Pa"}},
            {"command": "mesh.block", "params": {"body": "Bar", "kind": kind, "corners": corners, "divisions": divisions, "faces": faces}},
            {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Bar"}]}},
            {"command": "model.set.create", "params": face("f6")},
            {"command": "model.set.create", "params": face("f4")},
            {"command": "model.set.create", "params": {"name": "skin", "geometry": (1..=6).map(|i| json!({"body": "Bar", "face": format!("f{i}")})).collect::<Vec<_>>()}},
            {"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": E, "poisson": NU},
                "thermal": {"conductivity": 45, "expansion": ALPHA, "reference_temperature": T_REF}}},
            {"command": "model.section.create", "params": {"name": "solid", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {}}},
        ]))
        .unwrap();
    }

    /// The source: a hex20 bar conducting between 300 K at x = 0 and 400 K at x = L (a linear
    /// temperature), and the same bar clamped all round and warmed uniformly to 400 K (a
    /// uniform hydrostatic stress −E α ΔT / (1 − 2ν)), both solved and saved to `path`.
    fn source(path: &str) -> Session {
        let mut s = Session::default();
        bar(&mut s, "hex20", [4, 2, 2]);
        s.run_script(&json!([
            {"command": "model.boundary_condition.create", "params": {"name": "Cold", "type": "temperature", "target": {"set": "f6"}, "value": 300}},
            {"command": "model.boundary_condition.create", "params": {"name": "Hot", "type": "temperature", "target": {"set": "f4"}, "value": 400}},
            {"command": "model.step.create", "params": {"name": "Heat", "procedure": {"type": "heat_transfer"}, "boundary_conditions": ["Cold", "Hot"]}},
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "skin"}}},
            {"command": "model.load.create", "params": {"name": "Warm", "type": "temperature", "target": {"set": "all"}, "value": 400}},
            {"command": "model.step.create", "params": {"name": "Squeeze", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Warm"]}},
        ]))
        .unwrap();
        s.execute("solve.heat", &json!({"step": "Heat"})).unwrap();
        s.execute("solve.static", &json!({"step": "Squeeze"})).unwrap();
        s.execute("file.save", &json!({"path": path})).unwrap();
        s
    }

    fn squeeze() -> f64 {
        -E * ALPHA * (400.0 - T_REF) / (1.0 - 2.0 * NU)
    }

    #[test]
    fn uniform_and_linear_fields_map_exactly_onto_another_mesh() {
        let dir = scratch("exact");
        let file = dir.join("source.fieldcraft").to_string_lossy().into_owned();
        source(&file);
        // The target: the same bar in tet10, with other divisions (no node in common but the
        // corners).
        let mut t = Session::default();
        bar(&mut t, "tet10", [3, 3, 1]);
        let r = t.execute("model.initial_field.map", &json!({"name": "Ejected", "quantity": "temperature", "file": file, "step": "Heat"})).unwrap();
        let nodes = t.model.nodes().len();
        assert_eq!((r["inside"].as_u64(), r["nearby"].as_u64(), r["unmapped"].as_u64()), (Some(nodes as u64), Some(0), Some(0)), "{r}");
        let Some(fieldcraft_model::InitialField::Temperature { values, source }) = t.model.initial_field("Ejected").cloned() else { panic!() };
        assert!(source.ends_with("source.fieldcraft, step Heat"), "{source}");
        let mut worst = 0.0f64;
        for (label, v) in &values {
            let x = t.model.node(*label).unwrap().xyz[0];
            worst = worst.max((v - (300.0 + 100.0 * x / L)).abs());
        }
        assert!(worst < 1e-9, "linear temperature: off by {worst:e} K");
        t.execute("model.initial_field.map", &json!({"name": "Residual", "quantity": "stress", "file": file, "step": "Squeeze"})).unwrap();
        let Some(fieldcraft_model::InitialField::Stress { values, .. }) = t.model.initial_field("Residual").cloned() else { panic!() };
        assert_eq!(values.len(), nodes);
        let want = [squeeze(), squeeze(), squeeze(), 0.0, 0.0, 0.0];
        let worst = values.values().flat_map(|s| s.iter().zip(want).map(|(a, b)| (a - b).abs())).fold(0.0f64, f64::max);
        assert!(worst < 1e-7 * squeeze().abs(), "uniform stress: off by {worst:e} Pa");
        // A static step of the target starts from it: clamped all round, nothing moves and the
        // stress is the mapped one.
        t.run_script(&json!([
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "skin"}}},
            {"command": "model.step.create", "params": {"name": "Hold", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"]}},
            {"command": "model.initial_field.delete", "params": {"name": "Ejected"}},
        ]))
        .unwrap();
        t.execute("solve.static", &json!({"step": "Hold"})).unwrap();
        let Some(StepResult::Implicit(fieldcraft_solver::StepResult::Static(st))) = t.result("Hold").map(|r| r.as_ref().clone()) else { panic!() };
        let worst = st.stress.iter().flat_map(|s| s.iter().zip(want).map(|(a, b)| (a - b).abs())).fold(0.0f64, f64::max);
        assert!(worst < 1e-6 * squeeze().abs(), "solved stress: off by {worst:e} Pa");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_field_maps_from_this_model_and_nodes_outside_are_counted() {
        let dir = scratch("self");
        let file = dir.join("source.fieldcraft").to_string_lossy().into_owned();
        let mut s = source(&file);
        let r = s.execute("model.initial_field.map", &json!({"name": "T0", "quantity": "temperature", "step": "Heat"})).unwrap();
        assert_eq!((r["unmapped"].as_u64(), r["nearby"].as_u64()), (Some(0), Some(0)), "{r}");
        assert_eq!(s.model.initial_field("T0").map(|f| f.source()), Some("step Heat"));
        // The heat step takes no initial fields, so its result stays current.
        assert!(s.result("Heat").is_some());
        // A target 2 % longer than the source: its far end is 8 mm outside, within the default
        // tolerance (a quarter of the source's 100 mm elements) at the closest point's 400 K.
        let mut t = Session::default();
        bar_of(&mut t, "hex8", [5, 1, 1], 1.02 * L);
        let r = t.execute("model.initial_field.map", &json!({"name": "T0", "quantity": "temperature", "file": file, "step": "Heat"})).unwrap();
        assert_eq!((r["nearby"].as_u64(), r["unmapped"].as_u64()), (Some(4), Some(0)), "{r}");
        assert!((r["max_distance"].as_f64().unwrap() - 0.02 * L).abs() < 1e-12, "{r}");
        assert!((r["tolerance"].as_f64().unwrap() - 0.025).abs() < 1e-12, "{r}");
        let Some(fieldcraft_model::InitialField::Temperature { values, .. }) = t.model.initial_field("T0") else { panic!() };
        assert!(values.iter().all(|(n, v)| t.model.node(*n).unwrap().xyz[0] < L || (v - 400.0).abs() < 1e-9));
        // Within a tighter tolerance those nodes get no value, and the reply says so.
        let r = t
            .execute("model.initial_field.map", &json!({"name": "T0", "quantity": "temperature", "file": file, "step": "Heat", "tolerance": 0.001}))
            .unwrap();
        assert_eq!((r["nearby"].as_u64(), r["unmapped"].as_u64()), (Some(0), Some(4)), "{r}");
        assert!(r["warnings"][0].as_str().unwrap().contains("4 of the 24 target nodes"), "{r}");
        assert_eq!(t.model.initial_field("T0").map(|f| f.len()), Some(20));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn hostile_parameters_are_refused() {
        let dir = scratch("hostile");
        let file = dir.join("source.fieldcraft").to_string_lossy().into_owned();
        let mut s = source(&file);
        let missing = dir.join("missing.fieldcraft").to_string_lossy().into_owned();
        for (p, why) in [
            (json!({"name": "F", "step": "Heat"}), "quantity"),
            (json!({"name": "F", "quantity": "pressure", "step": "Heat"}), "quantity"),
            (json!({"name": "F", "quantity": "temperature"}), "step"),
            (json!({"name": "F", "quantity": "temperature", "step": "Nope"}), "no current result"),
            (json!({"name": "F", "quantity": "stress", "step": "Heat"}), "no nodal stresses"),
            (json!({"name": "F", "quantity": "temperature", "step": "Squeeze"}), "no nodal temperatures"),
            (json!({"name": "F", "quantity": "temperature", "step": "Heat", "frame": 5}), "frame"),
            (json!({"name": "F", "quantity": "stress", "step": "Squeeze", "frame": 0}), "frame"),
            (json!({"name": "F", "quantity": "temperature", "step": "Heat", "tolerance": -1}), "tolerance"),
            (json!({"name": "F", "quantity": "temperature", "step": "Heat", "set": "nope"}), "set"),
            (json!({"name": "F", "quantity": "temperature", "step": "Heat", "file": missing}), "missing.fieldcraft"),
            (json!({"name": "F", "quantity": "temperature", "step": "Nope", "file": file}), "keeps no results"),
            (json!({"name": "1 bad", "quantity": "temperature", "step": "Heat"}), "name"),
            (json!({"name": "F", "quantity": "temperature", "step": "Heat", "frame": -3}), "frame"),
            (json!({"name": "F", "quantity": "temperature", "step": "Heat", "tolerance": "x"}), "tolerance"),
        ] {
            let e = s.execute("model.initial_field.map", &p).unwrap_err().to_string();
            assert!(e.contains(why), "{p}: {e}");
        }
        // A target set far from the source (the edit leaves this model's results stale, so map
        // from the file).
        s.execute("mesh.node.add", &json!({"label": 99_999, "xyz": [5.0, 5.0, 5.0]})).unwrap();
        s.execute("model.set.create", &json!({"name": "far", "nodes": [99_999]})).unwrap();
        let p = json!({"name": "F", "quantity": "temperature", "file": file, "step": "Heat", "set": "far"});
        let e = s.execute("model.initial_field.map", &p).unwrap_err().to_string();
        assert!(e.contains("none of the 1 target nodes"), "{e}");
        assert!(s.model.initial_fields().is_empty());
        assert!(s.execute("model.initial_field.delete", &json!({"name": "F"})).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
