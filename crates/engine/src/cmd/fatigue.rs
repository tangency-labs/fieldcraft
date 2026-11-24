//! Stress-life fatigue: `model.fatigue.set` and `model.fatigue.delete` define how a linear static
//! step's load repeats, `results.fatigue` reports the life, damage and safety factor that gives
//! with the materials' S-N data (`model.material.fatigue`), computed from the step's result
//! ([`crate::fatigue`]). The fields themselves are results fields (`fatigue_log_life`,
//! `fatigue_damage`, `fatigue_safety_factor`) for results.contour, solve.probe and results.path.

use serde_json::{Value, json};

use super::path::current_result;
use super::results::{Want, step_for};
use super::{opt, req};
use crate::fatigue::{Evaluation, NodeFatigue, evaluate};
use crate::requirement::Requirement;
use crate::results::{Field, Quantity, shown, shown_all};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_model::FatigueAssessment;

pub const SET: CommandSpec = CommandSpec {
    id: "model.fatigue.set",
    label: "Fatigue Loading",
    params: r#"{"step": "Load", "loading": {"type": "constant_amplitude", "max": 1, "min": -1} | {"type": "stress_ratio", "ratio": 0.1, "scale"?: 1} | {"type": "history", "factors": [0, 1, -0.5, 0.8, …], "counting"?: "repeated" | "once"}, "measure": "signed_von_mises" | "max_abs_principal", "mean_stress": "none" | "goodman" | "gerber" | "soderberg" | "swt", "stress_factor"?: 1} — gives a linear static step its fatigue assessment (or replaces the one it has): the step's stress field is the response to its load, and the repeated load is that field times load factors (proportional loading: every stress component follows the same factor). `constant_amplitude` cycles between `min` and `max` times the step's load; `stress_ratio` from `scale` to R·scale (R = σ_min/σ_max, from -1 to 1); `history` is one block of load factors (2 to 100000), rainflow counted (the three-point procedure of ASTM E1049-85) as a block that repeats (`repeated`, the default: every cycle closes, Downing & Socie 1982; the life is in blocks) or happens once (`once`: the ranges left at the end count as half cycles; only its damage is given, no life), with Palmgren–Miner damage. `measure` is the scalar stress read at each node: the von Mises stress signed as the principal stress of largest magnitude (ductile metals), or that principal stress itself (brittle materials). `mean_stress` corrects the amplitude for a tensile mean (Goodman and Gerber need the materials' ultimate strength, Soderberg their yield strength; those three give a compressive mean no correction; `swt`, Smith, Watson and Topper 1970, reads the curve at √(σ_max σ_a), needs no strength, and does no damage where σ_max is not tensile). `stress_factor` (above 0, at most 100) multiplies the whole stress history, amplitude and mean alike: a stress concentration the mesh does not hold, or 1/k for a surface or size factor k (one number for both uses). Wherever a material has an ultimate strength, a peak stress at or above it is a static failure (life 0, safety factor at most 1) whatever the correction; a material without one gets no static check (a warning says so). The step's result stays current: no solve is needed. Fields: fatigue_log_life, fatigue_damage, fatigue_safety_factor (results.contour, solve.probe, results.path); summary: results.fatigue"#,
    form: Some(super::forms::FATIGUE_SET),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: set,
};

pub const DELETE: CommandSpec = CommandSpec {
    id: "model.fatigue.delete",
    label: "Delete Fatigue Loading",
    params: r#"{"step": "Load"} — removes a step's fatigue assessment"#,
    form: Some(super::forms::FATIGUE_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Step],
    run: delete,
};

pub const REPORT: CommandSpec = CommandSpec {
    id: "results.fatigue",
    label: "Fatigue",
    params: r#"{"step": "Load", "set"?: "fillet"} — the stress-life fatigue of a solved linear static step under its fatigue assessment (model.fatigue.set) and its materials' S-N data (model.material.fatigue), over every node with fatigue data or over a set's nodes: the counted load factor `cycles` (range, mean and count; `cycles_per_block`), the `life_unit` (cycles of a constant amplitude loading, blocks of a history that repeats, none for a history counted once, which has damage only) with its `life_note`, and `static_check` (every material in use has an ultimate strength), the node of least life (`min_life`: its life, log₁₀ life, damage per cycle or block, material, governing cycle's stress amplitude, mean and equivalent fully reversed amplitude, and safety factor), the node of least safety factor against the endurance limit (`min_safety_factor`, when a material has one), counts of nodes assessed, of infinite life, of static failures (an amplitude beyond the curve's end, a mean beyond the strength, or a peak stress at or above the ultimate strength), of those above the yield strength at the peak load (the stress-life method does not apply there) and of those read beyond a tabulated curve's points, and `warnings`. Stresses in the display units. A query: not journalled"#,
    form: Some(super::forms::RESULTS_FATIGUE),
    kind: CommandKind::Session,
    requires: &[Requirement::Result],
    run: report,
};

/// Most cycles a reply lists.
const LISTED_CYCLES: usize = 200;

fn set(s: &mut Session, p: &Value) -> Result<Value> {
    let step: String = req(p, "step")?;
    let mut fields = p.as_object().cloned().ok_or_else(|| EngineError::Params("the parameters are a JSON object".into()))?;
    fields.remove("step");
    let a: FatigueAssessment = serde_json::from_value(Value::Object(fields)).map_err(|e| EngineError::Params(format!("fatigue assessment: {e}")))?;
    s.model.set_fatigue(&step, a)?;
    Ok(json!({"step": step}))
}

fn delete(s: &mut Session, p: &Value) -> Result<Value> {
    let step: String = req(p, "step")?;
    s.model.remove_fatigue(&step)?;
    Ok(Value::Null)
}

/// A number for a reply: `"infinite"` when it is.
fn number(v: f64) -> Value {
    if v == f64::INFINITY { json!("infinite") } else { json!(v) }
}

/// A node's fatigue in a reply, its stresses as `stress` gives them. `life` is false for a
/// history counted once, which has damage but no life (null).
fn node_json(n: &NodeFatigue, stress_of: &dyn Fn(f64) -> f64, life: bool) -> Value {
    let stress = |v: f64| if v.is_finite() { json!(stress_of(v)) } else { number(v) };
    json!({
        "life": if life { number(n.point.life) } else { Value::Null },
        "log10_life": if life { number(n.point.log10_life()) } else { Value::Null },
        "damage": n.point.damage,
        "safety_factor": n.point.safety_factor.map(number),
        "material": n.material,
        "stress_amplitude": stress(n.point.amplitude),
        "mean_stress": stress(n.point.mean),
        "equivalent_amplitude": stress(n.point.equivalent),
        "static_failure": n.point.static_failure,
        "above_yield": n.above_yield,
        "extrapolated": n.point.extrapolated,
    })
}

fn report(s: &mut Session, p: &Value) -> Result<Value> {
    fieldcraft_model::progress::check()?;
    let step = step_for(s, p, Want::Any)?;
    let result = current_result(s, &step)?;
    let only: Option<Vec<usize>> = match opt::<String>(p, "set")? {
        Some(name) => {
            let set = s.model.resolve_set(&name)?;
            let nodes = result.nodes();
            Some(set.nodes.iter().filter_map(|n| nodes.binary_search(n).ok()).collect())
        }
        None => None,
    };
    let e = evaluate(&s.model, &result, only.as_deref())?;
    let display = s.units().clone();
    let a = crate::fatigue::assessment(&s.model, &step)?;
    Ok(summary(&s.model, &step, a, &e, &display))
}

/// Each probed node's fatigue (`fatigue`, its stresses in SI for the probe to show) in a probe
/// reply `out` of the static result `r`; why there is none in `fatigue_unavailable`.
pub(crate) fn add_to_probe(model: &fieldcraft_model::Model, r: &crate::StepResult, out: &mut Value) {
    let Some(entries) = out.get_mut("nodes").and_then(Value::as_array_mut) else { return };
    let labels = r.nodes();
    let at: Vec<Option<usize>> = entries
        .iter()
        .map(|e| e.get("node").and_then(Value::as_u64).and_then(|l| u32::try_from(l).ok()).and_then(|l| labels.binary_search(&l).ok()))
        .collect();
    let wanted: Vec<usize> = at.iter().flatten().copied().collect();
    match evaluate(model, r, Some(&wanted)) {
        Ok(e) => {
            for (entry, k) in entries.iter_mut().zip(at) {
                if let (Some(n), Some(o)) = (k.and_then(|k| e.nodes.get(k)).and_then(|(_, n)| n.as_ref()), entry.as_object_mut()) {
                    o.insert("fatigue".into(), node_json(n, &|x| x, !e.once));
                }
            }
        }
        Err(e) => {
            if let Some(o) = out.as_object_mut() {
                o.insert("fatigue_unavailable".into(), json!(e.to_string()));
            }
        }
    }
}

fn summary(model: &fieldcraft_model::Model, step: &str, a: &FatigueAssessment, e: &Evaluation, display: &fieldcraft_model::DisplaySystem) -> Value {
    let assessed: Vec<(u32, &NodeFatigue)> = e.nodes.iter().filter_map(|(l, n)| n.as_ref().map(|n| (*l, n))).collect();
    let count = |f: &dyn Fn(&NodeFatigue) -> bool| assessed.iter().filter(|(_, n)| f(n)).count();
    let xyz = |label: u32| model.node(label).map(|n| shown_all(Quantity::Length, n.xyz, display));
    // Least life (largest damage; the lowest label on a tie), then least safety factor.
    let worst = assessed.iter().filter(|(_, n)| n.point.damage > 0.0).fold(None, |best: Option<&(u32, &NodeFatigue)>, x| match best {
        Some(b) if b.1.point.damage >= x.1.point.damage => best,
        _ => Some(x),
    });
    let min_life = worst.map(|(l, n)| {
        let mut v = node_json(n, &|x| shown(Quantity::Pressure, x, display), !e.once);
        v["node"] = json!(l);
        v["xyz"] = json!(xyz(*l));
        v
    });
    let least_sf = assessed
        .iter()
        .filter_map(|(l, n)| n.point.safety_factor.map(|f| (*l, f)))
        .fold(None, |best: Option<(u32, f64)>, x| if best.is_some_and(|b| b.1 <= x.1) { best } else { Some(x) });
    let (infinite, failures, yielded, extrapolated) =
        (count(&|n| n.point.damage == 0.0), count(&|n| n.point.static_failure), count(&|n| n.above_yield), count(&|n| n.point.extrapolated));
    let mut warnings = e.warnings.clone();
    if yielded > 0 {
        warnings.push(format!(
            "{yielded} node(s) are above the yield strength at the peak load (von Mises stress times the largest load factor): the stress is not elastic there, and stress-life does not apply (strain-life would)"
        ));
    }
    if failures > 0 {
        warnings.push(format!(
            "{failures} node(s) are static failures (an amplitude at or beyond the curve's end, a mean at or past the strength of the correction, or a peak stress at or above the ultimate strength): life 0, safety factor at most 1, damage at least 2"
        ));
    }
    if extrapolated > 0 {
        warnings.push(format!("{extrapolated} node(s) read a tabulated S-N curve beyond its points (its end segment extended)"));
    }
    if assessed.is_empty() {
        warnings.push("no node of the selection has fatigue data".into());
    }
    let cycles: Vec<Value> = e.cycles.iter().take(LISTED_CYCLES).map(|c| json!({"range": c.range, "mean": c.mean, "count": c.count})).collect();
    let mut fields = vec![Field::FatigueDamage.id()];
    if !e.once {
        fields.insert(0, Field::FatigueLogLife.id());
    }
    if e.endurance {
        fields.push(Field::FatigueSafetyFactor.id());
    }
    json!({
        "step": step,
        "assessment": a,
        "life_unit": if e.constant { "cycles" } else if e.once { "none" } else { "blocks" },
        "life_note": if e.constant {
            "cycles of the constant amplitude loading"
        } else if e.once {
            "no life: the history is counted once, so only its damage (of one pass, with the ranges left at the end as half cycles) is given; use counting `repeated` for a life in blocks"
        } else {
            "blocks: repeats of the history as a block that repeats, every cycle closed (1 / damage)"
        },
        "static_check": e.static_check,
        "cycles": cycles,
        "distinct_cycles": e.cycles.len(),
        "cycles_per_block": e.cycles.iter().map(|c| c.count).sum::<f64>(),
        "nodes": {"assessed": assessed.len(), "infinite_life": infinite, "static_failures": failures, "above_yield": yielded, "extrapolated": extrapolated},
        "min_life": min_life,
        "min_safety_factor": least_sf.map(|(l, f)| json!({"node": l, "xyz": xyz(l), "value": number(f)})),
        "fields": fields,
        "warnings": warnings,
        "units": {"stress": crate::results::unit_of(Quantity::Pressure, display), "xyz": crate::results::unit_of(Quantity::Length, display)},
        "method": "stress-life (S-N) with proportional loading: the step's nodal stress times the load factors; linear (Palmgren–Miner) damage",
    })
}

#[cfg(test)]
mod tests;
