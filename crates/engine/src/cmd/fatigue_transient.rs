//! Fatigue of an exported, explicitly identified SI nodal stress history.
//! No current-session mesh is inferred from labels in an external probe reply.
use super::fatigue::field;
use crate::schema::{Field, json as field_json};
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};
use fieldcraft_solver::fatigue::{Basquin, MAX_HISTORY_POINTS, life::StressMean, miner, rainflow};
use serde_json::{Value, json};
use std::collections::BTreeSet;

const FORM: &[Field] = &[
    field_json("history", "Transient probe reply").req().help("Original solve.probe reply: SI seconds and Pa, global nodal stress, spectrum false."),
    field_json("source_id", "Source identity")
        .req()
        .help("Caller-verified source model and raw history identity; retained, not associated with the current mesh."),
    field_json("step", "Expected source step").req().help("Must equal the exported step name."),
    field_json("nodes", "Expected source nodes").req().help("Exact set of exported node labels, unique and positive."),
    field_json("max_time_step_seconds", "Maximum sample gap (s)")
        .req()
        .help("Explicit sampling contract; no interpolation or automatic peak recovery."),
    field_json("weights", "Stress projection").req().help("Fixed global [xx,yy,zz,xy,xz,yz] weights applied before rainflow."),
    field_json("curve", "Basquin curve (SI)").req().help("amplitude_ref Pa, cycles_ref cycles, exponent, endurance_limit Pa."),
    field_json("correction", "Mean stress correction").req().help("Tagged type none/goodman/gerber/soderberg with strengths in Pa."),
];
pub const TRANSIENT: CommandSpec = CommandSpec {
    id: "solve.fatigue.transient",
    label: "Transient Stress Fatigue",
    params: "{history,source_id,step,nodes,max_time_step_seconds,weights,curve,correction} — validate an exported SI transient nodal stress history, project, rainflow count and accumulate Miner damage; source labels never refer to the current session mesh",
    form: Some(FORM),
    kind: CommandKind::Model,
    requires: &[],
    run,
};
fn bad(why: &str) -> EngineError {
    EngineError::invalid("history", why)
}
fn run(_: &mut Session, p: &Value) -> Result<Value> {
    let history: Value = field(p, "history")?;
    let source_id: String = field(p, "source_id")?;
    let step: String = field(p, "step")?;
    let nodes: Vec<u32> = field(p, "nodes")?;
    let gap: f64 = field(p, "max_time_step_seconds")?;
    let weights: [f64; 6] = field(p, "weights")?;
    let curve: Basquin = field(p, "curve")?;
    let correction: StressMean = field(p, "correction")?;
    if source_id.trim().is_empty() || source_id.len() > 1024 {
        return Err(EngineError::invalid("source_id", "need a nonempty source identity of at most 1024 bytes"));
    }
    if step.is_empty() || history.get("step").and_then(Value::as_str) != Some(step.as_str()) {
        return Err(EngineError::invalid("step", "source step differs"));
    }
    let expected: BTreeSet<_> = nodes.iter().copied().collect();
    if nodes.is_empty() || nodes.len() > 64 || expected.len() != nodes.len() || expected.contains(&0) {
        return Err(EngineError::invalid("nodes", "need 1..64 unique positive source labels"));
    }
    if !gap.is_finite() || gap <= 0.0 {
        return Err(EngineError::invalid("max_time_step_seconds", "need a finite positive sample gap"));
    }
    if weights.iter().any(|x| !x.is_finite()) || weights.iter().all(|x| *x == 0.0) {
        return Err(EngineError::invalid("weights", "need a finite nonzero projection"));
    }
    curve.validate()?;
    correction.corrected_amplitude(0.0, 0.0)?;
    if history.get("spectrum").and_then(Value::as_bool) != Some(false) {
        return Err(bad("response spectra and missing spectrum declarations are not stress histories"));
    }
    if history.pointer("/units/time").and_then(Value::as_str) != Some("s") || history.pointer("/units/stress").and_then(Value::as_str) != Some("Pa") {
        return Err(bad("export must explicitly use seconds and Pa"));
    }
    if history.get("stress_components") != Some(&json!(["xx", "yy", "zz", "xy", "xz", "yz"])) {
        return Err(bad("missing or different stress component order"));
    }
    if history.get("stress_location").and_then(Value::as_str) != Some("global nodal average; shell top surface") {
        return Err(bad("unsupported or missing stress recovery location"));
    }
    let frames = history.get("frames").and_then(Value::as_array).ok_or_else(|| bad("missing frame time axis"))?;
    if frames.len() < 2 || frames.len() > MAX_HISTORY_POINTS / nodes.len() {
        return Err(bad("need at least two frames and at most one million node samples"));
    }
    let mut times = Vec::with_capacity(frames.len());
    for f in frames {
        let t = f.get("time").and_then(Value::as_f64).filter(|x| x.is_finite()).ok_or_else(|| bad("missing or nonfinite frame time"))?;
        if let Some(&previous) = times.last() {
            let dt: f64 = t - previous;
            if !dt.is_finite() || dt <= 0.0 || dt > gap * (1.0 + 1e-10) {
                return Err(bad("frame times must increase within the declared maximum gap"));
            }
        }
        times.push(t);
    }
    let duration = times.last().zip(times.first()).map(|(b, a)| b - a).ok_or_else(|| bad("missing duration"))?;
    if !duration.is_finite() {
        return Err(bad("duration overflow"));
    }
    let rows = history.get("nodes").and_then(Value::as_array).ok_or_else(|| bad("missing nodes"))?;
    if rows.len() != nodes.len() {
        return Err(bad("source node count differs"));
    }
    let mut seen = BTreeSet::new();
    let mut result = Vec::new();
    for row in rows {
        let label = row.get("node").and_then(Value::as_u64).and_then(|x| u32::try_from(x).ok()).ok_or_else(|| bad("invalid node label"))?;
        if !expected.contains(&label) || !seen.insert(label) {
            return Err(bad("source node identity differs or is duplicated"));
        }
        let samples = row.get("samples").and_then(Value::as_array).ok_or_else(|| bad("missing samples"))?;
        if samples.len() != times.len() {
            return Err(bad("node sample count differs from frame count"));
        }
        let mut projected = Vec::with_capacity(samples.len());
        for (sample, &time) in samples.iter().zip(&times) {
            if sample.get("time").and_then(Value::as_f64) != Some(time) {
                return Err(bad("node and frame times differ"));
            }
            let stress: [f64; 6] = serde_json::from_value(sample.get("stress").cloned().ok_or_else(|| bad("missing stress"))?)
                .map_err(|_| bad("need six finite stress components; null is not zero"))?;
            if stress.iter().any(|x| !x.is_finite()) {
                return Err(bad("nonfinite stress"));
            }
            let value: f64 = stress.iter().zip(weights).map(|(s, w)| s * w).sum();
            if !value.is_finite() {
                return Err(bad("stress projection overflow"));
            }
            projected.push(value);
        }
        let cycles = rainflow(&projected)?;
        let damage = miner(&cycles, &curve, correction)?;
        let mut cumulative = 0.0;
        let counted: Vec<_> = cycles.iter().map(|c| -> Result<Value> {
            let increment = miner(std::slice::from_ref(c), &curve, correction)?;
            cumulative += increment;
            Ok(json!({"range_pa":c.range,"amplitude_pa":c.range*0.5,"mean_pa":c.mean,"count":c.count,"damage":increment,"cumulative_damage":cumulative}))
        }).collect::<Result<_>>()?;
        let life = if damage > 0.0 {
            let x = duration / damage;
            x.is_finite().then_some(x)
        } else {
            None
        };
        result.push(json!({"node":label,"samples":samples.len(),"cycles":counted,"cycle_count":cycles.iter().map(|c|c.count).sum::<f64>(),"damage":damage,"duration_seconds":duration,"repeated_record_life_seconds":life}));
    }
    Ok(
        json!({"source_id":source_id,"source_step":step,"nodes":result,"weights":weights,"units":{"stress":"Pa","time":"s"},"method":"exported_transient_projection_rainflow_miner","limitations":"source identity is caller-verified; no current-mesh binding; global nodal stresses only; sampling adequacy not certified; finite record residual half cycles, no periodic closure; cumulative damage is in counting order, not chronological damage history"}),
    )
}
