//! Stress-life fatigue of a solved linear static step: its fatigue assessment
//! ([`fieldcraft_model::FatigueAssessment`], `model.fatigue.set`) read at every node with the
//! fatigue data of the materials around it (`fieldcraft_fem::fatigue` holds the method and its
//! sources).
//!
//! Nothing is stored: the fields are computed from the step's saved result and the saved
//! definitions whenever they are asked for (contours, probes, path plots, `results.fatigue`), so
//! they survive save and open with no change to the results encoding, and follow an edit of the
//! assessment or the materials' fatigue data at once, without solving again.
//!
//! At each node the reference stress is the step's nodal stress (averaged over the elements
//! around the node, as the step stores it). A node on the boundary of two materials with fatigue
//! data takes the worse of the two (larger damage, smaller safety factor); a node of no element
//! whose material has fatigue data has no value.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_fem::fatigue::{self as method, Cycle};
use fieldcraft_model::fatigue::Counting;
use fieldcraft_model::{Fatigue, FatigueAssessment, Loading, Model, SectionKind};
use fieldcraft_solver::StaticResult;

use crate::results::Kind;
use crate::{EngineError, Result, StepResult};

/// The most node-cycle evaluations one assessment makes (nodes times distinct cycles): a cap on
/// the work a hostile history can ask for.
pub const MAX_WORK: usize = 500_000_000;

/// The fatigue at one node.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeFatigue {
    /// The worst of the node's materials.
    pub point: method::Point,
    /// The material that point is of.
    pub material: String,
    /// The reference stress measure (the step's stress, before the load factors).
    pub reference: f64,
    /// The step's von Mises stress times the largest load factor exceeds the yield strength:
    /// the stress is not elastic there, and the stress-life method does not apply.
    pub above_yield: bool,
}

/// A fatigue assessment evaluated over a result.
#[derive(Debug, Clone)]
pub struct Evaluation {
    /// The step's result nodes ([`StepResult::nodes`] order), each with its fatigue or `None`.
    pub nodes: Vec<(u32, Option<NodeFatigue>)>,
    /// The load factor cycles counted, combined (range and mean in load factors).
    pub cycles: Vec<Cycle>,
    /// Whether the loading is one constant amplitude cycle (life in cycles, else in blocks).
    pub constant: bool,
    /// Whether some material in use has an endurance limit (the safety factor is defined).
    pub endurance: bool,
    /// The history is counted once: its damage is of one pass, and there is no life (a block
    /// that repeats is what a life counts).
    pub once: bool,
    /// Every material in use has an ultimate strength, so the peak stress was checked against it
    /// everywhere; otherwise (some or none) the warnings say where no static check was made.
    pub static_check: bool,
    /// What the user should know about where the method may not apply.
    pub warnings: Vec<String>,
}

/// The assessment of `step` in `model`, or why there is none.
pub fn assessment<'a>(model: &'a Model, step: &str) -> Result<&'a FatigueAssessment> {
    model.fatigue(step).ok_or_else(|| EngineError::Params(format!("step `{step}` has no fatigue assessment: define one with model.fatigue.set")))
}

/// The linear static result of `r`, or why fatigue cannot read it.
pub fn linear_static(r: &StepResult) -> Result<&StaticResult> {
    match r.kind() {
        Kind::Static(s) if s.nonlinear.is_none() => Ok(s),
        Kind::Static(s) => Err(EngineError::Params(format!(
            "step `{}` is nonlinear: its stress is not proportional to its load, so its cycles cannot be scaled from it",
            s.step
        ))),
        _ => Err(EngineError::Params(format!("step `{}` is not a linear static step: fatigue reads a linear static stress field", r.step()))),
    }
}

/// The load factor cycles of `loading`, combined.
pub fn cycles(loading: &Loading) -> Vec<Cycle> {
    match loading {
        Loading::History { factors, counting: Counting::Repeated } => method::rainflow_repeated(factors),
        Loading::History { factors, counting: Counting::Once } => method::rainflow(factors),
        _ => {
            let (max, min) = loading.constant().unwrap_or((0.0, 0.0));
            method::combine(vec![Cycle { range: max - min, mean: 0.5 * (max + min), count: 1.0 }])
        }
    }
}

/// The materials with fatigue data around each node of `model`'s mesh, by node label (each list
/// sorted and without repeats), and the warnings about the sections they come from.
fn node_materials(model: &Model) -> Result<(BTreeMap<u32, Vec<&str>>, Vec<String>)> {
    let mut out: BTreeMap<u32, BTreeSet<&str>> = BTreeMap::new();
    let mut warnings = Vec::new();
    for (name, section) in model.sections() {
        let Some(material) = model.material(&section.material) else { continue };
        if material.fatigue.is_none() {
            continue;
        }
        match section.kind {
            SectionKind::Shell { .. } => {
                warnings.push(format!("section `{name}` is a shell: its nodal stress is its top surface's, so the other surface is not assessed"))
            }
            SectionKind::Beam { .. } | SectionKind::Rod { .. } => {
                warnings.push(format!("section `{name}` is a beam or truss: its nodes carry no stress tensor, and are not assessed"));
                continue;
            }
            _ => {}
        }
        let set = model.resolve_set(&section.set)?;
        for e in set.elements {
            for n in model.element(e).map(|e| e.nodes.as_slice()).unwrap_or_default() {
                out.entry(*n).or_default().insert(section.material.as_str());
            }
        }
    }
    Ok((out.into_iter().map(|(n, m)| (n, m.into_iter().collect())).collect(), warnings))
}

/// The fatigue assessment of `r`'s step evaluated at its nodes (all of them, or those at the
/// indices `only` into [`StepResult::nodes`]; the others get `None`).
pub fn evaluate(model: &Model, r: &StepResult, only: Option<&[usize]>) -> Result<Evaluation> {
    let a = assessment(model, r.step())?;
    let s = linear_static(r)?;
    if s.stress.len() != s.nodes.len() {
        return Err(EngineError::Params(format!("step `{}` kept no nodal stresses to assess", s.step)));
    }
    let (around, mut warnings) = node_materials(model)?;
    let used: BTreeSet<&str> = around.values().flatten().copied().collect();
    if used.is_empty() {
        return Err(EngineError::Params("no section's material has fatigue data: give the material an S-N curve with model.material.fatigue".into()));
    }
    let data: BTreeMap<&str, &Fatigue> = used.iter().filter_map(|name| Some((*name, model.material(name)?.fatigue.as_ref()?))).collect();
    for (name, f) in &data {
        method::check(f, a.mean_stress).map_err(|e| EngineError::Params(format!("material `{name}`: {e} (model.material.fatigue)")))?;
    }
    for (name, f) in &data {
        if f.ultimate_strength.is_none() {
            warnings.push(format!("material `{name}` has no ultimate strength: no static check of the peak stress was made for it"));
        }
    }
    let static_check = data.values().all(|f| f.ultimate_strength.is_some());
    let cycles = cycles(&a.loading);
    let wanted: Vec<usize> = match only {
        Some(k) => k.to_vec(),
        None => (0..s.nodes.len()).collect(),
    };
    if wanted.len().saturating_mul(cycles.len().max(1)) > MAX_WORK {
        return Err(EngineError::Params(format!(
            "{} nodes times {} distinct cycles is more than {MAX_WORK} evaluations: shorten the history or assess fewer nodes",
            wanted.len(),
            cycles.len()
        )));
    }
    let peak = a.loading.peak();
    let extent = a.loading.extent();
    let mut nodes: Vec<(u32, Option<NodeFatigue>)> = s.nodes.iter().map(|n| (*n, None)).collect();
    for (done, i) in wanted.into_iter().enumerate() {
        if done % 2048 == 0 {
            fieldcraft_model::progress::check()?;
        }
        let (Some((label, slot)), Some(t)) = (nodes.get_mut(i).map(|(l, v)| (*l, v)), s.stress.get(i)) else { continue };
        let Some(materials) = around.get(&label) else { continue };
        if t.iter().any(|x| !x.is_finite()) {
            continue;
        }
        let m = method::measure(a.measure, t);
        let vm = fieldcraft_fem::von_mises(t) * peak;
        let mut worst: Option<NodeFatigue> = None;
        for name in materials {
            let Some(f) = data.get(name) else { continue };
            let p = method::assess(f, a.mean_stress, a.stress_factor, m, &cycles, extent)
                .map_err(|e| EngineError::Params(format!("material `{name}`: {e}")))?;
            let above_yield = f.yield_strength.is_some_and(|y| vm > y);
            let node = NodeFatigue { point: p, material: name.to_string(), reference: m, above_yield };
            worst = Some(match worst {
                Some(w) => {
                    let sf = match (w.point.safety_factor, p.safety_factor) {
                        (Some(a), Some(b)) => Some(a.min(b)),
                        (a, b) => a.or(b),
                    };
                    let mut keep = if p.damage > w.point.damage { node } else { w };
                    keep.point.safety_factor = sf;
                    keep.above_yield |= above_yield;
                    keep
                }
                None => node,
            });
        }
        *slot = worst;
    }
    let endurance = data.values().any(|f| f.endurance_limit.is_some());
    let once = matches!(a.loading, Loading::History { counting: Counting::Once, .. });
    Ok(Evaluation { nodes, constant: a.loading.constant().is_some(), cycles, endurance, once, static_check, warnings })
}

/// Why a history counted once has no life field.
pub const ONCE_HAS_NO_LIFE: &str = "a history counted once has damage (fatigue_damage) but no life: its ranges left at the end are counted as half cycles, so 1/damage is not a number of repeats of the history; use counting `repeated` for a life in blocks";

/// Which fatigue field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Quantity {
    /// log₁₀ of the life (cycles, or blocks of a history); infinite below the endurance limit.
    LogLife,
    /// Palmgren–Miner damage of one cycle or block.
    Damage,
    /// The safety factor against the endurance limit.
    SafetyFactor,
}

/// A fatigue field at every node of `r` ([`StepResult::nodes`] order): not a number where a node
/// has no fatigue data (or, for the safety factor, no endurance limit).
pub fn field(model: &Model, r: &StepResult, quantity: Quantity) -> Result<Vec<f64>> {
    let e = evaluate(model, r, None)?;
    if quantity == Quantity::LogLife && e.once {
        return Err(EngineError::Params(ONCE_HAS_NO_LIFE.into()));
    }
    if quantity == Quantity::SafetyFactor && !e.endurance {
        return Err(EngineError::Params("the safety factor needs an endurance limit: no material in use gives one (model.material.fatigue)".into()));
    }
    Ok(e.nodes
        .iter()
        .map(|(_, v)| match (v, quantity) {
            (None, _) => f64::NAN,
            (Some(n), Quantity::LogLife) => n.point.log10_life(),
            (Some(n), Quantity::Damage) => n.point.damage,
            (Some(n), Quantity::SafetyFactor) => n.point.safety_factor.unwrap_or(f64::NAN),
        })
        .collect())
}
