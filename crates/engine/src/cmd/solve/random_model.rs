//! Native modal random steps reuse harmonic FE recovery and the existing one-sided PSD
//! projection. A coherent spatial load is solved once per frequency before its square is
//! formed, retaining modal cross terms (Newland, Random Vibrations, 3rd ed.).
use super::random_frames;
use crate::{EngineError, Result, Session, StepResult};
use fieldcraft_model::harmonic::{Frequencies, Harmonic, Method};
use fieldcraft_model::random::Excitation;
use fieldcraft_model::{BoundaryCondition, Model, Procedure, Step, Target};
use fieldcraft_solver::{FrequencyOptions, StepResult as ImplicitResult};
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn constraints(model: &Model, step: &Step) -> Result<BTreeSet<(u32, usize)>> {
    let mut out = BTreeSet::new();
    for name in &step.boundary_conditions {
        let Some(BoundaryCondition::Displacement { target: Target::Set(set), values }) = model.boundary_conditions().get(name) else {
            return Err(EngineError::invalid("step", "modal random import currently requires explicit node-set zero displacement constraints"));
        };
        if values.iter().flatten().any(|v| *v != 0.0) {
            return Err(EngineError::invalid("step", "nonzero random prescribed motion is not supported"));
        }
        let nodes = model.set(set).ok_or_else(|| EngineError::invalid("step", "missing constraint set"))?;
        if !nodes.elements.is_empty() || !nodes.faces.is_empty() || !nodes.edges.is_empty() {
            return Err(EngineError::invalid("step", "random constraints require explicit node membership"));
        }
        for node in &nodes.nodes {
            for (dof, value) in values.iter().enumerate() {
                if value.is_some() {
                    out.insert((*node, dof));
                }
            }
        }
    }
    Ok(out)
}

pub(super) fn run(s: &mut Session, p: &Value, name: &str) -> Result<Value> {
    for key in [
        "input_psd",
        "base_direction",
        "frequencies",
        "transfer",
        "base_psd",
        "natural_frequency",
        "damping_ratio",
        "base_acceleration_psd",
        "stress_psd",
        "probe_indices",
    ] {
        if p.get(key).is_some() {
            return Err(EngineError::invalid(key, "native random steps use their stored excitation; external overrides are not supported"));
        }
    }
    let step = s.model.step(name).ok_or_else(|| EngineError::invalid("step", "missing random step"))?.clone();
    let Procedure::RandomVibration(random) = &step.procedure else {
        return Err(EngineError::invalid("step", "not a random step"));
    };
    random.validate()?;
    let modal = random.modal.as_ref().ok_or_else(|| EngineError::invalid("step", "random FE execution requires modal interval settings"))?;
    let [Excitation::CorrelatedLoads { spectrum }] = random.excitations.as_slice() else {
        return Err(EngineError::invalid("step", "native random execution supports one fully correlated spatial load pattern"));
    };
    if spectrum.windows(2).any(|w| w[0].psd != w[1].psd) {
        return Err(EngineError::invalid("step", "nonconstant multiplier PSD interpolation is not mapped"));
    }
    let source = s.model.step(&modal.source_step).ok_or_else(|| EngineError::invalid("step", "missing source frequency step"))?;
    if !matches!(source.procedure, Procedure::Frequency { modes, min_hz:None,max_hz:None } if modes==modal.modes)
        || s.model.steps().iter().position(|x| x.name == source.name) >= s.model.steps().iter().position(|x| x.name == step.name)
        || constraints(&s.model, source)? != constraints(&s.model, &step)?
    {
        return Err(EngineError::invalid(
            "step",
            "random step must follow its unwindowed frequency source with identical homogeneous constraints and mode count",
        ));
    }
    if step.loads.is_empty() {
        return Err(EngineError::invalid("step", "random spatial load pattern is empty"));
    }
    for load in &step.loads {
        if !matches!(
            s.model.loads().get(load),
            Some(fieldcraft_model::Load::Force { .. } | fieldcraft_model::Load::Moment { .. } | fieldcraft_model::Load::Pressure { .. })
        ) {
            return Err(EngineError::invalid("step", "only concentrated force/moment and pressure random patterns are mapped"));
        }
    }
    if s.model
        .steps()
        .iter()
        .take_while(|x| x.name != step.name)
        .any(|x| matches!(x.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_) | Procedure::ExplicitDynamic(_)))
    {
        return Err(EngineError::invalid("step", "random modal execution after a general load history needs a prestress mapping"));
    }
    let limits = super::limits(p)?;
    let threads = super::threads(s, p)?;
    let modes = match s.result(&modal.source_step).map(|r| r.as_ref()) {
        Some(StepResult::Implicit(ImplicitResult::Frequency(f))) => f.clone(),
        _ => super::with_threads(threads, || {
            fieldcraft_solver::solve_frequency(&s.model, Some(&modal.source_step), &FrequencyOptions::default(), &limits)
        })?,
    };
    let natural: Vec<f64> = modes.modes.iter().map(|m| m.frequency).collect();
    let frequencies = modal.frequencies(&natural)?;
    let mut harmonic_step = step.clone();
    harmonic_step.procedure = Procedure::Harmonic(Harmonic {
        frequencies: Frequencies::List(frequencies.clone()),
        method: Method::Modal { modes: modal.modes, damping_ratios: modal.damping_ratios.clone(), residual_vectors: 0 },
        damping: random.damping,
        load_factors: Default::default(),
    });
    let mut model = s.model.clone();
    model.replace_step(harmonic_step)?;
    let h = super::with_threads(threads, || fieldcraft_solver::solve_harmonic(&model, Some(name), &limits))?;
    if h.natural_frequencies.len() != natural.len()
        || h.natural_frequencies.iter().zip(&natural).any(|(a, b)| (a - b).abs() > 1e-8 * a.abs().max(b.abs()))
    {
        return Err(EngineError::invalid("step", "harmonic basis frequencies disagree with the source modal analysis"));
    }
    let table = json!(spectrum.iter().map(|p| [p.frequency, p.psd]).collect::<Vec<_>>());
    let mut out = random_frames::evaluate(&h, p, name, &table)?;
    random_frames::element_output(&model, &h, p, &mut out)?;
    out["method"] = json!("modal_correlated_force_psd");
    out["source_step"] = json!(modal.source_step);
    out["natural_frequencies"] = json!(natural);
    out["damping_ratios"] = json!(modal.damping_ratios);
    out["frequency_points"] = json!(frequencies.len());
    out["modal_full_residual_max"] = json!(h.frames.iter().map(|f| f.relative_residual).fold(0.0_f64, f64::max));
    out["modal_projected_residual_max"] = json!(h.frames.iter().filter_map(|f| f.projected_residual).fold(0.0_f64, f64::max));
    out["warnings"] = json!(h.warnings);
    out["result_storage"] = json!(
        "PSD/RMS returned by this command in SI; save this reply. The project retains the random definition and original input, not random result frames."
    );
    Ok(out)
}
