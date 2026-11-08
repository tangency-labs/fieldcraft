//! Bounded one-parameter Levenberg–Marquardt fits of measured frequencies.
//! Levenberg (1944), Quarterly of Applied Mathematics 2:164–168,
//! doi:10.1090/qam/10666. The scalar normal equation is damped in scaled
//! parameter coordinates; each trial is projected onto the supplied bounds.
//! This entry point intentionally reports its one-parameter scope. It does
//! not resolve multi-parameter identifiability or globally locate all minima.
use super::frequency_sensitivity;
use crate::{Error, FrequencyOptions, Limits, Result};
use fieldcraft_model::{Model, design::Parameter};
use serde::Serialize;

const MIN_INCREMENT_SCALE: f64 = 1e-3; // unit-factors: fraction of the normalized parameter scale, not a unit conversion
const INITIAL_DAMPING: f64 = 1e-3; // unit-factors: dimensionless relative LM diagonal damping
const STATIONARITY_TOLERANCE: f64 = 1e-9; // unit-factors: dimensionless normalized gradient tolerance
const MIN_DAMPING: f64 = 1e-12; // unit-factors: dimensionless lower LM damping bound
const MAX_DAMPING: f64 = 1e12; // unit-factors: dimensionless upper LM damping bound

#[derive(Debug, Clone)]
pub struct FrequencyFit {
    pub lower: f64,
    pub upper: f64,
    pub observed: Vec<f64>,
    pub mode_indices: Option<Vec<usize>>,
    /// Known independent Gaussian observation standard deviations, in Hz.
    /// None uses relative residuals and does not claim a confidence interval.
    pub standard_deviation: Option<Vec<f64>>,
    /// λ ((p-p_initial)/parameter_scale)^2.
    pub regularization: f64,
    pub max_iterations: usize,
    pub relative_increment: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct FrequencyEstimate {
    pub parameter: Parameter,
    pub estimate: f64,
    pub predicted: Vec<f64>,
    pub residuals: Vec<f64>,
    pub objective: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub model_evaluations: usize,
    pub evaluation_method: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub baseline: Option<super::FrequencySensitivity>,
    pub at_bound: bool,
    /// Local linearized 95% normal interval, only for known standard
    /// deviations, no regularization, an interior estimate and nonzero slope.
    pub confidence_interval: Option<[f64; 2]>,
    pub warnings: Vec<String>,
}

struct Evaluation {
    predicted: Vec<f64>,
    residuals: Vec<f64>,
    cost: f64,
    gradient: f64,
    information: f64,
    data_information: f64,
}

/// Fit one existing material/section/support property to the frequency step's
/// ordered simple modes. The input model is never changed. Each evaluation
/// recomputes actual FE modes and derivatives; there is no global E-scaling
/// surrogate. Clustered modes are refused by the derivative evaluator.
pub fn inverse_frequency(model: &Model, step: Option<&str>, parameter: &Parameter, fit: &FrequencyFit, limits: &Limits) -> Result<FrequencyEstimate> {
    if fit.observed.iter().any(|y| *y <= 0.0) {
        return Err(Error::Invalid("frequency observations must be positive".into()));
    }
    fit_scalar(model, parameter, fit, &fit.observed, "Simple ordered modes must retain their identity throughout the bounds.", |trial, increment| {
        let response = frequency_sensitivity(trial, step, parameter, increment, &FrequencyOptions::default(), limits)?;
        let indices = fit.mode_indices.as_deref().map_or_else(|| (0..fit.observed.len()).collect(), ToOwned::to_owned);
        indices
            .iter()
            .map(|&i| {
                response.modes.get(i).ok_or_else(|| Error::Invalid(format!("frequency observation mode index {i} is out of range"))).and_then(|m| {
                    if m.modal_mass == 0.0 {
                        Err(Error::NotSupported(format!("frequency observation points to near-zero mode index {i} (order {})", i + 1)))
                    } else {
                        Ok((m.frequency, m.frequency_derivative))
                    }
                })
            })
            .collect()
    })
}

pub(super) fn fit_scalar(
    model: &Model,
    parameter: &Parameter,
    fit: &FrequencyFit,
    relative_scales: &[f64],
    condition: &str,
    response: impl Fn(&Model, f64) -> Result<Vec<(f64, f64)>>,
) -> Result<FrequencyEstimate> {
    if let Some(indices) = &fit.mode_indices
        && (indices.len() != fit.observed.len() || indices.iter().enumerate().any(|(i, v)| indices[..i].contains(v)))
    {
        return Err(Error::Invalid("mode_indices must be unique zero-based indices, one per observation".into()));
    }
    let initial = parameter.value(model)?;
    if !fit.lower.is_finite() || !fit.upper.is_finite() || fit.lower >= fit.upper || initial < fit.lower || initial > fit.upper {
        return Err(Error::Invalid("frequency fit needs finite ordered bounds containing the initial parameter".into()));
    }
    if fit.observed.is_empty() || fit.observed.len() > 256 || fit.observed.iter().any(|y| !y.is_finite()) {
        return Err(Error::Invalid("frequency fit needs 1..=256 finite scalar observations".into()));
    }
    if !fit.regularization.is_finite()
        || fit.regularization < 0.0
        || !(1..=64).contains(&fit.max_iterations)
        || !fit.relative_increment.is_finite()
        || fit.relative_increment <= 0.0
        || fit.relative_increment > 0.1
    {
        return Err(Error::Invalid("invalid frequency fit regularization, iteration cap or relative increment".into()));
    }
    let sigma = fit.standard_deviation.as_deref().unwrap_or(relative_scales);
    if sigma.len() != fit.observed.len() || sigma.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(Error::Invalid("standard deviations must be positive finite and match observations".into()));
    }
    let scale = initial.abs().max((fit.upper - fit.lower) * 0.1);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(Error::Invalid("frequency fit parameter scale overflow".into()));
    }
    let evaluate = |value: f64| -> Result<Evaluation> {
        fieldcraft_model::progress::check()?;
        let mut trial = model.clone();
        parameter.set(&mut trial, value)?;
        let increment = value.abs().max(scale * MIN_INCREMENT_SCALE) * fit.relative_increment;
        let values = response(&trial, increment)?;
        if values.len() != fit.observed.len() {
            return Err(Error::Invalid("observation count must equal prediction count".into()));
        }
        let mut cost = 0.0;
        let mut gradient = 0.0;
        let mut data_information = 0.0;
        let mut predicted = Vec::new();
        let mut residuals = Vec::new();
        for (((prediction, derivative), obs), s) in values.iter().zip(&fit.observed).zip(sigma) {
            let r = (prediction - obs) / s;
            let j = derivative * scale / s;
            cost += r * r;
            gradient += j * r;
            data_information += j * j;
            predicted.push(*prediction);
            residuals.push(prediction - obs);
        }
        let prior = (value - initial) / scale;
        cost += fit.regularization * prior * prior;
        gradient += fit.regularization * prior;
        if !cost.is_finite()
            || !gradient.is_finite()
            || !data_information.is_finite()
            || !(data_information + fit.regularization).is_finite()
            || data_information <= f64::MIN_POSITIVE
        {
            return Err(Error::Invalid("frequency fit has non-finite values or no local parameter information".into()));
        }
        Ok(Evaluation { predicted, residuals, cost, gradient, information: data_information + fit.regularization, data_information })
    };
    let mut value = initial;
    let mut current = evaluate(value)?;
    let mut damping = INITIAL_DAMPING;
    for iteration in 0..=fit.max_iterations {
        let at_lower = value == fit.lower;
        let at_upper = value == fit.upper;
        let stationary = current.gradient.abs() / current.information.sqrt() <= STATIONARITY_TOLERANCE;
        if stationary || (at_lower && current.gradient >= 0.0) || (at_upper && current.gradient <= 0.0) {
            let at_bound = at_lower || at_upper;
            let confidence_interval = if fit.standard_deviation.is_some() && !at_bound && fit.regularization == 0.0 {
                let half = 1.959963984540054 * scale / current.data_information.sqrt();
                let interval = [value - half, value + half];
                interval.iter().all(|v| v.is_finite()).then_some(interval)
            } else {
                None
            };
            let mut warnings = vec![format!("One-parameter local fit. {condition}")];
            if confidence_interval.is_none() {
                warnings.push(
                    "No confidence interval: known observation standard deviations, an interior estimate and zero regularization are required."
                        .into(),
                );
            } else {
                warnings.push("Confidence interval is a local Gaussian linearization conditional on the model and known independent observation errors; it excludes model discrepancy.".into());
            }
            return Ok(FrequencyEstimate {
                parameter: parameter.clone(),
                estimate: value,
                predicted: current.predicted,
                residuals: current.residuals,
                objective: current.cost,
                iterations: iteration,
                evaluations: iteration + 1,
                model_evaluations: iteration + 1,
                evaluation_method: "finite_element_each_evaluation".into(),
                baseline: None,
                at_bound,
                confidence_interval,
                warnings,
            });
        }
        if iteration == fit.max_iterations {
            break;
        }
        let delta = -current.gradient / (current.information * (1.0 + damping));
        let trial = (value + scale * delta).clamp(fit.lower, fit.upper);
        if !trial.is_finite() {
            return Err(Error::Invalid("frequency fit step overflow".into()));
        }
        let candidate = evaluate(trial)?;
        if candidate.cost < current.cost {
            value = trial;
            current = candidate;
            damping = (damping / 3.0).max(MIN_DAMPING);
        } else {
            damping = (damping * 10.0).min(MAX_DAMPING);
        }
    }
    Err(Error::NotConverged(format!(
        "frequency inverse reached {} iterations, objective {} at parameter {}",
        fit.max_iterations, current.cost, value
    )))
}

/// Opt-in exact reuse for K(E)=E/E₀ K(E₀), with constant mass. This is
/// deliberately restricted to one isotropic material across every section,
/// with no springs or initial fields. Mixed-material/support fits must use
/// the general evaluator. One actual FE derivative solve supplies the basis.
pub fn inverse_frequency_uniform_young(
    model: &Model,
    step: Option<&str>,
    parameter: &Parameter,
    fit: &FrequencyFit,
    limits: &Limits,
) -> Result<FrequencyEstimate> {
    use fieldcraft_model::{Elastic, design::MaterialConstant};
    let Parameter::Material { name, constant: MaterialConstant::Young } = parameter else {
        return Err(Error::NotSupported("uniform Young scaling requires an isotropic Young parameter".into()));
    };
    if model.sections().is_empty()
        || model.sections().values().any(|s| s.material != *name)
        || model.sections().values().any(|s| matches!(s.kind, fieldcraft_model::SectionKind::CompositeShell { .. }))
        || !model.springs().is_empty()
        || !model.foundations().is_empty()
        || !model.initial_fields().is_empty()
        || !matches!(model.material(name).and_then(|m| m.elastic.as_ref()), Some(Elastic::Isotropic { .. }))
        || fit.lower <= 0.0
        || fit.observed.iter().any(|y| *y <= 0.0)
    {
        return Err(Error::NotSupported(
            "uniform Young scaling requires all sections in one isotropic material, no springs/initial fields and positive bounds/observations"
                .into(),
        ));
    }
    let initial = parameter.value(model)?;
    let baseline = frequency_sensitivity(model, step, parameter, initial * fit.relative_increment, &FrequencyOptions::default(), limits)?;
    let mut result = fit_scalar(
        model,
        parameter,
        fit,
        &fit.observed,
        "Exact uniform-isotropic Young scaling after one actual FE sensitivity solve; no repeated eigen solves. The model topology, Poisson ratio and mass are fixed.",
        |trial, _increment| {
            let value = parameter.value(trial)?;
            let ratio = (value / initial).sqrt();
            let indices = fit.mode_indices.as_deref().map_or_else(|| (0..fit.observed.len()).collect(), ToOwned::to_owned);
            indices
                .iter()
                .map(|&i| {
                    baseline.modes.get(i).ok_or_else(|| Error::Invalid(format!("frequency observation mode index {i} is out of range"))).and_then(
                        |m| {
                            if m.modal_mass == 0.0 {
                                Err(Error::NotSupported(format!("frequency observation points to near-zero mode index {i} (order {})", i + 1)))
                            } else {
                                let f = m.frequency * ratio;
                                Ok((f, f / (2.0 * value)))
                            }
                        },
                    )
                })
                .collect()
        },
    )?;
    result.model_evaluations = 1;
    result.evaluation_method = "exact_uniform_young_scaling".into();
    result.baseline = Some(baseline);
    Ok(result)
}
