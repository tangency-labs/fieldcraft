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
    let initial = parameter.value(model)?;
    if !fit.lower.is_finite() || !fit.upper.is_finite() || fit.lower >= fit.upper || initial < fit.lower || initial > fit.upper {
        return Err(Error::Invalid("frequency fit needs finite ordered bounds containing the initial parameter".into()));
    }
    if fit.observed.is_empty() || fit.observed.len() > 256 || fit.observed.iter().any(|y| !y.is_finite() || *y <= 0.0) {
        return Err(Error::Invalid("frequency fit needs 1..=256 positive finite observations".into()));
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
    let sigma = fit.standard_deviation.as_ref().unwrap_or(&fit.observed);
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
        let response = frequency_sensitivity(&trial, step, parameter, increment, &FrequencyOptions::default(), limits)?;
        if response.modes.len() != fit.observed.len() {
            return Err(Error::Invalid("observation count must equal the frequency step mode count".into()));
        }
        let mut cost = 0.0;
        let mut gradient = 0.0;
        let mut data_information = 0.0;
        let mut predicted = Vec::new();
        let mut residuals = Vec::new();
        for ((mode, obs), s) in response.modes.iter().zip(&fit.observed).zip(sigma) {
            let r = (mode.frequency - obs) / s;
            let j = mode.frequency_derivative * scale / s;
            cost += r * r;
            gradient += j * r;
            data_information += j * j;
            predicted.push(mode.frequency);
            residuals.push(mode.frequency - obs);
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
            let mut warnings = vec!["One-parameter local fit; mode identity must remain simple and ordered throughout the supplied bounds.".into()];
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
