//! Scaled bounded Levenberg–Marquardt with a data-only identifiability check.
//! Levenberg (1944), doi:10.1090/qam/10666. One-sided Jacobi SVD follows
//! the orthogonal column rotations of Golub & Van Loan, Matrix Computations.
//! Confidence ellipsoids use the local known-error Gaussian information matrix;
//! regularization never supplies missing observational information.
use super::{HarmonicObservation, StaticObservation};
use crate::{Error, Limits, Result};
use fieldcraft_model::{Model, design::Parameter};
use serde::{Deserialize, Serialize};

const RANK_RATIO: f64 = 1e-8; // Dimensionless singular-value ratio.
const GRADIENT_TOLERANCE: f64 = 1e-6; // unit-factors: whitened projected gradient, above FE roundoff
const MIN_DAMPING: f64 = 1e-12; // unit-factors: dimensionless LM damping
const MAX_DAMPING: f64 = 1e12; // unit-factors: dimensionless LM damping

const INITIAL_DAMPING: f64 = 0.001; // unit-factors: dimensionless initial LM damping
const INCREMENT_FLOOR: f64 = 0.001; // unit-factors: dimensionless parameter-scale fraction

#[derive(Clone, Debug, Deserialize)]
pub struct JointParameter {
    pub parameter: Parameter,
    pub lower_si: f64,
    pub upper_si: f64,
    pub scale_si: f64,
}
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum JointResponse {
    Static { step: String, node: u32, dof: u8 },
    MatchedFrequency { group: String, measured_mode: usize },
    ModeShape { group: String, measured_mode: usize, sensor: usize },
    Frequency { step: String, mode: usize },
    HarmonicReal { step: String, node: u32, dof: u8, frame: usize },
    HarmonicImaginary { step: String, node: u32, dof: u8, frame: usize },
}
#[derive(Clone, Debug, Deserialize)]
pub struct JointObservation {
    pub response: JointResponse,
    pub value_si: f64,
    pub standard_deviation_si: f64,
}
#[derive(Clone, Debug, Deserialize)]
pub struct JointFit {
    pub parameters: Vec<JointParameter>,
    pub observations: Vec<JointObservation>,
    pub regularization: f64,
    pub max_iterations: usize,
    pub relative_increment: f64,
}
#[derive(Debug, Serialize)]
pub struct JointEstimate {
    pub parameters: Vec<Parameter>,
    pub modal_pairings: Vec<super::ModalPairingReport>,
    pub estimate_si: Vec<f64>,
    pub predicted_si: Vec<f64>,
    pub residuals_si: Vec<f64>,
    pub jacobian_si: Vec<Vec<f64>>,
    pub parameter_scales_si: Vec<f64>,
    pub singular_values_scaled: Vec<f64>,
    pub rank: usize,
    pub condition_number: f64,
    pub covariance_si: Option<Vec<Vec<f64>>>,
    pub joint_projection_intervals_si: Option<Vec<[f64; 2]>>,
    pub confidence_level: f64,
    pub joint_chi_square: f64,
    pub objective: f64,
    pub evaluations: usize,
    pub iterations: usize,
    pub at_bound: bool,
    pub warnings: Vec<String>,
}

struct Evaluation {
    predicted: Vec<f64>,
    jac: Vec<Vec<f64>>,
    residual: Vec<f64>,
    cost: f64,
    gradient: Vec<f64>,
    info: Vec<Vec<f64>>,
    singular: Vec<f64>,
    covariance: Vec<Vec<f64>>,
    condition: f64,
}

/// At most three columns; rotate the whitened/scaled Jacobian itself, not its
/// normal matrix, so the rank check does not square its condition number.
fn information(j: &[Vec<f64>], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>, f64)> {
    let mut a = j.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..40 {
        let mut changed = false;
        for p in 0..n {
            for q in p + 1..n {
                let aa = a.iter().map(|r| r[p] * r[p]).sum::<f64>();
                let bb = a.iter().map(|r| r[q] * r[q]).sum::<f64>();
                let ab = a.iter().map(|r| r[p] * r[q]).sum::<f64>();
                if ab.abs() <= 8.0 * f64::EPSILON * aa.sqrt() * bb.sqrt() {
                    continue;
                }
                let tau = (bb - aa) / (2.0 * ab);
                let t = if tau >= 0.0 { 1.0 / (tau + libm::hypot(1.0, tau)) } else { -1.0 / (-tau + libm::hypot(1.0, tau)) };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for r in a.iter_mut().chain(v.iter_mut()) {
                    let x = r[p];
                    let y = r[q];
                    r[p] = c * x - s * y;
                    r[q] = s * x + c * y;
                }
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let singular: Vec<f64> = (0..n).map(|k| a.iter().map(|r| r[k] * r[k]).sum::<f64>().sqrt()).collect();
    let largest = singular.iter().copied().fold(0.0_f64, f64::max);
    let rank = singular.iter().filter(|s| s.is_finite() && **s > largest * RANK_RATIO && **s > f64::MIN_POSITIVE).count();
    let smallest = singular.iter().copied().fold(f64::INFINITY, f64::min);
    let condition = largest / smallest;
    if rank < n || !condition.is_finite() || condition > 1.0 / RANK_RATIO {
        let weakest = singular
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(i, _)| i)
            .ok_or_else(|| Error::Invalid("empty joint singular spectrum".into()))?;
        let direction: Vec<_> = v.iter().map(|row| row[weakest]).collect();
        return Err(Error::Invalid(format!(
            "joint inverse unidentifiable data Jacobian: rank {rank}/{n}, condition {condition}, weakest scaled parameter direction {direction:?} (input parameter order); regularization cannot repair missing information; change independent observations or parameters"
        )));
    }
    let mut cov = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..n {
            for (l, s) in singular.iter().enumerate() {
                cov[i][k] += v[i][l] * v[k][l] / (s * s);
            }
        }
    }
    Ok((singular, cov, condition))
}
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Result<Vec<f64>> {
    let n = b.len();
    for i in 0..n {
        let pivot = (i..n).max_by(|&p, &q| a[p][i].abs().total_cmp(&a[q][i].abs())).ok_or_else(|| Error::Invalid("empty joint system".into()))?;
        a.swap(i, pivot);
        b.swap(i, pivot);
        let d = a[i][i];
        if !d.is_finite() || d.abs() <= f64::MIN_POSITIVE {
            return Err(Error::Invalid("singular joint LM system".into()));
        }
        for k in i..n {
            a[i][k] /= d;
        }
        b[i] /= d;
        for r in 0..n {
            if r != i {
                let f = a[r][i];
                for k in i..n {
                    a[r][k] -= f * a[i][k];
                }
                b[r] -= f * b[i];
            }
        }
    }
    Ok(b)
}

fn validate_fit(initial: &[f64], fit: &JointFit) -> Result<()> {
    let n = fit.parameters.len();
    let m = fit.observations.len();
    if !(2..=3).contains(&n)
        || initial.len() != n
        || m < n
        || m > 256
        || !(1..=64).contains(&fit.max_iterations)
        || !fit.regularization.is_finite()
        || fit.regularization < 0.0
        || !fit.relative_increment.is_finite()
        || fit.relative_increment <= 0.0
        || fit.relative_increment > 0.1
    {
        return Err(Error::Invalid(
            "joint fit requires 2..=3 parameters, enough observations, finite regularization and bounded iterations/increment".into(),
        ));
    }
    for (i, p) in fit.parameters.iter().enumerate() {
        if !initial[i].is_finite()
            || !p.lower_si.is_finite()
            || !p.upper_si.is_finite()
            || p.lower_si >= p.upper_si
            || initial[i] < p.lower_si
            || initial[i] > p.upper_si
            || !p.scale_si.is_finite()
            || p.scale_si <= 0.0
            || fit.parameters[..i].iter().any(|q| q.parameter == p.parameter)
        {
            return Err(Error::Invalid("joint parameters need unique identities, valid bounds and positive SI scales".into()));
        }
    }
    if fit.observations.iter().any(|o| !o.value_si.is_finite() || !o.standard_deviation_si.is_finite() || o.standard_deviation_si <= 0.0) {
        return Err(Error::Invalid("joint observations require finite values and known positive SI standard deviations".into()));
    }
    Ok(())
}

/// Shared bounded vector fit. The evaluator returns predictions and an SI
/// Jacobian with one column per parameter. No synthetic truth is accepted.
pub fn fit_joint(initial: &[f64], fit: &JointFit, response: impl Fn(&[f64]) -> Result<(Vec<f64>, Vec<Vec<f64>>)>) -> Result<JointEstimate> {
    validate_fit(initial, fit)?;
    let n = fit.parameters.len();
    let m = fit.observations.len();
    let evaluate = |x: &[f64]| -> Result<Evaluation> {
        fieldcraft_model::progress::check()?;
        let (predicted, jac) = response(x)?;
        if predicted.len() != m
            || jac.len() != m
            || jac.iter().any(|r| r.len() != n)
            || predicted.iter().chain(jac.iter().flatten()).any(|x| !x.is_finite())
        {
            return Err(Error::Invalid("invalid joint predictions/Jacobian dimensions or values".into()));
        }
        let weighted: Vec<Vec<f64>> = jac
            .iter()
            .zip(&fit.observations)
            .map(|(row, o)| row.iter().zip(&fit.parameters).map(|(j, p)| j * p.scale_si / o.standard_deviation_si).collect())
            .collect();
        let (singular, covariance, condition) = information(&weighted, n)?;
        let residual: Vec<f64> = predicted.iter().zip(&fit.observations).map(|(y, o)| y - o.value_si).collect();
        let mut cost = 0.0;
        let mut gradient = vec![0.0; n];
        let mut info = vec![vec![0.0; n]; n];
        for ((row, r), o) in weighted.iter().zip(&residual).zip(&fit.observations) {
            let r = r / o.standard_deviation_si;
            cost += r * r;
            for i in 0..n {
                gradient[i] += row[i] * r;
                for k in 0..n {
                    info[i][k] += row[i] * row[k];
                }
            }
        }
        for i in 0..n {
            let prior = (x[i] - initial[i]) / fit.parameters[i].scale_si;
            cost += fit.regularization * prior * prior;
            gradient[i] += fit.regularization * prior;
            info[i][i] += fit.regularization;
        }
        if !cost.is_finite() || gradient.iter().chain(info.iter().flatten()).any(|x| !x.is_finite()) {
            return Err(Error::Invalid("joint scaled objective overflow".into()));
        }
        Ok(Evaluation { predicted, jac, residual, cost, gradient, info, singular, covariance, condition })
    };
    let mut x = initial.to_vec();
    let mut current = evaluate(&x)?;
    let mut damping = INITIAL_DAMPING;
    let mut floating_point_plateau = None;
    let mut normalized_gradient = f64::INFINITY;

    for iteration in 0..=fit.max_iterations {
        let mut projected = current.gradient.clone();
        for i in 0..n {
            let p = &fit.parameters[i];
            if (x[i] == p.lower_si && projected[i] > 0.0) || (x[i] == p.upper_si && projected[i] < 0.0) {
                projected[i] = 0.0;
            }
        }
        let diagonal_norm = projected.iter().enumerate().map(|(i, g)| g.abs() / current.info[i][i].sqrt()).fold(0.0_f64, f64::max);
        // The diagonal gradient alone hides unconverged weak directions when
        // identifiable parameters are correlated. Use the data-SVD covariance
        // to measure the projected gradient in the full information metric.
        // Keep the previous criterion for regularized fits (different Hessian).
        let norm = if fit.regularization == 0.0 {
            let decrement: f64 =
                projected.iter().enumerate().map(|(i, g)| g * current.covariance[i].iter().zip(&projected).map(|(c, g)| c * g).sum::<f64>()).sum();
            diagonal_norm.max(decrement.max(0.0).sqrt())
        } else {
            diagonal_norm
        };
        normalized_gradient = norm;
        if floating_point_plateau.is_some() || norm < GRADIENT_TOLERANCE {
            let at_bound = x.iter().zip(&fit.parameters).any(|(x, p)| *x == p.lower_si || *x == p.upper_si);
            let chi = if n == 2 { 5.991464547107979 } else { 7.814727903251179 };
            let covariance = (!at_bound && fit.regularization == 0.0).then(|| {
                current
                    .covariance
                    .iter()
                    .enumerate()
                    .map(|(i, r)| r.iter().enumerate().map(|(k, c)| c * fit.parameters[i].scale_si * fit.parameters[k].scale_si).collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            });
            let intervals = covariance.as_ref().map(|c| {
                (0..n)
                    .map(|i| {
                        let half = (chi * c[i][i]).sqrt();
                        [x[i] - half, x[i] + half]
                    })
                    .collect()
            });
            let mut warnings = vec!["Local data-only rank diagnosis; not global identifiability. Covariance is a known-independent-Gaussian-error linearization. Joint projection intervals are projections of one 95% ellipsoid, not marginal intervals. No covariance for regularization or bound-active estimates.".into()];
            if let Some(decrease) = floating_point_plateau {
                warnings.push(format!("Floating-point plateau after rejected step: gradient criterion not met, but the predicted relative decrease is at most sqrt(machine epsilon), below finite-element evaluation noise; normalized_gradient={norm}, predicted_decrease={decrease}, objective={}", current.cost));
            }
            return Ok(JointEstimate {
                modal_pairings: Vec::new(),
                parameters: fit.parameters.iter().map(|p| p.parameter.clone()).collect(),
                estimate_si: x,
                predicted_si: current.predicted,
                residuals_si: current.residual,
                jacobian_si: current.jac,
                parameter_scales_si: fit.parameters.iter().map(|p| p.scale_si).collect(),
                singular_values_scaled: current.singular,
                rank: n,
                condition_number: current.condition,
                covariance_si: covariance,
                joint_projection_intervals_si: intervals,
                confidence_level: 0.95,
                joint_chi_square: chi,
                objective: current.cost,
                evaluations: iteration + 1,
                iterations: iteration,
                at_bound,
                warnings,
            });
        }
        if iteration == fit.max_iterations {
            break;
        }
        let mut matrix = current.info.clone();
        for (i, row) in matrix.iter_mut().enumerate() {
            row[i] *= 1.0 + damping;
        }
        let delta = solve(matrix, current.gradient.iter().map(|g| -g).collect())?;
        let trial: Vec<_> = x.iter().zip(&delta).zip(&fit.parameters).map(|((x, d), p)| (x + p.scale_si * d).clamp(p.lower_si, p.upper_si)).collect();
        let candidate = evaluate(&trial)?;

        if candidate.cost < current.cost {
            x = trial;
            current = candidate;
            damping = (damping / 3.0).max(MIN_DAMPING);
        } else {
            // Use the undamped model: escalating LM damping must never manufacture convergence.
            // With a prior, data covariance is not the inverse of the full information matrix.
            let direction = if fit.regularization == 0.0 {
                current.covariance.iter().map(|row| row.iter().zip(&current.gradient).map(|(c, g)| c * g).sum()).collect()
            } else {
                solve(current.info.clone(), current.gradient.clone())?
            };
            let decrease = 0.5 * current.gradient.iter().zip(direction).map(|(g, d)| g * d).sum::<f64>();
            if decrease.is_finite() && decrease >= 0.0 && decrease <= f64::EPSILON.sqrt() * current.cost {
                floating_point_plateau = Some(decrease);
            }
            damping = (damping * 10.0).min(MAX_DAMPING);
        }
    }
    Err(Error::NotConverged(format!(
        "joint inverse iteration cap; parameters_si={x:?}, normalized_gradient={normalized_gradient}, damping={damping}, objective={}",
        current.cost
    )))
}

/// Assemble Jacobian columns from the existing model derivative implementations.
pub fn inverse_joint(model: &Model, fit: &JointFit, limits: &Limits) -> Result<JointEstimate> {
    inverse_joint_with_modes(model, fit, &[], limits)
}

/// Joint observations with fixed measured-mode groups; pair each trial afresh.
pub fn inverse_joint_with_modes(model: &Model, fit: &JointFit, groups: &[super::ModalObservationGroup], limits: &Limits) -> Result<JointEstimate> {
    super::modal_observation::validate(groups)?;
    for g in groups {
        let step = model.select_step(Some(&g.step), |p| matches!(p, fieldcraft_model::Procedure::Frequency { .. }), "frequency")?;
        if step.base_state.is_some() {
            return Err(Error::NotSupported(format!("sensitivity of prestressed step `{}`", step.name)));
        }
        if !matches!(step.procedure, fieldcraft_model::Procedure::Frequency { modes: 1..=64, min_hz: None, max_hz: None }) {
            return Err(Error::Invalid("modal observation groups require 1..=64 candidate modes and no frequency window".into()));
        }
    }
    let initial = fit.parameters.iter().map(|p| p.parameter.value(model)).collect::<std::result::Result<Vec<_>, _>>()?;
    validate_fit(&initial, fit)?;
    let evaluate = |x: &[f64]| {
        let mut trial = model.clone();
        for (p, x) in fit.parameters.iter().zip(x) {
            p.parameter.set(&mut trial, *x)?;
        }
        let mut predictions = vec![0.0; fit.observations.len()];
        let mut jac = vec![vec![0.0; x.len()]; fit.observations.len()];
        let mut reports = Vec::new();
        let mut static_steps: std::collections::BTreeMap<&str, Vec<(usize, StaticObservation)>> = std::collections::BTreeMap::new();
        for (row, o) in fit.observations.iter().enumerate() {
            if let JointResponse::Static { step, node, dof } = &o.response {
                static_steps.entry(step).or_default().push((row, StaticObservation { node: *node, dof: *dof, stress_element: None }));
            }
        }
        let parameters: Vec<_> = fit
            .parameters
            .iter()
            .zip(x)
            .map(|(p, x)| (p.parameter.clone(), x.abs().max(p.scale_si * INCREMENT_FLOOR) * fit.relative_increment))
            .collect();
        for (step, rows) in static_steps {
            let observations: Vec<_> = rows.iter().map(|(_, o)| *o).collect();
            let result = super::model::static_jacobian(&trial, step, &parameters, &observations, limits)?;
            for ((row, _), (value, derivatives)) in rows.into_iter().zip(result.values.into_iter().zip(result.derivatives)) {
                predictions[row] = value;
                jac[row] = derivatives;
            }
        }

        let mut active = vec![true; fit.observations.len()];
        for (col, p) in fit.parameters.iter().enumerate() {
            let mut modal_cache = std::collections::BTreeMap::new();
            let mut pair_cache = std::collections::BTreeMap::new();
            let mut frequency_cache = std::collections::BTreeMap::new();
            let mut harmonic_cache = std::collections::BTreeMap::new();
            let h = x[col].abs().max(p.scale_si * INCREMENT_FLOOR) * fit.relative_increment;
            for (row, o) in fit.observations.iter().enumerate() {
                let (value, derivative) = match &o.response {
                    JointResponse::MatchedFrequency { group, measured_mode } | JointResponse::ModeShape { group, measured_mode, .. } => {
                        let g = groups.iter().find(|g| &g.id == group).ok_or_else(|| Error::Invalid(format!("unknown modal group {group}")))?;
                        if !modal_cache.contains_key(&g.step) {
                            let result =
                                super::frequency_sensitivity_with_vectors(&trial, Some(&g.step), &p.parameter, h, &Default::default(), true, limits)?;
                            modal_cache.insert(g.step.clone(), result);
                        }
                        let result = modal_cache.get(&g.step).ok_or_else(|| Error::Internal("modal cache missing".into()))?;
                        if !pair_cache.contains_key(group) {
                            pair_cache.insert(group.clone(), super::modal_observation::pairing(g, result)?);
                        }
                        let paired = pair_cache.get(group).ok_or_else(|| Error::Internal("modal pairing missing".into()))?;
                        if paired.removed_modes.iter().any(|r| r.measured_mode == *measured_mode) {
                            active[row] = false;
                            continue;
                        }
                        if let JointResponse::ModeShape { sensor, .. } = &o.response {
                            super::modal_observation::component(g, result, paired, *measured_mode, *sensor)?
                        } else {
                            let index = paired.pair(*measured_mode).ok_or_else(|| Error::Invalid("measured mode out of range".into()))?.computed_mode;
                            let mode = result.modes.get(index).ok_or_else(|| Error::Internal("paired mode missing".into()))?;
                            (mode.frequency, mode.frequency_derivative)
                        }
                    }
                    JointResponse::Static { .. } => continue,
                    JointResponse::Frequency { step, mode } => {
                        if !frequency_cache.contains_key(step) {
                            frequency_cache.insert(
                                step.clone(),
                                super::frequency_sensitivity(&trial, Some(step), &p.parameter, h, &Default::default(), limits)?,
                            );
                        }
                        let r = frequency_cache.get(step).ok_or_else(|| Error::Internal("frequency cache missing".into()))?;
                        let r =
                            r.modes.get(*mode).ok_or_else(|| Error::Invalid("joint frequency mode index is zero-based and out of range".into()))?;
                        if r.modal_mass == 0.0 {
                            return Err(Error::NotSupported(format!("frequency observation points to near-zero mode at order {}", *mode + 1)));
                        }
                        (r.frequency, r.frequency_derivative)
                    }
                    JointResponse::HarmonicReal { step, node, dof, frame } | JointResponse::HarmonicImaginary { step, node, dof, frame } => {
                        let key = (step.clone(), *node, *dof);
                        if !harmonic_cache.contains_key(&key) {
                            harmonic_cache.insert(
                                key.clone(),
                                super::harmonic_sensitivity(
                                    &trial,
                                    Some(step),
                                    &p.parameter,
                                    h,
                                    HarmonicObservation { node: *node, dof: *dof },
                                    limits,
                                )?,
                            );
                        }
                        let r = harmonic_cache.get(&key).ok_or_else(|| Error::Internal("harmonic cache missing".into()))?;
                        let r = r.frames.get(*frame).ok_or_else(|| Error::Invalid("joint harmonic frame out of range".into()))?;
                        let i = usize::from(matches!(o.response, JointResponse::HarmonicImaginary { .. }));
                        (r.displacement[i], r.direct[i])
                    }
                };
                predictions[row] = value;
                jac[row][col] = derivative;
            }
            if col == 0 {
                reports = pair_cache.into_iter().map(|(group, pairing)| super::ModalPairingReport { group, pairing }).collect();
            }
        }
        Ok((predictions, jac, reports, active))
    };
    let active = if groups.iter().any(|g| g.allow_mac_fallback) { evaluate(&initial)?.3 } else { vec![true; fit.observations.len()] };
    let mut retained_fit = fit.clone();
    retained_fit.observations = fit.observations.iter().zip(&active).filter(|(_, keep)| **keep).map(|(o, _)| o.clone()).collect();
    if retained_fit.observations.len() < fit.parameters.len() {
        return Err(Error::Invalid("rank deficient after MAC removal: fewer retained observations than parameters".into()));
    }
    let mut estimate = fit_joint(&initial, &retained_fit, |x| {
        let (y, j, _, current_active) = evaluate(x)?;
        if current_active != active {
            return Err(Error::Invalid("MAC retained observation set changed during fit; review measured modes and restart".into()));
        }
        Ok((
            y.into_iter().zip(&active).filter(|(_, keep)| **keep).map(|(v, _)| v).collect(),
            j.into_iter().zip(&active).filter(|(_, keep)| **keep).map(|(v, _)| v).collect(),
        ))
    })?;
    if !groups.is_empty() {
        let (_, _, reports, _) = evaluate(&estimate.estimate_si)?;
        estimate.modal_pairings = reports;
        estimate.warnings.push("MAC assignments recomputed at every trial and at the returned estimate. Reference-normalized shape errors must have known independent sigma in normalized coordinates; noise in a shared measured reference is correlated and not represented by this covariance.".into());
    }
    Ok(estimate)
}

/// Exact homogeneous E/rho reuse for mixed static and frequency observations.
/// Baselines are actual FE solutions; their uniform scaling derivatives are exact.
/// No springs, point masses, nonzero prescribed motion or mass-dependent loads.
pub fn inverse_joint_uniform_material(model: &Model, fit: &JointFit, limits: &Limits) -> Result<JointEstimate> {
    use fieldcraft_model::{BoundaryCondition, Elastic, Load, design::MaterialConstant as C};
    if fit.parameters.len() != 2 {
        return Err(Error::Invalid("uniform joint reuse requires [Young, density]".into()));
    }
    let Parameter::Material { name, constant: C::Young } = &fit.parameters[0].parameter else {
        return Err(Error::Invalid("first uniform parameter must be Young".into()));
    };
    if fit.parameters[1].parameter != (Parameter::Material { name: name.clone(), constant: C::Density })
        || model.sections().is_empty()
        || model.sections().values().any(|s| s.material != *name)
        || model.sections().values().any(|s| matches!(s.kind, fieldcraft_model::SectionKind::CompositeShell { .. }))
        || !matches!(model.material(name).and_then(|m| m.elastic.as_ref()), Some(Elastic::Isotropic { .. }))
        || !model.springs().is_empty()
        || !model.foundations().is_empty()
        || !model.masses().is_empty()
        || !model.initial_fields().is_empty()
        || model
            .loads()
            .values()
            .any(|l| !matches!(l, Load::Force { .. } | Load::Moment { .. } | Load::Traction { .. } | Load::Pressure { .. } | Load::BodyForce { .. }))
        || model.boundary_conditions().values().any(|b| match b {
            BoundaryCondition::Fixed { .. } | BoundaryCondition::Symmetry { .. } => false,
            BoundaryCondition::Displacement { values, .. } => values.iter().flatten().any(|v| *v != 0.0),
            _ => true,
        })
        || fit.parameters.iter().any(|p| p.lower_si <= 0.0)
    {
        return Err(Error::NotSupported("uniform E/rho joint reuse requires homogeneous isotropic material, homogeneous constraints, no springs/point masses/initial fields or mass-dependent loads".into()));
    }
    let initial = fit.parameters.iter().map(|p| p.parameter.value(model)).collect::<std::result::Result<Vec<_>, _>>()?;
    validate_fit(&initial, fit)?;
    let mut modal = std::collections::BTreeMap::new();
    let mut statics = std::collections::BTreeMap::new();
    let mut baseline = Vec::new();
    for o in &fit.observations {
        let value = match &o.response {
            JointResponse::Static { step, node, dof } => {
                if !statics.contains_key(step) {
                    statics.insert(step.clone(), crate::solve_static(model, Some(step), limits)?);
                }
                let r = statics.get(step).ok_or_else(|| Error::Internal("uniform static baseline".into()))?;
                let i = r.nodes.binary_search(node).map_err(|_| Error::Invalid("uniform static observation node missing".into()))?;
                match dof {
                    1..=3 => r.displacement.get(i).and_then(|v| v.get(usize::from(*dof - 1))),
                    4..=6 => r.rotation.get(i).and_then(|v| v.get(usize::from(*dof - 4))),
                    _ => None,
                }
                .copied()
                .ok_or_else(|| Error::Invalid("uniform static observation dof must be 1..=6".into()))?
            }
            JointResponse::Frequency { step, mode } => {
                if !modal.contains_key(step) {
                    let selected = model.select_step(Some(step), |p| matches!(p, fieldcraft_model::Procedure::Frequency { .. }), "frequency")?;
                    if selected.base_state.is_some() {
                        return Err(Error::NotSupported(format!("sensitivity of prestressed step `{}`", selected.name)));
                    }
                    let fieldcraft_model::Procedure::Frequency { modes, min_hz: None, max_hz: None } = selected.procedure else {
                        return Err(Error::NotSupported("frequency windows in uniform joint reuse".into()));
                    };
                    let mut guarded = model.clone();
                    let mut selected = selected.clone();
                    selected.procedure = fieldcraft_model::Procedure::Frequency {
                        modes: modes.checked_add(1).ok_or_else(|| Error::Invalid("mode count overflow".into()))?,
                        min_hz: None,
                        max_hz: None,
                    };
                    guarded.replace_step(selected)?;
                    let r = crate::solve_frequency(&guarded, Some(step), &Default::default(), limits)?;
                    super::model::require_uncapped_modes(&r, limits)?;
                    let requested = modes as usize;
                    if r.modes.len() < requested || (r.modes.len() == requested && r.stats.equations > requested) {
                        return Err(Error::NotConverged("uniform joint baseline needs requested modes and an upper gap guard".into()));
                    }
                    for pair in r.modes.windows(2) {
                        if let [a, b] = pair
                            && (a.eigenvalue - b.eigenvalue).abs() <= 1e-5 * a.eigenvalue.abs().max(b.eigenvalue.abs())
                        {
                            return Err(Error::NotSupported("clustered modes in uniform joint baseline".into()));
                        }
                    }
                    modal.insert(step.clone(), r.modes.iter().take(requested).map(|m| m.frequency).collect::<Vec<_>>());
                }
                *modal.get(step).and_then(|m| m.get(*mode)).ok_or_else(|| Error::Invalid("uniform joint mode out of range".into()))?
            }
            _ => return Err(Error::NotSupported("uniform E/rho reuse supports static and simple frequency observations only".into())),
        };
        baseline.push(value);
    }
    let mut r = fit_joint(&initial, fit, |x| {
        let mut y = Vec::new();
        let mut j = Vec::new();
        for (o, b) in fit.observations.iter().zip(&baseline) {
            let (v, row) = if matches!(o.response, JointResponse::Static { .. }) {
                let v = b * initial[0] / x[0];
                (v, vec![-v / x[0], 0.0])
            } else {
                let v = b * ((x[0] / initial[0]) * (initial[1] / x[1])).sqrt();
                (v, vec![v / (2.0 * x[0]), -v / (2.0 * x[1])])
            };
            y.push(v);
            j.push(row);
        }
        Ok((y, j))
    })?;
    r.warnings.push(format!("Exact homogeneous E/rho identities after {} actual FE static baseline solves and {} actual modal baseline solves. Optimizer evaluations are not additional FE solves.",statics.len(),modal.len()));
    Ok(r)
}

#[cfg(test)]
mod platform_tests {
    use super::*;
    use fieldcraft_model::design::MaterialConstant;
    fn run(initial: &[f64], observations: &[f64], bias: f64, cap: usize) -> Result<JointEstimate> {
        let parameters = [MaterialConstant::Young, MaterialConstant::Density]
            .into_iter()
            .map(|constant| JointParameter {
                parameter: Parameter::Material { name: "m".into(), constant },
                lower_si: 0.5,
                upper_si: 2.5,
                scale_si: 1.0,
            })
            .collect();
        let fit = JointFit {
            parameters,
            observations: observations
                .iter()
                .enumerate()
                .map(|(mode, &value_si)| JointObservation {
                    response: JointResponse::Frequency { step: "F".into(), mode },
                    value_si,
                    standard_deviation_si: 1.0,
                })
                .collect(),
            regularization: 0.0,
            max_iterations: cap,
            relative_increment: 1e-5,
        };
        fit_joint(initial, &fit, |x| {
            Ok((vec![x[0], x[0], x[1], x[1]], vec![vec![1.0 + bias, 0.0], vec![1.0, 0.0], vec![0.0, 1.0 + bias], vec![0.0, 1.0]]))
        })
    }
    #[test]
    fn joint_zero_residual_snapshot() {
        let r = run(&[0.9, 1.8], &[1.0, 1.0, 2.0, 2.0], 0.0, 32).unwrap();
        assert_eq!(
            (r.estimate_si.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), r.objective.to_bits(), r.iterations, r.evaluations),
            (vec![4607182418500177320, 4611686018127547816], 4398033392809034205, 2, 3)
        );
        eprintln!(
            "JOINT_ZERO {:?}",
            (r.estimate_si.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), r.objective.to_bits(), r.iterations, r.evaluations)
        );
        assert!(!r.warnings.iter().any(|w| w.contains("Floating-point plateau")));
    }
    #[test]
    fn joint_derivative_bias_platform() {
        let r = run(&[1.0, 2.0], &[-49.0, 51.0, -48.0, 52.0], 5e-8, 1).unwrap();
        for (i, &center) in [1.0, 2.0].iter().enumerate() {
            let scanned = (-100..=100)
                .map(|k| center + f64::from(k) * 1e-4)
                .min_by(|a, b| {
                    let cost = |x: f64| (x - center + 50.0).powi(2) + (x - center - 50.0).powi(2);
                    cost(*a).total_cmp(&cost(*b))
                })
                .unwrap();
            assert_eq!(r.estimate_si[i].to_bits(), scanned.to_bits());
        }
        assert_eq!((r.iterations, r.evaluations), (1, 2));
        assert!(r.covariance_si.is_some());
        assert!(r.warnings.iter().any(|w| w.contains("Floating-point plateau") && w.contains("normalized_gradient=")));
        eprintln!("JOINT_PLATFORM estimate={:?} objective={} {:?}", r.estimate_si, r.objective, r.warnings);
        let error = run(&[0.9, 1.8], &[1.0, 1.0, 2.0, 2.0], 0.0, 1).unwrap_err().to_string();
        for key in ["parameters_si=", "normalized_gradient=", "damping=", "objective="] {
            assert!(error.contains(key), "{error}");
        }
        assert!(run(&[1.0, 2.0], &[-49.0, 51.0, -48.0, 52.0], 0.1, 32).is_err());
    }
}
