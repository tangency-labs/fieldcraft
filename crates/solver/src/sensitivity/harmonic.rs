//! Physical harmonic parameter derivatives at fixed excitation frequencies.
//! Central differences apply to assembled operators, including material damping
//! and prescribed-motion elimination, not independently solved responses.
use fieldcraft_fem::DofKind;
use fieldcraft_model::{
    Model, Procedure,
    design::{DampingParameter, Parameter},
    harmonic::Method,
};
use serde::Serialize;

use super::model::{cold, difference, matrix_difference, perturb};
use crate::analysis::{self, Setup};
use crate::complex::{self, Vector};
use crate::harmonic::{self, Matrices};
use crate::system::Mesh;
use crate::{Error, Limits, MassKind, Result};

const CANCELLATION_TOLERANCE: f64 = 1e-12; // unit-factors: dimensionless dynamic-operator cancellation threshold, not a unit conversion

#[derive(Debug, Clone, Copy)]
pub struct HarmonicObservation {
    pub node: u32,
    /// Structural DOF 1..=6, translations followed by rotations.
    pub dof: u8,
}

#[derive(Debug, Clone, Serialize)]
pub struct HarmonicDerivative {
    pub frequency: f64,
    pub displacement: [f64; 2],
    pub direct: [f64; 2],
    pub adjoint: [f64; 2],
    pub amplitude: f64,
    /// Principal phase in radians; its derivative is locally unwrapped.
    pub phase: f64,
    pub amplitude_derivative: f64,
    pub phase_derivative: f64,
    pub primal_relative_residual: f64,
    pub tangent_relative_residual: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct HarmonicSensitivity {
    pub parameter: Parameter,
    pub parameter_value: f64,
    pub increment: f64,
    pub frames: Vec<HarmonicDerivative>,
    pub warnings: Vec<String>,
}

/// One column/sweep of a bounded, invocation-local harmonic batch.
pub struct HarmonicSensitivityRequest {
    pub step: Option<String>,
    pub parameter: Parameter,
    pub increment: f64,
    pub observation: HarmonicObservation,
}
#[derive(Debug, Serialize)]
pub struct HarmonicSensitivityBatch {
    pub results: Vec<HarmonicSensitivity>,
    pub modal_basis_builds: usize,
}
struct CachedBasis {
    step: fieldcraft_model::Step,
    vectors: Vec<Vec<f64>>,
    eigenvalues: Vec<f64>,
}
/// Reuse modal/Ritz bases only within this immutable-model call. All assembly,
/// basis extraction and projections occur here; no cross-call cache is accepted.
pub fn harmonic_sensitivity_batch(model: &Model, requests: &[HarmonicSensitivityRequest], limits: &Limits) -> Result<HarmonicSensitivityBatch> {
    if requests.is_empty() || requests.len() > 16 {
        return Err(Error::Invalid("harmonic batch requires 1..=16 requests".into()));
    }
    let mut cache = Vec::new();
    let mut results = Vec::new();
    for r in requests {
        results.push(harmonic_sensitivity_cached(model, r.step.as_deref(), &r.parameter, r.increment, r.observation, limits, &mut cache)?);
    }
    Ok(HarmonicSensitivityBatch { results, modal_basis_builds: cache.len() })
}
/// Differentiate physical equilibrium in a direct or fixed modal/Ritz basis.
/// A zero observation refuses amplitude/phase derivatives rather than inventing them.
pub fn harmonic_sensitivity(
    model: &Model,
    step: Option<&str>,
    parameter: &Parameter,
    increment: f64,
    observation: HarmonicObservation,
    limits: &Limits,
) -> Result<HarmonicSensitivity> {
    harmonic_sensitivity_cached(model, step, parameter, increment, observation, limits, &mut Vec::new())
}
fn harmonic_sensitivity_cached(
    model: &Model,
    step: Option<&str>,
    parameter: &Parameter,
    increment: f64,
    observation: HarmonicObservation,
    limits: &Limits,
    cache: &mut Vec<CachedBasis>,
) -> Result<HarmonicSensitivity> {
    let damping_parameter = if let Parameter::Damping { step, component } = parameter { Some((step, *component)) } else { None };
    let perturbed = if damping_parameter.is_none() { Some(perturb(model, parameter, increment)?) } else { None };
    let resolved = model.resolve()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::Harmonic(_)), "harmonic")?;
    let Procedure::Harmonic(h) = &step.procedure else {
        return Err(Error::Internal("harmonic sensitivity step".into()));
    };
    if let Some((parameter_step, component)) = damping_parameter {
        if parameter_step != &step.name {
            return Err(Error::Invalid("damping parameter belongs to a different harmonic step".into()));
        }
        if matches!(component, DampingParameter::ModalRatio { .. }) && h.method == Method::Direct {
            return Err(Error::NotSupported("modal-ratio derivatives require modal response".into()));
        }
    }
    let mesh = Mesh::new(model, &resolved, limits.max_dofs)?;
    let setup = Setup { model, resolved: &resolved, mesh: &mesh, step };
    cold(&setup)?;
    let kind = observation
        .dof
        .checked_sub(1)
        .and_then(|i| DofKind::from_direction(i as usize))
        .ok_or_else(|| Error::Invalid("sensitivity dof must be 1..=6".into()))?;
    let node = mesh.node(observation.node).ok_or_else(|| Error::Invalid("sensitivity node does not exist".into()))?;
    let dof = mesh.dof(node, kind)?;
    let (graph, order) = analysis::pattern(&mesh);
    let mut warnings: Vec<String> = resolved.warnings.iter().map(ToString::to_string).chain(mesh.warnings.iter().cloned()).collect();
    warnings.extend(crate::system::inelastic_warnings(model));
    let red = analysis::reduction(&setup, true, &order, analysis::budget(limits, &mesh, &graph, 12), &mut warnings)?;
    if red.neq() == 0 {
        return Err(Error::Invalid("every harmonic degree of freedom is prescribed".into()));
    }
    let assemble = |m: &Model| -> Result<_> {
        let r = m.resolve()?;
        let other = Mesh::new(m, &r, limits.max_dofs)?;
        if other.labels != mesh.labels || other.dofs() != mesh.dofs() {
            return Err(Error::Invalid("parameter changed harmonic mesh layout".into()));
        }
        let a = analysis::assemble(&other, &graph, Some(MassKind::Consistent))?;
        let mut k = a.k;
        let mut mass = a.m.ok_or_else(|| Error::Internal("harmonic sensitivity mass".into()))?;
        k.mirror_lower();
        mass.mirror_lower();
        let (c, g) = harmonic::damping_operators(m, &r, &other, &graph)?;
        let f = harmonic::load_vector(m, &r, &other, step, &red.active)?;
        Ok((k, mass, c, g, f))
    };
    let original = assemble(model)?;
    let perturbed_operators = perturbed.as_ref().map(|(p, m)| -> Result<_> { Ok((assemble(p)?, assemble(m)?)) }).transpose()?;
    let reduced_modal = if matches!(h.method, Method::Modal { .. }) {
        if red.largest_prescribed() != 0.0 {
            return Err(Error::NotSupported("fixed-basis modal derivatives require homogeneous support motion".into()));
        }
        let [mut k, mut m] = red.reduce([&original.0, &original.1])?;
        k.mirror_lower();
        m.mirror_lower();
        let reduce_optional = |a: Option<&crate::sparse::Csr>| -> Result<_> {
            a.map(|a| -> Result<_> {
                let [mut r] = red.reduce([a])?;
                r.mirror_lower();
                Ok(r)
            })
            .transpose()
        };
        Some((k, m, reduce_optional(original.2.as_ref())?, reduce_optional(original.3.as_ref())?))
    } else {
        None
    };
    let matrices = reduced_modal.as_ref().map(|(k, m, c, g)| Matrices { stiffness: k, mass: m, viscous: c.as_ref(), loss: g.as_ref() });
    let modal = if let (Some(matrices), Method::Modal { modes, residual_vectors, .. }) = (&matrices, &h.method) {
        // Frequency and step damping do not enter the real K/M eigenproblem or
        // Ritz load. Preserve every other step field, including load factors,
        // constraints, mode count and residual cycles, in the exact cache key.
        let mut key = step.clone();
        key.name.clear();
        if let Procedure::Harmonic(h) = &mut key.procedure {
            h.frequencies = fieldcraft_model::harmonic::Frequencies::List(vec![0.0]);
            h.damping = Default::default();
            if let Method::Modal { damping_ratios, .. } = &mut h.method {
                damping_ratios.clear();
            }
        }
        let index = if let Some(index) = cache.iter().position(|entry| entry.step == key) {
            index
        } else {
            let force = Vector { real: red.restrict(&original.4.real), imaginary: red.restrict(&original.4.imaginary) };
            let modal = harmonic::Modal::with_residual_vectors(matrices, *modes as usize, &force, *residual_vectors as usize, limits)?;
            cache.push(CachedBasis { step: key, vectors: modal.basis_vectors().to_vec(), eigenvalues: modal.eigenvalues().to_vec() });
            cache.len() - 1
        };
        Some(cache.get(index).ok_or_else(|| Error::Internal("modal batch cache index".into()))?)
    } else {
        None
    };
    let ratios = if let Method::Modal { damping_ratios, .. } = &h.method { damping_ratios.as_slice() } else { &[] };
    if modal.is_some() && damping_parameter.is_none() && ratios.iter().any(|x| *x != 0.0) {
        return Err(Error::NotSupported(
            "physical-parameter fixed-basis derivatives with nonzero modal ratios need differentiated modal damping operators".into(),
        ));
    }
    if modal.is_some() {
        warnings.push("Fixed modal/Ritz basis derivative: differentiates projected physical operators, not the truncation basis. Full-basis limit agrees with direct sensitivity; validate truncation separately. Modal-ratio parameters hold the basis, mass and eigenfrequencies fixed.".into());
    }

    let mut weights = vec![0.0; mesh.dofs()];
    *weights.get_mut(dof).ok_or_else(|| Error::Internal("harmonic observation dof".into()))? = 1.0;
    let weights = Vector { real: red.restrict(&weights), imaginary: vec![0.0; red.neq()] };
    if let (Some(modal), Some(matrices)) = (&modal, matrices) {
        let derivative = if let Some((p, m)) = &perturbed_operators {
            let reduce = |a: &crate::sparse::Csr| -> Result<_> {
                let [mut r] = red.reduce([a])?;
                r.mirror_lower();
                Ok(r)
            };
            let optional = |p: Option<&crate::sparse::Csr>, m: Option<&crate::sparse::Csr>| -> Result<_> {
                match (p, m) {
                    (None, None) => Ok(None),
                    (Some(p), Some(m)) => Ok(Some(reduce(&matrix_difference(p, m, increment)?)?)),
                    _ => Err(Error::NotSupported("parameter changes modal damping sparsity".into())),
                }
            };
            Some(super::projected::Operators {
                k: reduce(&matrix_difference(&p.0, &m.0, increment)?)?,
                m: reduce(&matrix_difference(&p.1, &m.1, increment)?)?,
                c: optional(p.2.as_ref(), m.2.as_ref())?,
                g: optional(p.3.as_ref(), m.3.as_ref())?,
                f: Vector {
                    real: red.restrict(&difference(&p.4.real, &m.4.real, increment)?),
                    imaginary: red.restrict(&difference(&p.4.imaginary, &m.4.imaginary, increment)?),
                },
            })
        } else {
            None
        };
        let force = Vector { real: red.restrict(&original.4.real), imaginary: red.restrict(&original.4.imaginary) };
        warnings.push("Residual fields are projected-equation residuals, not full-space truncation-error bounds. All sparse operator projections occur once inside this timed sweep.".into());
        return super::projected::solve(
            super::projected::Sweep {
                basis: &modal.vectors,
                eigenvalues: &modal.eigenvalues,
                matrices,
                derivative,
                force,
                weights: &weights,
                harmonic: h,
                parameter,
                value: parameter.value(model)?,
                increment,
                warnings,
            },
            limits,
        );
    }
    let mut frames = Vec::new();
    for frequency in h.frequencies.points()? {
        fieldcraft_model::progress::check()?;
        let omega = 2.0 * std::f64::consts::PI * frequency;
        let dynamic = |(k, mass, c, g, f): &(_, _, Option<_>, Option<_>, Vector)| -> Result<_> {
            let (a, b) = Matrices { stiffness: k, mass, viscous: c.as_ref(), loss: g.as_ref() }.dynamic(omega, h.damping, limits)?;
            let rhs = Vector { real: red.rhs(&f.real, &a), imaginary: red.rhs(&f.imaginary, &b) };
            let [mut a, mut b] = red.reduce([&a, &b])?;
            a.mirror_lower();
            b.mirror_lower();
            Ok((a, b, rhs))
        };
        let (a, b, rhs) = dynamic(&original)?;
        // Catch complete cancellation even when round-off leaves a nonzero scalar
        // pivot. A factor pivot test scaled only by Z cannot detect this case.
        let [kr, mr] = red.reduce([&original.0, &original.1])?;
        let reference = kr.val.iter().zip(&mr.val).map(|(k, m)| k.abs() + omega * omega * m.abs()).fold(0.0_f64, f64::max);
        let dynamic_scale = a.val.iter().zip(&b.val).map(|(a, b)| a.abs() + b.abs()).fold(0.0_f64, f64::max);
        if dynamic_scale <= CANCELLATION_TOLERANCE * reference {
            return Err(Error::Singular(format!("harmonic sensitivity unresolved dynamic cancellation at {frequency} Hz")));
        }

        let (da, db, df) = if let Some((_, component)) = damping_parameter {
            let (operator, factor) = match component {
                DampingParameter::Alpha => (&original.1, omega),
                DampingParameter::Beta => (&original.0, omega),
                DampingParameter::LossFactor => (&original.0, 1.0),
                DampingParameter::ModalRatio { .. } => (&original.0, 0.0),
            };
            let mut full_db = operator.clone();
            for value in &mut full_db.val {
                *value *= factor;
            }
            let zero = vec![0.0; mesh.dofs()];
            let df = Vector { real: vec![0.0; red.neq()], imaginary: red.rhs(&zero, &full_db) };
            let [mut db] = red.reduce([&full_db])?;
            db.mirror_lower();
            let da = crate::sparse::Csr { val: vec![0.0; db.val.len()], ..db.clone() };
            (da, db, df)
        } else {
            let (p, m) = perturbed_operators.as_ref().ok_or_else(|| Error::Internal("missing perturbed harmonic operators".into()))?;
            let (ap, bp, fp) = dynamic(p)?;
            let (am, bm, fm) = dynamic(m)?;
            (
                matrix_difference(&ap, &am, increment)?,
                matrix_difference(&bp, &bm, increment)?,
                Vector { real: difference(&fp.real, &fm.real, increment)?, imaginary: difference(&fp.imaginary, &fm.imaginary, increment)? },
            )
        };
        let primal = complex::solve(&a, &b, &rhs, limits)?;
        let tangent = super::tangent_load(&da, &db, &primal.displacement, &df, limits)?;
        let derivative = super::direct(&a, &b, &tangent, limits)?;
        let adjoint = super::adjoint(&a, &b, &weights, &tangent, limits)?;
        let extract = |u: &Vector, prescribed: bool| -> Result<[f64; 2]> {
            let real = red.expand(&u.real, prescribed);
            let imaginary = red.expand(&u.imaginary, false);
            Ok([
                *real.get(dof).ok_or_else(|| Error::Internal("harmonic observation".into()))?,
                *imaginary.get(dof).ok_or_else(|| Error::Internal("harmonic observation".into()))?,
            ])
        };
        let value = extract(&primal.displacement, true)?;
        let direct = extract(&derivative.displacement, false)?;
        let [amplitude_derivative, phase_derivative] = super::amplitude_phase(value, direct)?;
        frames.push(HarmonicDerivative {
            frequency,
            displacement: value,
            direct,
            adjoint,
            amplitude: libm::hypot(value[0], value[1]),
            phase: libm::atan2(value[1], value[0]),
            amplitude_derivative,
            phase_derivative,
            primal_relative_residual: primal.residual.relative,
            tangent_relative_residual: derivative.residual.relative,
        });
    }
    Ok(HarmonicSensitivity { parameter: parameter.clone(), parameter_value: parameter.value(model)?, increment, frames, warnings })
}

/// Bounded fit to complex physical displacements, ordered as real, imaginary
/// pairs at the harmonic step's frequencies. Using Cartesian residuals avoids
/// phase wrapping and amplitude-only ambiguity. Known independent component
/// errors may provide a local conditional confidence interval.
pub fn inverse_harmonic(
    model: &Model,
    step: Option<&str>,
    parameter: &Parameter,
    fit: &super::FrequencyFit,
    observation: HarmonicObservation,
    limits: &Limits,
) -> Result<super::FrequencyEstimate> {
    if !fit.observed.len().is_multiple_of(2) {
        return Err(Error::Invalid("harmonic observations must contain real, imaginary pairs".into()));
    }
    let mut scales = Vec::with_capacity(fit.observed.len());
    for pair in fit.observed.as_chunks::<2>().0 {
        let magnitude = libm::hypot(pair[0], pair[1]);
        if magnitude == 0.0 && fit.standard_deviation.is_none() {
            return Err(Error::Invalid("zero complex observation requires known component standard deviations".into()));
        }
        scales.extend([magnitude, magnitude]);
    }
    super::inverse::fit_scalar(
        model,
        parameter,
        fit,
        &scales,
        "Direct harmonic Cartesian residuals; resonance may introduce multiple local minima.",
        |trial, increment| {
            let result = harmonic_sensitivity(trial, step, parameter, increment, observation, limits)?;
            Ok(result.frames.iter().flat_map(|f| [(f.displacement[0], f.direct[0]), (f.displacement[1], f.direct[1])]).collect())
        },
    )
}
