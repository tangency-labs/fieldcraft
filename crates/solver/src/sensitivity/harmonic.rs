//! Physical harmonic parameter derivatives at fixed excitation frequencies.
//! Central differences apply to assembled operators, including material damping
//! and prescribed-motion elimination, not independently solved responses.
use fieldcraft_fem::DofKind;
use fieldcraft_model::{Model, Procedure, design::Parameter, harmonic::Method};
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

/// Differentiate direct harmonic equilibrium with a fixed physical observation.
/// Modal truncation and parameter-dependent observation weights are unsupported.
/// A zero observation refuses amplitude/phase derivatives rather than inventing them.
pub fn harmonic_sensitivity(
    model: &Model,
    step: Option<&str>,
    parameter: &Parameter,
    increment: f64,
    observation: HarmonicObservation,
    limits: &Limits,
) -> Result<HarmonicSensitivity> {
    let (plus, minus) = perturb(model, parameter, increment)?;
    let resolved = model.resolve()?;
    let step = model.select_step(step, |p| matches!(p, Procedure::Harmonic(_)), "harmonic")?;
    let Procedure::Harmonic(h) = &step.procedure else {
        return Err(Error::Internal("harmonic sensitivity step".into()));
    };
    if h.method != Method::Direct {
        return Err(Error::NotSupported("harmonic sensitivity requires direct response; modal-basis derivatives are not implemented".into()));
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
    let p = assemble(&plus)?;
    let m = assemble(&minus)?;
    let mut weights = vec![0.0; mesh.dofs()];
    *weights.get_mut(dof).ok_or_else(|| Error::Internal("harmonic observation dof".into()))? = 1.0;
    let weights = Vector { real: red.restrict(&weights), imaginary: vec![0.0; red.neq()] };
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

        let (ap, bp, fp) = dynamic(&p)?;
        let (am, bm, fm) = dynamic(&m)?;
        let da = matrix_difference(&ap, &am, increment)?;
        let db = matrix_difference(&bp, &bm, increment)?;
        let df = Vector { real: difference(&fp.real, &fm.real, increment)?, imaginary: difference(&fp.imaginary, &fm.imaginary, increment)? };
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
