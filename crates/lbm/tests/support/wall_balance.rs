//! Shared test-only wall-stress diagnostics; no production API or numerical changes.
//!
//! Momentum exchange follows Mei et al., PRE 65:041203 (2002). The Musker request is
//! the inference of Malaspinas & Sagaut, JCP 275:25–40 (2014), distinct from the actual
//! impulse delivered by reconstruction and wall links. All inputs here use one consistent
//! system of units; `closed_step` takes lattice momenta/impulses and lattice face area.

use fieldcraft_lbm::{Error, Result};

const UNIT_NORMAL_TOLERANCE: f64 = 1e-12; // unit-factors: dimensionless squared-normal-length tolerance, not a conversion

/// One closed-box step. Positive traction is momentum transferred from fluid to walls.
/// Open-face fluxes, moving control volumes, and omitted body forces invalidate this budget.
#[derive(Debug, Clone, PartialEq)]
pub struct ClosedStep {
    pub requested: [f64; 3],
    pub applied: [f64; 3],
    /// Independently inferred from supplied body impulse minus stored-momentum change.
    pub carried: [f64; 3],
    pub momentum_residual: [f64; 3],
    pub momentum_scale: f64,
    pub absolute_roundoff_floor: f64,
    pub budget_ok: bool,
}

/// A documented scale and absolute floor avoid dividing by vanishing momentum at rest.
/// Scale = max of component magnitudes of both stored momenta, body and wall impulses.
/// Residual gate = relative_tolerance * scale + 64 ε max(mass,1) in lattice units.
/// `epsilon` is f64::EPSILON for the host sums on both backends; tolerances are 1e-12 / 1e-5.
#[allow(clippy::too_many_arguments)]
pub fn closed_step(
    before: [f64; 3],
    after: [f64; 3],
    body_impulse: [f64; 3],
    wall_impulse: [f64; 3],
    requested: [f64; 3],
    area: f64,
    mass: f64,
    relative_tolerance: f64,
    epsilon: f64,
) -> Result<ClosedStep> {
    if ![area, mass, relative_tolerance, epsilon].iter().all(|x| x.is_finite() && *x > 0.0)
        || !before.iter().chain(&after).chain(&body_impulse).chain(&wall_impulse).chain(&requested).all(|x| x.is_finite())
    {
        return Err(Error::Invalid("wall balance needs finite inputs and positive area, mass and tolerances".into()));
    }
    let carried = [0, 1, 2].map(|a| body_impulse[a] - (after[a] - before[a]));
    let residual = [0, 1, 2].map(|a| carried[a] - wall_impulse[a]);
    let scale = before.iter().chain(&after).chain(&body_impulse).chain(&wall_impulse).map(|x| x.abs()).fold(0.0, f64::max);
    let floor = 64.0 * epsilon * mass.max(1.0);
    let applied = wall_impulse.map(|x| x / area);
    let carried = carried.map(|x| x / area);
    let limit = relative_tolerance * scale + floor;
    if !limit.is_finite() || !applied.iter().chain(&carried).chain(&residual).all(|x| x.is_finite()) {
        return Err(Error::Invalid("wall balance overflow".into()));
    }
    Ok(ClosedStep {
        requested,
        applied,
        carried,
        momentum_residual: residual,
        momentum_scale: scale,
        absolute_roundoff_floor: floor,
        budget_ok: residual.iter().all(|r| r.abs() <= limit),
    })
}

/// Recover the model's requested kinematic shear from the target speed at its first cell.
/// Musker's profile at that cell uses the very same friction velocity as its matching cell.
/// This is an inference, never the stress the populations actually carry. Only flat, stationary
/// halfway walls with known normal/area are measured by the channel studies using this helper.
pub fn musker_request(target: [f64; 3], normal: [f64; 3], distance: f64, nu: f64) -> Result<[f64; 3]> {
    if !target.iter().chain(&normal).all(|x| x.is_finite()) || !distance.is_finite() || distance <= 0.0 || !nu.is_finite() || nu <= 0.0 {
        return Err(Error::Invalid("Musker request needs finite targets, positive distance and viscosity".into()));
    }
    let normal_length = normal.iter().map(|x| x * x).sum::<f64>();
    if (normal_length - 1.0).abs() > UNIT_NORMAL_TOLERANCE {
        return Err(Error::Invalid("wall normal must be a unit vector".into()));
    }
    let un = target.iter().zip(normal).map(|(u, n)| u * n).sum::<f64>();
    let tangent = [0, 1, 2].map(|a| target[a] - un * normal[a]);
    let speed = tangent.iter().map(|v| v * v).sum::<f64>().sqrt();
    if !speed.is_finite() {
        return Err(Error::Invalid("Musker request speed overflow".into()));
    }
    let friction = fieldcraft_lbm::wall::musker_friction(speed, distance, nu);
    if speed > 0.0 && !(friction > 0.0 && friction.is_finite()) {
        return Err(Error::Invalid("could not infer a finite positive Musker request".into()));
    }
    let stress = if speed > 0.0 { tangent.map(|v| v / speed * friction * friction) } else { [0.0; 3] };
    if !stress.iter().all(|x| x.is_finite()) {
        return Err(Error::Invalid("Musker request overflow".into()));
    }
    Ok(stress)
}

/// Carried drag is tested against its reference independently of the wall model's inference.
#[derive(Debug, Clone, PartialEq)]
pub struct DragValidation {
    pub requested_over_carried: f64,
    pub carried_relative_error: f64,
    pub carried_drag_ok: bool,
    // Proposed stress-transfer contract: request and carried stress within 1%.
    pub stress_agreement_ok: bool,
}

/// Compare carried drag with its independent reference using the supplied relative tolerance,
/// and requested/carried stress against the proposed 1% agreement gate. Both checks are
/// independent; passing the wall model's own inferred drag cannot satisfy the carried check.
/// Returns an error for non-positive/non-finite inputs or overflowing ratios.
pub fn validate_drag(requested: f64, carried: f64, reference: f64, tolerance: f64) -> Result<DragValidation> {
    if ![requested, carried, reference, tolerance].iter().all(|x| x.is_finite() && *x > 0.0) {
        return Err(Error::Invalid("drag validation needs positive finite stresses, reference and tolerance".into()));
    }
    let ratio = requested / carried;
    let error = carried / reference - 1.0;
    if !ratio.is_finite() || !error.is_finite() {
        return Err(Error::Invalid("drag validation overflow".into()));
    }
    Ok(DragValidation {
        requested_over_carried: ratio,
        carried_relative_error: error,
        carried_drag_ok: error.abs() <= tolerance,
        stress_agreement_ok: (ratio - 1.0).abs() <= 0.01,
    })
}

/// Stored momentum from the actual populations, without the velocity's half-force shift
/// or the wall-target substitution used for display. Fixed cell/direction/component order.
pub fn population_momentum(flow: &fieldcraft_lbm::Flow) -> [f64; 3] {
    let mut momentum = [0.0; 3];
    for cell in (0..flow.grid().len()).filter(|c| flow.is_fluid(*c)) {
        for (q, direction) in fieldcraft_lbm::lattice::C.iter().enumerate() {
            let f = flow.population(q, cell);
            for (value, c) in momentum.iter_mut().zip(direction) {
                *value += f64::from(*c) * f;
            }
        }
    }
    momentum
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn measured_n10_discrepancy_remains_a_failed_validation() {
        // FieldCraft cb976117, original Re_tau950/N10 audit published in #786.
        let ratio = 1.03713 / 0.83318;
        let r = validate_drag(0.0043803 * ratio, 0.0043803, 0.00509, 0.05).unwrap();
        assert!((r.requested_over_carried - 1.2448).abs() < 0.0001);
        assert!((r.carried_relative_error + 0.1394).abs() < 0.001);
        assert!(!r.carried_drag_ok);
        assert!(!r.stress_agreement_ok);
    }
    #[test]
    fn inverse_target_recovers_matching_cell_request() {
        for speed in [0.0, 0.001, 0.01, 0.1] {
            let (target, _) = fieldcraft_lbm::wall::wall_law(speed, 0.5, 1.5, 0.001);
            let r = musker_request([target, 0.0, 0.0], [0.0, 1.0, 0.0], 0.5, 0.001).unwrap();
            let ut = fieldcraft_lbm::wall::musker_friction(speed, 1.5, 0.001);
            assert!((r[0] - ut * ut).abs() < 1e-12, "{speed}: {r:?} vs {}", ut * ut);
        }
    }
    #[test]
    fn budget_detects_missing_impulse_and_handles_rest() {
        let report =
            |wall| closed_step([0.0; 3], [0.01, 0.0, 0.0], [0.03, 0.0, 0.0], wall, [0.012, 0.0, 0.0], 2.0, 10.0, 1e-12, f64::EPSILON).unwrap();
        let r = report([0.02, 0.0, 0.0]);
        assert!(r.budget_ok);
        assert!((r.carried[0] - 0.01).abs() < 1e-12);
        assert!(!report([0.01, 0.0, 0.0]).budget_ok);
        assert!(closed_step([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], 1.0, 1.0, 1e-12, f64::EPSILON).unwrap().budget_ok);
        assert!(closed_step([f64::NAN; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], 1.0, 1.0, 1e-12, f64::EPSILON).is_err());
    }
}
