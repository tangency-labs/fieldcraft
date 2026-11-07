//! Load-dependent residual Ritz vectors, formed with the unshifted static factor.
//!
//! Modal force truncation and its static correction follow SAND2025-03995,
//! Sierra/SD Theory Manual 5.24, §2.22. Successive K⁻¹M vectors follow Wilson,
//! Yuan & Dickens, *Earthquake Engineering & Structural Dynamics* 10 (1982),
//! 813–821, DOI 10.1002/eqe.4290100606. We retain the physical modes and twice
//! M-orthogonalize each additional vector against the entire retained basis.
//! The recurrence is also given in A. A. Abdallah (1990), *Dynamic Substructuring
//! by the Boundary Flexibility Vector Method of Component Mode Synthesis*,
//! §2.3, pp. 10–12, NASA NTRS 19900013497.
//! The two load quadratures seed separate chains. No frequency or damping is
//! changed, and the original-equation residual remains the convergence check.

use super::{Factor, Matrices, Vector, dot};
use crate::{Error, Result};

pub(super) fn augment(matrices: &Matrices<'_>, factor: &Factor, basis: &mut Vec<Vec<f64>>, force: &Vector, cycles: usize) -> Result<()> {
    let n = matrices.stiffness.n;
    if force.real.len() != n || force.imaginary.len() != n || force.real.iter().chain(&force.imaginary).any(|x| !x.is_finite()) {
        return Err(Error::Invalid("residual-vector load must match the modal equations and be finite".into()));
    }
    for load in [&force.real, &force.imaginary] {
        let scale = load.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        if scale == 0.0 {
            continue;
        }
        let mut next: Vec<f64> = load.iter().map(|x| x / scale).collect();
        let mut mass = vec![0.0; n];
        for _ in 0..cycles {
            fieldcraft_model::progress::check()?;
            if basis.len() == n {
                break;
            }
            factor.solve(&mut next)?;
            matrices.mass.mul(&next, &mut mass);
            let before = dot(&next, &mass);
            if !before.is_finite() || before < 0.0 {
                return Err(Error::NotConverged("residual vector has invalid mass norm".into()));
            }
            if before == 0.0 && next.iter().any(|x| *x != 0.0) {
                return Err(Error::NotSupported("massless static residual vectors cannot be mass-normalized".into()));
            }
            // Two full passes: only one sparse mass product per pass.
            for _ in 0..2 {
                matrices.mass.mul(&next, &mut mass);
                for phi in basis.iter() {
                    let coefficient = dot(phi, &mass);
                    for (x, p) in next.iter_mut().zip(phi) {
                        *x -= coefficient * p;
                    }
                }
            }
            matrices.mass.mul(&next, &mut mass);
            let norm = dot(&next, &mass);
            if !norm.is_finite() || norm < 0.0 {
                return Err(Error::NotConverged("residual vector lost a positive mass norm".into()));
            }
            // Relative linear-dependence test, not a physical response tolerance.
            if norm <= 1e-24 * before {
                break;
            }
            for x in &mut next {
                *x /= norm.sqrt();
            }
            matrices.mass.mul(&next, &mut mass);
            basis.push(next);
            next = mass.clone();
        }
    }
    Ok(())
}
