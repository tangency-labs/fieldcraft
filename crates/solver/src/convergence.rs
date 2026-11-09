//! Diagnostics in the actual, unshifted constrained frequency pencil, K φ = λ M φ.
//! A projected Ritz residual is not this physical-equation residual.
use crate::{Error, Result, sparse::Csr};

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ModeResidual {
    pub residual_norm: f64,
    pub stiffness_action_norm: f64,
    /// ||Kφ − λMφ||₂ / ||Kφ||₂; undefined at a zero denominator.
    pub relative_residual: Option<f64>,
    pub undefined_reason: Option<String>,
    /// ||Kφ|| <= sqrt(eps) ||K||_inf ||φ||: the quotient is ill-conditioned,
    /// particularly for rigid-body modes. This is a warning, not a convergence test.
    pub near_zero_stiffness_action: bool,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ModalConvergence {
    /// Actual applications of the Lanczos block recurrence, not the vector count.
    pub iterations: usize,
    pub basis_vectors: usize,
    pub basis_limit: usize,
    /// The effective basis limit (user/memory/dimension) was reached.
    /// Inspect full_space_reached to distinguish a truncated cap from the full space.
    pub basis_limit_reached: bool,
    pub full_space_reached: bool,
    /// Same order as the returned frequency modes; not buckling-pencil residuals.
    pub modes: Vec<ModeResidual>,
}

/// Evaluate a mode using the actual unshifted K and selected M after all constraint
/// reductions, including the terminal tangent in a prestressed solve. Do not pass
/// K−σM or expand φ into held equations (whose reactions are not residual equations).
pub fn frequency_residual(k: &Csr, m: &Csr, lambda: f64, phi: &[f64]) -> Result<ModeResidual> {
    if k.n != m.n || phi.len() != k.n || !lambda.is_finite() {
        return Err(Error::Invalid("frequency residual dimensions or eigenvalue".into()));
    }
    let mut kp = vec![0.0; k.n];
    let mut mp = vec![0.0; k.n];
    k.mul(phi, &mut kp);
    m.mul(phi, &mut mp);
    let norm = |v: &[f64]| v.iter().fold(0.0f64, |a, b| a.hypot(*b));
    let denominator = norm(&kp);
    let residual = kp.iter().zip(mp).fold(0.0f64, |a, (k, m)| a.hypot(k - lambda * m));
    let scale = (0..k.n).map(|r| k.row(r).1.iter().map(|v| v.abs()).sum::<f64>()).fold(0.0f64, f64::max);
    if !denominator.is_finite() || !residual.is_finite() || !scale.is_finite() {
        return Err(Error::Invalid("non-finite frequency residual".into()));
    }
    let quotient = residual / denominator;
    let (relative_residual, undefined_reason) = if denominator == 0.0 {
        (None, Some("zero_stiffness_action".into()))
    } else if !quotient.is_finite() {
        (None, Some("ratio_overflow".into()))
    } else {
        (Some(quotient), None)
    };
    Ok(ModeResidual {
        residual_norm: residual,
        stiffness_action_norm: denominator,
        relative_residual,
        undefined_reason,
        near_zero_stiffness_action: denominator <= f64::EPSILON.sqrt() * scale * norm(phi),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn diag(v: &[f64]) -> Csr {
        Csr { n: v.len(), ptr: (0..=v.len()).collect(), idx: (0..v.len()).collect(), val: v.to_vec() }
    }
    #[test]
    fn actual_pencil_residual_and_rigid_boundary() {
        let k = diag(&[0.0, 8.0, 12.0]);
        let m = diag(&[2.0, 2.0, 3.0]);
        let exact = frequency_residual(&k, &m, 4.0, &[0.0, 1.0, 1.0]).unwrap();
        assert_eq!(exact.relative_residual, Some(0.0));
        let wrong = frequency_residual(&k, &m, 3.0, &[0.0, 1.0, 1.0]).unwrap();
        assert!((wrong.relative_residual.unwrap() - 0.25).abs() < 1e-15);
        let rigid = frequency_residual(&k, &m, 0.0, &[1.0, 0.0, 0.0]).unwrap();
        assert_eq!(rigid.relative_residual, None);
        assert_eq!(rigid.undefined_reason.as_deref(), Some("zero_stiffness_action"));
        assert!(rigid.near_zero_stiffness_action);
        let inconsistent = frequency_residual(&k, &m, 1.0, &[1.0, 0.0, 0.0]).unwrap();
        assert_eq!(inconsistent.relative_residual, None);
        assert_eq!(inconsistent.residual_norm, 2.0);
        assert!(frequency_residual(&k, &m, 4.0, &[1.0]).is_err());
    }
}
