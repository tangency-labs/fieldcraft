//! Nelson's constrained particular solution followed by the mass-normalization
//! correction. R. B. Nelson (1976), AIAA Journal 14(9), 1201–1205,
//! doi:10.2514/3.7211. The simple-eigenvalue limitation and normalization term
//! are also discussed in NASA TM-89127 (1987), equations (5)–(8).
use crate::{Error, Limits, Result, complex, sparse::Csr};

pub(super) struct Operators<'a> {
    pub k: &'a Csr,
    pub m: &'a Csr,
    pub dk: &'a Csr,
    pub dm: Option<&'a Csr>,
}

/// φ must have unit modal mass. Returns its derivative, the differentiated
/// eigen-equation residual and the differentiated normalization residual.
pub(super) fn derivative(op: Operators<'_>, phi: &[f64], lambda: f64, dlambda: f64, limits: &Limits) -> Result<(Vec<f64>, f64, f64)> {
    let n = op.k.n;
    if phi.len() != n || n == 0 {
        return Err(Error::Invalid("Nelson mode must match nonempty equations".into()));
    }
    let mut mp = vec![0.0; n];
    let mut dmp = vec![0.0; n];
    let mut dkp = vec![0.0; n];
    op.m.mul(phi, &mut mp);
    op.dk.mul(phi, &mut dkp);
    if let Some(dm) = op.dm {
        dm.mul(phi, &mut dmp);
    }
    let rhs: Vec<_> = dkp.iter().zip(&dmp).zip(&mp).map(|((k, m), mass)| -k + lambda * m + dlambda * mass).collect();
    let z = op.k.add_scaled(-lambda, op.m)?;
    let pivot = phi
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))
        .map(|(i, _)| i)
        .ok_or_else(|| Error::Invalid("empty Nelson mode".into()))?;
    let reference = op.k.val.iter().zip(&op.m.val).map(|(k, m)| k.abs() + lambda.abs() * m.abs()).fold(0.0_f64, f64::max);
    if !reference.is_finite() || reference <= 0.0 || phi[pivot] == 0.0 {
        return Err(Error::Invalid("invalid Nelson mode scale".into()));
    }
    let mut constrained = z.clone();
    for row in 0..n {
        for at in constrained.ptr[row]..constrained.ptr[row + 1] {
            let col = constrained.idx[at];
            if row == pivot || col == pivot {
                constrained.val[at] = if row == col { reference } else { 0.0 };
            }
        }
    }
    let mut load = complex::Vector { real: rhs.clone(), imaginary: vec![0.0; n] };
    load.real[pivot] = 0.0;
    let zero = Csr { val: vec![0.0; constrained.val.len()], ..constrained.clone() };
    let mut derivative = complex::solve(&constrained, &zero, &load, limits)?.displacement.real;
    let mass_change: f64 = phi.iter().zip(&dmp).map(|(p, m)| p * m).sum();
    let projection: f64 = derivative.iter().zip(&mp).map(|(d, m)| d * m).sum();
    let alpha = -0.5 * mass_change - projection;
    for (d, p) in derivative.iter_mut().zip(phi) {
        *d += alpha * p;
    }
    let mut zd = vec![0.0; n];
    z.mul(&derivative, &mut zd);
    let error = zd.iter().zip(&rhs).map(|(a, b)| (a - b).abs()).fold(0.0_f64, f64::max);
    let scale = dkp.iter().zip(&dmp).zip(&mp).map(|((k, m), p)| k.abs() + lambda.abs() * m.abs() + dlambda.abs() * p.abs()).fold(0.0_f64, f64::max)
        + reference * derivative.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let residual = if scale > 0.0 { error / scale } else { error };
    let normalization = 2.0 * derivative.iter().zip(&mp).map(|(d, m)| d * m).sum::<f64>() + mass_change;
    if derivative.iter().any(|v| !v.is_finite()) || !residual.is_finite() || !normalization.is_finite() || residual > 1e-7 {
        return Err(Error::NotConverged(format!("Nelson differentiated eigen-equation residual {residual}")));
    }
    Ok((derivative, residual, normalization))
}
