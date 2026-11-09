//! Shared total-Lagrangian rod force and consistent material/stress tangent.
//! Derived from the St Venant–Kirchhoff strain energy, as in the nonlinear
//! transient formulation (Bathe & Cimento 1980, doi:10.1016/0045-7825(80)90051-1).
use crate::{Error, Result};
fn vector(v: &[f64], n: usize) -> Result<()> {
    if v.len() != n || v.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("finite rod vector".into()));
    }
    Ok(())
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

/// St Venant–Kirchhoff rod2, total Lagrangian: E_G=(d.d-L²)/(2L²).
/// f = EA E_G d/L; df/dd=EA/L [d⊗d/L²+E_G I].
pub(crate) fn evaluate(x: &[[f64; 3]], u: &[f64], ea: f64) -> Result<(Vec<f64>, Vec<f64>, f64)> {
    if x.len() != 2 || u.len() != 6 || !ea.is_finite() || ea <= 0. {
        return Err(Error::Invalid("finite rod2 data".into()));
    }
    vector(u, 6)?;
    let dx: [f64; 3] = std::array::from_fn(|i| x[1][i] - x[0][i]);
    let l2 = dot(&dx, &dx);
    let l = libm::sqrt(l2);
    let d: [f64; 3] = std::array::from_fn(|i| dx[i] + u[i + 3] - u[i]);
    let d2 = dot(&d, &d);
    if !l2.is_finite() || l2 <= 0. || !d2.is_finite() || d2 <= 1e-20 * l2 {
        return Err(Error::Invalid("collapsed finite rod2".into()));
    }
    let strain = (d2 - l2) / (2. * l2);
    let f: Vec<_> = (0..6).map(|i| if i < 3 { -1. } else { 1. } * ea / l * strain * d[i % 3]).collect();
    let k: Vec<_> = (0..36)
        .map(|i| {
            let r = i / 6;
            let c = i % 6;
            let sign = if (r < 3) == (c < 3) { 1. } else { -1. };
            sign * ea / l * (d[r % 3] * d[c % 3] / l2 + if r % 3 == c % 3 { strain } else { 0. })
        })
        .collect();
    let energy = 0.5 * ea * l * strain * strain;
    vector(&f, 6)?;
    vector(&k, 36)?;
    if !energy.is_finite() {
        return Err(Error::Invalid("rod energy overflow".into()));
    }
    Ok((f, k, energy))
}
