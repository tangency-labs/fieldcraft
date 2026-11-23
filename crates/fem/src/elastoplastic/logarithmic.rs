//! Lagrangian logarithmic plasticity (C. Miehe, N. Apel, M. Lambrecht, CMAME
//! 191:5383–5425, 2002, doi:10.1016/S0045-7825(02)00438-3): `E = ½ log C`,
//! `C = FᵀF`. The additive return acts on E, with conjugate stress T and tangent D.
//! The chain rule gives `S = L_log(C)[T]` and, for a Green strain increment H,
//! `C_mat[H] = L_log(C)[D[L_log(C)[H]]] + 2 D²log(C)[H,T]`.
//! Both derivatives use divided differences, including their coincident-node limits
//! (R. Bhatia, Matrix Analysis, Springer 1997, §V.3). No eigenvector derivative is needed.

use super::{HISTORY, log_divided_difference, small, value};
use crate::constitutive::{PAIRS, Response, Tangent, push_forward, sandwich};
use crate::dense::{det3, symmetric_eigen};
use crate::plasticity::J2;
use crate::rotation::{Mat3, mul, tensor, transpose, voigt};
use crate::{Error, Result};

/// `log[a,b]`, evaluated without cancellation at equal or nearby positive nodes.
fn first(a: f64, b: f64) -> f64 {
    let mean = 0.5 * a + 0.5 * b;
    let r = (0.5 * a - 0.5 * b) / mean;
    if r.abs() < 0.5 { log_divided_difference(r) / mean } else { (libm::log(a) - libm::log(b)) / (a - b) }
}

/// Symmetric `log[a,b,c]`. Sort so the recurrence divides by the widest gap;
/// near a triple root use the convergent log-series in complete homogeneous
/// polynomials h_k, including the exact limit −1/(2a²).
fn second(a: f64, b: f64, c: f64) -> f64 {
    let mut nodes = [a, b, c];
    nodes.sort_by(f64::total_cmp);
    let [a, b, c] = nodes;
    let m = 0.5 * a + 0.5 * c;
    if (c - a) / m >= 0.01 {
        return (first(a, b) - first(b, c)) / (a - c);
    }
    // Here each scaled node has magnitude below 0.005; twelve terms exceed f64 precision.
    let (u, v, w) = ((a - m) / m, (b - m) / m, (c - m) / m);
    let (s1, s2, s3) = (u + v + w, u * v + u * w + v * w, u * v * w);
    let (mut h0, mut h1, mut h2) = (0.0, 0.0, 1.0);
    let mut sum = -0.5;
    for k in 1..=12 {
        let h = s1 * h2 - s2 * h1 + s3 * h0;
        sum += if k % 2 == 0 { -h } else { h } / (k + 2) as f64;
        (h0, h1, h2) = (h1, h2, h);
    }
    sum / m / m
}

fn engineering(a: &Mat3) -> [f64; 6] {
    let v = voigt(a);
    [v[0], v[1], v[2], 2.0 * v[3], 2.0 * v[4], 2.0 * v[5]]
}

pub(super) fn isochoric(j2: &J2, f: &Mat3, temperature: Option<f64>, committed: &[f64], updated: &mut [f64]) -> Result<Response> {
    let invalid = || Error::Degenerate("non-finite or degenerate logarithmic plastic material point".into());
    let j = det3(f);
    if !j.is_finite() || j <= 0.0 || committed.iter().take(HISTORY).any(|x| !x.is_finite()) {
        return Err(invalid());
    }
    let (l, q) = symmetric_eigen(&mul(&transpose(f), f)).filter(|(l, _)| l.iter().all(|x| x.is_finite() && *x > 0.0)).ok_or_else(invalid)?;
    let qt = transpose(&q);
    let global = |a: &Mat3| mul(&mul(&q, a), &qt);
    let local = |a: &Mat3| mul(&mul(&qt, a), &q);
    let eps: Mat3 = std::array::from_fn(|a| std::array::from_fn(|b| if a == b { 0.5 * libm::log(l[a]) } else { 0.0 }));
    let mut history: [f64; HISTORY] = std::array::from_fn(|k| value(committed, k));
    let r = small(j2, &engineering(&global(&eps)), temperature, committed, &mut history, false)?;
    let t = local(&tensor(&r.stress));
    let d1: Mat3 = std::array::from_fn(|a| std::array::from_fn(|b| first(l[a], l[b])));
    let d2: [[[f64; 3]; 3]; 3] = std::array::from_fn(|a| std::array::from_fn(|b| std::array::from_fn(|c| second(l[a], l[b], l[c]))));
    let apply = |a: &Mat3| -> Mat3 { std::array::from_fn(|i| std::array::from_fn(|j| d1[i][j] * a[i][j])) };
    let s = global(&apply(&t));
    let mut material: Tangent = [[0.0; 6]; 6];
    for (col, &(i, j)) in PAIRS.iter().enumerate() {
        let mut h = [[0.0; 3]; 3];
        h[i][j] = if i == j { 1.0 } else { 0.5 };
        h[j][i] = h[i][j];
        let h = local(&h);
        let de = engineering(&global(&apply(&h)));
        let dt = std::array::from_fn(|a| (0..6).map(|b| r.tangent[a][b] * de[b]).sum());
        let dt = apply(&local(&tensor(&dt)));
        let ds = global(&std::array::from_fn(|a| {
            std::array::from_fn(|b| dt[a][b] + 2.0 * (0..3).map(|k| d2[a][k][b] * (h[a][k] * t[k][b] + t[a][k] * h[k][b])).sum::<f64>())
        }));
        let ds = voigt(&ds);
        for row in 0..6 {
            material[row][col] = ds[row];
        }
    }
    let response = Response { stress: voigt(&mul(&mul(f, &s), &transpose(f))), tangent: sandwich(&push_forward(f), &material) };
    if response.stress.iter().chain(response.tangent.iter().flatten()).chain(history.iter()).any(|x| !x.is_finite()) {
        return Err(invalid());
    }
    for (slot, v) in updated.iter_mut().zip(history) {
        *slot = v;
    }
    Ok(response)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn logarithm_divided_differences_include_repeated_nodes() {
        for a in [0.01, 1.0, 100.0] {
            assert_eq!(first(a, a), 1.0 / a);
            assert!((second(a, a, a) + 0.5 / a / a).abs() < 1e-14 / a / a);
            for delta in [0.0, 1e-12, 1e-7, 1e-3] {
                let b = a * (1.0 + delta);
                let c = a * (1.0 - delta);
                // For symmetric nodes log[a,a(1+δ),a(1−δ)] = ln(1−δ²)/(2a²δ²).
                // Include h₂ and higher terms in the reference, so the tolerance cannot hide them.
                let d2: f64 = delta * delta;
                let expected = if delta == 0.0 { -0.5 } else { (-d2).ln_1p() / (2.0 * d2) } / a / a;
                assert!((second(a, b, c) - expected).abs() < 1e-14 / a / a);
                assert_eq!(second(a, b, c), second(c, a, b));
            }
        }
        assert!((second(1.0, 1.0, 2.0) - (2.0f64.ln() - 1.0)).abs() < 1e-15);
        assert!((second(1.0, 2.0, 2.0) - (0.5 - 2.0f64.ln())).abs() < 1e-15);
    }

    #[test]
    fn logarithm_second_difference_resolves_the_quadratic_series_term() {
        for a in [0.015625, 1.0, 64.0] {
            for delta in [1.0 / 256.0, 1.0 / 512.0, 1.0 / 1024.0] {
                // All nodes are exactly representable. Removing the leading −½ isolates h₂,
                // whose coefficient is −¼, with the next term −δ⁴/6 retained in the reference.
                let value = second(a, a * (1.0 + delta), a * (1.0 - delta)) * a * a;
                let coefficient = (value + 0.5) / (delta * delta);
                let expected = -0.25 - delta * delta / 6.0 - delta.powi(4) / 8.0;
                assert!((coefficient - expected).abs() < 4.0 * f64::EPSILON / (delta * delta));
            }
        }
    }
}
