//! Eigenvalues: the dense symmetric eigenproblem of the Lanczos projection, and the shifted
//! block Lanczos method for the lowest modes of K φ = λ M φ.
//!
//! **Dense.** Householder reduction to tridiagonal form, then the implicit symmetric QR
//! algorithm with Wilkinson shifts, accumulating the transformations (G. H. Golub, C. F. Van
//! Loan, *Matrix Computations*, 4th ed., Johns Hopkins 2013, §8.3, algorithms 8.3.1 and 8.3.3).
//!
//! **Sparse.** The spectral transformation (T. Ericsson, A. Ruhe, *The spectral transformation
//! Lanczos method for the numerical solution of large sparse generalized symmetric eigenvalue
//! problems*, Math. Comp. 35(152):1251–1268, 1980): with a shift σ below the wanted
//! eigenvalues, the operator (K − σM)⁻¹M is symmetric in the M inner product and has eigenvalues
//! θ = 1/(λ − σ), so the lowest λ become the largest, well separated θ. A block of vectors
//! (R. G. Grimes, J. G. Lewis, H. D. Simon, *A shifted block Lanczos algorithm for solving
//! sparse symmetric generalized eigenproblems*, SIAM J. Matrix Anal. Appl. 15(1):228–272, 1994)
//! finds repeated eigenvalues up to the block size, such as the six rigid-body modes of a free
//! body. Every new block is reorthogonalized against all earlier ones, twice (full
//! reorthogonalization, B. N. Parlett, *The Symmetric Eigenvalue Problem*, SIAM 1998, ch. 13),
//! so no spurious copies appear. A Ritz pair has converged when its residual bound
//! ‖B_{j+1} s‖ is below a relative tolerance. The start block is pseudo-random from a fixed seed,
//! so results are reproducible.

use crate::{Error, Result};

/// Eigenvalues (ascending) and eigenvectors (column k of the row-major `vectors`) of a dense
/// symmetric matrix (row-major, n × n).
pub fn symmetric_eigen(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    if a.len() != n * n {
        return Err(Error::Internal("dense eigenproblem size".into()));
    }
    let mut t = a.to_vec();
    let mut q = vec![0.0; n * n];
    for i in 0..n {
        q[i * n + i] = 1.0;
    }
    // Householder tridiagonalization: for each column k, reflect A[k+1.., k] onto e1.
    let mut v = vec![0.0; n];
    let mut p = vec![0.0; n];
    for k in 0..n.saturating_sub(2) {
        let m = k + 1;
        let norm: f64 = (m..n).map(|i| t[i * n + k] * t[i * n + k]).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let x0 = t[m * n + k];
        let alpha = if x0 >= 0.0 { -norm } else { norm };
        for i in m..n {
            v[i] = t[i * n + k];
        }
        v[m] -= alpha;
        let vnorm: f64 = (m..n).map(|i| v[i] * v[i]).sum::<f64>().sqrt();
        if vnorm == 0.0 {
            continue;
        }
        for i in m..n {
            v[i] /= vnorm;
        }
        // p = A22 v, w = p − (vᵀp) v, A22 −= 2 (v wᵀ + w vᵀ).
        for i in m..n {
            p[i] = (m..n).map(|j| t[i * n + j] * v[j]).sum();
        }
        let kk: f64 = (m..n).map(|i| v[i] * p[i]).sum();
        for i in m..n {
            p[i] -= kk * v[i];
        }
        for i in m..n {
            for j in m..n {
                t[i * n + j] -= 2.0 * (v[i] * p[j] + p[i] * v[j]);
            }
        }
        t[m * n + k] = alpha;
        t[k * n + m] = alpha;
        for i in m + 1..n {
            t[i * n + k] = 0.0;
            t[k * n + i] = 0.0;
        }
        // Q = Q (I − 2 v vᵀ).
        for r in 0..n {
            let s: f64 = (m..n).map(|j| q[r * n + j] * v[j]).sum();
            for j in m..n {
                q[r * n + j] -= 2.0 * s * v[j];
            }
        }
    }
    let mut d: Vec<f64> = (0..n).map(|i| t[i * n + i]).collect();
    let mut e: Vec<f64> = (0..n.saturating_sub(1)).map(|i| t[(i + 1) * n + i]).collect();
    tridiagonal_qr(&mut d, &mut e, &mut q, n)?;
    // Sort ascending, carrying the vectors.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| d[a].total_cmp(&d[b]).then(a.cmp(&b)));
    let values = order.iter().map(|&i| d[i]).collect();
    let mut vectors = vec![0.0; n * n];
    for (new, &old) in order.iter().enumerate() {
        for r in 0..n {
            vectors[r * n + new] = q[r * n + old];
        }
    }
    Ok((values, vectors))
}

/// Diagonalize a symmetric tridiagonal matrix (diagonal d, sub-diagonal e) by implicit QR steps
/// with Wilkinson shifts, applying the rotations to the columns of q.
fn tridiagonal_qr(d: &mut [f64], e: &mut [f64], q: &mut [f64], n: usize) -> Result<()> {
    let small = |e: f64, a: f64, b: f64| e.abs() <= f64::EPSILON * (a.abs() + b.abs());
    let mut h = n.saturating_sub(1);
    let mut iterations = 0usize;
    let cap = 30 * n.max(1) + 30;
    while h > 0 {
        if small(e[h - 1], d[h - 1], d[h]) {
            e[h - 1] = 0.0;
            h -= 1;
            continue;
        }
        let mut l = h - 1;
        while l > 0 && !small(e[l - 1], d[l - 1], d[l]) {
            l -= 1;
        }
        if l > 0 {
            e[l - 1] = 0.0;
        }
        iterations += 1;
        if iterations > cap {
            return Err(Error::NotConverged("the tridiagonal QR iteration did not converge".into()));
        }
        // Wilkinson shift from the trailing 2 × 2 block.
        let dd = (d[h - 1] - d[h]) / 2.0;
        let eh = e[h - 1];
        let sign = if dd >= 0.0 { 1.0 } else { -1.0 };
        let mu = d[h] - eh * eh / (dd + sign * dd.hypot(eh));
        let mut x = d[l] - mu;
        let mut z = e[l];
        for k in l..h {
            let r = x.hypot(z);
            let (c, s) = if r == 0.0 { (1.0, 0.0) } else { (x / r, -z / r) };
            if k > l {
                e[k - 1] = c * e[k - 1] - s * z;
            }
            let (a, b, dq) = (d[k], e[k], d[k + 1]);
            d[k] = c * c * a - 2.0 * c * s * b + s * s * dq;
            d[k + 1] = s * s * a + 2.0 * c * s * b + c * c * dq;
            e[k] = c * s * (a - dq) + (c * c - s * s) * b;
            if k + 1 < h {
                z = -s * e[k + 1];
                e[k + 1] *= c;
                x = e[k];
            }
            for row in 0..n {
                let (qp, qq) = (q[row * n + k], q[row * n + k + 1]);
                q[row * n + k] = c * qp - s * qq;
                q[row * n + k + 1] = s * qp + c * qq;
            }
        }
    }
    Ok(())
}

/// A small deterministic pseudo-random generator (xorshift64*, S. Vigna, *An experimental
/// exploration of Marsaglia's xorshift generators, scrambled*, ACM TOMS 42(4), 2016) for
/// Lanczos start vectors.
pub struct Random(u64);

impl Random {
    pub fn new(seed: u64) -> Random {
        Random(seed.max(1))
    }

    /// Uniform in [−1, 1).
    pub fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let v = self.0.wrapping_mul(0x2545_f491_4f6c_dd1d);
        (v >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }
}

/// What the Lanczos iteration looks for, in terms of λ = σ + 1/θ.
#[derive(Debug, Clone, Copy)]
pub struct Wanted {
    /// Modes to return within the range.
    pub count: usize,
    /// Lower and upper bounds on λ (eigenvalues outside are skipped).
    pub min: f64,
    pub max: f64,
}

/// The operators the Lanczos iteration needs.
pub trait Operator: Sync {
    fn n(&self) -> usize;
    /// y = M x.
    fn mass(&self, x: &[f64], y: &mut [f64]);
    /// x ← (K − σM)⁻¹ x.
    fn solve(&self, x: &mut [f64]) -> Result<()>;
    /// The right-hand side the operator solves for, given a basis vector `x` and its M product
    /// `mx`: M x for the spectral transformation (the default); another matrix's product for
    /// other pencils (buckling, [`crate::buckle`]).
    fn operand(&self, _x: &[f64], mx: &[f64]) -> Vec<f64> {
        mx.to_vec()
    }
    /// The eigenvalue of a Ritz value θ of the operator: σ + 1/θ for the spectral
    /// transformation.
    fn eigenvalue(&self, shift: f64, theta: f64) -> f64 {
        shift + 1.0 / theta
    }
}

/// Converged eigenpairs: λ ascending, vectors M-orthonormal.
#[derive(Debug, Clone)]
pub struct Eigenpairs {
    pub values: Vec<f64>,
    pub vectors: Vec<Vec<f64>>,
    /// Converged eigenvalues of negative Ritz values θ (outside what the transformation looks
    /// for: negative buckling load factors), largest |θ| first.
    pub negative: Vec<f64>,
    /// Lanczos vectors used.
    pub basis: usize,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Orthonormalize a block against the basis and itself in the M inner product (classical
/// Gram–Schmidt, twice). A column that vanishes is replaced by a fresh random direction and its
/// diagonal coefficient stays zero. Returns the new vectors, their M products and the
/// upper-triangular coefficients `r[j][i]` of column i on new vector j.
fn orthonormalize(
    op: &impl Operator,
    w: Vec<Vec<f64>>,
    q: &[Vec<f64>],
    mq: &[Vec<f64>],
    rng: &mut Random,
    scratch: &mut [f64],
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = op.n();
    let k = w.len();
    let mut r = vec![vec![0.0; k]; k];
    let mut out: Vec<Vec<f64>> = Vec::with_capacity(k);
    let mut mout: Vec<Vec<f64>> = Vec::with_capacity(k);
    for (i, mut v) in w.into_iter().enumerate() {
        let mut fresh = false;
        for _attempt in 0..4 {
            op.mass(&v, scratch);
            let before = dot(&v, scratch).max(0.0).sqrt();
            for _ in 0..2 {
                for (qj, mqj) in q.iter().zip(mq) {
                    let c = dot(mqj, &v);
                    v.iter_mut().zip(qj).for_each(|(x, y)| *x -= c * y);
                }
                for (j, (qj, mqj)) in out.iter().zip(&mout).enumerate() {
                    let c = dot(mqj, &v);
                    v.iter_mut().zip(qj).for_each(|(x, y)| *x -= c * y);
                    if !fresh && let Some(x) = r.get_mut(j).and_then(|row| row.get_mut(i)) {
                        *x += c;
                    }
                }
            }
            op.mass(&v, scratch);
            let norm = dot(&v, scratch).max(0.0).sqrt();
            if norm > 1e-8 * before && norm > 0.0 {
                if !fresh && let Some(x) = r.get_mut(i).and_then(|row| row.get_mut(i)) {
                    *x = norm;
                }
                v.iter_mut().for_each(|x| *x /= norm);
                mout.push(scratch.iter().map(|x| x / norm).collect());
                out.push(v);
                break;
            }
            fresh = true;
            v = (0..n).map(|_| rng.uniform()).collect();
        }
    }
    (out, mout, r)
}

/// Shift-invert block Lanczos with full reorthogonalization.
pub fn lanczos(op: &impl Operator, shift: f64, wanted: Wanted, block: usize, max_basis: usize, tolerance: f64) -> Result<Eigenpairs> {
    let n = op.n();
    if n == 0 {
        return Err(Error::Invalid("there are no free degrees of freedom".into()));
    }
    let b = block.clamp(1, n);
    let max_basis = max_basis.min(n).max(b);
    let mut rng = Random::new(0x9e37_79b9_7f4a_7c15);
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut mq: Vec<Vec<f64>> = Vec::new();
    // The projected matrix T, grown block by block (row-major in a growing square).
    let mut t: Vec<Vec<f64>> = Vec::new();
    let mut scratch = vec![0.0; n];

    // Start block.
    let start: Vec<Vec<f64>> = (0..b).map(|_| (0..n).map(|_| rng.uniform()).collect()).collect();
    let (mut cur, mut mcur, _) = orthonormalize(op, start, &q, &mq, &mut rng, &mut scratch);
    if cur.len() < b {
        return Err(Error::NotConverged("could not build a Lanczos start block".into()));
    }
    let block_width = b;
    let mut prev_b: Option<Vec<Vec<f64>>> = None;
    let mut next_check = (wanted.count + b).max(2 * b);
    loop {
        fieldcraft_model::progress::check()?;
        let b = cur.len();
        let j0 = q.len();
        // W = Op Q_j.
        let mut w: Vec<Vec<f64>> = Vec::with_capacity(b);
        // The block's solves are independent: one task each.
        for y in fieldcraft_parallel::map(mcur.len(), |i| {
            let mut y = op.operand(cur.get(i).map(Vec::as_slice).unwrap_or_default(), mcur.get(i).map(Vec::as_slice).unwrap_or_default());
            op.solve(&mut y).map(|()| y)
        }) {
            w.push(y?);
        }
        // A_j = Q_jᵀ M W (symmetrized), W −= Q_j A_j (the earlier blocks go in the full
        // reorthogonalization below).
        let mut a = vec![vec![0.0; b]; b];
        for (i, mqi) in mcur.iter().enumerate() {
            for (j, wj) in w.iter().enumerate() {
                a[i][j] = dot(mqi, wj);
            }
        }
        for i in 0..b {
            for j in i + 1..b {
                let s = 0.5 * (a[i][j] + a[j][i]);
                a[i][j] = s;
                a[j][i] = s;
            }
        }
        q.append(&mut cur);
        mq.append(&mut mcur);
        let m = q.len();
        for row in &mut t {
            row.resize(m, 0.0);
        }
        t.resize(m, vec![0.0; m]);
        for i in 0..b {
            for j in 0..b {
                t[j0 + i][j0 + j] = a[i][j];
            }
        }
        if let Some(bp) = &prev_b {
            // T[block j][block j−1] = B_j, T[block j−1][block j] = B_jᵀ.
            let previous_width = bp.first().map_or(0, Vec::len);
            for (i, row) in bp.iter().enumerate() {
                for (j, value) in row.iter().enumerate() {
                    t[j0 + i][j0 - previous_width + j] = *value;
                    t[j0 - previous_width + j][j0 + i] = *value;
                }
            }
        }
        // Full reorthogonalization happens inside the QR against the whole basis.
        let full = m >= max_basis;
        let next_width = block_width.min(max_basis.saturating_sub(m));
        let (next, mnext, bnext) = if full {
            (Vec::new(), Vec::new(), Vec::new())
        } else if next_width == b {
            orthonormalize(op, w, &q, &mq, &mut rng, &mut scratch)
        } else {
            // A final partial block must reach the remaining dimensions. Its coupling
            // is rectangular: next_width rows by the previous block's b columns.
            let (next, mnext, _) = orthonormalize(op, w.iter().take(next_width).cloned().collect(), &q, &mq, &mut rng, &mut scratch);
            let coupling = mnext.iter().map(|v| w.iter().map(|column| dot(v, column)).collect()).collect();
            (next, mnext, coupling)
        };
        // Convergence check on a schedule (the dense eigenproblem costs m³).
        if m >= next_check || full || next.len() < next_width {
            next_check = m + (m / 4).max(b);
            let dense: Vec<f64> = t.iter().flatten().copied().collect();
            let (theta, s) = symmetric_eigen(&dense, m)?;
            // Exhausting a user/memory cap is not convergence. At a truncated cap,
            // check the actual transformed Ritz residual in the M norm. Only a basis
            // spanning all n dimensions may use the vanishing terminal residual.
            let residual_squared = |i: usize| -> Result<f64> {
                if (full && m < n) || (next_width < b && max_basis < n) {
                    let mut x = vec![0.0; n];
                    let mut mx = vec![0.0; n];
                    for (k, (qk, mqk)) in q.iter().zip(&mq).enumerate() {
                        let c = s[k * m + i];
                        x.iter_mut().zip(qk).for_each(|(a, b)| *a += c * b);
                        mx.iter_mut().zip(mqk).for_each(|(a, b)| *a += c * b);
                    }
                    let mut residual = op.operand(&x, &mx);
                    op.solve(&mut residual)?;
                    residual.iter_mut().zip(&x).for_each(|(a, b)| *a -= theta[i] * b);
                    op.mass(&residual, &mut mx);
                    return Ok(dot(&residual, &mx).max(0.0));
                }
                Ok(bnext.iter().map(|row| (0..b).map(|c| row[c] * s[(m - b + c) * m + i]).sum::<f64>().powi(2)).sum())
            };
            // Order by λ ascending = θ descending (θ > 0 below the spectrum).
            let mut idx: Vec<usize> = (0..m).collect();
            idx.sort_by(|&x, &y| theta[y].total_cmp(&theta[x]).then(x.cmp(&y)));
            let mut done = false;
            let mut accepted = Vec::new();
            let mut found = 0;
            for &i in &idx {
                let th = theta[i];
                if !(th > 0.0) {
                    break;
                }
                // Residual bound ‖B_{j+1} s_last‖ with s_last the last block of the Ritz vector.
                let res2 = residual_squared(i)?;
                if res2.sqrt() > tolerance * th {
                    break;
                }
                let lambda = op.eigenvalue(shift, th);
                if lambda > wanted.max {
                    done = true;
                    break;
                }
                if lambda >= wanted.min {
                    accepted.push(i);
                    found += 1;
                    if found == wanted.count {
                        done = true;
                        break;
                    }
                }
            }
            if full && !done {
                if m >= n {
                    done = true;
                } else {
                    return Err(Error::NotConverged(format!(
                        "found {found} of {} modes after {m} Lanczos vectors (Lanczos basis cap); raise the iteration limit or ask for fewer modes",
                        wanted.count
                    )));
                }
            }
            if done {
                // Negative Ritz values that have converged, largest magnitude first.
                let mut negative = Vec::new();
                for &i in idx.iter().rev() {
                    let th = theta[i];
                    if !(th < 0.0) {
                        break;
                    }
                    let res2 = residual_squared(i)?;
                    if res2.sqrt() > tolerance * th.abs() {
                        break;
                    }
                    negative.push(op.eigenvalue(shift, th));
                }
                let mut values = Vec::new();
                let mut vectors = Vec::new();
                for &i in &accepted {
                    values.push(op.eigenvalue(shift, theta[i]));
                    let mut x = vec![0.0; n];
                    for (k, qk) in q.iter().enumerate() {
                        let c = s[k * m + i];
                        x.iter_mut().zip(qk).for_each(|(a, b)| *a += c * b);
                    }
                    op.mass(&x, &mut scratch);
                    let norm = dot(&x, &scratch).max(0.0).sqrt();
                    if norm > 0.0 {
                        x.iter_mut().for_each(|a| *a /= norm);
                    }
                    // Deterministic sign: the largest component positive.
                    let big = x.iter().copied().fold(0.0f64, |m, v| if v.abs() > m.abs() { v } else { m });
                    if big < 0.0 {
                        x.iter_mut().for_each(|a| *a = -*a);
                    }
                    vectors.push(x);
                }
                return Ok(Eigenpairs { values, vectors, negative, basis: m });
            }
        }
        if full {
            return Err(Error::NotConverged("the Lanczos basis is exhausted".into()));
        }
        cur = next;
        mcur = mnext;
        if cur.len() < next_width {
            return Err(Error::NotConverged("the Lanczos block lost rank".into()));
        }
        prev_b = Some(bnext);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_symmetric_eigen() {
        let n = 9;
        let mut rng = Random::new(7);
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = rng.uniform();
                a[i * n + j] = v;
                a[j * n + i] = v;
            }
        }
        let (values, vectors) = symmetric_eigen(&a, n).unwrap();
        assert!(values.windows(2).all(|w| w[0] <= w[1]));
        for k in 0..n {
            for i in 0..n {
                let av: f64 = (0..n).map(|j| a[i * n + j] * vectors[j * n + k]).sum();
                assert!((av - values[k] * vectors[i * n + k]).abs() < 1e-12, "pair {k}");
            }
            let norm: f64 = (0..n).map(|i| vectors[i * n + k].powi(2)).sum();
            assert!((norm - 1.0).abs() < 1e-12);
        }
        assert_eq!(symmetric_eigen(&[3.0], 1).unwrap().0, vec![3.0]);
    }

    /// K and M diagonal.
    struct Diagonal {
        k: Vec<f64>,
        m: Vec<f64>,
        shift: f64,
    }

    impl Operator for Diagonal {
        fn n(&self) -> usize {
            self.k.len()
        }
        fn mass(&self, x: &[f64], y: &mut [f64]) {
            for i in 0..x.len() {
                y[i] = self.m[i] * x[i];
            }
        }
        fn solve(&self, x: &mut [f64]) -> Result<()> {
            for i in 0..x.len() {
                x[i] /= self.k[i] - self.shift * self.m[i];
            }
            Ok(())
        }
    }

    #[test]
    fn lanczos_boundary_cap_does_not_imply_convergence() {
        let op = Diagonal { k: (1..=60).map(f64::from).collect(), m: vec![1.0; 60], shift: 0.0 };
        let wanted = Wanted { count: 1, min: 0.0, max: f64::INFINITY };
        for cap in [6, 13] {
            assert!(lanczos(&op, 0.0, wanted, 6, cap, 1e-12).is_err());
        }
    }

    #[test]
    fn lanczos_boundary_partial_last_block_spans_the_full_space() {
        let op = Diagonal { k: (1..=25).map(f64::from).collect(), m: vec![1.0; 25], shift: 0.0 };
        let wanted = Wanted { count: 25, min: 0.0, max: f64::INFINITY };
        let pairs = lanczos(&op, 0.0, wanted, 6, 25, 1e-12).unwrap();
        assert_eq!(pairs.basis, 25);
        assert_eq!(pairs.values.len(), 25);
        for (i, (lambda, vector)) in pairs.values.iter().zip(&pairs.vectors).enumerate() {
            assert!((lambda - (i + 1) as f64).abs() < 1e-10);
            let residual = vector.iter().zip(&op.k).map(|(x, k)| ((k - lambda) * x).powi(2)).sum::<f64>().sqrt();
            assert!(residual / lambda < 1e-10);
        }
    }

    #[test]
    fn lanczos_finds_repeated_and_zero_eigenvalues() {
        // Six zero eigenvalues, then 1 twice, then distinct values; M not the identity.
        let n = 200;
        let lambdas: Vec<f64> = (0..n)
            .map(|i| {
                if i < 6 {
                    0.0
                } else if i < 8 {
                    1.0
                } else {
                    1.0 + (i - 7) as f64 * 0.5
                }
            })
            .collect();
        let m: Vec<f64> = (0..n).map(|i| 1.0 + (i % 7) as f64 * 0.1).collect();
        let k: Vec<f64> = lambdas.iter().zip(&m).map(|(l, m)| l * m).collect();
        let shift = -1e-3;
        let op = Diagonal { k, m, shift };
        let wanted = Wanted { count: 10, min: f64::NEG_INFINITY, max: f64::INFINITY };
        let pairs = lanczos(&op, shift, wanted, 6, n, 1e-9).unwrap();
        assert_eq!(pairs.values.len(), 10);
        let want = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.5, 2.0];
        for (got, w) in pairs.values.iter().zip(want) {
            assert!((got - w).abs() < 1e-8, "{:?}", pairs.values);
        }
        // M-orthonormal vectors.
        let mut y = vec![0.0; n];
        for a in &pairs.vectors {
            for b in &pairs.vectors {
                op.mass(b, &mut y);
                let d = dot(a, &y);
                assert!((d - if std::ptr::eq(a, b) { 1.0 } else { 0.0 }).abs() < 1e-8);
            }
        }
        // A range: values in [1.2, 3.2] (λ), at most 100.
        let ranged = lanczos(&op, shift, Wanted { count: 100, min: 1.2, max: 3.2 }, 6, n, 1e-9).unwrap();
        assert_eq!(ranged.values.len(), 4, "{:?}", ranged.values);
        assert!(ranged.values.iter().zip([1.5, 2.0, 2.5, 3.0]).all(|(g, w)| (g - w).abs() < 1e-8), "{:?}", ranged.values);
        // Caps: a basis too small to converge is an error, not a hang.
        assert!(lanczos(&op, shift, Wanted { count: 50, min: f64::NEG_INFINITY, max: f64::INFINITY }, 6, 24, 1e-9).is_err());
        // Tiny problems: the whole space.
        let small = Diagonal { k: vec![2.0, 1.0, 3.0], m: vec![1.0; 3], shift: 0.0 };
        let all = lanczos(&small, 0.0, Wanted { count: 5, min: f64::NEG_INFINITY, max: f64::INFINITY }, 6, 100, 1e-9).unwrap();
        assert_eq!(all.values.len(), 3);
        assert!((all.values[0] - 1.0).abs() < 1e-12 && (all.values[2] - 3.0).abs() < 1e-12, "{:?}", all.values);
    }
}
