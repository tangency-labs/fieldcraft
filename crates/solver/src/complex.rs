//! Complex symmetric linear systems, with native and real-equivalent direct factors.
//!
//! With `Re(u exp(i ω t))`, equilibrium is `(A+iB)u=f`, where
//! `A=K−ω²M`, `B=ωC+G`. Transpose symmetry is required, not Hermitian symmetry.
//! [`solve`] uses native complex supernodal LDLᵀ with 1×1/2×2 rook pivots
//! ([`native`], Bunch 1982 and Ashcraft et al. 1998). The existing Cholesky
//! symbolic structure and deterministic subtree schedule are reused. It stores
//! two real slots per original factor entry, compared with approximately four
//! in the doubled real formulation. Both components count toward the cap.
//!
//! [`solve_equivalent`] retains `[[A,−B],[−B,−A]] [Ur,Ui]=[Fr,−Fi]`, with
//! interleaved pairs, as a separate arithmetic check. Reference: Day and Heroux,
//! *Solving Complex-Valued Linear Systems via Equivalent Real Formulations*,
//! SIAM J. Sci. Comput. 23(2):480–498 (2001),
//! <https://doi.org/10.1137/S1064827500372262>. [`crate::ldlt`] supplies its real
//! pivoted factor. No normal equations, artificial damping or frequency shifts
//! are introduced. For B=0, both paths use one n-dimensional real factor.
//!
//! Symbolic capacity is checked before numeric allocation. A failed local
//! pivot is a refusal, not proof of physical resonance: neither factor delays
//! pivots to a later supernode. Refinement checks the original complex equations.
//!
//! Iterative refinement uses the original complex residual, with compensated
//! summation (N. J. Higham, *The Accuracy of Floating Point Summation*, SIAM J.
//! Sci. Comput. 14(4):783–799, 1993, <https://doi.org/10.1137/0914050>), and the
//! existing factor as correction solver (Higham 2002, §12.2). This also refines
//! real-only Cholesky solves; the equilibrium threshold is never relaxed.

use serde::{Deserialize, Serialize};

use crate::cholesky::{FactorError, UpperColumns};
use crate::ldlt::{Kind, Symmetric};
use crate::sparse::Csr;
use crate::{Error, Limits, Result};

mod native;

/// Both parts of a phasor, in the same degree-of-freedom order.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Vector {
    pub real: Vec<f64>,
    pub imaginary: Vec<f64>,
}

impl Vector {
    pub fn zeros(n: usize) -> Self {
        Self { real: vec![0.0; n], imaginary: vec![0.0; n] }
    }

    fn check(&self, n: usize) -> Result<()> {
        if self.real.len() != n || self.imaginary.len() != n {
            return Err(Error::Invalid("complex vector length does not match the matrix".into()));
        }
        if self.real.iter().chain(&self.imaginary).any(|v| !v.is_finite()) {
            return Err(Error::Invalid("complex vector contains a non-finite value".into()));
        }
        Ok(())
    }

    /// Component amplitudes (not the maximum length of the vector during a cycle).
    pub fn amplitudes(&self) -> Vec<f64> {
        self.real.iter().zip(&self.imaginary).map(|(a, b)| libm::hypot(*a, *b)).collect()
    }

    /// Component phases in radians, in [−π, π]; zero has the conventional phase 0.
    pub fn phases(&self) -> Vec<f64> {
        self.real.iter().zip(&self.imaginary).map(|(a, b)| if *a == 0.0 && *b == 0.0 { 0.0 } else { libm::atan2(*b, *a) }).collect()
    }
}

/// Residuals in the original complex equations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Residual {
    /// Euclidean norm of the algebraic imbalance in the supplied equation scaling.
    /// This has no single physical unit when translation and rotation are mixed.
    pub norm: f64,
    /// Infinity-norm ratio ‖f − Z u‖ / ‖f‖. A zero load and zero residual give zero.
    pub relative: f64,
    /// Infinity-norm backward error ‖f − Z u‖ / (‖Z‖ ‖u‖ + ‖f‖).
    pub backward: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    pub displacement: Vector,
    pub residual: Residual,
    /// Stored real entries in the numeric factor, excluding indices and workspace.
    pub factor_entries: usize,
    /// n for real/native complex factors, 2n for the real symmetric equivalent.
    pub factor_equations: usize,
}

/// Validate the CSR shape before it reaches numerical kernels. Input-derived indices
/// are checked even though model assembly normally constructs this representation.
pub(crate) fn validate(a: &Csr, limits: &Limits) -> Result<()> {
    if a.n == 0 || a.n > limits.max_dofs {
        return Err(Error::Limit(format!("complex matrix has {} equations; allowed 1..={}", a.n, limits.max_dofs)));
    }
    if a.n.checked_add(1) != Some(a.ptr.len())
        || a.ptr.first() != Some(&0)
        || a.ptr.last() != Some(&a.idx.len())
        || a.idx.len() != a.val.len()
        || a.ptr.windows(2).any(|w| w[0] > w[1])
    {
        return Err(Error::Invalid("invalid complex matrix CSR pointers or lengths".into()));
    }
    for r in 0..a.n {
        let (cols, vals) = a.row(r);
        if cols.iter().any(|c| *c >= a.n)
            || cols.windows(2).any(|w| w[0] >= w[1])
            || cols.binary_search(&r).is_err()
            || vals.iter().any(|v| !v.is_finite())
        {
            return Err(Error::Invalid(format!("invalid complex matrix row {r}: sorted columns, diagonal and finite entries required")));
        }
    }
    for r in 0..a.n {
        let (cols, vals) = a.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            let (mirror, values) = a.row(c);
            if mirror.binary_search(&r).ok().and_then(|at| values.get(at)).copied() != Some(v) {
                return Err(Error::Invalid(format!("complex matrix is not transpose symmetric at ({r}, {c})")));
            }
        }
    }
    Ok(())
}

fn pair(a: &Csr, b: &Csr, limits: &Limits) -> Result<()> {
    validate(a, limits)?;
    validate(b, limits)?;
    if a.n != b.n || a.ptr != b.ptr || a.idx != b.idx {
        return Err(Error::Invalid("complex matrix real and imaginary parts must share a pattern".into()));
    }
    Ok(())
}

/// Assemble the real equivalent in interleaved order. Explicit zero entries in
/// every 2 × 2 block preserve the symbolic pivot pairing.
fn equivalent(a: &Csr, b: &Csr, limits: &Limits) -> Result<Csr> {
    let limit = || Error::Limit("real equivalent complex matrix size overflow".into());
    let n = a.n.checked_mul(2).ok_or_else(limit)?;
    let entries = a.val.len().checked_mul(4).ok_or_else(limit)?;
    let lower_entries = a.val.len().checked_mul(2).and_then(|v| v.checked_add(a.n)).ok_or_else(limit)?;
    if lower_entries > limits.max_factor_entries {
        return Err(Error::factor_limit(lower_entries, limits.max_factor_entries));
    }
    let mut ptr = Vec::new();
    let mut idx = Vec::new();
    let mut val = Vec::new();
    ptr.try_reserve_exact(n.checked_add(1).ok_or_else(limit)?).map_err(|_| limit())?;
    idx.try_reserve_exact(entries).map_err(|_| limit())?;
    val.try_reserve_exact(entries).map_err(|_| limit())?;
    ptr.push(0);
    for r in 0..a.n {
        let (cols, ar) = a.row(r);
        let (_, br) = b.row(r);
        for imaginary in [false, true] {
            for ((&c, &x), &y) in cols.iter().zip(ar).zip(br) {
                idx.extend([2 * c, 2 * c + 1]);
                val.extend(if imaginary { [-y, -x] } else { [x, -y] });
            }
            ptr.push(idx.len());
        }
    }
    Ok(Csr { n, ptr, idx, val })
}

fn factor_error(e: FactorError, interleaved: bool) -> Error {
    match e {
        FactorError::NotPositiveDefinite { column } | FactorError::Singular { column } => {
            let equation = if interleaved { column / 2 } else { column };
            Error::Singular(format!(
                "harmonic factor has no acceptable pivot at equation {equation}; check restraints, undamped resonance or local pivoting"
            ))
        }
        FactorError::TooLarge { entries, limit } => Error::factor_limit(entries, limit),
        FactorError::OutOfMemory { entries } => Error::Limit(format!("cannot allocate harmonic factor with {entries} entries")),
        FactorError::Internal(message) => Error::Internal(message),
        FactorError::Cancelled => fieldcraft_model::Error::Cancelled.into(),
    }
}

/// Solve (A + i B) u = f, with a bounded native factorization and at most three
/// refinement steps. The caller supplies A and B at the requested frequency.
pub fn solve(a: &Csr, b: &Csr, rhs: &Vector, limits: &Limits) -> Result<Solution> {
    solve_using(a, b, rhs, limits, true)
}

/// The interleaved real equivalent, retained as an independent arithmetic path
/// for checking the native complex factor. It uses the same real-entry cap.
pub fn solve_equivalent(a: &Csr, b: &Csr, rhs: &Vector, limits: &Limits) -> Result<Solution> {
    solve_using(a, b, rhs, limits, false)
}

fn solve_using(a: &Csr, b: &Csr, rhs: &Vector, limits: &Limits, native: bool) -> Result<Solution> {
    fieldcraft_model::progress::check()?;
    pair(a, b, limits)?;
    rhs.check(a.n)?;
    let real_only = b.val.iter().all(|v| *v == 0.0);
    let (displacement, factor_entries, factor_equations) = if real_only {
        let f = Symmetric::new(&UpperColumns::from_symmetric(a), limits.max_factor_entries, Kind::Auto).map_err(|e| factor_error(e, false))?;
        let mut u = rhs.clone();
        f.solve(&mut u.real)?;
        f.solve(&mut u.imaginary)?;
        let u = refine(a, b, rhs, u, |r| {
            f.solve(&mut r.real)?;
            f.solve(&mut r.imaginary)
        })?;
        (u, f.entries(), a.n)
    } else if native {
        let f = native::Factor::new(a, b, limits.max_factor_entries).map_err(|e| factor_error(e, false))?;
        let mut u = rhs.clone();
        f.solve(&mut u)?;
        let u = refine(a, b, rhs, u, |r| f.solve(r))?;
        (u, f.entries(), a.n)
    } else {
        let matrix = equivalent(a, b, limits)?;
        let f =
            Symmetric::new(&UpperColumns::from_symmetric(&matrix), limits.max_factor_entries, Kind::Indefinite).map_err(|e| factor_error(e, true))?;
        let mut u: Vec<f64> = rhs.real.iter().zip(&rhs.imaginary).flat_map(|(r, i)| [*r, -*i]).collect();
        f.solve(&mut u)?;
        let mut out = Vector::zeros(a.n);
        for ((r, i), values) in out.real.iter_mut().zip(&mut out.imaginary).zip(u.as_chunks::<2>().0) {
            *r = values[0];
            *i = values[1];
        }
        let out = refine(a, b, rhs, out, |r| {
            let mut correction: Vec<f64> = r.real.iter().zip(&r.imaginary).flat_map(|(x, y)| [*x, -*y]).collect();
            f.solve(&mut correction)?;
            for ((x, y), pair) in r.real.iter_mut().zip(&mut r.imaginary).zip(correction.as_chunks::<2>().0) {
                *x = pair[0];
                *y = pair[1];
            }
            Ok(())
        })?;
        (out, f.entries(), matrix.n)
    };
    displacement.check(a.n)?;
    let residual = residual_unchecked(a, b, rhs, &displacement)?;
    // The load-relative certificate is intentionally independent of backward error:
    // a tiny backward error alone does not certify equilibrium near resonance.
    if residual.relative > 1e-8 {
        return Err(Error::NotConverged(format!(
            "harmonic original-equation residual {} exceeds 1e-8 (backward error {})",
            residual.relative, residual.backward
        )));
    }
    Ok(Solution { displacement, residual, factor_entries, factor_equations })
}

/// Certificate for any candidate (including a truncated modal solution), without
/// treating a small projected residual as a small full-equation residual.
pub fn residual(a: &Csr, b: &Csr, rhs: &Vector, u: &Vector, limits: &Limits) -> Result<Residual> {
    pair(a, b, limits)?;
    rhs.check(a.n)?;
    u.check(a.n)?;
    residual_unchecked(a, b, rhs, u)
}

fn imbalance(a: &Csr, b: &Csr, rhs: &Vector, u: &Vector) -> Result<Vector> {
    let mut rows = vec![[0.0; 2]; a.n];
    fieldcraft_parallel::chunks_mut(&mut rows, 4096, |block, rows| {
        for (offset, out) in rows.iter_mut().enumerate() {
            let r = block * 4096 + offset;
            let (cols, ar) = a.row(r);
            let (_, br) = b.row(r);
            let mut sum = [rhs.real[r], rhs.imaginary[r]];
            let mut error = [0.0; 2];
            for ((&c, &x), &y) in cols.iter().zip(ar).zip(br) {
                for values in [[-x * u.real[c], -x * u.imaginary[c]], [y * u.imaginary[c], -y * u.real[c]]] {
                    for i in 0..2 {
                        let next = sum[i] + values[i];
                        error[i] += if sum[i].abs() >= values[i].abs() { (sum[i] - next) + values[i] } else { (values[i] - next) + sum[i] };
                        sum[i] = next;
                    }
                }
            }
            *out = [sum[0] + error[0], sum[1] + error[1]];
        }
    });
    if rows.iter().flatten().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("overflow evaluating the complex equation imbalance".into()));
    }
    Ok(Vector { real: rows.iter().map(|v| v[0]).collect(), imaginary: rows.iter().map(|v| v[1]).collect() })
}

fn refine(a: &Csr, b: &Csr, rhs: &Vector, mut u: Vector, solve: impl Fn(&mut Vector) -> Result<()>) -> Result<Vector> {
    let norm = |v: &Vector| v.real.iter().zip(&v.imaginary).map(|(r, i)| libm::hypot(*r, *i)).fold(0.0_f64, f64::max);
    let mut residual = imbalance(a, b, rhs, &u)?;
    let mut size = norm(&residual);
    for _ in 0..3 {
        if size == 0.0 {
            break;
        }
        fieldcraft_model::progress::check()?;
        let mut correction = residual.clone();
        solve(&mut correction)?;
        let mut trial = u.clone();
        for (x, d) in trial.real.iter_mut().chain(&mut trial.imaginary).zip(correction.real.iter().chain(&correction.imaginary)) {
            *x += d;
        }
        trial.check(a.n)?;
        let next = imbalance(a, b, rhs, &trial)?;
        let next_size = norm(&next);
        if next_size >= size {
            break;
        }
        u = trial;
        residual = next;
        size = next_size;
    }
    Ok(u)
}

fn residual_unchecked(a: &Csr, b: &Csr, rhs: &Vector, u: &Vector) -> Result<Residual> {
    let imbalance = imbalance(a, b, rhs, u)?;
    let mut norm_r = 0.0_f64;
    let mut norm_f = 0.0_f64;
    let mut norm_u = 0.0_f64;
    let mut norm_z = 0.0_f64;
    let mut norm = 0.0_f64;
    for r in 0..a.n {
        let residual = libm::hypot(imbalance.real[r], imbalance.imaginary[r]);
        let load = libm::hypot(rhs.real[r], rhs.imaginary[r]);
        let displacement = libm::hypot(u.real[r], u.imaginary[r]);
        let (_, ar) = a.row(r);
        let (_, br) = b.row(r);
        let row: f64 = ar.iter().zip(br).map(|(a, b)| libm::hypot(*a, *b)).sum();
        if !residual.is_finite() || !load.is_finite() || !displacement.is_finite() || !row.is_finite() {
            return Err(Error::Invalid("overflow in harmonic residual evaluation".into()));
        }
        norm_r = norm_r.max(residual);
        norm = libm::hypot(norm, residual);
        norm_f = norm_f.max(load);
        norm_u = norm_u.max(displacement);
        norm_z = norm_z.max(row);
    }
    if norm_r == 0.0 {
        return Ok(Residual { norm: 0.0, relative: 0.0, backward: 0.0 });
    }
    if norm_f == 0.0 || !norm.is_finite() {
        return Err(Error::NotConverged("nonzero harmonic residual for zero force or unrepresentable residual norm".into()));
    }
    // Divide before multiplying, to avoid overflowing ‖Z‖ ‖u‖.
    let scale = norm_f.max(norm_z);
    let backward = (norm_r / scale) / ((norm_z / scale) * norm_u + norm_f / scale);
    let relative = norm_r / norm_f;
    if !relative.is_finite() || !backward.is_finite() {
        return Err(Error::Invalid("overflow in harmonic residual ratios".into()));
    }
    Ok(Residual { norm, relative, backward })
}
