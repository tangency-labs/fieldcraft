//! Steady-state harmonic equilibrium and real-mode superposition.
//!
//! The physical equations and damping convention are in
//! [`fieldcraft_model::harmonic`]. Real modes come from the existing shifted block
//! Lanczos implementation ([`crate::eigen`], Ericsson–Ruhe 1980 and
//! Grimes–Lewis–Simon 1994). The modal expansion follows J. K. Vandiver,
//! MIT 2.003SC (2011), lectures 24–25, and R. R. Craig / A. J. Kurdila,
//! *Fundamentals of Structural Dynamics*, 2nd ed. (2006), modal response analysis.
//! Material viscous and loss matrices are projected **in full**, Φᵀ C Φ and
//! Φᵀ G Φ; no diagonal approximation is made for nonproportional damping.
//! Global Rayleigh damping and modal ratios use the diagonal modal operator.
//! The reduced solve uses [`crate::complex`]. Both projected and full equation
//! residuals are returned: modal truncation is never certified by the former.
//! A spectral shift only extracts free-body modes: it cannot certify a static
//! inverse. Zero-frequency modal response therefore requires an unshifted
//! positive-definite stiffness. An undamped diagonal pole inside the relative
//! pivot resolution is refused, never moved or given artificial damping.

use fieldcraft_model::harmonic::Damping;

use crate::analysis::{BLOCK, TOLERANCE, factor_raw, lanczos_basis};
use crate::cholesky::{Factor, FactorError, PIVOT_TOLERANCE};
use crate::complex::{self, Vector};
use crate::eigen::{self, Operator, Wanted};
use crate::sparse::Csr;
use crate::{Error, Limits, Result};

mod model;
pub use model::{Frame, HarmonicResult, solve_harmonic};

/// Reduced physical matrices in one equation order and on the same CSR pattern.
/// C and G contain material and discrete contributions only; global step damping
/// is supplied at the frequency evaluation boundary.
#[derive(Clone, Copy)]
pub struct Matrices<'a> {
    pub stiffness: &'a Csr,
    pub mass: &'a Csr,
    pub viscous: Option<&'a Csr>,
    pub loss: Option<&'a Csr>,
}

impl Matrices<'_> {
    pub fn validate(&self, limits: &Limits) -> Result<()> {
        complex::validate(self.stiffness, limits)?;
        for a in [Some(self.mass), self.viscous, self.loss].into_iter().flatten() {
            complex::validate(a, limits)?;
            if a.n != self.stiffness.n || a.ptr != self.stiffness.ptr || a.idx != self.stiffness.idx {
                return Err(Error::Invalid("harmonic matrices must share dimensions and a CSR pattern".into()));
            }
        }
        Ok(())
    }

    /// K − ω²M and ω(αM + βK + C) + ηK + G, with ω in radians/second.
    pub fn dynamic(&self, omega: f64, damping: Damping, limits: &Limits) -> Result<(Csr, Csr)> {
        self.validate(limits)?;
        damping.validate()?;
        angular_frequency(omega)?;
        let a = self.stiffness.add_scaled(-omega * omega, self.mass)?;
        let mut b = Csr { val: vec![0.0; a.val.len()], ..a.clone() };
        for (i, value) in b.val.iter_mut().enumerate() {
            let k = self.stiffness.val[i];
            let m = self.mass.val[i];
            let c = self.viscous.and_then(|c| c.val.get(i)).copied().unwrap_or(0.0);
            let g = self.loss.and_then(|g| g.val.get(i)).copied().unwrap_or(0.0);
            *value = omega * (damping.alpha * m + damping.beta * k + c) + damping.loss_factor * k + g;
        }
        if a.val.iter().chain(&b.val).any(|v| !v.is_finite()) {
            return Err(Error::Invalid("overflow assembling the harmonic dynamic stiffness".into()));
        }
        Ok((a, b))
    }

    pub fn direct(&self, omega: f64, damping: Damping, force: &Vector, limits: &Limits) -> Result<complex::Solution> {
        let (a, b) = self.dynamic(omega, damping, limits)?;
        complex::solve(&a, &b, force, limits)
    }
}

fn angular_frequency(omega: f64) -> Result<()> {
    if !omega.is_finite() || omega < 0.0 || !(omega * omega).is_finite() {
        return Err(Error::Invalid("harmonic angular frequency must be finite, nonnegative and have a representable square".into()));
    }
    Ok(())
}

struct ShiftInvert<'a> {
    mass: &'a Csr,
    factor: &'a Factor,
}

impl Operator for ShiftInvert<'_> {
    fn n(&self) -> usize {
        self.mass.n
    }
    fn mass(&self, x: &[f64], y: &mut [f64]) {
        self.mass.mul(x, y);
    }
    fn solve(&self, x: &mut [f64]) -> Result<()> {
        self.factor.solve(x)
    }
}

fn factor_error(e: FactorError) -> Error {
    match e {
        FactorError::NotPositiveDefinite { column } | FactorError::Singular { column } => {
            Error::Singular(format!("harmonic mode extraction: no pivot at equation {column}"))
        }
        FactorError::TooLarge { entries, limit } => Error::factor_limit(entries, limit),
        FactorError::OutOfMemory { entries } => Error::Limit(format!("cannot allocate {entries} factor entries for harmonic modes")),
        FactorError::Internal(message) => Error::Internal(message),
        FactorError::Cancelled => fieldcraft_model::Error::Cancelled.into(),
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn dense(n: usize, values: Vec<f64>) -> Csr {
    Csr { n, ptr: (0..=n).map(|i| i * n).collect(), idx: (0..n).flat_map(|_| 0..n).collect(), val: values }
}

/// Project one symmetric operator in a fixed order; the same dot product is used
/// for each pair and its mirror, so rounding cannot break transpose symmetry.
fn project(matrix: Option<&Csr>, modes: &[Vec<f64>], n: usize) -> Vec<f64> {
    let count = modes.len();
    let mut out = vec![0.0; count * count];
    if let Some(matrix) = matrix {
        let mut y = vec![0.0; n];
        for j in 0..count {
            matrix.mul(&modes[j], &mut y);
            for i in 0..=j {
                let value = dot(&modes[i], &y);
                out[i * count + j] = value;
                out[j * count + i] = value;
            }
        }
    }
    out
}

/// A reusable modal basis and the full projections of material damping. No load
/// or frequency is retained, so one extraction serves an entire sweep.
pub struct Modal<'a> {
    eigenvalues: Vec<f64>,
    shapes: Vec<Vec<f64>>,
    pub factor_entries: usize,
    pub lanczos_vectors: usize,
    viscous: Vec<f64>,
    loss: Vec<f64>,
    n: usize,
    matrices: Matrices<'a>,
    unshifted_stiffness: bool,
}

pub struct ModalSolution {
    pub displacement: Vector,
    pub coordinates: Vector,
    pub projected_residual: complex::Residual,
    pub full_relative_residual: f64,
    pub full_residual_norm: f64,
    /// True only when the projected operator has no off-diagonal damping terms.
    pub diagonal: bool,
}

impl<'a> Modal<'a> {
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }
    pub fn shapes(&self) -> &[Vec<f64>] {
        &self.shapes
    }

    pub fn new(matrices: &Matrices<'a>, count: usize, limits: &Limits) -> Result<Self> {
        matrices.validate(limits)?;
        let n = matrices.stiffness.n;
        if count == 0 || count > fieldcraft_model::steps::MAX_MODES as usize {
            return Err(Error::Invalid("harmonic mode count is outside the model limit".into()));
        }
        let count = count.min(n);
        let entries = count.checked_mul(n).and_then(|x| count.checked_mul(count).and_then(|p| p.checked_mul(2)?.checked_add(x)));
        if entries.is_none_or(|e| e > limits.max_factor_entries) {
            return Err(Error::Limit("harmonic modal basis and projected damping exceed the entry budget".into()));
        }
        let k = matrices.stiffness;
        let m = matrices.mass;
        let (shift, factor) = match factor_raw(k, limits) {
            Ok(f) => (0.0, f),
            Err(FactorError::NotPositiveDefinite { .. }) => {
                let tk: f64 = k.diagonal().iter().sum();
                let tm: f64 = m.diagonal().iter().sum();
                if !(tk > 0.0 && tm > 0.0 && tk.is_finite() && tm.is_finite()) {
                    return Err(Error::Singular("harmonic mode extraction requires positive stiffness and mass traces".into()));
                }
                let shift = -1e-6 * tk / tm; // unit-factors: relative spectral shift, as frequency analysis
                (shift, factor_raw(&k.add_scaled(-shift, m)?, limits).map_err(factor_error)?)
            }
            Err(e) => return Err(factor_error(e)),
        };
        let op = ShiftInvert { mass: m, factor: &factor };
        let wanted = Wanted { count, min: f64::NEG_INFINITY, max: f64::INFINITY };
        let pairs = eigen::lanczos(&op, shift, wanted, BLOCK, lanczos_basis(limits, n), TOLERANCE)?;
        if pairs.values.len() != count {
            return Err(Error::NotConverged(format!("harmonic extraction found {} of {count} modes", pairs.values.len())));
        }
        // Verify against the original pencil, not just the Lanczos recurrence.
        let mut kphi = vec![0.0; n];
        let mut mphi = vec![0.0; n];
        let norm_k = (0..n).map(|i| k.row(i).1.iter().map(|v| v.abs()).sum::<f64>()).fold(0.0_f64, f64::max);
        let norm_m = (0..n).map(|i| m.row(i).1.iter().map(|v| v.abs()).sum::<f64>()).fold(0.0_f64, f64::max);
        for (j, (lambda, phi)) in pairs.values.iter().zip(&pairs.vectors).enumerate() {
            k.mul(phi, &mut kphi);
            m.mul(phi, &mut mphi);
            let scale = (norm_k + lambda.abs() * norm_m) * phi.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            let r = kphi.iter().zip(&mphi).map(|(k, m)| (k - lambda * m).abs()).fold(0.0_f64, f64::max);
            if !lambda.is_finite() || r > 1e-7 * scale || (dot(phi, &mphi) - 1.0).abs() > 1e-7 {
                return Err(Error::NotConverged("harmonic mode failed original-pencil residual or mass normalization".into()));
            }
            for (i, other) in pairs.vectors.iter().enumerate() {
                if i != j && dot(other, &mphi).abs() > 1e-7 {
                    return Err(Error::NotConverged("harmonic modes failed mass orthogonality".into()));
                }
            }
        }
        let viscous = project(matrices.viscous, &pairs.vectors, n);
        let loss = project(matrices.loss, &pairs.vectors, n);
        Ok(Self {
            eigenvalues: pairs.values,
            shapes: pairs.vectors,
            factor_entries: factor.entries(),
            lanczos_vectors: pairs.basis,
            viscous,
            loss,
            n,
            matrices: *matrices,
            unshifted_stiffness: shift == 0.0,
        })
    }

    pub fn solve(&self, omega: f64, damping: Damping, ratios: &[f64], force: &Vector, limits: &Limits) -> Result<ModalSolution> {
        let matrices = &self.matrices;
        angular_frequency(omega)?;
        damping.validate()?;
        let count = self.eigenvalues.len();
        if matrices.stiffness.n != self.n
            || force.real.len() != self.n
            || force.imaginary.len() != self.n
            || force.real.iter().chain(&force.imaginary).any(|v| !v.is_finite())
        {
            return Err(Error::Invalid("harmonic force or matrices do not match the modal basis".into()));
        }
        if ![0, 1, count].contains(&ratios.len()) || ratios.iter().any(|r| !r.is_finite() || *r < 0.0) {
            return Err(Error::Invalid("modal damping needs 0, 1 or one finite nonnegative ratio per extracted mode".into()));
        }
        if omega == 0.0 && !self.unshifted_stiffness {
            return Err(Error::NotSupported(
                "zero-frequency modal response requires a positive-definite unshifted stiffness; use the direct method to check the original system"
                    .into(),
            ));
        }
        let mut a = vec![0.0; count * count];
        let mut b: Vec<f64> = self.viscous.iter().zip(&self.loss).map(|(c, g)| omega * c + g).collect();
        let mut modal_c = Vec::with_capacity(count);
        let mut rhs = Vector::zeros(count);
        for (j, (lambda, phi)) in self.eigenvalues.iter().zip(&self.shapes).enumerate() {
            let ratio = ratios.get(if ratios.len() == 1 { 0 } else { j }).copied().unwrap_or(0.0);
            let c = 2.0 * ratio * lambda.max(0.0).sqrt();
            modal_c.push(c);
            a[j * count + j] = lambda - omega * omega;
            b[j * count + j] += omega * (damping.alpha + damping.beta * lambda + c) + damping.loss_factor * lambda;
            rhs.real[j] = dot(phi, &force.real);
            rhs.imaginary[j] = dot(phi, &force.imaginary);
        }
        let diagonal = b.iter().enumerate().all(|(i, v)| i / count == i % count || *v == 0.0);
        let a = dense(count, a);
        let b = dense(count, b);
        let q = if diagonal {
            let mut q = Vector::zeros(count);
            for j in 0..count {
                let ar = a.val[j * count + j];
                let bi = b.val[j * count + j];
                let scale = ar.abs().max(bi.abs());
                if scale == 0.0 {
                    return Err(Error::Singular(format!("undamped harmonic resonance in mode {}", j + 1)));
                }
                let pole_scale = self.eigenvalues[j].abs().max(omega * omega);
                if bi == 0.0 && ar.abs() <= PIVOT_TOLERANCE * pole_scale {
                    return Err(Error::NotConverged(format!(
                        "undamped forcing frequency is not numerically separated from mode {} within the relative pivot tolerance",
                        j + 1
                    )));
                }
                let (ar, bi) = (ar / scale, bi / scale);
                let denom = ar * ar + bi * bi;
                let (fr, fi) = (rhs.real[j] / scale, rhs.imaginary[j] / scale);
                q.real[j] = (ar * fr + bi * fi) / denom;
                q.imaginary[j] = (ar * fi - bi * fr) / denom;
            }
            q
        } else {
            complex::solve(&a, &b, &rhs, limits)?.displacement
        };
        let projected_residual = complex::residual(&a, &b, &rhs, &q, limits)?;
        let mut displacement = Vector::zeros(self.n);
        // Additional modal damping force is iω M Φ diag(2ζω_j) q.
        let mut damping_shape = Vector::zeros(self.n);
        for (j, phi) in self.shapes.iter().enumerate() {
            for (i, p) in phi.iter().enumerate() {
                displacement.real[i] += p * q.real[j];
                displacement.imaginary[i] += p * q.imaginary[j];
                damping_shape.real[i] -= omega * p * modal_c[j] * q.imaginary[j];
                damping_shape.imaginary[i] += omega * p * modal_c[j] * q.real[j];
            }
        }
        let (physical_a, physical_b) = matrices.dynamic(omega, damping, limits)?;
        let (full_relative_residual, full_residual_norm) =
            full_residual(&physical_a, &physical_b, matrices.mass, &displacement, force, &damping_shape)?;
        Ok(ModalSolution { displacement, coordinates: q, projected_residual, full_relative_residual, full_residual_norm, diagonal })
    }
}

fn full_residual(a: &Csr, b: &Csr, mass: &Csr, u: &Vector, force: &Vector, damping_shape: &Vector) -> Result<(f64, f64)> {
    let mut real = vec![0.0; a.n];
    let mut imaginary = vec![0.0; a.n];
    let mut scratch = vec![0.0; a.n];
    a.mul(&u.real, &mut real);
    b.mul(&u.imaginary, &mut scratch);
    for (r, b) in real.iter_mut().zip(&scratch) {
        *r -= b;
    }
    a.mul(&u.imaginary, &mut imaginary);
    b.mul(&u.real, &mut scratch);
    for (i, b) in imaginary.iter_mut().zip(&scratch) {
        *i += b;
    }
    mass.mul(&damping_shape.real, &mut scratch);
    for (r, b) in real.iter_mut().zip(&scratch) {
        *r += b;
    }
    mass.mul(&damping_shape.imaginary, &mut scratch);
    for (i, b) in imaginary.iter_mut().zip(&scratch) {
        *i += b;
    }
    let mut residual = 0.0_f64;
    let mut norm = 0.0_f64;
    let mut residual_norm = 0.0_f64;
    for i in 0..a.n {
        let r = (real[i] - force.real[i]).hypot(imaginary[i] - force.imaginary[i]);
        let f = force.real[i].hypot(force.imaginary[i]);
        if !r.is_finite() || !f.is_finite() {
            return Err(Error::Invalid("non-finite harmonic modal residual".into()));
        }
        residual = residual.max(r);
        residual_norm = residual_norm.hypot(r);
        norm = norm.max(f);
    }
    if residual == 0.0 {
        Ok((0.0, 0.0))
    } else if norm == 0.0 || !residual_norm.is_finite() || !(residual / norm).is_finite() {
        Err(Error::NotConverged("nonzero harmonic residual for zero force or unrepresentable residual ratio".into()))
    } else {
        Ok((residual / norm, residual_norm))
    }
}
