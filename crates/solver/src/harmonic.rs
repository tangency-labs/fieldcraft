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
mod residual;
pub use model::{
    ElementStressRecovery, Frame, HarmonicResult, MatrixExport, StressLocation, export_matrices, export_matrices_with_mass, solve_harmonic,
    solve_harmonic_with_mass, solve_harmonic_with_mass_cached,
};
pub(crate) use model::{damping_operators, damping_operators_with_mass, load_vector, stresses};

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
        self.dynamic_validated(omega, damping, limits)
    }

    // Only called after full validation while these immutable matrix borrows remain live.
    // Per-call limits and frequency-dependent overflow still need checking.
    fn dynamic_validated(&self, omega: f64, damping: Damping, limits: &Limits) -> Result<(Csr, Csr)> {
        if self.stiffness.n == 0 || self.stiffness.n > limits.max_dofs {
            return Err(Error::Limit(format!("complex matrix has {} equations; allowed 1..={}", self.stiffness.n, limits.max_dofs)));
        }
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
        FactorError::MemoryLimit { entries, estimated_bytes, limit_bytes } => Error::Limit(format!(
            "factorization needs {entries} entries; estimated factor storage including indices, matrices and workspace is {estimated_bytes} bytes, exceeding factor budget {limit_bytes} bytes"
        )),
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

/// A reusable modal basis and full damping projections. Optional residual vectors
/// are tailored to a load, but no frequency is retained: one basis serves a sweep.
pub struct Modal<'a> {
    data: std::sync::Arc<ModalData>,
    pub factor_entries: usize,
    pub lanczos_vectors: usize,
    n: usize,
    matrices: Matrices<'a>,
    unshifted_stiffness: bool,
}

#[derive(Debug)]
struct ModalData {
    eigenvalues: Vec<f64>,
    shapes: Vec<Vec<f64>>,
    viscous: Vec<f64>,
    loss: Vec<f64>,
    stiffness: Vec<f64>,
    mass: Vec<f64>,
    factor_entries: usize,
    lanczos_vectors: usize,
    unshifted_stiffness: bool,
}

mod cache;
pub use cache::ModalCache;

pub struct ModalSolution {
    pub displacement: Vector,
    pub coordinates: Vector,
    pub projected_residual: complex::Residual,
    pub full_relative_residual: f64,
    pub full_residual_norm: f64,
    /// True only when the projected dynamic operator is diagonal.
    pub diagonal: bool,
}

struct ModalCoordinates {
    q: Vector,
    coordinates: Vector,
    modal_c: Vec<f64>,
    projected_residual: complex::Residual,
    diagonal: bool,
    omega: f64,
}

impl<'a> Modal<'a> {
    pub fn eigenvalues(&self) -> &[f64] {
        &self.data.eigenvalues
    }
    pub fn shapes(&self) -> &[Vec<f64>] {
        &self.data.shapes[..self.data.eigenvalues.len()]
    }
    /// Complete fixed projection basis, including load-dependent Ritz corrections.
    pub fn basis_vectors(&self) -> &[Vec<f64>] {
        &self.data.shapes
    }
    pub fn residual_vectors(&self) -> usize {
        self.data.shapes.len() - self.data.eigenvalues.len()
    }

    pub fn new(matrices: &Matrices<'a>, count: usize, limits: &Limits) -> Result<Self> {
        Self::build(matrices, count, None, limits)
    }

    /// Add at most `cycles` static/Ritz vectors per load quadrature. The factor
    /// used by Lanczos is reused before it is released. A singular unshifted K
    /// is refused; the original modal method remains available for free bodies.
    pub fn with_residual_vectors(matrices: &Matrices<'a>, count: usize, force: &Vector, cycles: usize, limits: &Limits) -> Result<Self> {
        if cycles > 16 {
            return Err(Error::Invalid("harmonic residual-vector cycles must be 0 to 16".into()));
        }
        Self::build(matrices, count, Some((force, cycles)), limits)
    }

    fn build(matrices: &Matrices<'a>, count: usize, residual: Option<(&Vector, usize)>, limits: &Limits) -> Result<Self> {
        matrices.validate(limits)?;
        let n = matrices.stiffness.n;
        if count == 0 || count > fieldcraft_model::steps::MAX_MODES as usize {
            return Err(Error::Invalid("harmonic mode count is outside the model limit".into()));
        }
        let count = count.min(n);
        let budget_count = (count + 2 * residual.map_or(0, |(_, cycles)| cycles)).min(n);
        let entries = budget_count.checked_mul(n).and_then(|x| budget_count.checked_mul(budget_count).and_then(|p| p.checked_mul(4)?.checked_add(x)));
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
        if residual.is_some_and(|(_, cycles)| cycles > 0) && shift != 0.0 {
            return Err(Error::NotSupported("harmonic residual vectors require positive-definite unshifted stiffness".into()));
        }
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
        let mut shapes = pairs.vectors;
        if let Some((force, cycles)) = residual
            && cycles > 0
        {
            residual::augment(matrices, &factor, &mut shapes, force, cycles)?;
        }
        let basis_count = shapes.len();
        let (stiffness, mass) = if basis_count > count {
            (project(Some(k), &shapes, n), project(Some(m), &shapes, n))
        } else {
            let mut stiffness = vec![0.0; count * count];
            let mut mass = vec![0.0; count * count];
            for (j, lambda) in pairs.values.iter().enumerate() {
                stiffness[j * count + j] = *lambda;
                mass[j * count + j] = 1.0;
            }
            (stiffness, mass)
        };
        let viscous = project(matrices.viscous, &shapes, n);
        let loss = project(matrices.loss, &shapes, n);
        let data = std::sync::Arc::new(ModalData {
            eigenvalues: pairs.values,
            shapes,
            factor_entries: factor.entries(),
            lanczos_vectors: pairs.basis,
            viscous,
            loss,
            stiffness,
            mass,
            unshifted_stiffness: shift == 0.0,
        });
        Ok(Self::from_data(matrices, data))
    }

    fn from_data(matrices: &Matrices<'a>, data: std::sync::Arc<ModalData>) -> Self {
        Self {
            factor_entries: data.factor_entries,
            lanczos_vectors: data.lanczos_vectors,
            unshifted_stiffness: data.unshifted_stiffness,
            n: matrices.stiffness.n,
            matrices: *matrices,
            data,
        }
    }

    // Construction validates all matrices, whose shared borrows prohibit mutation.
    fn dynamic(&self, omega: f64, damping: Damping, limits: &Limits) -> Result<(Csr, Csr)> {
        self.matrices.dynamic_validated(omega, damping, limits)
    }

    pub fn solve(&self, omega: f64, damping: Damping, ratios: &[f64], force: &Vector, limits: &Limits) -> Result<ModalSolution> {
        self.solve_with_operator(omega, damping, ratios, force, limits).map(|(solution, _)| solution)
    }

    // The model-level power audit uses the same physical loss operator as the
    // full residual, avoiding a second construction without changing arithmetic.
    fn prepare_coordinates(
        &self,
        omega: f64,
        damping: Damping,
        ratios: &[f64],
        force: &Vector,
        limits: &Limits,
        projected_force: Option<&Vector>,
    ) -> Result<ModalCoordinates> {
        let matrices = &self.matrices;
        angular_frequency(omega)?;
        damping.validate()?;
        let count = self.data.shapes.len();
        let mode_count = self.data.eigenvalues.len();
        if matrices.stiffness.n != self.n
            || force.real.len() != self.n
            || force.imaginary.len() != self.n
            || force.real.iter().chain(&force.imaginary).any(|v| !v.is_finite())
        {
            return Err(Error::Invalid("harmonic force or matrices do not match the modal basis".into()));
        }
        if ![0, 1, mode_count].contains(&ratios.len()) || ratios.iter().any(|r| !r.is_finite() || *r < 0.0) {
            return Err(Error::Invalid("modal damping needs 0, 1 or one finite nonnegative ratio per extracted mode".into()));
        }
        if omega == 0.0 && !self.unshifted_stiffness {
            return Err(Error::NotSupported(
                "zero-frequency modal response requires a positive-definite unshifted stiffness; use the direct method to check the original system"
                    .into(),
            ));
        }
        let a: Vec<f64> = self.data.stiffness.iter().zip(&self.data.mass).map(|(k, m)| k - omega * omega * m).collect();
        let mut b: Vec<f64> = self.data.viscous.iter().zip(&self.data.loss).map(|(c, g)| omega * c + g).collect();
        if count > mode_count {
            for ((value, k), m) in b.iter_mut().zip(&self.data.stiffness).zip(&self.data.mass) {
                *value += omega * (damping.alpha * m + damping.beta * k) + damping.loss_factor * k;
            }
        }
        let mut modal_c = Vec::with_capacity(mode_count);
        let mut rhs = Vector::zeros(count);
        for (j, lambda) in self.data.eigenvalues.iter().enumerate() {
            let ratio = ratios.get(if ratios.len() == 1 { 0 } else { j }).copied().unwrap_or(0.0);
            let c = 2.0 * ratio * lambda.max(0.0).sqrt();
            modal_c.push(c);
            if count == mode_count {
                b[j * count + j] += omega * (damping.alpha + damping.beta * lambda + c) + damping.loss_factor * lambda;
                continue;
            }
            // Vᵀ M Φ diag(2ζω_j) Φᵀ M V, including roundoff in M-orthogonality.
            // Assemble one triangle and mirror it. Evaluating the two products
            // separately changes multiplication order and can break the exact
            // transpose symmetry required by the complex factorization.
            for row in 0..count {
                for col in 0..=row {
                    b[row * count + col] += omega * c * self.data.mass[row * count + j] * self.data.mass[j * count + col];
                    b[col * count + row] = b[row * count + col];
                }
            }
        }
        if let Some(projected) = projected_force {
            rhs = projected.clone();
        } else {
            for (j, phi) in self.data.shapes.iter().enumerate() {
                rhs.real[j] = dot(phi, &force.real);
                rhs.imaginary[j] = dot(phi, &force.imaginary);
            }
        }
        let diagonal = a.iter().zip(&b).enumerate().all(|(i, (a, b))| i / count == i % count || (*a == 0.0 && *b == 0.0));
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
                let pole_scale = self.data.stiffness[j * count + j].abs().max(omega * omega * self.data.mass[j * count + j].abs());
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
        let mut coordinates = Vector::zeros(mode_count);
        for j in 0..mode_count {
            let row = &self.data.mass[j * count..(j + 1) * count];
            coordinates.real[j] = dot(row, &q.real);
            coordinates.imaginary[j] = dot(row, &q.imaginary);
        }
        Ok(ModalCoordinates { q, coordinates, modal_c, projected_residual, diagonal, omega })
    }
    fn finish_response(
        &self,
        damping: Damping,
        force: &Vector,
        limits: &Limits,
        prepared: &ModalCoordinates,
        displacement: Vector,
        damping_shape: Vector,
    ) -> Result<(ModalSolution, Csr)> {
        let _timing = crate::timing::Span::new("harmonic_full_residual");
        let (physical_a, physical_b) = self.dynamic(prepared.omega, damping, limits)?;
        let (full_relative_residual, full_residual_norm) =
            full_residual(&physical_a, &physical_b, self.matrices.mass, &displacement, force, &damping_shape)?;
        Ok((
            ModalSolution {
                displacement,
                coordinates: prepared.coordinates.clone(),
                projected_residual: prepared.projected_residual,
                full_relative_residual,
                full_residual_norm,
                diagonal: prepared.diagonal,
            },
            physical_b,
        ))
    }

    fn solve_with_operator(&self, omega: f64, damping: Damping, ratios: &[f64], force: &Vector, limits: &Limits) -> Result<(ModalSolution, Csr)> {
        let prepared = self.prepare_coordinates(omega, damping, ratios, force, limits, None)?;
        let mut displacement = Vector::zeros(self.n);
        let mut damping_shape = Vector::zeros(self.n);
        let _timing = crate::timing::Span::new("harmonic_displacement_recovery");
        for (j, phi) in self.data.shapes.iter().enumerate() {
            for (i, p) in phi.iter().enumerate() {
                displacement.real[i] += p * prepared.q.real[j];
                displacement.imaginary[i] += p * prepared.q.imaginary[j];
                if j < self.data.eigenvalues.len() {
                    damping_shape.real[i] -= omega * p * prepared.modal_c[j] * prepared.coordinates.imaginary[j];
                    damping_shape.imaginary[i] += omega * p * prepared.modal_c[j] * prepared.coordinates.real[j];
                }
            }
        }
        drop(_timing);
        self.finish_response(damping, force, limits, &prepared, displacement, damping_shape)
    }

    // Phi Q in fixed DOF tiles and bounded frequency panels. Each output belongs to
    // one task; every sum follows ascending basis index, independent of scheduling.
    fn recover_batch(&self, coordinates: &[ModalCoordinates]) -> std::collections::VecDeque<(Vector, Vector)> {
        let _timing = crate::timing::Span::new("harmonic_displacement_recovery");
        let width = coordinates.len();
        let mut tile = vec![[0.0; 4]; self.n * width];
        fieldcraft_parallel::chunks_mut(&mut tile, 256 * width, |block, output| {
            let first = block * 256;
            for (j, phi) in self.data.shapes.iter().enumerate() {
                for (i, row) in output.chunks_mut(width).enumerate() {
                    let p = phi[first + i];
                    for (value, c) in row.iter_mut().zip(coordinates) {
                        value[0] += p * c.q.real[j];
                        value[1] += p * c.q.imaginary[j];
                        if j < self.data.eigenvalues.len() {
                            value[2] -= c.omega * p * c.modal_c[j] * c.coordinates.imaginary[j];
                            value[3] += c.omega * p * c.modal_c[j] * c.coordinates.real[j];
                        }
                    }
                }
            }
        });
        (0..width)
            .map(|f| {
                let mut u = Vector::zeros(self.n);
                let mut d = Vector::zeros(self.n);
                for (i, row) in tile.chunks_exact(width).enumerate() {
                    [u.real[i], u.imaginary[i], d.real[i], d.imaginary[i]] = row[f];
                }
                (u, d)
            })
            .collect()
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
        let r = libm::hypot(real[i] - force.real[i], imaginary[i] - force.imaginary[i]);
        let f = libm::hypot(force.real[i], force.imaginary[i]);
        if !r.is_finite() || !f.is_finite() {
            return Err(Error::Invalid("non-finite harmonic modal residual".into()));
        }
        residual = residual.max(r);
        residual_norm = libm::hypot(residual_norm, r);
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

#[cfg(test)]
mod batch_tests {
    use super::*;
    #[test]
    fn blocked_recovery_matches_scalar_and_closed_form_at_every_thread_count() -> Result<()> {
        let k = dense(2, vec![4.0, 0.0, 0.0, 9.0]);
        let m = dense(2, vec![1.0, 0.0, 0.0, 2.0]);
        let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
        let force = Vector { real: vec![1.0, 2.0], imaginary: vec![0.4, -0.3] };
        let limits = Limits::default();
        let damping = Damping { alpha: 0.1, beta: 0.02, loss_factor: 0.01 };
        let modal = Modal::new(&matrices, 2, &limits)?;
        let prepared = [0.0, 1.99, 2.0, 2.01, 3.0]
            .into_iter()
            .map(|w| modal.prepare_coordinates(w, damping, &[], &force, &limits, None))
            .collect::<Result<Vec<_>>>()?;
        let run = || -> Result<Vec<u64>> {
            let mut bits = Vec::new();
            for (p, (u, d)) in prepared.iter().zip(modal.recover_batch(&prepared)) {
                let (out, _) = modal.finish_response(damping, &force, &limits, p, u, d)?;
                let scalar = modal.solve(p.omega, damping, &[], &force, &limits)?;
                assert_eq!(out.displacement, scalar.displacement);
                assert!(out.full_relative_residual < 1e-12);
                for i in 0..2 {
                    let stiffness = k.val[i * 2 + i];
                    let mass = m.val[i * 2 + i];
                    let a = stiffness - p.omega * p.omega * mass;
                    let b = p.omega * (damping.alpha * mass + damping.beta * stiffness) + damping.loss_factor * stiffness;
                    let real = (a * force.real[i] + b * force.imaginary[i]) / (a * a + b * b);
                    let imag = (a * force.imaginary[i] - b * force.real[i]) / (a * a + b * b);
                    assert!((out.displacement.real[i] - real).abs() <= 1e-12 * real.abs().max(1.0));
                    assert!((out.displacement.imaginary[i] - imag).abs() <= 1e-12 * imag.abs().max(1.0));
                }
                bits.extend(out.displacement.real.iter().chain(&out.displacement.imaginary).map(|x| x.to_bits()));
            }
            Ok(bits)
        };
        let expected = fieldcraft_parallel::with_threads(1, run)?;
        for threads in [2, 7] {
            assert_eq!(expected, fieldcraft_parallel::with_threads(threads, run)?);
        }
        Ok(())
    }
}
