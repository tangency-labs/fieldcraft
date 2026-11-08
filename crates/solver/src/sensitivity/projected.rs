//! Fixed-basis Galerkin harmonic derivatives; projection is done once per sweep.
//! Z_r q=f_r; q'=Z_r^-1(f_r'-Z_r' q). Haug, Choi & Komkov (1986).
use super::{HarmonicDerivative, HarmonicSensitivity};
use crate::{
    Error, Limits, Result,
    complex::{self, Vector},
    harmonic::{Matrices, Modal},
    sparse::Csr,
};
use fieldcraft_model::{
    design::{DampingParameter, Parameter},
    harmonic::Harmonic,
};

pub(super) struct Operators {
    pub k: Csr,
    pub m: Csr,
    pub c: Option<Csr>,
    pub g: Option<Csr>,
    pub f: Vector,
}
fn project(a: Option<&Csr>, basis: &[Vec<f64>]) -> Vec<f64> {
    let n = basis.len();
    let mut out = vec![0.0; n * n];
    if let Some(a) = a {
        let mut product = vec![0.0; a.n];
        for (j, phi) in basis.iter().enumerate() {
            a.mul(phi, &mut product);
            for (i, psi) in basis.iter().enumerate().take(j + 1) {
                let value = psi.iter().zip(&product).map(|(x, y)| x * y).sum();
                out[i * n + j] = value;
                out[j * n + i] = value;
            }
        }
    }
    out
}
fn vector(x: &[f64], basis: &[Vec<f64>]) -> Vec<f64> {
    basis.iter().map(|phi| phi.iter().zip(x).map(|(p, x)| p * x).sum()).collect()
}
fn dense(n: usize, val: Vec<f64>) -> Csr {
    Csr { n, ptr: (0..=n).map(|i| i * n).collect(), idx: (0..n).flat_map(|_| 0..n).collect(), val }
}
pub(super) struct Sweep<'a, 'b> {
    pub modal: &'a Modal<'b>,
    pub matrices: Matrices<'a>,
    pub derivative: Option<Operators>,
    pub force: Vector,
    pub weights: &'a Vector,
    pub harmonic: &'a Harmonic,
    pub parameter: &'a Parameter,
    pub value: f64,
    pub increment: f64,
    pub warnings: Vec<String>,
}
pub(super) fn solve(s: Sweep<'_, '_>, limits: &Limits) -> Result<HarmonicSensitivity> {
    let basis = s.modal.basis_vectors();
    let n = basis.len();
    if n == 0 {
        return Err(Error::Invalid("empty modal derivative basis".into()));
    }
    let k = project(Some(s.matrices.stiffness), basis);
    let m = project(Some(s.matrices.mass), basis);
    let c = project(s.matrices.viscous, basis);
    let g = project(s.matrices.loss, basis);
    let f = Vector { real: vector(&s.force.real, basis), imaginary: vector(&s.force.imaginary, basis) };
    let weights = Vector { real: vector(&s.weights.real, basis), imaginary: vector(&s.weights.imaginary, basis) };
    let (dk, dm, dc, dg, df) = if let Some(d) = &s.derivative {
        (
            project(Some(&d.k), basis),
            project(Some(&d.m), basis),
            project(d.c.as_ref(), basis),
            project(d.g.as_ref(), basis),
            Vector { real: vector(&d.f.real, basis), imaginary: vector(&d.f.imaginary, basis) },
        )
    } else {
        (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n], Vector::zeros(n))
    };
    let ratios =
        if let fieldcraft_model::harmonic::Method::Modal { damping_ratios, .. } = &s.harmonic.method { damping_ratios.as_slice() } else { &[] };
    let damping = s.harmonic.damping;
    let mut frames = Vec::new();
    for frequency in s.harmonic.frequencies.points()? {
        fieldcraft_model::progress::check()?;
        let omega = std::f64::consts::TAU * frequency;
        let a = dense(n, k.iter().zip(&m).map(|(k, m)| k - omega * omega * m).collect());
        let mut b: Vec<_> = k
            .iter()
            .zip(&m)
            .zip(c.iter().zip(&g))
            .map(|((k, m), (c, g))| omega * (c + damping.alpha * m + damping.beta * k) + g + damping.loss_factor * k)
            .collect();
        let da = dense(n, dk.iter().zip(&dm).map(|(k, m)| k - omega * omega * m).collect());
        let mut db: Vec<_> = dk
            .iter()
            .zip(&dm)
            .zip(dc.iter().zip(&dg))
            .map(|((k, m), (c, g))| omega * (c + damping.alpha * m + damping.beta * k) + g + damping.loss_factor * k)
            .collect();
        for (j, lambda) in s.modal.eigenvalues().iter().enumerate() {
            let ratio = ratios.get(if ratios.len() == 1 { 0 } else { j }).copied().unwrap_or(0.0);
            let derivative = matches!(s.parameter,Parameter::Damping{component:DampingParameter::ModalRatio{mode},..} if *mode==j);
            for row in 0..n {
                for col in 0..n {
                    let term = 2.0 * omega * lambda.max(0.0).sqrt() * m[row * n + j] * m[j * n + col];
                    b[row * n + col] += ratio * term;
                    if derivative {
                        db[row * n + col] += term;
                    }
                }
            }
        }
        if let Parameter::Damping { component, .. } = s.parameter {
            match component {
                DampingParameter::Alpha => {
                    for (b, m) in db.iter_mut().zip(&m) {
                        *b = omega * m;
                    }
                }
                DampingParameter::Beta => {
                    for (b, k) in db.iter_mut().zip(&k) {
                        *b = omega * k;
                    }
                }
                DampingParameter::LossFactor => db.clone_from(&k),
                DampingParameter::ModalRatio { mode } => {
                    if *mode >= s.modal.eigenvalues().len() {
                        return Err(Error::Invalid("modal ratio index exceeds extracted modes".into()));
                    }
                }
            }
        }
        // Preserve transpose symmetry bit-for-bit also for modal damping outer products.
        for row in 0..n {
            for col in 0..row {
                b[row * n + col] = b[col * n + row];
                db[row * n + col] = db[col * n + row];
            }
        }
        let b = dense(n, b);
        let db = dense(n, db);
        let primal = complex::solve(&a, &b, &f, limits)?;
        let tangent = super::tangent_load(&da, &db, &primal.displacement, &df, limits)?;
        let derivative = complex::solve(&a, &b, &tangent, limits)?;
        let contract = |u: &Vector| -> [f64; 2] {
            [weights.real.iter().zip(&u.real).map(|(w, u)| w * u).sum(), weights.real.iter().zip(&u.imaginary).map(|(w, u)| w * u).sum()]
        };
        let value = contract(&primal.displacement);
        let direct = contract(&derivative.displacement);
        let adjoint = super::adjoint(&a, &b, &weights, &tangent, limits)?;
        let [amplitude_derivative, phase_derivative] = super::amplitude_phase(value, direct)?;
        frames.push(HarmonicDerivative {
            frequency,
            displacement: value,
            direct,
            adjoint,
            amplitude: libm::hypot(value[0], value[1]),
            phase: libm::atan2(value[1], value[0]),
            amplitude_derivative,
            phase_derivative,
            primal_relative_residual: primal.residual.relative,
            tangent_relative_residual: derivative.residual.relative,
        });
    }
    Ok(HarmonicSensitivity { parameter: s.parameter.clone(), parameter_value: s.value, increment: s.increment, frames, warnings: s.warnings })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::harmonic::{Damping, Frequencies, Method};
    #[test]
    fn full_and_ritz_completed_bases_match_two_coordinate_direct_derivatives() -> Result<()> {
        let limits = Limits::default();
        let k = dense(2, vec![4.0, -1.0, -1.0, 3.0]);
        let m = dense(2, vec![2.0, 0.1, 0.1, 1.0]);
        let dk = dense(2, vec![1.0, 0.0, 0.0, 0.0]);
        let zero = dense(2, vec![0.0; 4]);
        let force = Vector { real: vec![1.0, 0.0], imaginary: vec![0.0, 0.2] };
        let weights = Vector { real: vec![0.0, 1.0], imaginary: vec![0.0; 2] };
        let parameter = Parameter::SpringStiffness { name: "test".into() };
        for (count, cycles) in [(2, 0), (1, 1)] {
            let matrices = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None };
            let modal = Modal::with_residual_vectors(&matrices, count, &force, cycles, &limits)?;
            assert_eq!(modal.basis_vectors().len(), 2);
            let h = Harmonic {
                frequencies: Frequencies::List(vec![0.0, 0.15, 0.3]),
                method: Method::Modal { modes: count as u32, damping_ratios: vec![], residual_vectors: cycles as u32 },
                damping: Damping { alpha: 0.1, beta: 0.02, loss_factor: 0.01 },
                load_factors: Default::default(),
            };
            let r = solve(
                Sweep {
                    modal: &modal,
                    matrices,
                    derivative: Some(Operators { k: dk.clone(), m: zero.clone(), c: None, g: None, f: Vector::zeros(2) }),
                    force: force.clone(),
                    weights: &weights,
                    harmonic: &h,
                    parameter: &parameter,
                    value: 4.0,
                    increment: 1e-5,
                    warnings: vec![],
                },
                &limits,
            )?;
            for frame in r.frames {
                let omega = std::f64::consts::TAU * frame.frequency;
                let (a, b) = Matrices { stiffness: &k, mass: &m, viscous: None, loss: None }.dynamic(omega, h.damping, &limits)?;
                let primal = complex::solve(&a, &b, &force, &limits)?;
                let mut db = dk.clone();
                for v in &mut db.val {
                    *v *= omega * h.damping.beta + h.damping.loss_factor;
                }
                let tangent = crate::sensitivity::tangent_load(&dk, &db, &primal.displacement, &Vector::zeros(2), &limits)?;
                let derivative = complex::solve(&a, &b, &tangent, &limits)?;
                for (actual, expected) in
                    [(frame.direct[0], derivative.displacement.real[1]), (frame.direct[1], derivative.displacement.imaginary[1])]
                {
                    assert!((actual - expected).abs() < 1e-8 * expected.abs().max(1.0), "{actual} vs {expected}");
                }
            }
        }
        Ok(())
    }
}
