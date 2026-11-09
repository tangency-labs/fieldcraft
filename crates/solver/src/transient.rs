//! Linear fixed-step integration of reduced SI equations M a + C v + K u = f.
//!
//! Newmark (1959), doi:10.1061/JMCEA3.0000098; Hilber, Hughes and Taylor
//! (1977), doi:10.1002/eqe.4290050306. HHT uses alpha in [-1/3,0]:
//!
//! ```text
//! M a_next + (1+alpha)(C v_next + K u_next - f_next)
//! - alpha(C v + K u - f) = 0.
//! ```
//!
//! Initial acceleration uses physical equilibrium.
//! The same reduced-matrix integrator accepts projected modal operators, including
//! non-diagonal damping. It stores only the current state; callers budget outputs.

mod model;
pub use model::{Frame, Response, solve};

use crate::cholesky::{Factor, FactorError, UpperColumns};
use crate::sparse::Csr;
use crate::{Error, Limits, Result};

#[derive(Debug, Clone, Copy)]
pub enum Integration {
    AverageAcceleration,
    /// Conditionally stable; caller must check dt * highest omega <= sqrt(12).
    LinearAcceleration,
    Hht {
        alpha: f64,
    },
}

impl Integration {
    fn parameters(self) -> Result<(f64, f64, f64)> {
        match self {
            Self::AverageAcceleration => Ok((0.25, 0.5, 0.0)),
            Self::LinearAcceleration => Ok((1.0 / 6.0, 0.5, 0.0)),
            Self::Hht { alpha } if alpha.is_finite() && (-1.0 / 3.0..=0.0).contains(&alpha) => Ok(((1.0 - alpha).powi(2) / 4.0, 0.5 - alpha, alpha)),
            _ => Err(Error::Invalid("HHT alpha must be in [-1/3, 0]".into())),
        }
    }
}

#[derive(Debug, Clone)]
pub struct State {
    pub displacement: Vec<f64>,
    pub velocity: Vec<f64>,
    pub acceleration: Vec<f64>,
}

fn vector(values: &[f64], n: usize) -> Result<()> {
    if values.len() != n || values.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("transient vector dimension or finite-value check failed".into()));
    }
    Ok(())
}

fn factor(matrix: &Csr, limits: &Limits) -> Result<Factor> {
    let upper = UpperColumns::from_symmetric(matrix);
    let result = if let Some(bytes) = limits.max_factor_memory_bytes {
        Factor::with_memory_limit(&upper, limits.max_factor_entries, bytes).map(|(f, _)| f)
    } else {
        Factor::new(&upper, limits.max_factor_entries)
    };
    result.map_err(|e| match e {
        FactorError::NotPositiveDefinite { column } | FactorError::Singular { column } => {
            Error::Singular(format!("transient pivot at equation {column}"))
        }
        FactorError::TooLarge { entries, limit } => Error::factor_limit(entries, limit),
        FactorError::MemoryLimit { estimated_bytes, limit_bytes, .. } => {
            Error::Limit(format!("transient factor needs {estimated_bytes} bytes; budget {limit_bytes}"))
        }
        FactorError::OutOfMemory { entries } => Error::Limit(format!("transient factor cannot allocate {entries} entries")),
        FactorError::Internal(message) => Error::Internal(message),
        FactorError::Cancelled => fieldcraft_model::Error::Cancelled.into(),
    })
}

pub struct Integrator<'a> {
    stiffness: &'a Csr,
    damping: &'a Csr,
    effective: Factor,
    dt: f64,
    beta: f64,
    gamma: f64,
    alpha: f64,
}

impl<'a> Integrator<'a> {
    /// Matrices share a symmetric CSR pattern; M and effective operator must be SPD.
    /// No mass regularization is introduced for massless degrees of freedom.
    pub fn new(stiffness: &'a Csr, mass: &Csr, damping: &'a Csr, dt: f64, method: Integration, limits: &Limits) -> Result<Self> {
        crate::harmonic::Matrices { stiffness, mass, viscous: Some(damping), loss: None }.validate(limits)?;
        let (beta, gamma, alpha) = method.parameters()?;
        if !dt.is_finite() || dt <= 0.0 || !(beta * dt * dt).is_finite() || beta * dt * dt == 0.0 {
            return Err(Error::Invalid("transient time step must be positive and representable".into()));
        }
        let effective = mass.add_scaled((1.0 + alpha) * gamma * dt, damping)?.add_scaled((1.0 + alpha) * beta * dt * dt, stiffness)?;
        crate::complex::validate(&effective, limits)?;
        Ok(Self { stiffness, damping, effective: factor(&effective, limits)?, dt, beta, gamma, alpha })
    }

    /// Evaluate before creating the integrator, so the mass and effective factors
    /// do not coexist. f0 is the right-hand limit of the force at t=0.
    pub fn initial_state(
        stiffness: &Csr,
        mass: &Csr,
        damping: &Csr,
        displacement: Vec<f64>,
        velocity: Vec<f64>,
        force: &[f64],
        limits: &Limits,
    ) -> Result<State> {
        crate::harmonic::Matrices { stiffness, mass, viscous: Some(damping), loss: None }.validate(limits)?;
        let n = stiffness.n;
        vector(&displacement, n)?;
        vector(&velocity, n)?;
        vector(force, n)?;
        let mut ku = vec![0.0; n];
        let mut cv = vec![0.0; n];
        stiffness.mul(&displacement, &mut ku);
        damping.mul(&velocity, &mut cv);
        let mut acceleration: Vec<_> = force.iter().zip(ku).zip(cv).map(|((f, k), c)| f - k - c).collect();
        vector(&acceleration, n)?;
        factor(mass, limits)?.solve(&mut acceleration)?;
        vector(&acceleration, n)?;
        Ok(State { displacement, velocity, acceleration })
    }

    /// Transactional: errors leave state unchanged. Forces include prescribed-motion
    /// elimination or -M r a_ground, performed by the model adapter.
    pub fn advance(&self, state: &mut State, force: &[f64], next_force: &[f64]) -> Result<()> {
        fieldcraft_model::progress::check()?;
        let n = self.stiffness.n;
        for values in [&state.displacement[..], &state.velocity, &state.acceleration, force, next_force] {
            vector(values, n)?;
        }
        let u: Vec<_> = (0..n)
            .map(|i| state.displacement[i] + self.dt * state.velocity[i] + self.dt * self.dt * (0.5 - self.beta) * state.acceleration[i])
            .collect();
        let v: Vec<_> = (0..n).map(|i| state.velocity[i] + self.dt * (1.0 - self.gamma) * state.acceleration[i]).collect();
        let mut ku = vec![0.0; n];
        let mut cv = vec![0.0; n];
        let mut old_ku = vec![0.0; n];
        let mut old_cv = vec![0.0; n];
        self.stiffness.mul(&u, &mut ku);
        self.damping.mul(&v, &mut cv);
        self.stiffness.mul(&state.displacement, &mut old_ku);
        self.damping.mul(&state.velocity, &mut old_cv);
        let mut a: Vec<_> =
            (0..n).map(|i| (1.0 + self.alpha) * (next_force[i] - ku[i] - cv[i]) + self.alpha * (old_ku[i] + old_cv[i] - force[i])).collect();
        vector(&a, n)?;
        self.effective.solve(&mut a)?;
        let next_u: Vec<_> = u.iter().zip(&a).map(|(u, a)| u + self.beta * self.dt * self.dt * a).collect();
        let next_v: Vec<_> = v.iter().zip(&a).map(|(v, a)| v + self.gamma * self.dt * a).collect();
        for values in [&next_u, &next_v, &a] {
            vector(values, n)?;
        }
        *state = State { displacement: next_u, velocity: next_v, acceleration: a };
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn scalar(value: f64) -> Csr {
        Csr { n: 1, ptr: vec![0, 1], idx: vec![0], val: vec![value] }
    }

    #[test]
    fn step_force_matches_closed_form_and_physical_initial_acceleration() -> Result<()> {
        let (k, m, c) = (scalar(4.0), scalar(1.0), scalar(0.0));
        let limits = Limits::default();
        let mut state = Integrator::initial_state(&k, &m, &c, vec![0.0], vec![0.0], &[1.0], &limits)?;
        assert_eq!(state.acceleration, vec![1.0]);
        let integrator = Integrator::new(&k, &m, &c, 0.001, Integration::AverageAcceleration, &limits)?;
        for i in 1..=1000 {
            integrator.advance(&mut state, &[1.0], &[1.0])?;
            let exact = (1.0 - libm::cos(2.0 * i as f64 * 0.001)) / 4.0;
            assert!((state.displacement[0] - exact).abs() < 2e-7);
        }
        Ok(())
    }

    #[test]
    fn newmark_free_motion_follows_discrete_period_without_energy_decay() -> Result<()> {
        let (k, m, c) = (scalar(4.0), scalar(1.0), scalar(0.0));
        let limits = Limits::default();
        let mut state = Integrator::initial_state(&k, &m, &c, vec![1.0], vec![0.0], &[0.0], &limits)?;
        let integrator = Integrator::new(&k, &m, &c, 0.1, Integration::AverageAcceleration, &limits)?;
        let theta = 2.0 * libm::atan(0.1);
        for i in 1..=100 {
            integrator.advance(&mut state, &[0.0], &[0.0])?;
            assert!((state.displacement[0] - libm::cos(i as f64 * theta)).abs() < 1e-12);
            let energy = 0.5 * state.velocity[0].powi(2) + 2.0 * state.displacement[0].powi(2);
            assert!((energy - 2.0).abs() < 1e-12);
        }
        Ok(())
    }

    #[test]
    fn damped_free_response_and_nonzero_initial_velocity() -> Result<()> {
        let (k, m, c) = (scalar(4.0), scalar(1.0), scalar(0.4));
        let limits = Limits::default();
        let mut state = Integrator::initial_state(&k, &m, &c, vec![0.3], vec![0.7], &[0.0], &limits)?;
        assert!((state.acceleration[0] + 1.48).abs() < 1e-14);
        let integrator = Integrator::new(&k, &m, &c, 0.001, Integration::Hht { alpha: 0.0 }, &limits)?;
        let wd = libm::sqrt(3.96);
        for i in 1..=1000 {
            integrator.advance(&mut state, &[0.0], &[0.0])?;
            let t = i as f64 * 0.001;
            let exact = libm::exp(-0.2 * t) * (0.3 * libm::cos(wd * t) + 0.76 / wd * libm::sin(wd * t));
            assert!((state.displacement[0] - exact).abs() < 5e-7);
        }
        Ok(())
    }

    #[test]
    fn harmonic_force_matches_zero_initial_condition_solution() -> Result<()> {
        let (k, m, c) = (scalar(4.0), scalar(1.0), scalar(0.0));
        let limits = Limits::default();
        let mut state = Integrator::initial_state(&k, &m, &c, vec![0.0], vec![0.0], &[0.0], &limits)?;
        let integrator = Integrator::new(&k, &m, &c, 0.001, Integration::AverageAcceleration, &limits)?;
        for i in 1..=1000 {
            let t = i as f64 * 0.001;
            integrator.advance(&mut state, &[libm::sin(t - 0.001)], &[libm::sin(t)])?;
            let exact = (libm::sin(t) - 0.5 * libm::sin(2.0 * t)) / 3.0;
            assert!((state.displacement[0] - exact).abs() < 3e-7);
        }
        Ok(())
    }
}
