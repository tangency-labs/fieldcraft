//! Local displacement estimators for smooth, physically equilibrated Newmark states.
//!
//! Zienkiewicz & Xie (1991), doi:10.1002/eqe.4290200907, as reproduced
//! by Mayr/Wall/Gee, Adaptive time stepping for FSI solvers, equation (6).
//! The quadratic acceleration reconstruction follows the higher-order auxiliary
//! integration premise of Li/Zeng/Wiberg (1993), doi:10.1002/cnm.1640090402.
//! Its equal-step limit is independently given in Neumann/Schweizerhof (1999),
//! Analysis of Shell Structures under Transient Loading, equations (32)-(33).
//! We derive the unequal-step extension by integrating the unique quadratic
//! acceleration interpolant; no claim of verbatim reproduction of the 1993 paper.
//!
//! For h=current step and g=previous step:
//! c=((a1-a0)/h-(a0-am)/g)/(h+g),
//! u_aux=u0+h*v0+h*h*a0/2+h*h*(a1-a0)/6-c*h^4/12.
//! Thus u_Newmark-u_aux=(beta-1/6)*h*h*(a1-a0)+c*h^4/12.
//! The first step uses the two-point estimator. These are local smooth-solution
//! estimates, not bounds on global error; HHT weighted equilibrium is excluded.

use crate::{Error, Result};

use fieldcraft_model::dynamics::ErrorEstimator as Estimator;

pub fn displacement_error(estimator: Estimator, beta: f64, h: f64, a0: &[f64], a1: &[f64], previous: Option<(&[f64], f64)>) -> Result<Vec<f64>> {
    if !h.is_finite() || h <= 0.0 || !beta.is_finite() || a0.len() != a1.len() || a0.iter().chain(a1).any(|v| !v.is_finite()) {
        return Err(Error::Invalid("adaptive estimator requires finite equal-length states and positive step".into()));
    }
    if let Some((am, g)) = previous
        && (am.len() != a0.len() || !g.is_finite() || g <= 0.0 || am.iter().any(|v| !v.is_finite()))
    {
        return Err(Error::Invalid("adaptive estimator previous state or step is invalid".into()));
    }
    a0.iter()
        .zip(a1)
        .enumerate()
        .map(|(i, (a, b))| {
            let mut value = (beta - 1.0 / 6.0) * h * h * (b - a);
            if let Estimator::QuadraticAcceleration = estimator
                && let Some((am, g)) = previous
                && let Some(old) = am.get(i)
            {
                value += (((b - a) / h - (a - old) / g) / (h + g)) * h.powi(4) / 12.0;
            }
            if value.is_finite() { Ok(value) } else { Err(Error::Invalid("adaptive estimator overflow".into())) }
        })
        .collect()
}

/// Persistent controller; accepted history only. A rejected trial never mutates
/// the solution, forcing vector, acceleration history or work ledger.
pub struct Controller {
    pub accepted: u32,
    pub rejected: u32,
    pub minimum_used: f64,
    pub maximum_used: f64,
    next_step: f64,
    previous: Option<(Vec<f64>, f64)>,
}
impl Controller {
    pub fn new(step: f64) -> Self {
        Self { accepted: 0, rejected: 0, minimum_used: f64::INFINITY, maximum_used: 0.0, next_step: step, previous: None }
    }
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to(
        &mut self,
        state: &mut super::State,
        force: &mut Vec<f64>,
        start: f64,
        end: f64,
        stiffness: &crate::sparse::Csr,
        mass: &crate::sparse::Csr,
        damping: &crate::sparse::Csr,
        settings: &fieldcraft_model::dynamics::Adaptive,
        knots: &[f64],
        evaluate: &mut impl FnMut(f64) -> Result<Vec<f64>>,
        limits: &crate::Limits,
    ) -> Result<(f64, f64)> {
        let mut time = start;
        let mut dissipation = 0.0;
        let mut work = 0.0;
        while time < end {
            fieldcraft_model::progress::check()?;
            if self.accepted.saturating_add(self.rejected) >= settings.max_attempts {
                return Err(Error::Limit("adaptive transient maximum attempts exhausted".into()));
            }
            let boundary = knots.get(knots.partition_point(|t| *t <= time)).copied().filter(|t| *t < end).unwrap_or(end);
            let h = self.next_step.clamp(settings.min_step, settings.max_step).min(boundary - time);
            if !h.is_finite() || h <= 0.0 || time + h <= time {
                return Err(Error::Invalid("adaptive transient time increment is not representable".into()));
            }
            let next_time = if h == boundary - time { boundary } else { time + h };
            let next_force = evaluate(next_time)?;
            let mut trial = state.clone();
            super::Integrator::new(stiffness, mass, damping, h, super::Integration::AverageAcceleration, limits)?.advance(
                &mut trial,
                force,
                &next_force,
            )?;
            let estimate = displacement_error(
                settings.estimator,
                0.25,
                h,
                &state.acceleration,
                &trial.acceleration,
                self.previous.as_ref().map(|(a, g)| (a.as_slice(), *g)),
            )?;
            let mut ratio = 0.0_f64;
            for (e, (a, b)) in estimate.iter().zip(state.displacement.iter().zip(&trial.displacement)) {
                let scale = settings.absolute_tolerance + settings.relative_tolerance * a.abs().max(b.abs());
                let value = e.abs() / scale;
                if !scale.is_finite() || !value.is_finite() {
                    return Err(Error::Invalid("adaptive tolerance scaling overflow".into()));
                }
                ratio = ratio.max(value);
            }
            if !ratio.is_finite() {
                return Err(Error::Invalid("adaptive error ratio is not finite".into()));
            }
            let factor = if ratio == 0.0 { 2.0 } else { (0.9 / libm::cbrt(ratio)).clamp(0.2, 2.0) };
            self.next_step = (h * factor).clamp(settings.min_step, settings.max_step);
            if ratio > 1.0 {
                self.rejected += 1;
                if h <= settings.min_step || self.next_step >= h {
                    return Err(Error::Limit("adaptive transient tolerance cannot be met at min_step".into()));
                }
                continue;
            }
            let midpoint_v: Vec<_> = state.velocity.iter().zip(&trial.velocity).map(|(a, b)| 0.5 * (a + b)).collect();
            let mut cv = vec![0.0; midpoint_v.len()];
            damping.mul(&midpoint_v, &mut cv);
            dissipation += h * midpoint_v.iter().zip(cv).map(|(v, c)| v * c).sum::<f64>();
            work += force
                .iter()
                .zip(&next_force)
                .zip(state.displacement.iter().zip(&trial.displacement))
                .map(|((a, b), (u, v))| 0.5 * (a + b) * (v - u))
                .sum::<f64>();
            self.previous = if knots.binary_search_by(|v| v.total_cmp(&next_time)).is_ok() { None } else { Some((state.acceleration.clone(), h)) };
            *state = trial;
            *force = next_force;
            time = next_time;
            self.accepted += 1;
            self.minimum_used = self.minimum_used.min(h);
            self.maximum_used = self.maximum_used.max(h);
        }
        Ok((dissipation, work))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quadratic_acceleration_matches_integrated_polynomial_and_published_equal_step_limit() {
        let h = 0.2_f64;
        let g = 0.3_f64;
        let e = displacement_error(Estimator::QuadraticAcceleration, 0.25, h, &[0.0], &[h * h], Some((&[g * g], g))).unwrap();
        assert!((e[0] - (0.25 * h.powi(4) - h.powi(4) / 12.0)).abs() < 1e-16);
        let e = displacement_error(Estimator::QuadraticAcceleration, 0.25, h, &[0.4], &[0.6], Some((&[0.3], h))).unwrap();
        // Neumann/Schweizerhof (1999), equation 32, beta=1/4.
        let published = h * h / 24.0 * (0.3 - 4.0 * 0.4 + 3.0 * 0.6);
        assert!((e[0] - published).abs() < 1e-16);
        let zx = displacement_error(Estimator::ZienkiewiczXie, 0.25, h, &[0.4], &[0.6], None).unwrap();
        assert!((zx[0] - (0.25 - 1.0 / 6.0) * h * h * 0.2).abs() < 1e-16);
    }
}
