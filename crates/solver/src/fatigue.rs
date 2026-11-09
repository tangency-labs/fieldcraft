//! Finite-record fatigue counting, from ASTM E1049 §5.4.4 (three-point rainflow).
//!
//! Stress is in Pa. Ranges are peak-to-valley; amplitudes are half ranges.
//! Residual excursions count as half cycles: records are never implicitly periodic.
//! Basquin uses cycles, not reversals; Miner sums count / cycles-to-failure.
//! This kernel does not infer material fatigue data from elastic properties.

use crate::{Error, Result};

pub mod curve;
pub mod field;
pub mod history;
pub mod life;
pub mod spectral;

/// Bounded input and output storage, independent of user-supplied histories.
pub const MAX_HISTORY_POINTS: usize = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cycle {
    pub range: f64,
    pub mean: f64,
    pub count: f64,
}

fn cycle(a: f64, b: f64, count: f64) -> Result<Cycle> {
    let range = (b - a).abs();
    let mean = a * 0.5 + b * 0.5;
    if !range.is_finite() || !mean.is_finite() {
        return Err(Error::Invalid("fatigue cycle exceeds finite stress range".into()));
    }
    Ok(Cycle { range, mean, count })
}

/// Remove plateaus and monotonic interior points before counting. Every loop is
/// bounded by input size; every inner-loop reduction removes at least one point.
pub fn rainflow(samples: &[f64]) -> Result<Vec<Cycle>> {
    if samples.len() < 2 || samples.len() > MAX_HISTORY_POINTS {
        return Err(Error::Invalid("fatigue history needs 2..=1000000 samples".into()));
    }
    if samples.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("fatigue history contains non-finite stress".into()));
    }
    let mut turns: Vec<f64> = Vec::new();
    for &x in samples {
        if turns.last() == Some(&x) {
            continue;
        }
        if let [.., a, b] = turns.as_slice()
            && ((b > a && x > *b) || (b < a && x < *b))
        {
            turns.pop();
        }
        turns.push(x);
    }
    let mut stack: Vec<f64> = Vec::new();
    let mut cycles = Vec::new();
    for x in turns {
        stack.push(x);
        while let [.., a, b, c] = stack.as_slice() {
            let (a, b, c) = (*a, *b, *c);
            let previous = cycle(a, b, 1.0)?;
            let next = cycle(b, c, 1.0)?;
            if next.range < previous.range {
                break;
            }
            if stack.len() == 3 {
                cycles.push(Cycle { count: 0.5, ..previous });
                stack.remove(0);
            } else {
                cycles.push(previous);
                stack.truncate(stack.len() - 3);
                stack.push(c);
            }
        }
    }
    for pair in stack.windows(2) {
        if let [a, b] = pair {
            cycles.push(cycle(*a, *b, 0.5)?);
        }
    }
    Ok(cycles)
}

/// N(a) = cycles_ref * (amplitude_ref / a)^exponent (Basquin).
/// The reference amplitude and any endurance limit are fully reversed stresses.
#[derive(Debug, Clone, Copy, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Basquin {
    pub amplitude_ref: f64,
    pub cycles_ref: f64,
    pub exponent: f64,
    pub endurance_limit: f64,
}

impl Basquin {
    pub fn validate(&self) -> Result<()> {
        if [self.amplitude_ref, self.cycles_ref, self.exponent].iter().any(|x| !x.is_finite() || *x <= 0.0)
            || !self.endurance_limit.is_finite()
            || self.endurance_limit < 0.0
        {
            return Err(Error::Invalid("fatigue Basquin parameters must be finite and positive; endurance limit nonnegative".into()));
        }
        Ok(())
    }

    /// Damage per full cycle. Logarithms avoid overflowing an intermediate life.
    pub fn damage_per_cycle(&self, amplitude: f64) -> Result<f64> {
        self.validate()?;
        if !amplitude.is_finite() || amplitude < 0.0 {
            return Err(Error::Invalid("fatigue amplitude must be finite and nonnegative".into()));
        }
        if amplitude <= self.endurance_limit {
            return Ok(0.0);
        }
        let log_damage = self.exponent * (libm::log(amplitude) - libm::log(self.amplitude_ref)) - libm::log(self.cycles_ref);
        let damage = libm::exp(log_damage);
        if !damage.is_finite() || damage == 0.0 {
            return Err(Error::Invalid("fatigue damage outside representable positive range".into()));
        }
        Ok(damage)
    }
}

/// Palmgren–Miner damage without mean-stress correction. The caller must select
/// and record a mean-stress policy before using this for a physical assessment.
pub fn miner_uncorrected(cycles: &[Cycle], curve: &Basquin) -> Result<f64> {
    miner(cycles, curve, life::StressMean::None)
}

/// Miner damage with a caller-selected mean-stress law, applied per cycle.
pub fn miner(cycles: &[Cycle], curve: &impl curve::StressLife, correction: life::StressMean) -> Result<f64> {
    curve.validate()?;
    if cycles.len() > MAX_HISTORY_POINTS {
        return Err(Error::Invalid("fatigue cycle count exceeds storage limit".into()));
    }
    let mut sum = 0.0;
    for c in cycles {
        if !c.range.is_finite() || c.range < 0.0 || !c.mean.is_finite() || !c.count.is_finite() || c.count <= 0.0 {
            return Err(Error::Invalid("invalid fatigue cycle".into()));
        }
        sum += c.count * curve.damage_per_cycle(correction.corrected_amplitude(c.range * 0.5, c.mean)?)?;
        if !sum.is_finite() {
            return Err(Error::Invalid("fatigue Miner damage overflow".into()));
        }
    }
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve() -> Basquin {
        Basquin { amplitude_ref: 100_000_000.0, cycles_ref: 1_000_000.0, exponent: 3.0, endurance_limit: 0.0 }
    }

    #[test]
    fn finite_triangle_counts_two_halves_and_known_damage() {
        let cycles = rainflow(&[-100_000_000.0, 100_000_000.0, -100_000_000.0]).unwrap();
        assert_eq!(cycles, vec![Cycle { range: 200_000_000.0, mean: 0.0, count: 0.5 }; 2]);
        assert!((miner_uncorrected(&cycles, &curve()).unwrap() / 1e-6 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nested_excursion_and_residue_are_separate() {
        let cycles = rainflow(&[0.0, 4.0, 1.0, 3.0, 0.0]).unwrap();
        assert_eq!(
            cycles,
            vec![
                Cycle { range: 2.0, mean: 2.0, count: 1.0 },
                Cycle { range: 4.0, mean: 2.0, count: 0.5 },
                Cycle { range: 4.0, mean: 2.0, count: 0.5 },
            ]
        );
    }

    #[test]
    fn plateaus_and_constant_records() {
        assert!(rainflow(&[2.0, 2.0, 2.0]).unwrap().is_empty());
        assert_eq!(rainflow(&[0.0, 1.0, 1.0, 2.0, 1.0, 0.0]).unwrap(), rainflow(&[0.0, 2.0, 0.0]).unwrap());
        assert_eq!(rainflow(&[0.0, 1.0, 2.0]).unwrap(), vec![Cycle { range: 2.0, mean: 1.0, count: 0.5 }]);
    }

    #[test]
    fn hostile_inputs_are_refused() {
        assert!(rainflow(&[0.0, f64::NAN]).is_err());
        assert!(rainflow(&[-f64::MAX, f64::MAX]).is_err());
        assert!(rainflow(&[]).is_err());
        assert!(curve().damage_per_cycle(-1.0).is_err());
        assert_eq!(curve().damage_per_cycle(0.0).unwrap(), 0.0);
    }
}

#[cfg(test)]
mod published_counting {
    use super::*;
    /// ASTM E1049-85(2017), §5.4.4.2 / Fig.6 / Table X1.3.
    /// Numerical reversal/cycle facts only; no copied figure, prose or implementation.
    #[test]
    fn astm_e1049_fig6_each_cycle() {
        let counted = rainflow(&[-2., 1., -3., 5., -1., 3., -4., 4., -2.]).unwrap();
        let expected = [(3., -0.5, 0.5), (4., -1., 0.5), (4., 1., 1.), (8., 1., 0.5), (9., 0.5, 0.5), (8., 0., 0.5), (6., 1., 0.5)];
        assert_eq!(counted.len(), expected.len());
        for (c, (range, mean, count)) in counted.iter().zip(expected) {
            assert_eq!((c.range, c.mean, c.count), (range, mean, count));
        }
    }
}
