//! Dimensionless amplitudes evaluated in step-local SI seconds.
//!
//! Tables hold their endpoint values outside their domain. Smooth steps use the
//! quintic 10s³ − 15s⁴ + 6s⁵, with zero first and second derivatives at knots.
//! A periodic curve is a Fourier series in angular frequency (rad/s), not Hz.

use serde::{Deserialize, Serialize};

use crate::{Error, Result};

pub const MAX_AMPLITUDE_POINTS: usize = 100_000;
pub const MAX_HARMONICS: usize = 1024;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Amplitude {
    Tabular {
        points: Vec<[f64; 2]>,
    },
    SmoothStep {
        points: Vec<[f64; 2]>,
    },
    /// offset + sum(a_n cos(n ω (t-origin)) + b_n sin(n ω (t-origin))).
    Periodic {
        angular_frequency: f64,
        origin: f64,
        offset: f64,
        harmonics: Vec<[f64; 2]>,
    },
}

impl Amplitude {
    pub fn validate(&self) -> Result<()> {
        match self {
            Self::Tabular { points } | Self::SmoothStep { points } => {
                if !(2..=MAX_AMPLITUDE_POINTS).contains(&points.len())
                    || points.iter().flatten().any(|x| !x.is_finite())
                    || points.windows(2).any(|w| w[0][0] >= w[1][0] || !(w[1][0] - w[0][0]).is_finite())
                {
                    return Err(Error::Invalid("amplitude needs finite, strictly ascending time/value pairs (2 to 100000)".into()));
                }
            }
            Self::Periodic { angular_frequency, origin, offset, harmonics } => {
                if !angular_frequency.is_finite()
                    || *angular_frequency <= 0.0
                    || !origin.is_finite()
                    || !offset.is_finite()
                    || harmonics.is_empty()
                    || harmonics.len() > MAX_HARMONICS
                    || harmonics.iter().flatten().any(|x| !x.is_finite())
                {
                    return Err(Error::Invalid(
                        "periodic amplitude needs positive finite angular frequency and finite Fourier coefficients (1 to 1024)".into(),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Validate once before a time loop, then use the borrowed evaluator.
    pub fn evaluator(&self) -> Result<Evaluator<'_>> {
        self.validate()?;
        Ok(Evaluator(self))
    }
}

pub struct Evaluator<'a>(&'a Amplitude);

impl Evaluator<'_> {
    pub fn value(&self, time: f64) -> Result<f64> {
        if !time.is_finite() {
            return Err(Error::Invalid("amplitude time must be finite".into()));
        }
        let value = match self.0 {
            Amplitude::Tabular { points } | Amplitude::SmoothStep { points } => {
                let first = points.first().ok_or_else(|| Error::Invalid("empty amplitude".into()))?;
                let last = points.last().ok_or_else(|| Error::Invalid("empty amplitude".into()))?;
                if time <= first[0] {
                    first[1]
                } else if time >= last[0] {
                    last[1]
                } else {
                    let upper = points.partition_point(|p| p[0] <= time);
                    let a = points.get(upper - 1).ok_or_else(|| Error::Invalid("amplitude interval".into()))?;
                    let b = points.get(upper).ok_or_else(|| Error::Invalid("amplitude interval".into()))?;
                    let s = (time - a[0]) / (b[0] - a[0]);
                    let weight = if matches!(self.0, Amplitude::SmoothStep { .. }) { s * s * s * (10.0 + s * (-15.0 + 6.0 * s)) } else { s };
                    (1.0 - weight) * a[1] + weight * b[1]
                }
            }
            Amplitude::Periodic { angular_frequency, origin, offset, harmonics } => {
                let phase = angular_frequency * (time - origin);
                let mut value = *offset;
                for (i, [a, b]) in harmonics.iter().enumerate() {
                    let angle = (i + 1) as f64 * phase;
                    if !angle.is_finite() {
                        return Err(Error::Invalid("periodic amplitude phase overflow".into()));
                    }
                    value += a * libm::cos(angle) + b * libm::sin(angle);
                }
                value
            }
        };
        if !value.is_finite() {
            return Err(Error::Invalid("amplitude value overflow".into()));
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_hold_endpoints_and_allow_signed_values() -> Result<()> {
        let curve = Amplitude::Tabular { points: vec![[0.0, -2.0], [2.0, 4.0], [3.0, 1.0]] };
        let e = curve.evaluator()?;
        for (t, expected) in [(-1.0, -2.0), (0.0, -2.0), (1.0, 1.0), (2.0, 4.0), (4.0, 1.0)] {
            assert_eq!(e.value(t)?, expected);
        }
        Ok(())
    }

    #[test]
    fn smooth_step_and_fourier_conventions() -> Result<()> {
        let smooth = Amplitude::SmoothStep { points: vec![[0.0, 0.0], [1.0, 1.0]] };
        assert_eq!(smooth.evaluator()?.value(0.5)?, 0.5);
        assert!((smooth.evaluator()?.value(0.25)? - 0.103515625).abs() < 1e-15);
        let sine = Amplitude::Periodic { angular_frequency: std::f64::consts::PI, origin: 1.0, offset: 2.0, harmonics: vec![[0.0, 3.0]] };
        assert_eq!(sine.evaluator()?.value(1.0)?, 2.0);
        assert!((sine.evaluator()?.value(1.5)? - 5.0).abs() < 1e-14);
        Ok(())
    }

    #[test]
    fn invalid_tables_and_evaluation_are_errors() -> Result<()> {
        for points in [vec![], vec![[0.0, 1.0]], vec![[0.0, 1.0], [0.0, 2.0]], vec![[0.0, 1.0], [f64::NAN, 2.0]]] {
            assert!(Amplitude::Tabular { points }.validate().is_err());
        }
        let curve = Amplitude::Tabular { points: vec![[0.0, 1.0], [1.0, 2.0]] };
        assert!(curve.evaluator()?.value(f64::INFINITY).is_err());
        Ok(())
    }
}
