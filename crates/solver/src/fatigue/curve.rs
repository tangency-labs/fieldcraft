//! User-supplied S–N data: log-log interpolation (piecewise Basquin).
//! Surface and size factors reduce fully reversed strength, not cycle counts.
//! No empirical factor is inferred from elastic properties or geometry.
use super::Basquin;
use crate::{Error, Result};

pub trait StressLife {
    fn validate(&self) -> Result<()>;
    fn damage_per_cycle(&self, amplitude: f64) -> Result<f64>;
}
impl StressLife for Basquin {
    fn validate(&self) -> Result<()> {
        Basquin::validate(self)
    }
    fn damage_per_cycle(&self, amplitude: f64) -> Result<f64> {
        Basquin::damage_per_cycle(self, amplitude)
    }
}

/// Legacy Basquin JSON remains valid. A table explicitly supplies its reduction
/// factors and its extrapolation policy: outside the table use the end slopes.
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(untagged)]
pub enum SnCurve {
    Basquin(Basquin),
    Table(Table),
}
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Table {
    pub points: Vec<Point>,
    pub endurance_limit: f64,
    pub surface_factor: f64,
    pub size_factor: f64,
    pub extrapolation: Extrapolation,
}
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Extrapolation {
    EndSlopes,
    Reject,
}
#[derive(Debug, Clone, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Point {
    pub cycles: f64,
    pub amplitude: f64,
}
impl StressLife for SnCurve {
    fn validate(&self) -> Result<()> {
        match self {
            Self::Basquin(b) => b.validate(),
            Self::Table(t) => {
                if !(2..=1024).contains(&t.points.len())
                    || [t.surface_factor, t.size_factor].iter().any(|v| !v.is_finite() || *v <= 0.0 || *v > 1.0)
                    || !t.endurance_limit.is_finite()
                    || t.endurance_limit < 0.0
                    || t.points.iter().any(|p| !p.cycles.is_finite() || p.cycles <= 0.0 || !p.amplitude.is_finite() || p.amplitude <= 0.0)
                    || t.points.windows(2).any(|w| w[1].cycles <= w[0].cycles || w[1].amplitude >= w[0].amplitude)
                    || t.surface_factor * t.size_factor == 0.0
                {
                    return Err(Error::Invalid(
                        "S-N table needs 2..1024 increasing cycles, decreasing positive Pa amplitudes, nonnegative limit and factors in (0,1]".into(),
                    ));
                }
                Ok(())
            }
        }
    }
    fn damage_per_cycle(&self, amplitude: f64) -> Result<f64> {
        self.validate()?;
        if !amplitude.is_finite() || amplitude < 0.0 {
            return Err(Error::Invalid("invalid fatigue amplitude".into()));
        }
        let Self::Table(t) = self else {
            return match self {
                Self::Basquin(b) => b.damage_per_cycle(amplitude),
                _ => Err(Error::Invalid("invalid S-N law".into())),
            };
        };
        // Mean correction precedes this strength reduction; tensile/yield strengths
        // supplied to the mean law are not multiplied by these fatigue factors.
        let a = amplitude / (t.surface_factor * t.size_factor);
        if !a.is_finite() {
            return Err(Error::Invalid("adjusted fatigue amplitude overflow".into()));
        }
        if a <= t.endurance_limit {
            return Ok(0.0);
        }
        let first = t.points.first().ok_or_else(|| Error::Invalid("empty S-N table".into()))?;
        let last = t.points.last().ok_or_else(|| Error::Invalid("empty S-N table".into()))?;
        if matches!(t.extrapolation, Extrapolation::Reject) && (a > first.amplitude || a < last.amplitude) {
            return Err(Error::Invalid("amplitude outside S-N table; extrapolation rejected".into()));
        }
        let i = t.points.windows(2).position(|w| a >= w[1].amplitude).unwrap_or(t.points.len() - 2);
        let pair = t.points.get(i..i + 2).ok_or_else(|| Error::Invalid("invalid S-N segment".into()))?;
        let [p, q] = pair else {
            return Err(Error::Invalid("invalid S-N segment".into()));
        };
        let fraction = (libm::log(a) - libm::log(p.amplitude)) / (libm::log(q.amplitude) - libm::log(p.amplitude));
        let log_n = libm::log(p.cycles) + fraction * (libm::log(q.cycles) - libm::log(p.cycles));
        let d = libm::exp(-log_n);
        if !d.is_finite() || d == 0.0 {
            return Err(Error::Invalid("S-N damage outside representable range".into()));
        }
        Ok(d)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn table() -> SnCurve {
        SnCurve::Table(Table {
            points: vec![
                Point { cycles: 100., amplitude: 100. },
                Point { cycles: 10000., amplitude: 10. },
                Point { cycles: 10000000., amplitude: 1. },
            ],
            endurance_limit: 0.5,
            surface_factor: 0.8,
            size_factor: 0.5,
            extrapolation: Extrapolation::EndSlopes,
        })
    }
    #[test]
    fn segments_factors_and_limit() {
        let c = table();
        for (a, n) in [(40., 100.), (4., 10000.), (0.4, 10000000.), (4. * 10f64.sqrt(), 1000.), (0.4 * 10f64.sqrt(), 10000. * 1000f64.sqrt())] {
            assert!((c.damage_per_cycle(a).unwrap() * n - 1.).abs() < 1e-12);
        }
        assert_eq!(c.damage_per_cycle(0.2).unwrap(), 0.);
        assert!(c.damage_per_cycle(f64::NAN).is_err());
    }
    #[test]
    fn invalid_tables_and_extrapolation_are_rejected() {
        let SnCurve::Table(mut t) = table() else { panic!() };
        t.extrapolation = Extrapolation::Reject;
        assert!(SnCurve::Table(t.clone()).damage_per_cycle(41.).is_err());
        t.points[1].cycles = 100.;
        assert!(SnCurve::Table(t.clone()).validate().is_err());
        t.points[1].cycles = 10000.;
        t.size_factor = 1.1;
        assert!(SnCurve::Table(t).validate().is_err());
    }
}
