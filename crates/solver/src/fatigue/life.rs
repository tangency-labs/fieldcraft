//! Uniaxial mean-stress and strain-life equations: Stephens et al., Metal
//! Fatigue in Engineering, 2nd ed. (2000), strain-life approach; Smith,
//! Watson & Topper (1970), A stress-strain function for the fatigue of metals.
//! Public equation reference: Altair, Uniaxial Fatigue Analysis, strain-life
//! approach (Morrow and SWT). Stress is Pa; strain dimensionless; N is cycles.
use crate::{Error, Result};

#[derive(Debug, Clone, Copy, serde::Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum StressMean {
    None,
    Goodman { ultimate: f64 },
    Gerber { ultimate: f64 },
    Soderberg { yield_stress: f64 },
}

impl StressMean {
    /// Compressive mean stress receives no beneficial credit. This explicit
    /// conservative policy differs from extrapolating the tensile relations.
    pub fn corrected_amplitude(self, amplitude: f64, mean: f64) -> Result<f64> {
        if !amplitude.is_finite() || amplitude < 0.0 || !mean.is_finite() {
            return Err(Error::Invalid("fatigue amplitude/mean invalid".into()));
        }
        let (strength, quadratic) = match self {
            Self::None => return Ok(amplitude),
            Self::Goodman { ultimate } => (ultimate, false),
            Self::Gerber { ultimate } => (ultimate, true),
            Self::Soderberg { yield_stress } => (yield_stress, false),
        };
        if !strength.is_finite() || strength <= 0.0 || mean >= strength {
            return Err(Error::Invalid("mean-stress correction requires positive strength and mean below strength".into()));
        }
        let ratio = mean.max(0.0) / strength;
        let value = amplitude / (1.0 - if quadratic { ratio * ratio } else { ratio });
        if !value.is_finite() {
            return Err(Error::Invalid("corrected fatigue amplitude overflow".into()));
        }
        Ok(value)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct StrainLife {
    pub young: f64,
    pub strength: f64,
    pub ductility: f64,
    pub b: f64,
    pub c: f64,
}

#[derive(Debug, Clone, Copy)]
pub enum StrainMean {
    None,
    Morrow { mean_stress: f64 },
    Swt { maximum_stress: f64 },
}

impl StrainLife {
    /// Solve in log(2N), over 0.5..=1e16 cycles. Outside this declared
    /// numerical bracket is an error, not an endurance-limit assertion.
    /// None means zero strain amplitude only, not demonstrated infinite life.
    /// Total local strain and stress must describe the SAME physical cycle.
    /// One coefficient may be zero for the exact Basquin or Coffin-Manson limit;
    /// mean corrections require a positive fatigue strength coefficient.
    /// These limits are equation checks, not estimates of missing material data.
    pub fn cycles(self, amplitude: f64, correction: StrainMean) -> Result<Option<f64>> {
        if !self.young.is_finite()
            || self.young <= 0.0
            || [self.strength, self.ductility].iter().any(|x| !x.is_finite() || *x < 0.0)
            || (self.strength == 0.0 && self.ductility == 0.0)
            || [self.b, self.c].iter().any(|x| !x.is_finite() || *x >= 0.0)
            || !amplitude.is_finite()
            || amplitude < 0.0
        {
            return Err(Error::Invalid("invalid Coffin-Manson-Basquin parameters or total strain amplitude".into()));
        }
        let (a, b, c, d, target) = match correction {
            StrainMean::None => (self.strength / self.young, self.b, self.ductility, self.c, amplitude),
            StrainMean::Morrow { mean_stress } => {
                if self.strength == 0.0 || !mean_stress.is_finite() || mean_stress >= self.strength {
                    return Err(Error::Invalid("Morrow mean stress must be finite and below fatigue strength coefficient".into()));
                }
                ((self.strength - mean_stress.max(0.0)) / self.young, self.b, self.ductility, self.c, amplitude)
            }
            StrainMean::Swt { maximum_stress } => {
                if self.strength == 0.0 || !maximum_stress.is_finite() || maximum_stress <= 0.0 {
                    return Err(Error::Invalid(
                        "SWT requires positive fatigue strength and tensile maximum stress; compression-only cycles unsupported".into(),
                    ));
                }
                (
                    self.strength / self.young * self.strength,
                    2.0 * self.b,
                    self.strength * self.ductility,
                    self.b + self.c,
                    maximum_stress * amplitude,
                )
            }
        };
        if [a, b, c, d, target].iter().any(|x| !x.is_finite()) || a < 0.0 || c < 0.0 || (a == 0.0 && c == 0.0) {
            return Err(Error::Invalid("strain-life equation overflow or underflow".into()));
        }
        if amplitude == 0.0 {
            return Ok(None);
        }
        if target <= 0.0 {
            return Err(Error::Invalid("strain-life target underflow".into()));
        }
        let value = |x: f64| a * libm::exp(b * x) + c * libm::exp(d * x);
        let mut lo = 0.0;
        let mut hi = libm::log(2e16);
        if !value(lo).is_finite() || target > value(lo) || target < value(hi) {
            return Err(Error::Invalid("strain-life solution outside 0.5..=1e16 cycles".into()));
        }
        for _ in 0..128 {
            let mid = (lo + hi) * 0.5;
            if value(mid) > target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let x = (lo + hi) * 0.5;
        if (value(x) / target - 1.0).abs() > 1e-10 {
            return Err(Error::NotConverged("strain-life inversion residual".into()));
        }
        Ok(Some(0.5 * libm::exp(x)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mean_stress_known_factors_and_domain() {
        assert_eq!(StressMean::Goodman { ultimate: 400.0 }.corrected_amplitude(100.0, 200.0).unwrap(), 200.0);
        assert_eq!(StressMean::Gerber { ultimate: 400.0 }.corrected_amplitude(150.0, 200.0).unwrap(), 200.0);
        assert_eq!(StressMean::Soderberg { yield_stress: 200.0 }.corrected_amplitude(100.0, 100.0).unwrap(), 200.0);
        assert_eq!(StressMean::Goodman { ultimate: 400.0 }.corrected_amplitude(100.0, -200.0).unwrap(), 100.0);
        assert!(StressMean::Goodman { ultimate: 400.0 }.corrected_amplitude(0.0, 400.0).is_err());
    }
    #[test]
    fn independent_closed_form_strain_lives() {
        // b=c=-1/2: epsilon=(sigma_f/E+epsilon_f)/sqrt(2N).
        let m = StrainLife { young: 200e9, strength: 1e9, ductility: 0.5, b: -0.5, c: -0.5 };
        for (amplitude, correction) in
            [(0.00505, StrainMean::None), (0.005025, StrainMean::Morrow { mean_stress: 0.5e9 }), (0.00505, StrainMean::Swt { maximum_stress: 1e7 })]
        {
            let n = m.cycles(amplitude, correction).unwrap().unwrap();
            assert!((n / 5000.0 - 1.0).abs() < 1e-10);
        }
        assert!(m.cycles(0.01, StrainMean::Swt { maximum_stress: -1.0 }).is_err());
        assert!(m.cycles(1.0, StrainMean::None).is_err());
        assert_eq!(m.cycles(0.0, StrainMean::None).unwrap(), None);
    }
}
