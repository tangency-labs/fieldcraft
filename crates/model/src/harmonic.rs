//! Steady-state harmonic data, with SI frequencies and radians at the model boundary.
//!
//! The time convention is `Re(u exp(i ω t))`. Rayleigh coefficients define
//! C = α M + β K (α in s⁻¹, β in s); the structural loss factor η defines
//! G = η K, so Z = K − ω² M + i(ω C + G). Step and material contributions
//! are additive. Modal ratios add 2 ζ_j ω_j to each mass-normalized modal
//! damping diagonal; they are not a replacement for projected material damping.
//!
//! References: J. K. Vandiver, MIT 2.003SC (2011), lectures 24–25 (modal
//! orthogonality, Rayleigh damping, harmonic forcing); Abaqus 2026 public
//! Analysis/Theory manuals, *Material Damping*, *Damping in Dynamic Analysis*.
//! Constant η is a frequency-domain idealization: its low-frequency limit is
//! complex compliance, not the real static compliance. No hidden η/ω division
//! or discontinuous switch at zero frequency is part of this definition.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{Error, Result, finite};

pub const MAX_FREQUENCIES: usize = 10_000;

/// Additive viscous and structural damping, applicable to a material or step.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Damping {
    #[serde(default)]
    pub alpha: f64,
    #[serde(default)]
    pub beta: f64,
    #[serde(default)]
    pub loss_factor: f64,
}

impl Damping {
    pub fn validate(&self) -> Result<()> {
        for (name, value) in [("Rayleigh alpha", self.alpha), ("Rayleigh beta", self.beta), ("structural loss factor", self.loss_factor)] {
            if finite(name, value)? < 0.0 {
                return Err(Error::Invalid(format!("{name} must not be negative")));
            }
        }
        Ok(())
    }
}

/// Frequencies in cycles per second: an explicit strictly ascending list, or an
/// inclusive range. A list permits explicit refinement around resonances.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Frequencies {
    List(Vec<f64>),
    Range(FrequencyRange),
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FrequencyRange {
    pub start: f64,
    pub end: f64,
    pub count: usize,
    #[serde(default)]
    pub logarithmic: bool,
}

impl Frequencies {
    pub fn points(&self) -> Result<Vec<f64>> {
        let values = match self {
            Self::List(values) => {
                if values.is_empty() || values.len() > MAX_FREQUENCIES {
                    return Err(Error::Invalid(format!("a harmonic step needs 1 to {MAX_FREQUENCIES} frequencies")));
                }
                values.clone()
            }
            Self::Range(r) => {
                if !(2..=MAX_FREQUENCIES).contains(&r.count) {
                    return Err(Error::Invalid(format!("a harmonic range needs 2 to {MAX_FREQUENCIES} points")));
                }
                if !r.start.is_finite() || !r.end.is_finite() || r.start < 0.0 || r.start >= r.end || (r.logarithmic && r.start == 0.0) {
                    return Err(Error::Invalid(
                        "harmonic frequency range must be finite, ascending and nonnegative; logarithmic start must be positive".into(),
                    ));
                }
                (0..r.count)
                    .map(|i| {
                        if i == 0 {
                            return r.start;
                        }
                        if i + 1 == r.count {
                            return r.end;
                        }
                        let fraction = i as f64 / (r.count - 1) as f64;
                        if r.logarithmic {
                            ((1.0 - fraction) * r.start.ln() + fraction * r.end.ln()).exp()
                        } else {
                            (1.0 - fraction) * r.start + fraction * r.end
                        }
                    })
                    .collect()
            }
        };
        if values.iter().any(|f| !f.is_finite() || *f < 0.0 || !(2.0 * std::f64::consts::PI * f).powi(2).is_finite())
            || values.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(Error::Invalid(
                "harmonic frequencies must be finite, nonnegative, strictly ascending and have representable squared angular frequencies".into(),
            ));
        }
        Ok(values)
    }
}

/// The modal method uses existing real Lanczos modes; damping ratios are either
/// empty (zero), one shared value, or one per requested mode.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Method {
    Direct,
    Modal {
        modes: u32,
        #[serde(default)]
        damping_ratios: Vec<f64>,
    },
}

impl Method {
    pub fn validate(&self) -> Result<()> {
        if let Self::Modal { modes, damping_ratios } = self {
            if !(1..=crate::steps::MAX_MODES).contains(modes) {
                return Err(Error::Invalid(format!("harmonic modal count must be 1 to {}", crate::steps::MAX_MODES)));
            }
            if ![0, 1, *modes as usize].contains(&damping_ratios.len()) {
                return Err(Error::Invalid("modal damping ratios must have length 0, 1 or the requested mode count".into()));
            }
            for value in damping_ratios {
                if finite("modal damping ratio", *value)? < 0.0 {
                    return Err(Error::Invalid("modal damping ratio must not be negative".into()));
                }
            }
        }
        Ok(())
    }
}

fn one() -> f64 {
    1.0
}

/// A dimensionless amplitude multiplier and a phase in radians for an existing load.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoadFactor {
    #[serde(default = "one")]
    pub amplitude: f64,
    #[serde(default)]
    pub phase: f64,
}

impl Default for LoadFactor {
    fn default() -> Self {
        Self { amplitude: 1.0, phase: 0.0 }
    }
}

impl LoadFactor {
    pub fn validate(&self) -> Result<()> {
        if finite("harmonic load amplitude", self.amplitude)? < 0.0 {
            return Err(Error::Invalid("harmonic load amplitude must not be negative (use phase for its sign)".into()));
        }
        finite("harmonic load phase", self.phase)?;
        Ok(())
    }

    pub fn components(&self) -> (f64, f64) {
        (self.amplitude * self.phase.cos(), self.amplitude * self.phase.sin())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Harmonic {
    pub frequencies: Frequencies,
    pub method: Method,
    #[serde(default)]
    pub damping: Damping,
    #[serde(default)]
    pub load_factors: BTreeMap<String, LoadFactor>,
}

impl Harmonic {
    pub fn validate(&self) -> Result<()> {
        self.frequencies.points()?;
        self.method.validate()?;
        self.damping.validate()?;
        if self.load_factors.len() > crate::steps::MAX_STEP_ITEMS {
            return Err(Error::Invalid("too many harmonic load factors".into()));
        }
        for (name, factor) in &self.load_factors {
            crate::check_name("harmonic load", name)?;
            factor.validate()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frequency_ranges_keep_endpoints_and_refuse_collapsed_points() {
        let range = |start, end, count, logarithmic| Frequencies::Range(FrequencyRange { start, end, count, logarithmic }).points();
        assert_eq!(range(0.0, 2.0, 3, false).unwrap(), vec![0.0, 1.0, 2.0]);
        let log = range(1.0, 100.0, 3, true).unwrap();
        assert!((log[1] - 10.0).abs() < 1e-13);
        for (start, end, count, logarithmic) in
            [(0.0, 10.0, 5, true), (0.0, 1.0, usize::MAX, false), (2.0, 1.0, 3, false), (1.0, 1.0 + f64::EPSILON, 10, false)]
        {
            assert!(range(start, end, count, logarithmic).is_err());
        }
        for list in [vec![], vec![0.0, 0.0], vec![1.0, 0.0], vec![-1.0], vec![f64::NAN], vec![f64::MAX]] {
            assert!(Frequencies::List(list).points().is_err());
        }
    }

    #[test]
    fn invalid_damping_and_load_factors_are_refused() {
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(Damping { alpha: invalid, ..Damping::default() }.validate().is_err());
            assert!(Damping { beta: invalid, ..Damping::default() }.validate().is_err());
            assert!(Damping { loss_factor: invalid, ..Damping::default() }.validate().is_err());
            assert!(Method::Modal { modes: 2, damping_ratios: vec![invalid] }.validate().is_err());
        }
        assert!(Method::Modal { modes: 3, damping_ratios: vec![0.1; 2] }.validate().is_err());
        assert!(LoadFactor { amplitude: -1.0, phase: 0.0 }.validate().is_err());
        assert!(LoadFactor { amplitude: 1.0, phase: f64::INFINITY }.validate().is_err());
    }
}
