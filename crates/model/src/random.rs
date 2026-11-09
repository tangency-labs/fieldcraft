//! Random-vibration analysis-step data.  PSDs are one-sided and use Hz at the model boundary.

use serde::{Deserialize, Serialize};

use crate::{Error, Result, finite};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpectrumPoint {
    pub frequency: f64,
    pub psd: f64,
}

impl SpectrumPoint {
    pub fn validate(&self) -> Result<()> {
        if finite("random-vibration frequency", self.frequency)? < 0.0 || finite("random-vibration PSD", self.psd)? < 0.0 {
            return Err(Error::Invalid("random-vibration frequency and PSD must be non-negative".into()));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Excitation {
    Force {
        load: String,
        spectrum: Vec<SpectrumPoint>,
    },
    /// One scalar random multiplier on the sum of this step's spatial loads. Its PSD
    /// is per Hz; the physical load units are already carried by the spatial pattern.
    CorrelatedLoads {
        spectrum: Vec<SpectrumPoint>,
    },
    BaseAcceleration {
        dof: u8,
        spectrum: Vec<SpectrumPoint>,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RandomVibration {
    #[serde(default)]
    pub frequencies: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub modal: Option<ModalRandom>,
    pub excitations: Vec<Excitation>,
    #[serde(default)]
    pub damping: crate::harmonic::Damping,
}

impl RandomVibration {
    pub fn validate(&self) -> Result<()> {
        if let Some(modal) = &self.modal {
            if !self.frequencies.is_empty() {
                return Err(Error::Invalid("choose modal sampling or explicit frequencies, not both".into()));
            }
            modal.validate()?;
        }
        if self.modal.is_none()
            && (self.frequencies.len() < 2 || self.frequencies.windows(2).any(|w| !w[0].is_finite() || w[0] < 0.0 || w[1] <= w[0]))
        {
            return Err(Error::Invalid("random-vibration frequencies must be finite, non-negative and strictly increasing".into()));
        }
        if self.excitations.is_empty() {
            return Err(Error::Invalid("random-vibration needs at least one excitation".into()));
        }
        for excitation in &self.excitations {
            let spectrum = match excitation {
                Excitation::CorrelatedLoads { spectrum } | Excitation::Force { spectrum, .. } | Excitation::BaseAcceleration { spectrum, .. } => {
                    spectrum
                }
            };
            if spectrum.len() < 2 {
                return Err(Error::Invalid("random-vibration PSD needs at least two points".into()));
            }
            for point in spectrum {
                point.validate()?;
            }
            if spectrum.windows(2).any(|w| w[1].frequency <= w[0].frequency) {
                return Err(Error::Invalid("random-vibration PSD frequencies must be strictly increasing".into()));
            }
        }
        self.damping.validate()
    }
}

/// Modal interval sampling from the public Abaqus random-response analysis guide.
/// https://abaqus.uclouvain.be/English/SIMACAEANLRefMap/simaanl-c-randomresponse.htm
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModalRandom {
    pub source_step: String,
    pub modes: u32,
    pub damping_ratios: Vec<f64>,
    pub bands: Vec<ModalBand>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModalBand {
    pub lower: f64,
    pub upper: f64,
    pub points: usize,
    pub bias: f64,
    pub logarithmic: bool,
}
impl ModalRandom {
    pub fn validate(&self) -> Result<()> {
        crate::check_name("random modal source step", &self.source_step)?;
        crate::harmonic::Method::Modal { modes: self.modes, damping_ratios: self.damping_ratios.clone(), residual_vectors: 0 }.validate()?;
        if self.bands.is_empty() || self.bands.len() > 100 {
            return Err(Error::Invalid("random modal sampling needs 1..100 bands".into()));
        }
        for b in &self.bands {
            if !b.lower.is_finite()
                || !b.upper.is_finite()
                || b.lower < 0.0
                || b.upper <= b.lower
                || (b.logarithmic && b.lower == 0.0)
                || !(2..=crate::harmonic::MAX_FREQUENCIES).contains(&b.points)
                || !b.bias.is_finite()
                || b.bias <= 0.0
            {
                return Err(Error::Invalid("invalid random modal frequency band".into()));
            }
        }
        if self.bands.windows(2).any(|w| w[0].upper != w[1].lower) {
            return Err(Error::Invalid("random frequency bands must be contiguous and ascending".into()));
        }
        Ok(())
    }
    pub fn frequencies(&self, natural: &[f64]) -> Result<Vec<f64>> {
        self.validate()?;
        if natural.len() != self.modes as usize || natural.iter().any(|f| !f.is_finite() || *f <= 0.0) {
            return Err(Error::Invalid(
                "random response needs every requested positive natural frequency; rigid or missing modes are not mapped".into(),
            ));
        }
        let mut out = Vec::new();
        for b in &self.bands {
            let mut ends = vec![b.lower];
            ends.extend(natural.iter().copied().filter(|f| *f > b.lower && *f < b.upper));
            ends.push(b.upper);
            ends.sort_by(f64::total_cmp);
            ends.dedup();
            for pair in ends.windows(2) {
                crate::progress::check()?;
                let a = if b.logarithmic { libm::log(pair[0]) } else { pair[0] };
                let z = if b.logarithmic { libm::log(pair[1]) } else { pair[1] };
                for k in 0..b.points {
                    let y = -1.0 + 2.0 * k as f64 / (b.points - 1) as f64;
                    let fraction = if b.points < 4 { (y + 1.0) / 2.0 } else { (1.0 + libm::pow(y.abs(), 1.0 / b.bias) * y.signum()) / 2.0 };
                    // Stable monotone interpolation even for nearly repeated modes.
                    let x = a + fraction * (z - a);
                    let f = if k == 0 {
                        pair[0]
                    } else if k + 1 == b.points {
                        pair[1]
                    } else if b.logarithmic {
                        libm::exp(x)
                    } else {
                        x
                    };
                    let f = f.clamp(pair[0], pair[1]);
                    if out.last().is_some_and(|last| *last == f) {
                        continue;
                    }
                    if out.len() >= crate::harmonic::MAX_FREQUENCIES {
                        return Err(Error::Invalid("random sampling exceeds frequency limit".into()));
                    }
                    out.push(f);
                }
            }
        }
        crate::harmonic::Frequencies::List(out).points()
    }
}

#[cfg(test)]
mod sampling_tests {
    use super::*;
    #[test]
    fn nearly_repeated_modes_have_strictly_ascending_representable_samples() {
        let modal = ModalRandom {
            source_step: "Modes".into(),
            modes: 2,
            damping_ratios: vec![0.02; 2],
            bands: vec![ModalBand { lower: 0., upper: 1000., points: 200, bias: 3., logarithmic: false }],
        };
        let frequencies = modal.frequencies(&[100., 100. + 1e-12]).unwrap();
        assert!(frequencies.windows(2).all(|w| w[1] > w[0]));
        assert!(frequencies.contains(&100.));
        assert!(frequencies.contains(&(100. + 1e-12)));
        assert_eq!(frequencies.first(), Some(&0.));
        assert_eq!(frequencies.last(), Some(&1000.));
    }
}
