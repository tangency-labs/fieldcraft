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
    Force { load: String, spectrum: Vec<SpectrumPoint> },
    BaseAcceleration { dof: u8, spectrum: Vec<SpectrumPoint> },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RandomVibration {
    pub frequencies: Vec<f64>,
    pub excitations: Vec<Excitation>,
    #[serde(default)]
    pub damping: crate::harmonic::Damping,
}

impl RandomVibration {
    pub fn validate(&self) -> Result<()> {
        if self.frequencies.len() < 2 || self.frequencies.windows(2).any(|w| !w[0].is_finite() || w[0] < 0.0 || w[1] <= w[0]) {
            return Err(Error::Invalid("random-vibration frequencies must be finite, non-negative and strictly increasing".into()));
        }
        if self.excitations.is_empty() {
            return Err(Error::Invalid("random-vibration needs at least one excitation".into()));
        }
        for excitation in &self.excitations {
            let spectrum = match excitation {
                Excitation::Force { spectrum, .. } | Excitation::BaseAcceleration { spectrum, .. } => spectrum,
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
