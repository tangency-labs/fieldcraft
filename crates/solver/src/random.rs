//! Deterministic scalar random-vibration primitives.
//!
//! These are deliberately independent of the finite-element model.  The modal solver will use
//! them after projecting force or base-acceleration inputs onto the retained modes.  Frequencies
//! are in Hz, PSDs are one-sided (per Hz), and damping is the modal viscous ratio.

use crate::{Error, Result};

/// A one-sided PSD sample at a frequency in Hz.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RandomInput {
    pub frequency_hz: f64,
    pub psd: f64,
}

/// Result of Miles' narrow-band white-noise approximation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MilesResult {
    pub natural_frequency_hz: f64,
    pub quality_factor: f64,
    pub acceleration_psd: f64,
    pub acceleration_rms: f64,
}

/// Frequency-by-frequency scalar response and integrated RMS quantities.
#[derive(Debug, Clone, PartialEq)]
pub struct RandomResponse {
    pub frequencies_hz: Vec<f64>,
    pub displacement_psd: Vec<f64>,
    pub acceleration_psd: Vec<f64>,
    pub displacement_rms: f64,
    pub acceleration_rms: f64,
}

/// Combine a modal transfer magnitude with a one-sided base PSD and integrate displacement and
/// absolute-acceleration response PSDs. `transfer` is the dimensionless absolute displacement
/// transfer at each frequency; acceleration is obtained by `(2πf)^2`.
pub fn modal_base_response(frequencies_hz: &[f64], transfer: &[f64], base_psd: &[f64]) -> Result<RandomResponse> {
    if frequencies_hz.len() < 2 || frequencies_hz.len() != transfer.len() || transfer.len() != base_psd.len() {
        return Err(Error::Invalid("random response arrays must have equal length and at least two points".into()));
    }
    let mut displacement_psd = Vec::with_capacity(frequencies_hz.len());
    let mut acceleration_psd = Vec::with_capacity(frequencies_hz.len());
    for i in 0..frequencies_hz.len() {
        let f = frequencies_hz[i];
        if !f.is_finite() || f < 0.0 || !transfer[i].is_finite() || !base_psd[i].is_finite() || base_psd[i] < 0.0 {
            return Err(Error::Invalid("random response inputs must be finite and non-negative where applicable".into()));
        }
        if i > 0 && f <= frequencies_hz[i - 1] { return Err(Error::Invalid("random response frequencies must be strictly increasing".into())); }
        let u = transfer[i] * transfer[i] * base_psd[i];
        displacement_psd.push(u);
        acceleration_psd.push(u * (2.0 * std::f64::consts::PI * f).powi(4));
    }
    let u = integrate_psd(&frequencies_hz.iter().zip(&displacement_psd).map(|(&frequency_hz, &psd)| RandomInput { frequency_hz, psd }).collect::<Vec<_>>())?;
    let a = integrate_psd(&frequencies_hz.iter().zip(&acceleration_psd).map(|(&frequency_hz, &psd)| RandomInput { frequency_hz, psd }).collect::<Vec<_>>())?;
    Ok(RandomResponse { frequencies_hz: frequencies_hz.to_vec(), displacement_psd, acceleration_psd, displacement_rms: u.sqrt(), acceleration_rms: a.sqrt() })
}

/// Integrate a non-negative one-sided PSD by the trapezoidal rule.
pub fn integrate_psd(samples: &[RandomInput]) -> Result<f64> {
    if samples.len() < 2 {
        return Err(Error::Invalid("PSD integration needs at least two samples".into()));
    }
    let mut area = 0.0;
    for pair in samples.windows(2) {
        let a = pair[0];
        let b = pair[1];
        if !a.frequency_hz.is_finite() || !b.frequency_hz.is_finite() || b.frequency_hz <= a.frequency_hz {
            return Err(Error::Invalid("PSD frequencies must be finite and strictly increasing".into()));
        }
        if !a.psd.is_finite() || !b.psd.is_finite() || a.psd < 0.0 || b.psd < 0.0 {
            return Err(Error::Invalid("PSD values must be finite and non-negative".into()));
        }
        area += 0.5 * (b.frequency_hz - a.frequency_hz) * (a.psd + b.psd);
    }
    Ok(area)
}

/// Miles' white-noise estimate for a lightly damped SDOF under base acceleration.
///
/// `base_acceleration_psd` is the one-sided input PSD in acceleration²/Hz.  The formula is
/// `a_rms² = π/2 · f_n · Q · S_a(f_n)`, with `Q = 1/(2ζ)`.
pub fn miles_base_acceleration_rms(natural_frequency_hz: f64, damping_ratio: f64, base_acceleration_psd: f64) -> Result<MilesResult> {
    if !natural_frequency_hz.is_finite() || natural_frequency_hz <= 0.0 {
        return Err(Error::Invalid("natural frequency must be positive and finite".into()));
    }
    if !damping_ratio.is_finite() || damping_ratio <= 0.0 {
        return Err(Error::Invalid("damping ratio must be positive and finite".into()));
    }
    if !base_acceleration_psd.is_finite() || base_acceleration_psd < 0.0 {
        return Err(Error::Invalid("base acceleration PSD must be finite and non-negative".into()));
    }
    let q = 1.0 / (2.0 * damping_ratio);
    let variance = std::f64::consts::PI * 0.5 * natural_frequency_hz * q * base_acceleration_psd;
    Ok(MilesResult { natural_frequency_hz, quality_factor: q, acceleration_psd: base_acceleration_psd, acceleration_rms: variance.sqrt() })
}

/// Absolute acceleration transfer PSD for a base-excited SDOF oscillator.
pub fn sdof_base_acceleration_psd(frequency_hz: f64, natural_frequency_hz: f64, damping_ratio: f64, base_psd: f64) -> Result<f64> {
    if !frequency_hz.is_finite()
        || frequency_hz < 0.0
        || !natural_frequency_hz.is_finite()
        || natural_frequency_hz <= 0.0
        || !damping_ratio.is_finite()
        || damping_ratio <= 0.0
        || !base_psd.is_finite()
        || base_psd < 0.0
    {
        return Err(Error::Invalid("invalid SDOF random-vibration input".into()));
    }
    let r = frequency_hz / natural_frequency_hz;
    let denominator = (1.0 - r * r).powi(2) + (2.0 * damping_ratio * r).powi(2);
    Ok(base_psd * (1.0 + (2.0 * damping_ratio * r).powi(2) * r * r) / denominator)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_constant_psd() {
        let s = [RandomInput { frequency_hz: 0.0, psd: 4.0 }, RandomInput { frequency_hz: 10.0, psd: 4.0 }];
        assert_eq!(integrate_psd(&s).unwrap(), 40.0);
    }

    #[test]
    fn miles_formula_matches_reference() {
        let r = miles_base_acceleration_rms(10.0, 0.05, 2.0).unwrap();
        let expected = (std::f64::consts::PI * 0.5 * 10.0 * 10.0 * 2.0).sqrt();
        assert!((r.acceleration_rms - expected).abs() < 1e-14);
    }

    #[test]
    fn transfer_is_finite_at_resonance() {
        let p = sdof_base_acceleration_psd(10.0, 10.0, 0.05, 1.0).unwrap();
        assert!(p.is_finite() && p > 0.0);
    }

    #[test]
    fn modal_response_integrates_constant_transfer() {
        let r = modal_base_response(&[0.0, 1.0], &[1.0, 1.0], &[2.0, 2.0]).unwrap();
        assert!((r.displacement_rms - 2.0f64.sqrt()).abs() < 1e-14);
        assert!(r.acceleration_rms > 0.0);
    }
}
