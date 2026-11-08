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
    pub stress_psd: Option<Vec<f64>>,
    pub stress_rms: Option<f64>,
}

/// Stored field/probe contract for a completed random step. Values are ordered by model node or
/// element label; PSD curves remain frequency-resolved while RMS values are integrated scalars.
#[derive(Debug, Clone, PartialEq)]
pub struct RandomFieldOutput {
    pub frequencies_hz: Vec<f64>,
    pub node_displacement_psd: Vec<Vec<f64>>,
    pub node_acceleration_psd: Vec<Vec<f64>>,
    pub node_displacement_rms: Vec<f64>,
    pub node_acceleration_rms: Vec<f64>,
    pub element_stress_psd: Vec<Vec<f64>>,
    pub element_stress_rms: Vec<f64>,
    pub probe_indices: Vec<usize>,
}

/// Combine mass-normalized modal stress components with modal generalized-coordinate PSDs.
/// `modal_stress[mode][entity][component]` uses SI stress units; the returned PSD is stress²/Hz.
pub fn modal_stress_field_psd(modal_stress: &[Vec<[f64; 6]>], modal_coordinate_psd: &[f64]) -> Result<Vec<[f64; 6]>> {
    if modal_stress.is_empty() || modal_stress.len() != modal_coordinate_psd.len() {
        return Err(Error::Invalid("modal stress and coordinate PSD mode counts differ".into()));
    }
    let entities = modal_stress[0].len();
    if entities == 0 || modal_stress.iter().any(|m| m.len() != entities) || modal_coordinate_psd.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(Error::Invalid("modal stress field dimensions or PSD values are invalid".into()));
    }
    let mut out = vec![[0.0; 6]; entities];
    for (mode, field) in modal_stress.iter().enumerate() {
        for (entity, components) in field.iter().enumerate() {
            for c in 0..6 {
                if !components[c].is_finite() {
                    return Err(Error::Invalid("modal stress contains non-finite values".into()));
                }
                out[entity][c] += components[c] * components[c] * modal_coordinate_psd[mode];
            }
        }
    }
    Ok(out)
}

/// Generalized-force PSD for a fully correlated distributed force. `participation[mode]` is the
/// modal integral of the spatial load shape (force units per input-load unit).
pub fn distributed_force_modal_psd(participation: &[f64], load_psd: f64) -> Result<Vec<f64>> {
    if !load_psd.is_finite() || load_psd < 0.0 || participation.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("distributed force PSD inputs are invalid".into()));
    }
    Ok(participation.iter().map(|p| p * p * load_psd).collect())
}

/// Modal coordinate PSD from a scalar modal transfer magnitude and generalized-force PSD.
pub fn modal_coordinate_psd(transfer: &[f64], generalized_force_psd: &[f64]) -> Result<Vec<f64>> {
    if transfer.len() != generalized_force_psd.len()
        || transfer.iter().any(|v| !v.is_finite())
        || generalized_force_psd.iter().any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(Error::Invalid("modal transfer and generalized-force PSD lengths or values are invalid".into()));
    }
    Ok(transfer.iter().zip(generalized_force_psd).map(|(h, s)| h * h * s).collect())
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
        if i > 0 && f <= frequencies_hz[i - 1] {
            return Err(Error::Invalid("random response frequencies must be strictly increasing".into()));
        }
        let u = transfer[i] * transfer[i] * base_psd[i];
        displacement_psd.push(u);
        acceleration_psd.push(u * (2.0 * std::f64::consts::PI * f).powi(4));
    }
    let u = integrate_psd(
        &frequencies_hz.iter().zip(&displacement_psd).map(|(&frequency_hz, &psd)| RandomInput { frequency_hz, psd }).collect::<Vec<_>>(),
    )?;
    let a = integrate_psd(
        &frequencies_hz.iter().zip(&acceleration_psd).map(|(&frequency_hz, &psd)| RandomInput { frequency_hz, psd }).collect::<Vec<_>>(),
    )?;
    Ok(RandomResponse {
        frequencies_hz: frequencies_hz.to_vec(),
        displacement_psd,
        acceleration_psd,
        displacement_rms: u.sqrt(),
        acceleration_rms: a.sqrt(),
        stress_psd: None,
        stress_rms: None,
    })
}

/// Deterministic band-limited SDOF force-response variance by composite Simpson integration.
pub fn sdof_force_bandpass_variance(natural_frequency_hz: f64, damping_ratio: f64, force_psd: f64, low_hz: f64, high_hz: f64) -> Result<f64> {
    if !low_hz.is_finite() || !high_hz.is_finite() || low_hz < 0.0 || high_hz <= low_hz {
        return Err(Error::Invalid("invalid SDOF band limits".into()));
    }
    let n = 4096usize;
    let h = (high_hz - low_hz) / n as f64;
    let mut sum = 0.0;
    for i in 0..=n {
        let f = low_hz + i as f64 * h;
        let y = sdof_force_displacement_psd(f, natural_frequency_hz, damping_ratio, force_psd)?;
        sum += if i == 0 || i == n {
            y
        } else if i % 2 == 0 {
            2.0 * y
        } else {
            4.0 * y
        };
    }
    Ok(sum * h / 3.0)
}

/// Integrate independent modal stress transfer functions into one RMS value per mode/probe.
pub fn modal_stress_rms(frequencies_hz: &[f64], modal_stress_transfer: &[Vec<f64>], force_psd: &[f64]) -> Result<Vec<f64>> {
    if frequencies_hz.len() < 2 || force_psd.len() != frequencies_hz.len() || modal_stress_transfer.iter().any(|m| m.len() != frequencies_hz.len()) {
        return Err(Error::Invalid("modal stress arrays have inconsistent lengths".into()));
    }
    modal_stress_transfer
        .iter()
        .map(|transfer| {
            let samples = frequencies_hz
                .iter()
                .zip(transfer)
                .zip(force_psd)
                .map(|((&frequency_hz, &h), &psd)| RandomInput { frequency_hz, psd: h * h * psd })
                .collect::<Vec<_>>();
            Ok(integrate_psd(&samples)?.sqrt())
        })
        .collect()
}

/// Equivalent von-Mises RMS from six zero-mean stress-component RMS values. This is a second
/// moment indicator; without cross-spectral densities it is not a time-domain peak estimate.
pub fn von_mises_rms_from_components(rms: [f64; 6]) -> Result<f64> {
    if rms.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(Error::Invalid("stress RMS components must be finite and non-negative".into()));
    }
    let [xx, yy, zz, xy, xz, yz] = rms;
    Ok((0.5 * ((xx - yy).powi(2) + (yy - zz).powi(2) + (zz - xx).powi(2) + 6.0 * (xy * xy + xz * xz + yz * yz))).sqrt())
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

/// Absolute displacement PSD for a force-excited SDOF oscillator with unit mass.
pub fn sdof_force_displacement_psd(frequency_hz: f64, natural_frequency_hz: f64, damping_ratio: f64, force_psd: f64) -> Result<f64> {
    if !frequency_hz.is_finite()
        || frequency_hz < 0.0
        || !natural_frequency_hz.is_finite()
        || natural_frequency_hz <= 0.0
        || !damping_ratio.is_finite()
        || damping_ratio <= 0.0
        || !force_psd.is_finite()
        || force_psd < 0.0
    {
        return Err(Error::Invalid("invalid SDOF force random-vibration input".into()));
    }
    let w = 2.0 * std::f64::consts::PI * frequency_hz;
    let wn = 2.0 * std::f64::consts::PI * natural_frequency_hz;
    let den = (wn * wn - w * w).powi(2) + (2.0 * damping_ratio * wn * w).powi(2);
    Ok(force_psd / den)
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

    #[test]
    fn bandpass_variance_is_positive() {
        let v = sdof_force_bandpass_variance(10.0, 0.05, 2.0, 0.0, 100.0).unwrap();
        assert!(v.is_finite() && v > 0.0);
    }

    #[test]
    fn modal_stress_rms_integrates_modes() {
        let r = modal_stress_rms(&[0.0, 1.0], &[vec![1.0, 1.0], vec![2.0, 2.0]], &[2.0, 2.0]).unwrap();
        assert!((r[0] - 2.0f64.sqrt()).abs() < 1e-14);
        assert!((r[1] - 2.0 * 2.0f64.sqrt()).abs() < 1e-14);
    }

    #[test]
    fn von_mises_rms_is_a_second_moment_indicator() {
        let v = von_mises_rms_from_components([2.0, 0.0, 0.0, 0.0, 0.0, 0.0]).unwrap();
        assert!((v - 2.0).abs() < 1e-14);
    }

    #[test]
    fn modal_stress_field_psd_preserves_si_components_and_entity_indices() {
        let field = modal_stress_field_psd(&[vec![[2.0, 0.0, 0.0, 0.0, 0.0, 0.0]], vec![[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]], &[1.0, 4.0]).unwrap();
        assert_eq!(field, vec![[8.0, 0.0, 0.0, 0.0, 0.0, 0.0]]);
    }

    #[test]
    fn correlated_distributed_force_preserves_modal_participation_units() {
        let q = distributed_force_modal_psd(&[2.0, -3.0], 4.0).unwrap();
        assert_eq!(q, vec![16.0, 36.0]);
        let x = modal_coordinate_psd(&[0.5, 2.0], &q).unwrap();
        assert_eq!(x, vec![4.0, 144.0]);
    }
}
