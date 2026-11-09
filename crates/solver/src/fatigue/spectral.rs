//! Fixed linear stress projection and stationary Gaussian narrow-band damage.
//! Dirlik & Benasciutti (2021), Metals 11,1333, doi:10.3390/met11091333,
//! equations (1),(2),(5). Here moments use Hz, so crossing rates are in 1/s
//! WITHOUT a 2*pi divisor. PSD is one-sided. No multiaxial invariant PSD is inferred.
use super::{Basquin, life::StressMean};
use crate::{Error, HarmonicResult, RandomInput, Result};

const PROBABILITY_ROUNDOFF: f64 = 1e-12; // unit-factors: dimensionless probability/moment roundoff, not a conversion

/// Project six coherent stress phasors BEFORE squaring. Components are
/// xx,yy,zz,xy,xz,yz. Weights are dimensionless and fixed over frequency.
/// A plane normal n uses [nx²,ny²,nz²,2nxny,2nxnz,2nynz].
pub fn project(real: &[f64; 6], imaginary: &[f64; 6], weights: &[f64; 6], input_psd: f64) -> Result<f64> {
    if real.iter().chain(imaginary).chain(weights).any(|x| !x.is_finite())
        || weights.iter().all(|x| *x == 0.0)
        || !input_psd.is_finite()
        || input_psd < 0.0
    {
        return Err(Error::Invalid("invalid coherent stress projection".into()));
    }
    let r: f64 = real.iter().zip(weights).map(|(x, w)| x * w).sum();
    let i: f64 = imaginary.iter().zip(weights).map(|(x, w)| x * w).sum();
    let value = (r * r + i * i) * input_psd;
    if !value.is_finite() {
        return Err(Error::Invalid("projected stress PSD overflow".into()));
    }
    Ok(value)
}

/// Adapter to actual stored harmonic results. input_psd is already interpolated
/// on exactly these frames, in squared SI excitation units / Hz. The result must
/// represent one unit spatial excitation, or one explicitly coherent spatial
/// pattern driven by a COMMON scalar stochastic input. No independence inferred.
pub fn from_harmonic(h: &HarmonicResult, node: u32, weights: &[f64; 6], input_psd: &[f64]) -> Result<Vec<RandomInput>> {
    if h.frames.len() < 2 || h.frames.len() > 10000 || input_psd.len() != h.frames.len() {
        return Err(Error::Invalid("stress projection needs 2..10000 matching harmonic frames and input PSD values".into()));
    }
    let index = h.nodes.iter().position(|n| *n == node).ok_or_else(|| Error::Invalid("stress projection node absent".into()))?;
    let mut out = Vec::with_capacity(h.frames.len());
    for (frame, &g) in h.frames.iter().zip(input_psd) {
        let real = frame.stress_real.as_ref().and_then(|s| s.get(index)).ok_or_else(|| Error::Invalid("missing harmonic real stress".into()))?;
        let imaginary =
            frame.stress_imaginary.as_ref().and_then(|s| s.get(index)).ok_or_else(|| Error::Invalid("missing harmonic imaginary stress".into()))?;
        out.push(RandomInput { frequency_hz: frame.frequency, psd: project(real, imaginary, weights, g)? });
    }
    moments(&out)?;
    Ok(out)
}

/// m0,m1,m2,m4 for a piecewise-linear PSD, zero outside its stated band.
/// Three-point Gauss integrates f^j times each linear interval exactly (degree<=5).
/// Keeping nonnegative quadrature weights preserves moment inequalities.
pub fn moments(samples: &[RandomInput]) -> Result<[f64; 4]> {
    if samples.len() < 2
        || samples.len() > 10000
        || samples.iter().any(|s| !s.frequency_hz.is_finite() || s.frequency_hz < 0.0 || !s.psd.is_finite() || s.psd < 0.0)
        || samples.windows(2).any(|w| w[1].frequency_hz <= w[0].frequency_hz)
    {
        return Err(Error::Invalid("invalid one-sided stress spectrum".into()));
    }
    let mut m = [0.0; 4];
    let q = libm::sqrt(3.0 / 5.0);
    for pair in samples.windows(2) {
        if let [a, b] = pair {
            let width = b.frequency_hz - a.frequency_hz;
            for (x, w) in [(-q, 5.0 / 9.0), (0.0, 8.0 / 9.0), (q, 5.0 / 9.0)] {
                let t = (1.0 + x) * 0.5;
                let f = a.frequency_hz + t * width;
                let g = (1.0 - t) * a.psd + t * b.psd;
                let mass = (width * 0.5) * w * g;
                let f2 = f * f;
                for (dest, power) in m.iter_mut().zip([1.0, f, f2, f2 * f2]) {
                    *dest += mass * power;
                }
            }
        }
    }
    if m.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("spectral moment overflow".into()));
    }
    Ok(m)
}

#[derive(Debug, Clone, Copy)]
pub struct NarrowBand {
    pub moments_hz: [f64; 4],
    pub zero_crossing_rate: f64,
    pub irregularity: f64,
    pub damage: f64,
}

/// Rayleigh narrow-band approximation, not a broad-band rainflow predictor.
/// Caller must establish stationarity, Gaussianity, linear response, coverage,
/// and resolution. Returns bandwidth diagnostic, without silently certifying it.
/// Only a single Basquin slope with no endurance cutoff is supported here.
pub fn narrow_band(samples: &[RandomInput], curve: &Basquin, seconds: f64, mean: f64, correction: StressMean) -> Result<NarrowBand> {
    curve.validate()?;
    if curve.endurance_limit != 0.0 || !seconds.is_finite() || seconds <= 0.0 {
        return Err(Error::Invalid("narrow-band closed form needs no endurance cutoff and positive duration in seconds".into()));
    }
    let factor = correction.corrected_amplitude(1.0, mean)?;
    let m = moments(samples)?;
    let [m0, _, m2, m4] = m;
    if m0 == 0.0 {
        return Ok(NarrowBand { moments_hz: m, zero_crossing_rate: 0.0, irregularity: 0.0, damage: 0.0 });
    }
    if m2 <= 0.0 || m4 <= 0.0 {
        return Err(Error::Invalid("stress spectrum has no positive-frequency variance".into()));
    }
    let rate = libm::sqrt(m2 / m0);
    let peaks = libm::sqrt(m4 / m2);
    let irregularity = rate / peaks;
    let k = curve.exponent;
    let log_damage = libm::log(seconds) + libm::log(rate) - libm::log(curve.cycles_ref)
        + k * (0.5 * (libm::log(2.0) + libm::log(m0)) + libm::log(factor) - libm::log(curve.amplitude_ref))
        + libm::lgamma(1.0 + k * 0.5);
    let damage = libm::exp(log_damage);
    if !rate.is_finite()
        || rate <= 0.0
        || !irregularity.is_finite()
        || irregularity <= 0.0
        || irregularity > 1.0 + PROBABILITY_ROUNDOFF
        || !damage.is_finite()
        || damage == 0.0
    {
        return Err(Error::Invalid("narrow-band result outside representable domain".into()));
    }
    Ok(NarrowBand { moments_hz: m, zero_crossing_rate: rate, irregularity, damage })
}

/// Single-slope stationary Gaussian spectral approximations. TB is the 2005
/// four-moment w2 variant (review eqs.11,13), not the earlier w1 approximation.
#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Method {
    #[default]
    NarrowBand,
    Dirlik,
    BenasciuttiTovo,
}

/// Return D/D_NB. Coefficients follow Dirlik & Benasciutti (2021), eqs.8,9,11,13.
/// The exact/near-degenerate narrow-band limit is evaluated analytically rather
/// than dividing by vanishing coefficient denominators (threshold 1-alpha2<=1e-8).
pub fn broadband_ratio(m: [f64; 4], k: f64, method: Method) -> Result<f64> {
    let [m0, m1, m2, m4] = m;
    if m.iter().any(|x| !x.is_finite() || *x <= 0.0) || !k.is_finite() || k < 1.0 || k > 20.0 {
        return Err(Error::Invalid("broadband fatigue requires positive finite moments and 1<=Basquin exponent<=20".into()));
    }
    let a1 = (m1 / libm::sqrt(m0)) / libm::sqrt(m2);
    let a2 = libm::sqrt(m2 / m0) / libm::sqrt(m4 / m2);
    if !a1.is_finite() || !a2.is_finite() || a2 <= 0.0 || a1 > 1.0 + PROBABILITY_ROUNDOFF || a2 > a1 + PROBABILITY_ROUNDOFF {
        return Err(Error::Invalid("inconsistent spectral moment inequalities".into()));
    }
    if method == Method::NarrowBand || 1.0 - a2 <= 1e-8 {
        return Ok(1.0);
    }
    let ratio = match method {
        Method::NarrowBand => 1.0,
        Method::BenasciuttiTovo => {
            let delta = a1 - a2;
            // (1-a1)(1-a2) avoids cancellation in 1+a1*a2-a1-a2.
            let w = delta * (1.112 * (1.0 - a1) * (1.0 - a2) * libm::exp(2.11 * a2) + delta) / (1.0 - a2).powi(2);
            if !w.is_finite() || !(-PROBABILITY_ROUNDOFF..=1.0 + PROBABILITY_ROUNDOFF).contains(&w) {
                return Err(Error::Invalid("Benasciutti-Tovo weight outside probability domain".into()));
            }
            let w = w.clamp(0.0, 1.0);
            w + (1.0 - w) * libm::pow(a2, k - 1.0)
        }
        Method::Dirlik => {
            let xm = a1 * a2;
            let d1 = 2.0 * (xm - a2 * a2) / (1.0 + a2 * a2);
            let den = 1.0 - a2 - d1 + d1 * d1;
            if d1 <= 0.0 || den == 0.0 {
                return Err(Error::Invalid("degenerate Dirlik coefficients".into()));
            }
            let r = (a2 - xm - d1 * d1) / den;
            let d2 = den / (1.0 - r);
            let d3 = 1.0 - d1 - d2;
            // Eq.8 simplifies identically to Q=1.25*D1, avoiding cancellation.
            let q = 1.25 * d1;
            if [r, d1, d2, d3, q].iter().any(|x| !x.is_finite()) || d2 < -PROBABILITY_ROUNDOFF || d3 < -PROBABILITY_ROUNDOFF || q <= 0.0 {
                return Err(Error::Invalid("Dirlik mixture outside probability domain".into()));
            }
            let rayleigh_moment = libm::pow(2.0, k * 0.5) * libm::tgamma(1.0 + k * 0.5);
            (d1 * libm::pow(q, k) * libm::tgamma(1.0 + k) / rayleigh_moment + d2.max(0.0) * libm::pow(r.abs(), k) + d3.max(0.0)) / a2
        }
    };
    if !ratio.is_finite() || ratio <= 0.0 {
        return Err(Error::Invalid("invalid spectral damage ratio".into()));
    }
    Ok(ratio)
}

pub fn damage(samples: &[RandomInput], curve: &Basquin, seconds: f64, mean: f64, correction: StressMean, method: Method) -> Result<NarrowBand> {
    // Broad-band cycle means are distributed: a single global Goodman factor
    // is not a justified replacement for their joint amplitude/mean density.
    if method != Method::NarrowBand && (mean != 0.0 || !matches!(correction, StressMean::None)) {
        return Err(Error::NotSupported("broadband fatigue currently requires zero mean and correction=none".into()));
    }
    if method != Method::NarrowBand && !(1.0..=20.0).contains(&curve.exponent) {
        return Err(Error::Invalid("broadband fatigue exponent must be 1..=20".into()));
    }
    let mut result = narrow_band(samples, curve, seconds, mean, correction)?;
    if result.damage != 0.0 && method != Method::NarrowBand {
        result.damage *= broadband_ratio(result.moments_hz, curve.exponent, method)?;
        if !result.damage.is_finite() || result.damage <= 0.0 {
            return Err(Error::Invalid("broadband damage overflow or underflow".into()));
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn opposite_phases_cancel_before_square() {
        let w = [0.5, 0.5, 0.0, 0.0, 0.0, 0.0];
        assert_eq!(project(&[2.0, -2.0, 0.0, 0.0, 0.0, 0.0], &[0.0; 6], &w, 3.0).unwrap(), 0.0);
        assert_eq!(project(&[2.0, 2.0, 0.0, 0.0, 0.0, 0.0], &[0.0; 6], &w, 3.0).unwrap(), 12.0);
        // Both cases have identical six auto-spectra: phase is indispensable.
        assert!(project(&[0.0; 6], &[0.0; 6], &[0.0; 6], 1.0).is_err());
    }
    #[test]
    fn rectangular_spectrum_has_analytic_moments_and_k2_damage() {
        let s = [RandomInput { frequency_hz: 1.0, psd: 4.0 }, RandomInput { frequency_hz: 3.0, psd: 4.0 }];
        let m = moments(&s).unwrap();
        for (a, b) in m.into_iter().zip([8.0, 16.0, 104.0 / 3.0, 968.0 / 5.0]) {
            assert!((a / b - 1.0).abs() < 1e-12);
        }
        let curve = Basquin { amplitude_ref: 2.0, cycles_ref: 100.0, exponent: 2.0, endurance_limit: 0.0 };
        let r = narrow_band(&s, &curve, 10.0, 0.0, StressMean::None).unwrap();
        // Rayleigh E[a²]=2m0, so D=T*sqrt(13/3)*16/(100*4).
        assert!((r.damage / (0.4 * libm::sqrt(13.0 / 3.0)) - 1.0).abs() < 1e-12);
        assert!(r.irregularity < 1.0); // formula check, not a narrow-band validation
    }
}

#[cfg(test)]
mod broad_tests {
    use super::*;
    #[test]
    fn narrow_limit_zero_spectrum_and_invalid_moments() {
        for method in [Method::Dirlik, Method::BenasciuttiTovo] {
            assert_eq!(broadband_ratio([4.0, 8.0, 16.0, 64.0], 3.0, method).unwrap(), 1.0);
            assert!(broadband_ratio([1.0, 2.0, 1.0, 1.0], 3.0, method).is_err());
            let zero = [RandomInput { frequency_hz: 1.0, psd: 0.0 }, RandomInput { frequency_hz: 2.0, psd: 0.0 }];
            let curve = Basquin { amplitude_ref: 1.0, cycles_ref: 1.0, exponent: 3.0, endurance_limit: 0.0 };
            assert_eq!(damage(&zero, &curve, 10.0, 0.0, StressMean::None, method).unwrap().damage, 0.0);
            assert!(damage(&zero, &curve, 10.0, 1.0, StressMean::None, method).is_err());
        }
    }
    #[test]
    fn dirlik_closed_form_matches_independent_pdf_quadrature() {
        // Artificial realizable moments: alpha1=.6, alpha2=.5, xm=.3.
        // Eq.8 gives D1=2/25, R=242/533, D2=284089/363750, Q=1/10.
        let (d1, r, d2, q) = (2.0 / 25.0, 242.0 / 533.0, 284089.0 / 363750.0, 0.1);
        let d3 = 1.0 - d1 - d2;
        let k = 3.0;
        let mut integral = 0.0;
        let h = 20.0 / 20000.0;
        for j in 0..=20000 {
            let z = j as f64 * h;
            let pdf = d1 / q * libm::exp(-z / q) + d2 * z / (r * r) * libm::exp(-z * z / (2.0 * r * r)) + d3 * z * libm::exp(-z * z / 2.0);
            let weight = if j == 0 || j == 20000 {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * z.powi(3) * pdf;
        }
        integral *= h / 3.0;
        let expected = integral / (0.5 * libm::pow(2.0, k / 2.0) * libm::tgamma(1.0 + k / 2.0));
        assert!((broadband_ratio([1.0, 0.6, 1.0, 4.0], k, Method::Dirlik).unwrap() / expected - 1.0).abs() < 1e-8);
    }
    #[test]
    fn tb_2005_weight_and_scaling() {
        // alpha1=.6, alpha2=.5 -> w=.04+.08896 exp(1.055), k=3.
        let w = 0.04 + 0.08896 * libm::exp(1.055);
        let expected = w + (1.0 - w) * 0.25;
        let m = [1.0, 0.6, 1.0, 4.0];
        assert!((broadband_ratio(m, 3.0, Method::BenasciuttiTovo).unwrap() / expected - 1.0).abs() < 1e-12);
        for method in [Method::Dirlik, Method::BenasciuttiTovo] {
            let scale = m.map(|x| x * 9.0);
            assert!((broadband_ratio(scale, 3.0, method).unwrap() / broadband_ratio(m, 3.0, method).unwrap() - 1.0).abs() < 1e-12);
        }
    }
}

#[cfg(test)]
mod scale_tests {
    use super::*;
    #[test]
    fn amplitude_frequency_and_near_narrow_limits() {
        let curve = Basquin { amplitude_ref: 10.0, cycles_ref: 100.0, exponent: 3.0, endurance_limit: 0.0 };
        let s = [RandomInput { frequency_hz: 1.0, psd: 4.0 }, RandomInput { frequency_hz: 3.0, psd: 4.0 }];
        let louder = [RandomInput { frequency_hz: 1.0, psd: 16.0 }, RandomInput { frequency_hz: 3.0, psd: 16.0 }];
        let faster = [RandomInput { frequency_hz: 2.0, psd: 2.0 }, RandomInput { frequency_hz: 6.0, psd: 2.0 }];
        let thin = [RandomInput { frequency_hz: 10.0, psd: 1e6 }, RandomInput { frequency_hz: 10.000001, psd: 1e6 }];
        for method in [Method::NarrowBand, Method::Dirlik, Method::BenasciuttiTovo] {
            let d = damage(&s, &curve, 10.0, 0.0, StressMean::None, method).unwrap().damage;
            assert!((damage(&louder, &curve, 10.0, 0.0, StressMean::None, method).unwrap().damage / d - 8.0).abs() < 1e-10);
            assert!((damage(&faster, &curve, 10.0, 0.0, StressMean::None, method).unwrap().damage / d - 2.0).abs() < 1e-10);
            let nb = damage(&thin, &curve, 10.0, 0.0, StressMean::None, Method::NarrowBand).unwrap().damage;
            assert!((damage(&thin, &curve, 10.0, 0.0, StressMean::None, method).unwrap().damage / nb - 1.0).abs() < 1e-10);
        }
    }
}

#[cfg(test)]
mod published_benchmark {
    use super::*;
    /// Dirlik & Benasciutti 2021, Metals 11,1333, Table 1: band-limited
    /// rectangular spectra, f2/f1=5 and 20, slopes 3,5,7. Ratios are independent
    /// of PSD level and frequency scaling. 0.001 tolerance covers 3-digit printing.
    #[test]
    fn original_authors_table1_damage_ratios() {
        for (c, expected) in [(5., [0.980, 0.920, 0.867]), (20., [0.998, 0.932, 0.875])] {
            let psd = [RandomInput { frequency_hz: 1., psd: 1. }, RandomInput { frequency_hz: c, psd: 1. }];
            for (k, reference) in [3., 5., 7.].into_iter().zip(expected) {
                let curve = Basquin { amplitude_ref: 1., cycles_ref: 1., exponent: k, endurance_limit: 0. };
                let dk = damage(&psd, &curve, 1., 0., StressMean::None, Method::Dirlik).unwrap().damage;
                let tb = damage(&psd, &curve, 1., 0., StressMean::None, Method::BenasciuttiTovo).unwrap().damage;
                assert!((tb / dk - reference).abs() < 0.001, "c={c}, k={k}, ratio={}, reference={reference}", tb / dk);
            }
        }
    }
}
