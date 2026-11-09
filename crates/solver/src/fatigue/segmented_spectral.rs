//! Dirlik & Benasciutti (2021), Metals 11:1333, equations (4), (7), (8).
//! Integrate amplitude density / N(amplitude), splitting at every S-N knot and
//! endurance cutoff. No synthetic time history, range/mean pairing or TB extension.
use super::{
    curve::{Extrapolation, SnCurve, StressLife},
    life::StressMean,
    spectral::{self, Method, NarrowBand},
};
use crate::{Error, RandomInput, Result};

/// Log integral of t^(shape-1) exp(-t) over [lo,hi]. Substitute t=u² for a
/// smooth origin. Composite Simpson refinement is bounded; failure is an error.
/// For 1.5<=shape<=21, the omitted tail beyond mode+128+20 ln(1+lo)
/// is exponentially negligible relative to the retained positive integral.
fn log_gamma_interval(shape: f64, lo: f64, hi: f64) -> Result<f64> {
    if lo >= hi {
        return Ok(f64::NEG_INFINITY);
    }
    if !lo.is_finite() || lo > 1e8 {
        return Err(Error::Invalid("spectral S-N interval outside numerical domain".into()));
    }
    let end = hi.min(lo.max(shape - 1.) + 128. + 20. * lo.ln_1p());
    let a = lo.sqrt();
    let b = end.sqrt();
    let mode = (shape - 0.5).max(lo).min(end);
    let shift = (shape - 0.5) * mode.ln() - mode;
    let f = |u: f64| if u == 0. { 0. } else { ((2. * shape - 1.) * u.ln() - u * u - shift).exp() * 2. };
    let mut previous = 0.;
    for power in 5..=15 {
        let n = 1usize << power;
        let h = (b - a) / n as f64;
        let mut sum = f(a) + f(b);
        for i in 1..n {
            sum += (if i % 2 == 0 { 2. } else { 4. }) * f(a + i as f64 * h);
        }
        let integral = sum * h / 3.;
        if power > 5 && integral > 0. && (integral - previous).abs() <= 2e-9 * integral {
            return Ok(shift + integral.ln());
        }
        previous = integral;
    }
    Err(Error::Invalid("spectral S-N quadrature did not converge".into()))
}

pub fn damage(samples: &[RandomInput], curve: &SnCurve, seconds: f64, mean: f64, correction: StressMean, method: Method) -> Result<NarrowBand> {
    let SnCurve::Table(table) = curve else {
        return match curve {
            SnCurve::Basquin(b) => spectral::damage(samples, b, seconds, mean, correction, method),
            _ => Err(Error::Invalid("invalid spectral curve".into())),
        };
    };
    curve.validate()?;
    if method == Method::BenasciuttiTovo {
        return Err(Error::NotSupported(
            "benasciutti_tovo with segmented S-N: exponent-dependent correction has no implemented piecewise extension; use narrow_band or dirlik"
                .into(),
        ));
    }
    if matches!(table.extrapolation, Extrapolation::Reject) {
        return Err(Error::NotSupported("spectral segmented S-N with extrapolation=reject: unbounded amplitude density needs end_slopes".into()));
    }
    if method == Method::Dirlik && (mean != 0. || !matches!(correction, StressMean::None)) {
        return Err(Error::NotSupported("dirlik segmented S-N requires mean_stress_pa=0 and correction=none".into()));
    }
    if table.points.len() > 64 || !seconds.is_finite() || seconds <= 0. {
        return Err(Error::Invalid("spectral segmented S-N requires 2..64 points and positive duration".into()));
    }
    let mean_factor = correction.corrected_amplitude(1., mean)?;
    let strength = table.surface_factor * table.size_factor;
    let m = spectral::moments(samples)?;
    let [m0, m1, m2, m4] = m;
    // Validate all slopes even for an empty spectrum.
    let mut segments = Vec::new();
    for (i, pair) in table.points.windows(2).enumerate() {
        let [p, q] = pair else {
            return Err(Error::Invalid("missing S-N segment".into()));
        };
        let k = (q.cycles.ln() - p.cycles.ln()) / (p.amplitude.ln() - q.amplitude.ln());
        if !k.is_finite() || !(1. - 1e-10..=20. + 1e-10).contains(&k) {
            return Err(Error::NotSupported("spectral segmented S-N slopes must lie in 1..=20".into()));
        }
        let lo = if i + 2 == table.points.len() { 0. } else { q.amplitude };
        let hi = if i == 0 { f64::INFINITY } else { p.amplitude };
        segments.push((lo.max(table.endurance_limit), hi, k.clamp(1., 20.), p));
    }
    if m0 == 0. {
        return Ok(NarrowBand { moments_hz: m, zero_crossing_rate: 0., irregularity: 0., damage: 0. });
    }
    let rate = (m2 / m0).sqrt();
    let peaks = (m4 / m2).sqrt();
    let a2 = rate / peaks;
    if !rate.is_finite() || rate <= 0. || !a2.is_finite() || a2 <= 0. || a2 > 1. + 1e-10 {
        return Err(Error::Invalid("invalid spectral rates".into()));
    }
    // (weight, amplitude scale / sigma, exponential rather than Rayleigh)
    let (components, cycle_rate) = if method == Method::NarrowBand || 1. - a2 <= 1e-8 {
        (vec![(1., 1., false)], rate)
    } else {
        let a1 = (m1 / m0.sqrt()) / m2.sqrt();
        let xm = a1 * a2;
        let d1 = 2. * (xm - a2 * a2) / (1. + a2 * a2);
        let den = 1. - a2 - d1 + d1 * d1;
        let r = (a2 - xm - d1 * d1) / den;
        let d2 = den / (1. - r);
        let d3 = 1. - d1 - d2;
        let q = 1.25 * d1;
        if [d1, d2, d3, r, q].iter().any(|x| !x.is_finite()) || d1 <= 0. || d2 < 0. || d3 < 0. || q <= 0. || r == 0. {
            return Err(Error::Invalid("dirlik segmented S-N mixture outside probability domain".into()));
        }
        (vec![(d1, q, true), (d2, r.abs(), false), (d3, 1., false)], peaks)
    };
    let mut total = 0.;
    for (weight, scale, exponential) in components {
        if weight == 0. {
            continue;
        }
        let stress_scale = m0.sqrt() * scale * mean_factor / strength;
        if !stress_scale.is_finite() || stress_scale <= 0. {
            return Err(Error::Invalid("spectral S-N scale overflow".into()));
        }
        for &(lo, hi, k, p) in &segments {
            if lo >= hi {
                continue;
            }
            let (shape, lower, upper, log_scale) = if exponential {
                (k + 1., lo / stress_scale, hi / stress_scale, stress_scale.ln())
            } else {
                (1. + k / 2., (lo / stress_scale).powi(2) / 2., (hi / stress_scale).powi(2) / 2., stress_scale.ln() + 2f64.ln() / 2.)
            };
            let log_integral = log_gamma_interval(shape, lower, upper)?;
            total += (weight.ln() + k * (log_scale - p.amplitude.ln()) - p.cycles.ln() + log_integral).exp();
        }
    }
    let damage = seconds * cycle_rate * total;
    if !damage.is_finite() || damage <= 0. {
        return Err(Error::Invalid("spectral segmented S-N damage outside representable domain".into()));
    }
    Ok(NarrowBand { moments_hz: m, zero_crossing_rate: rate, irregularity: a2, damage })
}

#[cfg(test)]
mod tests {
    use super::super::{
        Basquin,
        curve::{Point, Table},
    };
    use super::*;
    fn curve() -> SnCurve {
        SnCurve::Table(Table {
            points: vec![Point { cycles: 100., amplitude: 10. }, Point { cycles: 100000., amplitude: 1. }],
            endurance_limit: 0.,
            surface_factor: 0.8,
            size_factor: 0.5,
            extrapolation: Extrapolation::EndSlopes,
        })
    }
    #[test]
    fn same_slope_matches_closed_forms_and_strength_direction() {
        let samples = [RandomInput { frequency_hz: 1., psd: 1. }, RandomInput { frequency_hz: 3., psd: 1. }];
        let reduced = Basquin { amplitude_ref: 4., cycles_ref: 100., exponent: 3., endurance_limit: 0. };
        for method in [Method::NarrowBand, Method::Dirlik] {
            let x = damage(&samples, &curve(), 100., 0., StressMean::None, method).unwrap();
            let exact = spectral::damage(&samples, &reduced, 100., 0., StressMean::None, method).unwrap();
            assert!((x.damage / exact.damage - 1.).abs() < 2e-8);
        }
    }
    #[test]
    fn unsupported_combinations_name_the_method_or_policy() {
        let samples = [RandomInput { frequency_hz: 1., psd: 1. }, RandomInput { frequency_hz: 2., psd: 1. }];
        assert!(damage(&samples, &curve(), 1., 0., StressMean::None, Method::BenasciuttiTovo).unwrap_err().to_string().contains("benasciutti_tovo"));
        assert!(damage(&samples, &curve(), 1., 1., StressMean::None, Method::Dirlik).unwrap_err().to_string().contains("dirlik"));
        let SnCurve::Table(mut t) = curve() else { panic!() };
        t.extrapolation = Extrapolation::Reject;
        assert!(
            damage(&samples, &SnCurve::Table(t), 1., 0., StressMean::None, Method::NarrowBand)
                .unwrap_err()
                .to_string()
                .contains("extrapolation=reject")
        );
    }
}

#[cfg(test)]
mod piecewise_reference {
    use super::super::curve::{Point, Table};
    use super::*;
    #[test]
    fn cross_knot_rayleigh_has_elementary_truncated_moments() {
        // Variance=1. N=1600/a² below a=1, N=1600/a^4 above a=1.
        // Integrating a^3 exp(-a²/2) and a^5 exp(-a²/2) by parts gives
        // E[1/N]=(2+10 exp(-1/2))/1600; cutoff .5 removes low-amplitude damage.
        let samples = [RandomInput { frequency_hz: 1., psd: 1. }, RandomInput { frequency_hz: 2., psd: 1. }];
        for limit in [0., 0.5] {
            let curve = SnCurve::Table(Table {
                points: vec![Point { cycles: 100., amplitude: 2. }, Point { cycles: 1600., amplitude: 1. }, Point { cycles: 6400., amplitude: 0.5 }],
                endurance_limit: limit,
                surface_factor: 1.,
                size_factor: 1.,
                extrapolation: Extrapolation::EndSlopes,
            });
            let r = damage(&samples, &curve, 10., 0., StressMean::None, Method::NarrowBand).unwrap();
            let low = if limit == 0. { 2. } else { 2.25 * (-0.125f64).exp() };
            let expected = 10. * (7f64 / 3.).sqrt() * (low + 10. * (-0.5f64).exp()) / 1600.;
            assert!((r.damage / expected - 1.).abs() < 2e-8);
        }
    }
}
