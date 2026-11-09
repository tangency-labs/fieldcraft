//! Peak component estimates, not simultaneous signed fields or time histories.
//! CQC: Der Kiureghian (1981), response-spectrum modal correlation for stationary
//! broad-band excitation. Unequal viscous damping form as documented in the
//! public Abaqus Analysis manual, Response spectrum analysis / CQC.
//! Modal component signs must be retained before the quadratic combination.

use crate::{Error, Result};

/// CQC correlation. Frequencies may be Hz or rad/s, but must share units.
pub fn correlation(f_i: f64, z_i: f64, f_j: f64, z_j: f64) -> Result<f64> {
    if [f_i, f_j].iter().any(|f| !f.is_finite() || *f <= 0.0) || [z_i, z_j].iter().any(|z| !z.is_finite() || *z <= 0.0 || *z >= 1.0) {
        return Err(Error::Invalid("CQC needs positive finite frequencies and modal damping in (0,1)".into()));
    }
    // r <= 1 avoids powers of large frequency ratios and preserves reciprocity.
    let (r, a, b) = if f_i >= f_j { (f_j / f_i, z_i, z_j) } else { (f_i / f_j, z_j, z_i) };
    let r2 = r * r;
    let numerator = 8.0 * libm::sqrt(a * b) * (a + r * b) * r * libm::sqrt(r);
    let denominator = (1.0 - r2).powi(2) + 4.0 * a * b * r * (1.0 + r2) + 4.0 * (a * a + b * b) * r2;
    let rho = numerator / denominator;
    if !rho.is_finite() {
        return Err(Error::Invalid("CQC correlation underflow or overflow".into()));
    }
    Ok(rho.clamp(0.0, 1.0))
}

pub fn srss(responses: &[f64]) -> Result<f64> {
    if responses.is_empty() || responses.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("SRSS needs finite modal components".into()));
    }
    let result = responses.iter().fold(0.0, |sum, x| libm::hypot(sum, *x));
    if !result.is_finite() {
        return Err(Error::Invalid("SRSS result overflow".into()));
    }
    Ok(result)
}

/// Bounded quadratic work. Component values include participation and spectrum
/// ordinate; do not take their absolute values before calling this function.
pub fn cqc(responses: &[f64], frequencies: &[f64], damping: &[f64]) -> Result<f64> {
    let n = responses.len();
    if n == 0 || n > 4096 || frequencies.len() != n || damping.len() != n || responses.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("CQC needs 1..4096 matching modal components, frequencies and damping ratios".into()));
    }
    let scale = responses.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let mut sum = 0.0;
    let mut magnitude = 0.0;
    for i in 0..n {
        correlation(frequencies[i], damping[i], frequencies[i], damping[i])?;
        let x = if scale > 0.0 { responses[i] / scale } else { 0.0 };
        sum += x * x;
        magnitude += x * x;
        for j in 0..i {
            let rho = correlation(frequencies[i], damping[i], frequencies[j], damping[j])?;
            let y = if scale > 0.0 { responses[j] / scale } else { 0.0 };
            let term = 2.0 * rho * x * y;
            sum += term;
            magnitude += term.abs();
        }
    }
    if sum < -64.0 * f64::EPSILON * magnitude {
        return Err(Error::Invalid("CQC quadratic form is negative beyond roundoff".into()));
    }
    let peak = scale * libm::sqrt(sum.max(0.0));
    if !peak.is_finite() {
        return Err(Error::Invalid("CQC peak overflow".into()));
    }
    Ok(peak)
}

/// Envelope over all permutations and signs of the 100-40-40 direction rule.
pub fn directional_100_40_40(peaks: [f64; 3]) -> Result<f64> {
    if peaks.iter().any(|x| !x.is_finite() || *x < 0.0) {
        return Err(Error::Invalid("directional peaks must be finite and nonnegative".into()));
    }
    let largest = peaks.into_iter().fold(0.0, f64::max);
    let result = 0.4 * peaks.iter().sum::<f64>() + 0.6 * largest;
    if !result.is_finite() {
        return Err(Error::Invalid("direction combination overflow".into()));
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coincident_modes_keep_signed_correlation() -> Result<()> {
        assert!((correlation(1.0, 0.05, 1.0, 0.05)? - 1.0).abs() < 1e-14);
        assert!((cqc(&[3.0, 4.0], &[1.0, 1.0], &[0.05, 0.05])? - 7.0).abs() < 1e-13);
        assert!((cqc(&[3.0, -4.0], &[1.0, 1.0], &[0.05, 0.05])? - 1.0).abs() < 1e-13);
        assert_eq!(srss(&[3.0, -4.0])?, 5.0);
        Ok(())
    }
    #[test]
    fn reciprocal_unequal_damping_and_invalid_inputs() -> Result<()> {
        assert_eq!(correlation(10.0, 0.02, 11.0, 0.05)?, correlation(11.0, 0.05, 10.0, 0.02)?);
        assert!(correlation(1.0, 0.0, 1.0, 0.0).is_err());
        assert!(cqc(&[1.0], &[], &[0.05]).is_err());
        assert!((directional_100_40_40([10.0, 5.0, 2.0])? - 12.8).abs() < 1e-14);
        Ok(())
    }
    #[test]
    fn wilson_published_five_mode_correlation_table() -> Result<()> {
        // E. L. Wilson, Chapter 15 (2014 revision), Table 15.1, p.15-11:
        // https://edwilson.neocities.org/bookshelf/2005/2005%20Chapter%2015%20MODIFIED%202014.pdf
        // Published rounded reference numbers only; no vendor input/model copied.
        let omega = [13.87, 13.93, 43.99, 44.19, 54.42];
        let published = [
            [1.0, 0.998, 0.006, 0.006, 0.004],
            [0.998, 1.0, 0.006, 0.006, 0.004],
            [0.006, 0.006, 1.0, 0.998, 0.180],
            [0.006, 0.006, 0.998, 1.0, 0.186],
            [0.004, 0.004, 0.180, 0.186, 1.0],
        ];
        // Both inputs (0.01 rad/s) and outputs (0.001) are rounded in the
        // publication. Propagate the printed frequency intervals; comparing
        // only with half an output digit incorrectly rejects modes 3 and 5.
        let equal_damping = |a: f64, b: f64| {
            let r = a / b;
            0.02 * (1.0 + r) * r.powf(1.5) / ((1.0 - r * r).powi(2) + 0.01 * r * (1.0 + r).powi(2))
        };
        for i in 0..5 {
            for j in 0..5 {
                let actual = correlation(omega[i], 0.05, omega[j], 0.05)?;
                assert!((actual - equal_damping(omega[i], omega[j])).abs() < 1e-13);
                if i == j {
                    assert!((actual - published[i][j]).abs() < 1e-13);
                    continue;
                }
                let mut low = f64::INFINITY;
                let mut high = f64::NEG_INFINITY;
                for da in [-0.005, 0.005] {
                    for db in [-0.005, 0.005] {
                        let rho = equal_damping(omega[i] + da, omega[j] + db);
                        low = low.min(rho);
                        high = high.max(rho);
                    }
                }
                assert!(
                    low <= published[i][j] + 0.0005 && high >= published[i][j] - 0.0005,
                    "pair {i},{j}: rounded-input interval [{low:.12},{high:.12}], published {}",
                    published[i][j]
                );
            }
        }
        // Original signed probe components, checked against the rounded table.
        let components = [3.0, -2.0, 1.0, -0.5, 0.2];
        let mut expected = 0.0;
        for i in 0..5 {
            for j in 0..5 {
                expected += components[i] * published[i][j] * components[j];
            }
        }
        assert!((cqc(&components, &omega, &[0.05; 5])? - libm::sqrt(expected)).abs() < 0.005);
        Ok(())
    }
}
