//! Bounded block-mean acceptance for LES studies (#796), not a physical-validation certificate.
//!
//! Student-t mean intervals: NIST/SEMATECH e-Handbook §1.3.5.2 (Snedecor & Cochran,
//! *Statistical Methods*, 1989), https://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm.
//! Critical values: §1.3.6.7.2, 0.975 column. Values are rounded upward; degrees of freedom
//! above 30 use 30 conservatively. Correlations use the fixed-mean sample autocorrelation
//! (NIST §1.3.5.12). The acceptance thresholds are the #785 study policy.
//!
//! Enough blocks for the lags tested: under independence a lag-k autocorrelation from N values
//! has a standard error near 1/√N (NIST §1.3.5.12's ±z/√N band, after Bartlett, J. R. Stat.
//! Soc. Suppl. 8:27–41, 1946), and Box & Jenkins (*Time Series Analysis: Forecasting and
//! Control*, 1976, §2.1.6) estimate lags only up to about N/4. Lags 1 to [`MAX_LAG`] on the
//! doubled-length series of N/2 values therefore need N/2 ≥ 4 · [`MAX_LAG`], so acceptance needs
//! at least [`MIN_BLOCKS`] = 8 · [`MAX_LAG`] base blocks. At 20 doubled values the |r| < 0.5 gate
//! sits 2.2 standard errors out; with fewer, independent data are rejected often.

use crate::{Error, Result};

/// Largest correlation lag checked, on both block lengths.
pub const MAX_LAG: usize = 5;
/// Fewest base blocks acceptance needs: the doubled series then has 4 · [`MAX_LAG`] values.
pub const MIN_BLOCKS: usize = 8 * MAX_LAG;

/// Maximum number of complete blocks retained by one stream.
pub const MAX_BLOCKS: usize = 65_536;
/// Maximum number of samples per block.
pub const MAX_BLOCK_LENGTH: usize = 1_000_000_000;

/// A bounded stream of equally spaced samples. Burn-in is a number of samples, discarded
/// before forming blocks. Incomplete blocks never enter any statistic.
#[derive(Debug, Clone)]
pub struct Blocks {
    len: usize,
    burn_in: u64,
    discarded: u64,
    filled: usize,
    sum: f64,
    means: Vec<f64>,
    exhausted: bool,
}

/// Sampling diagnostics for the largest even prefix of complete blocks.
/// Optional metrics are absent when their calculation is not finite.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Acceptance {
    /// Samples in each base block.
    pub block_samples: usize,
    /// Samples in the analysed even prefix.
    pub analysed_samples: u64,
    /// All complete blocks, including a possible unmatched last block.
    pub complete_blocks: usize,
    /// Base blocks in the analysed even prefix.
    pub analysed_blocks: usize,
    /// Adjacent block pairs in the same analysed window.
    pub doubled_blocks: usize,
    /// Samples in the incomplete trailing block, excluded from analysis.
    pub partial_samples: usize,
    /// Samples discarded during burn-in so far.
    pub discarded_burn_in: u64,
    /// Mean over the analysed window.
    pub mean: Option<f64>,
    /// Student-t 95% interval half-width in sample units.
    pub half_width_95: Option<f64>,
    /// Interval half-width divided by the absolute mean.
    pub relative_half_width: Option<f64>,
    /// Absolute difference of half-window means divided by the absolute mean.
    pub relative_half_window_drift: Option<f64>,
    /// Lags 1 through [`MAX_LAG`], all checked in absolute value.
    pub correlations: Vec<f64>,
    /// Lags 1 through [`MAX_LAG`] after combining adjacent base blocks.
    pub doubled_correlations: Vec<f64>,
    /// Relative interval half-width after combining adjacent base blocks.
    pub doubled_relative_half_width: Option<f64>,
    /// Absolute relative change in half-width when the block length doubles.
    pub doubled_half_width_change: Option<f64>,
    /// False includes insufficient evidence; reaching a run's cap never changes this to true.
    pub accepted: bool,
}

impl Blocks {
    /// Create a stream with `len` samples per block after `burn_in` samples.
    ///
    /// # Errors
    /// Returns [`Error::Invalid`] when `len` is zero or exceeds [`MAX_BLOCK_LENGTH`].
    pub fn new(len: usize, burn_in: u64) -> Result<Self> {
        if len == 0 || len > MAX_BLOCK_LENGTH {
            return Err(Error::Invalid(format!("block length must be in 1..={MAX_BLOCK_LENGTH}")));
        }
        Ok(Self { len, burn_in, discarded: 0, filled: 0, sum: 0.0, means: Vec::new(), exhausted: false })
    }

    /// Add one equally spaced sample, discarding it while burn-in remains.
    ///
    /// # Errors
    /// Returns [`Error::Invalid`] for a non-finite sample or magnitude above 1e150,
    /// or [`Error::Limit`] when the block cap or available allocation is exhausted.
    pub fn push(&mut self, value: f64) -> Result<()> {
        // This bound also keeps squared deviations and their bounded sums finite.
        if !value.is_finite() || value.abs() > 1e150 {
            return Err(Error::Invalid("statistical samples must be finite and at most 1e150 in magnitude".into()));
        }
        if self.discarded < self.burn_in {
            self.discarded += 1;
            return Ok(());
        }
        if self.means.len() >= MAX_BLOCKS {
            self.exhausted = true;
            return Err(Error::Limit("statistical block cap reached; acceptance remains inconclusive".into()));
        }
        if self.filled + 1 == self.len {
            self.means.try_reserve_exact(1).map_err(|_| Error::Limit("not enough memory for block means".into()))?;
        }
        self.sum += value / self.len as f64;
        self.filled += 1;
        if self.filled == self.len {
            self.means.push(self.sum);
            self.sum = 0.0;
            self.filled = 0;
        }
        Ok(())
    }

    /// All complete block means, including a possible unmatched last block.
    pub fn means(&self) -> &[f64] {
        &self.means
    }

    /// Descriptive mean over all complete blocks; may differ from the acceptance window.
    pub fn mean(&self) -> f64 {
        mean(&self.means)
    }
    /// Descriptive 95% interval half-width over all complete blocks.
    pub fn half_width(&self) -> f64 {
        half_width(&self.means)
    }
    /// Descriptive lag-one correlation over all complete blocks.
    pub fn lag1(&self) -> f64 {
        correlation(&self.means, 1)
    }

    /// Apply the acceptance gates to the largest even prefix, so both block lengths
    /// assess the same window. The descriptive accessors may include one additional block.
    pub fn acceptance(&self) -> Acceptance {
        self.assess(MIN_BLOCKS)
    }

    fn assess(&self, minimum: usize) -> Acceptance {
        let n = self.means.len() / 2 * 2;
        let base = self.means.get(..n).unwrap_or_default();
        let doubled: Vec<_> = base.as_chunks::<2>().0.iter().map(|[a, b]| (a + b) / 2.0).collect();
        let m = mean(base);
        let width = half_width(base);
        let long_width = half_width(&doubled);
        let correlations: Vec<_> = (1..=MAX_LAG).map(|lag| correlation(base, lag)).collect();
        let doubled_correlations: Vec<_> = (1..=MAX_LAG).map(|lag| correlation(&doubled, lag)).collect();
        let drift = (mean(base.get(n / 2..).unwrap_or_default()) - mean(base.get(..n / 2).unwrap_or_default())).abs() / m.abs();
        let relative = width / m.abs();
        let long_relative = long_width / m.abs();
        let change = if width == 0.0 && long_width == 0.0 { 0.0 } else { (long_width / width - 1.0).abs() };
        let accepted = !self.exhausted
            && n >= minimum
            && relative < 0.01
            && long_relative < 0.01
            && drift < 0.01
            && correlations.iter().chain(&doubled_correlations).all(|r| r.abs() < 0.5)
            && change < 0.2;
        let finite = |v: f64| v.is_finite().then_some(v);
        Acceptance {
            block_samples: self.len,
            analysed_samples: (n as u64).saturating_mul(self.len as u64),
            complete_blocks: self.means.len(),
            analysed_blocks: n,
            doubled_blocks: doubled.len(),
            partial_samples: self.filled,
            discarded_burn_in: self.discarded,
            mean: finite(m),
            half_width_95: finite(width),
            relative_half_width: finite(relative),
            relative_half_window_drift: finite(drift),
            correlations,
            doubled_correlations,
            doubled_relative_half_width: finite(long_relative),
            doubled_half_width_change: finite(change),
            accepted,
        }
    }
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let anchor = values.first().copied().unwrap_or(0.0);
    anchor + values.iter().map(|x| (x - anchor) / values.len() as f64).sum::<f64>()
}

fn half_width(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return f64::INFINITY;
    }
    let m = mean(values);
    let variance = values.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    critical(values.len() - 1) * (variance / values.len() as f64).sqrt()
}

fn correlation(values: &[f64], lag: usize) -> f64 {
    if values.len() <= lag {
        return 1.0;
    }
    let m = mean(values);
    let variance = values.iter().map(|x| (x - m).powi(2)).sum::<f64>();
    if variance == 0.0 {
        return 0.0;
    }
    values.iter().zip(values.iter().skip(lag)).map(|(a, b)| (a - m) * (b - m)).sum::<f64>() / variance
}

fn critical(df: usize) -> f64 {
    const T: [f64; 30] = [
        12.707, 4.304, 3.183, 2.777, 2.572, 2.448, 2.366, 2.307, 2.263, 2.229, 2.202, 2.180, 2.161, 2.146, 2.132, 2.121, 2.111, 2.102, 2.094, 2.087,
        2.081, 2.075, 2.070, 2.065, 2.061, 2.057, 2.053, 2.049, 2.046, 2.043,
    ];
    T.get(df.clamp(1, 30) - 1).copied().unwrap_or(f64::INFINITY)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn blocks(values: impl IntoIterator<Item = f64>) -> Blocks {
        let mut b = Blocks::new(1, 0).unwrap();
        for v in values {
            b.push(v).unwrap();
        }
        b
    }
    #[test]
    fn burn_in_partial_blocks_and_student_interval() {
        let mut b = Blocks::new(4, 2).unwrap();
        for x in [100.0, 100.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 100.0] {
            b.push(x).unwrap();
        }
        assert_eq!(b.means(), &[1.0, 3.0]);
        assert_eq!(b.mean(), 2.0);
        assert!((b.half_width() - 12.707).abs() < 1e-12);
        assert_eq!(b.acceptance().partial_samples, 1);
        assert!(!b.acceptance().accepted);
        assert!(blocks([1.0; MIN_BLOCKS]).acceptance().accepted);
        assert!(!blocks([1.0; MIN_BLOCKS - 1]).acceptance().accepted);
        assert!(!blocks([0.0; MIN_BLOCKS]).acceptance().accepted);
    }
    #[test]
    fn alternating_lag_two_and_drift_are_rejected() {
        let alternating = blocks((0..100).map(|i| 1.0 + if i % 2 == 0 { 1e-4 } else { -1e-4 }));
        let r = alternating.acceptance();
        assert!(r.relative_half_width.unwrap() < 0.01);
        assert!(r.correlations[0] < -0.9 && !r.accepted);
        let lag_two = blocks((0..100).map(|i| 1.0 + if i % 4 < 2 { 1e-4 } else { -1e-4 }));
        let r = lag_two.acceptance();
        assert!(r.correlations[0].abs() < 0.5);
        assert!(r.correlations[1] < -0.9 && !r.accepted);
        let drift = blocks((0..100).map(|i| 1.0 + i as f64 * 0.0003));
        let r = drift.acceptance();
        assert!(r.relative_half_width.unwrap() < 0.01);
        assert!(r.relative_half_window_drift.unwrap() > 0.01 && !r.accepted);
    }
    #[test]
    fn independent_noise_passes_and_pair_dependence_changes_interval() {
        // Original deterministic pseudo-random synthetic signal, no solver involved.
        let mut seed = 57_u64;
        let mut values = Vec::new();
        for _ in 0..2000 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            values.push(1.0 + 0.001 * ((seed >> 32) as f64 / u32::MAX as f64 - 0.5));
        }
        let independent = blocks(values.clone()).acceptance();
        assert!(independent.accepted, "{independent:?}");
        // Opposite neighbours reduce the doubled interval although lag1 remains <0.5.
        let mut b = Blocks::new(1, 0).unwrap();
        for p in values.as_chunks::<2>().0 {
            b.push(p[0]).unwrap();
            b.push(1.0 - 0.6 * (p[0] - 1.0) + 0.8 * (p[1] - 1.0)).unwrap();
        }
        let r = b.acceptance();
        assert!(r.correlations.iter().all(|v| v.abs() < 0.5), "{r:?}");
        assert!(r.doubled_half_width_change.unwrap() > 0.2 && !r.accepted, "{r:?}");
    }
    /// The #828 review: at 20 base blocks the lag tests on 10 doubled values rejected
    /// independent data often. At the minimum, most independent windows must pass.
    #[test]
    fn independent_windows_at_the_minimum_usually_pass() {
        let mut seed = 91_u64;
        let mut uniform = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1_u64 << 53) as f64
        };
        let (windows, mut passed, mut old_minimum) = (400, 0, 0);
        for _ in 0..windows {
            // Approximately normal (sum of twelve uniforms), original synthetic data.
            let values: Vec<f64> = (0..MIN_BLOCKS).map(|_| 1.0 + 1e-3 * ((0..12).map(|_| uniform()).sum::<f64>() - 6.0)).collect();
            passed += usize::from(blocks(values.iter().copied()).acceptance().accepted);
            old_minimum += usize::from(blocks(values.iter().copied().take(20)).assess(20).accepted);
        }
        let (rate, old_rate) = (passed as f64 / windows as f64, old_minimum as f64 / windows as f64);
        println!("independent windows accepted: {rate} at {MIN_BLOCKS} blocks, {old_rate} at 20");
        assert!(rate > 0.8 && old_rate < 0.65, "accepted {rate} at {MIN_BLOCKS} blocks, {old_rate} at 20");
    }
    #[test]
    fn hostile_and_capped_streams_refuse_without_losing_blocks() {
        for len in [0, usize::MAX] {
            assert!(Blocks::new(len, 0).is_err());
        }
        let mut b = Blocks::new(1, 0).unwrap();
        for value in [f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(b.push(value).is_err());
        }
        for _ in 0..MAX_BLOCKS {
            b.push(1.0).unwrap();
        }
        assert!(b.push(1.0).is_err());
        assert_eq!(b.means().len(), MAX_BLOCKS);
        assert!(!b.acceptance().accepted);
    }
}
