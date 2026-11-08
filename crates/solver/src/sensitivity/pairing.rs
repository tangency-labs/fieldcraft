//! Sensor-space real MAC and one-to-one maximum-weight mode assignment.
//! Allemang, Sound and Vibration 37(8), 14–23 (2003), MAC definition;
//! Kuhn, Naval Research Logistics Quarterly 2, 83–97 (1955), assignment.
//! Real modes only: sensor coordinates/directions and relative units must agree.
//! MAC does not resolve repeated eigenspaces or establish identifiability.
use crate::{Error, Result};
use serde::Serialize;

const MAX_MODES: usize = 64;
const MAX_SENSORS: usize = 4096;
const FORBIDDEN: f64 = 1e6; // unit-factors: dimensionless assignment penalty, above any feasible total

/// Pair order follows the measured vectors; indices into computed modes are zero-based.
#[derive(Debug, Serialize)]
pub struct ModePair {
    pub computed_mode: usize,
    pub mac: f64,
    /// Sign making the sensor dot product positive; not an amplitude calibration.
    pub sign: f64,
}

#[derive(Debug, Serialize)]
pub struct ModePairing {
    pub pairs: Vec<ModePair>,
    pub total_mac: f64,
    /// Best total MAC minus the best distinct feasible assignment; None if unique by dimensions.
    pub assignment_margin: Option<f64>,
}

fn unit_vector(x: &[f64]) -> Result<Vec<f64>> {
    if x.is_empty() || x.len() > MAX_SENSORS || x.iter().any(|v| !v.is_finite()) {
        return Err(Error::Invalid("MAC requires 1..=4096 finite sensor components".into()));
    }
    let scale = x.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Err(Error::Invalid("MAC cannot use a zero sensor vector".into()));
    }
    let norm = x.iter().map(|v| (v / scale).powi(2)).sum::<f64>().sqrt();
    Ok(x.iter().map(|v| (v / scale) / norm).collect())
}

/// Rectangular shortest augmenting-path Hungarian assignment, bounded by rows/columns.
/// All array indices below derive from validated dimensions or constructed assignments.
fn assign(cost: &[Vec<f64>], forbidden: Option<(usize, usize)>) -> Result<Option<(Vec<usize>, f64)>> {
    let n = cost.len();
    let m = cost.first().map_or(0, Vec::len);
    let mut u = vec![0.0; n + 1];
    let mut v = vec![0.0; m + 1];
    let mut p = vec![0; m + 1];
    let mut way = vec![0; m + 1];
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut minv = vec![f64::INFINITY; m + 1];
        let mut used = vec![false; m + 1];
        let mut found = false;
        for _ in 0..=m {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = f64::INFINITY;
            let mut j1 = 0;
            for j in 1..=m {
                if !used[j] {
                    let c = if forbidden == Some((i0 - 1, j - 1)) { FORBIDDEN } else { cost[i0 - 1][j - 1] };
                    let cur = c - u[i0] - v[j];
                    if cur < minv[j] {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if minv[j] < delta {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            if !delta.is_finite() || j1 == 0 {
                return Err(Error::Internal("MAC assignment has no augmenting column".into()));
            }
            for j in 0..=m {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                found = true;
                break;
            }
        }
        if !found {
            return Err(Error::NotConverged("MAC augmenting path cap".into()));
        }
        for _ in 0..=m {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
        if j0 != 0 {
            return Err(Error::NotConverged("MAC assignment reconstruction cap".into()));
        }
    }
    let mut indices = vec![0; n];
    for j in 1..=m {
        if p[j] != 0 {
            indices[p[j] - 1] = j - 1;
        }
    }
    if indices.iter().enumerate().any(|(i, j)| forbidden == Some((i, *j)) || cost[i][*j] == FORBIDDEN) {
        return Ok(None);
    }
    let total = indices.iter().enumerate().map(|(i, j)| 1.0 - cost[i][*j]).sum();
    Ok(Some((indices, total)))
}

/// Maximize the sum of MACs without reusing a computed mode. Reject weak pairs
/// and ambiguous *global* assignments. min_margin is an absolute total-MAC gap.
/// Thresholds are explicit caller policy, not universal statistical confidence.
/// Every candidate uses the identical ordered sensor set (no missing components).
pub fn pair_modes(measured: &[Vec<f64>], computed: &[Vec<f64>], min_mac: f64, min_margin: f64) -> Result<ModePairing> {
    if measured.is_empty()
        || measured.len() > computed.len()
        || computed.len() > MAX_MODES
        || !min_mac.is_finite()
        || !(0.0..=1.0).contains(&min_mac)
        || !min_margin.is_finite()
        || min_margin < 0.0
    {
        return Err(Error::Invalid("MAC pairing requires 1..=64 modes, enough candidates and finite nonnegative thresholds (MAC <=1)".into()));
    }
    let sensors = measured.first().map_or(0, Vec::len);
    if measured.iter().chain(computed).any(|x| x.len() != sensors) {
        return Err(Error::Invalid("MAC vectors must share the same sensor coordinates and dimensions".into()));
    }
    let measured = measured.iter().map(|x| unit_vector(x)).collect::<Result<Vec<_>>>()?;
    let computed = computed.iter().map(|x| unit_vector(x)).collect::<Result<Vec<_>>>()?;
    let dots: Vec<Vec<f64>> = measured.iter().map(|a| computed.iter().map(|b| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>()).collect()).collect();
    let cost: Vec<Vec<f64>> = dots
        .iter()
        .map(|row| {
            row.iter()
                .map(|d| {
                    let mac = (d * d).min(1.0);
                    if mac < min_mac { FORBIDDEN } else { 1.0 - mac }
                })
                .collect()
        })
        .collect();
    let (indices, total_mac) = assign(&cost, None)?.ok_or_else(|| Error::Invalid("no one-to-one modal assignment above minimum MAC".into()))?;
    let mut second: Option<f64> = None;
    // Any distinct assignment excludes at least one edge of the best assignment.
    for (i, j) in indices.iter().enumerate() {
        if let Some((_, total)) = assign(&cost, Some((i, *j)))? {
            second = Some(second.map_or(total, |old| old.max(total)));
        }
    }
    let margin = second.map(|x| (total_mac - x).max(0.0));
    let numerical_tie = 64.0 * f64::EPSILON * measured.len() as f64;
    if margin.is_some_and(|x| x <= min_margin.max(numerical_tie)) {
        return Err(Error::Invalid(format!(
            "ambiguous MAC assignment: total-MAC margin {margin:?}; add independent sensors or use a subspace observation"
        )));
    }
    let pairs = indices
        .iter()
        .enumerate()
        .map(|(i, j)| ModePair { computed_mode: *j, mac: 1.0 - cost[i][*j], sign: if dots[i][*j] < 0.0 { -1.0 } else { 1.0 } })
        .collect();
    Ok(ModePairing { pairs, total_mac, assignment_margin: margin })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn crossing_tracks_shapes_not_frequency_order() -> Result<()> {
        // M=I, K(p)=diag(1+p,1-p); the ordered modes exchange at p=0.
        // Both endpoints have simple, close eigenvalues; the exact crossing is excluded.
        let measured = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let computed = vec![vec![0.0, -2.0], vec![3.0, 0.0]];
        let r = pair_modes(&measured, &computed, 0.9, 0.01)?;
        assert_eq!(r.pairs.iter().map(|x| x.computed_mode).collect::<Vec<_>>(), vec![1, 0]);
        assert_eq!(r.pairs[1].sign, -1.0);
        let ordered_frequencies = [0.99_f64.sqrt(), 1.01_f64.sqrt()];
        assert_ne!(ordered_frequencies[0], ordered_frequencies[r.pairs[0].computed_mode]);
        Ok(())
    }
    #[test]
    fn global_assignment_beats_greedy_and_matches_enumeration() -> Result<()> {
        let cost = vec![vec![0.1, 0.2, 0.9], vec![0.11, 0.95, 0.8]];
        let Some((indices, total)) = assign(&cost, None)? else {
            return Err(Error::Internal("expected assignment".into()));
        };
        assert_eq!(indices, vec![1, 0]);
        for a in 0..3 {
            for b in 0..3 {
                if a != b {
                    assert!(total >= 2.0 - cost[0][a] - cost[1][b] - 1e-14);
                }
            }
        }
        Ok(())
    }
    #[test]
    fn rejects_ambiguous_incomplete_and_hostile_sensors() {
        let a = vec![vec![1.0, 0.0], vec![1.0, 0.0]];
        assert!(pair_modes(&a, &a, 0.0, 0.0).is_err());
        for bad in [vec![], vec![0.0, 0.0], vec![f64::NAN, 0.0], vec![f64::INFINITY, 0.0], vec![1.0]] {
            assert!(pair_modes(&[vec![1.0, 0.0]], &[bad], 0.0, 0.0).is_err());
        }
        assert!(pair_modes(&[vec![1.0, 0.0]], &[vec![0.0, 1.0]], 0.1, 0.0).is_err());
        assert!(pair_modes(&a, &a, f64::NAN, 0.0).is_err());
        assert!(pair_modes(&a, &a, 0.0, f64::NAN).is_err());
    }
    #[test]
    fn mac_survives_extreme_scaling_and_sign() -> Result<()> {
        let r = pair_modes(&[vec![1e300, -1e300]], &[vec![-1e-300, 1e-300]], 0.99, 0.0)?;
        assert!((r.pairs[0].mac - 1.0).abs() < 1e-14);
        assert_eq!(r.pairs[0].sign, -1.0);
        assert!(r.assignment_margin.is_none());
        Ok(())
    }
}
