//! Choosing among near-equal values the same way on every platform (#556).
//!
//! A decision that takes the smallest of some computed values (a distance, a cut's length, a
//! step's cost) reads their last bits when two of them are equal in exact arithmetic, as mirror
//! geometry makes them, and those bits differ between platforms' `sin`/`cos` and between
//! tessellations of the same part. Here values within a tolerance of each other form one group,
//! ordered by index inside it. The groups come from the sorted values (a gap above the
//! tolerance starts the next one), so the order is total; a comparison that calls values within
//! a tolerance equal is not, and a sort may refuse it.
//!
//! The tolerance follows what is compared: lengths take the merge tolerance
//! (`decompose::Settings::tol`), angles `region::ANGLE_TOLERANCE`, division counts a relative
//! `interval::GOAL_TIE`.

/// The indices of `keys` ordered by key, keys within `tol` of the one before them in that order
/// forming one group, kept in index order.
pub(crate) fn order(keys: &[f64], tol: f64) -> Vec<usize> {
    let key = |i: usize| keys.get(i).copied().unwrap_or(f64::INFINITY);
    let mut by_key: Vec<usize> = (0..keys.len()).collect();
    by_key.sort_by(|&a, &b| key(a).total_cmp(&key(b)).then(a.cmp(&b)));
    let mut grouped: Vec<(usize, usize)> = Vec::with_capacity(by_key.len());
    let (mut group, mut prev) = (0usize, None::<f64>);
    for i in by_key {
        let k = key(i);
        if prev.is_some_and(|p| k - p > tol) {
            group += 1;
        }
        prev = Some(k);
        grouped.push((group, i));
    }
    grouped.sort_unstable();
    grouped.into_iter().map(|(_, i)| i).collect()
}

/// The index of the smallest of `keys`: the first (lowest index) of those within `tol` of it.
pub(crate) fn min(keys: &[f64], tol: f64) -> Option<usize> {
    order(keys, tol).first().copied()
}

/// `v` moved by `ulps` ulps, up when positive: what another platform's rounding gives.
#[cfg(test)]
pub(crate) fn nudge(v: f64, ulps: i64) -> f64 {
    (0..ulps.unsigned_abs()).fold(v, |x, _| if ulps > 0 { x.next_up() } else { x.next_down() })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values equal but for their last bits keep their index order whichever way the bits fall;
    /// values clearly apart are ordered by value.
    #[test]
    fn near_equal_values_keep_their_index_order() {
        for (a, b) in [(0, 0), (1, -1), (-1, 1), (4, -3)] {
            let keys = [3.0, nudge(1.0, a), nudge(1.0, b), 2.0];
            assert_eq!(order(&keys, 1e-9), vec![1, 2, 3, 0], "{keys:?}");
            assert_eq!(min(&keys, 1e-9), Some(1));
        }
        assert_eq!(min(&[], 1e-9), None);
        assert_eq!(order(&[f64::NAN, 1.0], 1e-9).len(), 2);
    }
}
