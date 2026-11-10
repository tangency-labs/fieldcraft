//! Edge seeds: how many elements a mesher puts along a design edge, and how they are spaced.
//!
//! A seed names an edge of the geometry (`{body, edge}` with the edge's persistent name) and asks
//! for a number of divisions (`count`) or a target element length (`size`), optionally graded
//! geometrically toward one or both ends (`bias`). Its `constraint` says how firmly a mesher must
//! keep it: `exact` (the edge gets exactly that count, or the mesher refuses and names the edges
//! whose counts conflict), `min` (at least that count) or `free` (a goal the mesher may move to
//! satisfy its other constraints).
//!
//! **Spacing.** With `bias {ratio, toward}` the element lengths grow by a constant factor: the
//! largest element is `ratio` times the smallest. `toward: end` puts the small elements at the
//! edge's end, `start` at its start, `both` at both ends (largest in the middle). The edge's start
//! and end are those of its curve as the geometry gives it. [`Seed::spacing`] gives the node
//! positions as fractions of the arc length, which the meshers and the viewport both use.

use serde::{Deserialize, Serialize};

use crate::geometry::{Entity, GeometryRef};
use crate::{Error, Model, Result};

/// Most seeds a model holds.
pub const MAX_SEEDS: usize = 100_000;
/// Most divisions a seed asks for.
pub const MAX_SEED_COUNT: u32 = 100_000;
/// Largest bias ratio.
pub const MAX_BIAS_RATIO: f64 = 1e4;

/// How firmly a mesher keeps a seed's count.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SeedConstraint {
    /// Exactly the seeded count, or the mesher refuses.
    #[default]
    Exact,
    /// At least the seeded count.
    Min,
    /// A goal the mesher may change.
    Free,
}

impl SeedConstraint {
    pub fn name(self) -> &'static str {
        match self {
            SeedConstraint::Exact => "exact",
            SeedConstraint::Min => "min",
            SeedConstraint::Free => "free",
        }
    }
}

/// Which end of the edge gets the small elements.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BiasToward {
    Start,
    End,
    Both,
}

/// Geometric grading along the edge.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bias {
    /// The largest element over the smallest (at least 1).
    pub ratio: f64,
    pub toward: BiasToward,
}

/// The divisions a seed asks for.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Divisions {
    Count(u32),
    /// A target element length (model units); the count is the edge's length over it, rounded up.
    Size(f64),
}

/// A seed on one edge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "SeedRepr", into = "SeedRepr")]
pub struct Seed {
    pub edge: GeometryRef,
    pub divisions: Divisions,
    pub bias: Option<Bias>,
    pub constraint: SeedConstraint,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SeedRepr {
    edge: GeometryRef,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    count: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    size: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    bias: Option<Bias>,
    #[serde(default)]
    constraint: SeedConstraint,
}

impl TryFrom<SeedRepr> for Seed {
    type Error = Error;
    fn try_from(r: SeedRepr) -> Result<Seed> {
        let divisions = match (r.count, r.size) {
            (Some(n), None) => Divisions::Count(n),
            (None, Some(h)) => Divisions::Size(h),
            _ => return Err(Error::Invalid("a seed takes one of `count` or `size`".into())),
        };
        Seed::new(r.edge, divisions, r.bias, r.constraint)
    }
}

impl From<Seed> for SeedRepr {
    fn from(s: Seed) -> SeedRepr {
        let (count, size) = match s.divisions {
            Divisions::Count(n) => (Some(n), None),
            Divisions::Size(h) => (None, Some(h)),
        };
        SeedRepr { edge: s.edge, count, size, bias: s.bias, constraint: s.constraint }
    }
}

impl Seed {
    /// A checked seed: an edge reference, a count of 1 to [`MAX_SEED_COUNT`] or a positive size,
    /// and a bias ratio of 1 to [`MAX_BIAS_RATIO`].
    pub fn new(edge: GeometryRef, divisions: Divisions, bias: Option<Bias>, constraint: SeedConstraint) -> Result<Seed> {
        if !matches!(edge.entity, Entity::Edge(_)) {
            return Err(Error::Invalid(format!("a seed is on an edge, not the {}", edge.entity.kind())));
        }
        let what = format!("the seed on {edge}");
        match divisions {
            Divisions::Count(n) if !(1..=MAX_SEED_COUNT).contains(&n) => {
                return Err(Error::Invalid(format!("{what}: `count` must be 1 to {MAX_SEED_COUNT}, got {n}")));
            }
            Divisions::Size(h) if !(h > 0.0 && h.is_finite()) => {
                return Err(Error::Invalid(format!("{what}: `size` must be a positive number, got {h}")));
            }
            _ => {}
        }
        if let Some(b) = bias
            && !(b.ratio >= 1.0 && b.ratio <= MAX_BIAS_RATIO)
        {
            return Err(Error::Invalid(format!("{what}: the bias `ratio` must be 1 to {MAX_BIAS_RATIO}, got {}", b.ratio)));
        }
        Ok(Seed { edge, divisions, bias, constraint })
    }

    /// The count on an edge of `length`: the seeded count, or the length over the seeded size
    /// rounded up (at least 1, at most [`MAX_SEED_COUNT`]).
    pub fn count(&self, length: f64) -> u32 {
        match self.divisions {
            Divisions::Count(n) => n,
            Divisions::Size(h) => {
                let n = (length / h).ceil();
                if n.is_finite() { n.clamp(1.0, f64::from(MAX_SEED_COUNT)) as u32 } else { 1 }
            }
        }
    }

    /// The node positions of `n` divisions as fractions of the arc length from the edge's start:
    /// `n + 1` values from 0 to 1, graded by the bias.
    pub fn spacing(&self, n: u32) -> Vec<f64> {
        spacing(n, self.bias)
    }
}

/// `n + 1` node positions in [0, 1] for `n` divisions graded by `bias` (uniform without one).
pub fn spacing(n: u32, bias: Option<Bias>) -> Vec<f64> {
    let n = n.clamp(1, MAX_SEED_COUNT) as usize;
    // Element lengths, in order from the start. The growth factor is the portable `pow`, the
    // same on every platform (the system library's `powf` differs in the last bits between macOS
    // and Linux, #115); `powi` is repeated multiplication.
    let lengths: Vec<f64> = match bias {
        Some(b) if b.ratio > 1.0 && n > 1 => match b.toward {
            BiasToward::Start | BiasToward::End => {
                let q = libm::pow(b.ratio, 1.0 / (n - 1) as f64);
                let grow: Vec<f64> = (0..n).map(|i| q.powi(i as i32)).collect();
                if b.toward == BiasToward::Start { grow } else { grow.into_iter().rev().collect() }
            }
            BiasToward::Both => {
                // Small at both ends, growing to the middle: the ratio is reached over half the
                // elements on each side.
                let half = n.div_ceil(2);
                let q = if half > 1 { libm::pow(b.ratio, 1.0 / (half - 1) as f64) } else { 1.0 };
                (0..n).map(|i| q.powi(i.min(n - 1 - i) as i32)).collect()
            }
        },
        _ => vec![1.0; n],
    };
    let total: f64 = lengths.iter().sum();
    let mut out = Vec::with_capacity(n + 1);
    let mut at = 0.0;
    out.push(0.0);
    for l in lengths.iter().take(n - 1) {
        at += l / total;
        out.push(at);
    }
    out.push(1.0);
    out
}

impl Model {
    /// The edge seeds, sorted by edge.
    pub fn seeds(&self) -> &[Seed] {
        &self.seeds
    }

    /// The seed on `edge` (named exactly), if any.
    pub fn seed(&self, edge: &GeometryRef) -> Option<&Seed> {
        self.seeds.binary_search_by(|s| s.edge.cmp(edge)).ok().and_then(|i| self.seeds.get(i))
    }

    /// Seed an edge, replacing a seed it had.
    pub fn set_seed(&mut self, seed: Seed) -> Result<()> {
        let seed = Seed::new(seed.edge, seed.divisions, seed.bias, seed.constraint)?;
        match self.seeds.binary_search_by(|s| s.edge.cmp(&seed.edge)) {
            Ok(i) => {
                if let Some(slot) = self.seeds.get_mut(i) {
                    *slot = seed;
                }
            }
            Err(i) => {
                if self.seeds.len() >= MAX_SEEDS {
                    return Err(Error::Limit(format!("a model holds at most {MAX_SEEDS} edge seeds")));
                }
                self.seeds.insert(i, seed);
            }
        }
        Ok(())
    }

    /// Remove the seed on `edge`; returns whether there was one.
    pub fn remove_seed(&mut self, edge: &GeometryRef) -> bool {
        match self.seeds.binary_search_by(|s| s.edge.cmp(edge)) {
            Ok(i) => {
                self.seeds.remove(i);
                true
            }
            Err(_) => false,
        }
    }

    /// Remove every seed; returns how many there were.
    pub fn clear_seeds(&mut self) -> usize {
        std::mem::take(&mut self.seeds).len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge(n: &str) -> GeometryRef {
        GeometryRef::new("Body1", Entity::Edge(n.into())).unwrap()
    }

    #[test]
    fn spacing_is_uniform_or_geometric() {
        assert_eq!(spacing(4, None), vec![0.0, 0.25, 0.5, 0.75, 1.0]);
        // Ratio 8 over 4 elements: lengths 1, 2, 4, 8 (of 15) toward the start.
        let s = spacing(4, Some(Bias { ratio: 8.0, toward: BiasToward::Start }));
        let want = [0.0, 1.0 / 15.0, 3.0 / 15.0, 7.0 / 15.0, 1.0];
        assert!(s.iter().zip(want).all(|(a, b)| (a - b).abs() < 1e-12), "{s:?}");
        // Toward the end: the mirror image.
        let e = spacing(4, Some(Bias { ratio: 8.0, toward: BiasToward::End }));
        assert!(e.iter().zip(s.iter().rev()).all(|(a, b)| (a - (1.0 - b)).abs() < 1e-12), "{e:?}");
        // Both ends: symmetric, smallest at the ends.
        let b = spacing(5, Some(Bias { ratio: 4.0, toward: BiasToward::Both }));
        let l: Vec<f64> = b.windows(2).map(|w| w[1] - w[0]).collect();
        assert!((l[0] - l[4]).abs() < 1e-12 && (l[2] / l[0] - 4.0).abs() < 1e-9, "{l:?}");
        assert_eq!(spacing(1, Some(Bias { ratio: 5.0, toward: BiasToward::End })), vec![0.0, 1.0]);
    }

    /// Graded node positions are the same bit for bit on every platform (#115): the growth factor
    /// is the portable `pow`, never the system maths library's. The fingerprint is of every
    /// position's bits; it holds on macOS and Linux alike, and moves only when the spacing does.
    #[test]
    fn graded_spacing_is_the_same_on_every_platform() {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for (ratio, toward) in [(3.7, BiasToward::Start), (11.0, BiasToward::End), (5.3, BiasToward::Both)] {
            for n in [7, 13, 50] {
                for x in spacing(n, Some(Bias { ratio, toward })) {
                    h = (h ^ x.to_bits()).wrapping_mul(0x0100_0000_01b3);
                }
            }
        }
        assert_eq!(h, 0xaea6_85c2_3b0e_f37b, "{h:#018x}");
    }

    #[test]
    fn seeds_are_checked_replaced_and_round_trip() {
        let mut m = Model::default();
        let s = Seed::new(edge("A|B"), Divisions::Count(8), Some(Bias { ratio: 3.0, toward: BiasToward::End }), SeedConstraint::Exact).unwrap();
        m.set_seed(s.clone()).unwrap();
        m.set_seed(Seed::new(edge("A|C"), Divisions::Size(2.0), None, SeedConstraint::Min).unwrap()).unwrap();
        m.set_seed(Seed { divisions: Divisions::Count(4), ..s.clone() }).unwrap();
        assert_eq!(m.seeds().len(), 2);
        assert_eq!(m.seed(&edge("A|B")).unwrap().divisions, Divisions::Count(4));
        assert_eq!(m.seed(&edge("A|C")).unwrap().count(9.0), 5);

        let json = serde_json::to_value(m.seeds()).unwrap();
        assert_eq!(
            json[0],
            serde_json::json!({"edge": {"body": "Body1", "edge": "A|B"}, "count": 4, "bias": {"ratio": 3.0, "toward": "end"}, "constraint": "exact"})
        );
        let back: Vec<Seed> = serde_json::from_value(json).unwrap();
        assert_eq!(back, m.seeds());

        for bad in [
            serde_json::json!({"edge": {"body": "B", "edge": "e"}}),
            serde_json::json!({"edge": {"body": "B", "edge": "e"}, "count": 2, "size": 1.0}),
            serde_json::json!({"edge": {"body": "B", "edge": "e"}, "count": 0}),
            serde_json::json!({"edge": {"body": "B", "edge": "e"}, "count": 200000}),
            serde_json::json!({"edge": {"body": "B", "edge": "e"}, "size": -1.0}),
            serde_json::json!({"edge": {"body": "B", "face": "e"}, "count": 2}),
            serde_json::json!({"edge": {"body": "B", "edge": "e"}, "count": 2, "bias": {"ratio": 0.5, "toward": "end"}}),
            serde_json::json!({"edge": {"body": "B", "edge": "e"}, "count": 2, "constraint": "maybe"}),
        ] {
            assert!(serde_json::from_value::<Seed>(bad.clone()).is_err(), "{bad}");
        }
        assert!(m.remove_seed(&edge("A|B")));
        assert!(!m.remove_seed(&edge("A|B")));
        assert_eq!(m.clear_seeds(), 1);
    }
}
