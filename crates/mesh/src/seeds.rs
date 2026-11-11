//! Edge seeds as the meshers take them: a design edge's polyline with the count, bias and
//! constraint the model's seed gives it (`fieldcraft_model::seeds`).
//!
//! The meshers recover their own curves from the face-tagged boundary surface, so a seed is
//! matched to a mesher's curve by position: the curve's end points lie on the seed's end points
//! (either way round, which gives the orientation the bias is measured in) and each one's middle
//! lies on the other, within the tolerance of the two tessellations. The seed's node positions
//! are placed by arc length along the mesher's curve ([`place`]).

use fieldcraft_model::{Bias, SeedConstraint, seeds::spacing};

use crate::vec3::{P, dist, lerp, sub};
use crate::{Error, Result};

/// A seed on one design edge.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeSeed {
    /// The edge's name, for messages.
    pub name: String,
    /// The edge's polyline, start to end.
    pub points: Vec<[f64; 3]>,
    /// How far the polyline may lie from the exact edge (its tessellation's chordal tolerance).
    pub tolerance: f64,
    /// The divisions the seed asks for.
    pub count: usize,
    pub bias: Option<Bias>,
    pub constraint: SeedConstraint,
}

impl EdgeSeed {
    /// `n + 1` node positions as fractions of the arc length from the seed's start (reversed:
    /// from its end).
    pub fn fractions(&self, n: usize, reversed: bool) -> Vec<f64> {
        let f = spacing(u32::try_from(n).unwrap_or(u32::MAX), self.bias);
        if reversed { f.iter().rev().map(|x| 1.0 - x).collect() } else { f }
    }

    /// The seed in words: `F2:a|F2:b (8 exact)`.
    pub fn describe(&self) -> String {
        format!("{} ({} {})", self.name, self.count, self.constraint.name())
    }
}

fn length(poly: &[P]) -> f64 {
    poly.windows(2).map(|w| dist(w[0], w[1])).sum()
}

/// The point at fraction `t` of a polyline's arc length.
fn at_fraction(poly: &[P], t: f64) -> Option<P> {
    let total = length(poly);
    let target = t.clamp(0.0, 1.0) * total;
    let mut walked = 0.0;
    for w in poly.windows(2) {
        let l = dist(w[0], w[1]);
        if walked + l >= target && l > 0.0 {
            return Some(lerp(w[0], w[1], ((target - walked) / l).clamp(0.0, 1.0)));
        }
        walked += l;
    }
    poly.last().copied()
}

fn point_segment(p: P, a: P, b: P) -> f64 {
    let ab = sub(b, a);
    let l2 = crate::vec3::dot(ab, ab);
    let t = if l2 > 0.0 { (crate::vec3::dot(sub(p, a), ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    dist(p, lerp(a, b, t))
}

/// The distance from `p` to a polyline.
pub(crate) fn distance_to(poly: &[P], p: P) -> f64 {
    point_polyline(p, poly)
}

fn point_polyline(p: P, poly: &[P]) -> f64 {
    match poly {
        [] => f64::INFINITY,
        [a] => dist(p, *a),
        _ => poly.windows(2).map(|w| point_segment(p, w[0], w[1])).fold(f64::INFINITY, f64::min),
    }
}

/// Whether polyline `curve` is the seed's edge within `tol`: `Some(reversed)` when it runs from
/// the seed's end to its start.
pub(crate) fn matches(seed: &EdgeSeed, curve: &[P], tol: f64) -> Option<bool> {
    let tol = tol + seed.tolerance;
    let (Some(&a), Some(&b), Some(&s0), Some(&s1)) = (curve.first(), curve.last(), seed.points.first(), seed.points.last()) else {
        return None;
    };
    let forward = dist(a, s0) <= tol && dist(b, s1) <= tol;
    let backward = dist(a, s1) <= tol && dist(b, s0) <= tol;
    if !(forward || backward) {
        return None;
    }
    // A closed edge (both ends at one point) is oriented by its first quarter.
    let reversed = if forward && backward {
        let q = at_fraction(curve, 0.25)?;
        let (f, r) = (at_fraction(&seed.points, 0.25)?, at_fraction(&seed.points, 0.75)?);
        dist(q, r) < dist(q, f)
    } else {
        backward
    };
    let (cm, sm) = (at_fraction(curve, 0.5)?, at_fraction(&seed.points, 0.5)?);
    (point_polyline(cm, &seed.points) <= tol && point_polyline(sm, curve) <= tol).then_some(reversed)
}

/// The seed (index, reversed) of each curve, the first that matches; two seeds of different
/// counts on one curve are an error naming both.
pub(crate) fn match_curves(seeds: &[EdgeSeed], curves: &[Vec<P>], tol: f64) -> Result<Vec<Option<(usize, bool)>>> {
    let mut out = vec![None; curves.len()];
    if seeds.is_empty() {
        return Ok(out);
    }
    for (c, curve) in curves.iter().enumerate() {
        for (i, s) in seeds.iter().enumerate() {
            let Some(rev) = matches(s, curve, tol) else { continue };
            match out.get(c).copied().flatten() {
                None => {
                    if let Some(slot) = out.get_mut(c) {
                        *slot = Some((i, rev));
                    }
                }
                Some((j, _)) => {
                    if let Some(o) = seeds.get(j)
                        && (o.count, o.bias, o.constraint) != (s.count, s.bias, s.constraint)
                    {
                        return Err(Error::Invalid(format!("the edge seeds {} and {} are on the same edge", o.describe(), s.describe())));
                    }
                }
            }
        }
    }
    Ok(out)
}

/// Whether two seeds lie on one curve (each within the other's tolerance).
fn same_curve(a: &EdgeSeed, b: &EdgeSeed) -> bool {
    matches(a, &b.points, b.tolerance).is_some()
}

/// The seeds of `edges` (seeds on edges) and of `faces` (the seeds face seeds put on their
/// edges) the meshers take, matched by position as they match them: an edge two bodies share
/// (a partition's cut) has a name in each, so a seed on it in one body and a face seed on it in
/// the other meet only here. A seed on the edge wins over a face's; of face seeds on one curve,
/// the one with the most divisions wins (the first on a tie). Returns `edges` as they are, then
/// the face seeds kept, in order.
pub fn with_face_seeds(mut edges: Vec<EdgeSeed>, faces: Vec<EdgeSeed>) -> Vec<EdgeSeed> {
    let kept: Vec<bool> = faces
        .iter()
        .enumerate()
        .map(|(i, f)| {
            !edges.iter().any(|e| same_curve(f, e))
                && !faces.iter().enumerate().any(|(j, g)| j != i && (g.count > f.count || (g.count == f.count && j < i)) && same_curve(f, g))
        })
        .collect();
    edges.extend(faces.into_iter().zip(kept).filter_map(|(f, k)| k.then_some(f)));
    edges
}

/// The points at the given arc-length fractions along `poly` (its end points exact).
pub(crate) fn place(poly: &[P], fractions: &[f64]) -> Result<Vec<P>> {
    let (Some(&first), Some(&last)) = (poly.first(), poly.last()) else {
        return Err(Error::Invalid("empty polyline".into()));
    };
    let n = fractions.len();
    fractions
        .iter()
        .enumerate()
        .map(|(k, &t)| {
            if k == 0 {
                Ok(first)
            } else if k + 1 == n {
                Ok(last)
            } else {
                at_fraction(poly, t).ok_or_else(|| Error::Failed("seed spacing walked off the curve".into()))
            }
        })
        .collect()
}

/// A use of a curve: its index and whether it is traversed from its start.
pub(crate) type Use = (usize, bool);

/// The spacing of biased seeds carried to the curves interval assignment ties to them (#257).
///
/// `own[c]` is curve `c`'s seed (any identity `T`) and whether the curve runs against it; each
/// tie says two curve uses run side by side with the same count (opposite sides of a mapped face
/// or a sweep's linking face, both traversed the same way along it). A biased seed's spacing
/// spreads along the ties to every curve without a seed of its own, which it reaches through
/// unseeded curves only; the curve takes the seed, oriented as the tie gives, so its nodes are
/// placed as the seeded curve's are. Seeds are spread in curve order and the first one to reach a
/// curve keeps it. Returns the carried seed of each unseeded curve (`None` elsewhere).
pub(crate) fn carry<T: Clone>(own: &[Option<(T, bool)>], biased: impl Fn(&T) -> bool, ties: &[[Use; 2]]) -> Vec<Option<(T, bool)>> {
    let n = own.len();
    let mut next: Vec<Vec<(usize, bool, bool)>> = vec![Vec::new(); n];
    for &[(a, fa), (b, fb)] in ties {
        if a == b || a >= n || b >= n {
            continue;
        }
        // Traversed as tied, `a` runs forwards iff `fa`; `b`'s orientation follows from `a`'s.
        if let Some(l) = next.get_mut(a) {
            l.push((b, fa, fb));
        }
        if let Some(l) = next.get_mut(b) {
            l.push((a, fb, fa));
        }
    }
    let mut out: Vec<Option<(T, bool)>> = vec![None; n];
    for (start, seed) in own.iter().enumerate() {
        let Some((seed, rev)) = seed else { continue };
        if !biased(seed) {
            continue;
        }
        let mut stack = vec![(start, *rev)];
        while let Some((c, rev_c)) = stack.pop() {
            for &(d, fc, fd) in next.get(c).map(Vec::as_slice).unwrap_or(&[]) {
                if own.get(d).is_some_and(Option::is_some) || out.get(d).is_some_and(Option::is_some) {
                    continue;
                }
                // Traversed as tied, `c` runs against the seed iff it runs forwards and against
                // it, or backwards and with it; `d` runs the same way as `c`.
                let against = rev_c == fc;
                let rev_d = if against { fd } else { !fd };
                if let Some(slot) = out.get_mut(d) {
                    *slot = Some((seed.clone(), rev_d));
                }
                stack.push((d, rev_d));
            }
        }
    }
    out
}

/// The element size a seed asks for at each of its nodes (the mean length of the elements on
/// either side), for size fields.
pub fn node_sizes(seed: &EdgeSeed) -> Vec<(P, f64)> {
    let f = seed.fractions(seed.count, false);
    let total = length(&seed.points);
    let Ok(pts) = place(&seed.points, &f) else { return Vec::new() };
    pts.iter()
        .enumerate()
        .filter_map(|(k, p)| {
            let before = k.checked_sub(1).and_then(|j| Some(f.get(k)? - f.get(j)?));
            let after = f.get(k + 1).zip(f.get(k)).map(|(b, a)| b - a);
            let h = match (before, after) {
                (Some(x), Some(y)) => 0.5 * (x + y),
                (Some(x), None) | (None, Some(x)) => x,
                (None, None) => return None,
            } * total;
            (h > 0.0 && h.is_finite()).then_some((*p, h))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use fieldcraft_model::BiasToward;

    use super::*;

    fn seed(points: Vec<P>) -> EdgeSeed {
        EdgeSeed {
            name: "e".into(),
            points,
            tolerance: 0.0,
            count: 3,
            bias: Some(Bias { ratio: 4.0, toward: BiasToward::Start }),
            constraint: SeedConstraint::Exact,
        }
    }

    #[test]
    fn curves_match_either_way_round_and_carry_the_orientation() {
        let s = seed(vec![[0.0; 3], [5.0, 0.0, 0.0], [10.0, 0.0, 0.0]]);
        assert_eq!(matches(&s, &[[0.0; 3], [10.0, 0.0, 0.0]], 1e-6), Some(false));
        assert_eq!(matches(&s, &[[10.0, 0.0, 0.0], [0.0; 3]], 1e-6), Some(true));
        // Same ends, another path.
        assert_eq!(matches(&s, &[[0.0; 3], [5.0, 3.0, 0.0], [10.0, 0.0, 0.0]], 1e-6), None);
        assert_eq!(matches(&s, &[[0.0; 3], [9.0, 0.0, 0.0]], 1e-6), None);
        // Reversed fractions mirror the spacing.
        let f = s.fractions(3, false);
        let r = s.fractions(3, true);
        assert!(f.iter().zip(r.iter().rev()).all(|(a, b)| (a - (1.0 - b)).abs() < 1e-12));
        let pts = place(&[[0.0; 3], [10.0, 0.0, 0.0]], &f).unwrap();
        assert!((pts[1][0] - 10.0 * f[1]).abs() < 1e-12);
        // Sizes: the mean of the neighbouring elements.
        let h = node_sizes(&s);
        assert_eq!(h.len(), 4);
        assert!((h[0].1 - 10.0 * f[1]).abs() < 1e-9);
    }

    /// On a curve two bodies share, a seed on the edge wins over a face seed, and of two face
    /// seeds the one with more divisions; face seeds elsewhere are kept.
    #[test]
    fn face_seeds_yield_to_edge_seeds_and_finer_face_seeds_on_one_curve() {
        let named = |name: &str, points: Vec<P>, count| EdgeSeed { name: name.into(), count, bias: None, ..seed(points) };
        let cut = vec![[0.0; 3], [10.0, 0.0, 0.0]];
        let back = vec![[10.0, 0.0, 0.0], [0.0; 3]];
        let other = vec![[0.0; 3], [0.0, 10.0, 0.0]];
        let edges = vec![named("a:cut", cut.clone(), 4)];
        let faces = vec![named("b:cut", back.clone(), 9), named("b:other", other.clone(), 3), named("c:other", other.clone(), 6)];
        let names = |v: &[EdgeSeed]| v.iter().map(|s| s.name.clone()).collect::<Vec<_>>();
        assert_eq!(names(&with_face_seeds(edges, faces.clone())), ["a:cut", "c:other"]);
        // Without the edge seed, the face seed holds the cut; equal counts keep the first.
        let tie = vec![named("b:other", other.clone(), 6), named("c:other", other, 6)];
        assert_eq!(names(&with_face_seeds(Vec::new(), faces)), ["b:cut", "c:other"]);
        assert_eq!(names(&with_face_seeds(Vec::new(), tie)), ["b:other"]);
    }

    /// A biased seed's spacing carries along the ties to unseeded curves, oriented as each tie
    /// says, and stops at a curve with a seed of its own (#257).
    #[test]
    fn biased_spacing_carries_along_ties() {
        // 0 runs against its seed; 0 forwards beside 1 backwards; 1 forwards beside 2 forwards;
        // 3 has its own seed and 4 lies beyond it.
        let own = vec![Some((0, true)), None, None, Some((1, false)), None, None];
        let ties = [[(0, true), (1, false)], [(1, true), (2, true)], [(2, true), (3, true)], [(3, true), (4, true)]];
        let biased = [true, true];
        let c = carry(&own, |&i: &usize| biased[i], &ties);
        assert_eq!(c, vec![None, Some((0, false)), Some((0, false)), None, Some((1, false)), None]);
        // Unbiased seeds keep to their own curve.
        assert!(carry(&own, |_: &usize| false, &ties).iter().all(Option::is_none));
    }
}
