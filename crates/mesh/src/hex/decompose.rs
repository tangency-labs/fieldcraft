//! Automatic decomposition of a swept cross-section into mappable blocks.
//!
//! The regions of a cross-section ([`super::arrangement`]) are meshed by structured schemes
//! ([`super::region`]): mapped four-sided regions, the triangle primitive for three-cornered
//! ones, O-grids for rings and disks. A region none of them fits (an L, a U round a pocket, a
//! plate with several holes) is cut here by straight lines into pieces that fit, before interval
//! assignment sees the cross-section, so the cuts are curves like any other and their counts are
//! assigned with the rest. Nothing of this reaches the model: the cuts exist only inside the
//! mesher.
//!
//! The method follows medial-axis subdivision (Tam and Armstrong, *2D finite element mesh
//! generation by medial axis subdivision*, Advances in Engineering Software 13(5–6):313–324,
//! 1991): a region is split at its concave features, each by a line from the feature into the
//! region, and where it branches, until every piece is a simple primitive (Li, Armstrong,
//! McKeag, *Quad mesh generation for k-sided faces and hex mesh generation for trivalent
//! polyhedra*, Finite Elements in Analysis and Design 26(4):279–301, 1997). Corners are
//! classified by their interior angle as in submapping (Whiteley, White, Benzley, Blacker, *Two
//! and three-quarter dimensional meshing facilitators*, Engineering with Computers 12:144–154,
//! 1996). Each filled region is decomposed on its own, piece by piece:
//! 1. **Holes apart.** A piece with several holes is split by lines across the gap between
//!    neighbouring holes (perpendicular to the line between their centres), so that the pieces
//!    hold one hole each.
//! 2. **O-grids.** A piece with one hole that is not a ring of two smooth loops (the
//!    [`Scheme::Annulus`]) is cut into pieces round the hole: from the outer loop's corners towards
//!    the hole's centre or along their bisectors, from the hole's concave corners outwards, or
//!    from the hole along four directions (a long plate's layout, whose corners the next step cuts
//!    off). The layout whose pieces score best is taken.
//! 3. **Pieces without holes** are split by a search over candidate cuts: lines from concave
//!    corners (along the bisector or the extensions of the corner's edges), from the ends and the
//!    middle of strongly turning curves (fillets) along their normal and square to the curves
//!    beside them, and from convex corners along the bisector, each to the first boundary it
//!    meets (or a curve end near there). The cut taken is the best of a two-ply look-ahead: the
//!    fewest free pieces (counting their concave features), then the least quality shortfall,
//!    then the fewest pieces. Pieces are cut until they are good or no cut helps.
//!
//! **Scoring.** A piece is scored by meshing it alone with its scheme on the target size (coarser
//! for a fine size: the quality of a block depends on its shape): the smallest corner scaled
//! Jacobian of its quadrilaterals, capped by their aspect (shortest edge against the longest or
//! the target size), so thin strips and pieces much smaller than the elements score low. A piece
//! no structured scheme fits counts as free. The score takes the place of the medial axis's geometric rules in Tam and Armstrong: it prefers cuts
//! that cross the region where it narrows, square to its sides, and that leave blocks that do
//! not shear.
//!
//! **Do no harm.** A region's blocks are kept when they are fair (score at least 0.5) or no worse
//! than the free all-quad mesh of the whole region (scored alike); otherwise the region is left
//! to the free mesher. A region whose pieces cannot be given division counts together (mapped
//! pieces round a hole whose opposite sides chain back to a longer side) keeps no cuts either. Cuts never end inside a seeded curve, so edge seeds keep
//! their curves.

use super::arrangement::{Arrangement, Q};
use super::interval::Intervals;
use super::region::{self, Scheme};
use crate::quad::QuadOptions;
use crate::quality::planar_corner;
use crate::size::SizeField;
use crate::{Error, Result};

/// A ray crosses an edge when their cross product is at least this, relative to the edge.
const MIN_SINE: f64 = 1e-12; // unit-factors: relative to the edge's length
/// A crossing within this of a segment's ends is on it.
const ON_SEGMENT: f64 = 1e-9; // unit-factors: a fraction of the segment
/// Angles this close are one.
const SAME_ANGLE: f64 = 1e-6; // unit-factors: an angle in radians

/// The tag of cut lines (no wall behind them).
pub(crate) const CUT: usize = usize::MAX;

/// Steps (one cut, or one hole step) per region at most.
const MAX_STEPS: usize = 256;
/// Blocks this good (score) are taken whatever the free mesh would do.
const FAIR: f64 = 0.5;
/// How much worse than the region's free mesh (score) its blocks may be.
const HARM_MARGIN: f64 = 0.05;

/// A piece whose score is at least this is not cut further.
const GOOD: f64 = 0.8;
/// The quality the shortfall is measured from.
const TARGET: f64 = 0.9;
/// The shortfall a cut must save to be worth its extra piece.
const MIN_GAIN: f64 = 2e-3;
/// Candidate cuts tried in one search at most, and over one region (keeps a pathological region
/// bounded).
const SEARCH_BUDGET: usize = 2_000;
const REGION_BUDGET: usize = 20_000;
/// The edge ratio (longest over shortest edge) a piece's quadrilaterals may have before it lowers
/// the piece's score: the score is at most this over the ratio, so a strip twice as long per
/// element scores 0.5 and one at the analyst's limit of 10 scores 0.4.
const ASPECT_FREE: f64 = 4.0;
/// Boundary divisions a piece is scored with at most.
const SCORE_RESOLUTION: f64 = 96.0;
/// Candidate cuts compared per piece at most.
const MAX_CANDIDATES: usize = 40;
/// Boundary points of a piece the search takes at most (larger regions are left as they are).
const MAX_PIECE_POINTS: usize = 4096;
/// The smallest angle (radians) a cut may make with the boundary at its ends (about 25°).
const MIN_CUT_ANGLE: f64 = 0.436;
/// Interior angles beyond this (radians, about 200°) are concave features.
const CONCAVE: f64 = 3.49;
/// Turning (radians) along one curve from which it is a rounded corner (about 40°).
const ROUNDED: f64 = 0.70;

/// What the cutting needs to know: the target size in the cross-section plane, the merge
/// tolerance and the face mesher options.
pub(crate) struct Settings<'a> {
    pub field: &'a dyn SizeField,
    pub tol: f64,
    pub opts: &'a QuadOptions,
}

impl Settings<'_> {
    fn size(&self, p: Q) -> f64 {
        let h = self.field.size([p[0], p[1], 0.0]);
        if h > 0.0 && h.is_finite() { h } else { f64::INFINITY }
    }
}

/// One loop of a region as a polygon: its points (counter-clockwise for the outer loop, clockwise
/// for holes, so the region is always on the left), which of them are curve ends, and which
/// edges (`k` → `k + 1`) lie on curves a cut may not end inside (seeded ones).
#[derive(Debug, Clone)]
struct Poly {
    pts: Vec<Q>,
    end: Vec<bool>,
    locked: Vec<bool>,
}

impl Poly {
    fn len(&self) -> usize {
        self.pts.len()
    }

    fn at(&self, i: usize) -> Q {
        self.pts.get(i % self.len().max(1)).copied().unwrap_or([0.0; 2])
    }

    fn is_end(&self, i: usize) -> bool {
        self.end.get(i % self.len().max(1)).copied().unwrap_or(false)
    }

    /// The interior angle at point `i` (the region on the left of the loop).
    fn angle(&self, i: usize) -> f64 {
        let n = self.len();
        interior(self.at(i + n - 1), self.at(i), self.at(i + 1))
    }

    /// The curves of the loop, each the points from one curve end to the next (closed), starting
    /// at the first end; a loop without ends is one closed curve. With the position where each
    /// curve starts.
    fn curves(&self) -> (Vec<Vec<Q>>, Vec<usize>) {
        let n = self.len();
        let first = (0..n).find(|&i| self.is_end(i)).unwrap_or(0);
        let mut curves = Vec::new();
        let mut starts = Vec::new();
        let mut cur = vec![self.at(first)];
        starts.push(first);
        for k in 1..=n {
            let i = (first + k) % n;
            cur.push(self.at(i));
            if self.is_end(i) || k == n {
                curves.push(std::mem::take(&mut cur));
                if k < n {
                    cur.push(self.at(i));
                    starts.push(i);
                }
            }
        }
        (curves, starts)
    }
}

fn sub(a: Q, b: Q) -> Q {
    [a[0] - b[0], a[1] - b[1]]
}

fn len(a: Q) -> f64 {
    libm::hypot(a[0], a[1])
}

fn cross(a: Q, b: Q) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot(a: Q, b: Q) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn unit(a: Q) -> Option<Q> {
    let l = len(a);
    (l > 0.0 && l.is_finite()).then(|| [a[0] / l, a[1] / l])
}

/// The interior angle at `o` between the incoming edge from `a` and the outgoing edge to `b`, the
/// region on the left: π on a straight line, less at a convex corner, more at a concave one.
fn interior(a: Q, o: Q, b: Q) -> f64 {
    let (u, v) = (sub(o, a), sub(b, o));
    std::f64::consts::PI - crate::vec3::atan2(cross(u, v), dot(u, v))
}

/// The inward bisector at a corner of interior angle `theta` whose outgoing edge runs along unit
/// `w`: `w` turned left by half the angle.
fn bisector(w: Q, theta: f64) -> Q {
    let half = 0.5 * theta;
    {
        let (s, c) = (libm::sin(half), libm::cos(half));
        [w[0] * c - w[1] * s, w[0] * s + w[1] * c]
    }
}

fn centroid(pts: &[Q]) -> Q {
    let n = pts.len();
    let (mut cx, mut cy, mut a2) = (0.0, 0.0, 0.0);
    for k in 0..n {
        let (p, q) = (pts[k], pts[(k + 1) % n]);
        let c = cross(p, q);
        a2 += c;
        cx += (p[0] + q[0]) * c;
        cy += (p[1] + q[1]) * c;
    }
    if a2.abs() > 0.0 {
        [cx / (3.0 * a2), cy / (3.0 * a2)]
    } else {
        let s = pts.iter().fold([0.0, 0.0], |s, p| [s[0] + p[0], s[1] + p[1]]);
        [s[0] / n.max(1) as f64, s[1] / n.max(1) as f64]
    }
}

/// The region `r`'s loops as polygons (the region on the left of each).
fn region_polys(arr: &Arrangement, r: usize, locked_curves: &[bool]) -> Option<Vec<Poly>> {
    let region = arr.regions.get(r)?;
    let mut out = Vec::with_capacity(region.loops.len());
    for lp in &region.loops {
        let mut poly = Poly { pts: Vec::new(), end: Vec::new(), locked: Vec::new() };
        for &(c, fwd) in lp {
            let curve = arr.curves.get(c)?;
            let mut ids = curve.points.clone();
            if !fwd {
                ids.reverse();
            }
            ids.pop();
            let locked = locked_curves.get(c).copied().unwrap_or(false);
            for (k, id) in ids.into_iter().enumerate() {
                poly.pts.push(*arr.points.get(id)?);
                poly.end.push(k == 0);
                poly.locked.push(locked);
            }
        }
        if poly.len() < 2 {
            return None;
        }
        out.push(poly);
    }
    Some(out)
}

/// A point where a cut ends: a loop point, or a point inside edge `k` (`k` → `k + 1`).
#[derive(Debug, Clone, Copy, PartialEq)]
enum End {
    Point(usize, usize),
    Edge(usize, usize, Q),
}

impl End {
    fn loop_index(self) -> usize {
        match self {
            End::Point(l, _) | End::Edge(l, _, _) => l,
        }
    }

    fn at(self, loops: &[Poly]) -> Q {
        match self {
            End::Point(l, i) => loops.get(l).map_or([0.0; 2], |p| p.at(i)),
            End::Edge(_, _, p) => p,
        }
    }
}

/// The first boundary point the ray from `o` along unit `d` meets further than `min_t`, skipping
/// the edges at loop point `from` (where the ray starts), as an end snapped to a loop point
/// within `snap`, with its distance.
fn cast(loops: &[Poly], from: Option<(usize, usize)>, o: Q, d: Q, snap: f64, min_t: f64) -> Option<(f64, End)> {
    let mut best: Option<(f64, usize, usize)> = None;
    for (l, poly) in loops.iter().enumerate() {
        let n = poly.len();
        for k in 0..n {
            if let Some((fl, fi)) = from
                && fl == l
                && (k == fi || (k + 1) % n == fi)
            {
                continue;
            }
            let (a, b) = (poly.at(k), poly.at(k + 1));
            let e = sub(b, a);
            let den = cross(d, e);
            if den.abs() <= MIN_SINE * len(e) {
                continue;
            }
            let w = sub(a, o);
            let t = cross(w, e) / den;
            let u = cross(w, d) / den;
            if t > min_t && (-ON_SEGMENT..=1.0 + ON_SEGMENT).contains(&u) && best.is_none_or(|(bt, ..)| t < bt) {
                best = Some((t, l, k));
            }
        }
    }
    let (t, l, k) = best?;
    let poly = loops.get(l)?;
    let x = [o[0] + t * d[0], o[1] + t * d[1]];
    let n = poly.len();
    let (a, b) = (poly.at(k), poly.at(k + 1));
    let near = [(k, a), ((k + 1) % n, b)]
        .into_iter()
        .filter(|(_, p)| len(sub(*p, x)) <= snap)
        .min_by(|x1, x2| len(sub(x1.1, x)).total_cmp(&len(sub(x2.1, x))));
    Some((
        t,
        match near {
            Some((i, _)) => End::Point(l, i),
            None => End::Edge(l, k, x),
        },
    ))
}

/// The nearest curve end of loop `l` to `x` within `radius` (an alternative end for a cut).
fn nearby_end(loops: &[Poly], l: usize, x: Q, radius: f64) -> Option<End> {
    let poly = loops.get(l)?;
    (0..poly.len())
        .filter(|&i| poly.is_end(i) && len(sub(poly.at(i), x)) <= radius)
        .min_by(|&i, &j| len(sub(poly.at(i), x)).total_cmp(&len(sub(poly.at(j), x))))
        .map(|i| End::Point(l, i))
}

/// Do the open segments `ab` and `cd` cross or touch?
fn segments_meet(a: Q, b: Q, c: Q, d: Q) -> bool {
    let o = |p: Q, q: Q, r: Q| cross(sub(q, p), sub(r, p));
    let (d1, d2, d3, d4) = (o(a, b, c), o(a, b, d), o(c, d, a), o(c, d, b));
    (d1 > 0.0) != (d2 > 0.0) && (d3 > 0.0) != (d4 > 0.0) && d1 != 0.0 && d2 != 0.0 && d3 != 0.0 && d4 != 0.0
}

/// Distance from `p` to segment `ab`.
fn seg_dist(p: Q, a: Q, b: Q) -> f64 {
    super::arrangement::seg_dist2(p, a, b).0.sqrt()
}

/// Is `p` inside the region (inside the outer loop, outside the holes)?
fn inside(loops: &[Poly], p: Q) -> bool {
    let winding = |poly: &Poly| -> i32 {
        let mut wn = 0;
        for k in 0..poly.len() {
            let (a, b) = (poly.at(k), poly.at(k + 1));
            let side = cross(sub(b, a), sub(p, a));
            if a[1] <= p[1] {
                if b[1] > p[1] && side > 0.0 {
                    wn += 1;
                }
            } else if b[1] <= p[1] && side < 0.0 {
                wn -= 1;
            }
        }
        wn
    };
    loops.iter().map(winding).sum::<i32>() != 0
}

/// Whether the cut from `a` to `b` is admissible: inside the region, clear of other boundary
/// points, square enough to the boundary at both ends, not ending inside a locked edge, and not
/// along the boundary.
fn admissible(loops: &[Poly], a: End, b: End, h: f64) -> bool {
    let (pa, pb) = (a.at(loops), b.at(loops));
    let l = len(sub(pb, pa));
    if !(l > 0.25 * h) || a == b {
        return false;
    }
    for e in [a, b] {
        if let End::Edge(li, k, _) = e
            && loops.get(li).and_then(|p| p.locked.get(k)).copied().unwrap_or(true)
        {
            return false;
        }
    }
    // The angles the cut makes with the boundary at each end, on both sides.
    let Some(d) = unit(sub(pb, pa)) else { return false };
    for (e, dir) in [(a, d), (b, [-d[0], -d[1]])] {
        let (prev, next) = match e {
            End::Point(li, i) => {
                let Some(poly) = loops.get(li) else { return false };
                let n = poly.len();
                (poly.at(i + n - 1), poly.at(i + 1))
            }
            End::Edge(li, k, _) => {
                let Some(poly) = loops.get(li) else { return false };
                (poly.at(k), poly.at(k + 1))
            }
        };
        let o = e.at(loops);
        let tip = [o[0] + dir[0], o[1] + dir[1]];
        // The cut must leave into the region: between the outgoing and the incoming edge.
        let (a1, a2) = (interior(prev, o, tip), interior(tip, o, next));
        let total = interior(prev, o, next);
        if (a1 + a2 - total).abs() > SAME_ANGLE || a1 < MIN_CUT_ANGLE || a2 < MIN_CUT_ANGLE {
            return false;
        }
    }
    // No crossing, and clear of every other boundary point. An end on an edge touches that
    // edge by construction: its crossing test would read only the rounding of the end point (the
    // sign of a zero orientation), so a cut was taken or refused by noise at 1e-15 (#605).
    let clear = 0.05 * h.min(l);
    let on = |x: End, li: usize, k: usize| matches!(x, End::Edge(l2, k2, _) if l2 == li && k2 == k);
    for (li, poly) in loops.iter().enumerate() {
        for k in 0..poly.len() {
            let (c, e) = (poly.at(k), poly.at(k + 1));
            if !on(a, li, k) && !on(b, li, k) && segments_meet(pa, pb, c, e) {
                return false;
            }
            let is_end = |x: End| match x {
                End::Point(l2, i) => l2 == li && i == k,
                End::Edge(..) => false,
            };
            if !is_end(a) && !is_end(b) && seg_dist(c, pa, pb) < clear {
                return false;
            }
        }
    }
    inside(loops, [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])])
}

/// The region split by the cut from `a` to `b` (both on one loop) into two, each a list of loops
/// (holes go to the piece that contains them); `None` when the cut does not split it.
fn split(loops: &[Poly], a: End, b: End) -> Option<Vec<Vec<Poly>>> {
    // Rebuild the loops with the cut's ends inserted, and find their new positions.
    let mut rebuilt: Vec<Poly> = Vec::with_capacity(loops.len());
    let mut pos: [Option<(usize, usize)>; 2] = [None, None];
    for (l, poly) in loops.iter().enumerate() {
        let mut p = Poly { pts: Vec::with_capacity(poly.len() + 2), end: Vec::new(), locked: Vec::new() };
        for k in 0..poly.len() {
            let locked = poly.locked.get(k).copied().unwrap_or(false);
            p.pts.push(poly.at(k));
            p.end.push(poly.is_end(k));
            p.locked.push(locked);
            for (slot, e) in [a, b].into_iter().enumerate() {
                if e == End::Point(l, k)
                    && let (Some(x), Some(last)) = (pos.get_mut(slot), p.end.last_mut())
                {
                    *last = true;
                    *x = Some((l, p.pts.len() - 1));
                }
            }
            // Points inside this edge, in order along it.
            let mut inside: Vec<(f64, usize, Q)> = [a, b]
                .into_iter()
                .enumerate()
                .filter_map(|(slot, e)| match e {
                    End::Edge(le, ke, x) if le == l && ke == k => Some((len(sub(x, poly.at(k))), slot, x)),
                    _ => None,
                })
                .collect();
            inside.sort_by(|x, y| x.0.total_cmp(&y.0));
            for (_, slot, x) in inside {
                p.pts.push(x);
                p.end.push(true);
                p.locked.push(locked);
                if let Some(s) = pos.get_mut(slot) {
                    *s = Some((l, p.pts.len() - 1));
                }
            }
        }
        rebuilt.push(p);
    }
    let [Some((la, ia)), Some((lb, ib))] = pos else { return None };
    if (la, ia) == (lb, ib) {
        return None;
    }
    let loops = rebuilt;
    // A cut between two loops would join them (the hole steps handle holes).
    if la != lb {
        return None;
    }
    let poly = loops.get(la)?;
    let n = poly.len();
    let piece = |from: usize, to: usize| -> Poly {
        let mut p = Poly { pts: Vec::new(), end: Vec::new(), locked: Vec::new() };
        let count = (to + n - from) % n + 1;
        for k in 0..count {
            let i = (from + k) % n;
            p.pts.push(poly.at(i));
            p.end.push(poly.is_end(i));
            // The last edge of the piece is the cut.
            p.locked.push(if k + 1 == count { false } else { poly.locked.get(i).copied().unwrap_or(false) });
        }
        p
    };
    let (pa, pb) = (piece(ia, ib), piece(ib, ia));
    if pa.len() < 3 || pb.len() < 3 {
        return None;
    }
    let holes: Vec<Poly> = loops.iter().enumerate().filter(|(l, _)| *l != la).map(|(_, p)| p.clone()).collect();
    let mut out = vec![vec![pa], vec![pb]];
    for h in holes {
        let probe = h.at(0);
        let first = out.first().and_then(|p| p.first()).is_some_and(|p| inside(std::slice::from_ref(p), probe));
        if let Some(piece) = out.get_mut(usize::from(!first)) {
            piece.push(h);
        }
    }
    Some(out)
}

/// How good a piece is: `None` when no structured scheme fits it, else the worst
/// [`quad_score`] of its mesh on the target size.
fn score(piece: &[Poly], s: &Settings) -> Option<f64> {
    let [poly] = piece else { return ring_score(piece, s) };
    if poly.len() > MAX_PIECE_POINTS {
        return None;
    }
    let (curves, _) = poly.curves();
    let scheme = region::classify(std::slice::from_ref(&curves));
    if scheme == Scheme::Free {
        return None;
    }
    mesh_score(&[curves], scheme, s)
}

/// A piece with one hole: an annulus when its two loops are smooth, else not structured.
fn ring_score(piece: &[Poly], s: &Settings) -> Option<f64> {
    if piece.len() != 2 || piece.iter().map(Poly::len).sum::<usize>() > MAX_PIECE_POINTS {
        return None;
    }
    let curves: Vec<Vec<Vec<Q>>> = piece.iter().map(|p| p.curves().0).collect();
    let scheme = region::classify(&curves);
    if scheme != Scheme::Annulus {
        return None;
    }
    mesh_score(&curves, scheme, s)
}

/// Mesh loops of curves with a structured scheme on counts from their lengths, and measure.
fn mesh_score(loops: &[Vec<Vec<Q>>], scheme: Scheme, s: &Settings) -> Option<f64> {
    let mut goals: Vec<Vec<f64>> = Vec::with_capacity(loops.len());
    for lp in loops {
        let mut g = Vec::with_capacity(lp.len());
        for c in lp {
            let l: f64 = c.windows(2).map(|w| len(sub(w[1], w[0]))).sum();
            let mid = c.get(c.len() / 2).copied().unwrap_or([0.0; 2]);
            let goal = l / s.size(mid);
            if !(goal.is_finite() && goal <= 1e5) {
                return None;
            }
            g.push(goal);
        }
        goals.push(g);
    }
    // A fine mesh is scored on a coarser one of the same shape (the quality of a structured
    // block depends on its shape, the cost on its size).
    let total: f64 = goals.iter().flatten().sum();
    let coarsen = if total > SCORE_RESOLUTION { SCORE_RESOLUTION / total } else { 1.0 };
    let mut ia = Intervals::default();
    let vars: Vec<Vec<usize>> = goals.iter().map(|g| g.iter().map(|&x| ia.var(x * coarsen, 1)).collect()).collect();
    region::constrain(&mut ia, &vars, scheme);
    let counts = ia.solve().ok()?;
    let mut nodes: Vec<Vec<Q>> = Vec::with_capacity(loops.len());
    let mut starts: Vec<Vec<usize>> = Vec::with_capacity(loops.len());
    for (lp, v) in loops.iter().zip(&vars) {
        let mut pts = Vec::new();
        let mut st = Vec::new();
        for (c, &var) in lp.iter().zip(v) {
            st.push(pts.len());
            let n = counts.get(var).copied().unwrap_or(1);
            let mut d = divide(c, n);
            d.pop();
            pts.extend(d);
        }
        nodes.push(pts);
        starts.push(st);
    }
    let m = region::mesh_structured(&nodes, &starts, scheme, s.field, s.opts).ok()??;
    let worst = m.quads.iter().map(|q| quad_score(&q.map(|i| m.nodes.get(i as usize).copied().unwrap_or([0.0; 2])), s)).fold(f64::INFINITY, f64::min);
    worst.is_finite().then_some(worst)
}

/// A quadrilateral's quality for the scores: its smallest corner scaled Jacobian, and no more
/// than its aspect allows, its shortest edge against its longest or the target size (the layers
/// it is swept into are that thick): a thin strip or a small piece is square-cornered but
/// makes slivers.
fn quad_score(c: &[Q; 4], s: &Settings) -> f64 {
    let mut worst = f64::INFINITY;
    for i in 0..4 {
        worst = worst.min(planar_corner(c, i).1.unwrap_or(-1.0));
    }
    let edges = (0..4).map(|i| len(sub(c[(i + 1) % 4], c[i])));
    let (short, long) = edges.fold((f64::INFINITY, 0.0_f64), |(a, b), e| (a.min(e), b.max(e)));
    let long = long.max(s.size(centroid(c)));
    if long > 0.0 && long.is_finite() {
        worst = worst.min(ASPECT_FREE * short / long);
    }
    worst
}

/// `n + 1` points evenly along a polyline by length.
fn divide(c: &[Q], n: usize) -> Vec<Q> {
    let n = n.max(1);
    let lengths: Vec<f64> = c.windows(2).map(|w| len(sub(w[1], w[0]))).collect();
    let total: f64 = lengths.iter().sum();
    let mut out = Vec::with_capacity(n + 1);
    let (mut seg, mut acc) = (0usize, 0.0);
    for j in 0..=n {
        let target = total * j as f64 / n as f64;
        while seg + 1 < lengths.len() && acc + lengths.get(seg).copied().unwrap_or(0.0) < target {
            acc += lengths.get(seg).copied().unwrap_or(0.0);
            seg += 1;
        }
        let (a, b) = (c.get(seg).copied().unwrap_or([0.0; 2]), c.get(seg + 1).copied().unwrap_or([0.0; 2]));
        let l = lengths.get(seg).copied().unwrap_or(0.0);
        let t = if l > 0.0 { ((target - acc) / l).clamp(0.0, 1.0) } else { 0.0 };
        out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
    }
    if let (Some(first), Some(c0)) = (out.first_mut(), c.first()) {
        *first = *c0;
    }
    if let (Some(last), Some(c1)) = (out.last_mut(), c.last()) {
        *last = *c1;
    }
    out
}

/// The concave corners and concave rounded curves of a piece's loops.
fn concave_features(piece: &[Poly]) -> usize {
    let mut count = 0;
    for poly in piece {
        let n = poly.len();
        count += (0..n).filter(|&i| poly.is_end(i) && poly.angle(i) >= CONCAVE).count();
        let (curves, _) = poly.curves();
        for c in curves {
            let turning: f64 = (1..c.len().saturating_sub(1)).map(|k| std::f64::consts::PI - interior(c[k - 1], c[k], c[k + 1])).sum();
            if turning <= -ROUNDED {
                count += 1;
            }
        }
    }
    count
}

/// The value of a set of pieces: free pieces, quality shortfall, pieces.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Value {
    free: usize,
    shortfall: f64,
    pieces: usize,
}

impl Value {
    /// A piece's value from its score, or (free) from its concave features: each cut that takes
    /// one off counts as progress even while the piece is still free.
    fn of(score: Option<f64>, piece: &[Poly]) -> Value {
        match score {
            Some(q) => Value { free: 0, shortfall: (TARGET - q).max(0.0).powi(2), pieces: 1 },
            None => Value { free: 1 + concave_features(piece), shortfall: 0.0, pieces: 1 },
        }
    }

    fn add(self, o: Value) -> Value {
        Value { free: self.free + o.free, shortfall: self.shortfall + o.shortfall, pieces: self.pieces + o.pieces }
    }

    fn better(self, o: Value) -> bool {
        if self.free != o.free {
            return self.free < o.free;
        }
        if (self.shortfall - o.shortfall).abs() > MIN_GAIN * (self.pieces.max(o.pieces) - self.pieces.min(o.pieces)).max(1) as f64 {
            return self.shortfall < o.shortfall;
        }
        self.pieces < o.pieces
    }
}

/// The candidate cuts of a piece without holes (see the module docs).
fn candidates(loops: &[Poly], s: &Settings) -> Vec<(End, End)> {
    let mut out: Vec<(End, End, f64)> = Vec::new();
    let Some(poly) = loops.first() else { return Vec::new() };
    let n = poly.len();
    // Rays: (loop point, direction).
    let mut rays: Vec<(usize, Q)> = Vec::new();
    for i in 0..n {
        let (a, o, b) = (poly.at(i + n - 1), poly.at(i), poly.at(i + 1));
        let (Some(u), Some(w)) = (unit(sub(o, a)), unit(sub(b, o))) else { continue };
        let theta = poly.angle(i);
        // The inward bisector: the incoming direction turned left by half the interior angle
        // from the outgoing one.
        let bis = bisector(w, theta);
        if poly.is_end(i) {
            if theta >= CONCAVE {
                rays.push((i, u));
                rays.push((i, [-w[0], -w[1]]));
                rays.push((i, bis));
            } else if theta <= 0.83 * std::f64::consts::PI {
                rays.push((i, bis));
            }
        }
    }
    // Strongly turning curves: their ends and the point of half their turning, along the normal.
    let (curves, starts) = poly.curves();
    for (c, &st) in curves.iter().zip(&starts) {
        let m = c.len();
        if m < 3 {
            continue;
        }
        let turns: Vec<f64> = (1..m - 1).map(|k| std::f64::consts::PI - interior(c[k - 1], c[k], c[k + 1])).collect();
        let total: f64 = turns.iter().sum();
        if total.abs() < ROUNDED {
            continue;
        }
        let mut acc = 0.0;
        let mut half = 1;
        for (k, t) in turns.iter().enumerate() {
            acc += t;
            if acc.abs() >= 0.5 * total.abs() {
                half = k + 1;
                break;
            }
        }
        for k in [0, half, m - 1] {
            let i = (st + k) % n;
            let (a, o, b) = (poly.at(i + n - 1), poly.at(i), poly.at(i + 1));
            let (Some(u), Some(w)) = (unit(sub(o, a)), unit(sub(b, o))) else { continue };
            // The inward normal: left of the mean tangent.
            if let Some(t) = unit([u[0] + w[0], u[1] + w[1]]) {
                rays.push((i, [-t[1], t[0]]));
            }
            if k == 0 || k == m - 1 {
                // At a curve end, also square to the neighbouring curve.
                rays.push((i, [-u[1], u[0]]));
                rays.push((i, [-w[1], w[0]]));
            }
        }
    }
    for (i, d) in rays {
        let o = poly.at(i);
        let h = s.size(o);
        let Some((t, end)) = cast(loops, Some((0, i)), o, d, 0.3 * h, s.tol) else { continue };
        let a = End::Point(0, i);
        let x = [o[0] + t * d[0], o[1] + t * d[1]];
        let mut push = |b: End| {
            let same = |p: End, q: End| len(sub(p.at(loops), q.at(loops))) <= s.tol;
            let known = out.iter().any(|&(p, q, _)| (same(p, a) && same(q, b)) || (same(p, b) && same(q, a)));
            if !known && admissible(loops, a, b, h) {
                out.push((a, b, t));
            }
        };
        push(end);
        if let Some(alt) = nearby_end(loops, end.loop_index(), x, (0.75 * h).min(0.3 * t)) {
            push(alt);
        }
    }
    out.sort_by(|x, y| x.2.total_cmp(&y.2));
    out.truncate(MAX_CANDIDATES);
    out.into_iter().map(|(a, b, _)| (a, b)).collect()
}

/// The best value reachable from a piece with at most `depth` more cuts, and the first cut.
fn search(piece: &[Poly], depth: usize, s: &Settings, budget: &mut usize) -> (Value, Option<(End, End)>) {
    let here = score(piece, s);
    let leaf = Value::of(here, piece);
    // A disk is its own O-grid: cut into pieces round its centre it only gets worse, however a
    // coarse tessellation makes its curve ends look like corners.
    let disk = || matches!(piece, [poly] if region::classify(std::slice::from_ref(&poly.curves().0)) == Scheme::Disk);
    if depth == 0 || here.is_some_and(|q| q >= GOOD) || piece.len() != 1 || *budget == 0 || disk() {
        return (leaf, None);
    }
    let mut best = (leaf, None);
    for (a, b) in candidates(piece, s) {
        if *budget == 0 {
            break;
        }
        *budget -= 1;
        let Some(parts) = split(piece, a, b) else { continue };
        let mut total: Option<Value> = None;
        for part in &parts {
            let (v, _) = search(part, depth - 1, s, budget);
            total = Some(total.map_or(v, |t| t.add(v)));
        }
        if let Some(v) = total
            && v.better(best.0)
        {
            best = (v, Some((a, b)));
        }
    }
    best
}

/// Lines across the gaps between neighbouring holes of a region with several (step 1).
fn separate_holes(loops: &[Poly], s: &Settings) -> Vec<Vec<Q>> {
    let centres: Vec<Q> = loops.iter().map(|p| centroid(&p.pts)).collect();
    let mut out: Vec<Vec<Q>> = Vec::new();
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 1..loops.len() {
        let centre = |k: usize| centres.get(k).copied().unwrap_or([0.0; 2]);
        let nearest = (1..loops.len()).filter(|&j| j != i).min_by(|&a, &b| len(sub(centre(a), centre(i))).total_cmp(&len(sub(centre(b), centre(i)))));
        let Some(j) = nearest else { continue };
        let key = (i.min(j), i.max(j));
        if pairs.contains(&key) {
            continue;
        }
        pairs.push(key);
        let (ci, cj) = (centre(i), centre(j));
        let Some(axis) = unit(sub(cj, ci)) else { continue };
        // The middle of the gap: halfway between where the line of centres leaves each hole.
        let leave = |from: Q, d: Q, hole: usize| -> Option<Q> {
            let poly = loops.get(hole)?;
            let mut best: Option<f64> = None;
            for k in 0..poly.len() {
                let (a, b) = (poly.at(k), poly.at(k + 1));
                let e = sub(b, a);
                let den = cross(d, e);
                if den.abs() <= MIN_SINE * len(e) {
                    continue;
                }
                let w = sub(a, from);
                let (t, u) = (cross(w, e) / den, cross(w, d) / den);
                if t > 0.0 && (0.0..=1.0).contains(&u) && best.is_none_or(|bt| t > bt) {
                    best = Some(t);
                }
            }
            best.map(|t| [from[0] + t * d[0], from[1] + t * d[1]])
        };
        let (Some(pi), Some(pj)) = (leave(ci, axis, i), leave(cj, [-axis[0], -axis[1]], j)) else { continue };
        let m = [0.5 * (pi[0] + pj[0]), 0.5 * (pi[1] + pj[1])];
        let perp = [-axis[1], axis[0]];
        let h = s.size(m);
        let hit = |d: Q| cast(loops, None, m, d, 0.3 * h, s.tol);
        let (Some((_, e1)), Some((_, e2))) = (hit(perp), hit([-perp[0], -perp[1]])) else { continue };
        // A separator that ends on either of the two holes does not separate them.
        if [e1, e2].iter().any(|e| e.loop_index() == i || e.loop_index() == j) {
            continue;
        }
        let (p1, p2) = (e1.at(loops), e2.at(loops));
        if len(sub(p1, p2)) > 0.25 * h {
            out.push(vec![p1, p2]);
        }
    }
    out
}

/// The pieces a region falls into with `cuts` added (a new arrangement of its loops and the
/// cuts), each as its loops; `None` when the arrangement fails.
fn pieces_after(loops: &[Poly], cuts: &[Vec<Q>], tol: f64) -> Option<Vec<Vec<Poly>>> {
    let mut polylines: Vec<(Vec<Q>, usize)> = Vec::new();
    let mut vertices: Vec<Q> = Vec::new();
    for poly in loops {
        let mut pts = poly.pts.clone();
        pts.push(*poly.pts.first()?);
        polylines.push((pts, 0));
        vertices.extend((0..poly.len()).filter(|&i| poly.is_end(i)).map(|i| poly.at(i)));
        if !(0..poly.len()).any(|i| poly.is_end(i)) {
            vertices.push(poly.at(0));
        }
    }
    for c in cuts {
        vertices.extend([*c.first()?, *c.last()?]);
        polylines.push((c.clone(), CUT));
    }
    let arr = super::arrangement::build(&vertices, &polylines, tol, 0.0).ok()?;
    let mut out = Vec::new();
    for (r, region) in arr.regions.iter().enumerate() {
        let p = super::layered::interior_point(&arr, region).ok()?;
        if inside(loops, p) {
            let mut piece = region_polys(&arr, r, &[])?;
            // An edge on a locked edge of the region stays locked.
            for poly in &mut piece {
                for k in 0..poly.len() {
                    let (a, b) = (poly.at(k), poly.at(k + 1));
                    let m = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
                    let on_locked = loops
                        .iter()
                        .any(|l| (0..l.len()).any(|j| l.locked.get(j).copied().unwrap_or(false) && seg_dist(m, l.at(j), l.at(j + 1)) <= 4.0 * tol));
                    if let Some(x) = poly.locked.get_mut(k) {
                        *x = on_locked;
                    }
                }
            }
            out.push(piece);
        }
    }
    Some(out)
}

/// Whether a region cut by `cuts` can be given division counts that every piece's scheme
/// accepts: the interval assignment of its pieces together, as the cross-section's will be. Each
/// piece fitting on its own is not enough: mapped pieces round a hole chain their opposite sides,
/// and a chain that comes back to a side with one more curve on it has no solution (#397).
fn assignable(loops: &[Poly], cuts: &[Vec<Q>], s: &Settings) -> bool {
    let mut polylines: Vec<(Vec<Q>, usize)> = Vec::new();
    let mut vertices: Vec<Q> = Vec::new();
    for poly in loops {
        let (curves, _) = poly.curves();
        for c in curves {
            vertices.extend(c.first().into_iter().chain(c.last()).copied());
            polylines.push((c, 0));
        }
    }
    for c in cuts {
        vertices.extend(c.first().into_iter().chain(c.last()).copied());
        polylines.push((c.clone(), CUT));
    }
    let Ok(arr) = super::arrangement::build(&vertices, &polylines, s.tol, 0.0) else { return false };
    let mut used = Vec::with_capacity(arr.regions.len());
    for region in &arr.regions {
        let Ok(p) = super::layered::interior_point(&arr, region) else { return false };
        used.push(inside(loops, p));
    }
    let mut ia = Intervals::default();
    for curve in &arr.curves {
        let pts: Vec<Q> = curve.points.iter().filter_map(|&p| arr.points.get(p).copied()).collect();
        let length: f64 = pts.windows(2).map(|w| len(sub(w[1], w[0]))).sum();
        ia.var(length / s.size(pts.first().copied().unwrap_or([0.0; 2])), 1);
    }
    super::section::schemes(&arr, &used, &mut ia);
    ia.solve().is_ok()
}

/// O-grid cuts round the one hole of a region (step 2): of the layouts below, the one whose
/// pieces (each allowed one more cut) score best.
/// - From the outer loop's corners towards the hole's centre.
/// - From the outer loop's corners along their bisectors (towards the centre where a bisector
///   misses the hole).
/// - From the hole's concave corners away from its centre.
/// - From the hole along four directions (the axes, then the diagonals): the outer loop's
///   corners are then cut off by the pieces' own search, the layout of a long plate.
fn o_grid(loops: &[Poly], s: &Settings, budget: &mut usize) -> Vec<Vec<Q>> {
    let (Some(outer), Some(hole)) = (loops.first(), loops.get(1)) else { return Vec::new() };
    let c = centroid(&hole.pts);
    let corner = |p: &Poly, i: usize| p.is_end(i) && region::is_corner(p.angle(i));
    let concave = |p: &Poly, i: usize| p.is_end(i) && region::is_reflex(p.angle(i));
    // A cut from loop point (l, i) along d that must end on loop `target`.
    let ray = |l: usize, i: usize, d: Q, target: usize| -> Option<Vec<Q>> {
        let o = loops.get(l)?.at(i);
        let h = s.size(o);
        let (_, end) = cast(loops, Some((l, i)), o, d, 0.3 * h, s.tol)?;
        (end.loop_index() == target && admissible(loops, End::Point(l, i), end, h)).then(|| vec![o, end.at(loops)])
    };
    let outer_corners: Vec<usize> = (0..outer.len()).filter(|&i| corner(outer, i)).collect();
    let mut layouts: Vec<Vec<Vec<Q>>> = Vec::new();
    layouts.push(outer_corners.iter().filter_map(|&i| ray(0, i, unit(sub(c, outer.at(i)))?, 1)).collect());
    layouts.push(
        outer_corners
            .iter()
            .filter_map(|&i| {
                let o = outer.at(i);
                let w = unit(sub(outer.at(i + 1), o))?;
                let bis = bisector(w, outer.angle(i));
                ray(0, i, bis, 1).or_else(|| ray(0, i, unit(sub(c, o))?, 1))
            })
            .collect(),
    );
    layouts.push((0..hole.len()).filter(|&i| concave(hole, i)).filter_map(|i| ray(1, i, unit(sub(hole.at(i), c))?, 0)).collect());
    // From the hole's concave corners along both of its sides there, continued: a rectangular
    // opening becomes side strips and corner blocks whose lines carry on the opening's sides
    // (submapping, Whiteley et al.), so every corner stays near 90° where a diagonal would
    // leave 45° wedges.
    let n = hole.len();
    let extensions = (0..n).filter(|&i| concave(hole, i)).flat_map(|i| {
        let o = hole.at(i);
        [unit(sub(o, hole.at(i + n - 1))), unit(sub(o, hole.at(i + 1)))].into_iter().flatten().map(move |d| (i, d))
    });
    layouts.push(extensions.filter_map(|(i, d)| ray(1, i, d, 0)).collect());
    let s2 = std::f64::consts::FRAC_1_SQRT_2;
    for dirs in [[[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]], [[s2, s2], [-s2, s2], [-s2, -s2], [s2, -s2]]] {
        let mut cuts = Vec::new();
        for d in dirs {
            let Some((_, start)) = cast(std::slice::from_ref(hole), None, c, d, 0.0, 0.0) else { continue };
            let start = match start {
                End::Point(_, i) => End::Point(1, i),
                End::Edge(_, k, x) => End::Edge(1, k, x),
            };
            let o = start.at(loops);
            let h = s.size(o);
            // From just off the hole, so the ray starts inside the region.
            let o2 = [o[0] + d[0] * s.tol * 4.0, o[1] + d[1] * s.tol * 4.0];
            let Some((_, end)) = cast(loops, None, o2, d, 0.3 * h, s.tol) else { continue };
            if end.loop_index() == 0 && admissible(loops, start, end, h) {
                cuts.push(vec![o, end.at(loops)]);
            }
        }
        layouts.push(cuts);
    }
    let mut best: Option<(Value, Vec<Vec<Q>>)> = None;
    for cuts in layouts {
        if cuts.len() < 2 {
            continue;
        }
        let Some(pieces) = pieces_after(loops, &cuts, s.tol) else { continue };
        if !assignable(loops, &cuts, s) {
            continue;
        }
        let mut value: Option<Value> = None;
        for piece in &pieces {
            let (v, _) = search(piece, 1, s, budget);
            value = Some(value.map_or(v, |t| t.add(v)));
        }
        if let Some(v) = value
            && best.as_ref().is_none_or(|(b, _)| v.better(*b))
        {
            best = Some((v, cuts));
        }
    }
    best.map(|(_, c)| c).unwrap_or_default()
}

/// The cuts that decompose one region (its loops), and the pieces it ends in.
fn decompose_region(loops: &[Poly], s: &Settings) -> (Vec<Vec<Q>>, Vec<Vec<Poly>>) {
    let mut cuts: Vec<Vec<Q>> = Vec::new();
    let mut leaves: Vec<Vec<Poly>> = Vec::new();
    let mut work = vec![loops.to_vec()];
    let mut steps = 0;
    let mut remaining = REGION_BUDGET;
    while let Some(piece) = work.pop() {
        steps += 1;
        if steps > MAX_STEPS || remaining == 0 || piece.iter().map(Poly::len).sum::<usize>() > MAX_PIECE_POINTS {
            leaves.push(piece);
            continue;
        }
        let step: Vec<Vec<Q>> = match piece.len() {
            0 => continue,
            1 => {
                let mut budget = SEARCH_BUDGET.min(remaining);
                let (_, first) = search(&piece, 2, s, &mut budget);
                remaining -= SEARCH_BUDGET.min(remaining) - budget;
                if let Some((a, b)) = first
                    && let Some(parts) = split(&piece, a, b)
                {
                    cuts.push(vec![a.at(&piece), b.at(&piece)]);
                    work.extend(parts);
                    continue;
                }
                Vec::new()
            }
            2 if ring_score(&piece, s).is_some() => Vec::new(),
            2 => o_grid(&piece, s, &mut remaining),
            _ => separate_holes(&piece, s),
        };
        match (step.is_empty(), pieces_after(&piece, &step, s.tol)) {
            (false, Some(parts)) if parts.len() > 1 => {
                cuts.extend(step);
                work.extend(parts);
            }
            _ => leaves.push(piece),
        }
    }
    (cuts, leaves)
}

/// The score of the free all-quad mesh of a region on the target size, as a piece's (0 when it
/// cannot be meshed): what the region gets without cuts.
fn free_score(loops: &[Poly], s: &Settings) -> f64 {
    let min = if loops.len() > 1 { 6 } else { 4 };
    let mut nodes: Vec<Vec<Q>> = Vec::with_capacity(loops.len());
    for poly in loops {
        let (curves, _) = poly.curves();
        let goals: Vec<f64> =
            curves.iter().map(|c| c.windows(2).map(|w| len(sub(w[1], w[0]))).sum::<f64>() / s.size(c.first().copied().unwrap_or([0.0; 2]))).collect();
        if goals.iter().any(|g| !(g.is_finite() && *g <= 1e5)) {
            return 0.0;
        }
        // Counts near the goals, the loop even and long enough for the free mesher (as the
        // section's interval assignment makes them).
        let mut counts: Vec<usize> = goals.iter().map(|g| (g.round() as usize).max(1)).collect();
        let coarsest = |counts: &[usize]| {
            (0..counts.len()).max_by(|&i, &j| {
                let r = |k: usize| goals.get(k).copied().unwrap_or(0.0) / counts.get(k).copied().unwrap_or(1) as f64;
                r(i).total_cmp(&r(j)).then(j.cmp(&i))
            })
        };
        while counts.iter().sum::<usize>() < min || counts.iter().sum::<usize>() % 2 == 1 {
            let Some(k) = coarsest(&counts) else { return 0.0 };
            if let Some(c) = counts.get_mut(k) {
                *c += 1;
            }
        }
        let mut pts = Vec::new();
        for (c, &n) in curves.iter().zip(&counts) {
            let mut d = divide(c, n);
            d.pop();
            pts.extend(d);
        }
        nodes.push(pts);
    }
    let Ok(m) = crate::quad::free_quad_mesh(&nodes, s.field, s.opts) else { return 0.0 };
    let worst = m.quads.iter().map(|q| quad_score(&q.map(|i| m.nodes.get(i as usize).copied().unwrap_or([0.0; 2])), s)).fold(f64::INFINITY, f64::min);
    if worst.is_finite() { worst.max(0.0) } else { 0.0 }
}

/// The cuts that decompose the filled regions of `arr` (`used`) into blocks the structured
/// schemes fit, as polylines to add to the arrangement. Each region is decomposed on its own
/// (holes apart, O-grids, then the pieces), and kept when its blocks are fair or not worse than
/// the region's free mesh. `locked[c]`: curve `c` carries an edge seed, so no cut may end inside it.
pub(crate) fn cuts(arr: &Arrangement, used: &[bool], locked: &[bool], s: &Settings) -> Result<Vec<Vec<Q>>> {
    if !(s.tol > 0.0 && s.tol.is_finite()) {
        return Err(Error::Invalid("the cut tolerance must be a positive number".into()));
    }
    let mut out = Vec::new();
    for r in 0..arr.regions.len() {
        if !used.get(r).copied().unwrap_or(false) {
            continue;
        }
        let Some(loops) = region_polys(arr, r, locked) else { continue };
        // A structured region with a seeded curve is meshed as its seeds set it up: cuts would
        // add side equalities the exact seeds may not meet.
        let seeded = loops.iter().any(|p| p.locked.iter().any(|&l| l));
        if seeded && score(&loops, s).is_some() {
            continue;
        }
        let (cuts, leaves) = decompose_region(&loops, s);
        // A region whose pieces cannot be given counts together would fail the whole
        // cross-section: it stays as it is.
        if cuts.is_empty() || !assignable(&loops, &cuts, s) {
            continue;
        }
        // Do no harm: blocks of poor quality that are worse than the free mesh of the whole
        // region are not taken.
        let blocks = leaves.iter().filter_map(|l| score(l, s)).fold(f64::INFINITY, f64::min);
        if blocks.is_finite() && blocks < FAIR && blocks < free_score(&loops, s) - HARM_MARGIN {
            continue;
        }
        out.extend(cuts);
    }
    out.retain(|c| c.len() >= 2 && c.windows(2).all(|w| len(sub(w[1], w[0])) > s.tol) && c.iter().all(|p| p[0].is_finite() && p[1].is_finite()));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hex::arrangement;
    use crate::hex::section;
    use crate::size::UniformSize;
    use std::f64::consts::PI;

    /// A profile: loops of curves (polylines), the outer loop counter-clockwise, holes clockwise.
    type Profile = Vec<Vec<Vec<Q>>>;

    fn line(a: Q, b: Q) -> Vec<Q> {
        vec![a, b]
    }

    /// An arc about `c` of radius `r` from angle `a0` to `a1` (degrees), in steps of at most 5°.
    fn arc(c: Q, r: f64, a0: f64, a1: f64) -> Vec<Q> {
        let n = ((a1 - a0).abs() / 5.0).ceil().max(1.0) as usize;
        (0..=n)
            .map(|k| {
                // The mesher's portable trig, so every platform gets the same points (#556).
                let t = (a0 + (a1 - a0) * k as f64 / n as f64).to_radians();
                [c[0] + r * libm::cos(t), c[1] + r * libm::sin(t)]
            })
            .collect()
    }

    /// A closed loop through `pts` as straight curves.
    fn polygon(pts: &[Q]) -> Vec<Vec<Q>> {
        (0..pts.len()).map(|k| line(pts[k], pts[(k + 1) % pts.len()])).collect()
    }

    /// A circle as two half arcs, counter-clockwise or not.
    fn circle(c: Q, r: f64, ccw: bool) -> Vec<Vec<Q>> {
        if ccw { vec![arc(c, r, 0.0, 180.0), arc(c, r, 180.0, 360.0)] } else { vec![arc(c, r, 360.0, 180.0), arc(c, r, 180.0, 0.0)] }
    }

    /// The cross-section mesh of a profile and its extra (imprinted) curves at size `h`, with its
    /// regions decomposed: the mesh, which regions are filled and the arrangement.
    struct Meshed {
        section: section::Section,
        used: Vec<bool>,
        cuts: usize,
    }

    fn mesh_profile(profile: &Profile, extra: &[Vec<Q>], h: f64) -> Meshed {
        let mut polylines: Vec<(Vec<Q>, usize)> = Vec::new();
        let mut vertices: Vec<Q> = Vec::new();
        for lp in profile {
            for c in lp {
                vertices.push(c[0]);
                vertices.push(*c.last().unwrap());
                polylines.push((c.clone(), 0));
            }
        }
        for c in extra {
            vertices.push(c[0]);
            vertices.push(*c.last().unwrap());
            polylines.push((c.clone(), 1));
        }
        let polys: Vec<Poly> = profile
            .iter()
            .map(|lp| {
                let pts: Vec<Q> = lp.iter().flat_map(|c| c[..c.len() - 1].to_vec()).collect();
                Poly { end: vec![false; pts.len()], locked: vec![false; pts.len()], pts }
            })
            .collect();
        let tol = 1e-9;
        let field = UniformSize(h);
        let opts = QuadOptions::default();
        let settings = Settings { field: &field, tol, opts: &opts };
        let filled = |arr: &Arrangement| -> Vec<bool> {
            arr.regions.iter().map(|r| inside(&polys, super::super::layered::interior_point(arr, r).unwrap())).collect()
        };
        let mut arr = arrangement::build(&vertices, &polylines, tol, 0.0).unwrap();
        let mut used = filled(&arr);
        let c = cuts(&arr, &used, &vec![false; arr.curves.len()], &settings).unwrap();
        let cuts_added = c.len();
        if !c.is_empty() {
            polylines.extend(c.into_iter().map(|c| (c, CUT)));
            arr = arrangement::build(&vertices, &polylines, tol, 0.0).unwrap();
            used = filled(&arr);
        }
        let mut ia = Intervals::default();
        let pts = |c: usize| -> Vec<Q> { arr.curves[c].points.iter().map(|&p| arr.points[p]).collect() };
        for c in 0..arr.curves.len() {
            let l: f64 = pts(c).windows(2).map(|w| len(sub(w[1], w[0]))).sum();
            ia.var(l / h, 1);
        }
        let schemes = section::schemes(&arr, &used, &mut ia);
        let counts = ia.solve().unwrap();
        let place = |c: usize, n: usize| -> Result<Vec<Q>> { Ok(divide(&pts(c), n)) };
        let section = section::mesh(&arr, &used, schemes, &counts, &place, &field, &opts).unwrap();
        Meshed { section, used, cuts: cuts_added }
    }

    /// Every filled region structured, the quadrilaterals valid with corner scaled Jacobians of at
    /// least `min_sj`, conforming (each edge on one or two quadrilaterals, the ones on one along
    /// the profile, the others shared by two regions alike), and the area filled.
    fn check(m: &Meshed, profile: &Profile, min_sj: f64) -> f64 {
        let s = &m.section;
        for (r, scheme) in s.schemes.iter().enumerate() {
            if m.used[r] {
                assert_ne!(*scheme, Scheme::Free, "region {r} is free");
            }
        }
        let mut worst = f64::INFINITY;
        let mut area_sum = 0.0;
        let mut edges: std::collections::HashMap<(u32, u32), usize> = std::collections::HashMap::new();
        for q in s.quads.iter().flatten() {
            let c = q.map(|v| s.nodes[v as usize]);
            for i in 0..4 {
                worst = worst.min(planar_corner(&c, i).1.unwrap());
                let (a, b) = (q[i], q[(i + 1) % 4]);
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
            area_sum += 0.5 * (0..4).map(|k| cross(c[k], c[(k + 1) % 4])).sum::<f64>();
        }
        assert!(worst >= min_sj, "smallest corner scaled Jacobian {worst} below {min_sj}");
        let boundary: Vec<Q> = profile.iter().flatten().flatten().copied().collect();
        let on_profile = |p: Q| profile.iter().flatten().any(|c| c.windows(2).any(|w| seg_dist(p, w[0], w[1]) < 1e-6));
        for (&(a, b), &n) in &edges {
            assert!(n <= 2, "edge used {n} times");
            if n == 1 {
                let mid = [0.5 * (s.nodes[a as usize][0] + s.nodes[b as usize][0]), 0.5 * (s.nodes[a as usize][1] + s.nodes[b as usize][1])];
                assert!(on_profile(s.nodes[a as usize]) && on_profile(s.nodes[b as usize]), "an open edge inside the profile near {mid:?}");
            }
        }
        let _ = boundary;
        let profile_area: f64 = profile
            .iter()
            .map(|lp| {
                let pts: Vec<Q> = lp.iter().flat_map(|c| c[..c.len() - 1].to_vec()).collect();
                0.5 * (0..pts.len()).map(|k| cross(pts[k], pts[(k + 1) % pts.len()])).sum::<f64>()
            })
            .sum();
        // Curved boundaries are chords of the profile's polylines between the mesh nodes.
        assert!((area_sum - profile_area).abs() < 2e-3 * profile_area.abs(), "area {area_sum} against {profile_area}");
        worst
    }

    #[test]
    fn an_l_shape_is_cut_at_its_concave_corner() {
        let profile = vec![polygon(&[[0.0, 0.0], [40.0, 0.0], [40.0, 15.0], [15.0, 15.0], [15.0, 40.0], [0.0, 40.0]])];
        let m = mesh_profile(&profile, &[], 2.5);
        assert_eq!(m.cuts, 1);
        assert!(check(&m, &profile, 0.99) > 0.99);
    }

    #[test]
    fn a_plate_with_a_hole_gets_an_o_grid() {
        let mut profile = vec![polygon(&[[0.0, 0.0], [60.0, 0.0], [60.0, 40.0], [0.0, 40.0]])];
        profile.push(circle([30.0, 20.0], 8.0, false));
        let m = mesh_profile(&profile, &[], 2.5);
        assert_eq!(m.cuts, 4);
        // The O-grid's lines split the plate's corners (atan(20 / 30) = 34° of 90°).
        check(&m, &profile, 0.5);
    }

    /// A thin ring of air round a square pedestal (the Pi enclosure's, #397): the opening's sides
    /// carry on to the outer side, so the ring is side strips and corner blocks, all square,
    /// instead of a free mesh with slivers at the opening's corners. The outer sides carry curve
    /// ends off the opening's lines (cuts of the cells around), which the strips take in.
    #[test]
    fn a_ring_round_a_square_opening_continues_the_opening_sides() {
        let (x0, x1, y0, y1) = (29.94, 49.77, -49.0, -29.0);
        let outer = vec![
            line([x0, y0], [39.85, y0]),
            line([39.85, y0], [x1, y0]),
            line([x1, y0], [x1, -39.0]),
            line([x1, -39.0], [x1, y1]),
            line([x1, y1], [39.85, y1]),
            line([39.85, y1], [x0, y1]),
            line([x0, y1], [x0, -39.0]),
            line([x0, -39.0], [x0, y0]),
        ];
        let hole = polygon(&[[31.5, -30.5], [48.5, -30.5], [48.5, -47.5], [31.5, -47.5]]);
        let profile = vec![outer, hole];
        let m = mesh_profile(&profile, &[], 3.0);
        assert_eq!(m.cuts, 8);
        assert!(check(&m, &profile, 0.99) > 0.99);
    }

    #[test]
    fn a_plate_with_a_slot_gets_an_o_grid() {
        let mut profile = vec![polygon(&[[0.0, 0.0], [80.0, 0.0], [80.0, 40.0], [0.0, 40.0]])];
        // An obround slot: two R6 ends at x = 30 and 50 (clockwise).
        profile.push(vec![
            arc([30.0, 20.0], 6.0, 270.0, 90.0),
            line([30.0, 26.0], [50.0, 26.0]),
            arc([50.0, 20.0], 6.0, 90.0, -90.0),
            line([50.0, 14.0], [30.0, 14.0]),
        ]);
        let m = mesh_profile(&profile, &[], 2.5);
        check(&m, &profile, 0.5);
        // The intended O-grid on every platform (#556): Windows' last-bit-different arcs used to
        // give 10 regions and 2063 quads with a 0.375 corner.
        assert_eq!(m.section.schemes.len(), 9, "{:?}", m.section.schemes);
        assert_eq!(m.section.quads.iter().map(Vec::len).sum::<usize>(), 848);
    }

    #[test]
    fn a_plate_with_several_holes_is_split_between_them() {
        let mut profile = vec![polygon(&[[0.0, 0.0], [90.0, 0.0], [90.0, 40.0], [0.0, 40.0]])];
        for x in [20.0, 45.0, 70.0] {
            profile.push(circle([x, 20.0], 6.0, false));
        }
        let m = mesh_profile(&profile, &[], 2.0);
        check(&m, &profile, 0.5);
    }

    /// The half crank's cross-section (y ≤ 0) as the extrusion projects it: the outline, the
    /// bosses' circles where they rise above the arm, and the pocket's outline.
    fn crank_half() -> (Profile, Vec<Vec<Q>>) {
        let deg = |p: Q, c: Q| crate::vec3::atan2(p[1] - c[1], p[0] - c[0]).to_degrees();
        let spindle = [0.0, 0.0];
        let pedal = [170.0, 0.0];
        let j1 = [13.455063072993756, -18.65372021077189];
        let r1 = [15.471735240951828, -18.362928548901984];
        let r2 = [160.67227699330076, -12.384082712040556];
        let j2 = [161.58293237729853, -12.415835559271503];
        let c1 = [15.883151316213851, -28.354461805265395];
        let c2 = [160.91912663845798, -18.3790026658586];
        let a = |c: Q, r: f64, p: Q, q: Q, ccw: bool| {
            let (mut a0, mut a1) = (deg(p, c), deg(q, c));
            if ccw && a1 < a0 {
                a1 += 360.0;
            }
            if !ccw && a1 > a0 {
                a0 += 360.0;
            }
            let mut v = arc(c, r, a0, a1);
            v[0] = p;
            *v.last_mut().unwrap() = q;
            v
        };
        let outer = vec![
            a(spindle, 23.0, [-23.0, 0.0], j1, true),
            a(c1, libm::hypot(j1[0] - c1[0], j1[1] - c1[1]), j1, r1, false),
            line(r1, r2),
            a(c2, 6.0, r2, j2, false),
            a(pedal, 15.0, j2, [185.0, 0.0], true),
            line([185.0, 0.0], [177.0, 0.0]),
            a(pedal, 7.0, [177.0, 0.0], [163.0, 0.0], false),
            line([163.0, 0.0], [6.35, 0.0]),
            line([6.35, 0.0], [6.35, -4.85]),
            a([4.85, -4.85], 1.5, [6.35, -4.85], [4.85, -6.35], false),
            line([4.85, -6.35], [-4.85, -6.35]),
            a([-4.85, -4.85], 1.5, [-4.85, -6.35], [-6.35, -4.85], false),
            line([-6.35, -4.85], [-6.35, 0.0]),
            line([-6.35, 0.0], [-23.0, 0.0]),
        ];
        let extra = vec![
            a(spindle, 23.0, j1, [23.0, 0.0], true),
            a(pedal, 15.0, [155.0, 0.0], j2, true),
            line([48.0, 0.0], [48.0, -2.0]),
            a([53.0, -2.0], 5.0, [48.0, -2.0], [53.0, -7.0], true),
            line([53.0, -7.0], [125.0, -7.0]),
            a([125.0, -2.0], 5.0, [125.0, -7.0], [130.0, -2.0], true),
            line([130.0, -2.0], [130.0, 0.0]),
        ];
        (vec![outer], extra)
    }

    #[test]
    fn the_crank_half_section_is_all_blocks() {
        let (profile, extra) = crank_half();
        let m = mesh_profile(&profile, &extra, 3.0);
        // The arm meets the bosses about 20° off tangent: the corner there bounds the quality.
        let worst = check(&m, &profile, 0.37);
        assert!(m.cuts >= 4, "{}", m.cuts);
        let _ = worst;
    }

    /// One cut joining the T's two concave corners leaves two rectangles.
    #[test]
    fn a_t_shape_is_cut_between_its_concave_corners() {
        let profile = vec![polygon(&[[0.0, 0.0], [60.0, 0.0], [60.0, 10.0], [35.0, 10.0], [35.0, 40.0], [25.0, 40.0], [25.0, 10.0], [0.0, 10.0]])];
        let m = mesh_profile(&profile, &[], 2.0);
        assert_eq!(m.cuts, 1);
        check(&m, &profile, 0.99);
    }

    /// A cut never ends inside a curve that carries a seed: the L's long bottom side is locked,
    /// so its corner is cut the other way (along the bottom side's direction, ending on the
    /// unlocked left side).
    #[test]
    fn cuts_do_not_end_inside_seeded_curves() {
        let pts = [[0.0, 0.0], [40.0, 0.0], [40.0, 15.0], [15.0, 15.0], [15.0, 40.0], [0.0, 40.0]];
        let mut polylines: Vec<(Vec<Q>, usize)> = Vec::new();
        for k in 0..pts.len() {
            polylines.push((vec![pts[k], pts[(k + 1) % pts.len()]], 0));
        }
        let arr = arrangement::build(&pts, &polylines, 1e-9, 0.0).unwrap();
        let used: Vec<bool> = arr.regions.iter().map(|r| r.loops.len() == 1).collect();
        let field = UniformSize(2.5);
        let opts = QuadOptions::default();
        let settings = Settings { field: &field, tol: 1e-9, opts: &opts };
        let lies_on = |c: &arrangement::Curve, a: Q, b: Q| c.points.iter().all(|&p| seg_dist(arr.points[p], a, b) < 1e-9);
        for (a, b) in [([0.0, 0.0], [40.0, 0.0]), ([0.0, 40.0], [0.0, 0.0])] {
            let locked: Vec<bool> = arr.curves.iter().map(|c| lies_on(c, a, b)).collect();
            let c = cuts(&arr, &used, &locked, &settings).unwrap();
            assert_eq!(c.len(), 1);
            for p in [c[0][0], *c[0].last().unwrap()] {
                assert!(!(seg_dist(p, a, b) < 1e-9 && len(sub(p, a)) > 1e-9 && len(sub(p, b)) > 1e-9), "{c:?} ends inside the seeded side");
            }
        }
    }

    /// A loop of points, every one a curve end.
    fn corners(pts: &[Q]) -> Poly {
        Poly { pts: pts.to_vec(), end: vec![true; pts.len()], locked: vec![false; pts.len()] }
    }

    /// A cut that ends on an edge is admissible whichever side of the edge the rounding of its
    /// end point falls: the crank's O-grid lost a cut, and its pocket floor its blocks, when
    /// SolidCraft 0.16 moved the geometry by 1e-15 (#605).
    #[test]
    fn a_cut_ending_on_an_edge_does_not_read_the_rounding() {
        let square = [corners(&[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]])];
        for eps in [-1e-13, 0.0, 1e-13] {
            let a = End::Edge(0, 0, [5.0, eps]);
            let b = End::Edge(0, 2, [5.0, 10.0 - eps]);
            assert!(admissible(&square, a, b, 2.0), "{eps}");
        }
    }

    /// Mapped pieces round a hole chain their opposite sides; when one cut leaves the outer loop
    /// beside a corner instead of at it, the chain comes back to the side that also holds the
    /// short piece of outer side between them, which then needs no divisions: no counts exist,
    /// and the whole cross-section would fail (the Pi enclosure's base, #397).
    #[test]
    fn cuts_whose_pieces_cannot_be_given_counts_together_are_refused() {
        // A square plate with a round hole (16 points, clockwise; curve ends where cuts start).
        let at = |k: usize| {
            let a = -(k as f64) * PI / 8.0;
            [5.0 + 2.0 * a.cos(), 5.0 + 2.0 * a.sin()]
        };
        let hole: Vec<Q> = (0..16).map(at).collect();
        let mut end = vec![false; 16];
        for k in [2, 6, 7, 10, 14] {
            end[k] = true;
        }
        let outer = corners(&[[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [0.0, 0.5]]);
        let loops = vec![outer, Poly { locked: vec![false; 16], pts: hole, end }];
        let field = UniformSize(2.0);
        let opts = QuadOptions::default();
        let s = Settings { field: &field, tol: 1e-9, opts: &opts };
        let cut = |k: usize, b: Q| vec![at(k), b];
        // The O-grid: from the hole to the plate's corners.
        let o_grid = vec![cut(2, [10.0, 0.0]), cut(6, [0.0, 0.0]), cut(10, [0.0, 10.0]), cut(14, [10.0, 10.0])];
        assert!(assignable(&loops, &o_grid, &s));
        // One cut runs to the left side just above the corner, nearly along it: the piece below
        // takes it and the short piece of side as one side, opposite the cut beside it, which
        // the ring of pieces makes equal to the first.
        let beside = vec![cut(2, [10.0, 0.0]), cut(7, [0.0, 0.5]), cut(10, [0.0, 10.0]), cut(14, [10.0, 10.0])];
        assert!(!assignable(&loops, &beside, &s));
    }

    /// A strip one element across scores no better than its elements' aspect, and a piece much
    /// smaller than the target size no better than its elements against that size: their
    /// corners are square, their hexahedra slivers.
    #[test]
    fn thin_strips_and_small_pieces_score_low() {
        let field = UniformSize(2.5);
        let opts = QuadOptions::default();
        let s = Settings { field: &field, tol: 1e-9, opts: &opts };
        let rect = |w: f64, h: f64| vec![corners(&[[0.0, 0.0], [w, 0.0], [w, h], [0.0, h]])];
        let square = score(&rect(10.0, 10.0), &s).unwrap();
        let strip = score(&rect(10.0, 0.3), &s).unwrap();
        let small = score(&rect(0.3, 0.3), &s).unwrap();
        assert!(square > 0.99, "{square}");
        assert!(strip <= 4.0 * 0.3 / 2.5 + 1e-9, "{strip}");
        assert!(small <= 4.0 * 0.3 / 2.5 + 1e-9, "{small}");
    }

    /// A disk is meshed as its O-grid: a coarse tessellation whose curve ends turn like corners
    /// does not get it cut into wedges round its centre (they leave slivers at the rim).
    #[test]
    fn a_disk_is_not_cut() {
        let n = 10;
        let pts: Vec<Q> = (0..n).map(|k| k as f64 / n as f64 * 2.0 * PI).map(|a| [1.645 * a.cos(), 1.645 * a.sin()]).collect();
        let mut end = vec![false; n];
        for k in [0, 3, 5, 8] {
            end[k] = true;
        }
        let disk = vec![Poly { locked: vec![false; n], pts, end }];
        let field = UniformSize(2.5);
        let opts = QuadOptions::default();
        let s = Settings { field: &field, tol: 1e-9, opts: &opts };
        let (cuts, _) = decompose_region(&disk, &s);
        assert!(cuts.is_empty(), "{cuts:?}");
    }

    #[test]
    fn hostile_settings_are_an_error() {
        let arr = arrangement::build(&[], &[(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]], 0)], 1e-9, 0.0).unwrap();
        let field = UniformSize(1.0);
        let opts = QuadOptions::default();
        assert!(cuts(&arr, &[true], &[], &Settings { field: &field, tol: f64::NAN, opts: &opts }).is_err());
        let _ = PI;
    }
}
