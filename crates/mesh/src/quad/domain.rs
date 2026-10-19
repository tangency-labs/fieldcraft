//! Validation of a planar region given as boundary loops: finite, admissible coordinates, no
//! zero-length edges, no repeated points, no crossing or touching edges, holes inside the outer
//! loop and not inside each other. Every geometric decision is an exact `orient2d` or an exact
//! coordinate comparison.

use std::cmp::Ordering;

use crate::predicates::{Point2, Sign, orient2d, snap_to_domain};
use crate::{Error, Result};

/// The most boundary points of one region.
pub const MAX_BOUNDARY_POINTS: usize = 1 << 20;
/// The most loops (outer plus holes) of one region.
pub const MAX_LOOPS: usize = 1 << 14;
/// The most segment-cell insertions and pair tests of the intersection check.
const MAX_CHECK_WORK: usize = 1 << 28;

/// A validated region: the points in input order and the loops as point ids, the outer loop
/// first and counterclockwise, holes clockwise.
#[derive(Debug, Clone)]
pub(crate) struct Region {
    pub pts: Vec<[f64; 2]>,
    pub p2: Vec<Point2>,
    pub loops: Vec<Vec<u32>>,
}

impl Region {
    pub fn point(&self, v: u32) -> Result<&Point2> {
        self.p2.get(v as usize).ok_or_else(|| Error::Failed(format!("point {v} out of range")))
    }
}

/// The bounding box `(lo, hi)` of points (`lo > hi` when there are none).
pub(crate) fn bbox(pts: impl IntoIterator<Item = [f64; 2]>) -> ([f64; 2], [f64; 2]) {
    pts.into_iter()
        .fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |(lo, hi), p| ([lo[0].min(p[0]), lo[1].min(p[1])], [hi[0].max(p[0]), hi[1].max(p[1])]))
}

/// Admissible exact-predicate point from user coordinates (finite, within the model limit; tiny
/// magnitudes rounded to zero).
pub(crate) fn admissible(c: [f64; 2]) -> Result<[f64; 2]> {
    fieldcraft_model::check_coordinates(&c)?;
    Ok([snap_to_domain(c[0]), snap_to_domain(c[1])])
}

/// Validate loops given as point lists (the first is the outer boundary, the rest are holes, in
/// either orientation). Point ids number the points in input order, loop by loop.
pub(crate) fn region(input: &[Vec<[f64; 2]>]) -> Result<Region> {
    if input.is_empty() {
        return Err(Error::Invalid("a region needs an outer boundary loop".into()));
    }
    if input.len() > MAX_LOOPS {
        return Err(Error::Limit(format!("{} loops exceed the limit of {MAX_LOOPS}", input.len())));
    }
    let total: usize = input.iter().map(Vec::len).sum();
    if total > MAX_BOUNDARY_POINTS {
        return Err(Error::Limit(format!("{total} boundary points exceed the limit of {MAX_BOUNDARY_POINTS}")));
    }
    let mut pts = Vec::with_capacity(total);
    let mut loops = Vec::with_capacity(input.len());
    for (k, l) in input.iter().enumerate() {
        if l.len() < 3 {
            return Err(Error::Invalid(format!("loop {k} has {} points; a loop needs at least 3", l.len())));
        }
        let mut ids = Vec::with_capacity(l.len());
        for (i, c) in l.iter().enumerate() {
            let c = admissible(*c).map_err(|e| Error::Invalid(format!("loop {k} point {i}: {e}")))?;
            ids.push(pts.len() as u32);
            pts.push(c);
        }
        loops.push(ids);
    }
    let p2 = pts.iter().map(|c| Point2::new(c[0], c[1])).collect::<Result<Vec<_>>>()?;
    let mut r = Region { pts, p2, loops };
    check_distinct(&r)?;
    check_simple(&r)?;
    orient(&mut r)?;
    check_nesting(&r)?;
    Ok(r)
}

fn lex(a: &[f64; 2], b: &[f64; 2]) -> Ordering {
    a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1]))
}

/// Report a point by its loop and position.
fn locate(r: &Region, v: u32) -> String {
    for (k, l) in r.loops.iter().enumerate() {
        if let Some(i) = l.iter().position(|&x| x == v) {
            return format!("loop {k} point {i}");
        }
    }
    format!("point {v}")
}

/// No point appears twice (this also rejects zero-length edges).
fn check_distinct(r: &Region) -> Result<()> {
    let mut order: Vec<u32> = (0..r.pts.len() as u32).collect();
    order.sort_by(|&a, &b| match (r.pts.get(a as usize), r.pts.get(b as usize)) {
        (Some(p), Some(q)) => lex(p, q).then(a.cmp(&b)),
        _ => a.cmp(&b),
    });
    for w in order.windows(2) {
        if let [a, b] = *w
            && r.pts.get(a as usize) == r.pts.get(b as usize)
        {
            let what = if r.loops.iter().any(|l| l.windows(2).any(|e| e == [a, b] || e == [b, a]) || is_closing(l, a, b)) {
                "zero-length edge"
            } else {
                "repeated point"
            };
            let at = r.pts.get(a as usize).copied().map_or(crate::message::Part::from("an unknown point"), crate::message::point);
            return Err(Error::invalid(crate::message![format!("{what}: {} and {} coincide at ", locate(r, a), locate(r, b)), at]));
        }
    }
    Ok(())
}

fn is_closing(l: &[u32], a: u32, b: u32) -> bool {
    matches!((l.first(), l.last()), (Some(&f), Some(&z)) if (f == a && z == b) || (f == b && z == a))
}

/// The boundary segments `(a, b)` of every loop, in loop order.
pub(crate) fn segments(loops: &[Vec<u32>]) -> Vec<(u32, u32)> {
    let mut s = Vec::new();
    for l in loops {
        for (i, &a) in l.iter().enumerate() {
            if let Some(&b) = l.get((i + 1) % l.len()) {
                s.push((a, b));
            }
        }
    }
    s
}

/// `c` lies within the bounding box of `a`, `b` (for collinear points: on the segment).
fn on_segment(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    a[0].min(b[0]) <= c[0] && c[0] <= a[0].max(b[0]) && a[1].min(b[1]) <= c[1] && c[1] <= a[1].max(b[1])
}

/// Whether two boundary segments conflict: they cross or touch, or (when they share an end
/// point) they fold back on each other.
fn segments_conflict(r: &Region, s: (u32, u32), t: (u32, u32)) -> Result<bool> {
    let (a, b, c, d) = (r.point(s.0)?, r.point(s.1)?, r.point(t.0)?, r.point(t.1)?);
    let shared = [(s.0 == t.0, b, d, a), (s.0 == t.1, b, c, a), (s.1 == t.0, a, d, b), (s.1 == t.1, a, c, b)];
    if let Some(&(_, x, y, o)) = shared.iter().find(|e| e.0) {
        // Segments o–x and o–y: a conflict when collinear and on the same side of o.
        if orient2d(o, x, y) != Sign::Zero {
            return Ok(false);
        }
        let (o, x, y) = (o.coords(), x.coords(), y.coords());
        let dot = (x[0] - o[0]) * (y[0] - o[0]) + (x[1] - o[1]) * (y[1] - o[1]);
        return Ok(dot > 0.0);
    }
    let (o1, o2, o3, o4) = (orient2d(a, b, c), orient2d(a, b, d), orient2d(c, d, a), orient2d(c, d, b));
    let proper = (o1 * o2) == Sign::Negative && (o3 * o4) == Sign::Negative;
    let (pa, pb, pc, pd) = (a.coords(), b.coords(), c.coords(), d.coords());
    Ok(proper
        || (o1.is_zero() && on_segment(pa, pb, pc))
        || (o2.is_zero() && on_segment(pa, pb, pd))
        || (o3.is_zero() && on_segment(pc, pd, pa))
        || (o4.is_zero() && on_segment(pc, pd, pb)))
}

/// No two boundary segments cross, touch or overlap. Segments are bucketed on a uniform grid
/// over the bounding box; each pair is tested once, in the first cell both bounding boxes cover.
pub(crate) fn check_simple(r: &Region) -> Result<()> {
    let segs = segments(&r.loops);
    let n = segs.len();
    let (lo, hi) = bbox(r.pts.iter().copied());
    let g = ((n as f64).sqrt().ceil() as usize).clamp(1, 1024);
    let cell = |x: f64, axis: usize| -> usize {
        let (l, h) = (lo[axis], hi[axis]);
        if !(h > l) {
            return 0;
        }
        (((x - l) / (h - l) * g as f64) as usize).min(g - 1)
    };
    let mut boxes = Vec::with_capacity(n);
    let mut grid: Vec<Vec<u32>> = vec![Vec::new(); g * g];
    let mut work = 0usize;
    for (i, &(a, b)) in segs.iter().enumerate() {
        let (pa, pb) = (r.point(a)?.coords(), r.point(b)?.coords());
        let bx = [cell(pa[0].min(pb[0]), 0), cell(pa[1].min(pb[1]), 1), cell(pa[0].max(pb[0]), 0), cell(pa[1].max(pb[1]), 1)];
        for cy in bx[1]..=bx[3] {
            for cx in bx[0]..=bx[2] {
                work += 1;
                if let Some(c) = grid.get_mut(cy * g + cx) {
                    c.push(i as u32);
                }
            }
        }
        if work > MAX_CHECK_WORK {
            return Err(Error::Limit("boundary intersection check exceeded its work limit".into()));
        }
        boxes.push(bx);
    }
    for (ci, members) in grid.iter().enumerate() {
        let (cx, cy) = (ci % g, ci / g);
        for (k, &i) in members.iter().enumerate() {
            for &j in members.get(k + 1..).unwrap_or(&[]) {
                let (Some(bi), Some(bj), Some(&si), Some(&sj)) =
                    (boxes.get(i as usize), boxes.get(j as usize), segs.get(i as usize), segs.get(j as usize))
                else {
                    continue;
                };
                if bi[0].max(bj[0]) != cx || bi[1].max(bj[1]) != cy {
                    continue;
                }
                work += 1;
                if work > MAX_CHECK_WORK {
                    return Err(Error::Limit("boundary intersection check exceeded its work limit".into()));
                }
                if segments_conflict(r, si, sj)? {
                    return Err(Error::Invalid(format!(
                        "boundary edges cross, touch or overlap: {}–{} and {}–{}",
                        locate(r, si.0),
                        locate(r, si.1),
                        locate(r, sj.0),
                        locate(r, sj.1)
                    )));
                }
            }
        }
    }
    Ok(())
}

/// The exact orientation of a simple loop: the sign of the turn at its lexicographically
/// smallest point, a vertex of its convex hull.
pub(crate) fn loop_orientation(r: &Region, l: &[u32]) -> Result<Sign> {
    let (i, _) = l
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| r.pts.get(v as usize).map(|p| (i, p)))
        .min_by(|x, y| lex(x.1, y.1))
        .ok_or_else(|| Error::Invalid("empty loop".into()))?;
    let n = l.len();
    let at = |k: usize| l.get(k % n).copied().ok_or_else(|| Error::Failed("loop index".into()));
    Ok(orient2d(r.point(at(i + n - 1)?)?, r.point(at(i)?)?, r.point(at(i + 1)?)?))
}

/// The outer loop is counterclockwise and the holes clockwise (no repair).
pub(crate) fn check_orientation(r: &Region) -> Result<()> {
    for (k, l) in r.loops.iter().enumerate() {
        let want = if k == 0 { Sign::Positive } else { Sign::Negative };
        if loop_orientation(r, l)? != want {
            return Err(Error::Invalid(format!("loop {k} has the wrong orientation")));
        }
    }
    Ok(())
}

/// Outer loop counterclockwise, holes clockwise.
fn orient(r: &mut Region) -> Result<()> {
    for k in 0..r.loops.len() {
        let l = r.loops.get(k).cloned().unwrap_or_default();
        let want = if k == 0 { Sign::Positive } else { Sign::Negative };
        let s = loop_orientation(r, &l)?;
        if s == Sign::Zero {
            return Err(Error::Invalid(format!("loop {k} is degenerate (all points collinear)")));
        }
        if s != want
            && let Some(l) = r.loops.get_mut(k)
        {
            l.reverse();
        }
    }
    Ok(())
}

/// Winding number of loop `l` around `p` (exact: half-open crossing rule with `orient2d`).
pub(crate) fn winding(r: &Region, l: &[u32], p: &Point2) -> Result<i64> {
    let y = p.coords()[1];
    let mut w = 0i64;
    for (i, &a) in l.iter().enumerate() {
        let b = l.get((i + 1) % l.len()).copied().unwrap_or(a);
        let (pa, pb) = (r.point(a)?, r.point(b)?);
        let (ya, yb) = (pa.coords()[1], pb.coords()[1]);
        if ya <= y && yb > y && orient2d(pa, pb, p) == Sign::Positive {
            w += 1;
        } else if yb <= y && ya > y && orient2d(pa, pb, p) == Sign::Negative {
            w -= 1;
        }
    }
    Ok(w)
}

/// Every hole lies inside the outer loop and outside every other hole. The loops do not
/// intersect, so one point per hole decides.
pub(crate) fn check_nesting(r: &Region) -> Result<()> {
    let holes = r.loops.len().saturating_sub(1);
    if holes.saturating_mul(r.pts.len()) > MAX_CHECK_WORK {
        return Err(Error::Limit(format!("{holes} holes on {} boundary points exceed the nesting-check limit", r.pts.len())));
    }
    let outer = r.loops.first().cloned().unwrap_or_default();
    for (k, hole) in r.loops.iter().enumerate().skip(1) {
        let first = hole.first().copied().ok_or_else(|| Error::Invalid(format!("loop {k} is empty")))?;
        let p = r.point(first)?;
        if winding(r, &outer, p)? == 0 {
            return Err(Error::Invalid(format!("hole loop {k} lies outside the outer loop")));
        }
        for (m, other) in r.loops.iter().enumerate().skip(1) {
            if m != k && winding(r, other, p)? != 0 {
                return Err(Error::Invalid(format!("hole loop {k} lies inside hole loop {m}")));
            }
        }
    }
    Ok(())
}
