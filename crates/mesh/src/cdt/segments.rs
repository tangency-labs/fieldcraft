//! Segment recovery (Si & Gärtner 2005, as restated by Diazzi et al. 2023, §3.3): while an edge
//! of the input surface (a *segment*) is not an edge of the Delaunay tetrahedralization, split
//! it at a Steiner point and insert the point (Bowyer–Watson). Each Steiner point is the LNC
//! point `a + t (b − a)` of its original segment `a b` (Diazzi et al. §4.2–§4.3): it lies on the
//! segment exactly, and only `t` is computed in floating point, so the placement rule may round
//! but never moves a point off its segment.
//!
//! The placement follows Si's reference-point rule. The reference point `r` of a missing
//! subsegment `v₁ v₂` is the vertex inside its diametral sphere for which the circle through
//! `v₁, v₂, r` is largest. A vertex is *acute* when two segments meet there at less than 90°;
//! the rule depends on how many of the original segment's endpoints are acute:
//! - none: split where the smaller of the spheres centred at `v₁`, `v₂` through `r` cuts the
//!   segment, or at the midpoint when both radii exceed half the length;
//! - one (`w`): `p` is where the sphere centred at `w` through `r` cuts the segment; keep `p` if
//!   it is closer to `r` than to the far endpoint `v_o`, else shift it towards the near endpoint
//!   `v_w` by `|r p|` when that is under half of `|v_w p|`, else take the midpoint of `v_w p`;
//! - two: split at the midpoint first, after which each half has one acute endpoint.
//!
//! Any value that falls outside the open subsegment becomes the midpoint. Termination is Si's
//! (the spheres around acute vertices keep splits from cascading); the number of Steiner points
//! and of passes is still capped.

use crate::predicates::{Axis, GenericPoint, Lnc, Point3, generic_orient2d};
use crate::vec3::{cross, dot, norm, sub};
use crate::{Error, Result};

use super::Recovery;

/// Most tetrahedra searched for a reference point.
const NEAR_CAP: usize = 20_000;
/// Most passes over the segments.
const MAX_PASSES: usize = 10_000;

/// An input segment and its split points.
#[derive(Debug, Clone)]
pub(super) struct Segment {
    /// The original endpoints (explicit vertices), `ends[0] < ends[1]`.
    pub ends: [u32; 2],
    /// `(t, vertex)` along the segment from `ends[0]` (`t = 0`) to `ends[1]` (`t = 1`).
    pub splits: Vec<(f64, u32)>,
    /// Whether each original endpoint is acute.
    pub acute: [bool; 2],
}

impl Segment {
    pub(super) fn new(a: u32, b: u32) -> Self {
        let ends = [a.min(b), a.max(b)];
        Segment { ends, splits: vec![(0.0, ends[0]), (1.0, ends[1])], acute: [false; 2] }
    }

    /// The vertices from `from` to the other end, in order.
    pub(super) fn chain_from(&self, from: u32) -> Vec<u32> {
        let mut v: Vec<u32> = self.splits.iter().map(|&(_, v)| v).collect();
        if from != self.ends[0] {
            v.reverse();
        }
        v
    }
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(sub(a, b))
}

impl Recovery {
    /// Mark acute vertices: two segments at a vertex make an angle under 90° (floating point:
    /// this only steers where Steiner points go).
    pub(super) fn classify_acute(&mut self) -> Result<()> {
        let n = self.cx.points().len();
        let mut dirs: Vec<Vec<[f64; 3]>> = vec![Vec::new(); n];
        for s in &self.segments {
            let [a, b] = s.ends;
            let (pa, pb) = (self.cx.approx(a)?, self.cx.approx(b)?);
            if let Some(d) = dirs.get_mut(a as usize) {
                d.push(sub(pb, pa));
            }
            if let Some(d) = dirs.get_mut(b as usize) {
                d.push(sub(pa, pb));
            }
        }
        let acute: Vec<bool> = dirs.iter().map(|d| d.iter().enumerate().any(|(i, x)| d.iter().skip(i + 1).any(|y| dot(*x, *y) > 0.0))).collect();
        for s in &mut self.segments {
            s.acute = s.ends.map(|v| acute.get(v as usize).copied().unwrap_or(false));
        }
        Ok(())
    }

    /// Split missing segments until every one is an edge of the (Delaunay) tetrahedralization.
    /// Alg. 2 sweeps all segments until none is missing; an insertion only changes edges between
    /// vertices of the tetrahedra it deletes, so after the first sweep only the segments through
    /// those vertices are checked again.
    pub(super) fn recover_segments(&mut self) -> Result<()> {
        let mut scratch = Vec::new();
        let mut through: Vec<Vec<u32>> = vec![Vec::new(); self.cx.points().len()];
        for (s, seg) in self.segments.iter().enumerate() {
            for v in seg.ends {
                if let Some(l) = through.get_mut(v as usize) {
                    l.push(s as u32);
                }
            }
        }
        let mut queued = vec![false; self.segments.len()];
        let mut work: Vec<usize> = (0..self.segments.len()).collect();
        for _ in 0..MAX_PASSES {
            if work.is_empty() {
                return Ok(());
            }
            queued.fill(false);
            let mut next: Vec<usize> = Vec::new();
            let mut requeue = |s: usize, next: &mut Vec<usize>| {
                if let Some(q) = queued.get_mut(s)
                    && !*q
                {
                    *q = true;
                    next.push(s);
                }
            };
            for &s in &work {
                let mut i = 0;
                while let (Some(&(tu, u)), Some(&(tv, v))) = (self.segment(s)?.splits.get(i), self.segment(s)?.splits.get(i + 1)) {
                    if self.cx.has_edge(u, v, &mut scratch)? {
                        i += 1;
                        continue;
                    }
                    self.stats.missing_subsegments += 1;
                    let t = self.split_parameter(s, (tu, u), (tv, v))?;
                    let [a, b] = self.segment(s)?.ends;
                    let (pa, pb) = (explicit(self.cx.pt(a)?)?, explicit(self.cx.pt(b)?)?);
                    let w = self.cx.push_point(GenericPoint::Lnc(Lnc::new(&pa, &pb, t)?))?;
                    through.push(vec![s as u32]);
                    self.steiner.push((s, t));
                    if self.steiner.len() > self.options.max_steiner {
                        let at = self.cx.approx(u)?;
                        return Err(Error::limit(crate::message![
                            format!("more than {} Steiner points on boundary edges (the last near ", self.options.max_steiner),
                            crate::message::point(at),
                            "): the boundary has features much thinner than its triangles are long; raise the cap or tessellate coarser"
                        ]));
                    }
                    let touched = self.cx.insert(w, u)?;
                    if let Some(seg) = self.segments.get_mut(s) {
                        seg.splits.insert(i + 1, (t, w));
                    }
                    requeue(s, &mut next);
                    for x in touched {
                        for &r in through.get(x as usize).map(Vec::as_slice).unwrap_or_default() {
                            requeue(r as usize, &mut next);
                        }
                    }
                    // Both halves wait for the next sweep (Alg. 2 splits each missing segment
                    // once per sweep).
                    i += 2;
                }
            }
            work = next;
        }
        Err(Error::Limit(format!("segment recovery did not finish in {MAX_PASSES} passes")))
    }

    pub(super) fn segment(&self, s: usize) -> Result<&Segment> {
        self.segments.get(s).ok_or_else(|| super::complex::corrupt("segment id"))
    }

    /// Si's split parameter for the missing subsegment `(tu, u)`–`(tv, v)` of segment `s`.
    fn split_parameter(&mut self, s: usize, (tu, u): (f64, u32), (tv, v): (f64, u32)) -> Result<f64> {
        let (ends, acute) = {
            let seg = self.segment(s)?;
            (seg.ends, seg.acute)
        };
        let [a, b] = ends;
        let (pa, pb) = (self.cx.approx(a)?, self.cx.approx(b)?);
        let length = dist(pa, pb);
        let mid = 0.5 * (tu + tv);
        let r = self.reference_point(u, v)?;
        let chosen = match (r, acute) {
            (None, _) => mid,
            // Two acute ends: the first split is the midpoint.
            (_, [true, true]) if tu < 0.5 && tv > 0.5 => 0.5,
            (Some(r), acute) => {
                let pr = self.cx.approx(r)?;
                let (pu, pv) = (self.cx.approx(u)?, self.cx.approx(v)?);
                match acute {
                    [false, false] => {
                        let (ru, rv) = (dist(pu, pr), dist(pv, pr));
                        let half = 0.5 * dist(pu, pv);
                        if ru > half && rv > half {
                            mid
                        } else if ru <= rv {
                            tu + ru / length
                        } else {
                            tv - rv / length
                        }
                    }
                    _ => {
                        // The acute end w: the only one, or the end of the half this
                        // subsegment lies in.
                        let w_is_a = match acute {
                            [true, false] => true,
                            [false, true] => false,
                            _ => tv <= 0.5,
                        };
                        let pw = if w_is_a { pa } else { pb };
                        let radius = dist(pw, pr);
                        let tp = if w_is_a { radius / length } else { 1.0 - radius / length };
                        let (t_near, t_far) = if w_is_a { (tu, tv) } else { (tv, tu) };
                        let at = |t: f64| -> [f64; 3] { std::array::from_fn(|k| pa[k] + t * (pb[k] - pa[k])) };
                        let p = at(tp);
                        if dist(p, pr) < dist(p, at(t_far)) {
                            tp
                        } else {
                            let d = dist(pr, p);
                            let near = dist(at(t_near), p);
                            if d < 0.5 * near { if t_near < tp { tp - d / length } else { tp + d / length } } else { 0.5 * (t_near + tp) }
                        }
                    }
                }
            }
        };
        let t = if chosen.is_finite() && chosen > tu && chosen < tv { chosen } else { mid };
        if !(t > tu && t < tv) {
            let at = self.cx.approx(u)?;
            return Err(Error::Measured(
                crate::Kind::Failed,
                crate::message!["a missing segment near ", crate::message::point(at), " is too short to split further in double precision"],
            ));
        }
        Ok(t)
    }

    /// The reference point of the missing subsegment `u v`: the vertex strictly inside its
    /// diametral sphere with the largest circumcircle through `u, v` (floating point; `None` if
    /// rounding hides every encroaching vertex). Candidates are the vertices of the tetrahedra
    /// whose boxes meet the subsegment's box, connected to `u`'s star: they include every
    /// tetrahedron the subsegment crosses, whose vertices are the ones keeping it from being an
    /// edge. (On the corpus this picks the same point as searching the whole diametral ball, at a
    /// fraction of the cost.) A vertex exactly on the open subsegment makes the
    /// input invalid (a segment through a vertex), which is reported.
    fn reference_point(&mut self, u: u32, v: u32) -> Result<Option<u32>> {
        let (pu, pv) = (self.cx.approx(u)?, self.cx.approx(v)?);
        let c: [f64; 3] = std::array::from_fn(|k| 0.5 * (pu[k] + pv[k]));
        let r2 = 0.25 * dot(sub(pu, pv), sub(pu, pv));
        let lo: [f64; 3] = std::array::from_fn(|k| pu[k].min(pv[k]));
        let hi: [f64; 3] = std::array::from_fn(|k| pu[k].max(pv[k]));
        let near = self.cx.tets_near(u, lo, hi, NEAR_CAP)?;
        let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut best: Option<(f64, u32)> = None;
        for t in near {
            for x in self.cx.tet(t)? {
                if x == u || x == v || !seen.insert(x) {
                    continue;
                }
                let px = self.cx.approx(x)?;
                let d = sub(px, c);
                if dot(d, d) >= r2 {
                    continue;
                }
                let (ea, eb) = (sub(pu, px), sub(pv, px));
                let area2 = norm(cross(ea, eb));
                let radius = if area2 > 0.0 { norm(ea) * norm(eb) * norm(sub(ea, eb)) / (2.0 * area2) } else { f64::INFINITY };
                if !radius.is_finite() && self.on_segment(u, v, x)? {
                    let at = self.cx.approx(x)?;
                    return Err(Error::Measured(
                        crate::Kind::Invalid,
                        crate::message![
                            format!("vertex {x} at "),
                            crate::message::point(at),
                            " lies on a boundary edge: the surface is not a valid complex"
                        ],
                    ));
                }
                if best.is_none_or(|(b, _)| radius > b) {
                    best = Some((radius, x));
                }
            }
        }
        Ok(best.map(|(_, x)| x))
    }

    /// Whether `x` is exactly collinear with `u v` (it is inside their diametral ball, so then on
    /// the open segment).
    fn on_segment(&self, u: u32, v: u32, x: u32) -> Result<bool> {
        let (a, b, c) = (self.cx.pt(u)?, self.cx.pt(v)?, self.cx.pt(x)?);
        Ok(Axis::ALL.into_iter().all(|ax| generic_orient2d(a, b, c, ax).is_zero()))
    }
}

/// The explicit point of an original segment endpoint.
fn explicit(p: &GenericPoint) -> Result<Point3> {
    p.explicit().copied().ok_or_else(|| super::complex::corrupt("a segment endpoint is not an input vertex"))
}
