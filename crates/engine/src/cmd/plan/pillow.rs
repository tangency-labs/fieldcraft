//! The ring-and-pillow layout of a mesh plan (#455), for a body swept along the plan's direction
//! whose section the grid layout ([`super::blocks`]) does not fit: sides at an angle to the axes
//! and arcs where an arm runs into its bosses, as on a crank arm. The body is read in the local
//! frame of [`super::blocks`] (third axis along the sweep):
//! - **Levels.** The planar faces square to the sweep cut the layout into slabs. One section (the
//!   slab with the most material) is the reference: its outline, and every other slab's section
//!   is made of the same walls (a pocket's floor carves grid cells, a boss's top stands above the
//!   arm on its ring alone).
//! - **Rings.** Each boss (a butterfly ring of the plan) is an annulus of blocks from the hole
//!   inside it (the circle of a concentric O-grid ring, or a hole of the section inside the
//!   circle: a bore, with a block corner at the middle of each of its corner fillets) out to the
//!   boss's circle, through every level, so a boss's top stands on the blocks the arm joins.
//! - **Carved grid.** The rest of the section (the arm) is a grid of planes square to the in-plane
//!   axes, carved to the cells a margin inside the outline and outside every boss: the planar faces
//!   square to an axis (a pocket's sides, whose corner fillets the projection meets with the grid
//!   corner), the planes a margin off each boss's circle, and planes every few element sizes
//!   (closer across the arm than along it, so the carved boundary follows a taper).
//! - **Pillow.** A layer of blocks between the carved grid's boundary and the wall each boundary
//!   face faces: the outline, or a boss's circle where the arm runs into it (there the pillow's outer
//!   vertices are the annulus's). This is pillowing (Mitchell and Tautges 1995) between a carved grid
//!   and the geometry, as the grid-based methods project their boundary layer (Schneiders 1996).
//!   Each boundary vertex goes along the mean normal of its boundary faces to its wall, the
//!   landings kept in order along the wall and apart: a corner whose landing would pass or crowd
//!   a neighbour's goes between its neighbours' in proportion to its two faces, inside its fan. A
//!   grid corner whose faces face a boss and the outline lands on each and gets a corner block
//!   reaching the kink where the outline leaves the boss's circle: a junction that meets the
//!   circle at an angle (an arm's root fillet that is not tangent to the boss).
//!
//! Each slab keeps the blocks whose middle lies inside the body (generalized winding number, as
//! [`super::blocks`] carves). A body this layout does not fit is reported with the reason; the plan
//! then keeps its frame cuts and swept cells.
//!
//! References: S. A. Mitchell, T. J. Tautges, *Pillowing doublets: refining a mesh to ensure that
//! faces share at most one edge*, Proc. 4th International Meshing Roundtable, 231–240, 1995;
//! R. Schneiders, *A grid-based algorithm for the generation of hexahedral element meshes*,
//! Engineering with Computers 12:168–177, 1996.

use std::collections::BTreeMap;

use fieldcraft_geometry::FaceSurface;
use fieldcraft_mesh::project::Topology;
use fieldcraft_model::RingKind;
use fieldcraft_model::units::Quantity;

use super::blocks::{Input, Layout, PERM, corner_det, winding};

type V2 = [f64; 2];
type V3 = [f64; 3];

/// Grid spacing along the carved region's longer extent, in element sizes.
const SPACING_ALONG: f64 = 2.5;
/// Grid spacing across it, in element sizes: the carved boundary steps with a taper.
const SPACING_ACROSS: f64 = 1.25;
/// The smallest distance of a kept grid cell's corner from a wall, in element sizes: the pillow's
/// least thickness.
const MARGIN: f64 = 0.6;
/// The largest arc of a boss's circle one annulus block spans where the boss is the outline (30°);
/// otherwise the arc the arm's landings on it span on average.
const MAX_OWN_ARC: f64 = std::f64::consts::PI / 6.0;
/// The smallest such arc (8°).
const MIN_OWN_ARC: f64 = std::f64::consts::PI / 22.5;
/// The least angle at which the outline leaves a boss's circle: a tangent junction leaves a sliver
/// no block corner takes.
const MIN_KINK: f64 = 10.0 * std::f64::consts::PI / 180.0;
/// How much further than the nearest wall a boundary face's normal may meet the wall it faces.
const FACING_REACH: f64 = 2.0;
/// The least width of a grid cell beside a plane that bounds a corner fillet, in fillet radii.
const FILLET_CLEARANCE: f64 = 0.7;
/// Most grid cells before carving.
const MAX_CELLS: usize = 20_000;
/// Most blocks of the layout.
const MAX_BLOCKS: usize = 20_000;
/// Most points of one section loop.
const MAX_LOOP: usize = 200_000;
/// The least distance between neighbouring landings, as a fraction of the boundary edge between
/// their vertices.
const LANDING_GAP: f64 = 0.5;
/// How far inside its fan a corner's landing stays: its faces' normals tilted toward each other
/// by atan(0.2), 11°.
const FAN_INSET: f64 = 0.2;
/// Rounds of the ordering of the landings.
const ORDER_ROUNDS: usize = 500;
/// How much the layout's quads may miss a slab's section by, as a fraction of its area.
const AREA_TOLERANCE: f64 = 0.03;

fn sub(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}

fn add(a: V2, b: V2) -> V2 {
    [a[0] + b[0], a[1] + b[1]]
}

fn scale(a: V2, s: f64) -> V2 {
    [a[0] * s, a[1] * s]
}

fn dot(a: V2, b: V2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn cross(a: V2, b: V2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

fn norm(a: V2) -> f64 {
    a[0].hypot(a[1])
}

fn unit(a: V2) -> Option<V2> {
    let n = norm(a);
    (n.is_finite() && n > 1e-300).then(|| scale(a, 1.0 / n))
}

/// The closest point to `p` on segment `a b`.
fn near_segment(p: V2, a: V2, b: V2) -> V2 {
    let e = sub(b, a);
    let l = dot(e, e);
    let t = if l > 0.0 { (dot(sub(p, a), e) / l).clamp(0.0, 1.0) } else { 0.0 };
    add(a, scale(e, t))
}

/// A closed loop of a section: its points, counterclockwise for an outline and clockwise for a
/// hole (material on the left), and the B-rep face of each segment `i → i + 1`.
#[derive(Debug, Clone)]
struct Loop {
    pts: Vec<V2>,
    face: Vec<u32>,
}

impl Loop {
    fn segment(&self, i: usize) -> (V2, V2) {
        let n = self.pts.len();
        (self.pts.get(i % n).copied().unwrap_or_default(), self.pts.get((i + 1) % n).copied().unwrap_or_default())
    }

    fn area(&self) -> f64 {
        (0..self.pts.len()).map(|i| self.segment(i)).map(|(a, b)| cross(a, b)).sum::<f64>() / 2.0
    }

    /// Even-odd point in polygon.
    fn contains(&self, p: V2) -> bool {
        let mut odd = false;
        for i in 0..self.pts.len() {
            let (a, b) = self.segment(i);
            if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) {
                odd = !odd;
            }
        }
        odd
    }

    fn distance(&self, p: V2) -> f64 {
        (0..self.pts.len()).map(|i| self.segment(i)).map(|(a, b)| norm(sub(p, near_segment(p, a, b)))).fold(f64::INFINITY, f64::min)
    }

    /// The first point where the ray `p + t d` (t > `eps`) meets the loop: (t, point).
    fn ray(&self, p: V2, d: V2, eps: f64) -> Option<(f64, V2)> {
        let mut best: Option<(f64, V2)> = None;
        for i in 0..self.pts.len() {
            let (a, b) = self.segment(i);
            let e = sub(b, a);
            let den = cross(d, e);
            if den.abs() < 1e-300 {
                continue;
            }
            let w = sub(a, p);
            let (t, u) = (cross(w, e) / den, cross(w, d) / den);
            if t > eps && (-1e-12..=1.0 + 1e-12).contains(&u) && best.is_none_or(|x| t < x.0) {
                best = Some((t, add(p, scale(d, t))));
            }
        }
        best
    }

    /// The cumulative length at each point (the last entry the perimeter).
    fn lengths(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.pts.len() + 1);
        out.push(0.0);
        for i in 0..self.pts.len() {
            let (a, b) = self.segment(i);
            let last = out.last().copied().unwrap_or(0.0);
            out.push(last + norm(sub(b, a)));
        }
        out
    }
}

/// A boss of the plan in the local frame: its circle, the hole inside it (the inner loop of its
/// annulus) and the points of that hole that must be block corners (its corner fillets' middles).
#[derive(Debug, Clone)]
struct Boss {
    c: V2,
    r: f64,
    inner: Loop,
    corners: Vec<V2>,
}

/// A wall a boundary face of the carved grid faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Wall {
    Outline,
    Boss(usize),
}

/// Slice the tessellation at local height `z`: the section's loops. Within a face the triangles
/// share points, so its crossings chain by mesh edge; the pieces of different faces chain by their
/// ends (each face carries its own copy of the points it shares, within `join`).
fn slice(inp: &Input, local: &dyn Fn(V3) -> V3, z: f64, join: f64) -> Result<Vec<Loop>, String> {
    let length = |v: f64| (inp.show)(Quantity::Length, v);
    let m = inp.mesh;
    // Segments as (start key, end key, start, end, face); a key is a mesh edge (its two points).
    let mut segs: Vec<([u32; 2], [u32; 2], V2, V2, u32)> = Vec::new();
    for (t, &f) in m.triangles.iter().zip(&m.tri_face) {
        let Some(p) = t.iter().map(|&i| m.positions.get(i as usize).map(|q| local(*q))).collect::<Option<Vec<V3>>>() else { continue };
        let (Some(&a), Some(&b), Some(&c)) = (p.first(), p.get(1), p.get(2)) else { continue };
        let above = [a[2] > z, b[2] > z, c[2] > z];
        if above.iter().all(|x| *x) || above.iter().all(|x| !*x) {
            continue;
        }
        let mut cuts: Vec<([u32; 2], V2)> = Vec::with_capacity(2);
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            if above[i] != above[j] {
                let (pi, pj) = (p[i], p[j]);
                let s = (z - pi[2]) / (pj[2] - pi[2]);
                let key = [t[i].min(t[j]), t[i].max(t[j])];
                cuts.push((key, [pi[0] + s * (pj[0] - pi[0]), pi[1] + s * (pj[1] - pi[1])]));
            }
        }
        let [(ka, qa), (kb, qb)] = cuts.as_slice() else { continue };
        // Material on the left: the segment runs along the outward normal turned +90°.
        let n = {
            let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
            [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2]]
        };
        let along = [-n[1], n[0]];
        if dot(sub(*qb, *qa), along) >= 0.0 {
            segs.push((*ka, *kb, *qa, *qb, f));
        } else {
            segs.push((*kb, *ka, *qb, *qa, f));
        }
    }
    if segs.len() > MAX_LOOP {
        return Err(format!("the section at {} has {} segments", length(z), segs.len()));
    }
    // Pieces: chains within a face, by mesh edge.
    let mut by_start: BTreeMap<(u32, [u32; 2]), usize> = BTreeMap::new();
    for (i, s) in segs.iter().enumerate() {
        by_start.insert((s.4, s.0), i);
    }
    let mut used = vec![false; segs.len()];
    let mut pieces: Vec<Loop> = Vec::new();
    // Start from segments no other segment of the face ends at, then from the rest (closed loops).
    let mut ends: std::collections::BTreeSet<(u32, [u32; 2])> = std::collections::BTreeSet::new();
    for s in &segs {
        ends.insert((s.4, s.1));
    }
    let order: Vec<usize> = (0..segs.len()).filter(|&i| !ends.contains(&(segs[i].4, segs[i].0))).chain(0..segs.len()).collect();
    for i in order {
        if used.get(i).copied().unwrap_or(true) {
            continue;
        }
        let mut piece = Loop { pts: Vec::new(), face: Vec::new() };
        let mut k = i;
        while let Some(s) = segs.get(k) {
            if used[k] {
                break;
            }
            used[k] = true;
            piece.pts.push(s.2);
            piece.face.push(s.4);
            match by_start.get(&(s.4, s.1)) {
                Some(&next) if !used[next] => k = next,
                _ => {
                    piece.pts.push(s.3);
                    break;
                }
            }
        }
        pieces.push(piece);
    }
    // Loops: a piece's end joins the piece starting nearest it (within `join`).
    let mut loops = Vec::new();
    let mut taken = vec![false; pieces.len()];
    for start in 0..pieces.len() {
        if taken[start] {
            continue;
        }
        taken[start] = true;
        let mut lp = Loop { pts: Vec::new(), face: Vec::new() };
        let mut cur = start;
        for _ in 0..=pieces.len() {
            let Some(p) = pieces.get(cur) else { break };
            let n = p.face.len();
            lp.pts.extend(p.pts.iter().take(n));
            lp.face.extend(&p.face);
            let end = p.pts.last().copied().unwrap_or_default();
            let first = lp.pts.first().copied().unwrap_or_default();
            if norm(sub(end, first)) <= join {
                break;
            }
            let next = (0..pieces.len())
                .filter(|&j| !taken[j])
                .filter_map(|j| pieces.get(j).and_then(|q| q.pts.first()).map(|q| (norm(sub(*q, end)), j)))
                .filter(|(d, _)| *d <= join)
                .min_by(|a, b| a.0.total_cmp(&b.0));
            match next {
                Some((_, j)) => {
                    taken[j] = true;
                    cur = j;
                }
                None => return Err(format!("the section at {} does not close", length(z))),
            }
        }
        if lp.pts.len() >= 3 {
            loops.push(lp);
        }
    }
    Ok(loops)
}

/// The grid planes over `[lo, hi]`: the fixed ones inside, and between them (and the ends) planes
/// `step` apart at most; a fixed plane bounding a corner fillet of radius `r` keeps the planes
/// beside it `r` away.
fn planes(lo: f64, hi: f64, fixed: &[(f64, f64)], step: f64, eps: f64) -> Vec<f64> {
    let mut pts: Vec<(f64, f64)> = vec![(lo, 0.0), (hi, 0.0)];
    pts.extend(fixed.iter().copied().filter(|(v, _)| *v > lo + eps && *v < hi - eps));
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    pts.dedup_by(|a, b| {
        if (a.0 - b.0).abs() < eps {
            b.1 = b.1.max(a.1);
            true
        } else {
            false
        }
    });
    let mut out = Vec::new();
    for w in pts.windows(2) {
        let [(a, ra), (b, rb)] = [w[0], w[1]];
        let len = b - a;
        let mut n = (len / step).ceil().clamp(1.0, MAX_CELLS as f64) as usize;
        // A filleted corner's grid corner moves 0.29 of the radius along each axis onto the
        // fillet's middle: the cells beside it stay 0.7 of the radius wide.
        while n > 1 && len / n as f64 + 1e-12 < FILLET_CLEARANCE * ra.max(rb) {
            n -= 1;
        }
        out.extend((0..n).map(|k| a + len * k as f64 / n as f64));
    }
    out.push(hi);
    out
}

/// The parameter of point `q` along the closed polyline with cumulative `lengths`, and the point
/// at a parameter.
struct Param<'a> {
    lp: &'a Loop,
    cum: Vec<f64>,
}

impl Param<'_> {
    fn perimeter(&self) -> f64 {
        self.cum.last().copied().unwrap_or(0.0)
    }

    fn of(&self, q: V2) -> f64 {
        let mut best = (f64::INFINITY, 0.0);
        for i in 0..self.lp.pts.len() {
            let (a, b) = self.lp.segment(i);
            let h = near_segment(q, a, b);
            let d = norm(sub(q, h));
            if d < best.0 {
                best = (d, self.cum.get(i).copied().unwrap_or(0.0) + norm(sub(h, a)));
            }
        }
        best.1
    }

    fn at(&self, s: f64) -> V2 {
        let per = self.perimeter();
        if !(per > 0.0) {
            return self.lp.pts.first().copied().unwrap_or_default();
        }
        let s = s.rem_euclid(per);
        let i = self.cum.partition_point(|c| *c <= s).saturating_sub(1).min(self.lp.pts.len().saturating_sub(1));
        let (a, b) = self.lp.segment(i);
        let (c0, c1) = (self.cum.get(i).copied().unwrap_or(0.0), self.cum.get(i + 1).copied().unwrap_or(per));
        let t = if c1 > c0 { ((s - c0) / (c1 - c0)).clamp(0.0, 1.0) } else { 0.0 };
        add(a, scale(sub(b, a), t))
    }
}

/// The first point where the ray `p + t d` meets the circle about `c` of radius `r`.
fn ray_circle(p: V2, d: V2, c: V2, r: f64) -> Option<(f64, V2)> {
    let w = sub(p, c);
    let b = dot(w, d);
    let disc = b * b - (dot(w, w) - r * r);
    if disc < 0.0 {
        return None;
    }
    [-b - disc.sqrt(), -b + disc.sqrt()].into_iter().find(|t| *t > 0.0).map(|t| (t, add(p, scale(d, t))))
}

/// Whether the quad's corners turn one way (convex, not degenerate).
fn convex(q: &[V2; 4]) -> bool {
    let turns: Vec<f64> = (0..4).map(|t| cross(sub(q[(t + 1) % 4], q[t]), sub(q[(t + 2) % 4], q[(t + 1) % 4]))).collect();
    turns.iter().all(|x| *x > 0.0) || turns.iter().all(|x| *x < 0.0)
}

/// Plan the ring-and-pillow block topology of the body.
// Several per-axis and per-loop arrays are indexed together.
#[allow(clippy::needless_range_loop)]
pub(crate) fn layout(inp: &Input) -> Result<Layout, String> {
    let length = |v: f64| (inp.show)(Quantity::Length, v);
    let area_shown = |v: f64| (inp.show)(Quantity::Area, v);
    let d = inp.direction;
    let k = (0..3).find(|&k| (d[k].abs() - 1.0).abs() < 1e-9).ok_or("the sweep direction is not along an axis")?;
    let perm = PERM[k];
    let local = |p: V3| -> V3 { [p[perm[0]], p[perm[1]], p[perm[2]]] };
    let world = |l: V3| -> V3 {
        let mut w = [0.0; 3];
        for a in 0..3 {
            w[perm[a]] = l[a];
        }
        w
    };
    if inp.shell.is_some() {
        return Err("a shelled body takes the grid layout's wall layers".into());
    }
    if inp.mesh.triangles.is_empty() {
        return Err("the body has no tessellation".into());
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &inp.mesh.positions {
        let q = local(*p);
        for a in 0..3 {
            lo[a] = lo[a].min(q[a]);
            hi[a] = hi[a].max(q[a]);
        }
    }
    let diag = (0..3).map(|a| (hi[a] - lo[a]).powi(2)).sum::<f64>().sqrt();
    if !(diag.is_finite() && diag > 0.0 && inp.size.is_finite() && inp.size > 0.0) {
        return Err("the body has no extent".into());
    }
    let eps = 1e-6 * diag;
    let name = |f: usize| inp.shape.faces().get(f).map_or_else(|| format!("#{f}"), |x| x.name.clone());
    let size = inp.size;
    let margin = MARGIN * size;

    // Levels: the planar faces square to the sweep; every other face along it.
    let mut levels = vec![lo[2], hi[2]];
    // Planar faces square to an in-plane axis: (axis, offset, face).
    let mut sides: Vec<(usize, f64, usize)> = Vec::new();
    for (f, s) in inp.surfaces.iter().enumerate() {
        match s {
            Some(FaceSurface::Plane { normal, point }) => {
                let (n, q) = (local(*normal), local(*point));
                if (n[2].abs() - 1.0).abs() < 1e-6 {
                    levels.push(q[2]);
                } else if n[2].abs() > 1e-6 {
                    return Err(format!("face `{}` is a plane at an angle to the sweep", name(f)));
                } else if let Some(a) = (0..2).find(|&a| (n[a].abs() - 1.0).abs() < 1e-6) {
                    sides.push((a, q[a], f));
                }
            }
            Some(FaceSurface::Cylinder { axis, .. }) if (local(*axis)[2].abs() - 1.0).abs() < 1e-6 => {}
            _ => return Err(format!("face `{}` is not prismatic along the sweep", name(f))),
        }
    }
    levels.sort_by(f64::total_cmp);
    levels.dedup_by(|a, b| (*a - *b).abs() < eps);
    if levels.len() < 2 || levels.len() > 1001 {
        return Err(format!("{} levels", levels.len()));
    }

    // The sections of the slabs, and the reference: the one with the most material.
    let join = (4.0 * inp.mesh.weld_distance).max(1e-6 * diag);
    let mut sections: Vec<Vec<Loop>> = Vec::with_capacity(levels.len() - 1);
    for w in levels.windows(2) {
        sections.push(slice(inp, &local, 0.5 * (w[0] + w[1]), join)?);
    }
    let area = |s: &[Loop]| s.iter().map(Loop::area).sum::<f64>();
    let reference = (0..sections.len()).max_by(|a, b| area(&sections[*a]).total_cmp(&area(&sections[*b]))).ok_or("no slab")?;
    let refs = sections.get(reference).ok_or("no reference section")?;
    let outlines: Vec<&Loop> = refs.iter().filter(|l| l.area() > 0.0).collect();
    let [outline] = outlines.as_slice() else {
        return Err(format!("the widest section has {} outlines", outlines.len()));
    };
    let outline: &Loop = outline;
    let holes: Vec<&Loop> = refs.iter().filter(|l| l.area() < 0.0).collect();

    // The bosses: the plan's butterfly rings, each with the hole its annulus goes round.
    let ring_at = |r: &fieldcraft_model::PlanRing| -> Result<V2, String> {
        if (r.axis.iter().zip(&d).map(|(a, b)| a * b).sum::<f64>().abs() - 1.0).abs() > 1e-6 {
            return Err(format!("the ring of F{} is not along the sweep", r.feature));
        }
        let c = local(r.center);
        Ok([c[0], c[1]])
    };
    let mut bosses: Vec<Boss> = Vec::new();
    for r in inp.rings.iter().filter(|r| r.kind == RingKind::Butterfly) {
        let c = ring_at(r)?;
        if bosses.iter().any(|b| norm(sub(b.c, c)) < eps) {
            continue;
        }
        let inside: Vec<&&Loop> = holes.iter().filter(|h| h.pts.iter().all(|p| norm(sub(*p, c)) < r.radius - eps)).collect();
        let [inner] = inside.as_slice() else {
            return Err(format!("the boss of radius {} holds {} holes (this layout rings a boss round one)", length(r.radius), inside.len()));
        };
        // Its corners: the middles of the runs of blend faces (corner fillets).
        let mut corners = Vec::new();
        let n = inner.pts.len();
        let blend = |i: usize| inner.face.get(i % n).is_some_and(|f| inp.blends.contains(f));
        if (0..n).any(|i| !blend(i)) {
            let first = (0..n).find(|&i| !blend(i)).unwrap_or(0);
            let mut i = first + 1;
            while i <= first + n {
                if blend(i) {
                    let start = i;
                    while blend(i) && i < start + n {
                        i += 1;
                    }
                    let run: Vec<V2> = (start..=i).map(|j| inner.pts.get(j % n).copied().unwrap_or_default()).collect();
                    let cum: Vec<f64> = run
                        .windows(2)
                        .scan(0.0, |acc, w| {
                            *acc += norm(sub(w[1], w[0]));
                            Some(*acc)
                        })
                        .collect();
                    let half = cum.last().copied().unwrap_or(0.0) / 2.0;
                    let j = cum.partition_point(|c| *c < half);
                    corners.push(run.get(j + 1).copied().unwrap_or_default());
                }
                i += 1;
            }
        }
        bosses.push(Boss { c, r: r.radius, inner: (**inner).clone(), corners });
    }
    if bosses.is_empty() {
        return Err("no boss to ring".into());
    }
    for r in inp.rings.iter().filter(|r| r.kind == RingKind::OGrid) {
        let c = ring_at(r)?;
        if !bosses.iter().any(|b| norm(sub(b.c, c)) < 1e-6 * b.r.max(1.0) && r.radius < b.r) {
            return Err(format!("the hole of radius {} is not inside a boss", length(r.radius)));
        }
    }
    for h in &holes {
        if !bosses.iter().any(|b| b.inner.pts.first() == h.pts.first()) {
            return Err("a hole of the section is not inside a boss".into());
        }
    }
    let on_boss = |b: &Boss, p: V2| (norm(sub(p, b.c)) - b.r).abs() < 1e-4 * b.r.max(1.0) + join;
    // Each boss's circle runs along the outline (its own side).
    for b in &bosses {
        if !outline.pts.iter().any(|p| on_boss(b, *p)) {
            return Err(format!("the boss of radius {} is not on the outline", length(b.r)));
        }
    }

    // The kinks: where the outline leaves a boss's circle, at an angle.
    let n_out = outline.pts.len();
    // A segment of the outline on a boss: on a cylindrical face about the boss's axis, of its radius.
    let seg_on = |bi: usize, i: usize| {
        let (Some(b), Some(&f)) = (bosses.get(bi), outline.face.get(i)) else { return false };
        match inp.surfaces.get(f as usize) {
            Some(Some(FaceSurface::Cylinder { origin, radius, .. })) => {
                let o = local(*origin);
                (radius - b.r).abs() < 1e-6 * b.r.max(1.0) && norm(sub([o[0], o[1]], b.c)) < 1e-6 * b.r.max(1.0)
            }
            _ => false,
        }
    };
    let mut kinks: Vec<Vec<V2>> = vec![Vec::new(); bosses.len()];
    for (bi, b) in bosses.iter().enumerate() {
        for i in 0..n_out {
            let prev = (i + n_out - 1) % n_out;
            if seg_on(bi, prev) == seg_on(bi, i) {
                continue;
            }
            let p = outline.pts.get(i).copied().unwrap_or_default();
            // The outline's direction off the circle against the circle's tangent there.
            let (a0, a1) = if seg_on(bi, i) { outline.segment(prev) } else { outline.segment(i) };
            let dir = unit(sub(a1, a0)).ok_or("a degenerate outline segment")?;
            let radial = unit(sub(p, b.c)).ok_or("a kink at a boss's centre")?;
            let angle = dot(dir, radial).abs().clamp(0.0, 1.0).asin();
            if angle < MIN_KINK {
                return Err(format!("the outline meets the boss of radius {} tangentially (at {:.1}°)", length(b.r), angle.to_degrees()));
            }
            kinks[bi].push(p);
        }
    }

    // The carved region's extent: the outline outside the bosses.
    let (mut rlo, mut rhi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in &outline.pts {
        if bosses.iter().all(|b| norm(sub(*p, b.c)) > b.r + eps) {
            for a in 0..2 {
                rlo[a] = rlo[a].min(p[a]);
                rhi[a] = rhi[a].max(p[a]);
            }
        }
    }
    if !(0..2).all(|a| rhi[a] - rlo[a] > 2.0 * margin) {
        return Err("no arm between the bosses".into());
    }
    let along = if rhi[0] - rlo[0] >= rhi[1] - rlo[1] { 0 } else { 1 };
    // Fixed planes: the planar faces square to an in-plane axis in the region (not a boss's
    // hole), each with the largest corner fillet it bounds, and the planes a margin off each
    // boss's circle.
    let fillet_of = |a: usize, v: f64| -> f64 {
        inp.blends
            .iter()
            .filter_map(|&f| match inp.surfaces.get(f as usize) {
                Some(Some(FaceSurface::Cylinder { origin, radius, .. })) => {
                    let o = local(*origin);
                    ((o[a] - v).abs() - radius).abs().lt(&(1e-4 * radius.max(1.0))).then_some(*radius)
                }
                _ => None,
            })
            .fold(0.0, f64::max)
    };
    let mut fixed: [Vec<(f64, f64)>; 2] = [Vec::new(), Vec::new()];
    for &(a, v, f) in &sides {
        let pts: Vec<V2> = inp
            .mesh
            .triangles
            .iter()
            .zip(&inp.mesh.tri_face)
            .filter(|(_, tf)| **tf as usize == f)
            .flat_map(|(t, _)| t.iter().filter_map(|&i| inp.mesh.positions.get(i as usize).map(|p| local(*p))))
            .map(|p| [p[0], p[1]])
            .collect();
        if pts.iter().all(|p| bosses.iter().any(|b| norm(sub(*p, b.c)) < b.r)) {
            continue;
        }
        fixed[a].push((v, fillet_of(a, v)));
    }
    // Along the arm, where it runs into the bosses.
    for b in &bosses {
        for s in [-1.0, 1.0] {
            fixed[along].push((b.c[along] + s * (b.r + margin), 0.0));
        }
    }
    let spacing = |a: usize| if a == along { SPACING_ALONG * size } else { SPACING_ACROSS * size };
    let gp: [Vec<f64>; 2] = std::array::from_fn(|a| planes(rlo[a], rhi[a], &fixed[a], spacing(a), eps));
    let (nu, nv) = (gp[0].len() - 1, gp[1].len() - 1);
    if nu * nv > MAX_CELLS {
        return Err(format!("a grid of {nu}×{nv} cells is too large for one layout"));
    }
    let gpt = |i: usize, j: usize| -> V2 { [gp[0].get(i).copied().unwrap_or_default(), gp[1].get(j).copied().unwrap_or_default()] };

    // The kept cells: every corner inside the outline, `margin` off it and off every boss.
    let clear = |p: V2| outline.contains(p) && outline.distance(p) >= margin - eps && bosses.iter().all(|b| norm(sub(p, b.c)) >= b.r + margin - eps);
    let corner_ok: Vec<Vec<bool>> = (0..=nu).map(|i| (0..=nv).map(|j| clear(gpt(i, j))).collect()).collect();
    let ok = |i: usize, j: usize| corner_ok.get(i).and_then(|c| c.get(j)).copied().unwrap_or(false);
    let mut kept: std::collections::BTreeSet<(usize, usize)> = std::collections::BTreeSet::new();
    for i in 0..nu {
        for j in 0..nv {
            if ok(i, j) && ok(i + 1, j) && ok(i, j + 1) && ok(i + 1, j + 1) {
                kept.insert((i, j));
            }
        }
    }
    // Two cells meeting at a corner alone: the one with fewer neighbours goes.
    for _ in 0..=nu * nv {
        let mut changed = false;
        for i in 1..nu {
            for j in 1..nv {
                let q = [kept.contains(&(i - 1, j - 1)), kept.contains(&(i, j - 1)), kept.contains(&(i, j)), kept.contains(&(i - 1, j))];
                let pair = match q {
                    [true, false, true, false] => [(i - 1, j - 1), (i, j)],
                    [false, true, false, true] => [(i, j - 1), (i - 1, j)],
                    _ => continue,
                };
                let neighbours = |c: (usize, usize)| {
                    [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)]
                        .iter()
                        .filter(|(a, b)| kept.contains(&((c.0 as i64 + a) as usize, (c.1 as i64 + b) as usize)))
                        .count()
                };
                let gone = if neighbours(pair[0]) <= neighbours(pair[1]) { pair[0] } else { pair[1] };
                kept.remove(&gone);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    if kept.is_empty() {
        return Err("the arm has no grid cell a margin inside its walls".into());
    }

    // The carved boundary: one loop of grid edges, the cells on its left.
    let mut next: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
    for &(i, j) in &kept {
        let sides = [
            ((i, j), (i + 1, j), j.checked_sub(1).map(|jj| (i, jj))),
            ((i + 1, j), (i + 1, j + 1), Some((i + 1, j))),
            ((i + 1, j + 1), (i, j + 1), Some((i, j + 1))),
            ((i, j + 1), (i, j), i.checked_sub(1).map(|ii| (ii, j))),
        ];
        for (a, b, nb) in sides {
            if nb.is_none_or(|c| !kept.contains(&c)) && next.insert(a, b).is_some() {
                return Err("the carved grid pinches".into());
            }
        }
    }
    let start = next.keys().next().copied().ok_or("no carved boundary")?;
    let mut ring: Vec<(usize, usize)> = vec![start];
    for _ in 0..next.len() {
        let Some(&b) = ring.last().and_then(|v| next.get(v)) else { break };
        if b == start {
            break;
        }
        ring.push(b);
    }
    if ring.len() != next.len() {
        return Err("the carved grid has more than one boundary".into());
    }
    let n = ring.len();
    let pts: Vec<V2> = ring.iter().map(|&(i, j)| gpt(i, j)).collect();
    let at = |k: usize| pts.get(k % n).copied().unwrap_or_default();
    let normal = |k: usize| -> V2 {
        let e = unit(sub(at(k + 1), at(k))).unwrap_or([1.0, 0.0]);
        [e[1], -e[0]]
    };

    // A boss's circle as a wall: where the arm runs into it (its points inside the outline).
    let internal = |q: V2| outline.contains(q) && outline.distance(q) > join;
    let boss_distance = |bi: usize, p: V2| -> f64 {
        let Some(b) = bosses.get(bi) else { return f64::INFINITY };
        let q = add(b.c, scale(unit(sub(p, b.c)).unwrap_or([1.0, 0.0]), b.r));
        if internal(q) { (norm(sub(p, b.c)) - b.r).abs() } else { f64::INFINITY }
    };
    let wall_distance = |w: Wall, p: V2| match w {
        Wall::Outline => outline.distance(p),
        Wall::Boss(bi) => boss_distance(bi, p),
    };
    // The wall each boundary edge faces: the first its normal meets (a boss's circle where the
    // arm runs into it, or the outline), unless that is much further than the nearest wall (a
    // step of the carved boundary facing a slanted side): then the nearest.
    let walls: Vec<Wall> = (0..n)
        .map(|k| {
            let mid = scale(add(at(k), at(k + 1)), 0.5);
            let nk = normal(k);
            let mut candidates: Vec<(Wall, f64)> =
                (0..bosses.len()).map(Wall::Boss).chain([Wall::Outline]).map(|w| (w, wall_distance(w, mid))).collect();
            candidates.sort_by(|a, b| a.1.total_cmp(&b.1));
            let nearest = candidates.first().copied().unwrap_or((Wall::Outline, 0.0));
            let mut first: Option<(Wall, f64)> = outline.ray(mid, nk, eps).map(|h| (Wall::Outline, h.0));
            for (bi, b) in bosses.iter().enumerate() {
                if let Some((t, q)) = ray_circle(mid, nk, b.c, b.r)
                    && internal(q)
                    && first.is_none_or(|f| t < f.1)
                {
                    first = Some((Wall::Boss(bi), t));
                }
            }
            match first {
                Some((w, t)) if t <= FACING_REACH * nearest.1 => w,
                _ => nearest.0,
            }
        })
        .collect();
    let wall = |k: usize| walls.get(k % n).copied().unwrap_or(Wall::Outline);

    // The landings: parameters along the walls (arc length on the outline, angle on a boss).
    let along_outline = Param { lp: outline, cum: outline.lengths() };
    let param = |w: Wall, p: V2| -> f64 {
        match w {
            Wall::Outline => along_outline.of(p),
            Wall::Boss(b) => bosses.get(b).map_or(0.0, |b| (p[1] - b.c[1]).atan2(p[0] - b.c[0])),
        }
    };
    let point = |w: Wall, s: f64| -> V2 {
        match w {
            Wall::Outline => along_outline.at(s),
            Wall::Boss(b) => bosses.get(b).map_or([0.0; 2], |b| add(b.c, [b.r * s.cos(), b.r * s.sin()])),
        }
    };
    let period = |w: Wall| match w {
        Wall::Outline => along_outline.perimeter(),
        Wall::Boss(_) => std::f64::consts::TAU,
    };
    // A vertex's landing on wall `w` along `dir`; the nearest point of the wall when the ray
    // misses it or meets it much further than that.
    let hit = |w: Wall, p: V2, dir: V2| -> Option<V2> {
        let ray = match w {
            Wall::Outline => outline.ray(p, dir, eps).map(|h| h.1),
            Wall::Boss(b) => bosses.get(b).and_then(|b| ray_circle(p, dir, b.c, b.r)).map(|h| h.1),
        };
        let reach = FACING_REACH * wall_distance(w, p);
        match ray {
            Some(q) if norm(sub(q, p)) <= reach => Some(q),
            _ => match w {
                Wall::Outline => Some(along_outline.at(along_outline.of(p))),
                Wall::Boss(b) => bosses.get(b).map(|b| add(b.c, scale(unit(sub(p, b.c)).unwrap_or([1.0, 0.0]), b.r))),
            },
        }
    };
    let junction: Vec<bool> = (0..n).map(|k| wall(k + n - 1) != wall(k)).collect();
    // Per vertex: its landing's parameter on the wall of the edge before it and after it (one
    // landing, unless a junction).
    let (mut pin, mut pout) = (vec![0.0; n], vec![0.0; n]);
    let mut fixed_landing = vec![false; n];
    for k in 0..n {
        let (w0, w1) = (wall(k + n - 1), wall(k));
        let v = at(k);
        if junction[k] {
            let a = hit(w0, v, normal(k + n - 1)).ok_or("a corner's ray misses its wall")?;
            let b = hit(w1, v, normal(k)).ok_or("a corner's ray misses its wall")?;
            pin[k] = param(w0, a);
            pout[k] = param(w1, b);
            fixed_landing[k] = true;
        } else {
            let (n0, n1) = (normal(k + n - 1), normal(k));
            let dir = unit(add(n0, n1)).unwrap_or(n1);
            let q = hit(w1, v, dir).ok_or("a boundary vertex's ray misses its wall")?;
            pin[k] = param(w1, q);
            pout[k] = pin[k];
            fixed_landing[k] = cross(n0, n1).abs() < 1e-12;
        }
    }
    let unwrap = |x: f64, near: f64, w: Wall| {
        let per = period(w);
        x - per * ((x - near) / per).round()
    };
    // A corner's landing stays inside its corner's fan: between the rays along its two faces'
    // normals, each tilted toward the other by atan(FAN_INSET).
    let fan: Vec<Option<(f64, f64)>> = (0..n)
        .map(|k| {
            if fixed_landing[k] {
                return None;
            }
            let (n0, n1) = (normal(k + n - 1), normal(k));
            let w = wall(k);
            let turn = |a: V2, b: V2| unit(add(a, scale(b, FAN_INSET)));
            let (d0, d1) = (turn(n0, n1)?, turn(n1, n0)?);
            let ray = |dir: V2| match w {
                Wall::Outline => outline.ray(at(k), dir, eps).map(|h| h.1),
                Wall::Boss(b) => bosses.get(b).and_then(|b| ray_circle(at(k), dir, b.c, b.r)).map(|h| h.1),
            };
            let (p0, p1) = (unwrap(param(w, ray(d0)?), pout[k], w), unwrap(param(w, ray(d1)?), pout[k], w));
            Some((p0.min(p1), p0.max(p1)))
        })
        .collect();
    // A corner's landing goes between its neighbours' in proportion to its two boundary edges'
    // lengths (within its fan) when it passes a neighbour's, or comes within half an edge's
    // length of it (a step of the carved boundary, whose corners' bisectors land almost on the
    // next vertex's).
    let length_of = |w: Wall| match w {
        Wall::Outline => 1.0,
        Wall::Boss(b) => bosses.get(b).map_or(1.0, |b| b.r),
    };
    for _ in 0..ORDER_ROUNDS {
        let mut changed = false;
        for k in 0..n {
            if fixed_landing[k] {
                continue;
            }
            let w = wall(k);
            let a = unwrap(pout[(k + n - 1) % n], pout[k], w);
            let b = unwrap(pin[(k + 1) % n], pout[k], w);
            let (la, lb) = (norm(sub(at(k), at(k + n - 1))), norm(sub(at(k + 1), at(k))));
            let scale_w = length_of(w);
            let cramped = (pout[k] - a).abs() * scale_w < LANDING_GAP * la || (b - pout[k]).abs() * scale_w < LANDING_GAP * lb;
            let mut target = a + (b - a) * la / (la + lb).max(1e-300);
            if let Some(Some((lo, hi))) = fan.get(k) {
                let (lo, hi) = (unwrap(*lo, target, w), unwrap(*hi, target, w));
                target = target.clamp(lo.min(hi), lo.max(hi));
            }
            if (!(a.min(b) < pout[k] && pout[k] < a.max(b)) || cramped) && (target - pout[k]).abs() * scale_w > eps {
                pout[k] = target;
                pin[k] = target;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // The 2D vertices and quads.
    let mut v2: Vec<V2> = Vec::new();
    let mut grid_id: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut gid = |i: usize, j: usize, v2: &mut Vec<V2>| -> usize {
        *grid_id.entry((i, j)).or_insert_with(|| {
            v2.push(gpt(i, j));
            v2.len() - 1
        })
    };
    let mut quads: Vec<[usize; 4]> = Vec::new();
    for &(i, j) in &kept {
        let q = [gid(i, j, &mut v2), gid(i + 1, j, &mut v2), gid(i + 1, j + 1, &mut v2), gid(i, j + 1, &mut v2)];
        quads.push(q);
    }
    let base: Vec<usize> = ring.iter().map(|&(i, j)| gid(i, j, &mut v2)).collect();
    // Landings: one point per vertex, two at a junction.
    let mut land_in = vec![0usize; n];
    let mut land_out = vec![0usize; n];
    for k in 0..n {
        let w0 = wall(k + n - 1);
        v2.push(point(w0, pin[k]));
        land_in[k] = v2.len() - 1;
        if junction[k] {
            v2.push(point(wall(k), pout[k]));
            land_out[k] = v2.len() - 1;
        } else {
            land_out[k] = land_in[k];
        }
    }
    let mut checked: Vec<[usize; 4]> = Vec::new();
    for k in 0..n {
        checked.push([base[k], base[(k + 1) % n], land_in[(k + 1) % n], land_out[k]]);
    }
    // Corner blocks: a junction's two landings and the kink between them.
    let mut kink_id: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut junction_kink: Vec<Option<usize>> = vec![None; n];
    for k in 0..n {
        if !junction[k] {
            continue;
        }
        let bi = match (wall(k + n - 1), wall(k)) {
            (Wall::Boss(b), _) | (_, Wall::Boss(b)) => b,
            _ => return Err("a corner between two walls neither of which is a boss".into()),
        };
        let v = at(k);
        let (ki, _) = kinks
            .get(bi)
            .into_iter()
            .flatten()
            .enumerate()
            .map(|(i, p)| (i, norm(sub(*p, v))))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .ok_or("the outline does not leave a boss the arm runs into")?;
        let id = *kink_id.entry((bi, ki)).or_insert_with(|| {
            v2.push(kinks.get(bi).and_then(|kk| kk.get(ki)).copied().unwrap_or_default());
            v2.len() - 1
        });
        junction_kink[k] = Some(id);
        checked.push([base[k], land_in[k], id, land_out[k]]);
    }
    quads.extend(checked.iter().copied());

    // The annuli: per boss, the landings on its circle in order, the kinks, and its own side.
    let mut rings_made = 0;
    for (bi, b) in bosses.iter().enumerate() {
        let edges: Vec<usize> = (0..n).filter(|&k| wall(k) == Wall::Boss(bi)).collect();
        let starts: Vec<usize> = edges.iter().copied().filter(|&k| wall(k + n - 1) != Wall::Boss(bi)).collect();
        let [ja] = starts.as_slice() else {
            return Err(format!("the arm meets the boss of radius {} in {} places", length(b.r), starts.len()));
        };
        let mut run = vec![*ja];
        while run.len() < n && wall(run.last().copied().unwrap_or(0) + 1) == Wall::Boss(bi) {
            run.push((run.last().copied().unwrap_or(0) + 1) % n);
        }
        let jb = (run.last().copied().unwrap_or(0) + 1) % n;
        let (Some(ka), Some(kb)) = (junction_kink[*ja], junction_kink[jb]) else {
            return Err(format!("the arm's run along the boss of radius {} has no kinks at its ends", length(b.r)));
        };
        // Counterclockwise round the boss: the run backwards, then the kink at its start, the
        // boss's own side and the kink at its end.
        let mut outer: Vec<usize> = std::iter::once(land_out[*ja]).chain(run.iter().skip(1).map(|&k| land_out[k])).chain([land_in[jb]]).collect();
        outer.reverse();
        outer.push(ka);
        let ang = |p: V2| (p[1] - b.c[1]).atan2(p[0] - b.c[0]);
        let mut angles: Vec<f64> = Vec::with_capacity(outer.len());
        for &id in &outer {
            let mut x = ang(v2.get(id).copied().unwrap_or_default());
            if let Some(&last) = angles.last() {
                while x <= last {
                    x += std::f64::consts::TAU;
                }
            }
            angles.push(x);
        }
        let first = angles.first().copied().unwrap_or(0.0);
        let last = angles.last().copied().unwrap_or(0.0);
        let mut end = ang(v2.get(kb).copied().unwrap_or_default());
        while end <= last {
            end += std::f64::consts::TAU;
        }
        if end - first >= std::f64::consts::TAU - 1e-9 || angles.len() < 3 {
            return Err(format!("the arm's landings round the boss of radius {} are out of order", length(b.r)));
        }
        // The own side in steps of the run's mean step (within bounds).
        let step = ((angles[angles.len() - 2] - first) / (angles.len() - 2) as f64).clamp(MIN_OWN_ARC, MAX_OWN_ARC);
        let span = end - last;
        let count = ((span / step).round() as usize).max(1);
        for t in 1..count {
            let x = last + span * t as f64 / count as f64;
            v2.push(add(b.c, [b.r * x.cos(), b.r * x.sin()]));
            outer.push(v2.len() - 1);
            angles.push(x);
        }
        outer.push(kb);
        angles.push(end);
        // The inner loop: the hole's point on each outer point's ray, the corners on the nearest.
        let mut inner: Vec<V2> = Vec::with_capacity(angles.len());
        for &x in &angles {
            let h = b.inner.ray(b.c, [x.cos(), x.sin()], 0.0).ok_or("a ray from a boss's centre misses its hole")?;
            inner.push(h.1);
        }
        let wrap = |x: f64| (x + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI;
        let mut snapped = vec![false; inner.len()];
        for c in &b.corners {
            let ca = ang(*c);
            let best = (0..angles.len()).min_by(|p, q| wrap(angles[*p] - ca).abs().total_cmp(&wrap(angles[*q] - ca).abs())).ok_or("an empty ring")?;
            if snapped[best] {
                return Err(format!("the ring round the hole in the boss of radius {} has too few blocks for its corners", length(b.r)));
            }
            snapped[best] = true;
            inner[best] = *c;
        }
        let first_inner = v2.len();
        v2.extend(&inner);
        let m = outer.len();
        for t in 0..m {
            let q = [first_inner + t, first_inner + (t + 1) % m, outer[(t + 1) % m], outer[t]];
            checked.push(q);
            quads.push(q);
        }
        rings_made += 1;
    }
    for q in &checked {
        let pts = q.map(|i| v2.get(i).copied().unwrap_or_default());
        if !convex(&pts) {
            return Err(format!(
                "a block of the layout folds near ({:.4}, {:.4})",
                0.25 * pts.iter().map(|p| p[0]).sum::<f64>(),
                0.25 * pts.iter().map(|p| p[1]).sum::<f64>()
            ));
        }
    }

    // The slabs: each keeps the blocks whose middle lies inside the body.
    let mut t = Topology { vertices: Vec::new(), blocks: Vec::new() };
    let mut v3: BTreeMap<(usize, usize), u32> = BTreeMap::new();
    let quad_area = |q: &[usize; 4]| {
        let p = q.map(|i| v2.get(i).copied().unwrap_or_default());
        0.5 * (cross(p[0], p[1]) + cross(p[1], p[2]) + cross(p[2], p[3]) + cross(p[3], p[0])).abs()
    };
    for (s, w) in levels.windows(2).enumerate() {
        let (z0, z1) = (w[0], w[1]);
        let zm = 0.5 * (z0 + z1);
        let mut covered = 0.0;
        for q in &quads {
            let p = q.map(|i| v2.get(i).copied().unwrap_or_default());
            let c = scale(p.iter().fold([0.0; 2], |a, x| add(a, *x)), 0.25);
            if winding(inp.mesh, world([c[0], c[1], zm])) < 0.5 {
                continue;
            }
            covered += quad_area(q);
            let mut b = [0u32; 8];
            for (h, (lvl, z)) in [(s, z0), (s + 1, z1)].into_iter().enumerate() {
                for (c, &id) in q.iter().enumerate() {
                    let key = (id, lvl);
                    let vid = match v3.get(&key) {
                        Some(&x) => x,
                        None => {
                            let pt = v2.get(id).copied().unwrap_or_default();
                            t.vertices.push(world([pt[0], pt[1], z]));
                            let x = (t.vertices.len() - 1) as u32;
                            v3.insert(key, x);
                            x
                        }
                    };
                    b[4 * h + c] = vid;
                }
            }
            let pos = b.map(|v| t.vertices.get(v as usize).copied().unwrap_or_default());
            if corner_det(&pos) < 0.0 {
                b = [b[0], b[3], b[2], b[1], b[4], b[7], b[6], b[5]];
            }
            t.blocks.push(b);
            if t.blocks.len() > MAX_BLOCKS {
                return Err(format!("more than {MAX_BLOCKS} blocks: too many for one layout"));
            }
        }
        // The blocks cover the slab's section: no wall of it is left out of the layout.
        let section = sections.get(s).map_or(0.0, |l| area(l));
        if (covered - section).abs() > AREA_TOLERANCE * section.abs() {
            return Err(format!(
                "the layout covers {} of the {} section between {} and {}",
                area_shown(covered),
                area_shown(section),
                length(z0),
                length(z1)
            ));
        }
    }
    t.check().map_err(|e| format!("the layout's blocks do not fit: {e}"))?;
    Ok(Layout {
        topology: t,
        rings: rings_made,
        layers: Vec::new(),
        planes: [gp[0].len(), gp[1].len(), levels.len()],
        frames: Vec::new(),
        pillow: Some(margin),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_planes_keep_clear_of_a_filleted_corner() {
        let p = planes(0.0, 20.0, &[(10.0, 4.5 / FILLET_CLEARANCE)], 2.0, 1e-9);
        assert!(p.windows(2).all(|w| w[1] > w[0]), "{p:?}");
        let i = p.iter().position(|x| (x - 10.0).abs() < 1e-12).unwrap();
        assert!(p[i] - p[i - 1] >= 4.5 && p[i + 1] - p[i] >= 4.5, "{p:?}");
        assert!(planes(0.0, 20.0, &[], 2.0, 1e-9).windows(2).all(|w| (w[1] - w[0] - 2.0).abs() < 1e-12));
    }

    #[test]
    fn a_box_slices_into_its_rectangle() {
        let shape = fieldcraft_geometry::block([0.0; 3], [4.0, 3.0, 2.0]).unwrap();
        let mesh = shape.tessellate(0.01).unwrap();
        let surfaces = shape.face_surfaces().unwrap();
        let inp = Input {
            shape: &shape,
            surfaces: &surfaces,
            mesh: &mesh,
            direction: [0.0, 0.0, 1.0],
            rings: &[],
            blends: Vec::new(),
            shell: None,
            size: 1.0,
            show: &|q, v| crate::units::shown(crate::units::geometry(), crate::units::geometry(), q, v),
        };
        let loops = slice(&inp, &|p| p, 1.0, 1e-6).unwrap();
        assert_eq!(loops.len(), 1, "{loops:?}");
        assert!((loops[0].area() - 12.0).abs() < 1e-9, "{}", loops[0].area());
        assert!(loops[0].contains([2.0, 1.5]) && !loops[0].contains([5.0, 1.5]));
    }
}
