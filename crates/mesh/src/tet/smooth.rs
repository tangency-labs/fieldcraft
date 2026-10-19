//! Vertex smoothing: optimal Delaunay triangulation (ODT) moves, and optimization-based smoothing
//! of the worst tetrahedron around a vertex.
//!
//! Boundary vertices move on their B-rep face (in its tangent plane, then projected) or along
//! their feature curve; corners never move. Every move keeps every tetrahedron around the vertex
//! positive (exact `orient3d`), or is not made.
//!
//! Sources:
//! - L. Chen, J. Xu, *Optimal Delaunay triangulations*, J. Comput. Math. 22(2), 2004; L. Chen,
//!   M. Holst, *Efficient mesh optimization schemes based on Optimal Delaunay Triangulations*,
//!   CMAME 200, 2011; P. Alliez, D. Cohen-Steiner, M. Yvinec, M. Desbrun, *Variational
//!   tetrahedral meshing*, ACM TOG 24(3), 2005: the ODT update moves a vertex to the
//!   volume-weighted average of the circumcentres of its tetrahedra, weighted here by the density
//!   `h⁻³` of the size field (J. Tournois, C. Wormser, P. Alliez, M. Desbrun, *Interleaving
//!   Delaunay refinement and optimization for practical isotropic tetrahedron mesh generation*,
//!   ACM TOG 28(3), 2009, which interleaves these moves with refinement).
//! - L. Freitag, M. Jones, P. Plassmann, *An efficient parallel algorithm for mesh smoothing*,
//!   IMR 1995; L. Freitag, C. Ollivier-Gooch, *Tetrahedral mesh improvement using swapping and
//!   smoothing*, IJNME 40, 1997; Klingner and Shewchuk 2007, §3.1: nonsmooth optimization of the
//!   worst quality around a vertex (active set, steepest ascent along the smallest vector in the
//!   convex hull of the active gradients, line search), applied only when the worst improves.

use crate::Result;

use super::complex::{NONE, TetComplex, VertexKind, face_toward};
use super::measure::{self, P, add, dot, scale, sub};
use super::refine::Context;

/// The difference step, relative to the edge's length.
const STEP: f64 = 1e-6; // unit-factors: relative to the edge's length

/// The worst quality around `v` with `v` at `p` (the tetrahedra of `star`).
fn star_quality(cx: &TetComplex, v: u32, p: P, star: &[u32]) -> Result<f64> {
    let mut q = f64::INFINITY;
    for &t in star {
        let vs = cx.tet(t)?;
        let mut x = cx.coords(vs)?;
        for (i, &u) in vs.iter().enumerate() {
            if u == v
                && let Some(s) = x.get_mut(i)
            {
                *s = p;
            }
        }
        q = q.min(measure::quality(&x));
    }
    Ok(q)
}

/// How `v` may move: `None` for a corner, the identity for an interior vertex, the tangent
/// plane's normal for a face vertex, the tangent line for a curve vertex.
#[derive(Debug, Clone, Copy)]
enum Freedom {
    Free,
    Plane(P),
    Line(P),
}

fn freedom(cx: &mut TetComplex, v: u32, star: &[u32]) -> Result<Option<Freedom>> {
    Ok(match cx.kind(v) {
        VertexKind::Corner => None,
        VertexKind::Interior => Some(Freedom::Free),
        VertexKind::Face(_) => {
            let mut n = [0.0; 3];
            for &t in star {
                let vs = cx.tet(t)?;
                let ns = cx.adj(t)?;
                for k in 0..4 {
                    if ns[k] == NONE && vs[k] != v {
                        let f = face_toward(vs, k);
                        let (a, b, c) = (cx.xyz(f[0])?, cx.xyz(f[1])?, cx.xyz(f[2])?);
                        n = add(n, measure::cross(sub(b, a), sub(c, a)));
                    }
                }
            }
            let l = measure::norm(n);
            if !(l > 0.0) {
                return Ok(None);
            }
            Some(Freedom::Plane(scale(n, 1.0 / l)))
        }
        VertexKind::Curve(c) => {
            let ends: Vec<u32> = cx.neighbours(v)?.into_iter().filter(|&u| cx.curve(v, u) == Some(c)).collect();
            let [a, b] = ends.as_slice() else { return Ok(None) };
            let d = sub(cx.xyz(*b)?, cx.xyz(*a)?);
            let l = measure::norm(d);
            if !(l > 0.0) {
                return Ok(None);
            }
            Some(Freedom::Line(scale(d, 1.0 / l)))
        }
    })
}

fn restrict(d: P, f: Freedom) -> P {
    match f {
        Freedom::Free => d,
        Freedom::Plane(n) => sub(d, scale(n, dot(d, n))),
        Freedom::Line(t) => scale(t, dot(d, t)),
    }
}

/// Try moving `v` to `target` (restricted to its freedom and placed on its boundary entity),
/// at full, half and quarter steps; accept the first that keeps the star valid and its worst
/// quality at least `floor`. Returns whether it moved.
fn try_move(cx: &mut TetComplex, ctx: &Context, v: u32, star: &[u32], f: Freedom, target: P, floor: f64) -> Result<bool> {
    let x = cx.xyz(v)?;
    let d = restrict(sub(target, x), f);
    for s in [1.0, 0.5, 0.25] {
        let Some(p) = ctx.place(cx.kind(v), add(x, scale(d, s))) else { continue };
        if p.coords() == x || !cx.move_is_valid(v, &p, star)? {
            continue;
        }
        if star_quality(cx, v, p.coords(), star)? >= floor {
            cx.move_point(v, p)?;
            return Ok(true);
        }
    }
    Ok(false)
}

/// The quality a smoothing move may not go below (sin 17.5°) unless the star was already worse.
const ODT_FLOOR: f64 = 0.3;

/// One ODT pass over the vertices (ascending ids) that may move. Returns the number moved.
pub fn odt_pass(cx: &mut TetComplex, ctx: &Context) -> Result<usize> {
    let mut moved = 0;
    for v in 0..cx.point_count() as u32 {
        if !cx.is_used(v) || cx.kind(v) == VertexKind::Corner {
            continue;
        }
        let star = cx.star(v)?;
        let Some(f) = freedom(cx, v, &star)? else { continue };
        let (mut num, mut den) = ([0.0; 3], 0.0);
        for &t in &star {
            let x = cx.coords(cx.tet(t)?)?;
            let vol = measure::volume6(&x) / 6.0;
            let g = measure::centroid(&x);
            let (_, lmax) = measure::edge_range(&x);
            let c = measure::circumcenter(&x).filter(|c| measure::dist(*c, g) <= lmax).unwrap_or(g);
            let h = ctx.h(g, lmax);
            let w = vol / (h * h * h);
            if !(w > 0.0 && w.is_finite()) {
                continue;
            }
            num = add(num, scale(c, w));
            den += w;
        }
        if !(den > 0.0) {
            continue;
        }
        let target = scale(num, 1.0 / den);
        let before = star_quality(cx, v, cx.xyz(v)?, &star)?;
        if try_move(cx, ctx, v, &star, f, target, before.min(ODT_FLOOR))? {
            moved += 1;
        }
    }
    Ok(moved)
}

/// The smallest vector in the convex hull of `g` (Frank–Wolfe iterations).
fn min_norm(g: &[P]) -> P {
    let Some(&first) = g.first() else { return [0.0; 3] };
    let mut d = first;
    for _ in 0..64 {
        let mut s = first;
        let mut best = f64::INFINITY;
        for &gi in g {
            let v = dot(gi, d);
            if v < best {
                best = v;
                s = gi;
            }
        }
        let e = sub(d, s);
        let gap = dot(e, d);
        let ee = dot(e, e);
        if !(gap > 1e-30) || !(ee > 0.0) {
            break;
        }
        let gamma = (gap / ee).clamp(0.0, 1.0);
        d = sub(d, scale(e, gamma));
    }
    d
}

/// The star of a vertex as plain coordinates: per tetrahedron, its four corners and the slot of
/// the moving vertex (so a trial position costs no lookups).
struct LocalStar {
    tets: Vec<([P; 4], usize)>,
}

impl LocalStar {
    fn new(cx: &TetComplex, v: u32, star: &[u32]) -> Result<LocalStar> {
        let mut tets = Vec::with_capacity(star.len());
        for &t in star {
            let vs = cx.tet(t)?;
            let i = vs.iter().position(|&u| u == v).unwrap_or(0);
            tets.push((cx.coords(vs)?, i));
        }
        Ok(LocalStar { tets })
    }

    fn quality_at(&self, k: usize, p: P) -> f64 {
        let Some((x, i)) = self.tets.get(k) else { return f64::NEG_INFINITY };
        let mut y = *x;
        if let Some(s) = y.get_mut(*i) {
            *s = p;
        }
        measure::quality(&y)
    }

    fn worst_at(&self, p: P) -> f64 {
        (0..self.tets.len()).map(|k| self.quality_at(k, p)).fold(f64::INFINITY, f64::min)
    }
}

/// Optimization-based smoothing of `v`: raise the worst quality of its star (see the module
/// docs). Returns whether it moved.
pub fn smooth_vertex(cx: &mut TetComplex, ctx: &Context, v: u32) -> Result<bool> {
    if !cx.is_used(v) || cx.kind(v) == VertexKind::Corner {
        return Ok(false);
    }
    let star = cx.star(v)?;
    let Some(f) = freedom(cx, v, &star)? else { return Ok(false) };
    let local = LocalStar::new(cx, v, &star)?;
    let x0 = cx.xyz(v)?;
    // A length scale: the mean edge length from v.
    let mut len = 0.0;
    let mut count = 0usize;
    for (x, i) in &local.tets {
        for (k, p) in x.iter().enumerate() {
            if k != *i {
                len += measure::dist(*p, x0);
                count += 1;
            }
        }
    }
    let len = len / count.max(1) as f64;
    if !(len > 0.0 && len.is_finite()) {
        return Ok(false);
    }
    let delta = STEP * len;
    let mut x = x0;
    let mut worst = local.worst_at(x);
    let mut moved = false;
    for _ in 0..8 {
        // Qualities and gradients (central differences) of the tetrahedra near the worst.
        let mut qs = Vec::with_capacity(local.tets.len());
        for k in 0..local.tets.len() {
            let q = local.quality_at(k, x);
            let mut g = [0.0; 3];
            if q <= worst + 0.1 {
                for (c, gc) in g.iter_mut().enumerate() {
                    let (mut a, mut b) = (x, x);
                    a[c] += delta;
                    b[c] -= delta;
                    *gc = (local.quality_at(k, a) - local.quality_at(k, b)) / (2.0 * delta);
                }
            }
            qs.push((q, restrict(g, f)));
        }
        let active: Vec<P> = qs.iter().filter(|s| s.0 <= worst + 1e-4).map(|s| s.1).collect();
        let d = min_norm(&active);
        let dd = dot(d, d);
        if !(dd > 1e-20) {
            break;
        }
        // The step at which the linear model of an inactive tetrahedron meets the worst, capped
        // to a fifth of the edge length.
        let mut alpha = 0.2 * len / dd.sqrt();
        for &(q, g) in &qs {
            let gd = dot(g, d);
            if q > worst + 1e-4 && q <= worst + 0.1 && gd < dd {
                alpha = alpha.min((q - worst) / (dd - gd));
            }
        }
        let mut accepted = false;
        for _ in 0..6 {
            let Some(p) = ctx.place(cx.kind(v), add(x, scale(d, alpha))) else { break };
            let q = local.worst_at(p.coords());
            if p.coords() != x && q > worst && cx.move_is_valid(v, &p, &star)? {
                cx.move_point(v, p)?;
                x = p.coords();
                worst = q;
                accepted = true;
                moved = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    Ok(moved)
}
