//! Mortar segmentation and integration: the overlap of a secondary (non-mortar) face with a main
//! (mortar) face, and the integration points on it where both faces' shape functions and
//! positions are known.
//!
//! **Segmentation** (M. A. Puso, T. A. Laursen, *A mortar segment-to-segment contact method for
//! large deformation solid mechanics*, CMAME 193:601–629, 2004, §3; A. Popp, M. W. Gee, W. A.
//! Wall, *A finite deformation mortar contact formulation using a primal–dual active set
//! strategy*, IJNME 79:1354–1391, 2009, §3.3). Both faces are projected along the secondary
//! face's unit normal n onto its plane (through its centroid), and the main polygon is clipped by
//! the secondary one (Sutherland–Hodgman; a quadrilateral whose projection is not convex is
//! split into its two triangles first). Each overlap polygon is split into a fan of triangles
//! integrated with the seven-point degree-5 rule (J. Radon 1948; A. H. Stroud, *Approximate
//! Calculation of Multiple Integrals*, Prentice-Hall 1971, T2:5-1), which integrates the products
//! of two shape functions that the tie and the patch test need exactly on triangles and
//! parallelograms (on a general quadrilateral the inverse of its map is not a polynomial, and the
//! rule is accurate to its degree).
//!
//! **Points.** At an integration point y of the plane, each face's coordinates are found by
//! inverting its projected map: the linear one of a triangle, or the bilinear one of a
//! quadrilateral, x(u, v) = a + b u + c v + d u v, whose inverse is the root in [0, 1] of the
//! quadratic (b × d) u² + (b × c − h × d) u − h × c = 0 (h = y − a, × the scalar cross product of
//! the plane). The projection of a bilinear patch onto a plane is the bilinear patch of its
//! projected corners, so the main face's point found this way is the point of the main surface
//! along n through y, warped or not. The [`Sample`] then holds both faces' shape functions, both
//! points in space, the normal and the weight (the plane's area element times the rule's weight).
//!
//! The kernel is written over [`Real`], so the same code gives the values (`f64`), their
//! gradient ([`crate::dual::Dual`]) or their gradient and Hessian ([`crate::dual::Hyper`]) with
//! respect to the corner coordinates, variables `3 i + k` for secondary corner i and
//! `12 + 3 j + k` for main corner j: the exact linearisation of the segmentation, the normal and
//! the projection, which the consistent tangent of mortar contact needs (Popp et al. 2009, §4).
//!
//! **Edges of 2D models** ([`integrate_edges`]): plane strain, plane stress and axisymmetric
//! models touch along element edges, two-node segments in the XY plane. The main edge is
//! projected along the secondary edge's unit normal onto its line, which maps each main point to
//! one parameter of the secondary edge linearly, so the overlap is an interval of that parameter
//! (B. Yang, T. A. Laursen, X. Meng, *Two dimensional mortar contact methods for large
//! deformation frictional sliding*, IJNME 62:1183–1225, 2005, §3). It is integrated with the
//! three-point Gauss–Legendre rule, exact for the degree-3 products Φ_a N_b r that the
//! axisymmetric patch test needs (the edge's out-of-plane measure, the thickness or 2π r at the
//! secondary point, is part of the weight). Points stay three-component with z = 0, so the solver
//! takes edges through the same variables as faces.

use crate::dual::Real;
use fieldcraft_model::Extent;

/// A point in space over a scalar.
pub type Point<T> = [T; 3];

/// A point of the secondary face's plane.
type Flat<T> = [T; 2];

/// Radon's seven-point rule on a triangle, exact for degree 5: barycentric coordinates of the
/// second and third corners and weights summing to one.
fn rule() -> [([f64; 2], f64); 7] {
    let s = 15f64.sqrt();
    let (a, b) = ((6.0 - s) / 21.0, (9.0 + 2.0 * s) / 21.0);
    let (c, d) = ((6.0 + s) / 21.0, (9.0 - 2.0 * s) / 21.0);
    let (wa, wc) = ((155.0 - s) / 1200.0, (155.0 + s) / 1200.0);
    [([1.0 / 3.0, 1.0 / 3.0], 0.225), ([a, a], wa), ([b, a], wa), ([a, b], wa), ([c, c], wc), ([d, c], wc), ([c, d], wc)]
}

/// One integration point of a segment.
#[derive(Debug, Clone, Copy)]
pub struct Sample<T> {
    /// The rule's weight times the plane's area element (the secondary face's projected area).
    pub weight: T,
    /// The secondary face's shape functions at the point (0 past its corners).
    pub secondary: [T; 4],
    /// The main face's shape functions at the point.
    pub main: [T; 4],
    /// The point on the secondary face, and on the main face along the normal.
    pub xs: Point<T>,
    pub xm: Point<T>,
    /// The secondary face's unit normal (outward from its element).
    pub normal: Point<T>,
}

impl<T: Real> Sample<T> {
    /// The gap along the secondary normal, (x_m − x_s) · n: positive where the faces are apart.
    pub fn gap(&self) -> T {
        dot(sub(self.xm, self.xs), self.normal)
    }
}

fn sub<T: Real>(a: Point<T>, b: Point<T>) -> Point<T> {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot<T: Real>(a: Point<T>, b: Point<T>) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross<T: Real>(a: Point<T>, b: Point<T>) -> Point<T> {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn scale<T: Real>(a: Point<T>, s: T) -> Point<T> {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn cross2<T: Real>(a: Flat<T>, b: Flat<T>) -> T {
    a[0] * b[1] - a[1] * b[0]
}

fn sub2<T: Real>(a: Flat<T>, b: Flat<T>) -> Flat<T> {
    [a[0] - b[0], a[1] - b[1]]
}

/// A face's area vector: its outward normal times its area (flat triangle; for a bilinear
/// quadrilateral half the cross product of its diagonals, the area of its projection on the
/// plane normal to it).
pub fn area_vector<T: Real>(x: &[Point<T>]) -> Option<Point<T>> {
    let half = |v: Point<T>| scale(v, T::cst(0.5));
    match x {
        [a, b, c] => Some(half(cross(sub(*b, *a), sub(*c, *a)))),
        [a, b, c, d] => Some(half(cross(sub(*c, *a), sub(*d, *b)))),
        _ => None,
    }
}

/// A projected face as convex counter-clockwise pieces for the clipping: itself when it is
/// convex (no artificial diagonal, whose crossings with the other face's edges would change the
/// overlap's polygons, and the integration's round-off, as the faces move), else its two
/// triangles.
fn pieces<T: Real>(points: &[Flat<T>]) -> Vec<Vec<Flat<T>>> {
    let mut whole = points.to_vec();
    if signed_area(&whole).val() < 0.0 {
        whole.reverse();
    }
    let n = whole.len();
    let convex = (0..n).all(|i| {
        let (a, b, c) = (whole[i], whole[(i + 1) % n], whole[(i + 2) % n]);
        cross2(sub2(b, a), sub2(c, b)).val() > 0.0
    });
    if convex || n == 3 {
        return vec![whole];
    }
    [[0, 1, 2], [0, 2, 3]].iter().map(|t| counter_clockwise(t, points)).collect()
}

fn signed_area<T: Real>(p: &[Flat<T>]) -> T {
    let mut s = T::cst(0.0);
    for (i, a) in p.iter().enumerate() {
        let b = p.get((i + 1) % p.len()).copied().unwrap_or(*a);
        s = s + cross2(*a, b);
    }
    s * 0.5
}

/// The corners `t` of `points`, counter-clockwise.
fn counter_clockwise<T: Real>(t: &[usize; 3], points: &[Flat<T>]) -> Vec<Flat<T>> {
    let mut out: Vec<Flat<T>> = t.iter().filter_map(|i| points.get(*i).copied()).collect();
    if signed_area(&out).val() < 0.0 {
        out.reverse();
    }
    out
}

/// A vertex this close to a clip edge, relative to the edge's length, is on it.
const ON_EDGE: f64 = 1e-8;

/// A quadratic term below this fraction of the linear one is round-off: the equation is linear.
/// A span below it, along an edge of unit length, is no overlap.
const NEGLIGIBLE: f64 = 1e-12; // unit-factors: relative to the linear term, or to the edge's length

/// `subject` clipped by the convex counter-clockwise polygon `clip` (Sutherland–Hodgman), the
/// sides decided on the values. A vertex within [`ON_EDGE`] of a clip edge's line counts as on
/// it (kept, and never cut): the intersection of two nearly coincident edges is ill-conditioned,
/// its position and above all its derivatives swamped by round-off, while the vertex itself is
/// as good a corner of the overlap to that precision. Edges are cut only between a vertex
/// clearly inside and one clearly outside, where the intersection is well conditioned.
fn clip<T: Real>(subject: &[Flat<T>], clip: &[Flat<T>]) -> Vec<Flat<T>> {
    let mut out: Vec<Flat<T>> = subject.to_vec();
    for (i, a) in clip.iter().enumerate() {
        let Some(b) = clip.get((i + 1) % clip.len()).copied() else { break };
        let input = std::mem::take(&mut out);
        let edge = sub2(b, *a);
        let eps = ON_EDGE * (edge[0].val() * edge[0].val() + edge[1].val() * edge[1].val());
        let side = |p: Flat<T>| cross2(edge, sub2(p, *a));
        for (j, p) in input.iter().enumerate() {
            let Some(q) = input.get((j + 1) % input.len()).copied() else { break };
            let (sp, sq) = (side(*p), side(q));
            if sp.val() >= -eps {
                out.push(*p);
            }
            if (sp.val() > eps && sq.val() < -eps) || (sp.val() < -eps && sq.val() > eps) {
                let t = sp / (sp - sq);
                out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
            }
        }
        if out.is_empty() {
            break;
        }
    }
    out
}

/// The shape functions of a face (of 3 or 4 projected corners `p`) at plane point `y`: the
/// inverse of its linear or bilinear map. `None` when it has no inverse there.
fn shape_at<T: Real>(p: &[Flat<T>], y: Flat<T>) -> Option<[T; 4]> {
    let zero = T::cst(0.0);
    match p {
        [p0, p1, p2] => {
            let (e1, e2, h) = (sub2(*p1, *p0), sub2(*p2, *p0), sub2(y, *p0));
            let det = cross2(e1, e2);
            if !(det.val().abs() > 0.0) {
                return None;
            }
            let v = cross2(h, e2) / det;
            let w = cross2(e1, h) / det;
            Some([T::cst(1.0) - v - w, v, w, zero])
        }
        [p0, p1, p2, p3] => {
            let b = sub2(*p1, *p0);
            let c = sub2(*p3, *p0);
            let d = [p0[0] - p1[0] + p2[0] - p3[0], p0[1] - p1[1] + p2[1] - p3[1]];
            let h = sub2(y, *p0);
            let k2 = cross2(b, d);
            let k1 = cross2(b, c) - cross2(h, d);
            let k0 = -cross2(h, c);
            let u = if k2.val().abs() <= NEGLIGIBLE * k1.val().abs() {
                if !(k1.val().abs() > 0.0) {
                    return None;
                }
                -k0 / k1
            } else {
                // The roots of k2 u² + k1 u + k0, computed without cancellation; the one in [0, 1]
                // (the nearer when round-off puts both outside).
                let disc = k1 * k1 - k2 * k0 * 4.0;
                let root = if disc.val() > 0.0 { disc.sqrt() } else { T::cst(0.0) };
                let q = if k1.val() >= 0.0 { -(k1 + root) * 0.5 } else { (root - k1) * 0.5 };
                let mut candidates = [None, None];
                if q.val().abs() > 0.0 {
                    candidates[0] = Some(q / k2);
                    candidates[1] = Some(k0 / q);
                } else {
                    candidates[0] = Some(-k1 / (k2 * 2.0));
                }
                let distance = |u: f64| {
                    if u < 0.0 {
                        -u
                    } else if u > 1.0 {
                        u - 1.0
                    } else {
                        0.0
                    }
                };
                candidates.into_iter().flatten().filter(|u| u.val().is_finite()).min_by(|a, b| distance(a.val()).total_cmp(&distance(b.val())))?
            };
            // v from the better conditioned component of h − u b = v (c + d u).
            let m = [c[0] + d[0] * u, c[1] + d[1] * u];
            let k = if m[0].val().abs() >= m[1].val().abs() { 0 } else { 1 };
            if !(m[k].val().abs() > 0.0) {
                return None;
            }
            let v = (h[k] - b[k] * u) / m[k];
            let (s, t) = (u * 2.0 - 1.0, v * 2.0 - 1.0);
            let one = T::cst(1.0);
            Some([(one - s) * (one - t) * 0.25, (one + s) * (one - t) * 0.25, (one + s) * (one + t) * 0.25, (one - s) * (one + t) * 0.25])
        }
        _ => None,
    }
}

fn interpolate<T: Real>(x: &[Point<T>], n: &[T; 4]) -> Point<T> {
    let mut out = [T::cst(0.0); 3];
    for (xi, ni) in x.iter().zip(n) {
        for k in 0..3 {
            out[k] = out[k] + xi[k] * *ni;
        }
    }
    out
}

/// Integrate over the overlap of secondary face `secondary` with main face `main` (3 or 4
/// corners each, in their faces' node order): `f` is called at every integration point. `None`
/// when a face has another number of corners or the secondary face has no normal (collapsed);
/// faces that do not overlap give no point.
pub fn integrate<T: Real>(secondary: &[Point<T>], main: &[Point<T>], mut f: impl FnMut(&Sample<T>)) -> Option<()> {
    if !matches!(main.len(), 3 | 4) {
        return None;
    }
    let area = area_vector(secondary)?;
    let length = dot(area, area).sqrt();
    if !(length.val() > 0.0 && length.val().is_finite()) {
        return None;
    }
    let normal = scale(area, T::cst(1.0) / length);
    let count = T::cst(secondary.len() as f64);
    let mut center = [T::cst(0.0); 3];
    for x in secondary {
        for k in 0..3 {
            center[k] = center[k] + x[k] / count;
        }
    }
    // A basis of the plane, (e1, e2, n) right-handed, e1 towards the first corner.
    let r = sub(*secondary.first()?, center);
    let r = sub(r, scale(normal, dot(r, normal)));
    let rl = dot(r, r).sqrt();
    if !(rl.val() > 0.0) {
        return None;
    }
    let e1 = scale(r, T::cst(1.0) / rl);
    let e2 = cross(normal, e1);
    let flat = |x: &Point<T>| -> Flat<T> {
        let d = sub(*x, center);
        [dot(d, e1), dot(d, e2)]
    };
    let own: Vec<Flat<T>> = secondary.iter().map(flat).collect();
    let theirs: Vec<Flat<T>> = main.iter().map(flat).collect();
    let rule = rule();
    let (own_pieces, their_pieces) = (pieces(&own), pieces(&theirs));
    for s in &own_pieces {
        for m in &their_pieces {
            let overlap = clip(m, s);
            let Some(o) = overlap.first().copied() else { continue };
            for w in overlap.windows(2).skip(1) {
                let (p, q) = (w[0], w[1]);
                let doubled = cross2(sub2(p, o), sub2(q, o));
                if !(doubled.val() > 0.0) {
                    continue;
                }
                for (bary, weight) in rule {
                    let y = [0, 1].map(|k| o[k] + (p[k] - o[k]) * bary[0] + (q[k] - o[k]) * bary[1]);
                    let (Some(ns), Some(nm)) = (shape_at(&own, y), shape_at(&theirs, y)) else { continue };
                    f(&Sample {
                        weight: doubled * (0.5 * weight),
                        secondary: ns,
                        main: nm,
                        xs: interpolate(secondary, &ns),
                        xm: interpolate(main, &nm),
                        normal,
                    });
                }
            }
        }
    }
    Some(())
}

/// The out-of-plane measure of a 2D model at a point: the thickness, or 2π x.
fn measure<T: Real>(extent: Extent, x: T) -> T {
    match extent {
        Extent::Thickness(t) => T::cst(t),
        Extent::Revolved => x * (2.0 * std::f64::consts::PI),
    }
}

/// The unit tangent of a two-point edge in the XY plane and its length.
fn edge_frame<T: Real>(x: &[Point<T>]) -> Option<(Flat<T>, T)> {
    let [a, b] = x else { return None };
    let t = [b[0] - a[0], b[1] - a[1]];
    let length = (t[0] * t[0] + t[1] * t[1]).sqrt();
    if !(length.val() > 0.0 && length.val().is_finite()) {
        return None;
    }
    let inv = T::cst(1.0) / length;
    Some(([t[0] * inv, t[1] * inv], length))
}

/// An edge's area vector in a 2D model (two corners in the XY plane): its unit normal
/// (t_y, −t_x), outward when the edge runs counter-clockwise round its element, times its
/// length and its out-of-plane measure at its midpoint (∫ 2π r ds exactly for a straight edge).
pub fn edge_area_vector<T: Real>(x: &[Point<T>], extent: Extent) -> Option<Point<T>> {
    let (t, length) = edge_frame(x)?;
    let [a, b] = x else { return None };
    let size = length * measure(extent, (a[0] + b[0]) * 0.5);
    Some([t[1] * size, -(t[0] * size), T::cst(0.0)])
}

/// Integrate over the overlap of secondary edge `secondary` with main edge `main` of a 2D model
/// (two corners each, z ignored), as [`integrate`] does for faces: the samples' shape functions
/// are those of the edges' two corners (the last two are zero), the normal is the secondary
/// edge's (t_y, −t_x, 0) and the weight carries the out-of-plane measure `extent` at the
/// secondary point. `None` when an edge does not have two corners or the secondary one has no
/// length; edges that do not overlap, or a main edge along the normal, give no point.
pub fn integrate_edges<T: Real>(secondary: &[Point<T>], main: &[Point<T>], extent: Extent, mut f: impl FnMut(&Sample<T>)) -> Option<()> {
    let (t, length) = edge_frame(secondary)?;
    let ([s0, s1], [m0, m1]) = (secondary, main) else { return None };
    let inv = T::cst(1.0) / length;
    let along = |p: &Point<T>| ((p[0] - s0[0]) * t[0] + (p[1] - s0[1]) * t[1]) * inv;
    let (a, b) = (along(m0), along(m1));
    let span = b - a;
    if !(span.val().abs() > NEGLIGIBLE) || !span.val().is_finite() {
        return Some(());
    }
    let zero = T::cst(0.0);
    let one = T::cst(1.0);
    let (low, high) = if a.val() <= b.val() { (a, b) } else { (b, a) };
    let lo = if low.val() > 0.0 { low } else { zero };
    let hi = if high.val() < 1.0 { high } else { one };
    let width = hi - lo;
    if !(width.val() > 0.0) {
        return Some(());
    }
    let normal = [t[1], -t[0], zero];
    let size = width * length;
    // Three-point Gauss–Legendre mapped from [−1, 1] onto the overlap.
    for &(g, w) in crate::quadrature::gauss_legendre(3) {
        let xi = lo + width * (0.5 * (1.0 + g));
        let eta = (xi - a) / span;
        let xs = [s0[0] + (s1[0] - s0[0]) * xi, s0[1] + (s1[1] - s0[1]) * xi, zero];
        let xm = [m0[0] + (m1[0] - m0[0]) * eta, m0[1] + (m1[1] - m0[1]) * eta, zero];
        f(&Sample {
            weight: size * measure(extent, xs[0]) * (0.5 * w),
            secondary: [one - xi, xi, zero, zero],
            main: [one - eta, eta, zero, zero],
            xs,
            xm,
            normal,
        });
    }
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dual::{Dual, Hyper, VARS};

    fn square(x0: f64, y0: f64, side: f64, z: f64, up: bool) -> Vec<[f64; 3]> {
        let c = [[x0, y0, z], [x0 + side, y0, z], [x0 + side, y0 + side, z], [x0, y0 + side, z]];
        if up { c.to_vec() } else { vec![c[0], c[3], c[2], c[1]] }
    }

    /// ∫ N_a over the overlap and ∫ N_a g, per secondary corner.
    fn integrals(s: &[[f64; 3]], m: &[[f64; 3]]) -> ([f64; 4], [f64; 4], f64) {
        let (mut d, mut g, mut area) = ([0.0; 4], [0.0; 4], 0.0);
        integrate(s, m, |p| {
            area += p.weight;
            for a in 0..4 {
                d[a] += p.weight * p.secondary[a];
                g[a] += p.weight * p.secondary[a] * p.gap();
            }
        })
        .unwrap();
        (d, g, area)
    }

    #[test]
    fn a_covered_face_integrates_its_shape_functions_and_gap() {
        // Secondary face up (normal +z), main face facing it 0.1 above.
        let s = square(0.0, 0.0, 1.0, 0.0, true);
        let m = square(-0.5, -0.5, 2.0, 0.1, false);
        let (d, g, area) = integrals(&s, &m);
        assert!((area - 1.0).abs() < 1e-14);
        for a in 0..4 {
            assert!((d[a] - 0.25).abs() < 1e-14 && (g[a] - 0.025).abs() < 1e-14, "{d:?} {g:?}");
        }
        // Half covered: the overlap is the strip x > 0.5.
        let (d, _, area) = integrals(&s, &square(0.5, -0.5, 2.0, 0.1, false));
        assert!((area - 0.5).abs() < 1e-14);
        assert!((d[0] - 0.0625).abs() < 1e-14 && (d[1] - 0.1875).abs() < 1e-14, "{d:?}");
        // Apart: nothing.
        let (_, _, area) = integrals(&s, &square(3.0, 0.0, 1.0, 0.1, false));
        assert_eq!(area, 0.0);
    }

    /// The sums over a 3 × 3 grid of main faces (triangles in one corner) under secondary face
    /// `s`, against a tensor-product Gauss rule on the face itself.
    fn covered(s: &[[f64; 3]; 4]) -> ([f64; 4], [f64; 4], f64, f64) {
        let mut d = [0.0; 4];
        let mut area = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                let quad = square(-0.3 + 0.6 * i as f64, -0.3 + 0.6 * j as f64, 0.6, 0.2, false);
                let faces: Vec<Vec<[f64; 3]>> =
                    if i == 0 && j == 0 { vec![vec![quad[0], quad[1], quad[2]], vec![quad[0], quad[2], quad[3]]] } else { vec![quad] };
                for f in faces {
                    let (da, _, a) = integrals(s, &f);
                    area += a;
                    for k in 0..4 {
                        d[k] += da[k];
                    }
                }
            }
        }
        let mut reference = [0.0; 4];
        let mut exact = 0.0;
        let g = [-(0.6f64).sqrt(), 0.0, 0.6f64.sqrt()];
        let w = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        for (a, wa) in g.iter().zip(w) {
            for (b, wb) in g.iter().zip(w) {
                let n = [(1.0 - a) * (1.0 - b), (1.0 + a) * (1.0 - b), (1.0 + a) * (1.0 + b), (1.0 - a) * (1.0 + b)].map(|v| v / 4.0);
                let dxa: [f64; 3] =
                    std::array::from_fn(|k| 0.25 * (-(1.0 - b) * s[0][k] + (1.0 - b) * s[1][k] + (1.0 + b) * s[2][k] - (1.0 + b) * s[3][k]));
                let dxb: [f64; 3] =
                    std::array::from_fn(|k| 0.25 * (-(1.0 - a) * s[0][k] - (1.0 + a) * s[1][k] + (1.0 + a) * s[2][k] + (1.0 - a) * s[3][k]));
                let j = crate::frame::norm(crate::frame::cross(dxa, dxb));
                exact += wa * wb * j;
                for k in 0..4 {
                    reference[k] += wa * wb * j * n[k];
                }
            }
        }
        (d, reference, area, exact)
    }

    #[test]
    fn non_matching_main_faces_add_up_to_the_covered_integrals() {
        // A parallelogram: its shape functions are polynomials of the plane, integrated exactly.
        let (d, reference, area, exact) = covered(&[[0.0, 0.0, 0.0], [1.0, 0.2, 0.0], [1.1, 1.1, 0.0], [0.1, 0.9, 0.0]]);
        assert!((area - exact).abs() < 1e-13, "{area} {exact}");
        for k in 0..4 {
            assert!((d[k] - reference[k]).abs() < 1e-13, "{d:?} {reference:?}");
        }
        // A general quadrilateral: its inverse map is not polynomial; the rule is accurate to
        // its degree.
        let (d, reference, area, exact) = covered(&[[0.0, 0.0, 0.0], [1.1, 0.1, 0.0], [1.0, 0.9, 0.0], [0.05, 1.2, 0.0]]);
        assert!((area - exact).abs() < 1e-13, "{area} {exact}");
        for k in 0..4 {
            assert!((d[k] - reference[k]).abs() < 1e-5 * reference[k], "{d:?} {reference:?}");
        }
    }

    /// ∫ N_a g over a warped, tilted pair, as a function of all 24 coordinates.
    fn weighted_gap<T: Real>(x: &[T; VARS], a: usize) -> T {
        let s: Vec<Point<T>> = (0..4).map(|i| [x[3 * i], x[3 * i + 1], x[3 * i + 2]]).collect();
        let m: Vec<Point<T>> = (0..4).map(|j| [x[12 + 3 * j], x[13 + 3 * j], x[14 + 3 * j]]).collect();
        let mut g = T::cst(0.0);
        integrate(&s, &m, |p| g = g + p.weight * p.secondary[a] * p.gap()).unwrap();
        g
    }

    #[test]
    fn the_automatic_derivatives_are_those_of_the_integral() {
        let base: [f64; VARS] = [
            0.0, 0.0, 0.0, 1.0, 0.05, 0.02, 1.1, 0.9, -0.03, -0.05, 1.0, 0.01, // secondary
            -0.2, -0.3, 0.12, -0.25, 1.15, 0.08, 1.2, 1.1, 0.1, 0.9, -0.2, 0.15, // main, facing down
        ];
        let a = 2;
        let duals: [Dual; VARS] = std::array::from_fn(|k| Dual::var(base[k], k));
        let hypers: [Hyper; VARS] = std::array::from_fn(|k| Hyper::var(base[k], k));
        let d = weighted_gap(&duals, a);
        let e = weighted_gap(&hypers, a);
        let v = weighted_gap(&base, a);
        assert!((d.v - v).abs() < 1e-14 && (e.v - v).abs() < 1e-14);
        let h = 1e-6;
        let shifted = |k: usize, by: f64| {
            let mut x = base;
            x[k] += by;
            x
        };
        for k in 0..VARS {
            let fd = (weighted_gap(&shifted(k, h), a) - weighted_gap(&shifted(k, -h), a)) / (2.0 * h);
            assert!((d.g[k] - fd).abs() < 1e-7 && (e.g[k] - d.g[k]).abs() < 1e-13, "gradient {k}: {} {fd}", d.g[k]);
            // The Hessian's row k against differences of the gradient.
            let (plus, minus) = (shifted(k, h), shifted(k, -h));
            let gp = weighted_gap(&std::array::from_fn::<Dual, VARS, _>(|i| Dual::var(plus[i], i)), a);
            let gm = weighted_gap(&std::array::from_fn::<Dual, VARS, _>(|i| Dual::var(minus[i], i)), a);
            for j in 0..VARS {
                let fd = (gp.g[j] - gm.g[j]) / (2.0 * h);
                assert!((e.hessian(k, j) - fd).abs() < 1e-6, "hessian ({k}, {j}): {} {fd}", e.hessian(k, j));
            }
        }
        // A rigid translation changes nothing: the gradient sums to zero per direction.
        for c in 0..3 {
            let sum: f64 = (0..8).map(|n| d.g[3 * n + c]).sum();
            assert!(sum.abs() < 1e-13, "{sum}");
        }
    }

    #[test]
    fn edges_integrate_their_overlap_and_measure() {
        // Secondary edge along +x at y = 0 (normal −y), main edge 0.1 below, half overlapping.
        let s = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let m = [[2.5, -0.1, 0.0], [1.5, -0.1, 0.0]];
        for (extent, at) in [(Extent::Thickness(3.0), (3.0, 0.0)), (Extent::Revolved, (0.0, 2.0 * std::f64::consts::PI))] {
            let at = |x: f64| at.0 + at.1 * x;
            let (mut d, mut g, mut area) = ([0.0; 2], 0.0, 0.0);
            integrate_edges(&s, &m, extent, |p| {
                area += p.weight;
                g += p.weight * p.gap();
                for a in 0..2 {
                    d[a] += p.weight * p.secondary[a];
                }
                assert!((p.xm[1] + 0.1).abs() < 1e-15 && (p.xm[0] - p.xs[0]).abs() < 1e-15);
            })
            .unwrap();
            // ∫_{1.5}^{2} m(x) dx and ∫ (x − 1) m(x) dx (the second corner's share).
            let simpson = |f: &dyn Fn(f64) -> f64| (f(1.5) + 4.0 * f(1.75) + f(2.0)) * 0.5 / 6.0;
            assert!((area - simpson(&|x| at(x))).abs() < 1e-12, "{area}");
            assert!((d[1] - simpson(&|x| (x - 1.0) * at(x))).abs() < 1e-12);
            assert!((g - 0.1 * area).abs() < 1e-12);
        }
        // The area vector points along −y with the edge's measure.
        let v = edge_area_vector(&s, Extent::Revolved).unwrap();
        assert!(v[0].abs() < 1e-15 && (v[1] + 2.0 * std::f64::consts::PI * 1.5).abs() < 1e-12);
        // Apart, perpendicular or degenerate: nothing.
        let mut count = 0;
        integrate_edges(&s, &[[3.0, -0.1, 0.0], [4.0, -0.1, 0.0]], Extent::Revolved, |_| count += 1).unwrap();
        integrate_edges(&s, &[[1.5, -0.1, 0.0], [1.5, -1.0, 0.0]], Extent::Revolved, |_| count += 1).unwrap();
        assert_eq!(count, 0);
        assert!(integrate_edges(&[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], &m, Extent::Revolved, |_| {}).is_none());
    }

    #[test]
    fn edge_derivatives_are_those_of_the_integral() {
        let base: [f64; VARS] = std::array::from_fn(|k| match k {
            0 => 1.0,
            1 => 0.02,
            3 => 2.1,
            4 => -0.05,
            12 => 2.4,
            13 => -0.12,
            15 => 1.3,
            16 => -0.08,
            _ => 0.0,
        });
        let gap = |x: &[Hyper; VARS]| {
            let s = [[x[0], x[1], x[2]], [x[3], x[4], x[5]]];
            let m = [[x[12], x[13], x[14]], [x[15], x[16], x[17]]];
            let mut g = Hyper::cst(0.0);
            integrate_edges(&s, &m, Extent::Revolved, |p| g = g + p.weight * p.secondary[1] * p.gap()).unwrap();
            g
        };
        let value = |x: [f64; VARS]| gap(&std::array::from_fn(|k| Hyper::cst(x[k]))).v;
        let e = gap(&std::array::from_fn(|k| Hyper::var(base[k], k)));
        let h = 1e-6;
        for k in [0, 1, 3, 4, 12, 13, 15, 16] {
            let (mut p, mut m) = (base, base);
            p[k] += h;
            m[k] -= h;
            let fd = (value(p) - value(m)) / (2.0 * h);
            assert!((e.g[k] - fd).abs() < 1e-7, "gradient {k}: {} {fd}", e.g[k]);
            for j in [0, 1, 3, 4, 12, 13, 15, 16] {
                let (mut pp, mut pm, mut mp, mut mm) = (p, p, m, m);
                pp[j] += h;
                pm[j] -= h;
                mp[j] += h;
                mm[j] -= h;
                let fd = (value(pp) - value(pm) - value(mp) + value(mm)) / (4.0 * h * h);
                assert!((e.hessian(k, j) - fd).abs() < 2e-3, "hessian ({k}, {j}): {} {fd}", e.hessian(k, j));
            }
        }
    }

    #[test]
    fn degenerate_faces_give_nothing() {
        let s = square(0.0, 0.0, 1.0, 0.0, true);
        let collapsed = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        assert!(integrate(&collapsed, &s, |_| {}).is_none());
        assert!(integrate(&s, &s[..2], |_| {}).is_none());
        let mut count = 0;
        integrate(&s, &[[f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], |_| count += 1);
        assert_eq!(count, 0);
    }
}
