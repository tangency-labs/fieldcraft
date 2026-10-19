//! Where boundary vertices may go: projection onto the B-rep faces and the feature curves
//! between them.
//!
//! The quality stage inserts and moves boundary vertices; each one stays on its B-rep face (or
//! feature curve) by projection. [`BoundaryGeometry`] answers the projection. Two sources:
//! - [`AnalyticSurface`]: the exact surface of a face when the geometry kernel knows it (plane,
//!   cylinder, cone, sphere, torus), with closed-form closest points;
//! - [`Tessellation`]: the boundary surface's own triangles of that face, the fallback for
//!   free-form faces (the projection is then accurate to the tessellation's chordal tolerance).
//!
//! A point on a feature curve between faces `A` and `B` is projected onto the curve's polyline
//! in the tessellation, then refined by alternating projections onto `A` and `B` when both are
//! analytic (von Neumann's alternating projections converge to the intersection for surfaces
//! meeting transversally); the result is kept only when it stays close to the polyline.

use crate::size::{LocalSize, Region, SizeField};
use crate::surface::BoundarySurface;

use super::measure::{self, P, add, cross, dot, norm, scale, sub};

/// A cone's half-angle cosine is at least this.
const MIN_COSINE: f64 = 1e-12; // unit-factors: the cosine of an angle
/// A projection within this of a point, relative to its coordinates, is on it.
const ROUND_OFF: f64 = 1e-9; // unit-factors: relative to the coordinates' size
/// A curve turning less than this is straight.
const MIN_TURN: f64 = 1e-6; // unit-factors: an angle in radians

/// A face's exact surface. Axes are unit vectors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AnalyticSurface {
    Plane {
        normal: P,
        point: P,
    },
    Cylinder {
        origin: P,
        axis: P,
        radius: f64,
    },
    /// Apex, axis pointing into the cone, half-angle in radians.
    Cone {
        apex: P,
        axis: P,
        half_angle: f64,
    },
    Sphere {
        center: P,
        radius: f64,
    },
    Torus {
        center: P,
        axis: P,
        major: f64,
        minor: f64,
    },
}

fn unit(v: P) -> Option<P> {
    let l = norm(v);
    (l > 0.0 && l.is_finite()).then(|| scale(v, 1.0 / l))
}

/// `(sin x, cos x)` for `|x| ≤ π` from their Taylor series (to 20 terms, accurate to rounding),
/// so the result is the same on every platform (the C library's `sin` and `cos` are not
/// correctly rounded and differ between platforms).
pub fn sin_cos(x: f64) -> (f64, f64) {
    let x = x.clamp(-std::f64::consts::PI, std::f64::consts::PI);
    let (mut s, mut c) = (0.0, 0.0);
    let (mut ts, mut tc) = (x, 1.0);
    for k in 0..20 {
        s += ts;
        c += tc;
        let n = 2.0 * k as f64;
        ts *= -x * x / ((n + 2.0) * (n + 3.0));
        tc *= -x * x / ((n + 1.0) * (n + 2.0));
    }
    (s, c)
}

/// Any unit vector perpendicular to `a`.
fn perpendicular(a: P) -> P {
    let t = if a[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    unit(cross(a, t)).unwrap_or([0.0, 0.0, 1.0])
}

impl AnalyticSurface {
    /// The closest point on the (untrimmed) surface, or `None` where it is not unique (on an
    /// axis) or the surface's data is not finite.
    pub fn project(&self, p: P) -> Option<P> {
        let out = match *self {
            AnalyticSurface::Plane { normal, point } => {
                let n = unit(normal)?;
                sub(p, scale(n, dot(sub(p, point), n)))
            }
            AnalyticSurface::Cylinder { origin, axis, radius } => {
                let a = unit(axis)?;
                let d = sub(p, origin);
                let foot = add(origin, scale(a, dot(d, a)));
                let r = unit(sub(p, foot))?;
                add(foot, scale(r, radius))
            }
            AnalyticSurface::Sphere { center, radius } => add(center, scale(unit(sub(p, center))?, radius)),
            AnalyticSurface::Cone { apex, axis, half_angle } => {
                let a = unit(axis)?;
                let d = sub(p, apex);
                let along = dot(d, a);
                let radial = sub(d, scale(a, along));
                let r = unit(radial).unwrap_or_else(|| perpendicular(a));
                // The generator through p's half-plane: direction g = cos·a + sin·r.
                let (s, c) = sin_cos(half_angle);
                let g = add(scale(a, c), scale(r, s));
                let t = dot(d, g).max(0.0);
                add(apex, scale(g, t))
            }
            AnalyticSurface::Torus { center, axis, major, minor } => {
                let a = unit(axis)?;
                let d = sub(p, center);
                let radial = sub(d, scale(a, dot(d, a)));
                let r = unit(radial)?;
                let ring = add(center, scale(r, major));
                add(ring, scale(unit(sub(p, ring))?, minor))
            }
        };
        out.iter().all(|x| x.is_finite()).then_some(out)
    }

    /// Move a point of the surface by the displacement `d` in the surface's own parameters (UV):
    /// `p` is first put on the surface, `d` is resolved into the parameter steps along the
    /// surface's coordinate directions there (the tangential part of `d`; the normal part is
    /// dropped), and the surface is evaluated at the new parameters, so the result lies exactly
    /// on the surface. Parameters are local to `p` (angles measured from it), so the chart has no
    /// seam. A plane translates in its own plane; a cylinder, cone and torus turn about their
    /// axis (and a torus round its tube) by the arc length over the radius; a sphere follows the
    /// great circle in the direction of `d`. Angular steps are kept within half a turn. `None`
    /// where the parameters are singular (on an axis, at an apex) or the data is not finite.
    pub fn slide(&self, p: P, d: P) -> Option<P> {
        let q = self.project(p)?;
        let rotate = |r: P, t: P, angle: f64| {
            let (s, c) = sin_cos(angle.clamp(-std::f64::consts::PI, std::f64::consts::PI));
            add(scale(r, c), scale(t, s))
        };
        let out = match *self {
            AnalyticSurface::Plane { normal, .. } => {
                let n = unit(normal)?;
                add(q, sub(d, scale(n, dot(d, n))))
            }
            AnalyticSurface::Cylinder { origin, axis, radius } => {
                let a = unit(axis)?;
                let foot = add(origin, scale(a, dot(sub(q, origin), a)));
                let r = unit(sub(q, foot))?;
                let t = cross(a, r);
                let radial = rotate(r, t, dot(d, t) / radius);
                add(add(foot, scale(a, dot(d, a))), scale(radial, radius))
            }
            AnalyticSurface::Sphere { center, radius } => {
                let r = unit(sub(q, center))?;
                let tangent = sub(d, scale(r, dot(d, r)));
                let Some(t) = unit(tangent) else { return Some(q) };
                add(center, scale(rotate(r, t, norm(tangent) / radius), radius))
            }
            AnalyticSurface::Cone { apex, axis, half_angle } => {
                let a = unit(axis)?;
                let g0 = sub(q, apex);
                let along = norm(g0);
                let g = unit(g0)?;
                let r = unit(sub(g0, scale(a, dot(g0, a))))?;
                let t = cross(a, r);
                let rho = norm(sub(g0, scale(a, dot(g0, a))));
                let s = along + dot(d, g);
                if !(s > 0.0 && rho > 0.0) {
                    return None;
                }
                let r2 = rotate(r, t, dot(d, t) / rho);
                let (sa, ca) = sin_cos(half_angle);
                add(apex, scale(add(scale(a, ca), scale(r2, sa)), s))
            }
            AnalyticSurface::Torus { center, axis, major, minor } => {
                let a = unit(axis)?;
                let c = sub(q, center);
                let radial = sub(c, scale(a, dot(c, a)));
                let rho = norm(radial);
                let r = unit(radial)?;
                let ring = add(center, scale(r, major));
                let m = unit(sub(q, ring))?;
                let around = cross(a, r);
                let tube = cross(around, m);
                // Round the tube, then about the axis.
                let m2 = rotate(m, tube, dot(d, tube) / minor);
                let r2 = rotate(r, around, dot(d, around) / rho.max(f64::MIN_POSITIVE));
                let m3 = add(scale(r2, dot(m2, r)), scale(a, dot(m2, a)));
                add(add(center, scale(r2, major)), scale(m3, minor))
            }
        };
        out.iter().all(|x| x.is_finite()).then_some(out)
    }

    /// The unit normal (up to sign) at a point on the surface, if defined there.
    pub fn normal_at(&self, q: P) -> Option<P> {
        let n = match *self {
            AnalyticSurface::Plane { normal, .. } => normal,
            AnalyticSurface::Cylinder { origin, axis, .. } => {
                let a = unit(axis)?;
                let d = sub(q, origin);
                sub(d, scale(a, dot(d, a)))
            }
            AnalyticSurface::Sphere { center, .. } => sub(q, center),
            AnalyticSurface::Cone { apex, axis, .. } => {
                // Normal to both the generator through q and the circle through q.
                let a = unit(axis)?;
                let g = sub(q, apex);
                let tangent = cross(a, g);
                cross(g, tangent)
            }
            AnalyticSurface::Torus { center, axis, major, .. } => {
                let a = unit(axis)?;
                let d = sub(q, center);
                let r = unit(sub(d, scale(a, dot(d, a))))?;
                sub(q, add(center, scale(r, major)))
            }
        };
        unit(n)
    }

    /// The smallest radius of curvature at `p` (a cone's depends on the distance to its axis).
    pub fn radius_at(&self, p: P) -> f64 {
        match *self {
            AnalyticSurface::Cone { apex, axis, half_angle } => {
                let Some(a) = unit(axis) else { return f64::INFINITY };
                let d = sub(p, apex);
                let rho = norm(sub(d, scale(a, dot(d, a))));
                // Normal curvature of the circle of radius rho: cos(α)/rho.
                rho / sin_cos(half_angle).1.abs().max(MIN_COSINE)
            }
            AnalyticSurface::Plane { .. } => f64::INFINITY,
            AnalyticSurface::Cylinder { radius, .. } | AnalyticSurface::Sphere { radius, .. } => radius.abs(),
            AnalyticSurface::Torus { minor, .. } => minor.abs(),
        }
    }
}

/// sin 10°: faces meeting at a smaller angle are treated as tangent along their curve.
const TRANSVERSAL_SIN: f64 = 0.173_648_177_666_930_35;

/// The projections the quality stage needs.
pub trait BoundaryGeometry {
    /// The closest point to `p` on B-rep face `face`, near `p` (within a few element sizes).
    fn project_face(&self, face: u32, p: P) -> Option<P>;
    /// The closest point to `p` on the feature curve between faces `faces` (equal for a sharp
    /// edge inside one face).
    fn project_curve(&self, faces: [u32; 2], p: P) -> Option<P>;
}

/// A bounding-box tree over the boundary's triangles (or its feature segments) for
/// closest-point queries, filtered by tag.
#[derive(Debug, Clone)]
struct Tree {
    nodes: Vec<Node>,
    items: Vec<u32>,
}

#[derive(Debug, Clone, Copy)]
struct Node {
    lo: P,
    hi: P,
    /// Children (left; right is left + 1) or a leaf over `items[start..end]`.
    left: Option<usize>,
    start: usize,
    end: usize,
}

const LEAF: usize = 6;

impl Tree {
    fn build(boxes: &[(P, P)]) -> Tree {
        Tree::build_over(boxes, (0..boxes.len() as u32).collect())
    }

    /// A tree over the items `items` (indices into `boxes`) only.
    fn build_over(boxes: &[(P, P)], mut items: Vec<u32>) -> Tree {
        let empty = Node { lo: [0.0; 3], hi: [0.0; 3], left: None, start: 0, end: 0 };
        let mut nodes = vec![empty];
        let mut work = vec![(0usize, 0usize, items.len())];
        while let Some((ni, start, end)) = work.pop() {
            let Some(slice) = items.get_mut(start..end) else { continue };
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for &i in slice.iter() {
                if let Some((a, b)) = boxes.get(i as usize) {
                    for k in 0..3 {
                        lo[k] = lo[k].min(a[k]);
                        hi[k] = hi[k].max(b[k]);
                    }
                }
            }
            let mut node = Node { lo, hi, left: None, start, end };
            if end - start > LEAF {
                let ext = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
                let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
                    0
                } else if ext[1] >= ext[2] {
                    1
                } else {
                    2
                };
                let key = |i: &u32| boxes.get(*i as usize).map_or(0.0, |(a, b)| a[axis] + b[axis]);
                let mid = (end - start) / 2;
                slice.select_nth_unstable_by(mid, |a, b| key(a).total_cmp(&key(b)).then(a.cmp(b)));
                let left = nodes.len();
                nodes.push(empty);
                nodes.push(empty);
                node.left = Some(left);
                work.push((left, start, start + mid));
                work.push((left + 1, start + mid, end));
            }
            if let Some(n) = nodes.get_mut(ni) {
                *n = node;
            }
        }
        Tree { nodes, items }
    }

    /// The item nearest to `p` among those `accept` takes, by `distance` (squared, with the
    /// point and a fourth value carried along), searching boxes closer than the best so far.
    fn nearest4(&self, p: P, accept: &dyn Fn(u32) -> bool, distance: &dyn Fn(u32) -> (f64, [f64; 4])) -> Option<[f64; 4]> {
        let box_d2 = |n: &Node| -> f64 {
            (0..3)
                .map(|k| {
                    let d = (n.lo[k] - p[k]).max(0.0).max(p[k] - n.hi[k]);
                    d * d
                })
                .sum()
        };
        // Ties go to the lowest item, so the result does not depend on how the tree is built.
        let mut best: Option<(f64, u32, [f64; 4])> = None;
        let mut stack = vec![0usize];
        let mut steps = 0usize;
        while let Some(ni) = stack.pop() {
            steps += 1;
            if steps > 4 * self.nodes.len() + 16 {
                break;
            }
            let Some(n) = self.nodes.get(ni) else { continue };
            if best.is_some_and(|(d, _, _)| box_d2(n) > d) {
                continue;
            }
            match n.left {
                Some(l) => {
                    let (a, b) = (self.nodes.get(l), self.nodes.get(l + 1));
                    let (da, db) = (a.map_or(f64::INFINITY, box_d2), b.map_or(f64::INFINITY, box_d2));
                    if da <= db {
                        stack.push(l + 1);
                        stack.push(l);
                    } else {
                        stack.push(l);
                        stack.push(l + 1);
                    }
                }
                None => {
                    for &i in self.items.get(n.start..n.end).unwrap_or_default() {
                        if !accept(i) {
                            continue;
                        }
                        let (d, q) = distance(i);
                        if best.is_none_or(|(b, bi, _)| d < b || (d == b && i < bi)) {
                            best = Some((d, i, q));
                        }
                    }
                }
            }
        }
        best.map(|(_, _, q)| q)
    }
}

/// The closest point on triangle `abc` (Ericson, *Real-Time Collision Detection*, §5.1.5).
pub fn closest_on_triangle(p: P, a: P, b: P, c: P) -> P {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(a, scale(ab, v));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add(a, scale(ac, w));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, scale(sub(c, b), w));
    }
    let den = 1.0 / (va + vb + vc);
    let (v, w) = (vb * den, vc * den);
    let q = add(a, add(scale(ab, v), scale(ac, w)));
    if q.iter().all(|x| x.is_finite()) { q } else { a }
}

/// The closest point on segment `ab`.
pub fn closest_on_segment(p: P, a: P, b: P) -> P {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    if !(l2 > 0.0) {
        return a;
    }
    let t = (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0);
    measure::lerp(a, b, t)
}

/// The boundary surface as projection target, optionally with exact surfaces per face.
#[derive(Debug, Clone)]
pub struct Tessellation {
    points: Vec<P>,
    triangles: Vec<[u32; 3]>,
    tri_face: Vec<u32>,
    tri_tree: Tree,
    /// A tree over each face's own triangles: a query for one face searches only its triangles,
    /// not every triangle nearer the point (the same nearest triangle, far fewer boxes).
    face_trees: Vec<Tree>,
    /// Feature segments: vertices and the (sorted) face pair.
    segments: Vec<([u32; 2], [u32; 2])>,
    seg_tree: Tree,
    /// A tree over each face pair's segments.
    pair_trees: std::collections::BTreeMap<[u32; 2], Tree>,
    analytic: Vec<Option<AnalyticSurface>>,
    /// Faces made of several B-rep faces (virtual topology): each triangle's B-rep face (its
    /// piece) and each piece's exact surface, for a face with no one surface of its own. Empty
    /// when no face is merged.
    tri_piece: Vec<u32>,
    piece_analytic: Vec<Option<AnalyticSurface>>,
}

/// The B-rep faces of merged faces: per triangle of the surface a [`Tessellation`] is built on
/// (as given to [`BoundarySurface::new`]) the B-rep face it lies on, and each B-rep face's exact
/// surface when known.
#[derive(Debug, Clone, Default)]
pub struct Pieces {
    pub tri_piece: Vec<u32>,
    pub analytic: Vec<Option<AnalyticSurface>>,
}

impl Tessellation {
    /// The projection target of a boundary surface; `analytic[face]` gives a face's exact
    /// surface when known.
    pub fn new(surface: &BoundarySurface, analytic: Vec<Option<AnalyticSurface>>) -> Tessellation {
        let points = surface.points().to_vec();
        let triangles = surface.triangles().to_vec();
        let tri_face = surface.tri_face().to_vec();
        let pt = |v: u32| points.get(v as usize).copied().unwrap_or([0.0; 3]);
        let boxes: Vec<(P, P)> = triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (pt(t[0]), pt(t[1]), pt(t[2]));
                ([0, 1, 2].map(|k| a[k].min(b[k]).min(c[k])), [0, 1, 2].map(|k| a[k].max(b[k]).max(c[k])))
            })
            .collect();
        let tri_tree = Tree::build(&boxes);
        let mut by_face: Vec<Vec<u32>> = Vec::new();
        for (t, &f) in tri_face.iter().enumerate() {
            let f = f as usize;
            if f >= by_face.len() {
                by_face.resize_with(f + 1, Vec::new);
            }
            if let Some(list) = by_face.get_mut(f) {
                list.push(t as u32);
            }
        }
        let face_trees = by_face.into_iter().map(|items| Tree::build_over(&boxes, items)).collect();
        let face_of = |t: u32| tri_face.get(t as usize).copied().unwrap_or(u32::MAX);
        let segments: Vec<([u32; 2], [u32; 2])> = surface
            .feature_edges()
            .iter()
            .map(|e| {
                let (fa, fb) = (face_of(e.triangles[0]), face_of(e.triangles[1]));
                (e.vertices, [fa.min(fb), fa.max(fb)])
            })
            .collect();
        let seg_boxes: Vec<(P, P)> = segments
            .iter()
            .map(|(v, _)| {
                let (a, b) = (pt(v[0]), pt(v[1]));
                ([0, 1, 2].map(|k| a[k].min(b[k])), [0, 1, 2].map(|k| a[k].max(b[k])))
            })
            .collect();
        let seg_tree = Tree::build(&seg_boxes);
        let mut by_pair: std::collections::BTreeMap<[u32; 2], Vec<u32>> = std::collections::BTreeMap::new();
        for (i, (_, pair)) in segments.iter().enumerate() {
            by_pair.entry(*pair).or_default().push(i as u32);
        }
        let pair_trees = by_pair.into_iter().map(|(pair, items)| (pair, Tree::build_over(&seg_boxes, items))).collect();
        Tessellation {
            points,
            triangles,
            tri_face,
            tri_tree,
            face_trees,
            segments,
            seg_tree,
            pair_trees,
            analytic,
            tri_piece: Vec::new(),
            piece_analytic: Vec::new(),
        }
    }

    /// The same target with the B-rep pieces of its merged faces (`surface` the boundary
    /// surface it was built on): a point on a merged face that lies on no one surface goes onto
    /// the exact surface of the piece it is on, trimmed to that piece, instead of its triangles.
    pub fn with_pieces(mut self, surface: &BoundarySurface, pieces: Pieces) -> Tessellation {
        if pieces.tri_piece.is_empty() {
            return self;
        }
        self.tri_piece = surface.source_triangle().iter().map(|&t| pieces.tri_piece.get(t as usize).copied().unwrap_or(u32::MAX)).collect();
        self.piece_analytic = pieces.analytic;
        self
    }

    fn pt(&self, v: u32) -> P {
        self.points.get(v as usize).copied().unwrap_or([0.0; 3])
    }

    /// The exact surface of a face, if known.
    pub fn analytic(&self, face: u32) -> Option<&AnalyticSurface> {
        self.analytic.get(face as usize).and_then(|a| a.as_ref())
    }

    /// The exact surface of a point of face `face` on its triangle `t`: the face's own, else
    /// that of the B-rep piece of a merged face the triangle lies on.
    fn exact(&self, face: u32, t: u32) -> Option<&AnalyticSurface> {
        self.analytic(face).or_else(|| {
            let piece = *self.tri_piece.get(t as usize)?;
            self.piece_analytic.get(piece as usize).and_then(|a| a.as_ref())
        })
    }

    /// Whether triangle `u` lies on the same piece of face `face` as triangle `t` (always, for
    /// a face with a surface of its own or without pieces).
    fn same_piece(&self, face: u32, t: u32, u: u32) -> bool {
        self.analytic(face).is_some() || self.tri_piece.get(t as usize) == self.tri_piece.get(u as usize)
    }

    /// The nearest triangle of face `face` to `p` and the closest point on it.
    fn nearest_triangle(&self, face: u32, p: P) -> Option<(P, u32)> {
        let accept = |t: u32| self.tri_face.get(t as usize) == Some(&face);
        let dist = |t: u32| {
            let tri = self.triangles.get(t as usize).copied().unwrap_or([0; 3]);
            let q = closest_on_triangle(p, self.pt(tri[0]), self.pt(tri[1]), self.pt(tri[2]));
            let d = sub(q, p);
            // The triangle index travels through the f64 exactly (counts are capped far below
            // 2^53).
            (dot(d, d), [q[0], q[1], q[2], f64::from(t)])
        };
        self.face_tree(face).nearest4(p, &accept, &dist).map(|r| ([r[0], r[1], r[2]], r[3] as u32))
    }

    /// The closest point on the face's triangles, the longest edge of the triangle it is on
    /// (the local scale of the tessellation) and that triangle.
    fn on_tessellation(&self, face: u32, p: P) -> Option<(P, f64, u32)> {
        let (q, t) = self.nearest_triangle(face, p)?;
        let tri = self.triangles.get(t as usize).copied().unwrap_or([0; 3]);
        let (a, b, c) = (self.pt(tri[0]), self.pt(tri[1]), self.pt(tri[2]));
        Some((q, measure::dist(a, b).max(measure::dist(b, c)).max(measure::dist(c, a)), t))
    }

    /// The B-rep face whose triangles lie nearest to `p`.
    pub fn nearest_face(&self, p: P) -> Option<u32> {
        let dist = |t: u32| {
            let tri = self.triangles.get(t as usize).copied().unwrap_or([0; 3]);
            let d = sub(closest_on_triangle(p, self.pt(tri[0]), self.pt(tri[1]), self.pt(tri[2])), p);
            (dot(d, d), [0.0, 0.0, 0.0, f64::from(t)])
        };
        let r = self.tri_tree.nearest4(p, &|_| true, &dist)?;
        // The triangle index travels through the f64 exactly (triangle counts are capped far
        // below 2^53).
        self.tri_face.get(r[3] as usize).copied()
    }

    /// The tree over `face`'s triangles (the whole boundary's for a face out of range, which
    /// the queries' filter then finds nothing in).
    fn face_tree(&self, face: u32) -> &Tree {
        self.face_trees.get(face as usize).unwrap_or(&self.tri_tree)
    }

    /// The closest point on the curve's polyline, and the length of the segment it is on.
    fn on_polyline(&self, faces: [u32; 2], p: P) -> Option<(P, f64)> {
        let key = [faces[0].min(faces[1]), faces[0].max(faces[1])];
        let accept = |s: u32| self.segments.get(s as usize).is_some_and(|(_, f)| *f == key);
        let dist = |s: u32| {
            let v = self.segments.get(s as usize).map_or([0; 2], |(v, _)| *v);
            let (a, b) = (self.pt(v[0]), self.pt(v[1]));
            let q = closest_on_segment(p, a, b);
            let d = sub(q, p);
            (dot(d, d), [q[0], q[1], q[2], measure::dist(a, b)])
        };
        self.pair_trees.get(&key).unwrap_or(&self.seg_tree).nearest4(p, &accept, &dist).map(|r| ([r[0], r[1], r[2]], r[3]))
    }
}

/// How far off a face's triangles, in triangle sizes, a point of its surface may lie and still
/// be on the face, beyond the triangles' own gap to the surface (rounding only).
const PLANE_SLACK: f64 = 1e-6; // unit-factors: relative to the triangles' size

impl Tessellation {
    /// The closest point to `p` on the trimmed face `face`: the exact surface's closest point
    /// when it lies on the face (no further off the face's nearest triangle than twice that
    /// triangle's own chordal gap to the surface, which is rounding only for a plane), else the
    /// closest point of the face's triangles. Unlike [`BoundaryGeometry::project_face`] it never
    /// returns a point of the surface's extension beyond the face's boundary, whatever the
    /// distance from `p` — nor of the same surface's other faces: the pieces of one fillet split
    /// by a cut share their cylinder, and a point on one piece is off the others. On a merged
    /// face lying on no one surface, the surface is that of the piece its nearest triangle is
    /// on, trimmed to that piece.
    pub fn project_trimmed(&self, face: u32, p: P) -> Option<P> {
        let corners = |t: u32| {
            let tri = self.triangles.get(t as usize).copied().unwrap_or([0; 3]);
            (self.pt(tri[0]), self.pt(tri[1]), self.pt(tri[2]))
        };
        let (near, nearest) = self.nearest_triangle(face, p)?;
        let accept = |t: u32| self.tri_face.get(t as usize) == Some(&face) && self.same_piece(face, nearest, t);
        let Some(surface) = self.exact(face, nearest) else { return Some(near) };
        let Some(q) = surface.project(p) else { return Some(near) };
        // The nearest triangle to the surface's point, with its slack: the largest gap between
        // the triangle and the surface (at its edges' midpoints and centre: a chord's deepest
        // points), doubled, and a millionth of its size for rounding.
        let gap = |x: P| surface.project(x).map_or(0.0, |y| measure::dist(x, y));
        let dist = |t: u32| {
            let (a, b, c) = corners(t);
            let o = closest_on_triangle(q, a, b, c);
            let d = sub(o, q);
            let mid = |u: P, v: P| scale(add(u, v), 0.5);
            let deepest = [mid(a, b), mid(b, c), mid(c, a), scale(add(add(a, b), c), 1.0 / 3.0)].into_iter().map(gap).fold(0.0, f64::max);
            let size = measure::dist(a, b).max(measure::dist(b, c)).max(measure::dist(c, a));
            (dot(d, d), [o[0], o[1], o[2], 2.0 * deepest + PLANE_SLACK * size])
        };
        let r = self.face_tree(face).nearest4(q, &accept, &dist)?;
        Some(if measure::dist([r[0], r[1], r[2]], q) <= r[3] { q } else { near })
    }

    /// The outward unit normal of B-rep face `face` near `p`: the exact surface's (oriented like
    /// the nearest triangle; on a merged face, its piece's) when known, else the nearest
    /// triangle's.
    pub fn face_normal(&self, face: u32, p: P) -> Option<P> {
        let (_, t) = self.nearest_triangle(face, p)?;
        let tri = self.triangles.get(t as usize).copied().unwrap_or([0; 3]);
        let (a, b, c) = (self.pt(tri[0]), self.pt(tri[1]), self.pt(tri[2]));
        let flat = unit(cross(sub(b, a), sub(c, a)))?;
        match self.exact(face, t).and_then(|s| s.normal_at(p)) {
            Some(n) if dot(n, flat) < 0.0 => Some(scale(n, -1.0)),
            Some(n) => Some(n),
            None => Some(flat),
        }
    }
}

impl BoundaryGeometry for Tessellation {
    fn project_face(&self, face: u32, p: P) -> Option<P> {
        let (near, scale_len, t) = self.on_tessellation(face, p)?;
        match self.exact(face, t) {
            // The exact surface's closest point, kept when it is near the tessellation's: the
            // tessellation deviates from its surface by far less than half its triangles' size,
            // so a point further away is a projection onto another part of the surface (the far
            // side of a cylinder) and the tessellation's point is used.
            Some(s) => {
                let q = s.project(p)?;
                let tol = (0.5 * scale_len).max(2.0 * measure::dist(p, near));
                Some(if measure::dist(q, near) <= tol { q } else { near })
            }
            None => Some(near),
        }
    }

    fn project_curve(&self, faces: [u32; 2], p: P) -> Option<P> {
        let (near, seg_len) = self.on_polyline(faces, p)?;
        let (Some(a), Some(b)) = (self.analytic(faces[0]), self.analytic(faces[1])) else { return Some(near) };
        if faces[0] == faces[1] {
            return Some(near);
        }
        // Faces meeting tangentially (a blend's edge) define their curve badly: alternating
        // projections crawl along it. There the polyline lies on both faces to the
        // tessellation's accuracy and is kept.
        let (Some(na), Some(nb)) = (a.normal_at(near), b.normal_at(near)) else { return Some(near) };
        if norm(cross(na, nb)) < TRANSVERSAL_SIN {
            return Some(near);
        }
        let mut q = near;
        for _ in 0..40 {
            let Some(qa) = a.project(q) else { return Some(near) };
            let Some(qb) = b.project(qa) else { return Some(near) };
            let step = measure::dist(qb, q);
            q = qb;
            if step <= 1e-13 * (1.0 + norm(q)) {
                break;
            }
        }
        // Both surfaces must agree at the result, and it must stay near the polyline.
        let gap = a.project(q).map_or(f64::INFINITY, |qa| measure::dist(qa, q));
        let tol = (0.5 * seg_len).max(2.0 * measure::dist(p, near));
        if gap <= ROUND_OFF * (1.0 + norm(q)) && measure::dist(q, near) <= tol { Some(q) } else { Some(near) }
    }
}

impl Tessellation {
    /// Move a point of face `face` by the displacement `d` and keep it on the face: through the
    /// face's own parameters ([`AnalyticSurface::slide`]) when its exact surface is known and the
    /// result stays near the face's triangles (not on another part of the surface), else by
    /// projecting `p + d` onto the face's triangles (free-form faces: accurate to the
    /// tessellation's chordal tolerance, as [`BoundaryGeometry::project_face`]).
    pub fn slide_face(&self, face: u32, p: P, d: P) -> Option<P> {
        let target = add(p, d);
        let (near, scale_len, t) = self.on_tessellation(face, target)?;
        if let Some(q) = self.exact(face, t).and_then(|s| s.slide(p, d)) {
            let tol = (0.5 * scale_len).max(2.0 * measure::dist(target, near));
            if measure::dist(q, near) <= tol {
                return Some(q);
            }
        }
        self.project_face(face, target)
    }

    /// Move a point of the feature curve between faces `faces` by `d` and keep it on the curve
    /// (the curve's closest point to `p + d`: the move along the curve's parameter).
    pub fn slide_curve(&self, faces: [u32; 2], p: P, d: P) -> Option<P> {
        self.project_curve(faces, add(p, d))
    }
}

impl Tessellation {
    /// The local thickness at each vertex: the distance to the nearest triangle facing it from
    /// across the solid (a triangle that does not touch the vertex and whose outward normal
    /// points back against the vertex's, by more than 120°). Infinite where nothing faces it.
    /// Thin walls are this thick; elements larger than it would have to be flat.
    pub fn thickness(&self) -> Vec<f64> {
        let n = self.points.len();
        let mut normal = vec![[0.0; 3]; n];
        let tri_normal: Vec<P> = self
            .triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (self.pt(t[0]), self.pt(t[1]), self.pt(t[2]));
                cross(sub(b, a), sub(c, a))
            })
            .collect();
        for (t, m) in self.triangles.iter().zip(&tri_normal) {
            for &v in t {
                if let Some(x) = normal.get_mut(v as usize) {
                    *x = add(*x, *m);
                }
            }
        }
        (0..n as u32)
            .map(|v| {
                let p = self.pt(v);
                let Some(nv) = normal.get(v as usize).and_then(|m| unit(*m)) else { return f64::INFINITY };
                let accept = |t: u32| {
                    let Some(tri) = self.triangles.get(t as usize) else { return false };
                    !tri.contains(&v) && tri_normal.get(t as usize).and_then(|m| unit(*m)).is_some_and(|m| dot(m, nv) < -0.5)
                };
                let dist = |t: u32| {
                    let tri = self.triangles.get(t as usize).copied().unwrap_or([0; 3]);
                    let q = closest_on_triangle(p, self.pt(tri[0]), self.pt(tri[1]), self.pt(tri[2]));
                    let d = sub(q, p);
                    (dot(d, d), [q[0], q[1], q[2], 0.0])
                };
                // Only what lies behind the surface (inside the solid) counts.
                self.tri_tree.nearest4(p, &accept, &dist).map_or(f64::INFINITY, |r| {
                    let q = [r[0], r[1], r[2]];
                    if dot(sub(q, p), nv) < 0.0 { measure::dist(q, p) } else { f64::INFINITY }
                })
            })
            .collect()
    }
}

/// The curvature of the boundary, sampled at its vertices: for each edge inside one B-rep face,
/// the turn of the normal between its two triangles over the distance between their centroids
/// (a discrete estimate of the normal curvature across the edge), or the exact surface's radius.
/// Returns the smallest radius of curvature at each surface vertex (infinite where flat).
pub fn curvature_radii(surface: &BoundarySurface, geometry: &Tessellation) -> Vec<f64> {
    let pts = surface.points();
    let tris = surface.triangles();
    let faces = surface.tri_face();
    let mut radius = vec![f64::INFINITY; pts.len()];
    let pt = |v: u32| pts.get(v as usize).copied().unwrap_or([0.0; 3]);
    let centroid = |t: &[u32; 3]| scale(add(add(pt(t[0]), pt(t[1])), pt(t[2])), 1.0 / 3.0);
    for (i, (t, nb)) in tris.iter().zip(surface.neighbors()).enumerate() {
        let face = faces.get(i).copied().unwrap_or(u32::MAX);
        if let Some(s) = geometry.analytic(face) {
            for &v in t {
                if let Some(r) = radius.get_mut(v as usize) {
                    *r = r.min(s.radius_at(pt(v)));
                }
            }
            continue;
        }
        let n = measure::triangle_normal(pt(t[0]), pt(t[1]), pt(t[2]));
        for (e, &j) in nb.iter().enumerate() {
            let Some(u) = tris.get(j as usize) else { continue };
            if j as usize <= i || faces.get(j as usize) != Some(&face) {
                continue;
            }
            let m = measure::triangle_normal(pt(u[0]), pt(u[1]), pt(u[2]));
            // The chord between the unit normals: the turn angle to second order, and only
            // correctly rounded operations (the same on every platform).
            let turn = norm(sub(n, m));
            let d = measure::dist(centroid(t), centroid(u));
            if !(turn > MIN_TURN) || !(d > 0.0) {
                continue;
            }
            let r = d / turn;
            for v in [t[e], t[(e + 1) % 3]] {
                if let Some(x) = radius.get_mut(v as usize) {
                    *x = x.min(r);
                }
            }
        }
    }
    radius
}

/// A size field from a global size and the boundary's curvature, graded: the size at `x` is
/// `min(global, min_i (h_i + g |x − x_i|))` over curvature samples `h_i = θ r_i` (θ radians per
/// element along a curve of radius `r_i`), stored on a background grid relaxed by fast sweeping
/// and interpolated trilinearly.
#[derive(Debug, Clone)]
pub struct GradedSize {
    global: f64,
    lo: P,
    cell: f64,
    dims: [usize; 3],
    values: Vec<f64>,
}

/// The most nodes of the background grid.
pub const MAX_GRID_NODES: usize = 1 << 21;

impl GradedSize {
    /// `samples`: points with their local size. `grading`: the largest change of size per unit
    /// length.
    ///
    /// `locals`: regions with their size ([`crate::size::LocalSize`]), seeded at every grid node
    /// within a cell of the region (boxes and balls: inside too) at the local size plus the
    /// grading times the distance, then graded with the rest.
    pub fn new(global: f64, samples: &[(P, f64)], locals: &[LocalSize], lo: P, hi: P, grading: f64) -> GradedSize {
        let diag = measure::dist(lo, hi).max(f64::MIN_POSITIVE);
        let finest = samples.iter().map(|s| s.1).chain(locals.iter().map(|l| l.size)).filter(|h| *h > 0.0).fold(global, f64::min).max(diag * 1e-4);
        let mut cell = finest.max(diag / 64.0) * 0.5;
        let dims_for = |c: f64| [0, 1, 2].map(|k| (((hi[k] - lo[k]) / c).ceil().max(1.0) as usize) + 2);
        let mut dims = dims_for(cell);
        while dims.iter().product::<usize>() > MAX_GRID_NODES {
            cell *= 1.25;
            dims = dims_for(cell);
        }
        let lo = [lo[0] - cell, lo[1] - cell, lo[2] - cell];
        let n = dims.iter().product::<usize>();
        let mut values = vec![global; n];
        let idx = |i: usize, j: usize, k: usize| (k * dims[1] + j) * dims[0] + i;
        let g = grading.max(0.0);
        // Seed: each sample's cell corners.
        for &(x, h) in samples {
            if !(h > 0.0 && h.is_finite()) {
                continue;
            }
            let c = [0, 1, 2].map(|k| ((x[k] - lo[k]) / cell).floor().max(0.0) as usize);
            for dk in 0..2 {
                for dj in 0..2 {
                    for di in 0..2 {
                        let (i, j, k) = (c[0] + di, c[1] + dj, c[2] + dk);
                        if i >= dims[0] || j >= dims[1] || k >= dims[2] {
                            continue;
                        }
                        let node = [lo[0] + i as f64 * cell, lo[1] + j as f64 * cell, lo[2] + k as f64 * cell];
                        let v = h + g * measure::dist(node, x);
                        if let Some(s) = values.get_mut(idx(i, j, k)) {
                            *s = s.min(v);
                        }
                    }
                }
            }
        }
        // Seed: the grid nodes near each local region.
        let reach = cell * 3f64.sqrt();
        let range = |a: f64, b: f64, k: usize| -> (usize, usize) {
            let first = ((a - reach - lo[k]) / cell).ceil().max(0.0);
            let last = ((b + reach - lo[k]) / cell).floor().min((dims[k] - 1) as f64);
            if first.is_finite() && last.is_finite() && first <= last { (first as usize, last as usize) } else { (1, 0) }
        };
        let mut seed_near = |a: P, b: P, h: f64, distance: &dyn Fn(P) -> f64| {
            let r = [0, 1, 2].map(|k| range(a[k], b[k], k));
            for k in r[2].0..=r[2].1 {
                for j in r[1].0..=r[1].1 {
                    for i in r[0].0..=r[0].1 {
                        let node = [lo[0] + i as f64 * cell, lo[1] + j as f64 * cell, lo[2] + k as f64 * cell];
                        let d = distance(node);
                        if d <= reach
                            && let Some(s) = values.get_mut(idx(i, j, k))
                        {
                            *s = s.min(h + g * d);
                        }
                    }
                }
            }
        };
        for l in locals {
            if !(l.size > 0.0 && l.size.is_finite() && l.size < global) {
                continue;
            }
            match &l.region {
                Region::Triangles(tris) => {
                    for t in tris {
                        let a = [0, 1, 2].map(|k| t[0][k].min(t[1][k]).min(t[2][k]));
                        let b = [0, 1, 2].map(|k| t[0][k].max(t[1][k]).max(t[2][k]));
                        seed_near(a, b, l.size, &|p| measure::dist(p, closest_on_triangle(p, t[0], t[1], t[2])));
                    }
                }
                region => {
                    let (a, b) = region.bounds();
                    seed_near(a, b, l.size, &|p| region.distance(p));
                }
            }
        }
        // Fast sweeping: eight sweep orders, each relaxing against the three axis neighbours
        // behind the sweep direction (and diagonals through them over passes).
        if g > 0.0 {
            for _ in 0..2 {
                for dir in 0..8 {
                    let rev = [dir & 1 != 0, dir & 2 != 0, dir & 4 != 0];
                    for kk in 0..dims[2] {
                        let k = if rev[2] { dims[2] - 1 - kk } else { kk };
                        for jj in 0..dims[1] {
                            let j = if rev[1] { dims[1] - 1 - jj } else { jj };
                            for ii in 0..dims[0] {
                                let i = if rev[0] { dims[0] - 1 - ii } else { ii };
                                let mut best = values.get(idx(i, j, k)).copied().unwrap_or(global);
                                let back = |x: usize, r: bool, n: usize| -> Option<usize> {
                                    if r { (x + 1 < n).then_some(x + 1) } else { x.checked_sub(1) }
                                };
                                if let Some(pi) = back(i, rev[0], dims[0]) {
                                    best = best.min(values.get(idx(pi, j, k)).copied().unwrap_or(global) + g * cell);
                                }
                                if let Some(pj) = back(j, rev[1], dims[1]) {
                                    best = best.min(values.get(idx(i, pj, k)).copied().unwrap_or(global) + g * cell);
                                }
                                if let Some(pk) = back(k, rev[2], dims[2]) {
                                    best = best.min(values.get(idx(i, j, pk)).copied().unwrap_or(global) + g * cell);
                                }
                                if let Some(s) = values.get_mut(idx(i, j, k)) {
                                    *s = best;
                                }
                            }
                        }
                    }
                }
            }
        }
        GradedSize { global, lo, cell, dims, values }
    }

    fn at(&self, i: usize, j: usize, k: usize) -> f64 {
        let i = i.min(self.dims[0] - 1);
        let j = j.min(self.dims[1] - 1);
        let k = k.min(self.dims[2] - 1);
        self.values.get((k * self.dims[1] + j) * self.dims[0] + i).copied().unwrap_or(self.global)
    }
}

impl SizeField for GradedSize {
    fn size(&self, p: [f64; 3]) -> f64 {
        let f = [0, 1, 2].map(|k| ((p[k] - self.lo[k]) / self.cell).clamp(0.0, (self.dims[k] - 1) as f64));
        if f.iter().any(|x| !x.is_finite()) {
            return self.global;
        }
        let c = f.map(|x| (x.floor() as usize).min(usize::MAX - 1));
        let t = [f[0] - c[0] as f64, f[1] - c[1] as f64, f[2] - c[2] as f64];
        let mut h = 0.0;
        for dk in 0..2 {
            for dj in 0..2 {
                for di in 0..2 {
                    let w = (if di == 1 { t[0] } else { 1.0 - t[0] })
                        * (if dj == 1 { t[1] } else { 1.0 - t[1] })
                        * (if dk == 1 { t[2] } else { 1.0 - t[2] });
                    h += w * self.at(c[0] + di, c[1] + dj, c[2] + dk);
                }
            }
        }
        if h > 0.0 && h.is_finite() { h.min(self.global) } else { self.global }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analytic_projections() {
        let cyl = AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 2.0], radius: 3.0 };
        let q = cyl.project([1.0, 1.0, 5.0]).unwrap();
        assert!((norm([q[0], q[1], 0.0]) - 3.0).abs() < 1e-12 && (q[2] - 5.0).abs() < 1e-12);
        let s = AnalyticSurface::Sphere { center: [1.0, 0.0, 0.0], radius: 2.0 };
        assert!((measure::dist(s.project([5.0, 0.0, 0.0]).unwrap(), [3.0, 0.0, 0.0])) < 1e-12);
        let pl = AnalyticSurface::Plane { normal: [0.0, 0.0, 1.0], point: [0.0, 0.0, 2.0] };
        assert_eq!(pl.project([1.0, 2.0, 7.0]), Some([1.0, 2.0, 2.0]));
        let t = AnalyticSurface::Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major: 5.0, minor: 1.0 };
        let q = t.project([7.0, 0.0, 0.0]).unwrap();
        assert!(measure::dist(q, [6.0, 0.0, 0.0]) < 1e-12);
        let cone = AnalyticSurface::Cone { apex: [0.0; 3], axis: [0.0, 0.0, 1.0], half_angle: std::f64::consts::FRAC_PI_4 };
        let q = cone.project([2.0, 0.0, 0.0]).unwrap();
        assert!((q[0] - q[2]).abs() < 1e-12 && (q[0] - 1.0).abs() < 1e-12, "{q:?}");
        assert!(cyl.project([0.0, 0.0, 1.0]).is_none());
        for x in [-3.0, -1.0, 0.0, 0.3, 1.5, 3.1] {
            let (s, c) = sin_cos(x);
            assert!((s - f64::sin(x)).abs() < 1e-15 && (c - f64::cos(x)).abs() < 1e-15, "{x}");
        }
    }

    #[test]
    fn closest_points() {
        let (a, b, c) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(measure::dist(closest_on_triangle([0.2, 0.2, 3.0], a, b, c), [0.2, 0.2, 0.0]) < 1e-15);
        assert_eq!(closest_on_triangle([-1.0, -1.0, 0.0], a, b, c), a);
        let q = closest_on_triangle([1.0, 1.0, 0.0], a, b, c);
        assert!(measure::dist(q, [0.5, 0.5, 0.0]) < 1e-15);
        assert_eq!(closest_on_segment([0.5, 2.0, 0.0], a, b), [0.5, 0.0, 0.0]);
    }

    #[test]
    fn graded_size_takes_local_regions() {
        let boxed = LocalSize { region: Region::Box { lo: [8.0; 3], hi: [12.0; 3] }, size: 2.0 };
        let face = LocalSize { region: Region::Triangles(vec![[[0.0, 0.0, 20.0], [20.0, 0.0, 20.0], [0.0, 20.0, 20.0]]]), size: 1.0 };
        let f = GradedSize::new(8.0, &[], &[boxed, face], [0.0; 3], [20.0, 20.0, 20.0], 0.5);
        let inside = f.size([10.0, 10.0, 10.0]);
        assert!((inside - 2.0).abs() < 0.3, "{inside}");
        let on_face = f.size([5.0, 5.0, 20.0]);
        assert!((on_face - 1.0).abs() < 0.3, "{on_face}");
        // Four units off the box: about 2 + 0.5 × 4.
        let off = f.size([10.0, 10.0, 4.0]);
        assert!((off - 4.0).abs() < 0.6, "{off}");
        assert_eq!(f.size([20.0, 20.0, 0.0]), 8.0);
    }

    #[test]
    fn graded_size_grows_away_from_samples() {
        let f = GradedSize::new(10.0, &[([0.0; 3], 1.0)], &[], [0.0; 3], [20.0, 20.0, 20.0], 0.5);
        let near = f.size([0.0; 3]);
        let mid = f.size([4.0, 0.0, 0.0]);
        let far = f.size([20.0, 20.0, 20.0]);
        assert!((1.0..1.6).contains(&near), "{near}");
        assert!((mid - 3.0).abs() < 0.8, "{mid}");
        assert_eq!(far, 10.0);
        assert!(f.size([f64::NAN, 0.0, 0.0]) > 0.0);
    }

    /// A slide in UV keeps the point exactly on every kind of surface, moves it by the tangential
    /// part of the displacement (to first order), and turns a cylinder point by arc / radius.
    #[test]
    fn slides_stay_on_the_surface() {
        let surfaces = [
            AnalyticSurface::Plane { normal: [0.0, 0.0, 1.0], point: [0.0, 0.0, 2.0] },
            AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: 3.0 },
            AnalyticSurface::Sphere { center: [1.0, 0.0, 0.0], radius: 2.0 },
            AnalyticSurface::Cone { apex: [0.0; 3], axis: [0.0, 0.0, 1.0], half_angle: 0.4 },
            AnalyticSurface::Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major: 5.0, minor: 1.0 },
        ];
        let probe = [2.5, 1.5, 2.2];
        for s in surfaces {
            let q = s.project(probe).unwrap();
            for d in [[0.01, 0.0, 0.0], [0.0, 0.02, 0.0], [0.0, 0.0, -0.015], [0.3, -0.2, 0.1]] {
                let moved = s.slide(q, d).unwrap();
                let back = s.project(moved).unwrap();
                assert!(measure::dist(moved, back) < 1e-12, "{s:?} {d:?}: {moved:?} off by {}", measure::dist(moved, back));
                let n = s.normal_at(q).unwrap();
                let tangential = sub(d, scale(n, dot(d, n)));
                if norm(d) < 0.05 {
                    assert!(measure::dist(moved, add(q, tangential)) < 0.05 * norm(d) + 1e-12, "{s:?} {d:?}");
                }
            }
        }
        let cyl = AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: 2.0 };
        let q = cyl.slide([2.0, 0.0, 1.0], [0.0, std::f64::consts::PI, 0.5]).unwrap();
        assert!(measure::dist(q, [0.0, 2.0, 1.5]) < 1e-12, "{q:?}");
        // On the axis the parameters are singular.
        assert!(cyl.slide([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).is_none());
    }
}
