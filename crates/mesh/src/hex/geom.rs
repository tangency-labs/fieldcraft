//! Geometry the hex meshers place nodes on: the analytic surface of a B-rep face (exact closest
//! point), the triangles of the boundary surface (closest point through a bounding-volume tree,
//! for free-form faces), and the inside test by the generalized winding number.
//!
//! Closest points on the elementary surfaces are closed form. The point-to-triangle distance is
//! the Voronoi-region method of C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann 2005,
//! §5.1.5. The winding number is the sum of the solid angles the triangles subtend (A. Van
//! Oosterom, J. Strackee, *The solid angle of a plane triangle*, IEEE Trans. Biomed. Eng.
//! 30(2):125–126, 1983; A. Jacobson, L. Kavan, O. Sorkine-Hornung, *Robust inside-outside
//! segmentation using generalized winding numbers*, ACM TOG 32(4):33, 2013).
//!
//! These are floating-point measures for placing nodes and choosing tags; no topological decision
//! of the meshers rests on a bare sign here except the inside test, which is only asked at points
//! far from the boundary (region centres), where the winding number is close to 0 or 1.

use crate::surface::BoundarySurface;
use crate::vec3::{P, cross, dot, norm, sub};
use crate::{Error, Result};

pub(crate) fn add(a: P, b: P) -> P {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(crate) fn scale(a: P, s: f64) -> P {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(crate) fn dist(a: P, b: P) -> f64 {
    norm(sub(a, b))
}
/// The unit vector along `a`, or `None` when it has no direction.
pub(crate) fn unit(a: P) -> Option<P> {
    let n = norm(a);
    (n > 0.0 && n.is_finite()).then(|| scale(a, 1.0 / n))
}

/// The analytic surface a B-rep face lies on. Directions need not be unit length; they are
/// normalized when the surface is checked ([`Surface::checked`]). Lengths are model units and
/// angles radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Surface {
    Plane {
        point: [f64; 3],
        normal: [f64; 3],
    },
    Cylinder {
        origin: [f64; 3],
        axis: [f64; 3],
        radius: f64,
    },
    /// The apex, the axis pointing into the cone, and the half-angle (0, π/2).
    Cone {
        apex: [f64; 3],
        axis: [f64; 3],
        half_angle: f64,
    },
    Sphere {
        center: [f64; 3],
        radius: f64,
    },
    /// The centre of the ring, its axis, the ring (major) and tube (minor) radii.
    Torus {
        center: [f64; 3],
        axis: [f64; 3],
        major: f64,
        minor: f64,
    },
}

impl Surface {
    /// The surface with unit directions, or `None` when a number is not finite or a direction or
    /// radius is degenerate (such a face is treated as free-form).
    pub(crate) fn checked(&self) -> Option<Surface> {
        let fin = |p: &[f64; 3]| p.iter().all(|x| x.is_finite());
        let pos = |r: f64| r > 0.0 && r.is_finite();
        Some(match *self {
            Surface::Plane { point, normal } if fin(&point) => Surface::Plane { point, normal: unit(normal)? },
            Surface::Cylinder { origin, axis, radius } if fin(&origin) && pos(radius) => Surface::Cylinder { origin, axis: unit(axis)?, radius },
            Surface::Cone { apex, axis, half_angle } if fin(&apex) && half_angle > 0.0 && half_angle < std::f64::consts::FRAC_PI_2 => {
                Surface::Cone { apex, axis: unit(axis)?, half_angle }
            }
            Surface::Sphere { center, radius } if fin(&center) && pos(radius) => Surface::Sphere { center, radius },
            Surface::Torus { center, axis, major, minor } if fin(&center) && pos(major) && pos(minor) => {
                Surface::Torus { center, axis: unit(axis)?, major, minor }
            }
            _ => return None,
        })
    }

    /// The closest point on the (untrimmed) surface, or `None` where it is not unique (on a
    /// cylinder's axis, at a sphere's centre, …).
    pub(crate) fn project(&self, p: P) -> Option<P> {
        let q = match *self {
            Surface::Plane { point, normal } => sub(p, scale(normal, dot(sub(p, point), normal))),
            Surface::Cylinder { origin, axis, radius } => {
                let w = sub(p, origin);
                let h = dot(w, axis);
                let radial = unit(sub(w, scale(axis, h)))?;
                add(add(origin, scale(axis, h)), scale(radial, radius))
            }
            Surface::Cone { apex, axis, half_angle } => {
                let w = sub(p, apex);
                let h = dot(w, axis);
                let radial = unit(sub(w, scale(axis, h)))?;
                let rho = norm(sub(w, scale(axis, h)));
                let (c, s) = (libm::cos(half_angle), libm::sin(half_angle));
                let t = (h * c + rho * s).max(0.0);
                add(apex, add(scale(axis, t * c), scale(radial, t * s)))
            }
            Surface::Sphere { center, radius } => add(center, scale(unit(sub(p, center))?, radius)),
            Surface::Torus { center, axis, major, minor } => {
                let w = sub(p, center);
                let h = dot(w, axis);
                let ring = add(center, scale(unit(sub(w, scale(axis, h)))?, major));
                add(ring, scale(unit(sub(p, ring))?, minor))
            }
        };
        q.iter().all(|x| x.is_finite()).then_some(q)
    }

    /// The unit normal of the surface at its point closest to `p`, up to sign, or `None` where
    /// the closest point is not unique.
    pub(crate) fn normal(&self, p: P) -> Option<P> {
        match *self {
            Surface::Plane { normal, .. } => unit(normal),
            Surface::Cylinder { origin, axis, .. } => {
                let w = sub(p, origin);
                unit(sub(w, scale(axis, dot(w, axis))))
            }
            Surface::Cone { apex, axis, half_angle } => {
                let w = sub(p, apex);
                let radial = unit(sub(w, scale(axis, dot(w, axis))))?;
                unit(sub(scale(radial, libm::cos(half_angle)), scale(axis, libm::sin(half_angle))))
            }
            Surface::Sphere { center, .. } => unit(sub(p, center)),
            Surface::Torus { center, axis, major, .. } => {
                let w = sub(p, center);
                let ring = add(center, scale(unit(sub(w, scale(axis, dot(w, axis))))?, major));
                unit(sub(p, ring))
            }
        }
    }
}

/// The closest point to `p` on triangle `abc`.
pub(crate) fn closest_on_triangle(p: P, a: P, b: P, c: P) -> P {
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
        return add(a, scale(ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return add(a, scale(ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return add(b, scale(sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))));
    }
    let den = va + vb + vc;
    if !(den.abs() > 0.0) {
        return a;
    }
    add(a, add(scale(ab, vb / den), scale(ac, vc / den)))
}

#[derive(Debug, Clone, Copy)]
struct Bounds {
    lo: P,
    hi: P,
}

impl Bounds {
    const EMPTY: Bounds = Bounds { lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3] };
    fn add(&mut self, p: P) {
        for ((l, h), x) in self.lo.iter_mut().zip(self.hi.iter_mut()).zip(p) {
            *l = l.min(x);
            *h = h.max(x);
        }
    }
    /// Squared distance from `p` to the box (0 inside).
    fn dist2(&self, p: P) -> f64 {
        (0..3).map(|i| (self.lo[i] - p[i]).max(0.0).max(p[i] - self.hi[i]).powi(2)).sum()
    }
}

enum Node {
    Leaf { bounds: Bounds, first: usize, count: usize },
    Inner { bounds: Bounds, left: usize, right: usize },
}

impl Node {
    fn bounds(&self) -> &Bounds {
        match self {
            Node::Leaf { bounds, .. } | Node::Inner { bounds, .. } => bounds,
        }
    }
}

const LEAF_SIZE: usize = 8;
/// Nodes a closest-point query may visit before it settles for the best triangle found.
const MAX_QUERY_VISITS: usize = 1 << 20;

/// A bounding-volume tree over the boundary surface's triangles (median splits along the longest
/// axis), for closest-point queries restricted to a set of B-rep faces.
pub(crate) struct TriangleTree<'a> {
    surface: &'a BoundarySurface,
    order: Vec<u32>,
    nodes: Vec<Node>,
    normals: Vec<P>,
}

/// The closest triangle found by [`TriangleTree::closest`].
#[derive(Debug, Clone, Copy)]
pub(crate) struct Hit {
    pub triangle: u32,
    pub point: P,
    pub dist: f64,
}

impl<'a> TriangleTree<'a> {
    pub(crate) fn new(surface: &'a BoundarySurface) -> Result<Self> {
        let n = surface.triangles().len();
        let mut tree = TriangleTree { surface, order: (0..n as u32).collect(), nodes: Vec::new(), normals: Vec::with_capacity(n) };
        let mut centres = Vec::with_capacity(n);
        for t in 0..n as u32 {
            let [a, b, c] = tree.corners(t)?;
            centres.push(scale(add(add(a, b), c), 1.0 / 3.0));
            tree.normals.push(unit(cross(sub(b, a), sub(c, a))).unwrap_or([0.0; 3]));
        }
        if n > 0 {
            tree.build(0, n, &centres, 0)?;
        }
        Ok(tree)
    }

    pub(crate) fn corners(&self, t: u32) -> Result<[P; 3]> {
        let pts = self.surface.points();
        let tri = self.surface.triangles().get(t as usize).ok_or_else(|| Error::Failed(format!("no triangle {t}")))?;
        let p = |v: u32| pts.get(v as usize).copied().ok_or_else(|| Error::Failed(format!("no vertex {v}")));
        Ok([p(tri[0])?, p(tri[1])?, p(tri[2])?])
    }

    /// The unit normal of triangle `t` (outward), zero when degenerate.
    pub(crate) fn normal(&self, t: u32) -> P {
        self.normals.get(t as usize).copied().unwrap_or([0.0; 3])
    }

    fn build(&mut self, first: usize, count: usize, centres: &[P], depth: usize) -> Result<usize> {
        let mut bounds = Bounds::EMPTY;
        let mut cb = Bounds::EMPTY;
        for &t in self.order.get(first..first + count).unwrap_or(&[]) {
            for p in self.corners(t)? {
                bounds.add(p);
            }
            cb.add(centres.get(t as usize).copied().unwrap_or([0.0; 3]));
        }
        let id = self.nodes.len();
        // Depth is bounded by log2 of the count (median splits); the cap only guards the stack.
        if count <= LEAF_SIZE || depth > 64 {
            self.nodes.push(Node::Leaf { bounds, first, count });
            return Ok(id);
        }
        let ext = [cb.hi[0] - cb.lo[0], cb.hi[1] - cb.lo[1], cb.hi[2] - cb.lo[2]];
        let axis = if ext[0] >= ext[1] && ext[0] >= ext[2] {
            0
        } else if ext[1] >= ext[2] {
            1
        } else {
            2
        };
        let key = |t: &u32| centres.get(*t as usize).map_or(0.0, |c| c[axis]);
        if let Some(slice) = self.order.get_mut(first..first + count) {
            slice.sort_by(|a, b| key(a).total_cmp(&key(b)).then(a.cmp(b)));
        }
        self.nodes.push(Node::Leaf { bounds, first, count });
        let half = count / 2;
        let left = self.build(first, half, centres, depth + 1)?;
        let right = self.build(first + half, count - half, centres, depth + 1)?;
        if let Some(n) = self.nodes.get_mut(id) {
            *n = Node::Inner { bounds, left, right };
        }
        Ok(id)
    }

    /// The closest triangle to `p` among those `accept` admits (by triangle index).
    pub(crate) fn closest(&self, p: P, accept: &dyn Fn(u32) -> bool) -> Option<Hit> {
        self.search(p, f64::INFINITY, accept)
    }

    /// [`TriangleTree::closest`] among the triangles within distance `r` of `p` (the search
    /// skips everything farther).
    pub(crate) fn closest_within(&self, p: P, r: f64, accept: &dyn Fn(u32) -> bool) -> Option<Hit> {
        self.search(p, r * r, accept).filter(|h| h.dist <= r)
    }

    fn search(&self, p: P, limit2: f64, accept: &dyn Fn(u32) -> bool) -> Option<Hit> {
        let mut best: Option<Hit> = None;
        let mut best2 = limit2;
        let mut stack = vec![0usize];
        let mut visits = 0usize;
        while let Some(id) = stack.pop() {
            visits += 1;
            if visits > MAX_QUERY_VISITS {
                break;
            }
            let Some(node) = self.nodes.get(id) else { continue };
            if node.bounds().dist2(p) > best2 {
                continue;
            }
            match node {
                Node::Leaf { first, count, .. } => {
                    for &t in self.order.get(*first..*first + *count).unwrap_or(&[]) {
                        if !accept(t) {
                            continue;
                        }
                        let Ok([a, b, c]) = self.corners(t) else { continue };
                        let q = closest_on_triangle(p, a, b, c);
                        let d2 = dot(sub(p, q), sub(p, q));
                        if d2 < best2 || (d2 == best2 && best.is_some_and(|h| t < h.triangle)) {
                            best2 = d2;
                            best = Some(Hit { triangle: t, point: q, dist: d2.sqrt() });
                        }
                    }
                }
                Node::Inner { left, right, .. } => {
                    let (dl, dr) = (
                        self.nodes.get(*left).map_or(f64::INFINITY, |n| n.bounds().dist2(p)),
                        self.nodes.get(*right).map_or(f64::INFINITY, |n| n.bounds().dist2(p)),
                    );
                    // Visit the nearer child first (pushed last).
                    if dl <= dr {
                        stack.push(*right);
                        stack.push(*left);
                    } else {
                        stack.push(*left);
                        stack.push(*right);
                    }
                }
            }
        }
        best
    }
}

/// The most triangle evaluations the inside tests of one mesh may spend.
pub(crate) const MAX_INSIDE_WORK: u64 = 2_000_000_000;

/// Inside tests by the generalized winding number, with a shared work budget.
pub(crate) struct Inside<'a> {
    surface: &'a BoundarySurface,
    work: u64,
}

impl<'a> Inside<'a> {
    pub(crate) fn new(surface: &'a BoundarySurface) -> Self {
        Inside { surface, work: 0 }
    }

    /// The generalized winding number of the boundary at `x` (1 inside, 0 outside).
    pub(crate) fn winding(&mut self, x: P) -> Result<f64> {
        let n = self.surface.triangles().len() as u64;
        self.work = self.work.saturating_add(n);
        if self.work > MAX_INSIDE_WORK {
            return Err(Error::Limit(format!("inside tests exceed {MAX_INSIDE_WORK} triangle evaluations")));
        }
        let pts = self.surface.points();
        let mut sum = 0.0;
        for t in self.surface.triangles() {
            let (Some(a), Some(b), Some(c)) = (pts.get(t[0] as usize), pts.get(t[1] as usize), pts.get(t[2] as usize)) else { continue };
            let (a, b, c) = (sub(*a, x), sub(*b, x), sub(*c, x));
            let (la, lb, lc) = (norm(a), norm(b), norm(c));
            let den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la;
            sum += 2.0 * crate::vec3::atan2(dot(a, cross(b, c)), den);
        }
        Ok(sum / (4.0 * std::f64::consts::PI))
    }

    /// Whether `x` is inside the solid. A winding number away from both 0 and 1 means `x` is on
    /// or next to the boundary, which callers avoid; it is an error rather than a guess.
    pub(crate) fn contains(&mut self, x: P) -> Result<bool> {
        let w = self.winding(x)?;
        if (w - 1.0).abs() < 0.25 {
            Ok(true)
        } else if w.abs() < 0.25 {
            Ok(false)
        } else {
            Err(Error::failed(crate::message![
                "cannot tell whether ",
                crate::message::point(x),
                format!(" is inside the solid (winding number {w:.3})")
            ]))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: P, b: P) -> bool {
        dist(a, b) < 1e-12
    }

    #[test]
    fn analytic_projections() {
        let cyl = Surface::Cylinder { origin: [1.0, 0.0, 0.0], axis: [0.0, 0.0, 2.0], radius: 2.0 }.checked().unwrap();
        assert!(close(cyl.project([4.0, 0.0, 5.0]).unwrap(), [3.0, 0.0, 5.0]));
        assert!(cyl.project([1.0, 0.0, 7.0]).is_none());
        let sph = Surface::Sphere { center: [0.0; 3], radius: 3.0 }.checked().unwrap();
        assert!(close(sph.project([0.0, 0.0, 1.0]).unwrap(), [0.0, 0.0, 3.0]));
        let tor = Surface::Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major: 10.0, minor: 2.0 }.checked().unwrap();
        assert!(close(tor.project([13.0, 0.0, 0.0]).unwrap(), [12.0, 0.0, 0.0]));
        assert!(close(tor.project([10.0, 0.0, 5.0]).unwrap(), [10.0, 0.0, 2.0]));
        let cone = Surface::Cone { apex: [0.0; 3], axis: [0.0, 0.0, 1.0], half_angle: std::f64::consts::FRAC_PI_4 }.checked().unwrap();
        let q = cone.project([2.0, 0.0, 0.0]).unwrap();
        assert!(close(q, [1.0, 0.0, 1.0]), "{q:?}");
        let plane = Surface::Plane { point: [0.0, 0.0, 1.0], normal: [0.0, 0.0, 5.0] }.checked().unwrap();
        assert!(close(plane.project([3.0, 4.0, 9.0]).unwrap(), [3.0, 4.0, 1.0]));
        assert!(Surface::Sphere { center: [0.0; 3], radius: f64::NAN }.checked().is_none());
        assert!(Surface::Plane { point: [0.0; 3], normal: [0.0; 3] }.checked().is_none());
    }

    #[test]
    fn triangle_closest_point_regions() {
        let (a, b, c) = ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(close(closest_on_triangle([0.2, 0.2, 3.0], a, b, c), [0.2, 0.2, 0.0]));
        assert!(close(closest_on_triangle([-1.0, -1.0, 0.0], a, b, c), a));
        assert!(close(closest_on_triangle([2.0, -1.0, 0.0], a, b, c), b));
        assert!(close(closest_on_triangle([0.5, -1.0, 0.0], a, b, c), [0.5, 0.0, 0.0]));
        assert!(close(closest_on_triangle([1.0, 1.0, 0.0], a, b, c), [0.5, 0.5, 0.0]));
    }
}
