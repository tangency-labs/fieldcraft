//! Voxelisation: which cell centres of a grid lie inside a closed, outward-oriented triangle
//! surface, decided with exact predicates, and which triangle a lattice link meets.
//!
//! # Inside or outside
//! A cell centre `p` is inside when the winding number of the surface around it is positive. The
//! winding number is counted along the ray from `p` in +x: every triangle the ray crosses adds
//! the sign of its normal's x component (+1 leaving the solid, −1 entering it). All centres on one
//! grid line (fixed y, z) share the triangles their ray can meet, so each triangle is visited once
//! per line it covers:
//! - whether the line meets the triangle is the sign of three [`orient2d`] tests of the line's
//!   (y, z) against the triangle's projected edges;
//! - which centres lie before the triangle's plane is the sign of [`orient3d`] of the triangle and
//!   the centre, searched from a floating-point estimate of the crossing (the sign is monotone
//!   along the line, so the estimate only decides where the exact search starts).
//!
//! Degenerate positions (a line through a projected edge or vertex, a centre on a triangle's
//! plane) are resolved by one symbolic perturbation of the centre, `p + (ε, δ, δ²)` with
//! infinitesimals `1 ≫ ε ≫ δ` (Edelsbrunner & Mücke, *Simulation of Simplicity*, ACM TOG 9(1),
//! 1990): the line moves by `(δ, δ²)` in (y, z), so a line through a shared edge or vertex meets
//! exactly one of the triangles around it, and a centre on a plane the line crosses moves off it
//! by `ε` in +x. Every test sees the same perturbed point, so a watertight surface gives a
//! consistent count and the result never depends on rounding. Triangles parallel to x are met by
//! no perturbed line and are skipped.
//!
//! # Link faces
//! Momentum-exchange forces are reported per B-rep face, so every wall link (a fluid centre next
//! to a solid one) is given the face of the triangle nearest its midpoint ([`TriangleTree`]): a
//! bounding-volume tree with branch-and-bound nearest search and the closest point on a triangle
//! after C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann 2005, §5.1.5. Ties go to the
//! lower triangle index, so the assignment is deterministic. This is a classification for
//! reporting, not a meshing decision, so floating-point distances are adequate.

use fieldcraft_mesh::predicates::{Point2, Point3, Sign, orient2d, orient3d, snap_to_domain};
use fieldcraft_mesh::surface::BoundarySurface;

use crate::{Error, Result};

/// A uniform grid of cubic cells: cell `(i, j, k)` has its centre at
/// `origin + (i + ½, j + ½, k + ½) · spacing` and index `i + nx (j + ny k)`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Grid {
    pub origin: [f64; 3],
    pub spacing: f64,
    pub cells: [usize; 3],
}

impl Grid {
    /// The number of cells.
    pub fn len(&self) -> usize {
        self.cells[0] * self.cells[1] * self.cells[2]
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The index of cell `(i, j, k)`.
    #[inline(always)]
    pub fn index(&self, [i, j, k]: [usize; 3]) -> usize {
        i + self.cells[0] * (j + self.cells[1] * k)
    }

    /// The `(i, j, k)` of a cell index.
    #[inline(always)]
    pub fn coords(&self, c: usize) -> [usize; 3] {
        let [nx, ny, _] = self.cells;
        [c % nx, (c / nx) % ny, c / (nx * ny)]
    }

    /// The centre coordinate of cell `i` along `axis`, snapped into the exact predicates' domain.
    #[inline]
    pub fn centre_along(&self, axis: usize, i: usize) -> f64 {
        snap_to_domain(self.origin[axis] + (i as f64 + 0.5) * self.spacing)
    }

    /// The centre of a cell.
    pub fn centre(&self, c: usize) -> [f64; 3] {
        let ijk = self.coords(c);
        [self.centre_along(0, ijk[0]), self.centre_along(1, ijk[1]), self.centre_along(2, ijk[2])]
    }

    /// The cell containing a point, if it is in the grid.
    pub fn cell_at(&self, p: [f64; 3]) -> Option<usize> {
        let mut ijk = [0usize; 3];
        for a in 0..3 {
            let t = (p[a] - self.origin[a]) / self.spacing;
            if !(t >= 0.0 && t < self.cells[a] as f64) {
                return None;
            }
            ijk[a] = (t.floor() as usize).min(self.cells[a] - 1);
        }
        Some(self.index(ijk))
    }
}

/// The sign of `orient2d(a, b, p)` for `p` perturbed by `(δ, δ²)` in its two coordinates.
/// `orient2d` is affine in `p` with gradient `(a_v − b_v, b_u − a_u)`, so when it is zero the
/// sign is that of the first nonzero gradient component; both are zero only when `a = b`.
#[inline]
fn edge_sign(a: &Point2, b: &Point2, p: &Point2) -> Sign {
    let s = orient2d(a, b, p);
    if !s.is_zero() {
        return s;
    }
    let ([au, av], [bu, bv]) = (a.coords(), b.coords());
    if av != bv { Sign::of(av - bv) } else { Sign::of(bu - au) }
}

fn point3(surface: &BoundarySurface, v: u32) -> Result<Point3> {
    surface.point3(v).ok_or_else(|| Error::Invalid(format!("surface vertex {v} is missing")))
}

/// Mark every cell of `grid` whose centre is inside `surface` (a positive winding number):
/// `inside[c]` is set to true for those. Returns how many cells were marked that were not marked
/// before.
pub fn voxelise(surface: &BoundarySurface, grid: &Grid, inside: &mut [bool]) -> Result<usize> {
    let [nx, ny, nz] = grid.cells;
    if inside.len() != grid.len() {
        return Err(Error::Invalid(format!("{} flags for {} cells", inside.len(), grid.len())));
    }
    let lines = ny * nz;
    // Winding-number changes along each line: +s at the line's first cell and −s at the first
    // centre beyond the triangle's plane.
    let mut steps = vec![0i32; lines * (nx + 1)];
    let h = grid.spacing;
    let range = |axis: usize, lo: f64, hi: f64, n: usize| -> (usize, usize) {
        // Centres with lo ≤ centre ≤ hi, widened by one cell on each side: the exact tests decide.
        let first = ((lo - grid.origin[axis]) / h - 0.5).floor() - 1.0;
        let last = ((hi - grid.origin[axis]) / h - 0.5).ceil() + 1.0;
        let clamp = |t: f64| if t.is_nan() || t < 0.0 { 0usize } else { (t as usize).min(n - 1) };
        (clamp(first), clamp(last))
    };
    for tri in surface.triangles() {
        let [a, b, c] = [point3(surface, tri[0])?, point3(surface, tri[1])?, point3(surface, tri[2])?];
        let proj = |p: &Point3| -> Result<Point2> {
            let [_, y, z] = p.coords();
            Ok(Point2::new(y, z)?)
        };
        let (a2, b2, c2) = (proj(&a)?, proj(&b)?, proj(&c)?);
        let s = orient2d(&a2, &b2, &c2);
        if s.is_zero() {
            continue;
        }
        let ([ax, ay, az], [bx, by, bz], [cx, cy, cz]) = (a.coords(), b.coords(), c.coords());
        let (j0, j1) = range(1, ay.min(by).min(cy), ay.max(by).max(cy), ny);
        let (k0, k1) = range(2, az.min(bz).min(cz), az.max(bz).max(cz), nz);
        // The plane's normal, for the floating-point estimate of the crossing only.
        let (u, v) = ([bx - ax, by - ay, bz - az], [cx - ax, cy - ay, cz - az]);
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        for k in k0..=k1 {
            let z = grid.centre_along(2, k);
            for j in j0..=j1 {
                let y = grid.centre_along(1, j);
                let p2 = Point2::new(y, z)?;
                if edge_sign(&a2, &b2, &p2) != s || edge_sign(&b2, &c2, &p2) != s || edge_sign(&c2, &a2, &p2) != s {
                    continue;
                }
                // Centres before the plane (the +x ray from them crosses the triangle):
                // orient3d(a, b, c, p) = −(p − a)·n has the sign s there. A centre on the plane is
                // perturbed by +ε in x (n_x has the sign s) and counts as beyond it.
                let behind = |i: usize| -> Result<bool> { Ok(orient3d(&a, &b, &c, &Point3::new(grid.centre_along(0, i), y, z)?) == s) };
                let x = if n[0] != 0.0 { ax - (n[1] * (y - ay) + n[2] * (z - az)) / n[0] } else { ax };
                let t = ((x - grid.origin[0]) / h - 0.5).ceil();
                let mut split = if t.is_nan() || t < 0.0 { 0 } else { (t as usize).min(nx) };
                while split > 0 && !behind(split - 1)? {
                    split -= 1;
                }
                while split < nx && behind(split)? {
                    split += 1;
                }
                let line = (j + ny * k) * (nx + 1);
                let w = s as i32;
                if let Some(first) = steps.get_mut(line) {
                    *first += w;
                }
                if let Some(end) = steps.get_mut(line + split) {
                    *end -= w;
                }
            }
        }
    }
    let mut marked = 0;
    for line in 0..lines {
        let mut w = 0i32;
        for i in 0..nx {
            w += steps.get(line * (nx + 1) + i).copied().unwrap_or(0);
            if w > 0
                && let Some(flag) = inside.get_mut(line * nx + i)
                && !*flag
            {
                *flag = true;
                marked += 1;
            }
        }
    }
    Ok(marked)
}

/// A bounding-volume tree over a surface's triangles for nearest-triangle queries.
#[derive(Debug, Clone)]
pub struct TriangleTree {
    triangles: Vec<[[f64; 3]; 3]>,
    order: Vec<u32>,
    nodes: Vec<Node>,
}

#[derive(Debug, Clone, Copy)]
struct Node {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: the range `first..first + count` of `order`; inner: `count == 0` and the children
    /// are `first` and `first + 1`.
    first: u32,
    count: u32,
}

const LEAF: usize = 4;

fn bounds(tris: &[[[f64; 3]; 3]], ids: &[u32]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &t in ids {
        if let Some(tri) = tris.get(t as usize) {
            for p in tri {
                for a in 0..3 {
                    lo[a] = lo[a].min(p[a]);
                    hi[a] = hi[a].max(p[a]);
                }
            }
        }
    }
    (lo, hi)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The squared distance from `p` to triangle `abc` (Ericson 2005, §5.1.5: the closest point by
/// Voronoi regions of the vertices, edges and face).
fn distance2(p: [f64; 3], [a, b, c]: &[[f64; 3]; 3]) -> f64 {
    let (ab, ac, ap) = (sub(*b, *a), sub(*c, *a), sub(p, *a));
    let at = |s: f64, t: f64| [a[0] + s * ab[0] + t * ac[0], a[1] + s * ab[1] + t * ac[1], a[2] + s * ab[2] + t * ac[2]];
    let d2 = |q: [f64; 3]| {
        let d = sub(p, q);
        dot(d, d)
    };
    let (d1, d2v) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2v <= 0.0 {
        return d2(*a);
    }
    let bp = sub(p, *b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return d2(*b);
    }
    let vc = d1 * d4 - d3 * d2v;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return d2(at(d1 / (d1 - d3), 0.0));
    }
    let cp = sub(p, *c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return d2(*c);
    }
    let vb = d5 * d2v - d1 * d6;
    if vb <= 0.0 && d2v >= 0.0 && d6 <= 0.0 {
        return d2(at(0.0, d2v / (d2v - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return d2([b[0] + w * (c[0] - b[0]), b[1] + w * (c[1] - b[1]), b[2] + w * (c[2] - b[2])]);
    }
    let denom = 1.0 / (va + vb + vc);
    d2(at(vb * denom, vc * denom))
}

fn box_distance2(p: [f64; 3], n: &Node) -> f64 {
    let mut d = 0.0;
    for a in 0..3 {
        let e = (n.lo[a] - p[a]).max(0.0).max(p[a] - n.hi[a]);
        d += e * e;
    }
    d
}

impl TriangleTree {
    pub fn new(surface: &BoundarySurface) -> TriangleTree {
        let pts = surface.points();
        let at = |v: u32| pts.get(v as usize).copied().unwrap_or([0.0; 3]);
        let triangles: Vec<[[f64; 3]; 3]> = surface.triangles().iter().map(|t| [at(t[0]), at(t[1]), at(t[2])]).collect();
        let centroid = |t: u32, a: usize| triangles.get(t as usize).map_or(0.0, |tri| tri[0][a] + tri[1][a] + tri[2][a]);
        let mut order: Vec<u32> = (0..triangles.len() as u32).collect();
        let mut nodes = Vec::with_capacity(2 * triangles.len() / LEAF + 1);
        let (lo, hi) = bounds(&triangles, &order);
        nodes.push(Node { lo, hi, first: 0, count: order.len() as u32 });
        // Split nodes breadth-first; each split halves a range, so this ends after fewer than
        // 2 n / LEAF nodes.
        let mut next = 0;
        while next < nodes.len() {
            let Some(node) = nodes.get(next).copied() else { break };
            next += 1;
            let (first, count) = (node.first as usize, node.count as usize);
            if count <= LEAF {
                continue;
            }
            let Some(ids) = order.get_mut(first..first + count) else { continue };
            let extent = sub(node.hi, node.lo);
            let axis = if extent[0] >= extent[1] && extent[0] >= extent[2] {
                0
            } else if extent[1] >= extent[2] {
                1
            } else {
                2
            };
            ids.sort_by(|&p, &q| centroid(p, axis).total_cmp(&centroid(q, axis)).then(p.cmp(&q)));
            let half = count / 2;
            let (left, right) = ids.split_at(half);
            let (llo, lhi) = bounds(&triangles, left);
            let (rlo, rhi) = bounds(&triangles, right);
            let child = nodes.len() as u32;
            nodes.push(Node { lo: llo, hi: lhi, first: first as u32, count: half as u32 });
            nodes.push(Node { lo: rlo, hi: rhi, first: (first + half) as u32, count: (count - half) as u32 });
            if let Some(n) = nodes.get_mut(next - 1) {
                n.first = child;
                n.count = 0;
            }
        }
        TriangleTree { triangles, order, nodes }
    }

    /// The triangle nearest `p` (the lowest index among equally near ones).
    pub fn nearest(&self, p: [f64; 3]) -> Option<u32> {
        let mut best: Option<(f64, u32)> = None;
        let better = |d: f64, t: u32, best: &Option<(f64, u32)>| best.is_none_or(|(bd, bt)| d < bd || (d == bd && t < bt));
        let mut stack = vec![0u32];
        // Every node is pushed at most once.
        while let Some(i) = stack.pop() {
            let Some(node) = self.nodes.get(i as usize) else { continue };
            if best.is_some_and(|(bd, _)| box_distance2(p, node) > bd) {
                continue;
            }
            if node.count > 0 {
                let range = node.first as usize..(node.first + node.count) as usize;
                for &t in self.order.get(range).unwrap_or_default() {
                    if let Some(tri) = self.triangles.get(t as usize) {
                        let d = distance2(p, tri);
                        if better(d, t, &best) {
                            best = Some((d, t));
                        }
                    }
                }
            } else {
                let (l, r) = (node.first, node.first + 1);
                let dist = |c: u32| self.nodes.get(c as usize).map_or(f64::INFINITY, |n| box_distance2(p, n));
                // Visit the nearer child first (pushed last).
                if dist(l) <= dist(r) {
                    stack.extend([r, l]);
                } else {
                    stack.extend([l, r]);
                }
            }
        }
        best.map(|(_, t)| t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_mesh::surface::SurfaceOptions;

    /// The surface of the box [lo, hi], outward, two triangles per side, side k tagged k.
    pub(crate) fn cube(lo: [f64; 3], hi: [f64; 3]) -> BoundarySurface {
        let p = |i: usize| [if i & 1 == 0 { lo[0] } else { hi[0] }, if i & 2 == 0 { lo[1] } else { hi[1] }, if i & 4 == 0 { lo[2] } else { hi[2] }];
        let positions: Vec<[f64; 3]> = (0..8).map(p).collect();
        // Quads counter-clockwise seen from outside.
        let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
        let mut triangles = Vec::new();
        let mut tags = Vec::new();
        for (k, q) in quads.iter().enumerate() {
            triangles.push([q[0], q[1], q[2]]);
            triangles.push([q[0], q[2], q[3]]);
            tags.extend([k as u32; 2]);
        }
        BoundarySurface::new(&positions, &triangles, &tags, &SurfaceOptions::default()).unwrap()
    }

    #[test]
    fn projected_orientation_is_the_normal_x_sign() {
        let pt = |x, y, z| Point2::new(y, z).map(|p| (x, p)).unwrap().1;
        // Triangle (0,0,0), (0,1,0), (0,0,1): normal +x; its (y, z) projection is counter-clockwise.
        assert_eq!(orient2d(&pt(0.0, 0.0, 0.0), &pt(0.0, 1.0, 0.0), &pt(0.0, 0.0, 1.0)), Sign::Positive);
    }

    #[test]
    fn a_box_is_voxelised_exactly_even_with_faces_on_centres() {
        // Cells of size 1 from 0 to 8; the box's faces pass through cell centres (x = 2.5, 5.5)
        // and between them (y = 2, 6), and its edges lie on grid lines.
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [8, 8, 8] };
        let s = cube([2.5, 2.0, 2.5], [5.5, 6.0, 5.5]);
        let mut inside = vec![false; grid.len()];
        let n = voxelise(&s, &grid, &mut inside).unwrap();
        // Centres on a face are decided by the perturbation p + (ε, δ, δ²): beyond it in +x, y and z.
        let mut count = 0;
        for c in 0..grid.len() {
            let [x, y, z] = grid.centre(c);
            let strictly = x > 2.5 && x < 5.5 && y > 2.0 && y < 6.0 && z > 2.5 && z < 5.5;
            let maybe = (2.5..=5.5).contains(&x) && (2.0..=6.0).contains(&y) && (2.5..=5.5).contains(&z);
            if strictly {
                assert!(inside[c], "{x} {y} {z}");
            }
            if !maybe {
                assert!(!inside[c], "{x} {y} {z}");
            }
            count += usize::from(inside[c]);
        }
        assert_eq!(n, count);
        // x: 2.5 (moved inside), 3.5, 4.5 (5.5 moves outside); y: 2.5 … 5.5; z as x.
        assert_eq!(n, 3 * 4 * 3, "{n}");
    }

    #[test]
    fn nearest_triangle_and_its_face() {
        let s = cube([0.0; 3], [1.0; 3]);
        let tree = TriangleTree::new(&s);
        // Just outside the +y side, near the -x edge: side 3 (y = hi).
        let t = tree.nearest([0.1, 1.2, 0.5]).unwrap();
        assert_eq!(s.tri_face()[t as usize], 3);
        let t = tree.nearest([0.5, 0.5, -0.3]).unwrap();
        assert_eq!(s.tri_face()[t as usize], 4);
        // Against brute force on scattered points.
        for i in 0..200 {
            let f = |k: u32| ((i * 7919 + k * 104_729) % 1000) as f64 / 400.0 - 0.75;
            let p = [f(1), f(2), f(3)];
            let brute = (0..s.triangles().len() as u32)
                .min_by(|&a, &b| {
                    let tri = |t: u32| {
                        let v = s.triangles()[t as usize];
                        [s.points()[v[0] as usize], s.points()[v[1] as usize], s.points()[v[2] as usize]]
                    };
                    distance2(p, &tri(a)).total_cmp(&distance2(p, &tri(b))).then(a.cmp(&b))
                })
                .unwrap();
            assert_eq!(tree.nearest(p), Some(brute), "{p:?}");
        }
    }
}
