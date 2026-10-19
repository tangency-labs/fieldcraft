//! Contact geometry shared by the solvers: the closest point of an element face to a point (the
//! local search) and a bucket grid of face boxes (the global search).
//!
//! **Closest point.** A face is a flat triangle or a bilinear quadrilateral given by its corner
//! nodes, numbered so that their right-hand normal points out of the element. On a triangle the
//! closest point is found exactly by the region tests of C. Ericson, *Real-Time Collision
//! Detection*, Morgan Kaufmann 2005, §5.1.5. On a bilinear quadrilateral x(ξ, η) it is the
//! minimum of ½|x(ξ, η) − p|² over [−1, 1]², found by Newton's method on the stationarity
//! conditions x,ξ · (x − p) = x,η · (x − p) = 0 with the exact Hessian (the bilinear patch's
//! only second derivative is x,ξη), projected onto the square after each step: the
//! closest-point projection of node-to-segment contact (J. O. Hallquist, G. L. Goudreau, D. J.
//! Benson, *Sliding interfaces with contact-impact in large-scale Lagrangian computations*,
//! CMAME 51:107–137, 1985; P. Wriggers, *Computational Contact Mechanics*, 2nd ed., Springer
//! 2006, §4.1). The signed gap is measured along the face's unit normal at the closest point; the
//! lateral offset (how far the point lies beside the face rather than above or below it) tells a
//! point over the face from one past its edge.
//!
//! **Global search.** Each face's bounding box, grown by a margin, goes into every cell of a
//! uniform grid it overlaps (a bucket sort: D. J. Benson, J. O. Hallquist, *A single surface
//! contact algorithm for the post-buckling analysis of shell structures*, CMAME 78:141–163,
//! 1990); a point's candidates are the faces in its cell. The cell is twice the mean box size,
//! grown until the grid and its entries fit their caps, so the search is linear in the number of
//! faces and its memory bounded whatever the input. Entries are listed in face order, so the
//! candidates (and everything decided from them) do not depend on how the work is split.

pub use crate::frame::Vec3;
use crate::frame::{cross, dot, norm, sub};

/// Newton iterations of the bilinear closest-point search.
const MAX_NEWTON: usize = 30;

/// The closest point of a face to a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Projection {
    /// Face coordinates of the closest point: (ξ, η) ∈ [−1, 1]² on a quadrilateral, the
    /// barycentric coordinates of the second and third nodes on a triangle.
    pub xi: [f64; 2],
    /// The face's shape functions at the closest point, per node (0 past the face's nodes).
    pub weights: [f64; 4],
    pub point: Vec3,
    /// The unit outward normal at the closest point.
    pub normal: Vec3,
    /// (p − point) · normal: negative when the point is behind the face.
    pub gap: f64,
    /// |p − point|.
    pub distance: f64,
    /// The part of p − point across the normal: 0 when the point is over or under the face, its
    /// distance from the face's edge when it is beside it.
    pub lateral: f64,
}

/// The closest point of a face (three or four corner nodes, in its node order) to `p`; `None`
/// for another number of nodes, non-finite coordinates or a face collapsed at that point.
pub fn closest_point(face: &[Vec3], p: Vec3) -> Option<Projection> {
    if face.iter().chain(std::iter::once(&p)).flatten().any(|c| !c.is_finite()) {
        return None;
    }
    let (xi, weights, normal) = match face {
        [a, b, c] => triangle(*a, *b, *c, p)?,
        [a, b, c, d] => quadrilateral([*a, *b, *c, *d], p)?,
        _ => return None,
    };
    let mut point = [0.0; 3];
    for (x, w) in face.iter().zip(weights) {
        for k in 0..3 {
            point[k] += w * x[k];
        }
    }
    let d = sub(p, point);
    let gap = dot(d, normal);
    let distance = norm(d);
    let lateral = (distance * distance - gap * gap).max(0.0).sqrt();
    Some(Projection { xi, weights, point, normal, gap, distance, lateral })
}

/// A vector shorter than this fraction of its scale has no direction.
const NO_DIRECTION: f64 = 1e-12; // unit-factors: relative to the vector's scale

/// The unit vector along `v`, unless it is too short against `scale` (a squared length).
fn unit(v: Vec3, scale: f64) -> Option<Vec3> {
    let n = norm(v);
    (n > NO_DIRECTION * scale && n.is_finite()).then(|| v.map(|c| c / n))
}

/// Ericson's closest point on triangle abc: (barycentric v, w), weights, unit normal.
fn triangle(a: Vec3, b: Vec3, c: Vec3, p: Vec3) -> Option<([f64; 2], [f64; 4], Vec3)> {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let normal = unit(cross(ab, ac), dot(ab, ab).max(dot(ac, ac)))?;
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    let (va, vb, vc) = (d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2);
    let (v, w) = if d1 <= 0.0 && d2 <= 0.0 {
        (0.0, 0.0)
    } else if d3 >= 0.0 && d4 <= d3 {
        (1.0, 0.0)
    } else if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        (d1 / (d1 - d3), 0.0)
    } else if d6 >= 0.0 && d5 <= d6 {
        (0.0, 1.0)
    } else if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        (0.0, d2 / (d2 - d6))
    } else if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        (1.0 - t, t)
    } else {
        let s = va + vb + vc;
        (vb / s, vc / s)
    };
    if !(v.is_finite() && w.is_finite()) {
        return None;
    }
    Some(([v, w], [1.0 - v - w, v, w, 0.0], normal))
}

/// The bilinear shape functions and their derivatives at (ξ, η), nodes at (−1, −1), (1, −1),
/// (1, 1), (−1, 1).
fn bilinear(s: f64, t: f64) -> ([f64; 4], [f64; 4], [f64; 4]) {
    let n = [(1.0 - s) * (1.0 - t), (1.0 + s) * (1.0 - t), (1.0 + s) * (1.0 + t), (1.0 - s) * (1.0 + t)].map(|v| 0.25 * v);
    let ds = [-(1.0 - t), 1.0 - t, 1.0 + t, -(1.0 + t)].map(|v| 0.25 * v);
    let dt = [-(1.0 - s), -(1.0 + s), 1.0 + s, 1.0 - s].map(|v| 0.25 * v);
    (n, ds, dt)
}

fn combine(x: &[Vec3; 4], w: &[f64; 4]) -> Vec3 {
    let mut out = [0.0; 3];
    for (xi, wi) in x.iter().zip(w) {
        for k in 0..3 {
            out[k] += wi * xi[k];
        }
    }
    out
}

/// The closest point of a bilinear quadrilateral by projected Newton iterations.
fn quadrilateral(x: [Vec3; 4], p: Vec3) -> Option<([f64; 2], [f64; 4], Vec3)> {
    let size = dot(sub(x[2], x[0]), sub(x[2], x[0])).max(dot(sub(x[3], x[1]), sub(x[3], x[1])));
    if !(size > 0.0) {
        return None;
    }
    // x,ξη: constant over the patch.
    let twist = [0, 1, 2].map(|k| 0.25 * (x[0][k] - x[1][k] + x[2][k] - x[3][k]));
    let (mut s, mut t) = (0.0f64, 0.0f64);
    for _ in 0..MAX_NEWTON {
        let (n, ds, dt) = bilinear(s, t);
        let r = sub(combine(&x, &n), p);
        let (a, b) = (combine(&x, &ds), combine(&x, &dt));
        let g = [dot(a, r), dot(b, r)];
        let (haa, hbb) = (dot(a, a), dot(b, b));
        let mut hab = dot(a, b) + dot(r, twist);
        // Far from the patch the exact Hessian may be indefinite: Gauss–Newton then.
        if haa * hbb - hab * hab <= 1e-14 * haa * hbb {
            hab = dot(a, b);
        }
        let det = haa * hbb - hab * hab;
        if !(det > 0.0) {
            return None;
        }
        let step = [-(hbb * g[0] - hab * g[1]) / det, -(haa * g[1] - hab * g[0]) / det];
        let (s1, t1) = ((s + step[0]).clamp(-1.0, 1.0), (t + step[1]).clamp(-1.0, 1.0));
        let moved = (s1 - s).abs().max((t1 - t).abs());
        (s, t) = (s1, t1);
        if !moved.is_finite() {
            return None;
        }
        if moved < 1e-13 {
            break;
        }
    }
    let (n, ds, dt) = bilinear(s, t);
    let normal = unit(cross(combine(&x, &ds), combine(&x, &dt)), size)?;
    Some(([s, t], n, normal))
}

/// A box: its lowest and highest corners.
pub type Aabb = [Vec3; 2];

/// The bounding box of a face's nodes, grown by `margin` on every side.
pub fn bounds(face: &[Vec3], margin: f64) -> Aabb {
    let mut b = [[f64::INFINITY; 3], [f64::NEG_INFINITY; 3]];
    for x in face {
        for k in 0..3 {
            b[0][k] = b[0][k].min(x[k] - margin);
            b[1][k] = b[1][k].max(x[k] + margin);
        }
    }
    b
}

/// A uniform grid of boxes (a bucket sort), for the global contact search.
#[derive(Debug, Clone, Default)]
pub struct Grid {
    lowest: Vec3,
    cell: f64,
    dims: [usize; 3],
    /// Per cell, the range of its entries in `items` (compressed rows).
    start: Vec<usize>,
    items: Vec<u32>,
}

/// Most cells of a grid (whatever `max_cells` says).
pub const MAX_GRID_CELLS: usize = 1 << 24;

/// Most entries of a grid (box × cell overlaps).
pub const MAX_GRID_ITEMS: usize = 1 << 27;

/// The smallest grid cell, as a fraction of the boxes' span (point-like boxes).
const MIN_CELL: f64 = 1e-6; // unit-factors: relative to the boxes' span

impl Grid {
    /// The grid of `boxes` (their indices are the entries; boxes that are not finite are left
    /// out), with at most `max_cells` cells (and [`MAX_GRID_CELLS`]).
    pub fn new(boxes: &[Aabb], max_cells: usize) -> Grid {
        let finite = |b: &Aabb| b.iter().flatten().all(|c| c.is_finite()) && (0..3).all(|k| b[0][k] <= b[1][k]);
        let mut lowest = [f64::INFINITY; 3];
        let mut highest = [f64::NEG_INFINITY; 3];
        let (mut count, mut extent) = (0usize, 0.0);
        for b in boxes.iter().filter(|b| finite(b)) {
            for k in 0..3 {
                lowest[k] = lowest[k].min(b[0][k]);
                highest[k] = highest[k].max(b[1][k]);
            }
            extent += (0..3).map(|k| b[1][k] - b[0][k]).fold(0.0, f64::max);
            count += 1;
        }
        if count == 0 {
            return Grid::default();
        }
        let span = (0..3).map(|k| highest[k] - lowest[k]).fold(0.0, f64::max);
        let max_cells = max_cells.clamp(1, MAX_GRID_CELLS);
        // Twice the mean box size, but not below a millionth of the span (point-like boxes).
        let mut cell = (2.0 * extent / count as f64).max(MIN_CELL * span).max(f64::MIN_POSITIVE);
        for _ in 0..200 {
            let dims = [0, 1, 2].map(|k| ((highest[k] - lowest[k]) / cell).floor() + 1.0);
            let cells = dims[0] * dims[1] * dims[2];
            if !(cells <= max_cells as f64) {
                cell *= (cells / max_cells as f64).cbrt().max(1.01) * 1.01;
                continue;
            }
            let dims = dims.map(|d| d as usize);
            let mut grid = Grid { lowest, cell, dims, start: vec![0; dims[0] * dims[1] * dims[2] + 1], items: Vec::new() };
            // Count the entries per cell, then fill (compressed rows, boxes in order).
            let mut total = 0usize;
            for b in boxes.iter().filter(|b| finite(b)) {
                let (lo, hi) = (grid.coordinates(b[0]), grid.coordinates(b[1]));
                total = total.saturating_add((0..3).map(|k| hi[k] - lo[k] + 1).product::<usize>());
            }
            if total > MAX_GRID_ITEMS {
                cell *= 2.0;
                continue;
            }
            for b in boxes.iter().filter(|b| finite(b)) {
                grid.each_cell(b, |c, start| start[c + 1] += 1);
            }
            for c in 0..grid.start.len() - 1 {
                grid.start[c + 1] += grid.start[c];
            }
            let mut next = grid.start.clone();
            let mut items = vec![0u32; total];
            for (i, b) in boxes.iter().enumerate().filter(|(_, b)| finite(b)) {
                grid.each_cell(b, |c, _| {
                    if let (Some(slot), Some(item)) = (next.get_mut(c), u32::try_from(i).ok()) {
                        if let Some(e) = items.get_mut(*slot) {
                            *e = item;
                        }
                        *slot += 1;
                    }
                });
            }
            grid.items = items;
            return grid;
        }
        Grid::default()
    }

    /// The cell coordinates of a point inside the grid's span (clamped to it).
    fn coordinates(&self, x: Vec3) -> [usize; 3] {
        [0, 1, 2].map(|k| {
            let c = ((x[k] - self.lowest[k]) / self.cell).floor();
            if c > 0.0 { (c as usize).min(self.dims[k].saturating_sub(1)) } else { 0 }
        })
    }

    fn index(&self, c: [usize; 3]) -> usize {
        c[0] + self.dims[0] * (c[1] + self.dims[1] * c[2])
    }

    /// Call `f(cell, start)` for every cell a box overlaps.
    fn each_cell(&mut self, b: &Aabb, mut f: impl FnMut(usize, &mut Vec<usize>)) {
        let (lo, hi) = (self.coordinates(b[0]), self.coordinates(b[1]));
        for k in lo[2]..=hi[2] {
            for j in lo[1]..=hi[1] {
                for i in lo[0]..=hi[0] {
                    let c = self.index([i, j, k]);
                    f(c, &mut self.start);
                }
            }
        }
    }

    /// The boxes that may hold `p`: those listed in its cell, in box order (empty outside the
    /// grid).
    pub fn candidates(&self, p: Vec3) -> &[u32] {
        if self.items.is_empty() || p.iter().any(|c| !c.is_finite()) {
            return &[];
        }
        let outside = (0..3).any(|k| p[k] < self.lowest[k] || p[k] > self.lowest[k] + self.cell * self.dims[k] as f64);
        if outside {
            return &[];
        }
        let c = self.index(self.coordinates(p));
        match (self.start.get(c), self.start.get(c + 1)) {
            (Some(a), Some(b)) => self.items.get(*a..*b).unwrap_or_default(),
            _ => &[],
        }
    }

    /// The boxes that may overlap box `b`: those listed in the cells it overlaps, in box order
    /// without repeats (into `out`, cleared first).
    pub fn candidates_in(&self, b: &Aabb, out: &mut Vec<u32>) {
        out.clear();
        if self.items.is_empty() || b.iter().flatten().any(|c| !c.is_finite()) {
            return;
        }
        let top = [0, 1, 2].map(|k| self.lowest[k] + self.cell * self.dims[k] as f64);
        if (0..3).any(|k| b[1][k] < self.lowest[k] || b[0][k] > top[k]) {
            return;
        }
        let (lo, hi) = (self.coordinates(b[0]), self.coordinates(b[1]));
        for k in lo[2]..=hi[2] {
            for j in lo[1]..=hi[1] {
                for i in lo[0]..=hi[0] {
                    let c = self.index([i, j, k]);
                    if let (Some(a), Some(e)) = (self.start.get(c), self.start.get(c + 1)) {
                        out.extend_from_slice(self.items.get(*a..*e).unwrap_or_default());
                    }
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }

    /// The number of entries (box × cell overlaps).
    pub fn entries(&self) -> usize {
        self.items.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn triangle_projection_covers_the_face_its_edges_and_corners() {
        let t = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        // Above the inside: the foot of the perpendicular, positive gap along +z.
        let p = closest_point(&t, [0.5, 0.25, 0.3]).unwrap();
        assert!(close(p.gap, 0.3) && close(p.lateral, 0.0) && close(p.point[0], 0.5) && close(p.point[1], 0.25));
        assert!(close(p.weights.iter().sum::<f64>(), 1.0) && p.normal == [0.0, 0.0, 1.0]);
        // Below: a negative gap.
        assert!(close(closest_point(&t, [0.5, 0.5, -0.1]).unwrap().gap, -0.1));
        // Past the hypotenuse: the closest point is on it, the offset lateral.
        let p = closest_point(&t, [2.0, 2.0, 0.0]).unwrap();
        assert!(close(p.point[0], 1.0) && close(p.point[1], 1.0) && close(p.lateral, 2f64.sqrt()) && close(p.gap, 0.0));
        // Past a corner.
        let p = closest_point(&t, [-1.0, -1.0, 1.0]).unwrap();
        assert_eq!(p.weights, [1.0, 0.0, 0.0, 0.0]);
        assert!(close(p.distance, 3f64.sqrt()));
        // Degenerate and hostile faces give nothing.
        assert!(closest_point(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [0.0, 1.0, 0.0]).is_none());
        assert!(closest_point(&t, [f64::NAN, 0.0, 0.0]).is_none());
        assert!(closest_point(&t[..2], [0.0; 3]).is_none());
    }

    #[test]
    fn quadrilateral_projection_finds_the_foot_on_a_warped_patch() {
        let flat = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let p = closest_point(&flat, [1.5, 0.25, -0.2]).unwrap();
        assert!(close(p.xi[0], 0.5) && close(p.xi[1], -0.5) && close(p.gap, -0.2) && close(p.lateral, 0.0));
        assert!(close(p.weights.iter().sum::<f64>(), 1.0));
        // A warped (hyperbolic paraboloid) patch: the residual is normal to the surface.
        let warped = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2], [1.0, 1.0, 0.0], [0.0, 1.0, 0.2]];
        let q = [0.3, 0.6, 0.5];
        let p = closest_point(&warped, q).unwrap();
        let (_, ds, dt) = bilinear(p.xi[0], p.xi[1]);
        let r = sub(q, p.point);
        assert!(dot(r, combine(&warped, &ds)).abs() < 1e-12 && dot(r, combine(&warped, &dt)).abs() < 1e-12, "{p:?}");
        assert!(p.lateral < 1e-9 && p.gap > 0.0);
        // Beyond an edge: clamped, with the offset lateral.
        let p = closest_point(&flat, [3.0, 0.5, 0.0]).unwrap();
        assert!(close(p.xi[0], 1.0) && close(p.lateral, 1.0));
        // A quadrilateral collapsed to a segment has no normal.
        assert!(closest_point(&[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0; 3]], [0.5, 0.5, 0.0]).is_none());
    }

    #[test]
    fn the_grid_finds_every_box_holding_a_point_and_stays_bounded() {
        let boxes: Vec<Aabb> = (0..50).map(|i| [[i as f64, 0.0, 0.0], [i as f64 + 1.5, 1.0, 1.0]]).collect();
        let grid = Grid::new(&boxes, 1000);
        for k in 0..500 {
            let p = [k as f64 * 0.1 + 0.01, 0.5, 0.5];
            let found = grid.candidates(p);
            for (i, b) in boxes.iter().enumerate() {
                let inside = (0..3).all(|d| p[d] >= b[0][d] && p[d] <= b[1][d]);
                assert!(!inside || found.contains(&(i as u32)), "point {p:?} misses box {i}");
            }
            assert!(found.windows(2).all(|w| w[0] < w[1]), "candidates in box order");
        }
        assert!(grid.candidates([-5.0, 0.0, 0.0]).is_empty() && grid.candidates([f64::NAN; 3]).is_empty());
        // One huge box among small ones, a cap of one cell, boxes that are not finite: still
        // bounded and correct.
        let mut mixed = boxes.clone();
        mixed.push([[-1e6, -1e6, -1e6], [1e6, 1e6, 1e6]]);
        mixed.push([[f64::NAN, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        for cells in [1, 10, 1 << 30] {
            let g = Grid::new(&mixed, cells);
            assert!(g.candidates([10.2, 0.5, 0.5]).contains(&50) && g.candidates([10.2, 0.5, 0.5]).contains(&9));
            assert!(!g.candidates([10.2, 0.5, 0.5]).contains(&51));
            assert!(g.entries() <= MAX_GRID_ITEMS);
        }
        assert_eq!(Grid::new(&[], 10).candidates([0.0; 3]), &[] as &[u32]);
        // Box queries: every box that overlaps, once, in order.
        let mut found = Vec::new();
        grid.candidates_in(&[[10.2, 0.2, 0.2], [12.1, 0.4, 0.4]], &mut found);
        let expected: Vec<u32> = (0..50).filter(|i| boxes[*i as usize][0][0] <= 12.1 && boxes[*i as usize][1][0] >= 10.2).collect();
        assert!(expected.iter().all(|e| found.contains(e)) && found.windows(2).all(|w| w[0] < w[1]), "{found:?} {expected:?}");
        grid.candidates_in(&[[-9.0, 0.0, 0.0], [-8.0, 1.0, 1.0]], &mut found);
        assert!(found.is_empty());
    }
}
