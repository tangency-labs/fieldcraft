//! Display data of a flow field on a regular grid of cells (a lattice Boltzmann result): a slice
//! plane through the cells carrying their values, a field sampled at a surface beside the cells,
//! and streamlines traced through the velocity field.
//!
//! Values are cell-centred; a cell with no value (inside an obstacle) holds NaN. A slice joins the
//! cell centres of the layers either side of its plane (interpolated linearly between them) into
//! triangles carrying their corners' values, so [`crate::contour::split_triangle`] bands them as it
//! bands a nodal field; the outermost centres are moved out to the domain's faces, so the slice
//! covers the whole domain. A triangle touching cells with no value takes the mean of its other
//! corners there, and one with none is left out (the obstacle is drawn there).
//!
//! Streamlines follow the trilinear interpolation of the cell velocities between the cell centres
//! (solid cells at rest: the no-slip value), integrated along their arc length by the midpoint
//! (second-order Runge–Kutta) method with a step of half a cell. Each stops where it leaves the
//! domain, enters a solid cell or reaches still fluid, or at its step cap; all of them together
//! stop at a cap on points, so hostile counts cost a bounded time.

use crate::{Aabb3, Vec3};

/// Most quadrilaterals a slice is made of: a larger layer is sampled every few cells.
pub const MAX_SLICE_QUADS: usize = 1 << 20;
/// Most streamlines traced at once.
pub const MAX_STREAMLINES: usize = 2_000;
/// Most steps of one streamline.
pub const MAX_STREAMLINE_STEPS: usize = 100_000;
/// Most points of all the streamlines traced at once.
pub const MAX_STREAMLINE_POINTS: usize = 2_000_000;
/// A streamline's step, in cells.
const STEP: f64 = 0.5;
/// A speed below this fraction of the field's largest is still fluid, where a streamline ends.
const STILL: f64 = 1e-6; // unit-factors: relative to the field's largest speed
/// Interpolation weights summing to less than this leave no value to interpolate.
const NO_WEIGHT: f64 = 1e-9; // unit-factors: a sum of interpolation weights

/// A regular grid of cubic cells: cell `(i, j, k)` spans `origin + spacing · (i, j, k)` to one
/// spacing further, and its values sit at index `i + nx (j + ny k)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellGrid {
    pub origin: [f64; 3],
    pub spacing: f64,
    pub cells: [usize; 3],
}

impl CellGrid {
    /// The number of cells, or `None` for a grid that cannot hold values: no cells, a spacing
    /// that is not a positive number, an origin that is not finite, a count that overflows.
    pub fn cell_count(&self) -> Option<usize> {
        let ok = self.spacing.is_finite() && self.spacing > 0.0 && self.origin.iter().all(|c| c.is_finite()) && self.cells.iter().all(|&n| n > 0);
        if !ok {
            return None;
        }
        self.cells[0].checked_mul(self.cells[1])?.checked_mul(self.cells[2])
    }

    fn index(&self, [i, j, k]: [usize; 3]) -> usize {
        i + self.cells[0] * (j + self.cells[1] * k)
    }

    /// The `(i, j, k)` of cell index `c`.
    pub fn coords(&self, c: usize) -> [usize; 3] {
        let [nx, ny, _] = self.cells.map(|n| n.max(1));
        [c % nx, (c / nx) % ny, c / (nx * ny)]
    }

    /// The centre of cell `(i, j, k)`.
    pub fn centre(&self, ijk: [usize; 3]) -> Vec3 {
        Vec3::new(self.centre_along(0, ijk[0]), self.centre_along(1, ijk[1]), self.centre_along(2, ijk[2]))
    }

    fn centre_along(&self, axis: usize, i: usize) -> f64 {
        self.start(axis) + (i as f64 + 0.5) * self.spacing
    }

    /// The near end of the grid along `axis` (its origin there).
    pub fn start(&self, axis: usize) -> f64 {
        self.origin.get(axis).copied().unwrap_or(0.0)
    }

    /// The far end of the grid along `axis`.
    pub fn end(&self, axis: usize) -> f64 {
        self.start(axis) + self.cells.get(axis).copied().unwrap_or(0) as f64 * self.spacing
    }

    /// The box the cells fill.
    pub fn bounds(&self) -> Aabb3 {
        let mut b = Aabb3::EMPTY;
        b.add(Vec3::from(self.origin));
        b.add(Vec3::new(self.end(0), self.end(1), self.end(2)));
        b
    }

    /// The cell holding `p`, if it is in the grid.
    pub fn cell_at(&self, p: Vec3) -> Option<usize> {
        let mut ijk = [0usize; 3];
        for (axis, slot) in ijk.iter_mut().enumerate() {
            let n = self.cells.get(axis).copied().unwrap_or(0);
            let t = (coord(p, axis) - self.start(axis)) / self.spacing;
            if !(t >= 0.0 && t < n as f64) {
                return None;
            }
            *slot = (t.floor() as usize).min(n.saturating_sub(1));
        }
        Some(self.index(ijk))
    }

    /// The cells around `p` for trilinear interpolation between their centres: per axis the lower
    /// index, the upper one and the upper one's weight. `None` when `p` lies more than `margin`
    /// cells outside the grid.
    fn stencil(&self, p: Vec3, margin: f64) -> Option<[(usize, usize, f64); 3]> {
        let mut out = [(0, 0, 0.0); 3];
        for (axis, slot) in out.iter_mut().enumerate() {
            let n = self.cells.get(axis).copied().unwrap_or(0);
            let t = (coord(p, axis) - self.start(axis)) / self.spacing - 0.5;
            if !(t >= -0.5 - margin && t <= n as f64 - 0.5 + margin) {
                return None;
            }
            let t = t.clamp(0.0, n.saturating_sub(1) as f64);
            let lo = t.floor() as usize;
            let hi = (lo + 1).min(n.saturating_sub(1));
            *slot = (lo, hi, t - lo as f64);
        }
        Some(out)
    }
}

/// The eight cells of a trilinear stencil ([`CellGrid::stencil`]) with their weights.
fn corners(st: &[(usize, usize, f64); 3]) -> impl Iterator<Item = ([usize; 3], f64)> + '_ {
    (0..8).map(move |corner| {
        let mut ijk = [0; 3];
        let mut w = 1.0;
        for (axis, (lo, hi, f)) in st.iter().enumerate() {
            let up = (corner >> axis) & 1 == 1;
            ijk[axis] = if up { *hi } else { *lo };
            w *= if up { *f } else { 1.0 - f };
        }
        (ijk, w)
    })
}

/// The coordinate of `p` along `axis` (0, 1, 2: x, y, z).
fn coord(p: Vec3, axis: usize) -> f64 {
    match axis {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

/// A triangle carrying a value at each corner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValueTri {
    pub p: [Vec3; 3],
    pub v: [f64; 3],
}

/// A slice plane through a grid: its triangles with their corners' values, and its normal.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Slice {
    pub normal: Vec3,
    pub tris: Vec<ValueTri>,
}

/// The cell indices a slice samples along an axis of `n` cells, every `stride`-th and the last,
/// each with the coordinate its vertex is drawn at: the cell centre, the first and last moved out
/// to the grid's faces (`lo`, `hi`). One cell gives two vertices, one at each face.
fn samples(n: usize, stride: usize, lo: f64, hi: f64, centre: impl Fn(usize) -> f64) -> Vec<(usize, f64)> {
    let mut out: Vec<(usize, f64)> = (0..n).step_by(stride.max(1)).map(|i| (i, centre(i))).collect();
    if let Some(last) = n.checked_sub(1)
        && out.last().is_some_and(|(i, _)| *i != last)
    {
        out.push((last, centre(last)));
    }
    if out.len() == 1
        && let Some(&(i, _)) = out.first()
    {
        out.push((i, hi));
    }
    if let Some(first) = out.first_mut() {
        first.1 = lo;
    }
    if let Some(last) = out.last_mut() {
        last.1 = hi;
    }
    out
}

/// The slice of `values` (one per cell, NaN where there is none) on the plane normal to `axis`
/// (0, 1, 2: x, y, z) at coordinate `at`. `None` when the plane misses the grid, the axis is not
/// one, or the values do not match the grid.
pub fn slice(grid: &CellGrid, values: &[f64], axis: usize, at: f64) -> Option<Slice> {
    let n = grid.cell_count()?;
    if values.len() != n || axis > 2 || !at.is_finite() {
        return None;
    }
    let lo = grid.start(axis);
    if !(at >= lo && at <= grid.end(axis)) {
        return None;
    }
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    let (nw, nu, nv) = (grid.cells.get(axis).copied()?, grid.cells.get(u).copied()?, grid.cells.get(v).copied()?);
    // The layers either side of the plane, and the upper one's weight.
    let t = ((at - lo) / grid.spacing - 0.5).clamp(0.0, nw.saturating_sub(1) as f64);
    let k0 = t.floor() as usize;
    let k1 = (k0 + 1).min(nw.saturating_sub(1));
    let w = t - k0 as f64;
    let value = |a: usize, b: usize| -> f64 {
        let layer = |k: usize| {
            let mut ijk = [0; 3];
            ijk[axis] = k;
            ijk[u] = a;
            ijk[v] = b;
            values.get(grid.index(ijk)).copied().filter(|x| x.is_finite()).unwrap_or(f64::NAN)
        };
        let (x0, x1) = (layer(k0), layer(k1));
        match (x0.is_nan(), x1.is_nan()) {
            (false, false) => x0 + (x1 - x0) * w,
            // Beside an obstacle: the nearer layer's value, if it has one.
            (false, true) if w < 0.5 => x0,
            (true, false) if w >= 0.5 => x1,
            _ => f64::NAN,
        }
    };
    let stride = (nu.saturating_mul(nv) as f64 / MAX_SLICE_QUADS as f64).sqrt().ceil().max(1.0) as usize;
    let su = samples(nu, stride, grid.start(u), grid.end(u), |i| grid.centre_along(u, i));
    let sv = samples(nv, stride, grid.start(v), grid.end(v), |i| grid.centre_along(v, i));
    let grid_values: Vec<Vec<f64>> = su.iter().map(|&(a, _)| sv.iter().map(|&(b, _)| value(a, b)).collect()).collect();
    let at_vertex = |a: usize, b: usize| -> Option<(Vec3, f64)> {
        let (cu, cv) = (su.get(a)?.1, sv.get(b)?.1);
        let mut c = [0.0; 3];
        c[axis] = at;
        c[u] = cu;
        c[v] = cv;
        Some((Vec3::from(c), *grid_values.get(a)?.get(b)?))
    };
    let mut tris = Vec::new();
    for a in 0..su.len().saturating_sub(1) {
        for b in 0..sv.len().saturating_sub(1) {
            let (Some(c00), Some(c10), Some(c11), Some(c01)) = (at_vertex(a, b), at_vertex(a + 1, b), at_vertex(a + 1, b + 1), at_vertex(a, b + 1))
            else {
                continue;
            };
            let corners = [c00, c10, c11, c01];
            let (sum, count) = corners.iter().filter(|c| !c.1.is_nan()).fold((0.0, 0usize), |(s, n), c| (s + c.1, n + 1));
            if count == 0 {
                continue;
            }
            let mean = sum / count as f64;
            let [p0, p1, p2, p3] = corners.map(|c| c.0);
            let [v0, v1, v2, v3] = corners.map(|c| if c.1.is_nan() { mean } else { c.1 });
            tris.push(ValueTri { p: [p0, p1, p2], v: [v0, v1, v2] });
            tris.push(ValueTri { p: [p0, p2, p3], v: [v0, v2, v3] });
        }
    }
    let mut normal = [0.0; 3];
    normal[axis] = 1.0;
    Some(Slice { normal: Vec3::from(normal), tris })
}

/// The value of `values` (one per cell, NaN where there is none) at `p`, interpolated
/// trilinearly between the cell centres around it from those with a value; with none of them, the
/// nearest cell with a value within two cells. NaN outside the grid (more than a cell beyond it),
/// with nothing near, or when the values do not match the grid. For colouring an obstacle's
/// surface by the fluid beside it.
pub fn sample(grid: &CellGrid, values: &[f64], p: Vec3) -> f64 {
    if grid.cell_count() != Some(values.len()) {
        return f64::NAN;
    }
    sample_with(grid, |c| values.get(c).copied(), p)
}

/// [`sample`] reading each cell's value through `value` (`None` or not finite where there is
/// none), so a caller samples a few points without building a value for every cell.
pub fn sample_with(grid: &CellGrid, value: impl Fn(usize) -> Option<f64>, p: Vec3) -> f64 {
    let Some(st) = grid.stencil(p, 1.0) else { return f64::NAN };
    let at = |ijk: [usize; 3]| value(grid.index(ijk)).filter(|x| x.is_finite());
    let (mut sum, mut weight) = (0.0, 0.0);
    for (ijk, w) in corners(&st) {
        if let Some(x) = at(ijk).filter(|_| w > 0.0) {
            sum += w * x;
            weight += w;
        }
    }
    if weight > NO_WEIGHT {
        return sum / weight;
    }
    // The nearest cell with a value, two cells around the nearest centre at most.
    let near = st.map(|(lo, hi, f)| if f < 0.5 { lo } else { hi });
    let span = |axis: usize| {
        let n = near.get(axis).copied().unwrap_or(0);
        n.saturating_sub(2)..=(n + 2).min(grid.cells.get(axis).copied().unwrap_or(1).saturating_sub(1))
    };
    let mut best: Option<(f64, f64)> = None;
    for k in span(2) {
        for j in span(1) {
            for i in span(0) {
                let Some(x) = at([i, j, k]) else { continue };
                let d = grid.centre([i, j, k]).dist(p);
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, x));
                }
            }
        }
    }
    best.map_or(f64::NAN, |(_, x)| x)
}

/// A traced streamline: its points from the seed downstream and the speed at each.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Streamline {
    pub points: Vec<Vec3>,
    pub speed: Vec<f64>,
}

/// The velocity at `p` interpolated trilinearly between the cell centres (solid cells at rest);
/// `None` outside the grid.
fn velocity_at(grid: &CellGrid, velocity: &[[f32; 3]], solid: &[bool], p: Vec3) -> Option<Vec3> {
    let st = grid.stencil(p, 0.0)?;
    let mut u = Vec3::ZERO;
    for (ijk, w) in corners(&st) {
        let c = grid.index(ijk);
        if w == 0.0 || solid.get(c).copied().unwrap_or(false) {
            continue;
        }
        let Some(v) = velocity.get(c) else { continue };
        let v = Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2]));
        if v.is_finite() {
            u = u + v * w;
        }
    }
    Some(u)
}

/// Streamlines of `velocity` (one per cell; `solid` flags the cells inside obstacles, or is empty)
/// from each seed (at most [`MAX_STREAMLINES`]) downstream, each at most `max_steps` steps of half
/// a cell (capped at [`MAX_STREAMLINE_STEPS`]), all together at most [`MAX_STREAMLINE_POINTS`]
/// points. A seed outside the fluid gives no streamline; nothing is traced when the fields do not
/// match the grid.
pub fn streamlines(grid: &CellGrid, velocity: &[[f32; 3]], solid: &[bool], seeds: &[Vec3], max_steps: usize) -> Vec<Streamline> {
    let Some(n) = grid.cell_count() else { return Vec::new() };
    if velocity.len() != n || !(solid.is_empty() || solid.len() == n) {
        return Vec::new();
    }
    let fastest =
        velocity.iter().map(|v| Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2])).len()).filter(|s| s.is_finite()).fold(0.0, f64::max);
    if !(fastest > 0.0) {
        return Vec::new();
    }
    let still = STILL * fastest;
    let h = STEP * grid.spacing;
    // The flow's direction and speed at a point in the fluid.
    let dir = |p: Vec3| -> Option<(Vec3, f64)> {
        let cell = grid.cell_at(p)?;
        if solid.get(cell).copied().unwrap_or(false) {
            return None;
        }
        let u = velocity_at(grid, velocity, solid, p)?;
        let s = u.len();
        (s > still).then(|| (u * (1.0 / s), s))
    };
    let max_steps = max_steps.min(MAX_STREAMLINE_STEPS);
    let mut budget = MAX_STREAMLINE_POINTS;
    let mut out = Vec::new();
    for &seed in seeds.iter().take(MAX_STREAMLINES) {
        if budget == 0 {
            break;
        }
        let Some((mut d, s)) = dir(seed) else { continue };
        let mut line = Streamline { points: vec![seed], speed: vec![s] };
        budget -= 1;
        let mut p = seed;
        for _ in 0..max_steps {
            if budget == 0 {
                break;
            }
            let Some((mid, _)) = dir(p + d * (0.5 * h)) else { break };
            let q = p + mid * h;
            let Some((next, s)) = dir(q) else { break };
            line.points.push(q);
            line.speed.push(s);
            budget -= 1;
            (p, d) = (q, next);
        }
        if line.points.len() >= 2 {
            out.push(line);
        }
    }
    out
}

/// `count` seeds (at most [`MAX_STREAMLINES`]) on a lattice across the grid's face at the low
/// end of `axis` (`high`: the high end), half a cell inside it: over the solid cells seen along the
/// axis, the box around them grown by half its size on every side and kept on the face, so the
/// streamlines pass over and around the obstacles; over the whole face without solid cells.
/// Empty for an axis that is not one or a grid that holds no cells.
pub fn face_seeds(grid: &CellGrid, solid: &[bool], axis: usize, high: bool, count: usize) -> Vec<Vec3> {
    let Some(n) = grid.cell_count() else { return Vec::new() };
    let count = count.min(MAX_STREAMLINES);
    if axis > 2 || count == 0 {
        return Vec::new();
    }
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    let face = |a: usize| (grid.start(a), grid.end(a));
    let (mut ru, mut rv) = (face(u), face(v));
    // The solid cells' extent across the face, by cell index.
    let mut box_ = [usize::MAX, 0, usize::MAX, 0];
    for (c, _) in solid.iter().enumerate().take(n).filter(|(_, s)| **s) {
        let ijk = grid.coords(c);
        let (a, b) = (ijk[u], ijk[v]);
        box_ = [box_[0].min(a), box_[1].max(a), box_[2].min(b), box_[3].max(b)];
    }
    if box_[0] <= box_[1] {
        let grow = |(lo, hi): (f64, f64), first: usize, last: usize, axis: usize| {
            let a = grid.start(axis) + first as f64 * grid.spacing;
            let b = grid.start(axis) + (last + 1) as f64 * grid.spacing;
            let m = 0.5 * (b - a);
            ((a - m).max(lo), (b + m).min(hi))
        };
        ru = grow(ru, box_[0], box_[1], u);
        rv = grow(rv, box_[2], box_[3], v);
    }
    let (wu, wv) = (ru.1 - ru.0, rv.1 - rv.0);
    if !(wu > 0.0 && wv > 0.0) {
        return Vec::new();
    }
    let cu = ((count as f64 * wu / wv).sqrt().round() as usize).clamp(1, count);
    let cv = (count / cu).max(1);
    let at = if high { grid.end(axis) - 0.5 * grid.spacing } else { grid.start(axis) + 0.5 * grid.spacing };
    let mut out = Vec::with_capacity(cu * cv);
    for a in 0..cu {
        for b in 0..cv {
            let mut c = [0.0; 3];
            c[axis] = at;
            c[u] = ru.0 + (a as f64 + 0.5) / cu as f64 * wu;
            c[v] = rv.0 + (b as f64 + 0.5) / cv as f64 * wv;
            out.push(Vec3::from(c));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(n: [usize; 3]) -> CellGrid {
        CellGrid { origin: [0.0; 3], spacing: 1.0, cells: n }
    }

    /// A field linear in x (the value is the cell centre's x) slices to a plane with exact values
    /// at its corners, covering the whole domain.
    #[test]
    fn a_slice_covers_the_domain_and_interpolates() {
        let g = grid([4, 3, 2]);
        let values: Vec<f64> = (0..g.cell_count().unwrap()).map(|c| (c % 4) as f64 + 0.5).collect();
        let s = slice(&g, &values, 2, 1.0).unwrap();
        assert_eq!(s.normal, Vec3::Z);
        assert_eq!(s.tris.len(), 2 * 3 * 2, "4 × 3 samples, 3 × 2 quads");
        let mut b = Aabb3::EMPTY;
        for t in &s.tris {
            for (p, v) in t.p.iter().zip(t.v) {
                b.add(*p);
                assert_eq!(p.z, 1.0);
                // Inside, the value is x; at the faces, the outer cells' values.
                let expect = p.x.clamp(0.5, 3.5);
                assert!((v - expect).abs() < 1e-12, "{p:?} {v}");
            }
        }
        assert_eq!((b.min, b.max), (Vec3::new(0.0, 0.0, 1.0), Vec3::new(4.0, 3.0, 1.0)));
        // Between two layers the values are interpolated: z-dependent field, plane at z = 1.
        let vz: Vec<f64> = (0..g.cell_count().unwrap()).map(|c| (c / 12) as f64 * 10.0).collect();
        let s = slice(&g, &vz, 2, 1.0).unwrap();
        assert!(s.tris.iter().all(|t| t.v.iter().all(|v| (v - 5.0).abs() < 1e-12)));
        // Along x and y too; one cell along an axis still makes a quad across it.
        assert!(slice(&g, &values, 0, 2.0).is_some_and(|s| s.normal == Vec3::X && s.tris.len() == 2 * 2));
        assert!(slice(&grid([1, 1, 1]), &[3.0], 1, 0.5).is_some_and(|s| s.tris.len() == 2));
    }

    /// Cells with no value: a quad touching some takes its other corners' mean, one with none is
    /// left out; hostile planes and fields are refused.
    #[test]
    fn a_slice_skips_solid_cells_and_refuses_bad_input() {
        let g = grid([3, 3, 1]);
        let mut values = vec![1.0; 9];
        values[4] = f64::NAN;
        let s = slice(&g, &values, 2, 0.5).unwrap();
        assert_eq!(s.tris.len(), 8);
        assert!(s.tris.iter().all(|t| t.v == [1.0; 3]));
        let all_solid = slice(&g, &[f64::NAN; 9], 2, 0.5).unwrap();
        assert!(all_solid.tris.is_empty());
        for (axis, at) in [(3, 0.5), (2, -0.1), (2, 1.1), (2, f64::NAN), (2, f64::INFINITY)] {
            assert!(slice(&g, &values, axis, at).is_none(), "{axis} {at}");
        }
        assert!(slice(&g, &values[..8], 2, 0.5).is_none());
        let bad = [
            CellGrid { spacing: 0.0, ..g },
            CellGrid { spacing: f64::NAN, ..g },
            CellGrid { origin: [f64::INFINITY, 0.0, 0.0], ..g },
            CellGrid { cells: [0, 3, 1], ..g },
            CellGrid { cells: [usize::MAX, usize::MAX, 2], ..g },
        ];
        for b in bad {
            assert!(slice(&b, &values, 2, 0.5).is_none(), "{b:?}");
            assert!(sample(&b, &values, Vec3::ZERO).is_nan());
            assert!(streamlines(&b, &[[1.0, 0.0, 0.0]; 9], &[], &[Vec3::new(0.5, 0.5, 0.5)], 10).is_empty());
            assert!(face_seeds(&b, &[], 0, false, 10).is_empty());
        }
    }

    /// A huge layer is sampled every few cells, keeping its extent.
    #[test]
    fn a_large_slice_is_sampled() {
        let g = grid([2048, 1024, 1]);
        let values = vec![1.0; g.cell_count().unwrap()];
        let s = slice(&g, &values, 2, 0.5).unwrap();
        assert!(s.tris.len() <= 2 * MAX_SLICE_QUADS, "{}", s.tris.len());
        let xs = s.tris.iter().flat_map(|t| t.p).map(|p| p.x).fold(0.0, f64::max);
        assert_eq!(xs, 2048.0);
    }

    /// Sampling beside a solid block takes the fluid cells' values only.
    #[test]
    fn sampling_beside_solid_cells_reads_the_fluid() {
        let g = grid([4, 1, 1]);
        let values = [f64::NAN, f64::NAN, 2.0, 4.0];
        // On the block's face at x = 2: cells 1 (solid) and 2.
        assert_eq!(sample(&g, &values, Vec3::new(2.0, 0.5, 0.5)), 2.0);
        assert_eq!(sample(&g, &values, Vec3::new(3.0, 0.5, 0.5)), 3.0);
        // Deep inside the block, the nearest fluid cell; beyond the grid, nothing.
        assert_eq!(sample(&g, &values, Vec3::new(0.2, 0.5, 0.5)), 2.0);
        assert!(sample(&g, &values, Vec3::new(9.0, 0.5, 0.5)).is_nan());
        assert!(sample(&g, &[f64::NAN; 4], Vec3::new(1.0, 0.5, 0.5)).is_nan());
        assert!(sample(&g, &values, Vec3::new(f64::NAN, 0.5, 0.5)).is_nan());
    }

    /// A uniform flow along x carries a streamline straight to the outlet in half-cell steps; a
    /// solid cell stops it; still fluid and seeds outside give none.
    #[test]
    fn streamlines_follow_the_flow() {
        let g = grid([10, 3, 3]);
        let n = g.cell_count().unwrap();
        let u = vec![[2.0f32, 0.0, 0.0]; n];
        let seed = Vec3::new(0.5, 1.5, 1.5);
        let lines = streamlines(&g, &u, &[], &[seed], 1000);
        assert_eq!(lines.len(), 1);
        let l = &lines[0];
        assert!(l.points.iter().all(|p| (p.y - 1.5).abs() < 1e-12 && (p.z - 1.5).abs() < 1e-12));
        let last = *l.points.last().unwrap();
        assert!(last.x > 9.4 && last.x < 10.0, "{last:?}");
        assert_eq!(l.points.len(), 19);
        assert!(l.speed.iter().all(|s| (s - 2.0).abs() < 1e-6));
        // A step cap.
        assert_eq!(streamlines(&g, &u, &[], &[seed], 3)[0].points.len(), 4);
        // A solid cell across the way at x = 5.
        let mut solid = vec![false; n];
        solid[g.index([5, 1, 1])] = true;
        let blocked = streamlines(&g, &u, &solid, &[seed], 1000);
        assert!(blocked[0].points.iter().all(|p| p.x < 5.0), "{:?}", blocked[0].points.last());
        // Outside, in a solid cell, in still fluid: no streamline.
        let none = [Vec3::new(-1.0, 1.5, 1.5), Vec3::new(5.5, 1.5, 1.5), Vec3::new(f64::NAN, 0.0, 0.0)];
        assert!(streamlines(&g, &u, &solid, &none, 1000).is_empty());
        assert!(streamlines(&g, &vec![[0.0; 3]; n], &[], &[seed], 1000).is_empty());
        assert!(streamlines(&g, &u[..5], &[], &[seed], 1000).is_empty());
        assert!(streamlines(&g, &u, &solid[..5], &[seed], 1000).is_empty());
    }

    /// A rigid rotation about z: streamlines are circles, the midpoint rule keeps their radius.
    #[test]
    fn streamlines_of_a_vortex_close() {
        let g = CellGrid { origin: [-10.0, -10.0, 0.0], spacing: 0.5, cells: [40, 40, 2] };
        let u: Vec<[f32; 3]> = (0..g.cell_count().unwrap())
            .map(|c| {
                let p = g.centre([c % 40, (c / 40) % 40, c / 1600]);
                [-p.y as f32, p.x as f32, 0.0]
            })
            .collect();
        let lines = streamlines(&g, &u, &[], &[Vec3::new(5.0, 0.0, 0.5)], 2000);
        let l = &lines[0];
        assert_eq!(l.points.len(), 2001, "a closed streamline runs to its step cap");
        for p in &l.points {
            assert!(((p.x * p.x + p.y * p.y).sqrt() - 5.0).abs() < 0.05, "{p:?}");
        }
    }

    /// Hostile counts are capped: seeds, steps and points.
    #[test]
    fn streamline_work_is_bounded() {
        let g = grid([64, 64, 64]);
        let n = g.cell_count().unwrap();
        let u = vec![[0.0f32, 0.0, 1.0]; n];
        let seeds = face_seeds(&g, &[], 2, false, usize::MAX);
        assert_eq!(seeds.len(), 44 * 45, "capped at {MAX_STREAMLINES}");
        assert!(seeds.len() <= MAX_STREAMLINES);
        let lines = streamlines(&g, &u, &[], &seeds, usize::MAX);
        let points: usize = lines.iter().map(|l| l.points.len()).sum();
        assert!(points <= MAX_STREAMLINE_POINTS);
        assert!(lines.iter().all(|l| l.points.len() == 127), "each runs up the 64 cells in 126 steps");
    }

    /// Seeds cover the face, or the obstacle's shadow on it grown by half its size each way.
    #[test]
    fn seeds_lie_upstream_of_the_obstacle() {
        let g = grid([10, 10, 10]);
        let all = face_seeds(&g, &[], 0, false, 25);
        assert_eq!(all.len(), 25);
        assert!(all.iter().all(|p| p.x == 0.5 && (0.0..10.0).contains(&p.y) && (0.0..10.0).contains(&p.z)));
        assert_eq!(all[0], Vec3::new(0.5, 1.0, 1.0));
        // A 2 × 2 block of cells at y, z in [4, 6]: seeds over y, z in [3, 7].
        let mut solid = vec![false; 1000];
        for (j, k) in [(4, 4), (4, 5), (5, 4), (5, 5)] {
            solid[g.index([7, j, k])] = true;
        }
        let near = face_seeds(&g, &solid, 0, false, 16);
        assert_eq!(near.len(), 16);
        assert!(near.iter().all(|p| (3.0..=7.0).contains(&p.y) && (3.0..=7.0).contains(&p.z)), "{near:?}");
        let high = face_seeds(&g, &solid, 0, true, 4);
        assert!(high.iter().all(|p| p.x == 9.5));
        assert!(face_seeds(&g, &solid, 3, false, 4).is_empty());
        assert!(face_seeds(&g, &solid, 0, false, 0).is_empty());
        // A wide face gets more seeds across its width.
        let wide = face_seeds(&grid([1, 40, 10]), &[], 0, false, 16);
        assert_eq!(wide.len(), 8 * 2);
    }
}
