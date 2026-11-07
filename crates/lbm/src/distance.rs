//! Every fluid cell's distance to the nearest no-slip wall, for turbulence models that read it
//! (Spalart–Allmaras' destruction term).
//!
//! **The wall** is where the solver puts it: a halfway bounce-back link puts it at the link's
//! midpoint ([`WallSurface::from_links`]); Han, Ooka and Kikumoto's on-node wall-function bounce
//! puts it on the boundary node ([`WallSurface::from_nodes`]). Slip, inlet, outlet and far-field
//! faces are not walls. Every such point lies on the lattice of half cells, so the points are held
//! exactly, in doubled integer coordinates.
//!
//! **The distance** is the exact Euclidean one from each fluid cell's centre to the nearest wall
//! point (through solids too, across periodic faces: the nearest image), by the separable distance
//! transform of P. F. Felzenszwalb and D. P. Huttenlocher, *Distance transforms of sampled
//! functions*, Theory of Computing 8:415–428, 2012 (doi:10.4086/toc.2012.v008a019): lower envelopes
//! of parabolas, one axis at a time. The wall points split into eight classes by which coordinates
//! are half-integers, each taking the three passes on its own lattice; each cell keeps the least.
//! Integers throughout and breakpoints compared as exact fractions make the result exact (to the
//! final square root) and the same on any number of threads.
//!
//! **Work** is linear: at most 24 sweeps over about `(n_x + 2)(n_y + 2)(n_z + 2)` slots, no
//! iteration. Each axis is capped at [`MAX_AXIS`], the grid at `Limits::max_cells` and the working
//! memory at `Limits::max_memory`. With no wall point at all every fluid cell takes the box
//! diagonal ([`WallDistance::has_walls`]). Distances are in lattice units; solid cells hold 0.
//! Separate from Van Driest's distance in [`crate::wall`], which does not cross periodic faces.

use fieldcraft_parallel as parallel;

use crate::flow::{Link, Rule, filled};
use crate::lattice::C;
use crate::voxel::Grid;
use crate::{Error, Limits, Result};

/// The most cells along one axis.
pub const MAX_AXIS: usize = 1 << 24;

/// No sample: a squared distance larger than any.
const FAR: i64 = i64::MAX;

/// The no-slip wall points the solver uses, in doubled lattice coordinates (`2 p`), sorted and
/// without repeats.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct WallSurface {
    points: Vec<[i64; 3]>,
}

impl WallSurface {
    /// The wall points of a halfway bounce-back: each link of `links` (the flow's) that bounces
    /// back or reflects at a force group `g` with `no_slip[g]` puts a wall at its midpoint,
    /// `x − c_q / 2` for direction `q` of the cell at `x` (the population arrives from `−c_q`).
    pub fn from_links(grid: &Grid, links: &[Link], no_slip: &[bool]) -> Result<WallSurface> {
        let n = cell_count(grid)?;
        let mut points = Vec::new();
        for l in links {
            let (Rule::Bounce { group, .. } | Rule::Reflect { group, .. }) = l.rule else { continue };
            if !no_slip.get(group as usize).copied().unwrap_or(false) {
                continue;
            }
            let cell = l.cell as usize;
            let Some(cq) = C.get(usize::from(l.q)).filter(|_| cell < n) else {
                return Err(Error::Invalid(format!("a wall link (cell {}, direction {}) lies outside the grid", l.cell, l.q)));
            };
            let ijk = grid.coords(cell);
            points.push([0, 1, 2].map(|a| 2 * ijk[a] as i64 - i64::from(cq[a])));
        }
        Ok(WallSurface::sorted(points))
    }

    /// The wall points of the on-node wall-function bounce: `nodes` hold `(interior node, axis,
    /// inward, group)` as the flow lists them, the wall on the boundary node one cell from the
    /// interior node against `inward` (±1) along `axis`.
    pub fn from_nodes(grid: &Grid, nodes: &[(u32, usize, i32, u32)]) -> Result<WallSurface> {
        let n = cell_count(grid)?;
        let mut points = Vec::new();
        for &(cell, axis, inward, _) in nodes {
            if cell as usize >= n || axis > 2 || inward.abs() != 1 {
                return Err(Error::Invalid(format!("a wall node (cell {cell}, axis {axis}, inward {inward}) is not on the grid")));
            }
            let ijk = grid.coords(cell as usize);
            points.push([0, 1, 2].map(|a| 2 * ijk[a] as i64 - if a == axis { 2 * i64::from(inward) } else { 0 }));
        }
        Ok(WallSurface::sorted(points))
    }

    /// Wall points given in lattice coordinates, each coordinate a multiple of ½.
    pub fn from_points(points: &[[f64; 3]]) -> Result<WallSurface> {
        let limit = 4.0 * MAX_AXIS as f64;
        let mut doubled = Vec::with_capacity(points.len());
        for p in points {
            let d = p.map(|v| 2.0 * v);
            if d.iter().any(|v| !v.is_finite() || v.abs() > limit || v.fract() != 0.0) {
                return Err(Error::Invalid(format!("the wall point {p:?} is not on the lattice of half cells")));
            }
            doubled.push(d.map(|v| v as i64));
        }
        Ok(WallSurface::sorted(doubled))
    }

    fn sorted(mut points: Vec<[i64; 3]>) -> WallSurface {
        points.sort_unstable();
        points.dedup();
        WallSurface { points }
    }

    /// The number of wall points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// The wall points in lattice coordinates, in a fixed order.
    pub fn points(&self) -> impl Iterator<Item = [f64; 3]> + '_ {
        self.points.iter().map(|p| p.map(|v| v as f64 * 0.5))
    }
}

/// Every fluid cell's distance to the nearest wall point (lattice units), 0 in other cells.
#[derive(Debug, Clone, PartialEq)]
pub struct WallDistance {
    d: Vec<f64>,
    walls: bool,
}

impl WallDistance {
    /// The distance from every `fluid` cell of `grid` to the nearest of `walls`, across the faces
    /// of the `periodic` axes. Refused, with what to change, when the grid has no cells, an axis
    /// is longer than [`MAX_AXIS`], the grid is above `limits.max_cells`, the working memory is
    /// above `limits.max_memory`, or a wall point lies more than a cell outside the box.
    pub fn new(grid: &Grid, fluid: &[bool], walls: &WallSurface, periodic: [bool; 3], limits: &Limits) -> Result<WallDistance> {
        let cells = grid.cells;
        limits.check_cells(cells.map(|c| c as u64))?;
        let n = cell_count(grid)?;
        if fluid.len() != n {
            return Err(Error::Invalid(format!("{} fluid flags for a grid of {n} cells", fluid.len())));
        }
        // The working memory: both intermediate passes, the squared distances and the result.
        let [nx, ny, nz] = cells.map(|c| c as u64);
        let bytes = [nx.saturating_mul(ny + 2).saturating_mul(nz + 2), nx.saturating_mul(ny).saturating_mul(nz + 2), (n as u64).saturating_mul(2)]
            .iter()
            .fold(0u64, |s, v| s.saturating_add(*v))
            .saturating_mul(8);
        if bytes > limits.max_memory {
            return Err(Error::Limit(format!(
                "the wall distance needs {bytes} bytes, above the memory cap of {}: coarsen the spacing",
                limits.max_memory
            )));
        }

        let mut classes: [Vec<[i64; 3]>; 8] = Default::default();
        for p in &walls.points {
            let mut q = *p;
            for a in 0..3 {
                let span = 2 * cells[a] as i64;
                if periodic[a] {
                    q[a] = q[a].rem_euclid(span);
                } else if q[a] < -2 || q[a] > span {
                    return Err(Error::Invalid(format!("the wall point {:?} lies outside the box", p.map(|v| v as f64 * 0.5))));
                }
            }
            let class = (0..3).map(|a| ((q[a] & 1) as usize) << a).sum::<usize>();
            if let Some(c) = classes.get_mut(class) {
                c.push(q);
            }
        }

        let mut best = filled(n, FAR, "wall distances")?;
        for (class, points) in classes.iter().enumerate() {
            if !points.is_empty() {
                let parity = [0, 1, 2].map(|a| ((class >> a) & 1) as i64);
                transform(cells, periodic, parity, points, &mut best)?;
            }
        }

        let diagonal = cells.iter().map(|c| (*c as f64).powi(2)).sum::<f64>().sqrt();
        let mut d = filled(n, 0.0, "wall distances")?;
        for (cell, out) in d.iter_mut().enumerate() {
            if !fluid.get(cell).copied().unwrap_or(false) {
                continue;
            }
            let [x, y, z] = grid.coords(cell);
            *out = match best.get(z + nz as usize * (x + nx as usize * y)) {
                Some(&s) if s != FAR => (s as f64).sqrt() * 0.5,
                _ => diagonal,
            };
        }
        Ok(WallDistance { d, walls: !walls.is_empty() })
    }

    /// Cell `cell`'s distance to the nearest wall (lattice units); `None` outside the grid, for
    /// cells that are not fluid, and for a fluid cell at distance 0 (one centred on a wall point,
    /// such as an on-node boundary node passed as fluid).
    pub fn get(&self, cell: usize) -> Option<f64> {
        self.d.get(cell).copied().filter(|d| *d > 0.0)
    }

    /// Every cell's distance (lattice units), in cell order: 0 in cells that are not fluid.
    pub fn as_slice(&self) -> &[f64] {
        &self.d
    }

    /// Whether there was any wall point (without one every fluid cell holds the box diagonal).
    pub fn has_walls(&self) -> bool {
        self.walls
    }
}

/// The grid's cell count, refused when it is 0, an axis is longer than [`MAX_AXIS`] or the count
/// overflows.
fn cell_count(grid: &Grid) -> Result<usize> {
    let [nx, ny, nz] = grid.cells;
    if grid.cells.iter().any(|c| *c == 0 || *c > MAX_AXIS) {
        return Err(Error::Invalid(format!("a grid of {nx} × {ny} × {nz} cells: each axis needs 1 to {MAX_AXIS} cells for the wall distance")));
    }
    nx.checked_mul(ny).and_then(|v| v.checked_mul(nz)).ok_or_else(|| Error::Limit(format!("{nx} × {ny} × {nz} cells cannot be addressed")))
}

/// One class's three passes: `points` (doubled, wrapped along periodic axes) all have the
/// coordinate parities `parity`; each cell's squared doubled distance to the nearest of them is
/// kept in `best` where it is smaller (`best` laid out z fastest, then x, then y).
fn transform(cells: [usize; 3], periodic: [bool; 3], parity: [i64; 3], points: &[[i64; 3]], best: &mut [i64]) -> Result<()> {
    let [nx, ny, nz] = cells;
    // Slot m along an axis sits at doubled position 2 (m − 1) + parity: slots 0 ..= n + 1 hold
    // every point from one cell below the box to one above it.
    let (ly, lz) = (ny + 2, nz + 2);
    let slot = |p: i64, a: usize| ((p - parity[a]) / 2 + 1) as usize;
    let position = |m: usize, a: usize| 2 * (m as i64 - 1) + parity[a];

    let mut keyed: Vec<(usize, i64)> = points.iter().map(|p| (slot(p[1], 1) + ly * slot(p[2], 2), p[0])).collect();
    keyed.sort_unstable();
    let lines = ly * lz;
    let mut starts = filled(lines + 1, 0usize, "wall distance lines")?;
    for (line, _) in &keyed {
        if let Some(s) = starts.get_mut(line + 1) {
            *s += 1;
        }
    }
    for i in 0..lines {
        let before = starts.get(i).copied().unwrap_or(0);
        if let Some(s) = starts.get_mut(i + 1) {
            *s += before;
        }
    }
    let mut g1 = filled(nx * lines, FAR, "wall distance pass")?;
    parallel::chunks_mut(&mut g1, nx, |line, out| {
        let (a, b) = (starts.get(line).copied().unwrap_or(0), starts.get(line + 1).copied().unwrap_or(0));
        let samples: Vec<(i64, i64)> = keyed.get(a..b).unwrap_or(&[]).iter().map(|(_, x)| (*x, 0)).collect();
        envelope(&images(samples, periodic[0], nx), out);
    });

    let mut g2 = filled(ny * nx * lz, FAR, "wall distance pass")?;
    parallel::chunks_mut(&mut g2, ny, |line, out| {
        let (x, mz) = (line % nx, line / nx);
        let samples: Vec<(i64, i64)> =
            (0..ly).filter_map(|my| g1.get(x + nx * (my + ly * mz)).filter(|v| **v != FAR).map(|v| (position(my, 1), *v))).collect();
        envelope(&images(samples, periodic[1], ny), out);
    });
    drop(g1);

    parallel::chunks_mut(best, nz, |line, out| {
        let (x, y) = (line % nx, line / nx);
        let samples: Vec<(i64, i64)> =
            (0..lz).filter_map(|mz| g2.get(y + ny * (x + nx * mz)).filter(|v| **v != FAR).map(|v| (position(mz, 2), *v))).collect();
        let mut here = vec![FAR; out.len()];
        envelope(&images(samples, periodic[2], nz), &mut here);
        for (o, h) in out.iter_mut().zip(here) {
            *o = (*o).min(h);
        }
    });
    Ok(())
}

/// The samples of a line with, along a periodic axis of `n` cells (positions wrapped into
/// `0 .. 2n`), their images one period below and above: still sorted, and the nearest image of
/// every point is among them for any cell of the line.
fn images(samples: Vec<(i64, i64)>, periodic: bool, n: usize) -> Vec<(i64, i64)> {
    if !periodic || samples.is_empty() {
        return samples;
    }
    let span = 2 * n as i64;
    [-span, 0, span].iter().flat_map(|s| samples.iter().map(move |(p, f)| (p + s, *f))).collect()
}

/// The lower envelope of the parabolas `(X − p)² + f` of `samples` (doubled positions `p`
/// strictly increasing, `f` finite) at the cell centres `X = 2 i` of `out` (Felzenszwalb and
/// Huttenlocher's Algorithm 1). The breakpoint of two parabolas is the fraction
/// `((f_q + p_q²) − (f_v + p_v²)) / (2 (p_q − p_v))`, compared exactly in `i128`.
fn envelope(samples: &[(i64, i64)], out: &mut [i64]) {
    if samples.is_empty() {
        out.fill(FAR);
        return;
    }
    // The parabolas of the envelope, and each one's left breakpoint (unused for the first).
    let mut v: Vec<(i64, i64)> = Vec::with_capacity(samples.len());
    let mut z: Vec<(i128, i128)> = Vec::with_capacity(samples.len());
    for &(pq, fq) in samples {
        loop {
            let Some(&(pv, fv)) = v.last() else {
                v.push((pq, fq));
                z.push((0, 1));
                break;
            };
            let num = (i128::from(fq) + i128::from(pq) * i128::from(pq)) - (i128::from(fv) + i128::from(pv) * i128::from(pv));
            let den = 2 * (i128::from(pq) - i128::from(pv));
            if v.len() > 1
                && let Some(&(zn, zd)) = z.last()
                && num * zd <= zn * den
            {
                v.pop();
                z.pop();
                continue;
            }
            v.push((pq, fq));
            z.push((num, den));
            break;
        }
    }
    let mut k = 0;
    for (i, o) in out.iter_mut().enumerate() {
        let x = 2 * i as i64;
        while let Some(&(zn, zd)) = z.get(k + 1)
            && zn < i128::from(x) * zd
        {
            k += 1;
        }
        *o = match v.get(k) {
            Some(&(p, f)) => (x - p) * (x - p) + f,
            None => FAR,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(cells: [usize; 3]) -> Grid {
        Grid { origin: [0.0; 3], spacing: 1.0, cells }
    }

    /// The links a flow would build: every fluid cell pulls each direction from `x − c_q`; a solid
    /// cell there, or a face that is not periodic, is a bounce-back wall (group 0).
    fn links(g: &Grid, fluid: &[bool], periodic: [bool; 3]) -> Vec<Link> {
        let mut out = Vec::new();
        for cell in 0..g.len() {
            if !fluid[cell] {
                continue;
            }
            let ijk = g.coords(cell);
            for q in 1..C.len() {
                let mut src = [0usize; 3];
                let mut wall = false;
                for a in 0..3 {
                    let s = ijk[a] as i64 - i64::from(C[q][a]);
                    let n = g.cells[a] as i64;
                    if s < 0 || s >= n {
                        if periodic[a] {
                            src[a] = s.rem_euclid(n) as usize;
                        } else {
                            wall = true;
                        }
                    } else {
                        src[a] = s as usize;
                    }
                }
                if wall || !fluid[g.index(src)] {
                    out.push(Link { cell: cell as u32, q: q as u8, rule: Rule::Bounce { group: 0, du: 0.0 } });
                }
            }
        }
        out
    }

    /// The distance to the nearest point by checking every one, with the nearest periodic image.
    fn brute(g: &Grid, fluid: &[bool], walls: &WallSurface, periodic: [bool; 3]) -> Vec<f64> {
        let pts: Vec<[f64; 3]> = walls.points().collect();
        (0..g.len())
            .map(|cell| {
                if !fluid[cell] {
                    return 0.0;
                }
                let x = g.coords(cell).map(|v| v as f64);
                let mut best = f64::INFINITY;
                for p in &pts {
                    let mut s = 0.0;
                    for a in 0..3 {
                        let mut d = (x[a] - p[a]).abs();
                        if periodic[a] {
                            let n = g.cells[a] as f64;
                            d = d.rem_euclid(n);
                            d = d.min(n - d);
                        }
                        s += d * d;
                    }
                    best = best.min(s);
                }
                best.sqrt()
            })
            .collect()
    }

    /// A small deterministic generator (SplitMix64).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }

    /// Random boxes of solid cells in an `n³` grid, the periodic axes from `combo`'s bits.
    fn random_case(rng: &mut Rng, n: usize, combo: usize) -> (Grid, Vec<bool>, [bool; 3]) {
        let g = grid([n, n, n]);
        let periodic = [0, 1, 2].map(|a| combo >> a & 1 == 1);
        let mut fluid = vec![true; g.len()];
        for _ in 0..4 {
            let lo = [0; 3].map(|_| (rng.next() % n as u64) as usize);
            let size = [0; 3].map(|_| 1 + (rng.next() % (n as u64 / 3)) as usize);
            for k in lo[2]..(lo[2] + size[2]).min(n) {
                for j in lo[1]..(lo[1] + size[1]).min(n) {
                    for i in lo[0]..(lo[0] + size[0]).min(n) {
                        fluid[g.index([i, j, k])] = false;
                    }
                }
            }
        }
        (g, fluid, periodic)
    }

    fn check_random(n: usize, cases: usize) {
        let mut rng = Rng(7);
        for case in 0..cases {
            let (g, fluid, periodic) = random_case(&mut rng, n, case % 8);
            let walls = WallSurface::from_links(&g, &links(&g, &fluid, periodic), &[true]).unwrap();
            let d = WallDistance::new(&g, &fluid, &walls, periodic, &Limits::default()).unwrap();
            let expected = brute(&g, &fluid, &walls, periodic);
            if walls.is_empty() {
                continue;
            }
            for (cell, (a, b)) in d.as_slice().iter().zip(&expected).enumerate() {
                assert!((a - b).abs() <= 1e-12, "case {case} periodic {periodic:?} cell {cell}: {a} vs brute force {b}");
            }
        }
    }

    #[test]
    fn equals_brute_force_on_random_obstacles() {
        check_random(10, 16);
    }

    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn equals_brute_force_on_random_obstacles_24() {
        check_random(24, 20);
    }

    #[test]
    fn halfway_plane_walls_are_half_a_cell_off() {
        // A channel between y = −½ and y = n_y − ½, periodic along x and z.
        let g = grid([5, 8, 3]);
        let fluid = vec![true; g.len()];
        let periodic = [true, false, true];
        let walls = WallSurface::from_links(&g, &links(&g, &fluid, periodic), &[true]).unwrap();
        let d = WallDistance::new(&g, &fluid, &walls, periodic, &Limits::default()).unwrap();
        for cell in 0..g.len() {
            let j = g.coords(cell)[1] as f64;
            assert_eq!(d.get(cell), Some((j + 0.5).min(7.5 - j)), "cell {cell}");
        }
        // One wall: d = j + ½ exactly.
        let one: Vec<Link> = links(&g, &fluid, periodic).into_iter().filter(|l| C[usize::from(l.q)][1] == 1).collect();
        let walls = WallSurface::from_links(&g, &one, &[true]).unwrap();
        let d = WallDistance::new(&g, &fluid, &walls, periodic, &Limits::default()).unwrap();
        for cell in 0..g.len() {
            assert_eq!(d.get(cell), Some(g.coords(cell)[1] as f64 + 0.5));
        }
    }

    #[test]
    fn slip_groups_are_not_walls_and_on_node_walls_sit_on_the_node() {
        let g = grid([4, 6, 4]);
        let fluid = vec![true; g.len()];
        let all = links(&g, &fluid, [true, false, true]);
        assert!(WallSurface::from_links(&g, &all, &[false]).unwrap().is_empty());
        assert!(WallSurface::from_links(&g, &all, &[]).unwrap().is_empty());
        // Interior nodes j = 1 whose boundary node is j = 0: the wall is at y = 0.
        let nodes: Vec<(u32, usize, i32, u32)> = (0..g.len()).filter(|c| g.coords(*c)[1] == 1).map(|c| (c as u32, 1, 1, 0)).collect();
        let walls = WallSurface::from_nodes(&g, &nodes).unwrap();
        assert!(walls.points().all(|p| p[1] == 0.0));
        // The boundary nodes do not collide: they are left out of the fluid, as the flow does.
        let bulk: Vec<bool> = (0..g.len()).map(|c| g.coords(c)[1] > 0).collect();
        let d = WallDistance::new(&g, &bulk, &walls, [true, false, true], &Limits::default()).unwrap();
        for c in 0..g.len() {
            let j = g.coords(c)[1];
            assert_eq!(d.get(c), (j > 0).then_some(j as f64), "cell {c}");
        }
    }

    /// A sphere of radius `r` cells at the centre of an `n³` periodic box: the largest error of the
    /// distance against the sphere's, `|x − c| − r`, over the fluid cells (cells).
    fn sphere_error(n: usize, r: f64) -> f64 {
        let g = grid([n, n, n]);
        let c = (n as f64 - 1.0) / 2.0;
        let dist = |cell: usize| g.coords(cell).iter().map(|v| (*v as f64 - c).powi(2)).sum::<f64>().sqrt();
        let fluid: Vec<bool> = (0..g.len()).map(|cell| dist(cell) > r).collect();
        let periodic = [true; 3];
        let walls = WallSurface::from_links(&g, &links(&g, &fluid, periodic), &[true]).unwrap();
        let d = WallDistance::new(&g, &fluid, &walls, periodic, &Limits::default()).unwrap();
        // Within the box's inscribed shell, where the sphere is the nearest wall in every image.
        (0..g.len()).filter(|cell| fluid[*cell] && dist(*cell) < c).map(|cell| (d.as_slice()[cell] - (dist(cell) - r)).abs()).fold(0.0, f64::max)
    }

    #[test]
    fn a_sphere_is_within_a_cell_and_converges() {
        let coarse = sphere_error(16, 4.3);
        let fine = sphere_error(32, 8.6);
        assert!(coarse < 1.0 && fine < 1.0, "errors {coarse} and {fine} cells");
        // Relative to the radius the error shrinks with the grid.
        assert!(fine / 8.6 < 0.75 * coarse / 4.3, "relative errors {} then {}", coarse / 4.3, fine / 8.6);
    }

    #[test]
    fn crosses_periodic_faces() {
        // One wall plane at x = −½ ≡ n − ½: the last cell is half a cell from it.
        let g = grid([6, 2, 2]);
        let fluid = vec![true; g.len()];
        let walls = WallSurface::from_points(&[[-0.5, 0.0, 0.0], [-0.5, 1.0, 0.0], [-0.5, 0.0, 1.0], [-0.5, 1.0, 1.0]]).unwrap();
        let d = WallDistance::new(&g, &fluid, &walls, [true; 3], &Limits::default()).unwrap();
        for cell in 0..g.len() {
            let i = g.coords(cell)[0] as f64;
            assert_eq!(d.get(cell), Some((i + 0.5).min(5.5 - i)));
        }
    }

    #[test]
    fn same_on_any_number_of_threads() {
        let mut rng = Rng(11);
        let (g, fluid, periodic) = random_case(&mut rng, 12, 5);
        let walls = WallSurface::from_links(&g, &links(&g, &fluid, periodic), &[true]).unwrap();
        let run = |t| parallel::with_threads(t, || WallDistance::new(&g, &fluid, &walls, periodic, &Limits::default()).unwrap());
        let one = run(1);
        assert_eq!(one, run(2));
        assert_eq!(one, run(7));
    }

    #[test]
    fn hostile_input_is_refused_without_panic() {
        let limits = Limits::default();
        let g = grid([4, 4, 4]);
        let fluid = vec![true; 64];
        let none = WallSurface::default();
        // No walls: the box diagonal everywhere, finite.
        let d = WallDistance::new(&g, &fluid, &none, [true; 3], &limits).unwrap();
        assert!(!d.has_walls());
        assert!(d.as_slice().iter().all(|v| *v == 48f64.sqrt()));
        // Solid cells hold 0 and read None.
        let mut some = fluid.clone();
        some[3] = false;
        let d = WallDistance::new(&g, &some, &none, [true; 3], &limits).unwrap();
        assert_eq!((d.as_slice()[3], d.get(3), d.get(64)), (0.0, None, None));
        // Empty or oversized grids, a wrong mask, a memory cap.
        for cells in [[0, 4, 4], [MAX_AXIS + 1, 1, 1], [usize::MAX, usize::MAX, 2]] {
            assert!(WallDistance::new(&grid(cells), &[], &none, [false; 3], &limits).is_err(), "{cells:?}");
        }
        assert!(WallDistance::new(&g, &fluid[..10], &none, [false; 3], &limits).is_err());
        let tight = Limits { max_memory: 100, ..limits };
        assert!(matches!(WallDistance::new(&g, &fluid, &none, [false; 3], &tight), Err(Error::Limit(_))));
        let small = Limits { max_cells: 10, ..limits };
        assert!(matches!(WallDistance::new(&g, &fluid, &none, [false; 3], &small), Err(Error::Limit(_))));
        // Points off the half-cell lattice, not finite, or outside a box that is not periodic.
        for p in [[0.25, 0.0, 0.0], [f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0], [1e300, 0.0, 0.0]] {
            assert!(WallSurface::from_points(&[p]).is_err(), "{p:?}");
        }
        let far = WallSurface::from_points(&[[9.0, 0.0, 0.0]]).unwrap();
        assert!(WallDistance::new(&g, &fluid, &far, [false; 3], &limits).is_err());
        assert!(WallDistance::new(&g, &fluid, &far, [true; 3], &limits).is_ok());
        // Links and nodes off the grid.
        let bad = [Link { cell: 64, q: 1, rule: Rule::Bounce { group: 0, du: 0.0 } }];
        assert!(WallSurface::from_links(&g, &bad, &[true]).is_err());
        let bad = [Link { cell: 0, q: 40, rule: Rule::Bounce { group: 0, du: 0.0 } }];
        assert!(WallSurface::from_links(&g, &bad, &[true]).is_err());
        for node in [(64, 0, 1, 0), (0, 3, 1, 0), (0, 0, 2, 0)] {
            assert!(WallSurface::from_nodes(&g, &[node]).is_err(), "{node:?}");
        }
    }

    /// The work grows linearly: eight times the cells take no more than twice eight times as long.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn work_is_linear_in_the_cells() {
        let time = |n: usize| {
            let g = grid([n, n, n]);
            let c = n as f64 / 2.0;
            let fluid: Vec<bool> =
                (0..g.len()).map(|cell| g.coords(cell).iter().map(|v| (*v as f64 - c).powi(2)).sum::<f64>() > (n as f64 / 5.0).powi(2)).collect();
            let periodic = [true, false, true];
            let walls = WallSurface::from_links(&g, &links(&g, &fluid, periodic), &[true]).unwrap();
            let start = std::time::Instant::now();
            for _ in 0..3 {
                std::hint::black_box(WallDistance::new(&g, &fluid, &walls, periodic, &Limits::default()).unwrap());
            }
            start.elapsed().as_secs_f64() / 3.0
        };
        let (small, large) = (time(64), time(128));
        eprintln!("wall distance: 64³ {small:.4} s, 128³ {large:.4} s, ratio {:.2}", large / small);
        assert!(large <= 2.0 * 8.0 * small, "64³ {small} s, 128³ {large} s");
    }
}
