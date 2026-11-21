//! The wall nodes of the half-way wall: the first solid layer.
//!
//! With `spalding_bounce_surface` every cell whose centre is in the air stays open (fluid); the wall
//! is the face half way between the last open centre and the first solid centre, and that first
//! solid cell is the boundary node of Han, Ooka and Kikumoto's on-node wall-function bounce (their
//! eq. 14, a node with the wall through it), so the wall sits half a cell inside the surface-true
//! geometry the voxeliser gives. The node code of [`crate::flow`] takes any cell flagged fluid that
//! has a solid face neighbour for a node, so the rule is a change of that mask before it runs:
//! [`promote`] flags the node cells fluid and the loop does the rest.
//!
//! A solid cell becomes a node along an axis when
//! 1. the cell one step towards the air along that axis is open,
//! 2. the cell one step the other way is an obstacle's solid cell (a wall behind the node: without
//!    one the node has no wall direction and would be plain fluid, a leak), and
//! 3. the cell two steps behind is not open (a wall of two solid cells, open on both sides, would
//!    otherwise turn both cells into nodes with no solid between them).
//!
//! A cell that is not a node stays solid and the open cell beside it keeps plain half-way
//! bounce-back, which seals the wall and gives it no wall law. This happens to walls of fewer than
//! three solid centres across, and to cells in the layer of an inlet, outlet or slip face of the box
//! (the node loop gives those cells the face's rule, a velocity inlet inside the body).
//!
//! The decision reads the voxeliser's masks as they stand before any cell is promoted, so nodes do
//! not chain, and it visits cells in index order: the result does not depend on threads.

use crate::flow::{NONE, filled};
use crate::voxel::Grid;
use crate::{Error, Result};
use fieldcraft_mesh::surface::BoundarySurface;

/// What [`promote`] did, per grid.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WallLayerReport {
    /// Solid cells flagged as wall nodes.
    pub nodes: usize,
    /// Wall faces: pairs of an open cell and an obstacle's solid cell that share a face.
    pub faces: usize,
    /// Of those, the faces whose solid cell is not a node (a wall of fewer than three solid
    /// centres, or the layer of an open box face): plain half-way bounce-back, no wall law.
    pub plain_faces: usize,
    /// The share of the obstacles' surface area (inside the box) that no solid cell covers: a wall
    /// thinner than the spacing has no centre inside it, and the flow leaks through. 0 where not
    /// measured.
    pub unsealed_share: f64,
    /// Of the nodes, those with fewer than [`NARROW_CELLS`] open cells between the wall and the
    /// opposite wall (or an obstacle) along the normal: a passage too narrow for the law's
    /// second cell to lie in the log layer. 0 where not measured.
    pub narrow_nodes: usize,
}

impl WallLayerReport {
    /// The share of the wall area that runs plain bounce-back, 0 without walls.
    pub fn plain_share(&self) -> f64 {
        if self.faces == 0 { 0.0 } else { self.plain_faces as f64 / self.faces as f64 }
    }

    /// The share of the nodes in a passage narrower than [`NARROW_CELLS`] cells, 0 without nodes.
    pub fn narrow_share(&self) -> f64 {
        if self.nodes == 0 { 0.0 } else { self.narrow_nodes as f64 / self.nodes as f64 }
    }
}

/// The open cells to the opposite wall below which setup warns that a wall's passage is narrow: the
/// law reads the second open cell, which needs room above it to lie in the log layer and not in the
/// passage's core (a choice, not a measurement).
pub const NARROW_CELLS: usize = 12;

/// How many of `starts` (an interior cell, its wall's axis and the sign of its inward normal) have
/// fewer than [`NARROW_CELLS`] open cells (`open`) in a row along the normal before a cell that is
/// not open. A run that leaves the box, or that is still open at [`NARROW_CELLS`], is not narrow.
pub fn narrow_nodes(grid: &Grid, periodic: [bool; 3], open: &[bool], starts: impl IntoIterator<Item = (usize, usize, i32)>) -> usize {
    let mut narrow = 0;
    for (cell, axis, sign) in starts {
        let Some(n) = grid.cells.get(axis).map(|&n| n as i64) else { continue };
        let mut p = grid.coords(cell);
        let mut run = 0;
        let blocked = loop {
            if !open.get(grid.index(p)).copied().unwrap_or(false) {
                break true;
            }
            run += 1;
            if run >= NARROW_CELLS {
                break false;
            }
            let Some(next) = p.get(axis).map(|&v| v as i64 + i64::from(sign.signum())) else { break false };
            let next = if periodic.get(axis).copied().unwrap_or(false) { next.rem_euclid(n) } else { next };
            if !(0..n).contains(&next) {
                break false;
            }
            if let Some(v) = p.get_mut(axis) {
                *v = next as usize;
            }
        };
        narrow += usize::from(blocked && run < NARROW_CELLS);
    }
    narrow
}

/// The share of the wall area above which setup warns that walls run plain bounce-back.
pub const PLAIN_SHARE_WARNING: f64 = 0.10;

/// The share of the surface area above which setup warns of a wall with no centre inside.
pub const UNSEALED_WARNING: f64 = 0.002;

/// The share above which a wall with no centre inside is refused (a choice, not a measurement).
pub const UNSEALED_LIMIT: f64 = 0.02;

/// Most surface samples [`unsealed_area`] takes.
const MAX_SAMPLES: usize = 40_000_000;

/// Samples a surface at about one per cell face: its area, and the area of the samples whose cell
/// has no solid cell (`solid`) in the 3 x 3 x 3 block around it. A body at least two cells thick
/// always has one next to its surface; a wall with no centre inside it has none, and the flow leaks
/// through it. Triangles are visited in order, each cut into `m x m` equal sub-triangles sampled at
/// their centroids, so the result does not depend on threads. Returns `(unsealed, total)` areas
/// inside the box.
pub fn unsealed_area(surface: &BoundarySurface, grid: &Grid, solid: &[bool]) -> Result<(f64, f64)> {
    if solid.len() != grid.len() {
        return Err(Error::Invalid("the wall layer was given a mask that does not fit the grid".into()));
    }
    let (points, triangles) = (surface.points(), surface.triangles());
    let h = grid.spacing;
    if !(h.is_finite() && h > 0.0) {
        return Err(Error::Invalid("the grid spacing is not usable".into()));
    }
    let mut samples = 0usize;
    let (mut bad, mut total) = (0.0f64, 0.0f64);
    for t in triangles {
        let (Some(a), Some(b), Some(c)) = (points.get(t[0] as usize), points.get(t[1] as usize), points.get(t[2] as usize)) else { continue };
        let sub = |u: [f64; 3], v: [f64; 3]| [0, 1, 2].map(|i| v[i] - u[i]);
        let (e1, e2) = (sub(*a, *b), sub(*a, *c));
        let cross = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        let area = 0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        let edge = [e1, e2, sub(*b, *c)].iter().map(|e| (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt()).fold(0.0, f64::max);
        if !(area.is_finite() && edge.is_finite()) || area <= 0.0 {
            continue;
        }
        let m = (edge / h).ceil().clamp(1.0, 256.0) as usize;
        samples = samples.saturating_add(m * m);
        if samples > MAX_SAMPLES {
            return Err(Error::Limit(format!("the obstacles' surface needs more than {MAX_SAMPLES} samples at spacing {h}")));
        }
        let each = area / (m * m) as f64;
        // Sub-triangles of the m x m split: upright (i, j) and inverted ones.
        for i in 0..m {
            for j in 0..m - i {
                for inverted in [false, true] {
                    if inverted && i + j + 1 >= m {
                        continue;
                    }
                    let (s, r) = if inverted {
                        ((i as f64 + 2.0 / 3.0) / m as f64, (j as f64 + 2.0 / 3.0) / m as f64)
                    } else {
                        ((i as f64 + 1.0 / 3.0) / m as f64, (j as f64 + 1.0 / 3.0) / m as f64)
                    };
                    let p = [0, 1, 2].map(|k| a[k] + s * e1[k] + r * e2[k]);
                    let cell = [0, 1, 2].map(|k| ((p[k] - grid.origin[k]) / h).floor());
                    if cell.iter().zip(grid.cells).any(|(c, n)| !c.is_finite() || *c < 0.0 || *c >= n as f64) {
                        continue;
                    }
                    total += each;
                    let at = cell.map(|c| c as i64);
                    let mut found = false;
                    for dz in -1..=1i64 {
                        for dy in -1..=1i64 {
                            for dx in -1..=1i64 {
                                let q = [at[0] + dx, at[1] + dy, at[2] + dz];
                                if q.iter().zip(grid.cells).all(|(v, n)| *v >= 0 && *v < n as i64) {
                                    found |= solid.get(grid.index([q[0] as usize, q[1] as usize, q[2] as usize])).copied().unwrap_or(false);
                                }
                            }
                        }
                    }
                    if !found {
                        bad += each;
                    }
                }
            }
        }
    }
    Ok((bad, total))
}

/// The cell one step from `p` along `axis` (`up`: the positive way) on a grid that is periodic where
/// `periodic` says so, `None` past a closed face.
fn step(grid: &Grid, periodic: [bool; 3], p: [usize; 3], axis: usize, up: bool, by: usize) -> Option<usize> {
    let len = *grid.cells.get(axis)?;
    let here = *p.get(axis)?;
    let there = if up {
        let t = here + by;
        if t < len {
            t
        } else if periodic[axis] {
            t % len
        } else {
            return None;
        }
    } else if here >= by {
        here - by
    } else if periodic[axis] {
        (here + len - by % len) % len
    } else {
        return None;
    };
    let mut q = p;
    *q.get_mut(axis)? = there;
    Some(grid.index(q))
}

/// Flag the wall nodes of the half-way wall as fluid: `fluid` and `owner` are the voxeliser's masks
/// (a cell is solid when `fluid` is false, and an obstacle's when `owner` is not [`NONE`]); the node
/// cells become fluid with no owner. `closed(p)` is true for a cell that may not be a node (in the
/// layer of an open face of the box).
pub fn promote(
    grid: &Grid,
    periodic: [bool; 3],
    fluid: &mut [bool],
    owner: &mut [u32],
    closed: &dyn Fn([usize; 3]) -> bool,
) -> Result<WallLayerReport> {
    let n = grid.len();
    if fluid.len() != n || owner.len() != n {
        return Err(Error::Invalid("the wall layer was given masks that do not fit the grid".into()));
    }
    let open: Vec<bool> = fluid.to_vec();
    let solid = |c: Option<usize>| c.is_some_and(|c| !open.get(c).copied().unwrap_or(true) && owner.get(c).copied().unwrap_or(NONE) != NONE);
    let is_open = |c: Option<usize>| c.is_some_and(|c| open.get(c).copied().unwrap_or(false));
    let mut report = WallLayerReport::default();
    let mut nodes = filled(n, false, "wall nodes")?;
    for c in 0..n {
        if open.get(c).copied().unwrap_or(true) || owner.get(c).copied().unwrap_or(NONE) == NONE {
            continue;
        }
        let p = grid.coords(c);
        let mut faces = 0;
        let mut node = false;
        for axis in 0..3 {
            for up in [false, true] {
                // `up`: the open cell is on the positive side, the solid behind it on the negative one.
                if !is_open(step(grid, periodic, p, axis, up, 1)) {
                    continue;
                }
                faces += 1;
                let behind = step(grid, periodic, p, axis, !up, 1);
                let beyond = step(grid, periodic, p, axis, !up, 2);
                if solid(behind) && !is_open(beyond) && !closed(p) {
                    node = true;
                }
            }
        }
        report.faces += faces;
        if node {
            report.nodes += 1;
            if let Some(x) = nodes.get_mut(c) {
                *x = true;
            }
        } else {
            report.plain_faces += faces;
        }
    }
    for (c, _) in nodes.iter().enumerate().filter(|(_, x)| **x) {
        if let Some(f) = fluid.get_mut(c) {
            *f = true;
        }
        if let Some(o) = owner.get_mut(c) {
            *o = NONE;
        }
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_passage_narrower_than_the_limit_is_counted() {
        // Open runs of 5 and 14 cells between solid cells, periodic along z.
        let pattern = format!("s{}s{}", "o".repeat(5), "o".repeat(14));
        let grid = Grid { cells: [1, 1, pattern.len()], spacing: 1.0, origin: [0.0; 3] };
        let open: Vec<bool> = pattern.chars().map(|c| c == 'o').collect();
        // From the first open cell of each run, along +z.
        assert_eq!(narrow_nodes(&grid, [true; 3], &open, [(1, 2, 1), (7, 2, 1)]), 1);
        // The box edge ends a run on a non-periodic axis without making it narrow.
        assert_eq!(narrow_nodes(&grid, [true, true, false], &open, [(grid.cells[2] - 1, 2, 1)]), 0);
        // A cell index outside the grid, or a bad axis, is skipped.
        assert_eq!(narrow_nodes(&grid, [true; 3], &open, [(1, 7, 1)]), 0);
    }

    /// A column of `len` cells along z, periodic in x and y, with `pattern` ('o' open, 's' solid).
    fn column(pattern: &str) -> (Grid, Vec<bool>, Vec<u32>) {
        let grid = Grid { cells: [1, 1, pattern.len()], spacing: 1.0, origin: [0.0; 3] };
        let fluid: Vec<bool> = pattern.chars().map(|c| c == 'o').collect();
        let owner: Vec<u32> = pattern.chars().map(|c| if c == 's' { 0 } else { NONE }).collect();
        (grid, fluid, owner)
    }

    fn nodes(pattern: &str) -> (String, WallLayerReport) {
        let (grid, mut fluid, mut owner) = column(pattern);
        let before = fluid.clone();
        let report = promote(&grid, [true, true, false], &mut fluid, &mut owner, &|_| false).expect("promote");
        let out = fluid
            .iter()
            .zip(&before)
            .map(|(now, was)| {
                if *was {
                    'o'
                } else if *now {
                    'n'
                } else {
                    's'
                }
            })
            .collect();
        (out, report)
    }

    #[test]
    fn a_wall_of_three_or_more_centres_gets_a_node_on_each_open_side() {
        let (out, r) = nodes("ooosssooo");
        assert_eq!(out, "ooonsnooo");
        assert_eq!((r.nodes, r.faces, r.plain_faces), (2, 2, 0));
        let (out, r) = nodes("ooossssooo");
        assert_eq!(out, "ooonssnooo");
        assert_eq!(r.plain_faces, 0);
    }

    #[test]
    fn thinner_walls_stay_solid_and_sealed() {
        for pattern in ["ooosooo", "ooossooo"] {
            let (out, r) = nodes(pattern);
            assert_eq!(out, pattern, "{pattern}: no node");
            assert_eq!(r.nodes, 0);
            assert_eq!(r.plain_share(), 1.0);
        }
    }

    #[test]
    fn a_slab_against_the_box_needs_only_two_cells() {
        // The box face behind is a wall: nothing is open beyond the cell behind the node.
        let (out, r) = nodes("sssooo");
        assert_eq!(out, "ssnooo");
        assert_eq!(r.plain_faces, 0);
        let (out, _) = nodes("ssooo");
        assert_eq!(out, "snooo");
        let (out, _) = nodes("sooo");
        assert_eq!(out, "sooo", "a single cell against the box has nothing behind it");
    }

    #[test]
    fn a_cell_in_an_open_face_layer_is_no_node() {
        let (grid, mut fluid, mut owner) = column("ooosssooo");
        let r = promote(&grid, [true, true, false], &mut fluid, &mut owner, &|p| p[2] == 3).expect("promote");
        assert_eq!(r.nodes, 1, "the other side still has its node");
        assert_eq!(r.plain_faces, 1);
    }

    #[test]
    fn masks_that_do_not_fit_are_refused() {
        let (grid, mut fluid, mut owner) = column("oosso");
        fluid.pop();
        assert!(promote(&grid, [true, true, false], &mut fluid, &mut owner, &|_| false).is_err());
    }
}
