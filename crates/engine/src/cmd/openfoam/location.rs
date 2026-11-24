//! A point in the fluid for `snappyHexMesh`'s `locationInMesh`: the side of the obstacles the case
//! keeps. Inside and outside come from the lattice voxeliser ([`fieldcraft_lbm::voxel`]), which
//! decides cell centres with exact predicates, the same test that marks the obstacles' solid
//! cells, so no second, approximate test is used. The point is the centre of a finest cell in a
//! lattice cell of the fluid that is a whole spacing from every wall and box face, reached by a
//! breadth-first flood fill from the open faces (from the box's faces when none is open), so a closed pocket is never chosen.

use std::collections::VecDeque;

use fieldcraft_decks::openfoam::{MeshOptions, cell_sizes, finest_cell_of, finest_centre};
use fieldcraft_lbm::voxel::{Grid, voxelise};
use fieldcraft_lbm::{Limits, Obstacle};
use fieldcraft_model::{BoxFace, FaceCondition, FluidDomain, progress};

use crate::{EngineError, Result};

/// A point in the fluid of `domain` around `obstacles`, at a finest cell's centre.
pub(super) fn find(domain: &FluidDomain, obstacles: &[Obstacle], mesh: &MeshOptions) -> Result<[f64; 3]> {
    let cells = domain.cells()?;
    let total = cells.iter().try_fold(1u64, |acc, &n| acc.checked_mul(n));
    match total {
        Some(t) if t <= Limits::MAX_CELLS => {}
        _ => {
            let shown = total.map_or_else(|| "more than u64::MAX".to_string(), |t| t.to_string());
            return Err(EngineError::invalid(
                "location",
                format!("a lattice of {shown} cells is above {} to search for a point in the fluid: give `location`", Limits::MAX_CELLS),
            ));
        }
    }
    let cells = cells.map(|n| n as usize);
    let grid = Grid { origin: domain.min, spacing: domain.spacing, cells };
    let solid = solid_cells(&grid, obstacles)?;
    let (_, finest) = cell_sizes(domain, mesh)?;
    // Finest cells per lattice spacing: a power of two, at least one.
    let per_spacing = (domain.spacing / finest).round() as u64;

    let mut visited = vec![false; grid.len()];
    let mut queue = VecDeque::new();
    for seed in seeds(domain, &grid, cells, &solid) {
        if first_visit(&mut visited, seed) {
            queue.push_back(seed);
        }
    }
    while let Some(c) = queue.pop_front() {
        let ijk = grid.coords(c);
        if clear(&grid, &solid, ijk)
            && let Some(p) = point_in_fluid(domain, obstacles, mesh, ijk, per_spacing, finest)?
        {
            return Ok(p);
        }
        for n in neighbours(cells, ijk) {
            let q = grid.index(n);
            if solid.get(q) == Some(&false) && first_visit(&mut visited, q) {
                queue.push_back(q);
            }
        }
    }
    Err(EngineError::invalid("location", "no point in the fluid a whole spacing away from every wall: give `location`"))
}

/// Which cells of the lattice are inside an obstacle, decided by the voxeliser.
fn solid_cells(grid: &Grid, obstacles: &[Obstacle]) -> Result<Vec<bool>> {
    let mut solid = vec![false; grid.len()];
    for o in obstacles {
        progress::check()?;
        voxelise(&o.surface, grid, &mut solid)?;
    }
    Ok(solid)
}

/// Whether a box face admits fluid: velocity and pressure faces are open, the others walls.
fn is_open(condition: &FaceCondition) -> bool {
    matches!(condition, FaceCondition::Velocity { .. } | FaceCondition::Pressure { .. })
}

/// The fluid cells the search starts from: the layer next to each open face, faces in
/// [`BoxFace::ALL`] order and cells by increasing index; with no open face, the layer next to each
/// of the six box faces, so the fluid touching the box and never a closed cavity inside an obstacle.
fn seeds(domain: &FluidDomain, grid: &Grid, cells: [usize; 3], solid: &[bool]) -> Vec<usize> {
    let fluid = |c: usize| solid.get(c) == Some(&false);
    let mut open: Vec<BoxFace> = BoxFace::ALL.into_iter().filter(|&f| is_open(domain.faces.get(f))).collect();
    if open.is_empty() {
        open = BoxFace::ALL.to_vec();
    }
    let mut out = Vec::new();
    for face in open {
        let axis = face.axis();
        let len = cells.get(axis).copied().unwrap_or(0);
        let layer = if face.is_max() { len.saturating_sub(1) } else { 0 };
        out.extend((0..grid.len()).filter(|&c| grid.coords(c).get(axis) == Some(&layer) && fluid(c)));
    }
    out
}

/// Marks `c` visited; true when it was not visited before.
fn first_visit(visited: &mut [bool], c: usize) -> bool {
    match visited.get_mut(c) {
        Some(v) if !*v => {
            *v = true;
            true
        }
        _ => false,
    }
}

/// Whether cell `[x, y, z]` and its 26 neighbours all lie in the lattice and are fluid: the
/// point is then a whole spacing from every wall and box face.
fn clear(grid: &Grid, solid: &[bool], [x, y, z]: [usize; 3]) -> bool {
    let [nx, ny, nz] = grid.cells;
    if x == 0 || y == 0 || z == 0 || x + 1 >= nx || y + 1 >= ny || z + 1 >= nz {
        return false;
    }
    (0..3).all(|dk| {
        (0..3).all(|dj| {
            (0..3).all(|di| {
                let c = grid.index([x + di - 1, y + dj - 1, z + dk - 1]);
                solid.get(c) == Some(&false)
            })
        })
    })
}

/// The lattice cell reached from `ijk` by one step along `axis`, up or down, inside `cells`.
fn step(cells: [usize; 3], ijk: [usize; 3], axis: usize, up: bool) -> Option<[usize; 3]> {
    let v = *ijk.get(axis)?;
    let w = if up { v.checked_add(1)? } else { v.checked_sub(1)? };
    if w >= *cells.get(axis)? {
        return None;
    }
    let mut out = ijk;
    *out.get_mut(axis)? = w;
    Some(out)
}

/// The face neighbours of `ijk` in the order −x, +x, −y, +y, −z, +z.
fn neighbours(cells: [usize; 3], ijk: [usize; 3]) -> impl Iterator<Item = [usize; 3]> {
    (0..3usize).flat_map(move |axis| [false, true].into_iter().filter_map(move |up| step(cells, ijk, axis, up)))
}

/// A point the user gave, as the case writes it: the centre of the finest cell holding it (a
/// point on no cell face of any level), which the voxeliser must find outside every obstacle.
pub(super) fn given(domain: &FluidDomain, obstacles: &[Obstacle], mesh: &MeshOptions, p: [f64; 3]) -> Result<[f64; 3]> {
    let (_, finest) = cell_sizes(domain, mesh)?;
    let centre = finest_centre(domain, mesh, finest_cell_of(domain, mesh, p)?)?;
    outside_all(obstacles, centre, finest)?.ok_or_else(|| {
        EngineError::invalid("location", "the centre of the finest cell holding it lies inside an obstacle: give a point in the fluid")
    })
}

/// `p`, as the voxeliser places it, when it is outside every obstacle: a one-cell grid centred on
/// it, decided by the same exact test as the lattice's cells.
fn outside_all(obstacles: &[Obstacle], p: [f64; 3], finest: f64) -> Result<Option<[f64; 3]>> {
    progress::check()?;
    let one = Grid { origin: p.map(|v| v - finest / 2.0), spacing: finest, cells: [1, 1, 1] };
    let mut inside = [false];
    for o in obstacles {
        voxelise(&o.surface, &one, &mut inside)?;
    }
    let [hit] = inside;
    Ok((!hit).then(|| one.centre(0)))
}

/// The point of an accepted lattice cell: the centre of its finest cell nearest the cell's centre,
/// when the voxeliser finds that point outside every obstacle. `None` when it is inside one.
fn point_in_fluid(
    domain: &FluidDomain,
    obstacles: &[Obstacle],
    mesh: &MeshOptions,
    ijk: [usize; 3],
    per_spacing: u64,
    finest: f64,
) -> Result<Option<[f64; 3]>> {
    let k = ijk.map(|v| v as u64 * per_spacing + per_spacing / 2);
    outside_all(obstacles, finest_centre(domain, mesh, k)?, finest)
}

#[cfg(test)]
mod tests {
    use fieldcraft_decks::openfoam::{RefineBox, cell_sizes};
    use fieldcraft_lbm::Obstacle;
    use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
    use fieldcraft_model::FluidDomain;
    use serde_json::{Value, json};

    use super::*;

    /// A closed box from `lo` to `hi`: eight corners, twelve outward triangles, one tag per face.
    fn cube_obstacle(lo: [f64; 3], hi: [f64; 3]) -> Obstacle {
        let (positions, triangles) = cube_mesh(lo, hi);
        let tri_face: Vec<u32> = (0..triangles.len() as u32).map(|t| t / 2).collect();
        let surface = BoundarySurface::new(&positions, &triangles, &tri_face, &SurfaceOptions::default()).unwrap();
        Obstacle { body: "Cube".into(), surface, faces: Vec::new(), exact: Vec::new() }
    }

    /// The corners and outward triangles of the box from `lo` to `hi`.
    fn cube_mesh(lo: [f64; 3], hi: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        // Corner k takes `hi` on the axes flagged in its row: the order the quads below expect.
        const HIGH: [[usize; 3]; 8] = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        let positions: Vec<[f64; 3]> = HIGH.iter().map(|flags| std::array::from_fn(|a| if flags[a] == 1 { hi[a] } else { lo[a] })).collect();
        let quads: [[u32; 4]; 6] = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]];
        let mut triangles: Vec<[u32; 3]> = Vec::new();
        for q in &quads {
            triangles.push([q[0], q[1], q[2]]);
            triangles.push([q[0], q[2], q[3]]);
        }
        (positions, triangles)
    }

    /// A 3D box from the origin to `max` with lattice `spacing`; `faces` (when given) replaces the
    /// default walls.
    fn domain(spacing: f64, max: [f64; 3], faces: Option<Value>) -> FluidDomain {
        let mut value = json!({
            "material": "Air", "min": [0.0, 0.0, 0.0], "max": max, "spacing": spacing, "obstacles": ["Cube"]
        });
        if let Some(faces) = faces {
            value["faces"] = faces;
        }
        serde_json::from_value(value).unwrap()
    }

    /// Velocity in at the x-min face, pressure out at the x-max face.
    fn channel_faces() -> Option<Value> {
        Some(json!({"x_min": {"type": "velocity", "velocity": [1, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}))
    }

    #[test]
    fn the_point_is_in_the_fluid_and_off_every_face() {
        let max = [0.01, 0.004, 0.004];
        let (lo, hi) = ([0.004, 0.001, 0.001], [0.006, 0.003, 0.003]);
        let domain = domain(0.0005, max, channel_faces());
        let obstacles = [cube_obstacle(lo, hi)];
        let refined = MeshOptions { cell_factor: 2, refine_boxes: vec![RefineBox { min: [0.0; 3], max, spacing: Some(0.000_125) }] };
        for mesh in [MeshOptions::default(), refined] {
            let (_, finest) = cell_sizes(&domain, &mesh).unwrap();
            let p = find(&domain, &obstacles, &mesh).unwrap();
            for a in 0..3 {
                assert!(p[a] >= domain.min[a] && p[a] <= max[a], "{p:?} is outside the box");
                let whole = (p[a] - domain.min[a]) / finest - 0.5;
                assert!((whole - whole.round()).abs() < 1e-9, "{p:?} is not a finest cell's centre");
            }
            let in_cube = (0..3).all(|a| p[a] >= lo[a] && p[a] <= hi[a]);
            assert!(!in_cube, "{p:?} is inside the cube");
        }
    }

    /// A hollow shell: the outer box with outward triangles and the inner box with its triangles
    /// reversed, in one surface (its cavity has winding number 0, so it reads as fluid).
    fn hollow_obstacle(outer: ([f64; 3], [f64; 3]), inner: ([f64; 3], [f64; 3])) -> Obstacle {
        let ((pa, ta), (pb, tb)) = (cube_mesh(outer.0, outer.1), cube_mesh(inner.0, inner.1));
        let shift = pa.len() as u32;
        let mut positions = pa;
        positions.extend(pb);
        let mut triangles = ta;
        triangles.extend(tb.iter().map(|t| [t[0] + shift, t[2] + shift, t[1] + shift]));
        let tri_face = vec![0; triangles.len()];
        let surface = BoundarySurface::new(&positions, &triangles, &tri_face, &SurfaceOptions::default()).unwrap();
        Obstacle { body: "Shell".into(), surface, faces: Vec::new(), exact: Vec::new() }
    }

    #[test]
    fn a_closed_pocket_is_never_chosen() {
        // No open face: the search starts at the box's faces, so the cavity of a hollow obstacle
        // (clear cells before any above the shell) is not chosen; the answer is the same every call.
        let domain = domain(0.0005, [0.004, 0.004, 0.01], None);
        let obstacles = [hollow_obstacle(([-0.001, -0.001, -0.001], [0.005, 0.005, 0.004]), ([0.001; 3], [0.003; 3]))];
        let first = find(&domain, &obstacles, &MeshOptions::default()).unwrap();
        let second = find(&domain, &obstacles, &MeshOptions::default()).unwrap();
        assert_eq!(first, second);
        assert!(first[2] > 0.004, "{first:?} is in the cavity");
    }

    #[test]
    fn a_box_full_of_obstacle_is_refused() {
        let domain = domain(0.0005, [0.01, 0.004, 0.004], channel_faces());
        let obstacles = [cube_obstacle([-0.001; 3], [0.011, 0.005, 0.005])];
        let err = find(&domain, &obstacles, &MeshOptions::default()).unwrap_err();
        assert!(err.to_string().contains("give `location`"), "{err}");
    }

    #[test]
    fn a_given_point_is_its_finest_cells_centre_and_must_be_in_the_fluid() {
        let max = [0.01, 0.004, 0.004];
        let domain = domain(0.0005, max, channel_faces());
        let obstacles = [cube_obstacle([0.004, 0.001, 0.001], [0.006, 0.003, 0.003])];
        let mesh = MeshOptions::default();
        // On the face x = 2 mm between two cells (0.5 mm): moved into the upper one's centre.
        let p = given(&domain, &obstacles, &mesh, [0.002, 0.0021, 0.0021]).unwrap();
        let expected = [0.00225, 0.00225, 0.00225];
        for a in 0..3 {
            assert!((p[a] - expected[a]).abs() <= 1e-15, "{p:?} is not {expected:?}");
        }
        // Outside the box: refused, not clamped onto it.
        let err = given(&domain, &obstacles, &mesh, [100.0, 100.0, 100.0]).unwrap_err();
        assert!(err.to_string().contains("not inside the domain box"), "{err}");
        // Inside the cube: refused, never written.
        let err = given(&domain, &obstacles, &mesh, [0.005, 0.002, 0.002]).unwrap_err();
        assert!(err.to_string().contains("inside an obstacle"), "{err}");
    }

    #[test]
    fn too_large_a_lattice_is_refused() {
        let domain = domain(1e-6, [1.0, 1.0, 1.0], channel_faces());
        let err = find(&domain, &[], &MeshOptions::default()).unwrap_err();
        assert!(err.to_string().contains("give `location`"), "{err}");
    }
}
