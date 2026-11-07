use super::*;
use crate::octree::geom::tests::box_surface;
use crate::size::UniformSize;

#[test]
fn a_box_meshes_with_exact_volume_and_tagged_faces() {
    let s = box_surface([0.0; 3], [10.0, 8.0, 6.0]);
    let m = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    let r = &m.report;
    assert!((r.volume - 480.0).abs() < 1e-9 * 480.0, "{r:?}");
    assert!(r.min_scaled_jacobian > 0.3, "{r:?}");
    let faces: BTreeSet<u32> = m.boundary.iter().map(|b| b.brep_face).collect();
    assert_eq!(faces.len(), 6);
    // Each boundary quad lies in the plane of its face.
    for b in &m.boundary {
        let hx = m.hexes[b.element as usize];
        let q = topo::FACES[(b.face - 1) as usize].map(|i| m.nodes[hx[i] as usize]);
        let (axis, value) = match b.brep_face {
            0 => (2, 0.0),
            1 => (2, 6.0),
            2 => (1, 0.0),
            3 => (1, 8.0),
            4 => (0, 0.0),
            _ => (0, 10.0),
        };
        assert!(q.iter().all(|p| (p[axis] - value).abs() < 1e-9), "{b:?} {q:?}");
    }
}

/// A closed cylinder of radius `r` and height `h` on the z axis, `n` segments: faces 0 (bottom
/// cap), 1 (top cap), 2 (side).
pub(crate) fn cylinder_surface(r: f64, h: f64, n: usize) -> BoundarySurface {
    let mut p = vec![[0.0, 0.0, 0.0], [0.0, 0.0, h]];
    for z in [0.0, h] {
        for i in 0..n {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            p.push([r * a.cos(), r * a.sin(), z]);
        }
    }
    let (b, t) = (|i: usize| (2 + i % n) as u32, |i: usize| (2 + n + i % n) as u32);
    let mut tris = Vec::new();
    let mut face = Vec::new();
    for i in 0..n {
        tris.push([0, b(i + 1), b(i)]);
        face.push(0);
        tris.push([1, t(i), t(i + 1)]);
        face.push(1);
        tris.push([b(i), b(i + 1), t(i + 1)]);
        tris.push([b(i), t(i + 1), t(i)]);
        face.extend([2, 2]);
    }
    BoundarySurface::new(&p, &tris, &face, &crate::surface::SurfaceOptions::default()).unwrap()
}

/// A UV sphere of radius `r`, one face.
pub(crate) fn sphere_surface(r: f64, n: usize) -> BoundarySurface {
    let rings = n / 2;
    let mut p = vec![[0.0, 0.0, -r], [0.0, 0.0, r]];
    for j in 1..rings {
        let th = std::f64::consts::PI * j as f64 / rings as f64 - std::f64::consts::FRAC_PI_2;
        for i in 0..n {
            let a = std::f64::consts::TAU * i as f64 / n as f64;
            p.push([r * th.cos() * a.cos(), r * th.cos() * a.sin(), r * th.sin()]);
        }
    }
    let id = |j: usize, i: usize| (2 + (j - 1) * n + i % n) as u32;
    let mut tris = Vec::new();
    for i in 0..n {
        tris.push([0, id(1, i + 1), id(1, i)]);
        tris.push([1, id(rings - 1, i), id(rings - 1, i + 1)]);
    }
    for j in 1..rings - 1 {
        for i in 0..n {
            tris.push([id(j, i), id(j, i + 1), id(j + 1, i + 1)]);
            tris.push([id(j, i), id(j + 1, i + 1), id(j + 1, i)]);
        }
    }
    let face = vec![0; tris.len()];
    BoundarySurface::new(&p, &tris, &face, &crate::surface::SurfaceOptions::default()).unwrap()
}

#[test]
fn a_cylinder_meshes_with_its_rims_captured() {
    let s = cylinder_surface(5.0, 8.0, 64);
    let m = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    let r = &m.report;
    assert!(((r.volume - r.surface_volume) / r.surface_volume).abs() < 0.01, "{r:?}");
    let faces: BTreeSet<u32> = m.boundary.iter().map(|b| b.brep_face).collect();
    assert_eq!(faces.len(), 3);
}

#[test]
fn a_sphere_meshes() {
    let s = sphere_surface(5.0, 48);
    let m = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    let r = &m.report;
    assert!(((r.volume - r.surface_volume) / r.surface_volume).abs() < 0.01, "{r:?}");
}

#[test]
fn a_slab_gets_the_levels_its_thickness_needs_up_front() {
    // 1 mm thick: three cells across need cells of 1/3, two levels below a size of 3.
    let s = box_surface([0.0; 3], [12.0, 9.0, 1.0]);
    let m = octree_hex_mesh(&s, &UniformSize(3.0), &OctreeOptions::default()).unwrap();
    let r = &m.report;
    assert_eq!(r.thinnest_wall, 1.0, "{r:?}");
    assert_eq!(r.root_size, 3.0, "{r:?}");
    assert_eq!(r.levels, 2, "{r:?}");
    assert_eq!(r.attempts, 1, "{r:?}");
    assert!((r.volume - 108.0).abs() < 1e-9 * 108.0, "{r:?}");
}

#[test]
fn walls_beyond_the_limits_fail_before_anything_is_built() {
    let s = box_surface([0.0; 3], [12.0, 9.0, 1.0]);
    // Too few levels: the requested size needs 2, the shrunk root 1.
    let shallow = OctreeOptions { max_depth: 0, curvature_levels: 0, ..OctreeOptions::default() };
    let e = octree_hex_mesh(&s, &UniformSize(3.0), &shallow).unwrap_err().to_string();
    assert!(e.contains("walls 1 thick") && e.contains("levels"), "{e}");
    // Too few cells: refused with the count, at once.
    let small = OctreeOptions { max_elements: 50, ..OctreeOptions::default() };
    let start = std::time::Instant::now();
    let e = octree_hex_mesh(&s, &UniformSize(3.0), &small).unwrap_err().to_string();
    assert!(e.contains("walls 1 thick") && e.contains("more than"), "{e}");
    assert!(start.elapsed().as_secs_f64() < 5.0);
}

#[test]
fn the_root_shrinks_when_the_requested_size_would_need_too_many_cells() {
    // A size of 8 needs 3 levels for a 1 mm wall (cells of 0.30); a root of 3 needs 2 levels
    // with cells of exactly 1/3. With 2 levels allowed only the shrunk root fits.
    let s = box_surface([0.0; 3], [12.0, 9.0, 1.0]);
    let o = OctreeOptions { max_depth: 2, ..OctreeOptions::default() };
    let m = octree_hex_mesh(&s, &UniformSize(8.0), &o).unwrap();
    assert!((m.report.root_size - 3.0).abs() < 1e-12, "{:?}", m.report);
    assert_eq!(m.report.levels, 2);
}

#[test]
fn hostile_options_are_errors() {
    let s = box_surface([0.0; 3], [12.0, 9.0, 1.0]);
    for o in [
        OctreeOptions { max_depth: 9, ..OctreeOptions::default() },
        OctreeOptions { max_elements: 0, ..OctreeOptions::default() },
        OctreeOptions { max_elements: 1, ..OctreeOptions::default() },
        OctreeOptions { thickness_cells: f64::NAN, ..OctreeOptions::default() },
    ] {
        assert!(octree_hex_mesh(&s, &UniformSize(3.0), &o).is_err(), "{o:?}");
    }
    // The deepest grid allowed on a size far below the part: refused by the cell limit.
    let deep = OctreeOptions { max_depth: 5, max_elements: 10_000, ..OctreeOptions::default() };
    assert!(octree_hex_mesh(&s, &UniformSize(0.05), &deep).is_err());
}

#[test]
fn meshing_is_deterministic() {
    let s = cylinder_surface(5.0, 8.0, 64);
    let a = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    let b = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    assert_eq!(a, b);
}

/// A size box refines the grid inside it to the level its size needs (closure may carry that
/// level further), and the root cells keep the global size; a size deeper than `max_depth` allows is refused before meshing.
#[test]
fn a_size_box_refines_the_grid_inside_it() {
    use crate::size::{LocalSize, Region, RegionSize};
    let s = box_surface([0.0; 3], [18.0, 18.0, 18.0]);
    let local = LocalSize { region: Region::Box { lo: [0.0; 3], hi: [6.0; 3] }, size: 2.0 };
    let field = RegionSize::new(6.0, 0.3, vec![local.clone()]).unwrap();
    let m = octree_hex_mesh(&s, &field, &OctreeOptions::default()).unwrap();
    let edge = |h: &[u32; 8]| {
        let (a, b) = (m.nodes[h[0] as usize], m.nodes[h[1] as usize]);
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    };
    let centre = |h: &[u32; 8]| [0, 1, 2].map(|k| h.iter().map(|&v| m.nodes[v as usize][k]).sum::<f64>() / 8.0);
    let inside: Vec<f64> = m.hexes.iter().filter(|h| centre(h).iter().all(|&c| (1.0..5.0).contains(&c))).map(edge).collect();
    assert!(!inside.is_empty() && inside.iter().all(|&e| e < 2.0 * 1.2), "{inside:?}");
    assert_eq!((m.report.root_size, m.report.levels), (6.0, 1), "{:?}", m.report);
    let shallow = OctreeOptions { max_depth: 1, curvature_levels: 1, ..OctreeOptions::default() };
    let deep = RegionSize::new(6.0, 0.3, vec![LocalSize { size: 0.5, ..local }]).unwrap();
    let e = octree_hex_mesh(&s, &deep, &shallow).unwrap_err().to_string();
    assert!(e.contains("max_depth"), "{e}");
}

/// The thinnest boundary element of `m` against its face on the boundary: the distance from the
/// opposite face's centre to the boundary face's plane over the boundary face's mean edge.
fn thinnest_layer(m: &OctreeHexMesh) -> f64 {
    let mut thinnest = f64::INFINITY;
    for b in &m.boundary {
        let hx = m.hexes[b.element as usize];
        let f = (b.face - 1) as usize;
        let q = topo::FACES[f].map(|i| m.nodes[hx[i] as usize]);
        let o = topo::FACES[[1, 0, 4, 5, 2, 3][f]].map(|i| m.nodes[hx[i] as usize]);
        let n = cross(sub(q[2], q[0]), sub(q[3], q[1]));
        let n = n.map(|x| x / norm(n));
        let c = |p: [P; 4]| [0, 1, 2].map(|k| p.iter().map(|x| x[k]).sum::<f64>() / 4.0);
        let height = dot(sub(c(o), c(q)), n).abs();
        let edge = (0..4).map(|i| norm(sub(q[(i + 1) % 4], q[i]))).sum::<f64>() / 4.0;
        thinnest = thinnest.min(height / edge);
    }
    thinnest
}

#[test]
fn the_core_recedes_to_half_a_cell_from_the_boundary() {
    // A block of 3 × 3 × 3 unit cells 0.3 of a cell inside a box on every side: the padding and
    // the pillow would share a 0.3 gap. The faces recede to half a cell; edge and corner
    // vertices move along their diagonal normals, by the same step, and the core's elements stay
    // well shaped.
    let s = box_surface([0.0; 3], [3.6; 3]);
    let g = Geometry::new(&s, 30.0, 160.0).unwrap();
    let (points, hexes) = topo::tests::block(3);
    let points: Vec<P> = points.iter().map(|p| p.map(|x| x + 0.3)).collect();
    let mut dist: Vec<f64> = points.iter().map(|&p| g.nearest(p, |_| true).unwrap().2.sqrt()).collect();
    let mut w = Work { size: vec![1.0; hexes.len()], role: Vec::new(), core_count: hexes.len(), points, hexes };
    let q: Vec<[u32; 4]> = topo::boundary(&w.hexes).iter().map(|x| x.1).collect();
    let boundary: BTreeSet<u32> = q.iter().flatten().copied().collect();
    assert!(boundary.iter().all(|&v| (dist[v as usize] - 0.3).abs() < 1e-12));
    // Vertices inside a face: one coordinate on the core's boundary.
    let faces: Vec<u32> = boundary
        .iter()
        .copied()
        .filter(|&v| w.points[v as usize].iter().filter(|&&x| (x - 0.3).abs() < 1e-12 || (x - 3.3).abs() < 1e-12).count() == 1)
        .collect();
    recede(&g, &mut w, &q, &mut dist);
    for &v in &boundary {
        let d = g.nearest(w.points[v as usize], |_| true).unwrap().2.sqrt();
        assert!((d - dist[v as usize]).abs() < 1e-12 && d > 0.4, "vertex {v} at {:?}: {d}", w.points[v as usize]);
    }
    assert_eq!(faces.len(), 24);
    assert!(faces.iter().all(|&v| dist[v as usize] > 0.5 - 1e-9), "{dist:?}");
    let worst = w.hexes.iter().map(|hx| smooth::min_scaled_jacobian(&hx.map(|v| w.points[v as usize]))).fold(f64::INFINITY, f64::min);
    assert!(worst > 0.9, "{worst}");
    // On a block, every boundary vertex can still be pushed outward.
    let pts = &w.points;
    assert!(topo::hidden_vertices(&w.hexes, &vec![true; w.hexes.len()], &|v| pts[v as usize], BUFFER_MARGIN).is_empty());
}

#[test]
fn a_core_one_cell_thick_keeps_half_its_thickness() {
    // A slab of 3 × 3 × 1 unit cells 0.3 inside a box, its elements sized 1.4 (as a transition
    // element's edges are shorter than its cell): each face would recede 0.4, which leaves a
    // slab 0.2 thick whose elements are still perfect boxes, so only the edge lengths stop it.
    let s = box_surface([-0.3, -0.3, 0.7], [3.3, 3.3, 2.3]);
    let g = Geometry::new(&s, 30.0, 160.0).unwrap();
    let (points, hexes) = topo::tests::block(3);
    let hexes: Vec<[u32; 8]> = hexes.into_iter().filter(|hx| hx.iter().all(|&v| (1.0..=2.0).contains(&points[v as usize][2]))).collect();
    assert_eq!(hexes.len(), 9);
    let mut dist: Vec<f64> = points.iter().map(|&p| g.nearest(p, |_| true).unwrap().2.sqrt()).collect();
    let mut w = Work { size: vec![1.4; hexes.len()], role: Vec::new(), core_count: hexes.len(), points, hexes };
    let q: Vec<[u32; 4]> = topo::boundary(&w.hexes).iter().map(|x| x.1).collect();
    let start = w.points.clone();
    recede(&g, &mut w, &q, &mut dist);
    // The edges across the slab, as long as it is thick.
    let across: Vec<f64> = w
        .hexes
        .iter()
        .flat_map(|hx| crate::quality::jacobian::HEX_EDGES.iter().map(move |&(a, b)| (hx[a] as usize, hx[b] as usize)))
        .filter(|&(a, b)| (start[a][2] - start[b][2]).abs() > 0.5)
        .map(|(a, b)| norm(sub(w.points[a], w.points[b])))
        .collect();
    assert_eq!(across.len(), 9 * 4);
    let thinnest = across.iter().copied().fold(f64::INFINITY, f64::min);
    assert!((0.5 - 1e-9..0.7).contains(&thinnest), "{thinnest}");
}

/// A stand-in for a build's mesh: only its smallest scaled Jacobian matters to [`rebuilt`].
fn mesh_with(min_scaled_jacobian: f64) -> Box<OctreeHexMesh> {
    let report = OctreeReport {
        attempts: 0,
        grid_cells: 0,
        root_size: 1.0,
        levels: 0,
        thinnest_wall: f64::INFINITY,
        narrowest_patch: f64::INFINITY,
        refined_cells: Vec::new(),
        closed_cells: Vec::new(),
        core_elements: 0,
        layer_elements: 0,
        frame: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        corners: 0,
        curves: 0,
        patches: 0,
        sweeps: 0,
        min_scaled_jacobian,
        mean_scaled_jacobian: min_scaled_jacobian,
        volume: 0.0,
        surface_volume: 0.0,
        boundary_distance: 0.0,
        retries: Vec::new(),
    };
    Box::new(OctreeHexMesh { nodes: Vec::new(), hexes: Vec::new(), boundary: Vec::new(), report })
}

/// What [`rebuilt`] returns when its builds return `builds` in turn, with four attempts: the
/// kept mesh's smallest scaled Jacobian, builds and reasons, or the error.
fn rebuild_with(builds: Vec<Result<Attempt>>) -> std::result::Result<(f64, usize, Vec<String>), Error> {
    let mut builds = std::collections::VecDeque::from(builds);
    let mut build = |_: &[P], _: usize| builds.pop_front().unwrap_or_else(|| Err(Error::Failed("no more builds".into())));
    let m = rebuilt(4, &|attempt| 10 * (attempt + 1), &mut build)?;
    Ok((m.report.min_scaled_jacobian, m.report.attempts, m.report.retries.iter().map(|r| r.to_string()).collect()))
}

#[test]
fn a_rebuild_is_kept_only_while_it_raises_the_worst_element() {
    let done = |sj: f64| Ok(Attempt::Done(mesh_with(sj), Vec::new()));
    let poor = |sj: f64| Ok(Attempt::Done(mesh_with(sj), vec![[0.0; 3]]));
    let retry = || Ok(Attempt::Retry(Retry { at: vec![[0.0; 3]], why: message!["2 elements stay inverted"] }));
    let failed = || Err(Error::Failed("the pad layer folds".into()));
    let limit = || Err(Error::Limit("the grid needs 99 cells".into()));
    let cancelled = || Err(Error::Model(fieldcraft_model::Error::Cancelled));

    // A mesh without poor elements is returned at once.
    assert_eq!(rebuild_with(vec![done(0.5)]).unwrap(), (0.5, 1, vec![]));
    // A rebuild that raises the worst element is kept, one that does not returns the earlier mesh.
    let (sj, builds, why) = rebuild_with(vec![poor(0.05), done(0.3)]).unwrap();
    assert_eq!((sj, builds, why.len()), (0.3, 2, 1));
    assert!(why[0].contains("1 elements have a scaled Jacobian under 0.2"), "{why:?}");
    let (sj, builds, _) = rebuild_with(vec![poor(0.05), poor(0.04)]).unwrap();
    assert_eq!((sj, builds), (0.05, 2));
    // Rebuilds continue while each raises it, up to the last attempt, whose mesh is returned.
    let (sj, builds, why) = rebuild_with(vec![poor(0.05), poor(0.06), poor(0.07), poor(0.08)]).unwrap();
    assert_eq!((sj, builds, why.len()), (0.08, 4, 3));

    // Once a valid mesh is in hand, a rebuild that fails keeps it and says why: over the grid
    // cap, any other error, or after a failure it retried.
    let (sj, builds, why) = rebuild_with(vec![poor(0.05), limit()]).unwrap();
    assert_eq!((sj, builds), (0.05, 2));
    assert!(why[1].contains("needs more than 20 grid cells") && why[1].contains("99 cells"), "{why:?}");
    let (sj, builds, why) = rebuild_with(vec![poor(0.05), failed()]).unwrap();
    assert_eq!((sj, builds), (0.05, 2));
    assert!(why[1].contains("the rebuild failed") && why[1].contains("the pad layer folds"), "{why:?}");
    let (sj, builds, why) = rebuild_with(vec![poor(0.05), retry(), failed()]).unwrap();
    assert_eq!((sj, builds, why.len()), (0.05, 3, 3));
    assert!(why[1].contains("inverted"), "{why:?}");
    let (sj, builds, _) = rebuild_with(vec![poor(0.05), retry(), retry(), retry()]).unwrap();
    assert_eq!((sj, builds), (0.05, 4));
    // Cancellation is never swallowed.
    assert!(matches!(rebuild_with(vec![poor(0.05), cancelled()]), Err(Error::Model(fieldcraft_model::Error::Cancelled))));

    // Without a valid mesh: a failure propagates, a grid over its cap ends the retries, and
    // retries that never mesh fail with every reason.
    assert!(matches!(rebuild_with(vec![failed()]), Err(Error::Failed(_))));
    assert!(matches!(rebuild_with(vec![retry(), failed()]), Err(Error::Failed(m)) if m == "the pad layer folds"));
    let e = rebuild_with(vec![retry(), limit()]).unwrap_err().to_string();
    assert!(e.contains("failed after 2 builds") && e.contains("needs more than 20 grid cells"), "{e}");
    let e = rebuild_with(vec![retry(), retry(), retry(), retry()]).unwrap_err().to_string();
    assert!(e.contains("failed after 4 builds") && e.matches("inverted").count() == 4, "{e}");
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn boxes_whose_core_is_a_fraction_of_a_cell_from_a_face_keep_their_layers() {
    // At a size of 2, a 7 or 7.5 high box leaves its core a fraction of a cell from its top and
    // bottom faces. Without the core receding, the 7 box's boundary elements came out 0.012 of
    // their width thick (scaled Jacobian 1e-4) and the 7.5 box failed after four builds.
    for height in [7.0, 7.5] {
        let s = box_surface([0.0; 3], [10.0, 8.0, height]);
        let m = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
        let r = &m.report;
        assert!(thinnest_layer(&m) > 0.4, "{height}: {} {r:?}", thinnest_layer(&m));
        assert!(r.min_scaled_jacobian > 0.3, "{height}: {r:?}");
        assert!((r.volume - 80.0 * height).abs() < 1e-9 * 80.0 * height, "{height}: {r:?}");
    }
}
