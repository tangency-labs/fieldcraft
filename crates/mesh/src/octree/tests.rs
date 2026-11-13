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

/// With the exact surfaces given, the boundary nodes lie on the cylinder and on its rims, not on
/// the chords of its 24 facets (which stray up to 0.043 from it).
#[test]
fn boundary_nodes_lie_on_the_exact_surfaces_when_given() {
    use crate::tet::geometry::{AnalyticSurface, Tessellation};
    let (r, h) = (5.0, 8.0);
    let s = cylinder_surface(r, h, 24);
    let analytic = vec![
        Some(AnalyticSurface::Plane { normal: [0.0, 0.0, -1.0], point: [0.0; 3] }),
        Some(AnalyticSurface::Plane { normal: [0.0, 0.0, 1.0], point: [0.0, 0.0, h] }),
        Some(AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: r }),
    ];
    let exact = Tessellation::new(&s, analytic);
    let side = |m: &OctreeHexMesh| -> Vec<f64> {
        let mut nodes: BTreeSet<u32> = BTreeSet::new();
        for b in m.boundary.iter().filter(|b| b.brep_face == 2) {
            nodes.extend(topo::FACES[(b.face - 1) as usize].map(|i| m.hexes[b.element as usize][i]));
        }
        nodes.iter().map(|&v| (m.nodes[v as usize][0].hypot(m.nodes[v as usize][1]) - r).abs()).collect()
    };
    let chords = octree_hex_mesh(&s, &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    assert!(side(&chords).iter().fold(0.0f64, |a, &d| a.max(d)) > 1e-3);
    let m = octree_hex_mesh_exact(&s, Some(&exact), &UniformSize(2.0), &OctreeOptions::default()).unwrap();
    let off = side(&m);
    assert!(off.len() > 50, "{} side nodes", off.len());
    let worst = off.iter().fold(0.0f64, |a, &d| a.max(d));
    assert!(worst < 1e-9 * r, "a side node {worst} off the cylinder");
    assert!(m.report.min_scaled_jacobian > 0.3, "{:?}", m.report);
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

/// A prism along x of length `l` over the right triangle (0, 0), (a, 0), (0, a) in (y, z), its
/// two 45° wedges along x; every side divided into `n` and the length into `m`, five faces.
fn wedge_prism(a: f64, l: f64, n: usize, m: usize) -> BoundarySurface {
    let mut ids: BTreeMap<(usize, usize, usize), u32> = BTreeMap::new();
    let mut pts: Vec<[f64; 3]> = Vec::new();
    let mut id = |i: usize, j: usize, k: usize| -> u32 {
        *ids.entry((i, j, k)).or_insert_with(|| {
            pts.push([l * k as f64 / m as f64, a * i as f64 / n as f64, a * j as f64 / n as f64]);
            (pts.len() - 1) as u32
        })
    };
    let (mut tris, mut tags): (Vec<[u32; 3]>, Vec<u32>) = (Vec::new(), Vec::new());
    // The ends: triangular grids, (i, j) with i + j ≤ n.
    for (k, tag) in [(0, 0), (m, 1)] {
        for i in 0..n {
            for j in 0..n - i {
                tris.push([id(i, j, k), id(i + 1, j, k), id(i, j + 1, k)]);
                tags.push(tag);
                if i + j + 1 < n {
                    tris.push([id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k)]);
                    tags.push(tag);
                }
            }
        }
    }
    // The sides, each a strip of quads along x: z = 0 (j = 0), y = 0 (i = 0), the slope (i + j = n).
    let sides: [(Box<dyn Fn(usize) -> (usize, usize)>, u32); 3] =
        [(Box::new(|s| (s, 0)), 2), (Box::new(|s| (0, s)), 3), (Box::new(move |s| (n - s, s)), 4)];
    for (at, tag) in &sides {
        for s in 0..n {
            for k in 0..m {
                let ((i0, j0), (i1, j1)) = (at(s), at(s + 1));
                let q = [id(i0, j0, k), id(i1, j1, k), id(i1, j1, k + 1), id(i0, j0, k + 1)];
                tris.extend([[q[0], q[1], q[2]], [q[0], q[2], q[3]]]);
                tags.extend([*tag, *tag]);
            }
        }
    }
    // Every triangle facing outward: away from the prism's centroid.
    let centre = [l / 2.0, a / 3.0, a / 3.0];
    for t in &mut tris {
        let [p, q, r] = t.map(|v| pts[v as usize]);
        let nrm = crate::vec3::cross(crate::vec3::sub(q, p), crate::vec3::sub(r, p));
        let mid = [(p[0] + q[0] + r[0]) / 3.0, (p[1] + q[1] + r[1]) / 3.0, (p[2] + q[2] + r[2]) / 3.0];
        if crate::vec3::dot(nrm, crate::vec3::sub(mid, centre)) < 0.0 {
            t.swap(1, 2);
        }
    }
    BoundarySurface::new(&pts, &tris, &tags, &crate::surface::SurfaceOptions::default()).unwrap()
}

/// A wedge's wall is measured no nearer its edge than the local size: the prism's 45° wedges are
/// at least `h tan 45° = h` thick there, however finely the faces are tessellated. Measured at
/// every triangle centre, the fine tessellation's thinnest wall was its triangles' distance from
/// the edge, a twentieth of the coarse one's, and the grid it asked for refined with it (the
/// partitioned bracket's octree, meshed from a finer surface, then failed).
#[test]
fn a_wedges_wall_does_not_depend_on_its_tessellation() {
    let (a, l, h) = (10.0, 40.0, 2.0);
    let options = OctreeOptions::default();
    let mut thinnest = Vec::new();
    for (n, m) in [(1, 1), (20, 40)] {
        let s = wedge_prism(a, l, n, m);
        let g = Geometry::new(&s, options.sharp_angle, options.blade_angle).unwrap();
        let t = triangle_thickness(&g, &UniformSize(h));
        thinnest.push(t.iter().copied().fold(f64::INFINITY, f64::min));
    }
    for (k, w) in thinnest.iter().enumerate() {
        assert!(*w >= h * (1.0 - 1e-9), "tessellation {k}: the thinnest wall is {w}, under {h}: {thinnest:?}");
    }
}
