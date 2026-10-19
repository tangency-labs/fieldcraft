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
