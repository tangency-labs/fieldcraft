use std::f64::consts::PI;

use super::intersect::{segment_meets_triangle, surface_pair_intersects, triangles_meet};
use super::*;

type Mesh = (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<u32>);

/// An axis-aligned cube, counter-clockwise from outside, one face tag per side.
fn cube(lo: [f64; 3], size: f64) -> Mesh {
    let p =
        (0..8).map(|i| [lo[0] + size * f64::from(i & 1), lo[1] + size * f64::from((i >> 1) & 1), lo[2] + size * f64::from((i >> 2) & 1)]).collect();
    let t = vec![[0, 2, 3], [0, 3, 1], [4, 5, 7], [4, 7, 6], [0, 1, 5], [0, 5, 4], [2, 6, 7], [2, 7, 3], [0, 4, 6], [0, 6, 2], [1, 3, 7], [1, 7, 5]];
    (p, t, vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5])
}

fn flip(m: &mut Mesh) {
    for t in &mut m.1 {
        t.swap(1, 2);
    }
}

fn merge(a: Mesh, b: Mesh) -> Mesh {
    let (mut p, mut t, mut f) = a;
    let n = p.len() as u32;
    let faces = f.iter().max().map_or(0, |m| m + 1);
    p.extend(b.0);
    t.extend(b.1.iter().map(|t| t.map(|v| v + n)));
    f.extend(b.2.iter().map(|x| x + faces));
    (p, t, f)
}

/// Every triangle with its own three vertices, as SolidCraft's tessellator emits them.
fn soup(m: &Mesh) -> Mesh {
    let mut p = Vec::new();
    let mut t = Vec::new();
    for tri in &m.1 {
        let base = p.len() as u32;
        p.extend(tri.map(|v| m.0[v as usize]));
        t.push([base, base + 1, base + 2]);
    }
    (p, t, m.2.clone())
}

/// A closed parametric surface on an `n × m` grid, (u, v) ↦ point, periodic in u; `poles` caps
/// v = 0 and v = 1 with single vertices (a sphere), otherwise v is periodic too (a torus).
fn grid(n: u32, m: u32, poles: bool, f: impl Fn(f64, f64) -> [f64; 3]) -> Mesh {
    let mut p = Vec::new();
    let rows = if poles { m - 1 } else { m };
    for j in 0..rows {
        for i in 0..n {
            let v = if poles { f64::from(j + 1) / f64::from(m) } else { f64::from(j) / f64::from(m) };
            p.push(f(f64::from(i) / f64::from(n), v));
        }
    }
    let id = |i: u32, j: u32| (j % rows) * n + i % n;
    let mut t = Vec::new();
    let quads = if poles { rows - 1 } else { rows };
    for j in 0..quads {
        for i in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            t.push([a, b, c]);
            t.push([a, c, d]);
        }
    }
    if poles {
        let (south, north) = (p.len() as u32, p.len() as u32 + 1);
        p.push(f(0.0, 0.0));
        p.push(f(0.0, 1.0));
        for i in 0..n {
            t.push([south, id(i + 1, 0), id(i, 0)]);
            t.push([north, id(i, rows - 1), id(i + 1, rows - 1)]);
        }
    }
    let faces = vec![0; t.len()];
    let mut mesh = (p, t, faces);
    // Orient outwards.
    let six_v: f64 = mesh.1.iter().map(|t| crate::vec3::det(mesh.0[t[0] as usize], mesh.0[t[1] as usize], mesh.0[t[2] as usize])).sum();
    if six_v < 0.0 {
        flip(&mut mesh);
    }
    mesh
}

fn torus(n: u32, m: u32) -> Mesh {
    grid(n, m, false, |u, v| {
        let (a, b) = (2.0 * PI * u, 2.0 * PI * v);
        [(3.0 + b.cos()) * a.cos(), (3.0 + b.cos()) * a.sin(), b.sin()]
    })
}

fn sphere(n: u32, m: u32) -> Mesh {
    grid(n, m, true, |u, v| {
        let (a, b) = (2.0 * PI * u, PI * v);
        [b.sin() * a.cos(), b.sin() * a.sin(), -b.cos()]
    })
}

fn moved(m: Mesh, d: [f64; 3]) -> Mesh {
    (m.0.iter().map(|p| [p[0] + d[0], p[1] + d[1], p[2] + d[2]]).collect(), m.1, m.2)
}

fn prepare(m: &Mesh) -> Result<BoundarySurface> {
    BoundarySurface::new(&m.0, &m.1, &m.2, &SurfaceOptions::default())
}

fn defects(m: &Mesh) -> Defects {
    match prepare(m) {
        Err(Error::Surface(d)) => d,
        other => panic!("expected defects, got {other:?}"),
    }
}

fn kinds(d: &Defects) -> Vec<&'static str> {
    let mut k: Vec<_> = d.listed.iter().map(Defect::kind).collect();
    k.dedup();
    k
}

#[test]
fn cube_is_a_valid_boundary_with_its_edges_as_features() {
    let s = prepare(&cube([0.0; 3], 2.0)).unwrap();
    assert_eq!((s.points().len(), s.triangles().len()), (8, 12));
    assert_eq!(s.shells().len(), 1);
    let sh = s.shells()[0];
    assert_eq!((sh.genus, sh.triangles, sh.vertices, sh.void), (0, 12, 8, false));
    assert!((s.volume() - 8.0).abs() < 1e-12);
    assert_eq!(s.feature_edges().len(), 12);
    assert!(s.feature_edges().iter().all(|e| e.kind == FeatureKind::FaceBoundary && (e.angle - 90.0).abs() < 1e-9));
    assert_eq!(s.feature_corners(), (0..8).collect::<Vec<u32>>());
    for (t, n) in s.neighbors().iter().enumerate() {
        for (k, &m) in n.iter().enumerate() {
            let (a, b) = (s.triangles()[t][k], s.triangles()[t][(k + 1) % 3]);
            let other = s.triangles()[m as usize];
            assert!((0..3).any(|j| other[j] == b && other[(j + 1) % 3] == a), "triangle {t} edge {k}");
        }
    }
    // Feature edges name first the triangle that traverses them low → high.
    for e in s.feature_edges() {
        let t = s.triangles()[e.triangles[0] as usize];
        assert!((0..3).any(|j| t[j] == e.vertices[0] && t[(j + 1) % 3] == e.vertices[1]));
    }
}

#[test]
fn a_triangle_soup_welds_into_the_solid() {
    let c = cube([1.0, -2.0, 0.5], 3.0);
    let s = prepare(&soup(&c)).unwrap();
    assert_eq!((s.points().len(), s.triangles().len()), (8, 12));
    assert_eq!(s.source_vertex().len(), 8);
    assert_eq!(s.source_triangle(), (0..12).collect::<Vec<u32>>());
    assert_eq!(s.tri_face(), c.2.as_slice());
    assert!((s.volume() - 27.0).abs() < 1e-9);

    // Tessellation noise below the tolerance welds; with exact welding it leaves holes.
    let mut noisy = soup(&c);
    for (i, p) in noisy.0.iter_mut().enumerate() {
        p[0] += 1e-12 * (i % 3) as f64;
    }
    assert_eq!(prepare(&noisy).unwrap().points().len(), 8);
    let exact = SurfaceOptions { weld_tolerance: 0.0, ..SurfaceOptions::default() };
    assert!(matches!(BoundarySurface::new(&noisy.0, &noisy.1, &noisy.2, &exact), Err(Error::Surface(_))));
    // Or an absolute distance (a CAD kernel's tolerance).
    let mut noisy = soup(&c);
    for (i, p) in noisy.0.iter_mut().enumerate() {
        p[1] += 1.5e-6 * (i % 2) as f64;
    }
    assert!(matches!(prepare(&noisy), Err(Error::Surface(_))));
    let kernel = SurfaceOptions { weld_distance: 2e-6, ..SurfaceOptions::default() };
    assert_eq!(BoundarySurface::new(&noisy.0, &noisy.1, &noisy.2, &kernel).unwrap().points().len(), 8);
}

#[test]
fn untagged_faces_use_the_dihedral_angle() {
    let c = cube([0.0; 3], 1.0);
    let s = BoundarySurface::new(&c.0, &c.1, &[], &SurfaceOptions::default()).unwrap();
    assert!(s.tri_face().iter().all(|&f| f == 0));
    assert_eq!(s.feature_edges().len(), 12);
    assert!(s.feature_edges().iter().all(|e| e.kind == FeatureKind::Sharp));
    let blunt = SurfaceOptions { feature_angle: 90.0, ..SurfaceOptions::default() };
    assert!(BoundarySurface::new(&c.0, &c.1, &[], &blunt).unwrap().feature_edges().is_empty());
}

#[test]
fn torus_and_sphere_have_their_genus() {
    let t = prepare(&torus(16, 8)).unwrap();
    assert_eq!((t.shells().len(), t.shells()[0].genus), (1, 1));
    let s = prepare(&sphere(24, 12)).unwrap();
    assert_eq!((s.shells().len(), s.shells()[0].genus), (1, 0));
    // Inscribed in the unit sphere: a little less than its volume.
    let ball = 4.0 / 3.0 * PI;
    assert!(s.volume() < ball && s.volume() > 0.9 * ball, "{}", s.volume());
    // A smooth surface tagged as one face has no features at a fine tessellation.
    assert!(prepare(&torus(64, 32)).unwrap().feature_edges().is_empty());
}

#[test]
fn a_void_is_a_negative_shell_inside_another() {
    let mut inner = cube([1.0; 3], 1.0);
    flip(&mut inner);
    let s = prepare(&merge(cube([0.0; 3], 3.0), inner)).unwrap();
    assert_eq!(s.shells().len(), 2);
    assert!(!s.shells()[0].void && s.shells()[1].void);
    assert!((s.volume() - 26.0).abs() < 1e-9);
    assert!(s.tri_shell()[..12].iter().all(|&x| x == 0) && s.tri_shell()[12..].iter().all(|&x| x == 1));

    // Two separate solids are two outer shells.
    let two = prepare(&merge(cube([0.0; 3], 1.0), cube([5.0, 0.0, 0.0], 1.0))).unwrap();
    assert!(two.shells().iter().all(|s| !s.void && s.volume > 0.0));

    // A void facing outwards, and a solid facing inwards, are inverted.
    let d = defects(&merge(cube([0.0; 3], 3.0), cube([1.0; 3], 1.0)));
    assert_eq!(d.listed, vec![Defect::Inverted { shell: 1, triangle: 12, volume: 1.0 }]);
    let mut c = cube([0.0; 3], 1.0);
    flip(&mut c);
    assert_eq!(kinds(&defects(&c)), vec!["inverted"]);
}

#[test]
fn holes_and_bad_edges_are_located() {
    let mut open = cube([0.0; 3], 1.0);
    open.1.remove(3);
    open.2.remove(3);
    let d = defects(&open);
    assert_eq!(d.total, 3);
    assert!(d.listed.iter().all(|x| x.kind() == "open_edge" && x.at().is_some()));
    assert!(d.to_string().starts_with("3 defects; open edge"), "{d}");

    let mut flipped = cube([0.0; 3], 1.0);
    flipped.1[5].swap(0, 1);
    let d = defects(&flipped);
    assert_eq!(kinds(&d), vec!["flipped"]);
    assert_eq!(d.total, 3);
    assert!(d.listed.iter().all(|x| matches!(x, Defect::Flipped { triangles, .. } if triangles.contains(&5))));

    // A fin on an edge: three triangles there.
    let mut fin = cube([0.0; 3], 1.0);
    fin.0.push([0.5, -1.0, -1.0]);
    fin.1.push([0, 1, 8]);
    fin.2.push(9);
    assert!(kinds(&defects(&fin)).contains(&"non_manifold_edge"));
}

#[test]
fn cubes_touching_at_a_corner_are_not_a_manifold() {
    let d = defects(&merge(cube([0.0; 3], 1.0), cube([1.0; 3], 1.0)));
    assert_eq!(d.listed, vec![Defect::NonManifoldVertex { vertex: 7, at: [1.0; 3] }]);
}

#[test]
fn intersecting_and_touching_solids_are_self_intersections() {
    let d = defects(&merge(cube([0.0; 3], 2.0), cube([1.0; 3], 2.0)));
    assert_eq!(kinds(&d), vec!["self_intersection"]);
    // A smaller block resting on a face: coplanar contact without shared vertices.
    let d = defects(&merge(cube([0.0; 3], 2.0), cube([2.0, 0.5, 0.5], 1.0)));
    assert_eq!(kinds(&d), vec!["self_intersection"]);
    // A void whose corner is pulled onto the outer surface (onto the bottom face's diagonal).
    let mut inner = cube([1.0; 3], 1.0);
    inner.0[0] = [1.0, 1.0, 0.0];
    flip(&mut inner);
    let d = defects(&merge(cube([0.0; 3], 3.0), inner));
    assert_eq!(kinds(&d), vec!["self_intersection"]);
}

#[test]
fn many_intersections_are_capped() {
    // Small cubes all crossing the top of a large one.
    let mut m = cube([0.0; 3], 100.0);
    for i in 0..60 {
        m = merge(m, cube([1.5 * f64::from(i) + 0.1, 0.1, 99.5], 1.0));
    }
    let d = defects(&m);
    assert_eq!(d.listed.len(), MAX_LISTED_DEFECTS);
    assert!(d.listed.windows(2).all(|w| match (&w[0], &w[1]) {
        (Defect::SelfIntersection { triangles: a, .. }, Defect::SelfIntersection { triangles: b, .. }) => a <= b,
        _ => false,
    }));
}

#[test]
fn hostile_input_is_an_error() {
    let c = cube([0.0; 3], 1.0);
    let mut nan = c.clone();
    nan.0[3][1] = f64::NAN;
    assert_eq!(defects(&nan).listed, vec![Defect::BadCoordinate { vertex: 3 }]);
    let mut far = c.clone();
    far.0[0][0] = 1e10;
    assert_eq!(kinds(&defects(&far)), vec!["bad_coordinate"]);
    let mut bad = c.clone();
    bad.1[2][1] = 99;
    assert_eq!(defects(&bad).listed, vec![Defect::BadIndex { triangle: 2 }]);
    let mut degenerate = c.clone();
    degenerate.0.push([2.0, 0.0, 0.0]);
    degenerate.1.push([0, 1, 8]);
    degenerate.2.push(0);
    assert_eq!(kinds(&defects(&degenerate)), vec!["degenerate"]);

    let o = SurfaceOptions::default();
    assert!(matches!(BoundarySurface::new(&[], &[], &[], &o), Err(ref e) if e.kind() == Some(crate::Kind::Invalid)));
    assert!(matches!(BoundarySurface::new(&c.0, &c.1, &[0, 1], &o), Err(ref e) if e.kind() == Some(crate::Kind::Invalid)));
    assert!(matches!(BoundarySurface::new(&[[1.0; 3]; 3], &[[0, 1, 2]], &[], &o), Err(ref e) if e.kind() == Some(crate::Kind::Invalid)));
    for w in [-1.0, f64::NAN, 1e-20, 0.5, f64::INFINITY] {
        assert!(BoundarySurface::new(&c.0, &c.1, &c.2, &SurfaceOptions { weld_tolerance: w, ..o }).is_err(), "{w}");
    }
    for d in [-1.0, f64::NAN, f64::INFINITY, 0.01] {
        assert!(BoundarySurface::new(&c.0, &c.1, &c.2, &SurfaceOptions { weld_distance: d, ..o }).is_err(), "{d}");
    }
    for a in [0.0, -1.0, 181.0, f64::NAN] {
        assert!(BoundarySurface::new(&c.0, &c.1, &c.2, &SurfaceOptions { feature_angle: a, ..o }).is_err(), "{a}");
    }
}

#[test]
fn tessellation_noise_snaps_into_the_predicates_domain() {
    let mut c = cube([0.0; 3], 1.0);
    c.0[0] = [1e-40, -1e-300, 0.0];
    let s = BoundarySurface::new(&c.0, &c.1, &c.2, &SurfaceOptions { weld_tolerance: 0.0, ..SurfaceOptions::default() }).unwrap();
    assert_eq!(s.points()[0], [0.0; 3]);
    assert!((0..8).all(|v| s.point3(v).is_some()));
}

#[test]
fn output_is_deterministic() {
    let m = merge(merge(soup(&torus(20, 10)), moved(sphere(16, 8), [0.0, 0.0, 10.0])), cube([20.0; 3], 1.0));
    let a = prepare(&m).unwrap();
    assert_eq!(a, prepare(&m).unwrap());
    assert_eq!(a.shells().iter().map(|s| s.genus).collect::<Vec<_>>(), vec![1, 0, 0]);
}

#[test]
fn random_soups_never_panic() {
    let mut seed = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let specials = [0.0, -0.0, 1.0, -1.0, 1e-40, 1e9, f64::NAN, f64::INFINITY, 0.5];
    for _ in 0..2000 {
        let n = (next() % 9 + 1) as usize;
        let p: Vec<[f64; 3]> =
            (0..n).map(|_| [0, 1, 2].map(|_| if next() % 4 == 0 { specials[(next() % 9) as usize] } else { (next() % 5) as f64 })).collect();
        let t: Vec<[u32; 3]> = (0..next() % 12).map(|_| [0, 1, 2].map(|_| (next() % (n as u64 + 1)) as u32)).collect();
        let _ = BoundarySurface::new(&p, &t, &[], &SurfaceOptions::default());
    }
}

fn p3(c: [f64; 3]) -> Point3 {
    Point3::try_from(c).unwrap()
}

#[test]
fn exact_pair_tests() {
    let [o, x, y, z] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]].map(p3);
    let t = [&o, &x, &y];
    // Segments: piercing, touching a vertex, coplanar crossing, coplanar outside, parallel.
    assert!(segment_meets_triangle(&p3([0.2, 0.2, -1.0]), &p3([0.2, 0.2, 1.0]), t));
    assert!(segment_meets_triangle(&p3([1.0, 0.0, 0.0]), &p3([1.0, 0.0, 5.0]), t));
    assert!(segment_meets_triangle(&p3([-1.0, 0.3, 0.0]), &p3([2.0, 0.3, 0.0]), t));
    assert!(!segment_meets_triangle(&p3([0.6, 0.6, 0.0]), &p3([2.0, 2.0, 0.0]), t));
    assert!(!segment_meets_triangle(&p3([0.0, 0.0, 1.0]), &p3([1.0, 0.0, 1.0]), t));
    // On the hypotenuse x + y = 1 of a triangle whose corners are not exact in binary.
    let (a, b) = (p3([0.1, 0.9, 0.0]), p3([0.9, 0.1, 0.0]));
    let mid = p3([0.5, 0.5, 0.0]);
    assert!(
        segment_meets_triangle(&p3([0.5, 0.5, -1.0]), &p3([0.5, 0.5, 1.0]), [&o, &a, &b])
            == (orient2d_sign(&a, &b, &mid) * orient2d_sign(&a, &b, &o) >= 0)
    );

    assert!(triangles_meet(t, [&p3([0.1, 0.1, -1.0]), &p3([0.1, 0.1, 1.0]), &p3([0.5, 0.1, 0.0])]));
    assert!(!triangles_meet(t, [&p3([0.1, 0.1, 0.5]), &p3([0.6, 0.1, 0.5]), &p3([0.1, 0.6, 0.5])]));

    // Shared edge o–x: folded flat onto t, flat on the other side, and bent out of plane.
    let ids = |s: [u32; 3], t: [u32; 3]| [s, t];
    let (w, w2, w3) = (p3([0.5, 0.5, 0.0]), p3([0.5, -0.5, 0.0]), p3([0.5, 0.5, 1.0]));
    assert!(surface_pair_intersects(ids([0, 1, 2], [1, 0, 3]), t, [&x, &o, &w]));
    assert!(!surface_pair_intersects(ids([0, 1, 2], [1, 0, 3]), t, [&x, &o, &w2]));
    assert!(!surface_pair_intersects(ids([0, 1, 2], [1, 0, 3]), t, [&x, &o, &w3]));

    // Shared vertex o: a coplanar triangle inside t's corner, one outside it, one crossing t, one
    // out of plane touching only at o, and one whose edge runs along t's edge o–x.
    let at_o = ids([0, 1, 2], [0, 3, 4]);
    assert!(surface_pair_intersects(at_o, t, [&o, &p3([0.3, 0.1, 0.0]), &p3([0.1, 0.3, 0.0])]));
    assert!(!surface_pair_intersects(at_o, t, [&o, &p3([-0.3, -0.1, 0.0]), &p3([-0.1, -0.3, 0.0])]));
    assert!(surface_pair_intersects(at_o, t, [&o, &p3([0.3, 0.3, -1.0]), &p3([0.3, 0.3, 1.0])]));
    assert!(!surface_pair_intersects(at_o, t, [&o, &z, &p3([-1.0, -1.0, 1.0])]));
    assert!(surface_pair_intersects(at_o, t, [&o, &p3([0.5, 0.0, 0.0]), &p3([0.5, -1.0, 0.0])]));
}

/// The exact orientation of `c` against the directed line `ab` seen from +z, as -1, 0, 1.
fn orient2d_sign(a: &Point3, b: &Point3, c: &Point3) -> i32 {
    let p = |q: &Point3| crate::predicates::Point2::new(q.coords()[0], q.coords()[1]).unwrap();
    crate::predicates::orient2d(&p(a), &p(b), &p(c)) as i32
}

/// A vertex of one face on an undivided edge of its neighbour (a T-junction, as independent
/// face tessellations leave where a face is split) is closed by splitting that edge.
#[test]
fn t_junctions_are_closed() {
    let (mut p, mut t, mut f) = cube([0.0; 3], 1.0);
    // Edge 0–2 (x = 0, z = 0) gets a vertex at its middle on the x = 0 side only.
    let m = p.len() as u32;
    p.push([0.0, 0.5, 0.0]);
    let i = t.iter().position(|x| *x == [0, 6, 2]).unwrap();
    t[i] = [0, 6, m];
    t.push([m, 6, 2]);
    f.push(f[i]);
    let (sp, st, sf) = soup(&(p, t, f));
    let s = BoundarySurface::new(&sp, &st, &sf, &SurfaceOptions { weld_tolerance: 0.0, weld_distance: 1e-9, feature_angle: 30.0 }).unwrap();
    assert!((s.volume() - 1.0).abs() < 1e-12, "{}", s.volume());
    assert_eq!(s.triangles().len(), 14);
    // Without a weld distance nothing is split, and the open edge is reported.
    assert!(BoundarySurface::new(&sp, &st, &sf, &SurfaceOptions { weld_tolerance: 0.0, weld_distance: 0.0, feature_angle: 30.0 }).is_err());
}
