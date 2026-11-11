use super::*;
use crate::size::UniformSize;
use crate::surface::{BoundarySurface, SurfaceOptions};
use crate::vec3::{cross, dot, sub};

/// Ear clipping of a simple planar polygon (counter-clockwise about `normal`).
fn ears(points: &[[f64; 3]], poly: &[u32]) -> Vec<[u32; 3]> {
    let p = |i: u32| points[i as usize];
    let normal = (0..poly.len()).fold([0.0; 3], |n, k| {
        let (a, b) = (p(poly[k]), p(poly[(k + 1) % poly.len()]));
        [n[0] + (a[1] - b[1]) * (a[2] + b[2]), n[1] + (a[2] - b[2]) * (a[0] + b[0]), n[2] + (a[0] - b[0]) * (a[1] + b[1])]
    });
    let turn = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| dot(cross(sub(b, a), sub(c, a)), normal);
    let mut v: Vec<u32> = poly.to_vec();
    let mut out = Vec::new();
    while v.len() > 3 {
        let n = v.len();
        let k = (0..n)
            .find(|&k| {
                let (a, b, c) = (p(v[(k + n - 1) % n]), p(v[k]), p(v[(k + 1) % n]));
                turn(a, b, c) > 0.0
                    && v.iter().all(|&o| {
                        let q = p(o);
                        o == v[(k + n - 1) % n]
                            || o == v[k]
                            || o == v[(k + 1) % n]
                            || !(turn(a, b, q) >= 0.0 && turn(b, c, q) >= 0.0 && turn(c, a, q) >= 0.0)
                    })
            })
            .unwrap();
        out.push([v[(k + n - 1) % n], v[k], v[(k + 1) % n]]);
        v.remove(k);
    }
    out.push([v[0], v[1], v[2]]);
    out
}

/// A closed surface from planar polygonal faces (counter-clockwise from outside), tagged by face.
pub(crate) fn polyhedron(points: &[[f64; 3]], faces: &[Vec<u32>]) -> BoundarySurface {
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    for (f, poly) in faces.iter().enumerate() {
        for t in ears(points, poly) {
            tris.push(t);
            tags.push(f as u32);
        }
    }
    BoundarySurface::new(points, &tris, &tags, &SurfaceOptions::default()).unwrap()
}

/// A prism over a counter-clockwise polygon from z0 to z1.
pub(crate) fn prism(poly: &[[f64; 2]], z0: f64, z1: f64) -> BoundarySurface {
    let n = poly.len() as u32;
    let mut pts: Vec<[f64; 3]> = poly.iter().map(|p| [p[0], p[1], z0]).collect();
    pts.extend(poly.iter().map(|p| [p[0], p[1], z1]));
    let mut faces = vec![(0..n).rev().collect::<Vec<_>>(), (n..2 * n).collect()];
    for k in 0..n {
        let j = (k + 1) % n;
        faces.push(vec![k, j, n + j, n + k]);
    }
    polyhedron(&pts, &faces)
}

/// The mesh has no crack ([`crate::quality::conformity`]): its free faces close up, no two nodes
/// coincide, one mid-side node per element edge.
fn assert_conforming(m: &HexMesh) {
    let nodes: Vec<(u32, [f64; 3])> = m.nodes.iter().enumerate().map(|(i, p)| (i as u32, *p)).collect();
    let elements = (0..m.element_count()).filter_map(|e| Some((e as u32, m.kind, m.element(e)?)));
    let c = crate::quality::conformity(&nodes, elements, 1e-9).unwrap();
    assert!(c.is_conforming(), "{:?} {:?} {:?}", c.open_edges.first(), c.coincident.first(), c.midside_mismatches.first());
    assert_eq!(c.free_faces.len(), m.boundary.len());
}

fn check(m: &HexMesh, volume: f64) {
    assert!(m.report.gate.passed());
    assert_conforming(m);
    assert!((m.report.volume - volume).abs() < 1e-9 * volume, "{} vs {volume}", m.report.volume);
    assert!(m.report.min_scaled_jacobian > 0.0);
    // Every boundary face is tagged and the boundary faces are exactly the hex faces used once.
    assert!(!m.boundary.is_empty());
}

#[test]
fn box_extrudes_into_a_grid() {
    let s = prism(&[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], 0.0, 3.0);
    let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions::default()).unwrap();
    check(&m, 24.0);
    assert_eq!(m.element_count(), 24, "{}", m.report.description);
    assert!((m.report.min_scaled_jacobian - 1.0).abs() < 1e-12);
    // Boundary: 2·(4·2 + 4·3 + 2·3) = 52 quads.
    assert_eq!(m.boundary.len(), 52);
}

/// Finer near the face x = 4: 0.25 within 0.5 of it, 1 elsewhere.
struct FineAtEnd;

impl crate::size::SizeField for FineAtEnd {
    fn size(&self, p: [f64; 3]) -> f64 {
        if p[0] > 3.5 { 0.25 } else { 1.0 }
    }
}

/// A layer spans the whole cross-section, so the layers follow a size field that is fine near
/// one side face only, wherever the surface lists its points: the 4 × 4 × 3 box (extruded along
/// z, its largest faces' normal) with its points in their own order and rotated so a corner on
/// the fine side comes first extrudes in the same 12 layers (3 / 0.25) either way.
/// The layer count once followed a line through the surface's first point alone, so it changed
/// with the order the kernel listed the faces in (SolidCraft v0.40).
#[test]
fn layers_follow_the_size_field_anywhere_on_the_cross_section() {
    let poly = [[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]];
    let s = prism(&poly, 0.0, 3.0);
    let n = poly.len() as u32;
    let mut pts: Vec<[f64; 3]> = poly.iter().map(|p| [p[0], p[1], 0.0]).collect();
    pts.extend(poly.iter().map(|p| [p[0], p[1], 3.0]));
    let mut faces = vec![(0..n).rev().collect::<Vec<_>>(), (n..2 * n).collect()];
    for k in 0..n {
        let j = (k + 1) % n;
        faces.push(vec![k, j, n + j, n + k]);
    }
    // Point k of the rotated list is point (k + 5) mod 8 of the original: the first is (4, 0, 3).
    let rotate = |v: u32| (v + 2 * n - 5) % (2 * n);
    let rotated_pts: Vec<[f64; 3]> = (0..2 * n).map(|k| pts[((k + 5) % (2 * n)) as usize]).collect();
    let rotated_faces: Vec<Vec<u32>> = faces.iter().map(|f| f.iter().map(|&v| rotate(v)).collect()).collect();
    let r = polyhedron(&rotated_pts, &rotated_faces);
    assert_eq!(r.points()[0], [4.0, 0.0, 3.0]);
    let opts = HexOptions { method: HexMethod::Extrude, ..HexOptions::default() };
    for (what, surface) in [("listed", &s), ("rotated", &r)] {
        let m = hex_mesh(surface, &[], &FineAtEnd, &opts).unwrap();
        check(&m, 48.0);
        assert!(m.report.description.to_string().contains("12 layers"), "{what}: {}", m.report.description);
        let z = along(&m, 2, |p| p[0].abs() < 1e-9 && p[1].abs() < 1e-9);
        assert_eq!(z.len(), 13, "{what}: {} {z:?}", m.report.description);
    }
}

#[test]
fn l_shape_extrudes_with_a_free_or_mapped_section() {
    let l = [[0.0, 0.0], [6.0, 0.0], [6.0, 2.0], [2.0, 2.0], [2.0, 6.0], [0.0, 6.0]];
    let s = prism(&l, 0.0, 2.0);
    let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions::default()).unwrap();
    check(&m, 40.0);
}

#[test]
fn quadratic_box() {
    let s = prism(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]], 0.0, 2.0);
    let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { quadratic: true, ..HexOptions::default() }).unwrap();
    check(&m, 8.0);
    assert_eq!(m.kind, ElementKind::Hex20);
    // 3³ corner nodes + 3·(2·3·3) = 54 edge nodes.
    assert_eq!(m.nodes.len(), 27 + 54);
}

/// A block whose top is the curved surface z = 2 + 0.3 sin(πx/4), tessellated finely.
fn curved_block() -> BoundarySurface {
    let n = 16usize;
    let top = |x: f64| 2.0 + 0.3 * (std::f64::consts::PI * x / 4.0).sin();
    let mut pts: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.0, 2.0, 0.0], [0.0, 2.0, 0.0]];
    // Top grid rows at y = 0 and y = 2: ids 4..4+n+1 and 5+n..
    for y in [0.0, 2.0] {
        for k in 0..=n {
            let x = 4.0 * k as f64 / n as f64;
            pts.push([x, y, top(x)]);
        }
    }
    let t0 = |k: usize| (4 + k) as u32;
    let t2 = |k: usize| (5 + n + k) as u32;
    let mut faces: Vec<Vec<u32>> = vec![vec![0, 3, 2, 1]];
    // y = 0 side: 0, 1, then the top edge from x = 4 back to 0.
    faces.push([0u32, 1].into_iter().chain((0..=n).rev().map(t0)).collect());
    faces.push([2u32, 3].into_iter().chain((0..=n).map(t2)).collect());
    faces.push(vec![1, 2, t2(n), t0(n)]);
    faces.push(vec![3, 0, t0(0), t2(0)]);
    // The top: strips of two triangles, one face.
    let mut tris: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();
    for (f, poly) in faces.iter().enumerate() {
        for t in ears(&pts, poly) {
            tris.push(t);
            tags.push(f as u32);
        }
    }
    for k in 0..n {
        tris.push([t0(k), t0(k + 1), t2(k + 1)]);
        tris.push([t0(k), t2(k + 1), t2(k)]);
        tags.extend([5, 5]);
    }
    BoundarySurface::new(&pts, &tris, &tags, &SurfaceOptions::default()).unwrap()
}

#[test]
fn mapped_block_with_a_curved_top() {
    let s = curved_block();
    let opts = HexOptions { method: HexMethod::Mapped, ..HexOptions::default() };
    let m = hex_mesh(&s, &[], &UniformSize(0.5), &opts).unwrap();
    assert!(m.report.gate.passed());
    assert!(m.report.min_scaled_jacobian > 0.8, "{}", m.report.min_scaled_jacobian);
    assert!(m.report.description.to_string().contains("mapped block of 8 × 4 × 4"), "{}", m.report.description);
    // Top nodes lie on the tessellated top.
    let quadratic = hex_mesh(&s, &[], &UniformSize(0.5), &HexOptions { quadratic: true, ..opts }).unwrap();
    assert!(quadratic.report.gate.passed());
    assert_conforming(&quadratic);
    assert!((quadratic.report.volume - s.volume()).abs() < 1e-3 * s.volume(), "{} vs {}", quadratic.report.volume, s.volume());
}

/// A thick-walled pipe along z (radii `ri` < `ro`, height `h`, `n` segments around, `n` even):
/// eight faces, each cylinder and annulus in two halves, with their analytic surfaces.
fn pipe(ri: f64, ro: f64, h: f64, n: usize) -> (BoundarySurface, Vec<Option<Surface>>) {
    let ring = |r: f64, z: f64| -> Vec<[f64; 3]> {
        (0..n)
            .map(|k| {
                let t = std::f64::consts::TAU * k as f64 / n as f64;
                [r * t.cos(), r * t.sin(), z]
            })
            .collect()
    };
    let mut pts = Vec::new();
    for (r, z) in [(ro, 0.0), (ri, 0.0), (ro, h), (ri, h)] {
        pts.extend(ring(r, z));
    }
    let id = |ring: usize, k: usize| (ring * n + k % n) as u32;
    let (bo, bi, to, ti) = (0, 1, 2, 3);
    let mut tris: Vec<[u32; 3]> = Vec::new();
    let mut tags: Vec<u32> = Vec::new();
    // A quad as two triangles, turned to face `outward`.
    let mut quad = |q: [u32; 4], outward: [f64; 3], tag: u32| {
        let p = |i: u32| pts[i as usize];
        let n0 = cross(sub(p(q[1]), p(q[0])), sub(p(q[2]), p(q[0])));
        let q = if dot(n0, outward) < 0.0 { [q[0], q[3], q[2], q[1]] } else { q };
        tris.push([q[0], q[1], q[2]]);
        tris.push([q[0], q[2], q[3]]);
        tags.extend([tag, tag]);
    };
    for k in 0..n {
        let half = u32::from(k >= n / 2);
        let t = std::f64::consts::TAU * (k as f64 + 0.5) / n as f64;
        let radial = [t.cos(), t.sin(), 0.0];
        quad([id(bo, k), id(bo, k + 1), id(to, k + 1), id(to, k)], radial, half);
        quad([id(bi, k), id(bi, k + 1), id(ti, k + 1), id(ti, k)], [-radial[0], -radial[1], 0.0], 2 + half);
        quad([id(bo, k), id(bo, k + 1), id(bi, k + 1), id(bi, k)], [0.0, 0.0, -1.0], 4 + half);
        quad([id(to, k), id(to, k + 1), id(ti, k + 1), id(ti, k)], [0.0, 0.0, 1.0], 6 + half);
    }
    let surface = BoundarySurface::new(&pts, &tris, &tags, &SurfaceOptions::default()).unwrap();
    let cyl = |r: f64| Some(Surface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: r });
    let plane = |z: f64, nz: f64| Some(Surface::Plane { point: [0.0, 0.0, z], normal: [0.0, 0.0, nz] });
    (surface, vec![cyl(ro), cyl(ro), cyl(ri), cyl(ri), plane(0.0, -1.0), plane(0.0, -1.0), plane(h, 1.0), plane(h, 1.0)])
}

/// Every node of the boundary faces tagged with `faces` lies at `radius` from the z axis.
fn on_radius(m: &HexMesh, faces: &[u32], radius: f64) -> bool {
    m.boundary.iter().filter(|b| faces.contains(&b.brep_face)).all(|b| {
        let local = m.kind.faces()[b.face as usize - 1];
        let e = m.element(b.element as usize).unwrap();
        local.iter().all(|&i| {
            let p = m.nodes[e[i] as usize];
            ((p[0] * p[0] + p[1] * p[1]).sqrt() - radius).abs() < 1e-9 * radius
        })
    })
}

#[test]
fn pipe_extrudes_and_revolves_with_nodes_on_the_cylinders() {
    let (s, faces) = pipe(2.0, 3.0, 4.0, 32);
    for method in [HexMethod::Auto, HexMethod::Revolve] {
        for quadratic in [false, true] {
            let opts =
                HexOptions { method, quadratic, chordal_tolerance: 3.0 * (1.0 - (std::f64::consts::PI / 32.0).cos()), ..HexOptions::default() };
            let m = hex_mesh(&s, &faces, &UniformSize(0.5), &opts).unwrap();
            let expected = if method == HexMethod::Auto { HexMethod::Extrude } else { HexMethod::Revolve };
            assert_eq!(m.report.method, expected, "{}", m.report.description);
            assert!(m.report.gate.passed());
            assert!(m.report.min_scaled_jacobian > 0.5, "{method:?} {quadratic}: {} {}", m.report.min_scaled_jacobian, m.report.description);
            assert!(on_radius(&m, &[0, 1], 3.0) && on_radius(&m, &[2, 3], 2.0), "{method:?} {quadratic}");
            let exact = std::f64::consts::PI * (9.0 - 4.0) * 4.0;
            let tolerance = if quadratic { 1e-3 } else { 2e-2 };
            assert!((m.report.volume - exact).abs() < tolerance * exact, "{method:?} {quadratic}: {} vs {exact}", m.report.volume);
            // Every boundary face is tagged with one of the pipe's faces.
            assert!(m.boundary.iter().all(|b| b.brep_face < 8));
        }
    }
}

#[test]
fn box_faces_are_tagged_and_meshes_are_deterministic() {
    let s = prism(&[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], 0.0, 3.0);
    let a = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions::default()).unwrap();
    let b = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions::default()).unwrap();
    assert_eq!(a, b);
    // prism(): face 0 bottom, 1 top, then the sides along the polygon's edges.
    let count = |f: u32| a.boundary.iter().filter(|x| x.brep_face == f).count();
    assert_eq!([count(0), count(1), count(2), count(3), count(4), count(5)], [8, 8, 12, 6, 12, 6]);
    // Each boundary face lies in its B-rep face's plane.
    let plane = |f: u32| -> (usize, f64) {
        match f {
            0 => (2, 0.0),
            1 => (2, 3.0),
            2 => (1, 0.0),
            3 => (0, 4.0),
            4 => (1, 2.0),
            _ => (0, 0.0),
        }
    };
    for b in &a.boundary {
        let (axis, value) = plane(b.brep_face);
        let e = a.element(b.element as usize).unwrap();
        for &i in ElementKind::Hex8.faces()[b.face as usize - 1] {
            assert!((a.nodes[e[i] as usize][axis] - value).abs() < 1e-12, "{b:?}");
        }
    }
}

#[test]
fn an_l_prism_sweeps_one_to_one() {
    let l = [[0.0, 0.0], [6.0, 0.0], [6.0, 2.0], [2.0, 2.0], [2.0, 6.0], [0.0, 6.0]];
    let s = prism(&l, 0.0, 2.0);
    let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { method: HexMethod::Sweep, ..HexOptions::default() }).unwrap();
    assert_eq!(m.report.method, HexMethod::Sweep);
    check(&m, 40.0);
}

/// Hex20 meshes are conforming whichever method makes them (#260), on a curved body too: the
/// mid-side nodes put on the curved top are shared by the elements on both sides of an edge.
#[test]
fn hex20_meshes_are_conforming_by_every_method() {
    let l = [[0.0, 0.0], [6.0, 0.0], [6.0, 2.0], [2.0, 2.0], [2.0, 6.0], [0.0, 6.0]];
    let s = prism(&l, 0.0, 2.0);
    for method in [HexMethod::Extrude, HexMethod::Sweep] {
        let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { method, quadratic: true, ..HexOptions::default() }).unwrap();
        assert_eq!(m.kind, ElementKind::Hex20);
        assert_conforming(&m);
    }
    let curved = curved_block();
    for method in [HexMethod::Extrude, HexMethod::Mapped] {
        let m = hex_mesh(&curved, &[], &UniformSize(0.5), &HexOptions { method, quadratic: true, ..HexOptions::default() }).unwrap();
        assert_conforming(&m);
    }
    let (pipe, faces) = pipe(1.0, 2.0, 2.0, 16);
    for method in [HexMethod::Extrude, HexMethod::Revolve] {
        let m = hex_mesh(&pipe, &faces, &UniformSize(0.5), &HexOptions { method, quadratic: true, ..HexOptions::default() }).unwrap();
        assert_conforming(&m);
    }
}

#[test]
fn unsweepable_bodies_say_why_and_hostile_input_is_an_error() {
    // A tetrahedron: every face is slanted to every other.
    let tet = polyhedron(
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &[vec![0, 2, 1], vec![0, 1, 3], vec![1, 2, 3], vec![0, 3, 2]],
    );
    match hex_mesh(&tet, &[], &UniformSize(0.2), &HexOptions::default()) {
        Err(e) if e.kind() == Some(crate::Kind::Unsupported) => {
            let why = e.to_string();
            for m in ["extrude", "sweep", "revolve", "mapped"] {
                assert!(why.contains(m), "{why}");
            }
        }
        other => panic!("{other:?}"),
    }
    let s = prism(&[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], 0.0, 3.0);
    let bad_sizes: [&dyn crate::size::SizeField; 3] = [&UniformSize(f64::NAN), &UniformSize(0.0), &UniformSize(-1.0)];
    for f in bad_sizes {
        assert!(hex_mesh(&s, &[], f, &HexOptions::default()).is_err());
    }
    assert!(matches!(hex_mesh(&s, &[], &UniformSize(1e-4), &HexOptions::default()), Err(ref e) if e.kind() == Some(crate::Kind::Limit)));
    assert!(
        matches!(hex_mesh(&s, &[], &UniformSize(0.1), &HexOptions { max_elements: 100, ..HexOptions::default() }), Err(ref e) if e.kind() == Some(crate::Kind::Limit))
    );
    assert!(hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { chordal_tolerance: f64::NAN, ..HexOptions::default() }).is_err());
    // Analytic surfaces with garbage numbers are treated as free-form, not trusted.
    let garbage = vec![Some(Surface::Cylinder { origin: [f64::NAN; 3], axis: [0.0; 3], radius: -1.0 }); 6];
    let m = hex_mesh(&s, &garbage, &UniformSize(1.0), &HexOptions::default()).unwrap();
    assert_eq!(m.element_count(), 24);
}

// ---- Conforming cells ----

/// A prism over a counter-clockwise polygon from z0 to z1 whose side from vertex k to k + 1 gets
/// tag `side_tags[k]` (bottom 0, top 1); quads of one tag make one face.
fn tagged_prism(poly: &[[f64; 2]], side_tags: &[u32], z0: f64, z1: f64) -> BoundarySurface {
    let n = poly.len();
    let mut pts: Vec<[f64; 3]> = poly.iter().map(|p| [p[0], p[1], z0]).collect();
    pts.extend(poly.iter().map(|p| [p[0], p[1], z1]));
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    for t in ears(&pts, &(0..n as u32).rev().collect::<Vec<_>>()) {
        tris.push(t);
        tags.push(0);
    }
    for t in ears(&pts, &(n as u32..2 * n as u32).collect::<Vec<_>>()) {
        tris.push(t);
        tags.push(1);
    }
    for (k, &tag) in side_tags.iter().enumerate().take(n) {
        let (a, b) = (k as u32, ((k + 1) % n) as u32);
        let (c, d) = (b + n as u32, a + n as u32);
        tris.push([a, b, c]);
        tris.push([a, c, d]);
        tags.extend([tag, tag]);
    }
    BoundarySurface::new(&pts, &tris, &tags, &SurfaceOptions::default()).unwrap()
}

fn cells_of<'a>(surfaces: &'a [BoundarySurface], faces: &'a [Vec<Option<Surface>>]) -> Vec<Cell<'a>> {
    surfaces.iter().zip(faces).map(|(surface, faces)| Cell { name: "cell", surface, faces }).collect()
}

/// Boxes `[x0, x1] × [y0, y1] × [0, h]` as prisms (faces 0 bottom, 1 top, 2..6 the sides).
fn boxes(spans: &[([f64; 2], [f64; 2])], h: f64) -> Vec<BoundarySurface> {
    spans.iter().map(|(x, y)| prism(&[[x[0], y[0]], [x[1], y[0]], [x[1], y[1]], [x[0], y[1]]], 0.0, h)).collect()
}

/// Every node of the mesh at a distinct position: the cells' nodes were merged.
fn distinct_positions(m: &HexMesh) -> bool {
    let mut keys: Vec<[u64; 3]> = m.nodes.iter().map(|p| p.map(|x| (x + 0.0).to_bits())).collect();
    keys.sort_unstable();
    keys.windows(2).all(|w| w[0] != w[1])
}

#[test]
fn boxes_side_by_side_conform() {
    // Two boxes sharing x = 2, then a 2 × 2 arrangement whose middle edge four cells share.
    let pairs = boxes(&[([0.0, 2.0], [0.0, 1.0]), ([2.0, 3.0], [0.0, 1.0])], 1.0);
    let none = vec![Vec::new(); 2];
    let m = hex_mesh_cells(&cells_of(&pairs, &none), &UniformSize(0.5), &HexOptions::default()).unwrap();
    assert_eq!(m.interfaces, 1);
    assert_eq!(m.blocks.iter().map(|b| b[0] * b[1] * b[2]).sum::<usize>(), m.mesh.element_count());
    assert_eq!(m.mesh.element_count(), 6 * 2 * 2);
    // One structured grid: 7 × 3 × 3 nodes, none duplicated on the interface.
    assert_eq!(m.mesh.nodes.len(), 7 * 3 * 3);
    assert!(distinct_positions(&m.mesh));
    assert!((m.mesh.report.volume - 3.0).abs() < 1e-12);
    assert!(m.mesh.report.gate.passed());
    // The boundary is the outside only: 2·(6·2 + 6·2 + 2·2) quads.
    assert_eq!(m.mesh.boundary.len(), 56);
    assert_eq!(m.element_cell.iter().filter(|&&c| c == 0).count(), 16);

    let grid = boxes(&[([0.0, 1.0], [0.0, 1.0]), ([1.0, 2.0], [0.0, 1.0]), ([0.0, 1.0], [1.0, 2.0]), ([1.0, 2.0], [1.0, 2.0])], 1.0);
    let none = vec![Vec::new(); 4];
    let a = hex_mesh_cells(&cells_of(&grid, &none), &UniformSize(0.25), &HexOptions::default()).unwrap();
    assert_eq!(a.interfaces, 4);
    assert_eq!(a.mesh.nodes.len(), 9 * 9 * 5);
    assert!(distinct_positions(&a.mesh));
    // Deterministic, and hex20 shares its mid-edge nodes across the interfaces too.
    let b = hex_mesh_cells(&cells_of(&grid, &none), &UniformSize(0.25), &HexOptions::default()).unwrap();
    assert_eq!(a, b);
    let q = hex_mesh_cells(&cells_of(&grid, &none), &UniformSize(0.25), &HexOptions { quadratic: true, ..HexOptions::default() }).unwrap();
    assert_eq!(q.mesh.kind, ElementKind::Hex20);
    // Corner nodes plus one node per lattice edge: 9·9·5 + (8·9·5)·2 + 9·9·4.
    assert_eq!(q.mesh.nodes.len(), 405 + 720 + 324);
    assert!(distinct_positions(&q.mesh));
    assert_conforming(&q.mesh);
}

#[test]
fn cells_that_do_not_meet_face_to_face_are_refused() {
    // One box against two half-height boxes: its face at x = 1 meets two faces (a T-junction).
    let t = vec![
        prism(&[[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 2.0]], 0.0, 1.0),
        prism(&[[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]], 0.0, 1.0),
        prism(&[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]], 0.0, 1.0),
    ];
    let none = vec![Vec::new(); 3];
    match hex_mesh_cells(&cells_of(&t, &none), &UniformSize(0.5), &HexOptions::default()) {
        Err(Error::Failed(why)) => assert!(why.contains("interface") || why.contains("conform"), "{why}"),
        other => panic!("{other:?}"),
    }
    // A cell that is no logical cube says which.
    let l = vec![prism(&[[0.0, 0.0], [6.0, 0.0], [6.0, 2.0], [2.0, 2.0], [2.0, 6.0], [0.0, 6.0]], 0.0, 2.0)];
    match hex_mesh_cells(&cells_of(&l, &[Vec::new()]), &UniformSize(1.0), &HexOptions::default()) {
        Err(Error::Unsupported(why)) => assert!(why.starts_with("cell `cell`: not a logical cube"), "{why}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn hostile_cell_input_is_an_error() {
    let two = boxes(&[([0.0, 1.0], [0.0, 1.0]), ([1.0, 2.0], [0.0, 1.0])], 1.0);
    let none = vec![Vec::new(); 2];
    let cells = cells_of(&two, &none);
    assert!(hex_mesh_cells(&[], &UniformSize(1.0), &HexOptions::default()).is_err());
    for method in [HexMethod::Sweep, HexMethod::Extrude, HexMethod::Revolve, HexMethod::Swept] {
        assert!(
            matches!(hex_mesh_cells(&cells, &UniformSize(1.0), &HexOptions { method, ..HexOptions::default() }), Err(ref e) if e.kind() == Some(crate::Kind::Invalid))
        );
    }
    let bad_sizes: [&dyn crate::size::SizeField; 3] = [&UniformSize(f64::NAN), &UniformSize(0.0), &UniformSize(-1.0)];
    for f in bad_sizes {
        assert!(hex_mesh_cells(&cells, f, &HexOptions::default()).is_err());
    }
    assert!(matches!(hex_mesh_cells(&cells, &UniformSize(1e-4), &HexOptions::default()), Err(ref e) if e.kind() == Some(crate::Kind::Limit)));
    assert!(
        matches!(hex_mesh_cells(&cells, &UniformSize(0.1), &HexOptions { max_elements: 100, ..HexOptions::default() }), Err(ref e) if e.kind() == Some(crate::Kind::Limit))
    );
    assert!(hex_mesh_cells(&cells, &UniformSize(1.0), &HexOptions { chordal_tolerance: -1.0, ..HexOptions::default() }).is_err());
    // The same cell twice: two cells on the same side of every face.
    let same = [two[0].clone(), two[0].clone()];
    assert!(hex_mesh_cells(&cells_of(&same, &none), &UniformSize(0.5), &HexOptions::default()).is_err());
    // Garbage analytic surfaces are not trusted.
    let garbage = vec![vec![Some(Surface::Cylinder { origin: [f64::NAN; 3], axis: [0.0; 3], radius: -1.0 }); 6]; 2];
    let m = hex_mesh_cells(&cells_of(&two, &garbage), &UniformSize(0.5), &HexOptions::default()).unwrap();
    assert_eq!(m.mesh.element_count(), 16);
}

/// The four cells of an O-grid round a hole of radius 1 in the square [−2, 2]², height `h`: each
/// cell the hole's quarter between two diagonals and a side of the square. The hole's quarter is
/// one face (tag 5) or two halves of the cylinder (tags 5 and 6: virtual topology joins them).
fn o_grid(h: f64, arc_segments: usize, split_arcs: bool) -> (Vec<BoundarySurface>, Vec<Vec<Option<Surface>>>) {
    use std::f64::consts::FRAC_PI_4;
    let mut surfaces = Vec::new();
    let mut faces = Vec::new();
    for q in 0..4 {
        let turn = q as f64 * 2.0 * FRAC_PI_4;
        let rot = |p: [f64; 2]| [p[0] * turn.cos() - p[1] * turn.sin(), p[0] * turn.sin() + p[1] * turn.cos()];
        let mut poly = vec![rot([FRAC_PI_4.cos(), -FRAC_PI_4.sin()]), rot([2.0, -2.0]), rot([2.0, 2.0])];
        let mut tags = vec![2, 3, 4];
        for k in 0..arc_segments {
            let t = FRAC_PI_4 - 2.0 * FRAC_PI_4 * k as f64 / arc_segments as f64;
            poly.push(rot([t.cos(), t.sin()]));
            tags.push(if split_arcs && 2 * k >= arc_segments { 6 } else { 5 });
        }
        // Exact corners where the cells meet.
        let snap = |p: [f64; 2]| p.map(|x| if (x.abs() - 2.0).abs() < 1e-12 { 2.0 * x.signum() } else { x });
        let poly: Vec<[f64; 2]> = poly.into_iter().map(snap).collect();
        surfaces.push(tagged_prism(&poly, &tags, 0.0, h));
        let cyl = Some(Surface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: 1.0 });
        faces.push(vec![None, None, None, None, None, cyl, cyl]);
    }
    (surfaces, faces)
}

#[test]
fn an_o_grid_round_a_hole_conforms_with_graded_radial_lines() {
    for split_arcs in [false, true] {
        let (surfaces, faces) = o_grid(0.5, 24, split_arcs);
        let opts = HexOptions { chordal_tolerance: 1.0 - (std::f64::consts::PI / 48.0).cos(), ..HexOptions::default() };
        let m = hex_mesh_cells(&cells_of(&surfaces, &faces), &UniformSize(0.5), &opts).unwrap();
        assert_eq!(m.interfaces, 4, "{split_arcs}");
        assert!(m.mesh.report.gate.passed());
        assert!(m.mesh.report.min_scaled_jacobian > 0.3, "{}", m.mesh.report.min_scaled_jacobian);
        assert!(distinct_positions(&m.mesh));
        // A ring of quads round the hole: every block is (around) × (radial) × (through).
        let sorted = |b: [usize; 3]| {
            let mut b = b;
            b.sort_unstable();
            b
        };
        assert!(m.blocks.iter().all(|&b| sorted(b) == sorted(m.blocks[0])), "{:?}", m.blocks);
        let per_block = m.blocks[0].iter().product::<usize>();
        assert_eq!(m.mesh.element_count(), 4 * per_block);
        // Nodes on the hole lie on it.
        let hole: Vec<u32> = if split_arcs { vec![5, 6] } else { vec![5] };
        assert!(on_radius(&m.mesh, &hole, 1.0), "{split_arcs}");
        let on_hole = m.mesh.boundary.iter().filter(|b| hole.contains(&b.brep_face)).count();
        assert!(on_hole > 0);
        // Radial grading: the first spacing off the hole is smaller than the last at the square.
        let radii: std::collections::BTreeSet<u64> =
            m.mesh.nodes.iter().filter(|p| p[2] == 0.0 && p[1].abs() < 1e-12 && p[0] > 0.0).map(|p| p[0].to_bits()).collect();
        let r: Vec<f64> = radii.into_iter().map(f64::from_bits).collect();
        assert!(r.len() >= 3, "{r:?}");
        let first = r[1] - r[0];
        let last = r[r.len() - 1] - r[r.len() - 2];
        assert!(first < 0.75 * last, "{r:?}");
        // Spacing grows by at most the growth rate plus the rounding of the count.
        assert!(r.windows(3).all(|w| (w[2] - w[1]) <= 1.35 * (w[1] - w[0])), "{r:?}");
    }
}

fn seed(
    name: &str,
    points: Vec<[f64; 3]>,
    count: usize,
    bias: Option<fieldcraft_model::Bias>,
    constraint: fieldcraft_model::SeedConstraint,
) -> crate::seeds::EdgeSeed {
    crate::seeds::EdgeSeed { name: name.into(), points, tolerance: 0.0, count, bias, constraint }
}

/// The distinct values of coordinate `k` of the nodes `on` accepts, sorted.
fn along(m: &HexMesh, k: usize, on: impl Fn([f64; 3]) -> bool) -> Vec<f64> {
    let mut v: Vec<f64> = m.nodes.iter().filter(|p| on(**p)).map(|p| p[k]).collect();
    v.sort_by(f64::total_cmp);
    v.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    v
}

#[test]
fn seeds_give_a_swept_block_exact_counts_and_geometric_spacing() {
    use fieldcraft_model::{Bias, BiasToward, SeedConstraint};
    let s = prism(&[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], 0.0, 3.0);
    let bias = Some(Bias { ratio: 4.0, toward: BiasToward::Start });
    let opts = HexOptions {
        seeds: vec![
            // Through the thickness, 7 layers small at z = 0 (the seed runs up).
            seed("thickness", vec![[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]], 7, bias, SeedConstraint::Exact),
            // Along x on the top: 5 divisions, the seed runs from x = 4 to 0.
            seed("length", vec![[4.0, 2.0, 3.0], [0.0, 2.0, 3.0]], 5, None, SeedConstraint::Exact),
        ],
        ..HexOptions::default()
    };
    let m = hex_mesh(&s, &[], &UniformSize(1.0), &opts).unwrap();
    check(&m, 24.0);
    assert_eq!(m.element_count(), 5 * 2 * 7, "{}", m.report.description);
    // The nodes on the seeded edges.
    assert_eq!(along(&m, 0, |p| (p[1] - 2.0).abs() < 1e-9 && (p[2] - 3.0).abs() < 1e-9).len(), 6);
    let z = along(&m, 2, |p| p[0].abs() < 1e-9 && p[1].abs() < 1e-9);
    assert_eq!(z.len(), 8, "{z:?}");
    // Element heights grow geometrically from the bottom, the last 4 times the first.
    let dz: Vec<f64> = z.windows(2).map(|w| w[1] - w[0]).collect();
    let q = 4f64.powf(1.0 / 6.0);
    assert!(dz.windows(2).all(|w| (w[1] / w[0] - q).abs() < 1e-9), "{dz:?}");
    // Seeded the other way round: small at the top.
    let down = HexOptions { seeds: vec![seed("t", vec![[0.0, 0.0, 3.0], [0.0, 0.0, 0.0]], 7, bias, SeedConstraint::Exact)], ..opts.clone() };
    let m = hex_mesh(&s, &[], &UniformSize(1.0), &down).unwrap();
    let z = along(&m, 2, |p| p[0].abs() < 1e-9 && p[1].abs() < 1e-9);
    assert!(z[7] - z[6] < z[1] - z[0], "{z:?}");
    // The one-to-one sweep and the mapped mesher keep the same seeds.
    for method in [HexMethod::Sweep, HexMethod::Mapped] {
        let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { method, ..opts.clone() }).unwrap();
        let z = along(&m, 2, |p| p[0].abs() < 1e-9 && p[1].abs() < 1e-9);
        let dz: Vec<f64> = z.windows(2).map(|w| w[1] - w[0]).collect();
        assert!(dz.len() == 7 && dz.windows(2).all(|w| (w[1] / w[0] - q).abs() < 1e-9), "{method:?} {dz:?}");
        assert_eq!(m.element_count(), 5 * 2 * 7, "{method:?}");
    }
    // A minimum raises a count; a free seed is a goal the other constraints may move (its
    // opposite edges want 4, so 6 holds and 2 would not).
    let min = HexOptions { seeds: vec![seed("x", vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]], 2, None, SeedConstraint::Min)], ..opts.clone() };
    let front = |p: [f64; 3]| p[1].abs() < 1e-9 && p[2].abs() < 1e-9;
    assert_eq!(along(&hex_mesh(&s, &[], &UniformSize(1.0), &min).unwrap(), 0, front).len(), 5);
    let free = HexOptions { seeds: vec![seed("x", vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]], 6, None, SeedConstraint::Free)], ..opts.clone() };
    assert_eq!(along(&hex_mesh(&s, &[], &UniformSize(1.0), &free).unwrap(), 0, front).len(), 7);
}

#[test]
fn conflicting_exact_seeds_are_refused_with_their_names() {
    use fieldcraft_model::SeedConstraint;
    let s = prism(&[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], 0.0, 3.0);
    // Opposite edges of the mapped block (and of the swept cross-section) with different counts.
    let opts = HexOptions {
        method: HexMethod::Mapped,
        seeds: vec![
            seed("front", vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]], 5, None, SeedConstraint::Exact),
            seed("back", vec![[0.0, 2.0, 3.0], [4.0, 2.0, 3.0]], 3, None, SeedConstraint::Exact),
        ],
        ..HexOptions::default()
    };
    let e = hex_mesh(&s, &[], &UniformSize(1.0), &opts).unwrap_err().to_string();
    assert!(e.contains("front (5 exact)") && e.contains("back (3 exact)"), "{e}");
    let e = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { method: HexMethod::Auto, ..opts }).unwrap_err().to_string();
    assert!(e.contains("front (5 exact)") && e.contains("back (3 exact)"), "{e}");
}

#[test]
fn a_seed_on_a_shared_edge_holds_in_every_cell() {
    use fieldcraft_model::{Bias, BiasToward, SeedConstraint};
    let pairs = boxes(&[([0.0, 2.0], [0.0, 1.0]), ([2.0, 3.0], [0.0, 1.0])], 1.0);
    let none = vec![Vec::new(); 2];
    // The edge both cells share, at x = 2, along z.
    let bias = Some(Bias { ratio: 3.0, toward: BiasToward::Both });
    let opts =
        HexOptions { seeds: vec![seed("shared", vec![[2.0, 0.0, 0.0], [2.0, 0.0, 1.0]], 6, bias, SeedConstraint::Exact)], ..HexOptions::default() };
    let m = hex_mesh_cells(&cells_of(&pairs, &none), &UniformSize(0.5), &opts).unwrap();
    assert!(m.blocks.iter().all(|b| b[2] == 6), "{:?}", m.blocks);
    assert!(distinct_positions(&m.mesh));
    let z = along(&m.mesh, 2, |p| (p[0] - 2.0).abs() < 1e-9 && p[1].abs() < 1e-9);
    assert_eq!(z.len(), 7, "{z:?}");
    // Small at both ends, symmetric.
    assert!((z[1] - z[0] - (z[6] - z[5])).abs() < 1e-9 && z[3] - z[2] > z[1] - z[0], "{z:?}");
}

/// A biased seed on one edge carries its spacing to the edges interval assignment ties to it
/// (#257): the parallel edges of a mapped block, the opposite curves of an extruded
/// cross-section's mapped region, the side edges of a sweep's linking faces, and the parallel
/// edges of every cell they reach. The grid stays a grid: every node lies on one of the seed's
/// planes.
#[test]
fn a_biased_seed_grades_the_edges_tied_to_it() {
    use fieldcraft_model::{Bias, BiasToward, SeedConstraint};
    let s = prism(&[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]], 0.0, 3.0);
    let bias = Some(Bias { ratio: 3.0, toward: BiasToward::Start });
    let all = |_: [f64; 3]| true;
    // Along x on the front bottom edge, small at x = 4 (the seed runs from 4 to 0).
    let x_seed = seed("x", vec![[4.0, 0.0, 0.0], [0.0, 0.0, 0.0]], 6, bias, SeedConstraint::Exact);
    for method in [HexMethod::Mapped, HexMethod::Extrude] {
        let m = hex_mesh(&s, &[], &UniformSize(1.0), &HexOptions { method, seeds: vec![x_seed.clone()], ..HexOptions::default() }).unwrap();
        let x = along(&m, 0, all);
        assert_eq!(x.len(), 7, "{method:?} {x:?}");
        assert!(x[6] - x[5] < x[1] - x[0], "{method:?} {x:?}");
    }
    // Through the thickness at x = y = 0, small at z = 0: the block mapped, and an L swept
    // (up z, the one way it sweeps).
    let l = prism(&[[0.0, 0.0], [6.0, 0.0], [6.0, 2.0], [2.0, 2.0], [2.0, 6.0], [0.0, 6.0]], 0.0, 3.0);
    let z_seed = seed("z", vec![[0.0, 0.0, 0.0], [0.0, 0.0, 3.0]], 6, bias, SeedConstraint::Exact);
    for (body, method) in [(&s, HexMethod::Mapped), (&l, HexMethod::Sweep)] {
        let m = hex_mesh(body, &[], &UniformSize(1.0), &HexOptions { method, seeds: vec![z_seed.clone()], ..HexOptions::default() }).unwrap();
        let z = along(&m, 2, all);
        assert_eq!(z.len(), 7, "{method:?} {z:?}");
        assert!(z[1] - z[0] < z[6] - z[5], "{method:?} {z:?}");
    }
    // Two cells: a seed on the far side of the first carries through their interface to the
    // far side of the second.
    let pairs = boxes(&[([0.0, 2.0], [0.0, 1.0]), ([2.0, 3.0], [0.0, 1.0])], 1.0);
    let none = vec![Vec::new(); 2];
    let opts = HexOptions { seeds: vec![seed("z", vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], 6, bias, SeedConstraint::Exact)], ..HexOptions::default() };
    let m = hex_mesh_cells(&cells_of(&pairs, &none), &UniformSize(0.5), &opts).unwrap();
    let z = along(&m.mesh, 2, all);
    assert_eq!(z.len(), 7, "{z:?}");
    assert!(z[1] - z[0] < z[6] - z[5], "{z:?}");
}

/// Closed faces in one piece, as a kernel gives them that does not split them at their seams
/// (SolidCraft from v0.40): built from a grid round the z axis, `n` segments round, the seam at
/// angle 0. Triangles are turned to face away from `inside` (a point or the axis).
struct Revolved {
    pts: Vec<[f64; 3]>,
    tris: Vec<[u32; 3]>,
    tags: Vec<u32>,
}

impl Revolved {
    fn tri(&mut self, t: [u32; 3], outward: [f64; 3]) {
        let p = |i: u32| self.pts[i as usize];
        let n = cross(sub(p(t[1]), p(t[0])), sub(p(t[2]), p(t[0])));
        self.tris.push(if dot(n, outward) < 0.0 { [t[0], t[2], t[1]] } else { t });
    }

    fn surface(&self, seams: &[Vec<u32>]) -> BoundarySurface {
        let s = BoundarySurface::new(&self.pts, &self.tris, &self.tags, &SurfaceOptions::default()).unwrap();
        let seams: Vec<Vec<[f64; 3]>> = seams.iter().map(|s| s.iter().map(|&v| self.pts[v as usize]).collect()).collect();
        s.with_seams(&seams, 0.0).unwrap()
    }
}

fn angle(k: usize, n: usize) -> f64 {
    std::f64::consts::TAU * k as f64 / n as f64
}

/// A cylinder of radius `r` and height `h` along z: its side one face (tag 0) with a seam at
/// angle 0, the bottom (1) and top (2) fans round their centres. Returns the surface with and
/// without its seam marked, and the faces' surfaces.
fn seamed_cylinder(r: f64, h: f64, n: usize, m: usize) -> (BoundarySurface, BoundarySurface, Vec<Option<Surface>>) {
    let mut g = Revolved { pts: Vec::new(), tris: Vec::new(), tags: Vec::new() };
    for j in 0..=m {
        for k in 0..n {
            let t = angle(k, n);
            g.pts.push([r * t.cos(), r * t.sin(), h * j as f64 / m as f64]);
        }
    }
    let id = |j: usize, k: usize| (j * n + k % n) as u32;
    let (bottom, top) = (g.pts.len() as u32, g.pts.len() as u32 + 1);
    g.pts.push([0.0, 0.0, 0.0]);
    g.pts.push([0.0, 0.0, h]);
    for j in 0..m {
        for k in 0..n {
            let t = angle(2 * k + 1, 2 * n);
            let out = [t.cos(), t.sin(), 0.0];
            g.tri([id(j, k), id(j, k + 1), id(j + 1, k + 1)], out);
            g.tri([id(j, k), id(j + 1, k + 1), id(j + 1, k)], out);
            g.tags.extend([0, 0]);
        }
    }
    for k in 0..n {
        g.tri([bottom, id(0, k), id(0, k + 1)], [0.0, 0.0, -1.0]);
        g.tri([top, id(m, k), id(m, k + 1)], [0.0, 0.0, 1.0]);
        g.tags.extend([1, 2]);
    }
    let seam: Vec<u32> = (0..=m).map(|j| id(j, 0)).collect();
    let surfaces = vec![
        Some(Surface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: r }),
        Some(Surface::Plane { point: [0.0; 3], normal: [0.0, 0.0, -1.0] }),
        Some(Surface::Plane { point: [0.0, 0.0, h], normal: [0.0, 0.0, 1.0] }),
    ];
    (g.surface(&[seam]), g.surface(&[]), surfaces)
}

/// A tube along z (radii `ri` < `ro`, height `h`), its walls each one face with a seam at angle
/// 0 (tags 0 outer, 1 inner) and its ends annuli (2 bottom, 3 top), `n` segments round.
fn seamed_tube(ri: f64, ro: f64, h: f64, n: usize) -> (BoundarySurface, Vec<Option<Surface>>) {
    let mut g = Revolved { pts: Vec::new(), tris: Vec::new(), tags: Vec::new() };
    for (r, z) in [(ro, 0.0), (ri, 0.0), (ro, h), (ri, h)] {
        for k in 0..n {
            let t = angle(k, n);
            g.pts.push([r * t.cos(), r * t.sin(), z]);
        }
    }
    let id = |ring: usize, k: usize| (ring * n + k % n) as u32;
    for k in 0..n {
        let t = angle(2 * k + 1, 2 * n);
        let out = [t.cos(), t.sin(), 0.0];
        for (q, outward, tag) in [
            ([id(0, k), id(0, k + 1), id(2, k + 1), id(2, k)], out, 0),
            ([id(1, k), id(1, k + 1), id(3, k + 1), id(3, k)], [-out[0], -out[1], 0.0], 1),
            ([id(0, k), id(0, k + 1), id(1, k + 1), id(1, k)], [0.0, 0.0, -1.0], 2),
            ([id(2, k), id(2, k + 1), id(3, k + 1), id(3, k)], [0.0, 0.0, 1.0], 3),
        ] {
            g.tri([q[0], q[1], q[2]], outward);
            g.tri([q[0], q[2], q[3]], outward);
            g.tags.extend([tag, tag]);
        }
    }
    let surfaces = vec![
        Some(Surface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: ro }),
        Some(Surface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: ri }),
        Some(Surface::Plane { point: [0.0; 3], normal: [0.0, 0.0, -1.0] }),
        Some(Surface::Plane { point: [0.0, 0.0, h], normal: [0.0, 0.0, 1.0] }),
    ];
    (g.surface(&[vec![id(0, 0), id(2, 0)], vec![id(1, 0), id(3, 0)]]), surfaces)
}

/// A whole circle in one edge, as a kernel gives a seamed tube's ends, is divided into at least
/// 16 pieces (22.5° each) however coarse the size: at a size of a sixth of the circumference it
/// got 6, a hexagon whose prism lost 17 % of the tube. Swept or extruded, the tube's volume is now
/// the exact one within the 16-gon's chords.
#[test]
fn a_whole_circle_gets_sixteen_divisions_at_a_coarse_size() {
    use std::f64::consts::PI;
    let (ri, ro, h) = (0.5, 1.0, 1.0);
    let (s, faces) = seamed_tube(ri, ro, h, 64);
    let exact = PI * (ro * ro - ri * ri) * h;
    // The 16-gon keeps (16 / 2π) sin(2π / 16) of the circle's area.
    let chords = 1.0 - 16.0 / (2.0 * PI) * (2.0 * PI / 16.0).sin();
    for method in [HexMethod::Sweep, HexMethod::Extrude] {
        let opts = HexOptions { method, chordal_tolerance: ro * (1.0 - (PI / 64.0).cos()), ..HexOptions::default() };
        let m = hex_mesh(&s, &faces, &UniformSize(2.0 * PI * 0.75 / 6.0), &opts).unwrap();
        assert!(m.report.gate.passed(), "{method:?}");
        let outer: std::collections::BTreeSet<[u64; 2]> = m
            .nodes
            .iter()
            .filter(|p| p[2] == 0.0 && ((p[0] * p[0] + p[1] * p[1]).sqrt() - ro).abs() < 1e-9)
            .map(|p| [p[0].to_bits(), p[1].to_bits()])
            .collect();
        assert!(outer.len() >= 16, "{method:?}: {} nodes round the outer circle", outer.len());
        let error = (exact - m.report.volume) / exact;
        assert!((0.0..=chords + 1e-3).contains(&error), "{method:?}: {} against {exact} ({error})", m.report.volume);
    }
}

/// A torus about z (radii `major` > `minor`), one face with its two seams: the meridian at
/// angle 0 round the tube and the parallel on the outer equator round the axis.
fn seamed_torus(major: f64, minor: f64, nu: usize, nv: usize) -> (BoundarySurface, BoundarySurface, Vec<Option<Surface>>) {
    let mut g = Revolved { pts: Vec::new(), tris: Vec::new(), tags: Vec::new() };
    for i in 0..nu {
        for j in 0..nv {
            let (u, v) = (angle(i, nu), angle(j, nv));
            let rho = major + minor * v.cos();
            g.pts.push([rho * u.cos(), rho * u.sin(), minor * v.sin()]);
        }
    }
    let id = |i: usize, j: usize| ((i % nu) * nv + j % nv) as u32;
    for i in 0..nu {
        for j in 0..nv {
            let (u, v) = (angle(2 * i + 1, 2 * nu), angle(2 * j + 1, 2 * nv));
            let out = [v.cos() * u.cos(), v.cos() * u.sin(), v.sin()];
            g.tri([id(i, j), id(i + 1, j), id(i + 1, j + 1)], out);
            g.tri([id(i, j), id(i + 1, j + 1), id(i, j + 1)], out);
            g.tags.extend([0, 0]);
        }
    }
    let meridian: Vec<u32> = (0..=nv).map(|j| id(0, j)).collect();
    let parallel: Vec<u32> = (0..=nu).map(|i| id(i, 0)).collect();
    let surfaces = vec![Some(Surface::Torus { center: [0.0; 3], axis: [0.0, 0.0, 1.0], major, minor })];
    (g.surface(&[meridian, parallel]), g.surface(&[]), surfaces)
}

/// A sphere of radius `r` about the origin, one face with its seam from pole to pole at angle 0.
fn seamed_sphere(r: f64, nu: usize, nv: usize) -> (BoundarySurface, Vec<Option<Surface>>) {
    let mut g = Revolved { pts: Vec::new(), tris: Vec::new(), tags: Vec::new() };
    let (north, south) = (0u32, 1u32);
    g.pts.push([0.0, 0.0, r]);
    g.pts.push([0.0, 0.0, -r]);
    // Rings 1 .. nv - 1 from the north pole.
    for j in 1..nv {
        let phi = std::f64::consts::PI * j as f64 / nv as f64;
        for i in 0..nu {
            let u = angle(i, nu);
            g.pts.push([r * phi.sin() * u.cos(), r * phi.sin() * u.sin(), r * phi.cos()]);
        }
    }
    let id = |j: usize, i: usize| (2 + (j - 1) * nu + i % nu) as u32;
    let centre = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0];
    for i in 0..nu {
        let (a, b) = (id(1, i), id(1, i + 1));
        let out = centre(g.pts[0], g.pts[a as usize], g.pts[b as usize]);
        g.tri([north, a, b], out);
        let (a, b) = (id(nv - 1, i), id(nv - 1, i + 1));
        let out = centre(g.pts[1], g.pts[a as usize], g.pts[b as usize]);
        g.tri([south, a, b], out);
        g.tags.extend([0, 0]);
    }
    for j in 1..nv - 1 {
        for i in 0..nu {
            let q = [id(j, i), id(j, i + 1), id(j + 1, i + 1), id(j + 1, i)];
            let out = centre(g.pts[q[0] as usize], g.pts[q[1] as usize], g.pts[q[2] as usize]);
            g.tri([q[0], q[1], q[2]], out);
            g.tri([q[0], q[2], q[3]], out);
            g.tags.extend([0, 0]);
        }
    }
    let seam: Vec<u32> = std::iter::once(north).chain((1..nv).map(|j| id(j, 0))).chain(std::iter::once(south)).collect();
    (g.surface(&[seam]), vec![Some(Surface::Sphere { center: [0.0; 3], radius: r })])
}

/// A whole cylinder's side, one face with a seam, is a linking face (Blacker 1996: its loop is
/// the bottom circle, the seam, the top circle and the seam back): the cylinder sweeps one to one
/// from its bottom to its top, its nodes on the cylinder. Without the seam the side has two
/// loops and is no linking face.
#[test]
fn a_cylinder_whose_side_is_one_face_sweeps_across_its_seam() {
    let (s, plain, faces) = seamed_cylinder(1.0, 2.0, 24, 4);
    let brep = Brep::new(&s, &faces).unwrap();
    let side = &brep.faces[0];
    assert_eq!(side.loops.len(), 1, "{:?}", side.loops);
    assert_eq!(side.loops[0].len(), 4, "bottom, seam, top, seam back: {:?}", side.loops[0]);
    let seam = side.loops[0].iter().find(|&&(e, _)| brep.is_seam(e)).unwrap().0;
    assert_eq!(side.loops[0].iter().filter(|&&(e, _)| e == seam).count(), 2);
    let opts = HexOptions { method: HexMethod::Sweep, ..HexOptions::default() };
    let m = hex_mesh(&s, &faces, &UniformSize(0.25), &opts).unwrap();
    check_curved(&m, &s);
    assert!(on_radius(&m, &[0], 1.0), "the side's nodes lie on the cylinder");
    let e = hex_mesh(&plain, &faces, &UniformSize(0.25), &opts).unwrap_err().to_string();
    assert!(e.contains("boundary loops"), "{e}");
}

/// A torus in one face, with its seams: the meridian seam is a level of the revolution, so the
/// torus revolves; without the seams it has no level to revolve from.
#[test]
fn a_torus_in_one_face_revolves_from_its_meridian_seam() {
    let (s, plain, faces) = seamed_torus(3.0, 1.0, 32, 16);
    let opts = HexOptions { method: HexMethod::Revolve, ..HexOptions::default() };
    let m = hex_mesh(&s, &faces, &UniformSize(0.4), &opts).unwrap();
    check_curved(&m, &s);
    let e = hex_mesh(&plain, &faces, &UniformSize(0.4), &opts).unwrap_err().to_string();
    assert!(e.contains("no source and target faces"), "{e}");
}

/// A sphere in one face with a seam from pole to pole: its loop is the seam there and back, and
/// the structured meshers refuse it saying why (it reaches its axis of revolution), never
/// failing on the seam.
#[test]
fn a_sphere_with_a_seam_is_refused_saying_why() {
    let (s, faces) = seamed_sphere(1.0, 24, 12);
    let brep = Brep::new(&s, &faces).unwrap();
    assert_eq!(brep.faces.len(), 1);
    assert_eq!(brep.faces[0].loops.len(), 1, "{:?}", brep.faces[0].loops);
    assert_eq!(brep.faces[0].loops[0].len(), 2);
    let e = hex_mesh(&s, &faces, &UniformSize(0.3), &HexOptions::default()).unwrap_err();
    assert!(matches!(e.kind(), Some(crate::Kind::Unsupported)), "{e}");
}

/// A hex mesh of a curved body: valid, conforming, its volume the faceted boundary's within the
/// chords' share.
fn check_curved(m: &HexMesh, s: &BoundarySurface) {
    assert!(m.report.gate.passed(), "{}", m.report.description);
    assert_conforming(m);
    assert!(m.report.min_scaled_jacobian > 0.2, "{}", m.report.min_scaled_jacobian);
    assert!((m.report.volume - s.volume()).abs() < 2e-2 * s.volume(), "{} vs {}", m.report.volume, s.volume());
}
