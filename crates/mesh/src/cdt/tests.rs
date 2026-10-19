//! Tests of the constrained Delaunay tetrahedralization: exact conformity to every input
//! triangle, the enclosed volume, the validity gate, and the classic polyhedra that have no
//! tetrahedralization without Steiner points.

use super::meet::{in_closed_tet, meet_properly};
use super::*;
use crate::predicates::{GenericPoint, Lnc, Point3};
use crate::surface::{BoundarySurface, SurfaceOptions};
use crate::vec3::{cross, dot, norm, sub};

/// A closed surface from points and triangles, each triangle oriented outward from `centre`
/// (the shapes here are star-shaped about it).
pub(crate) fn surface(points: &[[f64; 3]], triangles: &[[u32; 3]], centre: [f64; 3]) -> BoundarySurface {
    let oriented: Vec<[u32; 3]> = triangles
        .iter()
        .map(|&[a, b, c]| {
            let (pa, pb, pc) = (points[a as usize], points[b as usize], points[c as usize]);
            let n = cross(sub(pb, pa), sub(pc, pa));
            if dot(n, sub(pa, centre)) > 0.0 { [a, b, c] } else { [a, c, b] }
        })
        .collect();
    let faces: Vec<u32> = (0..oriented.len() as u32).collect();
    BoundarySurface::new(points, &oriented, &faces, &SurfaceOptions::default()).unwrap()
}

/// The checks every result must pass: the gate (validity, volume), every input triangle covered
/// by boundary faces of exactly its area, each boundary face on its triangle's plane.
fn check(s: &BoundarySurface, m: &TetMesh) {
    let report = m.gate(s.volume()).unwrap();
    assert!(report.passed(), "{report:?}");
    let mut area = vec![0.0; s.triangles().len()];
    for b in &m.boundary {
        let t = m.tets[b.tet as usize];
        let face = face_toward(t, b.opposite as usize).map(|v| m.points[v as usize]);
        area[b.triangle as usize] += 0.5 * norm(cross(sub(face[1], face[0]), sub(face[2], face[0])));
        let tri = s.triangles()[b.triangle as usize].map(|v| s.points()[v as usize]);
        let tn = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
        for p in face {
            let d = dot(tn, sub(p, tri[0])) / norm(tn);
            assert!(d.abs() < 1e-9, "{d}");
        }
    }
    for (i, t) in s.triangles().iter().enumerate() {
        let p = t.map(|v| s.points()[v as usize]);
        let want = 0.5 * norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
        assert!((area[i] - want).abs() <= 1e-9 * want.max(1e-300), "triangle {i}: {} vs {want}", area[i]);
    }
}

pub(crate) fn cube() -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let p: Vec<[f64; 3]> = (0..8).map(|i| [(i & 1) as f64, (i >> 1 & 1) as f64, (i >> 2 & 1) as f64]).collect();
    let quads = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]];
    let t = quads.iter().flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]).collect();
    (p, t)
}

#[test]
fn a_cube_conforms() {
    let (p, t) = cube();
    let s = surface(&p, &t, [0.5; 3]);
    let m = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    check(&s, &m);
    // A square's diagonal is cocircular with the other two corners, so it is not strongly
    // Delaunay; at most one split per face (its midpoint) recovers it.
    assert!(m.steiner.len() <= 6, "{:?}", m.steiner);
    assert_eq!(m.boundary.len(), 12 + 2 * m.steiner.len());
}

/// Schönhardt's twisted prism (Schönhardt 1928): no tetrahedralization without Steiner points.
fn schonhardt(twist_deg: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut p = Vec::new();
    for z in 0..2 {
        for i in 0..3 {
            let a = (i as f64) * 2.0 * std::f64::consts::PI / 3.0 + if z == 1 { twist_deg.to_radians() } else { 0.0 };
            p.push([a.cos(), a.sin(), z as f64 * 1.5]);
        }
    }
    let mut t = vec![[0, 1, 2], [3, 4, 5]];
    for i in 0..3u32 {
        let j = (i + 1) % 3;
        // The reflex diagonal: from bottom i to top j (the top is turned towards j).
        t.push([i, j, 3 + j]);
        t.push([i, 3 + j, 3 + i]);
    }
    (p, t)
}

#[test]
fn schonhardt_prism_gets_steiner_points_and_conforms() {
    let (p, t) = schonhardt(30.0);
    let s = surface(&p, &t, [0.0, 0.0, 0.75]);
    let m = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    check(&s, &m);
    assert!(!m.steiner.is_empty(), "Schönhardt's prism has no tetrahedralization without Steiner points");
    for sp in &m.steiner {
        assert!(sp.t > 0.0 && sp.t < 1.0);
    }
}

/// A UV sphere: many cospherical vertices.
pub(crate) fn uv_sphere(rings: usize, sectors: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut p = vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
    for r in 1..rings {
        let th = std::f64::consts::PI * r as f64 / rings as f64;
        for s in 0..sectors {
            let ph = 2.0 * std::f64::consts::PI * s as f64 / sectors as f64;
            p.push([th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()]);
        }
    }
    let id = |r: usize, s: usize| (2 + (r - 1) * sectors + s % sectors) as u32;
    let mut t = Vec::new();
    for s in 0..sectors {
        t.push([0, id(1, s), id(1, s + 1)]);
        t.push([1, id(rings - 1, s + 1), id(rings - 1, s)]);
        for r in 1..rings - 1 {
            t.push([id(r, s), id(r + 1, s), id(r + 1, s + 1)]);
            t.push([id(r, s), id(r + 1, s + 1), id(r, s + 1)]);
        }
    }
    (p, t)
}

#[test]
fn a_sphere_conforms_and_is_deterministic() {
    let (p, t) = uv_sphere(12, 24);
    let s = surface(&p, &t, [0.0; 3]);
    let a = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    check(&s, &a);
    let b = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    assert_eq!(a, b);
}

#[test]
fn caps_are_reported() {
    let (p, t) = schonhardt(30.0);
    let s = surface(&p, &t, [0.0, 0.0, 0.75]);
    let none = CdtOptions { max_steiner: 0, max_tets: 1_000_000 };
    assert!(matches!(tetrahedralize(&s, &none), Err(ref e) if e.kind() == Some(crate::Kind::Limit)));
    let few = CdtOptions { max_steiner: 1000, max_tets: 10 };
    assert!(matches!(tetrahedralize(&s, &few), Err(ref e) if e.kind() == Some(crate::Kind::Limit)));
}

fn e(x: f64, y: f64, z: f64) -> GenericPoint {
    GenericPoint::Explicit(Point3::new(x, y, z).unwrap())
}

#[test]
fn proper_meeting_handles_coplanar_and_shared_cases() {
    let p = vec![
        e(0.0, 0.0, 0.0),  // 0
        e(1.0, 0.0, 0.0),  // 1
        e(0.0, 1.0, 0.0),  // 2
        e(0.0, 0.0, 1.0),  // 3: tet 0 1 2 3
        e(1.0, 1.0, 0.0),  // 4: beside face 0 1 2, in its plane
        e(0.2, 0.2, 0.0),  // 5: inside face 0 1 2
        e(0.0, 0.0, -1.0), // 6: below
        e(0.5, 0.5, 2.0),  // 7: above, off the tet
        e(2.0, 2.0, 2.0),  // 8: far
        e(0.1, 0.1, 0.1),  // 9: inside the tet
    ];
    let tet = [0, 1, 2, 3];
    let yes = |b: &[u32]| assert!(meet_properly(&p, &tet, b).unwrap(), "{b:?}");
    let no = |b: &[u32]| assert!(!meet_properly(&p, &tet, b).unwrap(), "{b:?}");
    yes(&[0, 1, 2]); // a face
    yes(&[1, 2, 4]); // coplanar, sharing edge 1 2, on the other side
    yes(&[0, 1, 6]); // shares an edge, below
    yes(&[3, 7, 8]); // shares a vertex, away
    yes(&[6, 8, 4]); // disjoint
    no(&[0, 1, 5]); // coplanar, overlapping the face 0 1 2
    no(&[1, 2, 5]); // coplanar, overlapping
    no(&[6, 7, 4]); // passes through the tet
    no(&[3, 9, 7]); // shares a vertex, enters the tet
    no(&[5, 6, 8]); // touches the face 0 1 2 at an interior point
    assert!(in_closed_tet(&p, tet, 9).unwrap());
    assert!(in_closed_tet(&p, tet, 5).unwrap());
    assert!(!in_closed_tet(&p, tet, 8).unwrap());
    // An LNC point on edge 0 3 is on the tet's boundary.
    let mut q = p.clone();
    q.push(GenericPoint::Lnc(Lnc::new(&Point3::new(0.0, 0.0, 0.0).unwrap(), &Point3::new(0.0, 0.0, 1.0).unwrap(), 0.3).unwrap()));
    assert!(in_closed_tet(&q, tet, 10).unwrap());
    assert!(!meet_properly(&q, &tet, &[10, 6, 8]).unwrap(), "touches edge 0 3 at a non-shared point");
    // Tetrahedra sharing a face meet properly; overlapping ones do not.
    assert!(meet_properly(&p, &tet, &[0, 1, 2, 6]).unwrap());
    assert!(!meet_properly(&p, &tet, &[0, 1, 2, 9]).unwrap());
}

/// Two solids: a wide tetrahedron standing on the big triangle z = 0, and a small one whose top
/// vertex is just under that triangle's centre. The Delaunay edge from the big apex to the small
/// top pierces the big triangle while all its edges are present: only facet recovery (the
/// gift-wrapping of a cavity) restores it.
fn pierced() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<[f64; 3]>) {
    let p = vec![
        [-10.0, -10.0, 0.0],
        [10.0, -10.0, 0.0],
        [0.0, 10.0, 0.0],
        [0.0, 0.0, 5.0],
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [0.0, 1.0, -1.0],
        [0.0, 0.0, -0.01],
    ];
    let t = vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [2, 0, 3], [4, 5, 6], [4, 5, 7], [5, 6, 7], [6, 4, 7]];
    (p, t, vec![[0.0, -3.0, 1.0], [0.0, -0.3, -0.7]])
}

/// Orients each triangle outward from the centre of the shell it belongs to.
fn shells(points: &[[f64; 3]], triangles: &[[u32; 3]], centres: &[[f64; 3]], per_shell: usize) -> BoundarySurface {
    let oriented: Vec<[u32; 3]> = triangles
        .iter()
        .enumerate()
        .map(|(i, &[a, b, c])| {
            let centre = centres[i / per_shell];
            let (pa, pb, pc) = (points[a as usize], points[b as usize], points[c as usize]);
            let n = cross(sub(pb, pa), sub(pc, pa));
            if dot(n, sub(pa, centre)) > 0.0 { [a, b, c] } else { [a, c, b] }
        })
        .collect();
    let faces: Vec<u32> = (0..oriented.len() as u32).collect();
    BoundarySurface::new(points, &oriented, &faces, &SurfaceOptions::default()).unwrap()
}

#[test]
fn a_pierced_facet_is_recovered_by_gift_wrapping() {
    let (p, t, centres) = pierced();
    let s = shells(&p, &t, &centres, 4);
    assert_eq!(s.shells().len(), 2);
    let m = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    check(&s, &m);
    assert!(m.stats.cavities > 0, "{:?}", m.stats);
    assert!(m.stats.gift_wrapped_tets > 0);
}

/// A Chazelle-style polyhedron (Chazelle 1984): a block whose top is cut by V-notches along `x`
/// and whose bottom is cut by V-notches along `y`, the notch edges almost touching. Its
/// tetrahedralizations need many Steiner points.
fn chazelle(notches: usize, gap: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let n = 2 * notches;
    let at = |i: usize| i as f64 / n as f64;
    // Top height by y index: ridges at even, valleys at odd. Bottom height by x index.
    let top = |j: usize| if j.is_multiple_of(2) { 1.0 } else { 0.5 + gap / 2.0 };
    let bottom = |i: usize| if i.is_multiple_of(2) { 0.0 } else { 0.5 - gap / 2.0 };
    let mut p = Vec::new();
    let mut id = std::collections::HashMap::new();
    let mut vertex = |x: usize, y: usize, up: bool, p: &mut Vec<[f64; 3]>| -> u32 {
        *id.entry((x, y, up)).or_insert_with(|| {
            p.push([at(x), at(y), if up { top(y) } else { bottom(x) }]);
            (p.len() - 1) as u32
        })
    };
    let mut t = Vec::new();
    let quad = |q: [u32; 4], out: [f64; 3], p: &Vec<[f64; 3]>, t: &mut Vec<[u32; 3]>| {
        for tri in [[q[0], q[1], q[2]], [q[0], q[2], q[3]]] {
            let (a, b, c) = (p[tri[0] as usize], p[tri[1] as usize], p[tri[2] as usize]);
            let nrm = cross(sub(b, a), sub(c, a));
            t.push(if dot(nrm, out) > 0.0 { tri } else { [tri[0], tri[2], tri[1]] });
        }
    };
    for i in 0..n {
        for j in 0..n {
            let q = [vertex(i, j, true, &mut p), vertex(i + 1, j, true, &mut p), vertex(i + 1, j + 1, true, &mut p), vertex(i, j + 1, true, &mut p)];
            quad(q, [0.0, 0.0, 1.0], &p, &mut t);
            let q =
                [vertex(i, j, false, &mut p), vertex(i + 1, j, false, &mut p), vertex(i + 1, j + 1, false, &mut p), vertex(i, j + 1, false, &mut p)];
            quad(q, [0.0, 0.0, -1.0], &p, &mut t);
        }
    }
    for k in 0..n {
        for (fixed, out) in [(0, -1.0), (n, 1.0)] {
            // Walls x = fixed (along y) and y = fixed (along x).
            let q = [
                vertex(fixed, k, false, &mut p),
                vertex(fixed, k + 1, false, &mut p),
                vertex(fixed, k + 1, true, &mut p),
                vertex(fixed, k, true, &mut p),
            ];
            quad(q, [out, 0.0, 0.0], &p, &mut t);
            let q = [
                vertex(k, fixed, false, &mut p),
                vertex(k + 1, fixed, false, &mut p),
                vertex(k + 1, fixed, true, &mut p),
                vertex(k, fixed, true, &mut p),
            ];
            quad(q, [0.0, out, 0.0], &p, &mut t);
        }
    }
    (p, t)
}

#[test]
fn chazelle_polyhedron_conforms() {
    let (p, t) = chazelle(3, 0.02);
    let faces: Vec<u32> = (0..t.len() as u32).collect();
    let s = BoundarySurface::new(&p, &t, &faces, &SurfaceOptions::default()).unwrap();
    assert_eq!(s.shells().len(), 1);
    let m = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    check(&s, &m);
    assert!(!m.steiner.is_empty());
}

#[test]
fn hostile_and_degenerate_surfaces_are_errors() {
    // A flat "solid" (all points coplanar) bounds no volume.
    let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
    let t = vec![[0, 1, 2], [1, 3, 2], [0, 2, 1], [1, 2, 3]];
    let faces = vec![0; 4];
    if let Ok(s) = BoundarySurface::new(&p, &t, &faces, &SurfaceOptions::default()) {
        assert!(tetrahedralize(&s, &CdtOptions::default()).is_err());
    }
}

#[test]
fn the_gate_reports_an_inverted_tetrahedron() {
    let (p, t) = cube();
    let s = surface(&p, &t, [0.5; 3]);
    let mut m = tetrahedralize(&s, &CdtOptions::default()).unwrap();
    m.tets[0].swap(0, 1);
    let report = m.gate(s.volume()).unwrap();
    assert_eq!(report.invalid, vec![1]);
    assert!(report.into_result().is_err());
}
