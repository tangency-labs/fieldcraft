use super::complex::*;
use super::measure;
use super::ops::*;
use crate::predicates::Point3;

/// A cube [0, 1]³ split into six tetrahedra around the diagonal 0–6, boundary faces tagged by
/// cube side (0..6).
pub(crate) fn cube() -> TetComplex {
    let p: Vec<[f64; 3]> = (0..8).map(|i| [(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64]).collect();
    // Vertex i = (x, y, z) bits; diagonal 0 → 7.
    let raw = [[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7], [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]];
    let pts: Vec<Point3> = p.iter().map(|c| Point3::try_from(*c).unwrap()).collect();
    let tets: Vec<[u32; 4]> = raw
        .iter()
        .map(|t| {
            let t = t.map(|v: usize| v as u32);
            let x = t.map(|v| pts[v as usize].coords());
            if measure::volume6(&x) > 0.0 { t } else { [t[1], t[0], t[2], t[3]] }
        })
        .collect();
    let side = |f: [u32; 3]| -> Option<u32> {
        let c = f.map(|v| p[v as usize]);
        for axis in 0..3 {
            for (s, val) in [(0, 0.0), (1, 1.0)] {
                if c.iter().all(|q| q[axis] == val) {
                    return Some((axis * 2 + s) as u32);
                }
            }
        }
        None
    };
    let mut boundary = Vec::new();
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            if let Some(tag) = side(face_toward(*vs, k)) {
                boundary.push((t as u32, k, tag));
            }
        }
    }
    let kinds = vec![VertexKind::Corner; 8];
    TetComplex::new(&p, kinds, &tets, &boundary, &[], 1000, 10_000).unwrap()
}

fn volume(cx: &TetComplex) -> f64 {
    cx.live_tets().map(|(_, vs)| measure::volume6(&cx.coords(vs).unwrap()) / 6.0).sum()
}

#[test]
fn cube_builds_and_rings() {
    let mut cx = cube();
    assert_eq!(cx.live_count(), 6);
    assert!((volume(&cx) - 1.0).abs() < 1e-15);
    let ring = cx.edge_ring(0, 7).unwrap().unwrap();
    assert!(ring.closed);
    assert_eq!(ring.verts.len(), 6);
    assert_eq!(ring.tets.len(), 6);
    let ring = cx.edge_ring(0, 1).unwrap().unwrap();
    assert!(!ring.closed);
    assert_eq!(ring.verts.len(), ring.tets.len() + 1);
}

#[test]
fn interior_insertion_and_rollback() {
    let mut cx = cube();
    cx.journal.recording = true;
    let start = cx.journal.steps.len();
    let p = Point3::try_from([0.5, 0.5, 0.5]).unwrap();
    let t = match cx.walk(0, &p, 100).unwrap().unwrap() {
        Walk::Inside(t) => t,
        w => panic!("{w:?}"),
    };
    let v = insert_delaunay(&mut cx, p, VertexKind::Interior, Site::Interior(t), 1000).unwrap().unwrap();
    assert!(cx.is_used(v));
    assert!((volume(&cx) - 1.0).abs() < 1e-15);
    assert!(cx.live_count() > 6);
    cx.rollback(start).unwrap();
    cx.validate().unwrap();
    assert_eq!(cx.live_count(), 6);
    assert_eq!(cx.point_count(), 8);
    assert!((volume(&cx) - 1.0).abs() < 1e-15);
}

#[test]
fn boundary_face_insertion_keeps_tags() {
    let mut cx = cube();
    // A point on the z = 0 side.
    let p = Point3::try_from([0.3, 0.6, 0.0]).unwrap();
    let mut site = None;
    for (t, _) in cx.live_tets().collect::<Vec<_>>() {
        let vs = cx.tet(t).unwrap();
        for k in 0..4 {
            if cx.adj(t).unwrap()[k] == NONE && cx.tags(t).unwrap()[k] == 4 {
                let f = face_toward(vs, k).map(|v| cx.xyz(v).unwrap());
                // inside the triangle in xy
                let s = |a: [f64; 3], b: [f64; 3]| (b[0] - a[0]) * (0.6 - a[1]) - (b[1] - a[1]) * (0.3 - a[0]);
                let signs = [s(f[0], f[1]), s(f[1], f[2]), s(f[2], f[0])];
                if signs.iter().all(|x| *x > 0.0) || signs.iter().all(|x| *x < 0.0) {
                    site = Some(Site::Face(t, k));
                }
            }
        }
    }
    let v = insert_delaunay(&mut cx, p, VertexKind::Face(4), site.unwrap(), 1000).unwrap().unwrap();
    assert!(cx.is_used(v));
    assert!((volume(&cx) - 1.0).abs() < 1e-15);
    let c = cx.compact().unwrap();
    let area: f64 = c
        .boundary
        .iter()
        .filter(|b| b.2 == 4)
        .map(|&(t, k, _)| {
            let f = face_toward(c.tets[t as usize], k).map(|v| c.points[v as usize]);
            0.5 * measure::norm(measure::cross(measure::sub(f[1], f[0]), measure::sub(f[2], f[0])))
        })
        .sum();
    assert!((area - 1.0).abs() < 1e-15, "{area}");
}

#[test]
fn edge_and_face_removal_keep_the_mesh_valid() {
    let mut cx = cube();
    let p = Point3::try_from([0.5, 0.45, 0.55]).unwrap();
    let Walk::Inside(t) = cx.walk(0, &p, 100).unwrap().unwrap() else { panic!() };
    insert_delaunay(&mut cx, p, VertexKind::Interior, Site::Interior(t), 1000).unwrap().unwrap();
    let mut changed = 0;
    for _ in 0..3 {
        let edges = edges_of(&cx, &cx.live_tets().map(|(t, _)| t).collect::<Vec<_>>()).unwrap();
        for (a, b) in edges {
            if cx.edge_tet(a, b).unwrap().is_some() && remove_edge(&mut cx, a, b).unwrap() {
                changed += 1;
            }
        }
        let tets: Vec<u32> = cx.live_tets().map(|(t, _)| t).collect();
        for t in tets {
            for k in 0..4 {
                if cx.is_live(t) && remove_face(&mut cx, t, k).unwrap() {
                    changed += 1;
                }
            }
        }
    }
    assert!(changed > 0);
    cx.validate().unwrap();
    assert!((volume(&cx) - 1.0).abs() < 1e-14);
    for (_, vs) in cx.live_tets() {
        assert_eq!(cx.orient(vs).unwrap(), crate::predicates::Sign::Positive);
    }
}

mod pipeline {
    use super::super::geometry::AnalyticSurface;
    use super::super::{Order, QualityOptions, quality_mesh};
    use crate::cdt::tests::{cube, surface, uv_sphere};
    use crate::cdt::{CdtOptions, tetrahedralize};
    use crate::surface::{BoundarySurface, SurfaceOptions};

    /// A 10 mm cube with one B-rep face per side.
    fn box_surface() -> BoundarySurface {
        let (p, t) = cube();
        let p: Vec<[f64; 3]> = p.iter().map(|c| c.map(|x| 10.0 * x)).collect();
        let s = surface(&p, &t, [5.0; 3]);
        // One face tag per side: the side's axis and end.
        let faces: Vec<u32> = s
            .triangles()
            .iter()
            .map(|t| {
                let c = t.map(|v| s.points()[v as usize]);
                (0..3).find_map(|k| [0.0, 10.0].iter().position(|&x| c.iter().all(|q| q[k] == x)).map(|e| (2 * k + e) as u32)).unwrap()
            })
            .collect();
        BoundarySurface::new(s.points(), s.triangles(), &faces, &SurfaceOptions::default()).unwrap()
    }

    /// A point of a side's plane past the side's edge is off the face, however large the side's
    /// triangles: the fillet box's nodes landed on the untrimmed planes, keeping the sharp
    /// corners (#402, #418).
    #[test]
    fn a_plane_is_trimmed_to_its_face() {
        let s = box_surface();
        let analytic = (0..6)
            .map(|f: usize| {
                let (k, e) = (f / 2, f % 2);
                let mut normal = [0.0; 3];
                normal[k] = if e == 0 { -1.0 } else { 1.0 };
                let mut point = [0.0; 3];
                point[k] = 10.0 * e as f64;
                Some(crate::tet::geometry::AnalyticSurface::Plane { normal, point })
            })
            .collect();
        let t = crate::tet::geometry::Tessellation::new(&s, analytic);
        let q = t.project_trimmed(0, [0.5, 11.0, 5.0]).unwrap();
        assert!(crate::vec3::dist(q, [0.0, 10.0, 5.0]) < 1e-9, "{q:?}");
        let q = t.project_trimmed(0, [0.5, 9.0, 5.0]).unwrap();
        assert!(crate::vec3::dist(q, [0.0, 9.0, 5.0]) < 1e-9, "{q:?}");
    }

    #[test]
    fn a_box_is_refined_improved_and_deterministic() {
        let s = box_surface();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let o = QualityOptions { size: Some(2.5), ..QualityOptions::default() };
        let a = quality_mesh(&s, &cdt, Vec::new(), &o).unwrap();
        let r = &a.report;
        assert!(r.output.tets > 100, "{r:?}");
        assert_eq!(r.output.inadmissible, 0);
        assert_eq!(r.output.slivers, 0, "{r:?}");
        assert!(r.output.min_dihedral > 25.0 && r.output.max_dihedral < 155.0, "{r:?}");
        // Flat faces: the volume is the box's, and each side keeps its tag on its own plane.
        assert!((r.volume - 1000.0).abs() < 1e-9, "{}", r.volume);
        for b in &a.boundary {
            let t: [u32; 4] = a.elements[b.tet as usize].clone().try_into().unwrap();
            let f = super::super::complex::face_toward(t, b.opposite as usize).map(|v| a.points[v as usize]);
            let (k, end) = ((b.face / 2) as usize, if b.face % 2 == 0 { 0.0 } else { 10.0 });
            assert!(f.iter().all(|p| p[k] == end), "face {} at {f:?}", b.face);
        }
        // The same input gives the same mesh, bit for bit.
        let b = quality_mesh(&s, &cdt, Vec::new(), &o).unwrap();
        assert_eq!(a, b);
    }

    /// Refinement splitting a face edge whose midpoint lies inside a feature edge's diametral
    /// sphere splits the feature edge instead (#303): otherwise each split halves the distance
    /// of the face vertices to the curve and they pile up against it.
    #[test]
    fn a_face_edge_split_that_encroaches_a_curve_splits_the_curve() {
        use super::super::complex::{NONE, VertexKind, face_toward};
        use super::super::geometry::Tessellation;
        use super::super::ops::{Site, insert_delaunay};
        use super::super::refine::{Context, MAX_CAVITY, split_edge};
        use super::super::{complex_from_cdt, measure};
        use crate::predicates::Point3;
        use crate::size::UniformSize;
        // Whether `p` lies in triangle `f` of the z = 0 plane.
        let contains = |f: [[f64; 3]; 3], p: [f64; 3]| {
            let side = |a: [f64; 3], b: [f64; 3]| (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
            let s = [side(f[0], f[1]), side(f[1], f[2]), side(f[2], f[0])];
            s.iter().all(|&x| x >= 0.0) || s.iter().all(|&x| x <= 0.0)
        };
        let s = box_surface();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let (mut cx, curve_faces) = complex_from_cdt(&s, &cdt, &QualityOptions::default()).unwrap();
        let (geometry, size) = (Tessellation::new(&s, Vec::new()), UniformSize(5.0));
        let ctx = Context { size: &size, geometry: &geometry, curve_faces: &curve_faces };
        // A curve edge on the z = 0 side, and a face vertex beside it on that side.
        let ((a, b), _) = cx.curve_edges().into_iter().find(|&((a, b), _)| [a, b].iter().all(|&v| cx.xyz(v).unwrap()[2] == 0.0)).unwrap();
        let (xa, xb) = (cx.xyz(a).unwrap(), cx.xyz(b).unwrap());
        let len = measure::dist(xa, xb);
        let centre = [5.0, 5.0, 0.0];
        let mid = measure::mid(xa, xb);
        let inward = measure::scale(measure::sub(centre, mid), 0.3 * len / measure::dist(centre, mid));
        let q = Point3::try_from(measure::add(mid, inward)).unwrap();
        let (t, k) = cx
            .live_tets()
            .flat_map(|(t, vs)| (0..4).map(move |k| (t, vs, k)))
            .find(|&(t, vs, k)| {
                cx.adj(t).unwrap()[k] == NONE && cx.tags(t).unwrap()[k] == 4 && contains(face_toward(vs, k).map(|v| cx.xyz(v).unwrap()), q.coords())
            })
            .map(|(t, _, k)| (t, k))
            .unwrap();
        let v = insert_delaunay(&mut cx, q, VertexKind::Face(4), Site::Face(t, k), MAX_CAVITY).unwrap().unwrap();
        assert!(cx.edge_tet(v, a).unwrap().is_some());
        let w = split_edge(&mut cx, &ctx, v, a, 4, 0.0).unwrap().unwrap();
        assert!(matches!(cx.kind(w), VertexKind::Curve(_)), "{:?}", cx.kind(w));
        assert!(measure::dist(cx.xyz(w).unwrap(), mid) < 1e-9);
        cx.validate().unwrap();
    }

    #[test]
    fn a_seeded_edge_gets_the_seeded_spacing_and_the_size_grows_from_it() {
        use crate::seeds::EdgeSeed;
        use fieldcraft_model::SeedConstraint;
        let s = box_surface();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let edge = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let seed = EdgeSeed { name: "front".into(), points: edge, tolerance: 0.0, count: 16, bias: None, constraint: SeedConstraint::Exact };
        let plain = QualityOptions { size: Some(5.0), seeds: vec![EdgeSeed { count: 1, ..seed.clone() }], ..QualityOptions::default() };
        let seeded = QualityOptions { seeds: vec![seed], ..plain.clone() };
        let a = quality_mesh(&s, &cdt, Vec::new(), &plain).unwrap();
        let b = quality_mesh(&s, &cdt, Vec::new(), &seeded).unwrap();
        let (na, nb) = (a.report.seeds[0].2, b.report.seeds[0].2);
        // Unseeded, the 10 mm edge gets about two 5 mm elements; seeded, about sixteen.
        assert!(na <= 4, "{na}");
        assert_eq!(nb, 16);
        assert_eq!(b.report.seeds[0].1, 16);
        assert!(b.report.output.tets > a.report.output.tets && b.report.output.inadmissible == 0);
        // Biased: the mesh's nodes on the edge are the seed's.
        let bias = Some(fieldcraft_model::Bias { ratio: 4.0, toward: fieldcraft_model::BiasToward::Start });
        let biased = QualityOptions { seeds: vec![EdgeSeed { count: 8, bias, ..seeded.seeds[0].clone() }], ..plain.clone() };
        let c = quality_mesh(&s, &cdt, Vec::new(), &biased).unwrap();
        let mut x: Vec<f64> = c.points.iter().filter(|p| p[1].abs() < 1e-9 && p[2].abs() < 1e-9).map(|p| p[0]).collect();
        x.sort_by(f64::total_cmp);
        let want: Vec<f64> = biased.seeds[0].fractions(8, false).iter().map(|f| 10.0 * f).collect();
        assert_eq!(x.len(), want.len(), "{x:?}");
        assert!(x.iter().zip(&want).all(|(a, b)| (a - b).abs() < 1e-9), "{x:?} {want:?}");
    }

    #[test]
    fn a_sphere_keeps_its_boundary_on_the_sphere_to_second_order() {
        let (p, t) = uv_sphere(10, 20);
        let p: Vec<[f64; 3]> = p.iter().map(|c| c.map(|x| 5.0 * x)).collect();
        let s = surface(&p, &t, [0.0; 3]);
        let faces = vec![0u32; s.triangles().len()];
        let s = BoundarySurface::new(s.points(), s.triangles(), &faces, &SurfaceOptions::default()).unwrap();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let sphere = vec![Some(AnalyticSurface::Sphere { center: [0.0; 3], radius: 5.0 })];
        let o = QualityOptions { size: Some(1.5), order: Order::Quadratic, ..QualityOptions::default() };
        let m = quality_mesh(&s, &cdt, sphere, &o).unwrap();
        assert_eq!(m.kind, fieldcraft_model::ElementKind::Tet10);
        assert_eq!(m.report.straightened, 0);
        // Every node of a boundary face, corners and mid-edge nodes, lies on the sphere.
        let mut on = 0;
        for b in &m.boundary {
            let e = &m.elements[b.tet as usize];
            let k = b.opposite as usize;
            let corners: Vec<usize> = (0..4).filter(|&i| i != k).collect();
            let mids =
                [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)].iter().enumerate().filter(|(_, (i, j))| *i != k && *j != k).map(|(n, _)| 4 + n);
            for i in corners.into_iter().chain(mids) {
                let q = m.points[e[i] as usize];
                let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
                assert!((r - 5.0).abs() < 1e-9, "a boundary node {i} of {e:?} at {q:?}, radius {r}; report {:?}", m.report);
                on += 1;
            }
        }
        assert!(on > 0);
        assert!(m.report.output.min_dihedral > 15.0, "{:?}", m.report.output);
        // Every element passes the certified check of the gate.
        for e in &m.elements {
            let xyz: Vec<[f64; 3]> = e.iter().map(|&v| m.points[v as usize]).collect();
            let v = crate::quality::element_validity(fieldcraft_model::ElementKind::Tet10, &xyz).unwrap();
            assert_eq!(v.validity, crate::quality::Validity::Valid);
        }
    }

    #[test]
    fn caps_stop_the_stage_without_panicking() {
        let s = box_surface();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        for o in [
            QualityOptions { size: Some(1e-3), max_points: 200, ..QualityOptions::default() },
            QualityOptions { size: Some(0.5), max_tets: 500, ..QualityOptions::default() },
            QualityOptions { size: Some(f64::NAN), ..QualityOptions::default() },
            QualityOptions { size: Some(-1.0), ..QualityOptions::default() },
        ] {
            match quality_mesh(&s, &cdt, Vec::new(), &o) {
                Ok(m) => assert!(!m.report.stopped.is_empty(), "{:?}", m.report),
                Err(e) => assert!(matches!(e.kind(), Some(crate::Kind::Limit | crate::Kind::Invalid)), "{e}"),
            }
        }
    }
}

/// Determinism (#115): the quality mesh is a function of the geometry, the same bit for bit on
/// every platform, run and listing order of the input.
mod determinism {
    use crate::cdt::{CdtOptions, tetrahedralize};
    use crate::surface::{BoundarySurface, SurfaceOptions};
    use crate::tet::canonical::CanonicalSoup;
    use crate::tet::{Order, QualityMesh, QualityOptions, quality_mesh};

    /// cos and sin of k · 15° from square roots only (correctly rounded everywhere), so the
    /// test input itself is the same on every platform.
    fn unit_24(k: usize) -> [f64; 2] {
        let (s2, s3, s6) = (2f64.sqrt(), 3f64.sqrt(), 6f64.sqrt());
        let c = [1.0, (s6 + s2) / 4.0, s3 / 2.0, s2 / 2.0, 0.5, (s6 - s2) / 4.0, 0.0];
        // Quadrant q, step j: cos(90° q + 15° j).
        let (q, j) = (k / 6 % 4, k % 6);
        let (x, y) = (c[j], c[6 - j]);
        match q {
            0 => [x, y],
            1 => [-y, x],
            2 => [-x, -y],
            _ => [y, -x],
        }
    }

    /// A washer (inner radius 6, outer 15, height 4) as a 24-sided prism ring: faces bottom 0,
    /// top 1, inside 2, outside 3; triangles oriented outward.
    pub(crate) fn washer() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<u32>) {
        const N: usize = 24;
        let mut p = Vec::new();
        for (r, z) in [(6.0, 0.0), (15.0, 0.0), (6.0, 4.0), (15.0, 4.0)] {
            for k in 0..N {
                let [c, s] = unit_24(k);
                p.push([r * c, r * s, z]);
            }
        }
        let id = |ring: usize, k: usize| (ring * N + k % N) as u32;
        let (mut t, mut f) = (Vec::new(), Vec::new());
        for k in 0..N {
            // (ring a, ring b): quad a_k b_k b_k+1 a_k+1, outward when listed this way.
            for (a, b, tag) in [(1, 0, 0), (2, 3, 1), (0, 2, 2), (3, 1, 3)] {
                let q = [id(a, k), id(b, k), id(b, k + 1), id(a, k + 1)];
                t.push([q[0], q[1], q[2]]);
                t.push([q[0], q[2], q[3]]);
                f.extend([tag, tag]);
            }
        }
        (p, t, f)
    }

    fn mesh(p: &[[f64; 3]], t: &[[u32; 3]], f: &[u32]) -> QualityMesh {
        let soup = CanonicalSoup::new(p, t, f).unwrap();
        let s = BoundarySurface::new(&soup.positions, &soup.triangles, &soup.tri_face, &SurfaceOptions::default()).unwrap();
        assert!(s.volume() > 0.0);
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let mut m =
            quality_mesh(&s, &cdt, Vec::new(), &QualityOptions { size: Some(3.0), order: Order::Quadratic, ..QualityOptions::default() }).unwrap();
        soup.restore_tags(&mut m);
        m
    }

    /// FNV-1a over every coordinate's bits, the elements and the tagged boundary faces.
    pub(crate) fn fingerprint(m: &QualityMesh) -> u64 {
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        let mut eat = |x: u64| {
            for b in x.to_le_bytes() {
                h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
            }
        };
        m.points.iter().flatten().for_each(|c| eat(c.to_bits()));
        m.elements.iter().flatten().for_each(|&v| eat(u64::from(v)));
        m.boundary.iter().for_each(|b| eat(u64::from(b.tet) << 40 | u64::from(b.opposite) << 32 | u64::from(b.face)));
        h
    }

    /// The washer's tet10 mesh has one fingerprint on every platform (CI checks it on Linux
    /// x86-64, macOS arm64 and Windows): only correctly rounded arithmetic and the crate's own
    /// `atan2` and `cbrt` make its decisions, and no hash map is iterated into it. A change that
    /// moves the tet mesh on purpose updates the constant.
    #[test]
    fn the_quality_mesh_is_the_same_on_every_platform() {
        let (p, t, f) = washer();
        let m = mesh(&p, &t, &f);
        assert!(m.elements.len() > 100, "{:?}", m.report);
        // No needles either: the scaled Jacobian analysts check is above the skill's limit.
        let sj = m
            .elements
            .iter()
            .map(|e| {
                let xyz: Vec<[f64; 3]> = e.iter().map(|&v| m.points[v as usize]).collect();
                crate::quality::element_quality(m.kind, &xyz).unwrap().get(crate::quality::Metric::ScaledJacobian).unwrap()
            })
            .fold(f64::INFINITY, f64::min);
        assert!(sj >= 0.3, "min scaled Jacobian {sj}");
        assert_eq!(fingerprint(&m), 11_252_377_805_886_396_373, "{} tets, {} nodes", m.elements.len(), m.points.len());
    }

    /// Listing the triangles, their corners and the vertices in another order, and numbering
    /// the faces differently, gives the same mesh with the faces' own tags.
    #[test]
    fn the_quality_mesh_does_not_depend_on_the_input_order() {
        let (p, t, f) = washer();
        let a = mesh(&p, &t, &f);
        let n = p.len() as u32;
        let q: Vec<[f64; 3]> = p.iter().rev().copied().collect();
        let u: Vec<[u32; 3]> = t.iter().rev().map(|x| [n - 1 - x[1], n - 1 - x[2], n - 1 - x[0]]).collect();
        let g: Vec<u32> = f.iter().rev().map(|&x| 10 + 3 - x).collect();
        let b = mesh(&q, &u, &g);
        assert_eq!((&a.points, &a.elements), (&b.points, &b.elements));
        assert_eq!(a.boundary.len(), b.boundary.len());
        for (x, y) in a.boundary.iter().zip(&b.boundary) {
            assert_eq!((x.tet, x.opposite, 10 + 3 - x.face), (y.tet, y.opposite, y.face));
        }
    }
}

/// A cap tetrahedron (two boundary faces of one B-rep face) whose two other neighbours each have
/// three boundary faces: deleting it would leave them joined to nothing by a face, floating free
/// (a singular stiffness: the dome-on-a-plate mesh on Linux, #438), so it stays.
#[test]
fn deleting_a_cap_never_leaves_a_neighbour_floating() {
    let p: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.3, 0.3, -1.0], [1.0, 1.0, 1.0]];
    let pts: Vec<Point3> = p.iter().map(|c| Point3::try_from(*c).unwrap()).collect();
    let positive = |t: [u32; 4]| if measure::volume6(&t.map(|v| pts[v as usize].coords())) > 0.0 { t } else { [t[1], t[0], t[2], t[3]] };
    // The cap, and a tetrahedron on each of its faces 012 and 123.
    let tets = [positive([0, 1, 2, 3]), positive([0, 1, 2, 4]), positive([1, 2, 3, 5])];
    let key = |f: [u32; 3]| {
        let mut f = f;
        f.sort_unstable();
        f
    };
    let mut boundary = Vec::new();
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            let f = key(face_toward(*vs, k));
            let shared = tets.iter().enumerate().any(|(u, ws)| u != t && (0..4).any(|j| key(face_toward(*ws, j)) == f));
            if !shared {
                boundary.push((t as u32, k, 0));
            }
        }
    }
    let mut cx = TetComplex::new(&p, vec![VertexKind::Corner; 6], &tets, &boundary, &[], 100, 100).unwrap();
    assert!(!cx.delete_boundary_tet(0).unwrap(), "the cap was deleted, leaving its neighbours floating");
    assert_eq!(cx.live_tets().count(), 3);
    cx.validate().unwrap();
}

/// A vertex inserted on a curved boundary face whose removed faces would ring a kept one: the
/// improvement's insertion on the rounded box's fillet in the tet10 run on Linux (SolidCraft
/// 0.18.1's tessellation). Vertex 10 tops a tetrahedron with three boundary faces; two of them
/// hold the new vertex in their circumcircle, the third does not, so the faces removed around it
/// form an annulus. Fanned from the new vertex, the inner loop and the kept face would enclose a
/// tetrahedron with four boundary faces, joined to the mesh by its vertices alone. The insertion
/// is refused instead.
#[test]
fn a_boundary_insertion_never_rings_a_kept_face() {
    let p: Vec<[f64; 3]> = vec![
        [37.62925187046595, 27.000000045728676, 19.925947437850475],
        [37.346388050558545, 26.568637357839815, 19.96821980678724],
        [37.864854899432075, 26.37493832376346, 19.858206230145555],
        [38.13849211461335, 26.78943127788964, 19.77164302800985],
        [37.141787936642764, 26.70315679040015, 19.514142391752234],
        [37.20584345508358, 26.832275662291323, 19.42626772210176],
        [37.48739155401909, 26.567563748691178, 19.955283105884796],
        [37.92769917669904, 26.41671699947128, 19.453486681697488],
        [37.43800406478175, 26.208143082650103, 19.55074124168729],
        [37.07336596672176, 26.83524701335969, 19.993268865411512],
        [37.58907678923575, 26.604654649192828, 19.937498683663613],
    ];
    let tets = [[3, 8, 6, 9], [3, 9, 7, 8], [8, 7, 5, 9], [3, 9, 0, 7], [9, 5, 0, 7], [8, 5, 4, 9], [1, 0, 3, 9], [1, 3, 0, 2], [0, 1, 2, 10]];
    // The faces on the fillet (tag 20); the others were shared with the rest of the mesh.
    let fillet: [&[usize]; 9] = [&[1], &[], &[], &[], &[], &[], &[1, 2], &[0, 2], &[0, 1, 2]];
    let key = |f: [u32; 3]| {
        let mut f = f;
        f.sort_unstable();
        f
    };
    let mut boundary = Vec::new();
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            let f = key(face_toward(*vs, k));
            let shared = tets.iter().enumerate().any(|(u, ws)| u != t && (0..4).any(|j| key(face_toward(*ws, j)) == f));
            if !shared {
                boundary.push((t as u32, k, if fillet[t].contains(&k) { 20 } else { 0 }));
            }
        }
    }
    let on_fillet = [false, true, true, true, false, false, true, false, false, true, true];
    let kinds: Vec<VertexKind> = on_fillet.iter().map(|&f| if f { VertexKind::Face(20) } else { VertexKind::Corner }).collect();
    let mut cx = TetComplex::new(&p, kinds, &tets, &boundary, &[], 100, 100).unwrap();
    // On the fillet face of tetrahedron 0, as far from the other vertices as the improvement asks.
    let q = Point3::try_from([37.57631755601779, 26.730747346632697, 19.941167252265036]).unwrap();
    insert_delaunay_apart(&mut cx, q, VertexKind::Face(20), Site::Face(0, 1), 4096, 0.06539086641130744).unwrap();
    cx.validate().unwrap();
    let floating: Vec<[u32; 4]> = cx.live_tets().filter(|(t, _)| cx.adj(*t).unwrap().iter().all(|&n| n == NONE)).map(|(_, vs)| vs).collect();
    assert!(floating.is_empty(), "tetrahedra with every face on the boundary: {floating:?}");
}

/// The validity gate refuses a tetrahedron that shares no face with the rest of the mesh, even
/// when the volume adds up: it says why instead of leaving the solver to find it singular.
#[test]
fn the_gate_refuses_a_detached_tetrahedron() {
    let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]];
    let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
    let v = 2.0 / 6.0;
    let e = super::linear_gate(&p, &tets, v).unwrap_err().to_string();
    assert!(e.contains("share no face") && e.contains("(1, 2)"), "{e}");
    assert!(super::linear_gate(&p, &tets[..1], v / 2.0).is_ok());
}
