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

/// A complex of these tetrahedra (each turned positive) and vertex kinds, every face not shared on
/// the boundary with the tag `tag` gives it (by its sorted vertices).
fn complex_of(p: &[[f64; 3]], raw: &[[u32; 4]], kinds: Vec<VertexKind>, tag: impl Fn([u32; 3]) -> u32) -> TetComplex {
    let tets: Vec<[u32; 4]> = raw
        .iter()
        .map(|t| {
            let x = t.map(|v| p[v as usize]);
            if measure::volume6(&x) > 0.0 { *t } else { [t[1], t[0], t[2], t[3]] }
        })
        .collect();
    let mut boundary = Vec::new();
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            let f = sorted3(face_toward(*vs, k));
            if !tets.iter().enumerate().any(|(u, ws)| u != t && (0..4).any(|j| sorted3(face_toward(*ws, j)) == f)) {
                boundary.push((t as u32, k, tag(f)));
            }
        }
    }
    TetComplex::new(p, kinds, &tets, &boundary, &[], 100, 100).unwrap()
}

/// The link condition (Dey, Edelsbrunner, Guha, Nekhayev 1999). The cube's boundary edge 0–1 may
/// be contracted; its inner diagonal 0–7 may not, as both ends see the boundary and the edge does
/// not. Round an edge of three tetrahedra, the triangle of the ring is empty (its edges are in the
/// mesh, the face is not): contracting one of its edges would lay that edge on another, so the
/// third vertex, seen by both ends off the edge's ring, refuses it.
#[test]
fn the_link_condition_refuses_contractions_that_pinch() {
    let mut cx = cube();
    assert!(link_condition(&mut cx, 0, 1).unwrap());
    assert!(!link_condition(&mut cx, 0, 7).unwrap());
    let p = [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [-0.5, 0.8, 0.0], [-0.5, -0.8, 0.0]];
    let mut cx = complex_of(&p, &[[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 2]], vec![VertexKind::Corner; 5], |_| 0);
    assert!(!link_condition(&mut cx, 2, 3).unwrap());
    assert!(link_condition(&mut cx, 0, 2).unwrap());
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
        let (mut cx, curve_faces) = complex_from_cdt(&s, &cdt, &[], &QualityOptions::default()).unwrap();
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

    /// An exact seed keeps its count beside a finer seed, which grades the size field down at
    /// their shared corner (refinement split the coarse seed's end segments there, #894); a
    /// `min` seed there may gain divisions but keeps its nodes.
    #[test]
    fn an_exact_seed_keeps_its_count_beside_a_finer_one() {
        use crate::seeds::EdgeSeed;
        use fieldcraft_model::SeedConstraint;
        let s = box_surface();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let seed = |name: &str, to: [f64; 3], count, constraint| EdgeSeed {
            name: name.into(),
            points: vec![[0.0; 3], to],
            tolerance: 0.0,
            count,
            bias: None,
            constraint,
        };
        let fine = seed("y", [0.0, 10.0, 0.0], 20, SeedConstraint::Exact);
        let exact = QualityOptions {
            size: Some(5.0),
            seeds: vec![seed("x", [10.0, 0.0, 0.0], 3, SeedConstraint::Exact), fine.clone()],
            ..QualityOptions::default()
        };
        let m = quality_mesh(&s, &cdt, Vec::new(), &exact).unwrap();
        let counts: Vec<usize> = m.report.seeds.iter().map(|s| s.2).collect();
        assert_eq!(counts, [3, 20]);
        let mut x: Vec<f64> = m.points.iter().filter(|p| p[1].abs() < 1e-9 && p[2].abs() < 1e-9).map(|p| p[0]).collect();
        x.sort_by(f64::total_cmp);
        assert_eq!(x.len(), 4, "{x:?}");
        assert!(x.iter().zip([0.0, 10.0 / 3.0, 20.0 / 3.0, 10.0]).all(|(a, b)| (a - b).abs() < 1e-9), "{x:?}");
        // The seeds given the other way round: the same mesh.
        let swapped = QualityOptions { seeds: exact.seeds.iter().rev().cloned().collect(), ..exact.clone() };
        let n = quality_mesh(&s, &cdt, Vec::new(), &swapped).unwrap();
        assert_eq!((&n.points, &n.elements), (&m.points, &m.elements));
        // At least three: the nodes are kept, and refinement may add more beside the fine seed.
        let min = QualityOptions { seeds: vec![seed("x", [10.0, 0.0, 0.0], 3, SeedConstraint::Min), fine], ..exact };
        let m = quality_mesh(&s, &cdt, Vec::new(), &min).unwrap();
        assert!(m.report.seeds[0].2 >= 3, "{:?}", m.report.seeds);
        for want in [10.0 / 3.0, 20.0 / 3.0] {
            assert!(m.points.iter().any(|p| (p[0] - want).abs() < 1e-9 && p[1].abs() < 1e-9 && p[2].abs() < 1e-9), "{want}");
        }
    }

    /// The turns between the triangles of a face with a known exact surface are the
    /// tessellation's, never feature curves: SolidCraft tessellates the trimmed sphere of the
    /// sphere-on-plate part with needles along the trim that face the other way, and the curves
    /// their folds made pinned needles the refinement could not remove (min dihedral 6.5°). A
    /// face without its exact surface (an imported mesh) keeps its creases.
    #[test]
    fn a_face_with_an_exact_surface_has_no_creases_inside() {
        use super::super::complex_from_cdt;
        let (p, t) = uv_sphere(4, 8);
        let p: Vec<[f64; 3]> = p.iter().map(|c| c.map(|x| 5.0 * x)).collect();
        let s = surface(&p, &t, [0.0; 3]);
        let faces = vec![0u32; s.triangles().len()];
        let s = BoundarySurface::new(s.points(), s.triangles(), &faces, &SurfaceOptions::default()).unwrap();
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let o = QualityOptions { feature_angle: 30.0, ..QualityOptions::default() };
        let (imported, curves) = complex_from_cdt(&s, &cdt, &[], &o).unwrap();
        assert!(!imported.curve_edges().is_empty() && !curves.is_empty());
        let sphere = [Some(AnalyticSurface::Sphere { center: [0.0; 3], radius: 5.0 })];
        let (exact, curves) = complex_from_cdt(&s, &cdt, &sphere, &o).unwrap();
        assert!(exact.curve_edges().is_empty() && curves.is_empty(), "{curves:?}");
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

    /// Refined far below its tessellation's facets (8 facets round, sag 0.23, elements 0.5), a
    /// cylinder's boundary nodes lie on the exact cylinder and caps, mid-side nodes included:
    /// the nodes refinement adds on a facet move out together (#674).
    #[test]
    fn nodes_added_below_the_facet_size_lie_on_the_exact_faces() {
        let (r, h) = (3.0, 1.5);
        let s = crate::octree::tests::cylinder_surface(r, h, 8);
        let cdt = tetrahedralize(&s, &CdtOptions::default()).unwrap();
        let analytic = vec![
            Some(AnalyticSurface::Plane { normal: [0.0, 0.0, -1.0], point: [0.0; 3] }),
            Some(AnalyticSurface::Plane { normal: [0.0, 0.0, 1.0], point: [0.0, 0.0, h] }),
            Some(AnalyticSurface::Cylinder { origin: [0.0; 3], axis: [0.0, 0.0, 1.0], radius: r }),
        ];
        let o = QualityOptions { size: Some(0.5), order: Order::Quadratic, ..QualityOptions::default() };
        let m = quality_mesh(&s, &cdt, analytic, &o).unwrap();
        assert_eq!((m.report.unsnapped, m.report.snap_gap), (0, 0.0), "{:?}", m.report);
        let mut side = 0;
        for b in &m.boundary {
            let e = &m.elements[b.tet as usize];
            let k = b.opposite as usize;
            let mids = crate::quality::jacobian::TET_EDGES.iter().enumerate().filter(|(_, (i, j))| *i != k && *j != k);
            for i in (0..4).filter(|&i| i != k).chain(mids.map(|(n, _)| 4 + n)) {
                let q = m.points[e[i] as usize];
                let gap = match b.face {
                    0 => q[2].abs(),
                    1 => (q[2] - h).abs(),
                    _ => (q[0].hypot(q[1]) - r).abs(),
                };
                assert!(gap < 1e-9, "node {i} of {e:?} at {q:?} is {gap} off face {}; report {:?}", b.face, m.report);
                side += usize::from(b.face == 2);
            }
        }
        assert!(side > 500, "{side}");
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

/// A cap tetrahedron whose two other faces meet at an edge that is on the boundary already:
/// deleting it would put that edge in four boundary faces, the groups of tetrahedra on either
/// side joined by the edge alone (#671), so it stays. Each group has two tetrahedra, so neither
/// neighbour would float on its own.
#[test]
fn deleting_a_cap_never_pinches_the_mesh_at_an_edge() {
    // The cap 0-1-2-3: boundary faces 1-2-3 and 0-2-3, inner faces 0-1-2 and 0-1-3 over the
    // edge 0-1. Beyond each inner face a pair of tetrahedra, open to each other round 0-1.
    let p =
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5, 1.0], [0.5, 0.5, -1.0], [0.5, -0.3, -1.0], [0.5, -0.3, 1.0], [1.5, 0.1, -1.0], [1.5, 0.1, 1.0]];
    let tets = [[0, 1, 2, 3], [0, 1, 3, 4], [1, 3, 4, 6], [0, 1, 2, 5], [1, 2, 5, 7]];
    let mut cx = complex_of(&p, &tets, vec![VertexKind::Face(0); 8], |_| 0);
    assert!(cx.edge_ring(0, 1).unwrap().is_some_and(|r| !r.closed));
    assert!(!cx.delete_boundary_tet(0).unwrap(), "the cap was deleted, pinching the mesh at the edge 0-1");
    assert_eq!(cx.live_tets().count(), 5);
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

/// A move keeps each boundary face round the vertex within 30° of its orientation: positive
/// tetrahedra alone let a boundary vertex fold its faces over the ones beside it. Corner 1 of
/// the cube moved out to x = 3 leaves both its tetrahedra positive but turns its faces on the
/// side x = 1 by 63°; moved to x = 1.1, by 6°.
#[test]
fn a_move_never_turns_a_boundary_face_over() {
    let mut cx = cube();
    let to = |p: [f64; 3]| Point3::try_from(p).unwrap();
    let before = format!("{cx:?}");
    assert!(!cx.try_move(1, to([3.0, 0.0, 0.0])).unwrap());
    assert_eq!(format!("{cx:?}"), before);
    assert!(cx.try_move(1, to([1.1, 0.0, 0.0])).unwrap());
}

#[test]
fn checked_moves_record_only_commits_and_can_be_rolled_back() {
    let mut cx = cube();
    cx.index_skin(1.0).unwrap();
    cx.journal.recording = true;
    let before = format!("{cx:?}");
    let mesh = format!("{:?}", cx.compact().unwrap());
    let from = cx.point(1).unwrap();
    assert!(!cx.try_move(1, from).unwrap());
    assert!(cx.try_move(NONE, from).is_err());
    assert_eq!(format!("{cx:?}"), before);

    let mut unused = cx.clone();
    let v = unused.add_point(Point3::try_from([2.0; 3]).unwrap(), VertexKind::Interior).unwrap();
    let unchanged = format!("{unused:?}");
    assert!(unused.try_move(v, Point3::try_from([2.1; 3]).unwrap()).is_err());
    assert!(unused.try_move(v, unused.point(v).unwrap()).is_err());
    assert_eq!(format!("{unused:?}"), unchanged);

    let to = Point3::try_from([1.1, 0.0, 0.0]).unwrap();
    assert!(cx.try_move(1, to).unwrap());
    assert_eq!(cx.xyz(1).unwrap(), to.coords());
    assert!(matches!(cx.journal.steps.as_slice(), [Step::Moved { v: 1, from: p }] if p.coords() == from.coords()));
    cx.validate().unwrap();
    cx.rollback(0).unwrap();
    cx.validate().unwrap();
    assert!(cx.journal.steps.is_empty() && cx.journal.recording);
    assert_eq!(format!("{:?}", cx.compact().unwrap()), mesh);
    let restored = format!("{cx:?}");
    assert!(!cx.try_move(1, Point3::try_from([3.0, 0.0, 0.0]).unwrap()).unwrap());
    assert_eq!(format!("{cx:?}"), restored);
}

/// A supplied partial star could miss the small tetrahedron on the other side of vertex 0.
/// The committing API derives its own full star, even when the pieces meet only at the vertex.
#[test]
fn a_checked_move_cannot_omit_a_piece_of_the_star() {
    let mut cx = complex_of(
        &[[0.0; 3], [0.125, 0.0, 0.0], [0.0, 0.125, 0.0], [0.0, 0.0, 0.125], [-10.0, 0.0, 0.0], [0.0, -10.0, 0.0], [0.0, 0.0, -10.0]],
        &[[0, 1, 2, 3], [0, 4, 5, 6]],
        vec![VertexKind::Corner; 7],
        |_| 0,
    );
    cx.journal.recording = true;
    let p = Point3::try_from([0.0625; 3]).unwrap();
    assert!(cx.move_keeps_star(0, &p, &[1]).unwrap());
    let star = cx.star(0).unwrap();
    assert!(!cx.move_keeps_star(0, &p, &star).unwrap());
    let before = format!("{cx:?}");
    assert!(!cx.try_move(0, p).unwrap());
    assert_eq!(format!("{cx:?}"), before);
}

/// A snap may satisfy every geometric guard while reducing an already poor star's quality.
/// Reject it before committing: moving back leaves spurious journal/index entries behind.
#[test]
fn a_rejected_snap_never_commits_or_journals() {
    use super::geometry::BoundaryGeometry;
    use super::refine::Context;
    use crate::size::UniformSize;

    struct Plane;
    impl BoundaryGeometry for Plane {
        fn project_face(&self, _: u32, p: [f64; 3]) -> Option<[f64; 3]> {
            Some([p[0], p[1], 0.005])
        }
        fn project_curve(&self, _: [u32; 2], _: [f64; 3]) -> Option<[f64; 3]> {
            None
        }
    }
    let mut cx = complex_of(
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.25, 0.25, 0.01]],
        &[[0, 1, 2, 3]],
        vec![VertexKind::Corner, VertexKind::Corner, VertexKind::Corner, VertexKind::Face(0)],
        |_| 0,
    );
    cx.index_skin(1.0).unwrap();
    cx.journal.recording = true;
    let star = cx.star(3).unwrap();
    let p = Point3::try_from([0.25, 0.25, 0.005]).unwrap();
    assert!(cx.move_keeps_star(3, &p, &star).unwrap());
    assert!(cx.skin_allows(3, &p, &star).unwrap());
    let vs = cx.tet(0).unwrap();
    let trial = vs.map(|v| if v == 3 { p.coords() } else { cx.xyz(v).unwrap() });
    assert!(measure::quality(&trial) < worst(&cx, &star).unwrap().min(0.3));
    let before = format!("{cx:?}");
    let ctx = Context { size: &UniformSize(1.0), geometry: &Plane, curve_faces: &[] };
    assert_eq!(super::snap_to_geometry(&mut cx, &ctx, true).unwrap().off, 1);
    assert_eq!(format!("{cx:?}"), before);
}

/// Sweeps that still move vertices when they run out say so in the report's reasons, as refine
/// and improve do. Here the face is the vertex's mirror image through its tetrahedron's base, so
/// the whole move and half of it (onto the base) are refused and every sweep moves it a quarter of
/// the way: it halves its height each time and is still off after the last.
#[test]
fn snapping_reports_running_out_of_sweeps() {
    use super::geometry::BoundaryGeometry;
    use super::refine::Context;
    use crate::size::UniformSize;

    struct Mirror;
    impl BoundaryGeometry for Mirror {
        fn project_face(&self, _: u32, p: [f64; 3]) -> Option<[f64; 3]> {
            Some([p[0], p[1], -p[2]])
        }
        fn project_curve(&self, _: [u32; 2], _: [f64; 3]) -> Option<[f64; 3]> {
            None
        }
    }
    let mut cx = complex_of(
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.25, 0.25, 0.5]],
        &[[0, 1, 2, 3]],
        vec![VertexKind::Corner, VertexKind::Corner, VertexKind::Corner, VertexKind::Face(0)],
        |_| 0,
    );
    cx.index_skin(1.0).unwrap();
    let ctx = Context { size: &UniformSize(1.0), geometry: &Mirror, curve_faces: &[] };
    let snapped = super::snap_to_geometry(&mut cx, &ctx, false).unwrap();
    assert_eq!(snapped.off, 1);
    assert!(snapped.stopped.as_ref().is_some_and(|s| s.contains(&format!("after {} sweeps", super::SNAP_SWEEPS))));
    let z = 0.5 * 0.5_f64.powi(super::SNAP_SWEEPS as i32);
    assert!((cx.xyz(3).unwrap()[2] - z).abs() < 1e-9 * z, "{:?}", cx.xyz(3));
    assert!((snapped.gap - 2.0 * z).abs() < 1e-9 * z, "{snapped:?}");
}

/// With the skin indexed, a move or a change whose boundary faces would pass through another part
/// of the skin is refused, however positive its tetrahedra. A square pyramid, its apex 0 in its
/// top face, under a separate tetrahedron whose bottom face lies 0.1 above: lifting the apex 0.2
/// keeps the pyramid positive and turns its top faces by 11°, but pushes them through that face,
/// whether the apex moves or a change puts a new vertex there; lifting it 0.05 does not.
#[test]
fn the_skin_never_passes_through_itself() {
    let p = [
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, -1.0, 0.0],
        [0.0, 0.0, -1.0],
        [-1.0, -1.0, 0.1],
        [1.0, -1.0, 0.1],
        [0.0, 1.0, 0.1],
        [0.0, 0.0, 1.0],
    ];
    let pyramid = [[0, 1, 2, 5], [0, 2, 3, 5], [0, 3, 4, 5], [0, 4, 1, 5]];
    let mut tets = pyramid.to_vec();
    tets.push([6, 7, 8, 9]);
    let plain = complex_of(&p, &tets, vec![VertexKind::Corner; p.len()], |_| 0);
    let mut indexed = plain.clone();
    indexed.index_skin(1.0).unwrap();
    let to = |z: f64| Point3::try_from([0.0, 0.0, z]).unwrap();
    for (mut cx, skin) in [(plain, false), (indexed, true)] {
        assert!(cx.clone().try_move(0, to(0.05)).unwrap());
        let mut trial = cx.clone();
        trial.journal.recording = true;
        let before = format!("{trial:?}");
        assert_eq!(trial.try_move(0, to(0.2)).unwrap(), !skin);
        if skin {
            assert_eq!(format!("{trial:?}"), before);
        }
        // The same lift as a change: a new apex 10, the old one gone.
        let v = cx.add_point(to(0.2), VertexKind::Face(0)).unwrap();
        let mut change = Change { old: (0..4).collect(), ..Change::default() };
        for t in pyramid {
            let x = [v, t[1], t[2], t[3]];
            let x = if measure::volume6(&cx.coords(x).unwrap()) > 0.0 { x } else { [x[1], x[0], x[2], x[3]] };
            let top = x.iter().position(|&u| u == 5).unwrap();
            change.new_boundary.push((face_toward(x, top), 0));
            change.new.push(x);
        }
        assert_eq!(cx.replace(&change).unwrap().is_some(), !skin);
    }
}

/// The skin index's cell sets only how much a check reads, never what it finds: the same moves on
/// two blocks of cubes stacked with a narrow gap (a vertex of the lower block's top lifted into
/// the upper block passes through its skin without inverting anything) pass or fail alike with
/// one cell for everything, cells of the median boundary face, and cells so small that every face
/// is read by every check.
#[test]
fn the_skin_index_cell_does_not_change_a_check() {
    let (nx, nz) = (3, 2);
    let mut p = Vec::new();
    let mut raw = Vec::new();
    for z0 in [0.0, 0.55] {
        let base = p.len();
        let id = |x: usize, y: usize, z: usize| (base + x + (nx + 1) * (y + (nx + 1) * z)) as u32;
        for z in 0..=nz {
            for y in 0..=nx {
                for x in 0..=nx {
                    p.push([x as f64 / nx as f64, y as f64 / nx as f64, z0 + 0.25 * z as f64]);
                }
            }
        }
        // Six tetrahedra round each cube's diagonal from its corner 0 to 7.
        for z in 0..nz {
            for y in 0..nx {
                for x in 0..nx {
                    let c = |i: usize| id(x + (i & 1), y + ((i >> 1) & 1), z + (i >> 2));
                    for [a, b] in [[1, 3], [3, 2], [2, 6], [6, 4], [4, 5], [5, 1]] {
                        raw.push([c(0), c(a), c(b), c(7)]);
                    }
                }
            }
        }
    }
    let tets: Vec<[u32; 4]> =
        raw.iter().map(|t| if measure::volume6(&t.map(|v| p[v as usize])) > 0.0 { *t } else { [t[1], t[0], t[2], t[3]] }).collect();
    let mut faces: std::collections::HashMap<[u32; 3], usize> = std::collections::HashMap::new();
    for vs in &tets {
        for k in 0..4 {
            *faces.entry(sorted3(face_toward(*vs, k))).or_default() += 1;
        }
    }
    let mut boundary = Vec::new();
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            if faces.get(&sorted3(face_toward(*vs, k))) == Some(&1) {
                boundary.push((t as u32, k, 0));
            }
        }
    }
    let plain = TetComplex::new(&p, vec![VertexKind::Corner; p.len()], &tets, &boundary, &[], 1000, 10_000).unwrap();
    // The lower block's top vertices, lifted by 0.05 to 0.3: the higher ones pass into the upper
    // block, and only the skin check stops them.
    let top: Vec<u32> = (0..p.len() as u32).filter(|&v| p[v as usize][2] == 0.5).collect();
    let mut runs = Vec::new();
    for cell in [1e6, 1.0, 1e-6] {
        let mut cx = plain.clone();
        cx.index_skin(cell).unwrap();
        let mut moved = Vec::new();
        for (i, &v) in top.iter().enumerate() {
            for lift in [0.05, 0.1, 0.2, 0.3] {
                let at = cx.point(v).unwrap().coords();
                let to = Point3::try_from([at[0] + 0.01 * i as f64, at[1], at[2] + lift]).unwrap();
                moved.push(cx.clone().try_move(v, to).unwrap());
            }
        }
        runs.push(moved);
    }
    assert!(runs[0].iter().any(|&m| m) && runs[0].iter().any(|&m| !m), "{:?}", runs[0]);
    assert_eq!(runs[0], runs[1]);
    assert_eq!(runs[0], runs[2]);
}

/// The fillet of a plate on a rib (our own recipe: a 60 × 30 × 6 mm plate on a 6 mm triangular
/// rib, its top edges rounded at 2 mm, meshed at 3 mm), where vertex 1 has three boundary faces
/// round it: (0, 1, 2), (0, 1, 5) and (1, 2, 5). Removing edge 0–1 ends in a 2-2 flip of the
/// boundary onto the edge 2–5, which the mesh has already: one of the new boundary faces is
/// (1, 2, 5), held by tetrahedron 0 too. The fillet's curvature lets every orientation pass. The
/// removal must make no face the mesh has.
#[test]
fn an_edge_removal_never_flips_the_boundary_onto_a_face_it_has() {
    let p = [
        [42.66838317202295, 11.015049414715271, 39.676406229007135],
        [42.18976625943718, 10.738523431897063, 39.858043495916334],
        [42.144070440995236, 11.092938909632773, 39.62524418327976],
        [41.984755325992694, 10.464231377518912, 39.91145084982469],
        [42.37942326300024, 11.234853867339881, 40.483532518490776],
        [42.52174059791378, 10.687602908819951, 39.86884416654581],
    ];
    let tets = [[1, 3, 2, 5], [4, 1, 0, 5], [1, 0, 2, 4]];
    // The faces on the fillet (tag 1); the others were shared with the rest of the mesh.
    let fillet = [[1, 2, 5], [2, 3, 5], [0, 1, 5], [0, 1, 2]];
    let kinds = [1, 1, 1, 1, 0, 1].map(|f| if f == 1 { VertexKind::Face(1) } else { VertexKind::Interior }).to_vec();
    let mut cx = complex_of(&p, &tets, kinds, |f| u32::from(fillet.contains(&f)));
    remove_edge_at(&mut cx, 2, 0, 1).unwrap();
    cx.validate().unwrap();
    let mut held: std::collections::BTreeMap<[u32; 3], (u32, bool)> = std::collections::BTreeMap::new();
    for (t, vs) in cx.live_tets() {
        let adj = cx.adj(t).unwrap();
        for k in 0..4 {
            let e = held.entry(sorted3(face_toward(vs, k))).or_default();
            e.0 += 1;
            e.1 |= adj[k] == NONE;
        }
    }
    let folded: Vec<_> = held.into_iter().filter(|(_, (n, on_boundary))| *n > 2 || (*n == 2 && *on_boundary)).collect();
    assert!(folded.is_empty(), "faces held twice over: {folded:?}");
}

/// The validity gate refuses a tetrahedron that shares no face with the rest of the mesh, even
/// when the volume adds up: it says why instead of leaving the solver to find it singular.
#[test]
fn the_gate_refuses_a_detached_tetrahedron() {
    let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]];
    let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
    let v = 2.0 / 6.0;
    let e = super::linear_gate(&p, &tets, &faces_but(&tets, &[]), v).unwrap_err().to_string();
    assert!(e.contains("share no face") && e.contains("(1, 2)"), "{e}");
    assert!(super::linear_gate(&p, &tets[..1], &faces_but(&tets[..1], &[]), v / 2.0).is_ok());
}

/// `tets` with each one turned positive (two corners swapped where its volume is negative).
fn positive(p: &[[f64; 3]], tets: &[[u32; 4]]) -> Vec<[u32; 4]> {
    tets.iter().map(|&t| if measure::volume6(&t.map(|v| p[v as usize])) > 0.0 { t } else { [t[1], t[0], t[2], t[3]] }).collect()
}

/// The validity gate refuses two tetrahedra that share a face with each other but meet the rest
/// of the mesh at one edge alone (#671): neither is detached on its own, and the volume adds up,
/// so only the pinched boundary shows it. The same pair meeting the body at one vertex is refused
/// too.
#[test]
fn the_gate_refuses_tetrahedra_joined_by_an_edge_or_a_vertex() {
    // The body: two tetrahedra on the face 0-1-2. The pair beyond the edge 1-2 (x + y > 1).
    let p = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 1.0, 1.0],
        [1.0, 1.0, -1.0],
        [0.7, 0.7, -1.0],
        [2.0, 0.0, 1.0],
        [2.0, 0.0, -1.0],
    ];
    let body = [[0, 1, 2, 3], [0, 1, 2, 4]];
    // The faces two tetrahedra hold.
    let inner = [[0, 1, 2], [1, 2, 6], [1, 8, 9]];
    let gate = |tets: &[[u32; 4]]| {
        let boundary = faces_but(tets, &inner);
        super::linear_gate(&p, tets, &boundary, super::boundary_volume(&p, tets, &boundary).unwrap())
    };
    let body = positive(&p, &body);
    assert!(gate(&body).is_ok());
    // Joined by the edge 1-2: the two share the face 1-2-6.
    let edge = positive(&p, &[body[0], body[1], [1, 2, 5, 6], [1, 2, 6, 7]]);
    let e = gate(&edge).unwrap_err().to_string();
    assert!(e.contains("1 boundary edges are not in exactly two boundary faces (the first, 1–2,"), "{e}");
    assert!(super::detached_tets(&edge, &super::face_use(&edge)).is_empty());
    // Joined by the vertex 1: the two share the face 1-8-9.
    let vertex = positive(&p, &[body[0], body[1], [1, 8, 9, 5], [1, 8, 9, 6]]);
    let e = gate(&vertex).unwrap_err().to_string();
    assert!(e.contains("the mesh is pinched at 1 vertices (the first, 1,"), "{e}");
    assert!(super::detached_tets(&vertex, &super::face_use(&vertex)).is_empty());
}

/// The validity gate refuses a mesh folded over itself even when the fold's two layers add up to
/// the volume the boundary encloses: a face held by three tetrahedra.
#[test]
fn the_gate_refuses_a_face_of_three_tetrahedra() {
    // Two tetrahedra on either side of the base triangle, and a third over it again.
    let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.2, 0.2, 1.0], [0.2, 0.2, -1.0], [0.3, 0.3, 0.5]];
    let tets = [[0, 1, 2, 3], [0, 2, 1, 4], [0, 1, 2, 5]];
    let v = (1.0 + 1.0 + 0.5) / 6.0;
    let e = super::linear_gate(&p, &tets, &faces_but(&tets, &[[0, 1, 2]]), v).unwrap_err().to_string();
    assert!(e.contains("held by more than two"), "{e}");
    assert!(super::linear_gate(&p, &tets[..2], &faces_but(&tets[..2], &[[0, 1, 2]]), 2.0 / 6.0).is_ok());
}

/// The validity gate refuses a face listed as boundary that two tetrahedra hold, one on either
/// side: the mesh closed onto its own boundary, as across a narrow gap. The volume adds up and no
/// face is held by three, so only the listed boundary shows it.
#[test]
fn the_gate_refuses_a_boundary_face_of_two_tetrahedra() {
    let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.2, 0.2, 1.0], [0.2, 0.2, -1.0]];
    let tets = [[0, 1, 2, 3], [0, 2, 1, 4]];
    let v = 2.0 / 6.0;
    let e = super::linear_gate(&p, &tets, &faces_but(&tets, &[]), v).unwrap_err().to_string();
    assert!(e.contains("1 boundary faces are held by two tetrahedra"), "{e}");
    assert!(super::linear_gate(&p, &tets, &faces_but(&tets, &[[0, 1, 2]]), v).is_ok());
    // A face held once but left out of the list, or listed twice.
    let mut short = faces_but(&tets, &[[0, 1, 2]]);
    let last = short.pop().unwrap();
    assert!(super::linear_gate(&p, &tets, &short, v).is_err());
    short.extend([last, last]);
    assert!(super::linear_gate(&p, &tets, &short, v).is_err());
}

/// Every face of `tets` except `inner`, as a boundary list.
fn faces_but(tets: &[[u32; 4]], inner: &[[u32; 3]]) -> Vec<super::TaggedFace> {
    let mut out = Vec::new();
    for (t, vs) in tets.iter().enumerate() {
        for k in 0..4 {
            if !inner.contains(&sorted3(face_toward(*vs, k))) {
                out.push(super::TaggedFace { tet: t as u32, opposite: k as u8, face: 0 });
            }
        }
    }
    out
}

/// A needle of a sphere's tessellation, its corners nearly along a meridian and ordered to face
/// into the sphere, is reported with its face; the triangles that face out are not, and a face
/// without an exact surface is passed over.
#[test]
fn a_facet_facing_into_its_face_is_reported() {
    use super::geometry::{AnalyticSurface, reversed_facets};
    let r = 10.0;
    let on = |polar: f64, azimuth: f64| {
        let (p, a) = (polar.to_radians(), azimuth.to_radians());
        [r * p.sin() * a.cos(), r * p.sin() * a.sin(), r * p.cos()]
    };
    // A cap facing out (counter-clockwise seen from above), and the needle.
    let mut points = vec![on(15.0, 0.0), on(15.0, 120.0), on(15.0, 240.0)];
    points.extend([on(20.0, 0.0), on(25.0, 2.0), on(30.0, 0.0)]);
    let outward = |t: [u32; 3]| {
        let [a, b, c] = t.map(|v| points[v as usize]);
        let n = super::measure::cross(super::measure::sub(b, a), super::measure::sub(c, a));
        let m = super::measure::scale(super::measure::add(super::measure::add(a, b), c), 1.0 / 3.0);
        super::measure::dot(n, m) > 0.0
    };
    let needle = if outward([3, 4, 5]) { [3, 5, 4] } else { [3, 4, 5] };
    let triangles = vec![[0, 1, 2], needle];
    assert!(outward([0, 1, 2]) && !outward(needle));
    let sphere = [Some(AnalyticSurface::Sphere { center: [0.0; 3], radius: r })];
    let found = reversed_facets(&points, &triangles, &[0, 0], &sphere);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!((found[0].triangle, found[0].face), (1, 0));
    assert!(found[0].angle > 90.0 && found[0].angle <= 180.0, "{found:?}");
    // Without the exact surface there is nothing to compare with.
    assert!(reversed_facets(&points, &triangles, &[0, 0], &[None]).is_empty());
}
