use super::*;
use fieldcraft_model::Model;

/// The plate (40 × 30 × 10) with its Ø10 hole round (15, 15).
const HOLE: ([f64; 2], f64) = ([15.0, 15.0], 5.0);

fn plate() -> Session {
    let mut s = super::super::datum::tests::plate().0;
    s.execute("mesh.hex.generate", &json!({"size": 2.0})).unwrap();
    s
}

/// The model in millimetres, as the plate is drawn (it is stored in metres).
fn mm(m: &Model) -> Model {
    m.converted(&|q, v| crate::units::geometry().express(q, v).ok()).unwrap()
}

/// Node `n`'s position in millimetres.
fn at(s: &Session, n: u32) -> [f64; 3] {
    crate::units::to_geometry_point(s.model.node(n).unwrap().xyz)
}

fn hole_radius(p: [f64; 3]) -> f64 {
    ((p[0] - HOLE.0[0]).powi(2) + (p[1] - HOLE.0[1]).powi(2)).sqrt()
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
}

/// The node nearest `p` among those `keep` accepts.
fn nearest(m: &Model, p: [f64; 3], keep: impl Fn([f64; 3]) -> bool) -> u32 {
    let mut best = (f64::INFINITY, 0);
    for (n, x) in m.nodes() {
        if keep(x.xyz) && dist(x.xyz, p) < best.0 {
            best = (dist(x.xyz, p), *n);
        }
    }
    best.1
}

/// An interior node (off every face of the plate).
fn interior(m: &Model) -> u32 {
    nearest(m, [30.0, 15.0, 5.0], |p| p[2] > 1e-6 && p[2] < 10.0 - 1e-6 && p[0] > 1e-6 && p[0] < 40.0 - 1e-6 && p[1] > 1e-6 && p[1] < 30.0 - 1e-6)
}

fn invalid(s: &mut Session) -> u64 {
    s.execute("mesh.quality", &json!({})).unwrap()["invalid"]["count"].as_u64().unwrap_or(0)
}

/// `by` moves the grabbed node by the vector, the reply carries the quality before and after,
/// one undo step takes it back and redo does it again.
#[test]
fn by_moves_a_node_in_one_undo_step() {
    let mut s = plate();
    let start = s.model.clone();
    let n = interior(&mm(&s.model));
    let p0 = at(&s, n);
    let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0.3, -0.2, 0.1]})).unwrap();
    assert_eq!((r["moved"].as_u64(), r["scale"].as_f64(), r["grabbed"].as_u64()), (Some(1), Some(1.0), Some(1)), "{r}");
    assert!(dist(at(&s, n), [p0[0] + 0.3, p0[1] - 0.2, p0[2] + 0.1]) < 1e-12);
    assert!(r["before"]["min_scaled_jacobian"].as_f64().unwrap() > r["after"]["min_scaled_jacobian"].as_f64().unwrap(), "{r}");
    assert!(r["worse"].as_u64().unwrap() > 0 && r["before"]["mean_scaled_jacobian"].is_number(), "{r}");
    let moved = s.model.clone();
    s.execute("edit.undo", &json!({})).unwrap();
    assert_eq!(s.model, start);
    s.execute("edit.redo", &json!({})).unwrap();
    assert_eq!(s.model, moved);
}

/// `to` puts the grabbed nodes' centroid (or the sphere's centre) at the point.
#[test]
fn to_puts_the_centroid_at_the_point() {
    let mut s = plate();
    let n = interior(&mm(&s.model));
    let p0 = at(&s, n);
    let to = [p0[0] + 0.25, p0[1], p0[2] - 0.25];
    s.execute("mesh.node.move", &json!({"nodes": [n], "to": to})).unwrap();
    assert!(dist(at(&s, n), to) < 1e-12);
    let c = at(&s, n);
    let r = s.execute("mesh.node.move", &json!({"center": c, "radius": 3.0, "to": [c[0], c[1] + 0.2, c[2]]})).unwrap();
    assert!(r["moved"].as_u64().unwrap() > 1, "{r}");
    assert!((at(&s, n)[1] - c[1] - 0.2).abs() < 1e-9, "the node at the centre takes all of it");
}

/// Each falloff moves the nodes within the radius, none beyond, the grabbed node by all of it;
/// they move the followers differently.
#[test]
fn each_falloff_moves_only_within_the_radius() {
    let mut followers = Vec::new();
    for falloff in ["smooth", "linear", "none"] {
        let mut s = plate();
        let start = mm(&s.model);
        let n = interior(&mm(&s.model));
        let p0 = start.node(n).unwrap().xyz;
        let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0.0, 0.0, 0.2], "radius": 4.0, "falloff": falloff})).unwrap();
        assert!(r["moved"].as_u64().unwrap() > 1, "{falloff}: {r}");
        let mut sum = 0.0;
        let now = mm(&s.model);
        for (k, x) in now.nodes() {
            let old = start.node(*k).unwrap().xyz;
            if old != x.xyz {
                assert!(dist(old, p0) < 4.0 + 1e-9, "{falloff}: node {k} beyond the radius moved");
                sum += dist(old, x.xyz);
            }
        }
        assert!((at(&s, n)[2] - p0[2] - 0.2).abs() < 1e-12);
        followers.push(sum);
    }
    assert!(followers[0] < followers[1] && followers[1] < followers[2], "{followers:?}");
}

/// Harmonic: the nodes within the radius follow, the move falls to nothing at the radius and no
/// node beyond it moves; falloff and strength do not go with it.
#[test]
fn harmonic_propagation_is_zero_outside_the_radius() {
    let mut s = plate();
    let start = mm(&s.model);
    let n = interior(&mm(&s.model));
    let p0 = start.node(n).unwrap().xyz;
    let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0.2, 0.0, 0.0], "radius": 5.0, "propagate": "harmonic"})).unwrap();
    assert!(r["moved"].as_u64().unwrap() > 5, "{r}");
    let now = mm(&s.model);
    for (k, x) in now.nodes() {
        let old = start.node(*k).unwrap().xyz;
        if old != x.xyz {
            let d = dist(old, p0);
            assert!(d < 5.0 + 1e-9, "node {k} beyond the radius moved");
            assert!(dist(old, x.xyz) <= 0.2 + 1e-12);
        }
    }
    assert_eq!(invalid(&mut s), 0);
    let e = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0.1, 0.0, 0.0], "radius": 5.0, "propagate": "harmonic", "falloff": "linear"}));
    assert!(e.is_err());
    assert!(s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0.1, 0.0, 0.0], "propagate": "harmonic"})).is_err(), "needs a radius");
}

/// On the geometry: a node on the hole slides round it (stays on the cylinder), its followers on
/// the caps stay in their planes; a node at a corner of the plate is held (an error says how to
/// move it), and moves with the geometry off.
#[test]
fn on_geometry_keeps_nodes_on_faces_and_edges_and_holds_vertices() {
    let mut s = plate();
    let start = mm(&s.model);
    let on_hole = |p: [f64; 3]| (hole_radius(p) - HOLE.1).abs() < 1e-6;
    let n = nearest(&mm(&s.model), [20.0, 15.0, 5.0], |p| on_hole(p) && p[2] > 1e-6 && p[2] < 10.0 - 1e-6);
    let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0.0, 0.6, 0.4], "radius": 3.0, "propagate": "harmonic"})).unwrap();
    assert!(r["moved"].as_u64().unwrap() > 1, "{r}");
    assert_ne!(at(&s, n), start.node(n).unwrap().xyz);
    let now = mm(&s.model);
    for (k, x) in now.nodes() {
        let old = start.node(*k).unwrap().xyz;
        if on_hole(old) {
            assert!(on_hole(x.xyz), "node {k} left the hole: r = {}", hole_radius(x.xyz));
        }
        for z in [0.0, 10.0] {
            if (old[2] - z).abs() < 1e-9 {
                assert!((x.xyz[2] - z).abs() < 1e-9, "node {k} left its cap");
            }
        }
    }
    // A node on the hole's rim (the edge of the hole and a cap) slides along the rim.
    let rim = nearest(&start, [20.0, 15.0, 10.0], |p| on_hole(p) && (p[2] - 10.0).abs() < 1e-9);
    let mut s2 = plate();
    s2.execute("mesh.node.move", &json!({"nodes": [rim], "by": [0.0, 0.8, 0.5]})).unwrap();
    let x = at(&s2, rim);
    assert!(on_hole(x) && (x[2] - 10.0).abs() < 1e-9 && x != start.node(rim).unwrap().xyz, "{x:?}");
    // The plate's corner is a vertex.
    let corner = nearest(&start, [0.0, 0.0, 0.0], |_| true);
    let e = s2.execute("mesh.node.move", &json!({"nodes": [corner], "by": [0.3, 0.3, 0.3]})).unwrap_err();
    assert!(e.to_string().contains("on_geometry"), "{e}");
    s2.execute("mesh.node.move", &json!({"nodes": [corner], "by": [0.3, 0.3, 0.3], "on_geometry": false})).unwrap();
    assert!(dist(at(&s2, corner), [0.3, 0.3, 0.3]) < 1e-9);
}

/// The guard scales a move that would invert elements back so that none does; without it they
/// invert.
#[test]
fn the_guard_scales_back_to_no_inversion() {
    let mut s = plate();
    let n = interior(&mm(&s.model));
    let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [3.0, 0.0, 0.0]})).unwrap();
    let scale = r["scale"].as_f64().unwrap();
    assert!(scale > 0.0 && scale < 1.0, "{r}");
    assert_eq!(r["after"]["invalid"], 0, "{r}");
    assert_eq!(invalid(&mut s), 0);
    let mut s = plate();
    let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [3.0, 0.0, 0.0], "guard": "none"})).unwrap();
    assert_eq!(r["scale"], 1.0);
    assert!(r["after"]["invalid"].as_u64().unwrap() > 0, "{r}");
}

/// The preview of the interactive tool is the command's move, without doing it.
#[test]
fn the_preview_is_the_move_without_doing_it() {
    let mut s = plate();
    let n = interior(&mm(&s.model));
    let p = json!({"nodes": [n], "by": [0.2, 0.1, 0.0], "radius": 4.0});
    let prep = prepare(&s, &p).unwrap();
    assert!(prep.nodes() > 1 && prep.handles().contains(&n) && prep.sphere.is_some());
    let start = s.model.clone();
    let preview = prep.preview(&s, prep.displacement(&p).unwrap()).unwrap();
    assert_eq!(s.model, start);
    s.execute("mesh.node.move", &p).unwrap();
    for (k, x) in &preview.positions {
        assert_eq!(s.model.node(*k).unwrap().xyz, *x);
    }
    assert!(preview.elements.iter().any(|e| e.after != e.before));
}

#[test]
fn hostile_move_params_are_refused() {
    let mut s = plate();
    let n = interior(&mm(&s.model));
    let start = s.model.clone();
    for p in [
        json!({}),
        json!({"nodes": [n]}),
        json!({"by": [1, 0, 0]}),
        json!({"nodes": [n], "by": [1, 0]}),
        json!({"nodes": [n], "by": ["x", 0, 0]}),
        json!({"nodes": [n], "by": [1e308, 1e308, 1e308]}),
        json!({"nodes": [n], "by": [1, 0, 0], "to": [0, 0, 0]}),
        json!({"nodes": [n], "center": [0, 0, 0], "radius": 1, "by": [1, 0, 0]}),
        json!({"center": [0, 0, 0], "by": [1, 0, 0]}),
        json!({"nodes": [n], "by": [1, 0, 0], "radius": -1}),
        json!({"nodes": [n], "by": [1, 0, 0], "radius": "big"}),
        json!({"nodes": [n], "by": [1, 0, 0], "strength": 0}),
        json!({"nodes": [n], "by": [1, 0, 0], "strength": 2}),
        json!({"nodes": [n], "by": [1, 0, 0], "falloff": "cubic"}),
        json!({"nodes": [n], "by": [1, 0, 0], "propagate": "magic"}),
        json!({"nodes": [n], "by": [1, 0, 0], "guard": "maybe"}),
        json!({"nodes": [n], "by": [1, 0, 0], "on_geometry": "yes"}),
        json!({"nodes": [99_999_999], "by": [1, 0, 0]}),
        json!({"nodes": [-1], "by": [1, 0, 0]}),
        json!({"nodes": [], "by": [1, 0, 0]}),
        json!({"faces": [], "by": [1, 0, 0]}),
        json!({"faces": [{"body": "Nope", "face": "F1"}], "by": [1, 0, 0]}),
        json!({"center": [1e9, 0, 0], "radius": 1, "by": [1, 0, 0]}),
        json!({"nodes": [n], "by": [1, 0, 0], "radius": 1e-300}),
        json!({"center": [-1e300, -1e300, 0], "radius": 1e-300, "by": [1, 0, 0]}),
        json!({"nodes": [n], "by": [1, 0, 0], "radius": 1e-300, "propagate": "harmonic"}),
        json!({"nodes": [n], "by": [1, 0, 0], "min_quality": 1.5}),
        json!({"nodes": [n], "by": [1, 0, 0], "min_quality": -0.1}),
        json!({"nodes": [n], "by": [1, 0, 0], "guard": "none", "min_quality": 0.2}),
    ] {
        assert!(s.execute("mesh.node.move", &p).is_err(), "{p} should be refused");
        assert_eq!(s.model, start, "{p} changed the model");
    }
}

/// A zero move, or one the guard scales back to nothing, changes nothing and records no undo
/// step; `min_quality` sets how far the guard lets a move go.
#[test]
fn a_move_of_nothing_records_no_undo_step_and_min_quality_is_adjustable() {
    let mut s = plate();
    let n = interior(&mm(&s.model));
    let (journal, start) = (s.journal.len(), s.model.clone());
    let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [0, 0, 0], "radius": 3.0})).unwrap();
    assert_eq!((r["moved"].as_u64(), r["nothing_moved"].as_bool()), (Some(0), Some(true)), "{r}");
    assert_eq!((s.journal.len(), &s.model), (journal, &start));
    let scale = |s: &mut Session, q: f64| {
        let r = s.execute("mesh.node.move", &json!({"nodes": [n], "by": [3.0, 0.0, 0.0], "min_quality": q})).unwrap();
        s.execute("edit.undo", &json!({})).unwrap();
        (r["scale"].as_f64().unwrap(), r["after"]["min_scaled_jacobian"].as_f64().unwrap())
    };
    let (loose, strict) = (scale(&mut s, 0.0), scale(&mut s, 0.5));
    assert!(strict.0 < loose.0 && strict.1 >= 0.5 - 1e-9, "{loose:?} {strict:?}");
}

/// The live preview's cap: more followers than it allows are counted, not prepared; the live
/// preview is the command's move.
#[test]
fn a_capped_prepare_counts_a_large_region_and_the_live_preview_matches() {
    let s = plate();
    let n = interior(&mm(&s.model));
    let p = json!({"nodes": [n], "by": [0.2, 0.1, 0.0], "radius": 4.0, "propagate": "harmonic"});
    assert!(matches!(prepare_capped(&s, &p, Some(3)).unwrap(), Capped::TooLarge(k) if k > 3));
    let Capped::Ready(prep) = prepare_capped(&s, &p, None).unwrap() else { panic!("not prepared") };
    let (exact, first, again) = (
        prep.preview(&s, [0.2, 0.1, 0.0]).unwrap(),
        prep.preview_live(&s, [0.2, 0.1, 0.0]).unwrap(),
        prep.preview_live(&s, [0.2, 0.1, 0.0]).unwrap(),
    );
    assert_eq!(exact.positions, first.positions);
    assert_eq!(first, again, "the reused baseline gives the same move");
}
