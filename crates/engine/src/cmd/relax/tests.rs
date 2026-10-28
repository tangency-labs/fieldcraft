use super::*;

/// The plate (40 × 30 × 10) with its Ø10 hole round (15, 15): a point on the hole's wall.
const ON_HOLE: [f64; 3] = [20.0, 15.0, 5.0];
const HOLE: ([f64; 2], f64) = ([15.0, 15.0], 5.0);

fn plate() -> Session {
    super::super::datum::tests::plate().0
}

/// The model in millimetres, as the plate is drawn (it is stored in metres).
fn mm(m: &Model) -> Model {
    m.converted(&|q, v| crate::units::geometry().express(q, v).ok()).unwrap()
}

fn hole_radius(p: [f64; 3]) -> f64 {
    ((p[0] - HOLE.0[0]).powi(2) + (p[1] - HOLE.0[1]).powi(2)).sqrt()
}

/// Which analytic face a node of the plate lies on, if any.
fn on_faces(p: [f64; 3]) -> Vec<&'static str> {
    let mut f = Vec::new();
    if (hole_radius(p) - HOLE.1).abs() < 1e-6 {
        f.push("hole");
    }
    if p[2].abs() < 1e-6 || (p[2] - 10.0).abs() < 1e-6 {
        f.push("cap");
    }
    f
}

/// Spikes near the hole: interior nodes pushed about; nodes on the hole turned round it and
/// moved along its axis (still on the cylinder), nodes on the caps moved in their plane.
fn perturb(s: &mut Session, centre: [f64; 3], radius: f64, amount: f64) {
    let nodes: Vec<(u32, [f64; 3])> = mm(&s.model).nodes().iter().map(|(n, p)| (*n, p.xyz)).collect();
    for (k, (n, p)) in nodes.into_iter().enumerate() {
        let d2: f64 = (0..3).map(|i| (p[i] - centre[i]).powi(2)).sum();
        if d2 > radius * radius {
            continue;
        }
        let wiggle = |m: usize| amount * ((k * m % 7) as f64 - 3.0) / 3.0;
        let faces = on_faces(p);
        let q = match faces.as_slice() {
            [] => [p[0] + wiggle(3), p[1] + wiggle(5), p[2] + wiggle(11)],
            ["hole"] => {
                let a = wiggle(3) / HOLE.1;
                let (x, y) = (p[0] - HOLE.0[0], p[1] - HOLE.0[1]);
                [HOLE.0[0] + x * a.cos() - y * a.sin(), HOLE.0[1] + x * a.sin() + y * a.cos(), p[2] + wiggle(5)]
            }
            ["cap"] if hole_radius(p) > HOLE.1 + 1e-6 => [p[0] + wiggle(3), p[1] + wiggle(5), p[2]],
            _ => continue,
        };
        // Through the command (millimetres, as typed), as develop moved them; a mid-side node,
        // which the command does not move, directly (its millimetres stored as metres).
        if s.execute("mesh.node.move", &json!({"nodes": [n], "to": q, "on_geometry": false, "guard": "none"})).is_err() {
            s.model.set_node(n, crate::units::geometry_point(q)).unwrap();
        }
    }
}

fn check_on_geometry(before: &Model, after: &Model) {
    for (n, p) in after.nodes() {
        let old = before.node(*n).unwrap().xyz;
        for f in on_faces(old) {
            match f {
                "hole" => assert!((hole_radius(p.xyz) - HOLE.1).abs() < 1e-6, "node {n} left the hole: r = {}", hole_radius(p.xyz)),
                _ => assert!((p.xyz[2] - old[2]).abs() < 1e-6, "node {n} left its cap"),
            }
        }
    }
}

fn invalid(s: &mut Session) -> u64 {
    s.execute("mesh.quality", &json!({})).unwrap()["invalid"]["count"].as_u64().unwrap_or(0)
}

/// The issue's hex test on a design: spikes near the hole relax, every surface node stays on its
/// face (the hole's cylinder exactly), quality rises, nothing inverts, nothing outside the sphere
/// moves; one undo step; the field is kept and saved with the model.
#[test]
fn hex_spikes_near_a_hole_relax_on_the_geometry() {
    let mut s = plate();
    s.execute("mesh.hex.generate", &json!({"size": 2.0})).unwrap();
    perturb(&mut s, ON_HOLE, 3.0, 0.35);
    let perturbed = s.model.clone();
    let r = s.execute("mesh.node.relax", &json!({"center": ON_HOLE, "radius": 4.0, "iterations": 10})).unwrap();
    let (before, after) = (r["before"]["min_scaled_jacobian"].as_f64().unwrap(), r["after"]["min_scaled_jacobian"].as_f64().unwrap());
    assert!(after > before + 0.05, "{r}");
    assert_eq!(r["after"]["invalid"], 0, "{r}");
    assert!(r["moved"].as_u64().unwrap() > 0 && r["field"] == "Relax1", "{r}");
    check_on_geometry(&mm(&perturbed), &mm(&s.model));
    let (perturbed_mm, now) = (mm(&perturbed), mm(&s.model));
    for (n, p) in now.nodes() {
        let old = perturbed_mm.node(*n).unwrap().xyz;
        if old != p.xyz {
            let d2: f64 = (0..3).map(|i| (old[i] - ON_HOLE[i]).powi(2)).sum();
            assert!(d2 < 16.0, "node {n} outside the sphere moved");
        }
    }
    assert_eq!(invalid(&mut s), 0);
    // The field is saved with the model.
    let text = serde_json::to_string(&s.model).unwrap();
    let back: Model = serde_json::from_str(&text).unwrap();
    assert_eq!(back.relaxation("Relax1"), s.model.relaxation("Relax1"));
    s.execute("edit.undo", &json!({})).unwrap();
    assert_eq!(s.model, perturbed);
}

/// The same on tet10 with curved mid-side nodes on the hole: they stay on the cylinder.
#[test]
fn tet10_relaxes_with_curved_mid_side_nodes_on_the_hole() {
    let mut s = plate();
    s.execute("mesh.tet.generate", &json!({"size": 2.5, "order": "quadratic"})).unwrap();
    perturb(&mut s, ON_HOLE, 3.0, 0.2);
    // Put the mid-side nodes back on their edges' middles, on the hole where the edge is.
    let mids: Vec<(u32, u32, u32)> = s
        .model
        .elements()
        .values()
        .flat_map(|e| e.kind.edges().iter().map(|ed| (e.nodes[ed[0]], e.nodes[ed[1]], e.nodes[ed[2]])).collect::<Vec<_>>())
        .collect();
    let original = mm(&s.model);
    for (a, b, m) in mids {
        let (pa, pb, pm) = (original.node(a).unwrap().xyz, original.node(b).unwrap().xyz, original.node(m).unwrap().xyz);
        let mut q = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0];
        if on_faces(pm).contains(&"hole") {
            let r = hole_radius(q);
            q = [HOLE.0[0] + (q[0] - HOLE.0[0]) * HOLE.1 / r, HOLE.0[1] + (q[1] - HOLE.0[1]) * HOLE.1 / r, q[2]];
        }
        s.model.set_node(m, crate::units::geometry_point(q)).unwrap();
    }
    let perturbed = s.model.clone();
    let r = s.execute("mesh.node.relax", &json!({"center": ON_HOLE, "radius": 4.0, "iterations": 8})).unwrap();
    assert!(r["after"]["min_scaled_jacobian"].as_f64().unwrap() >= r["before"]["min_scaled_jacobian"].as_f64().unwrap(), "{r}");
    assert_eq!(r["after"]["invalid"], 0, "{r}");
    check_on_geometry(&mm(&perturbed), &mm(&s.model));
}

/// A partition cut is not geometry: nodes on it move freely, shared by the cells on both sides.
#[test]
fn nodes_on_a_partition_cut_move_freely() {
    let mut s = plate();
    s.execute("geometry.datum", &json!({"name": "Cut", "by": {"type": "offset", "plane": "yz", "offset": 27}})).unwrap();
    s.execute("geometry.partition", &json!({"name": "P", "datum": "Cut"})).unwrap();
    s.execute("mesh.cells", &json!({"size": 2.0, "default": "hex"})).unwrap();
    let centre = [27.0, 15.0, 5.0];
    perturb(&mut s, centre, 2.5, 0.3);
    let before = s.model.clone();
    s.execute("mesh.node.relax", &json!({"center": centre, "radius": 3.5, "iterations": 10})).unwrap();
    let (before, now) = (mm(&before), mm(&s.model));
    let moved_off = now
        .nodes()
        .iter()
        .filter(|(n, p)| {
            let old = before.node(**n).unwrap().xyz;
            (old[0] - 27.0).abs() < 1e-6 && on_faces(old).is_empty() && (p.xyz[0] - 27.0).abs() > 1e-6
        })
        .count();
    assert!(moved_off > 0, "nodes on the cut are free");
    assert_eq!(invalid(&mut s), 0);
}

/// A field marked `reapply` is applied again after a remesh; `mesh.node.relax.apply` applies it
/// on demand; around the worst element works.
#[test]
fn a_kept_field_is_applied_again_after_a_remesh() {
    let mut s = plate();
    s.execute("mesh.hex.generate", &json!({"size": 3.0})).unwrap();
    perturb(&mut s, ON_HOLE, 4.0, 0.3);
    let r = s.execute("mesh.node.relax", &json!({"center": ON_HOLE, "radius": 5.0, "field": "Spike", "reapply": true})).unwrap();
    assert_eq!(r["field"], "Spike");
    let r = s.execute("mesh.hex.generate", &json!({"size": 1.5})).unwrap();
    assert_eq!(r["relaxations_applied"][0]["field"], "Spike", "{r}");
    assert_eq!(invalid(&mut s), 0);
    let r = s.execute("mesh.node.relax.apply", &json!({"field": "Spike"})).unwrap();
    assert!(r["summary"].as_str().unwrap().contains("Spike"));
    let worst = s.execute("mesh.quality", &json!({})).unwrap()["metrics"]["scaled_jacobian"]["worst"].as_u64().unwrap();
    let r = s.execute("mesh.node.relax", &json!({"around_element": worst})).unwrap();
    assert!(r["after"]["min_scaled_jacobian"].as_f64() >= r["before"]["min_scaled_jacobian"].as_f64(), "{r}");
}

#[test]
fn hostile_parameters_never_panic() {
    let mut s = plate();
    s.execute("mesh.hex.generate", &json!({"size": 4.0})).unwrap();
    for p in [
        json!({}),
        json!({"center": [0, 0, 0]}),
        json!({"center": [1e300, 0, 0], "radius": 1}),
        json!({"center": [0, 0, 0], "radius": -1}),
        json!({"center": [0, 0, 0], "radius": 1, "nodes": [1]}),
        json!({"nodes": [1], "strength": 2}),
        json!({"nodes": [1], "iterations": 0}),
        json!({"nodes": [1], "iterations": 1000}),
        json!({"nodes": [1], "falloff": "cubic"}),
        json!({"around_element": 99999999}),
        json!({"set": "nope"}),
        json!({"nodes": [1], "bodies": ["nope"]}),
    ] {
        assert!(s.execute("mesh.node.relax", &p).is_err(), "{p}");
    }
    assert!(s.execute("mesh.node.relax.apply", &json!({"field": "nope"})).is_err());
    assert!(s.execute("mesh.node.relax.apply", &json!({})).is_err());
}

/// A 20 mm box with an R2 fillet along its edge at x = z = 20 (axis x = z = 18, along y), tet
/// meshed; its fillet faces' persistent names.
fn filleted_box() -> (crate::tests::TempDir, Session, Vec<String>) {
    let script = json!({"commands": [
        {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 20, "body_name": "B"}},
        {"command": "solid.fillet", "params": {"edges": [[20, 10, 20]], "radius": 2}},
    ]});
    let bytes = fieldcraft_geometry::script_design(&script.to_string()).unwrap();
    let (dir, path) = crate::tests::design_file(&bytes);
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": path})).unwrap();
    let features = s.execute("geometry.features", &json!({"body": "B"})).unwrap();
    let fillet = features["bodies"][0]["features"].as_array().unwrap().iter().find(|f| f["type"] == "fillet").unwrap().clone();
    let faces: Vec<String> = fillet["faces"].as_array().unwrap().iter().map(|f| f.as_str().unwrap().to_string()).collect();
    assert!(!faces.is_empty(), "{fillet}");
    s.execute("mesh.tet.generate", &json!({"size": 1.5})).unwrap();
    (dir, s, faces)
}

fn fillet_radius(p: [f64; 3]) -> f64 {
    ((p[0] - 18.0).powi(2) + (p[2] - 18.0).powi(2)).sqrt()
}

/// #471: relaxing by CAD face. Relaxing the fillet faces of the tet mesh raises the worst
/// element there, inverts none, keeps the fillet's nodes on its cylinder
/// and moves nothing away from the fillet's elements. Faces that are not meshed are refused.
#[test]
fn relaxing_fillet_faces_improves_them() {
    let (_d, mut s, faces) = filleted_box();
    let perturbed = s.model.clone();
    let refs: Vec<Value> = faces.iter().map(|f| json!({"body": "B", "face": f})).collect();
    let r = s.execute("mesh.node.relax", &json!({"faces": refs, "iterations": 10})).unwrap();
    let (before, after) = (r["before"]["min_scaled_jacobian"].as_f64().unwrap(), r["after"]["min_scaled_jacobian"].as_f64().unwrap());
    assert!(after > before, "{r}");
    assert_eq!(r["after"]["invalid"], 0, "{r}");
    assert_eq!(invalid(&mut s), 0);
    let (perturbed, now) = (mm(&perturbed), mm(&s.model));
    for (n, p) in now.nodes() {
        let old = perturbed.node(*n).unwrap().xyz;
        if old[0] > 18.0 && old[2] > 18.0 && (fillet_radius(old) - 2.0).abs() < 1e-6 {
            assert!((fillet_radius(p.xyz) - 2.0).abs() < 1e-6, "node {n} left the fillet");
        }
        if old != p.xyz {
            assert!(fillet_radius(old) < 6.0, "node {n} far from the fillet moved");
        }
    }
    // Faces: required objects, known and meshed; not with a sphere.
    for bad in [
        json!({"faces": []}),
        json!({"faces": "F1"}),
        json!({"faces": [{"body": "B", "face": "Nope"}]}),
        json!({"faces": [{"body": "B", "face": "Nope"}, {"body": "B", "face": "Nope2"}]}),
        json!({"faces": [{"body": "Nope", "face": faces[0]}]}),
        json!({"faces": [{"body": "B"}]}),
        json!({"faces": refs, "center": [0, 0, 0], "radius": 1}),
    ] {
        assert!(s.execute("mesh.node.relax", &bad).is_err(), "{bad}");
    }
    // A listed face with no mesh is named, the others relaxed.
    let mut some = refs.clone();
    some.push(json!({"body": "B", "face": "Nope"}));
    let r = s.execute("mesh.node.relax", &json!({"faces": some})).unwrap();
    assert_eq!(r["faces_without_mesh"], json!([{"body": "B", "face": "Nope"}]), "{r}");
    // Elements by label are accepted too.
    let e = *s.model.elements().keys().next().unwrap();
    s.execute("mesh.node.relax", &json!({"elements": [e]})).unwrap();
    assert!(s.execute("mesh.node.relax", &json!({"set": "All"})).is_err());
}

/// Whether a node of the 20 mm box lies on two of its planes (an edge or a vertex of the box,
/// or a tangent line of its fillet).
fn on_box_edge(p: [f64; 3]) -> bool {
    p.iter().filter(|c| c.abs() < 1e-6 || (**c - 20.0).abs() < 1e-6).count() >= 2
}

/// Strong mode on the whole filleted box: the smallest scaled Jacobian does not drop, the mean rises,
/// nothing inverts, the fillet's nodes stay on its cylinder, the plane faces' nodes on their
/// planes, and the nodes on the box's edges and the fillet's tangent lines do not move. The
/// reply carries each round's region and quality.
#[test]
fn strong_relaxation_of_a_filleted_box_keeps_faces_and_edges() {
    let (_d, mut s, _) = filleted_box();
    let before = s.model.clone();
    // The tet mesher leaves no element below the skill's 0.3, so relax those below 0.45.
    let r = s.execute("mesh.node.relax", &json!({"mode": "strong", "bodies": ["B"], "below": 0.45})).unwrap();
    let (b, a) = (&r["before"], &r["after"]);
    // The worst tets may not improve: their corners can all be held on the box's edges.
    assert!(a["min_scaled_jacobian"].as_f64().unwrap() >= b["min_scaled_jacobian"].as_f64().unwrap(), "{r}");
    assert!(a["mean_scaled_jacobian"].as_f64().unwrap() > b["mean_scaled_jacobian"].as_f64().unwrap(), "{r}");
    assert_eq!(a["invalid"], 0, "{r}");
    assert_eq!(invalid(&mut s), 0);
    let rounds = r["rounds"].as_array().unwrap();
    assert!(!rounds.is_empty() && rounds.iter().all(|k| k["region_nodes"].as_u64() >= k["bad_nodes"].as_u64() && k["min"].is_number()), "{r}");
    let (before, now) = (mm(&before), mm(&s.model));
    for (n, p) in now.nodes() {
        let (old, new) = (before.node(*n).unwrap().xyz, p.xyz);
        if on_box_edge(old) {
            assert_eq!(old, new, "node {n} on an edge moved");
        }
        if old[0] > 18.0 && old[2] > 18.0 && (fillet_radius(old) - 2.0).abs() < 1e-6 {
            assert!((fillet_radius(new) - 2.0).abs() < 1e-6, "node {n} left the fillet");
        }
        for k in 0..3 {
            for plane in [0.0, 20.0] {
                if (old[k] - plane).abs() < 1e-9 {
                    assert!((new[k] - plane).abs() < 1e-6, "node {n} left the plane {k} = {plane}");
                }
            }
        }
    }
}

/// Strong mode on a distorted hex plate: over the guarded mode's sphere it lifts the smallest
/// scaled Jacobian past the threshold, close to the guarded mode (which also slides the edge nodes) (with nodes on edges held, where the guarded mode slides
/// them); by quality threshold (the elements below 0.5, morphed 4 mm round) it untangles
/// every element. Nothing inverts, nodes on the hole stay on its cylinder and the caps' on their
/// planes, nodes on the hole's edges do not move.
#[test]
fn strong_relaxation_untangles_a_distorted_plate() {
    let mut s = plate();
    s.execute("mesh.hex.generate", &json!({"size": 2.0})).unwrap();
    perturb(&mut s, ON_HOLE, 6.0, 0.35);
    let perturbed = s.model.clone();
    let sphere = json!({"center": ON_HOLE, "radius": 7.0, "iterations": 20});
    let guarded = s.execute("mesh.node.relax", &sphere).unwrap();
    let min = |r: &Value, at: &str| r[at]["min_scaled_jacobian"].as_f64().unwrap();
    for strong in [json!({"mode": "strong", "center": ON_HOLE, "radius": 7.0}), json!({"mode": "strong", "below": 0.5, "reach": 4.0})] {
        s.model = perturbed.clone();
        let r = s.execute("mesh.node.relax", &strong).unwrap();
        assert!(min(&r, "after") > strong["below"].as_f64().unwrap_or(0.7) && min(&r, "before") < 0.0, "{strong}: {r}");
        // Close to the guarded mode, which also slides the hole's edge nodes.
        assert!(strong.get("below").is_some() || min(&r, "after") > min(&guarded, "after") - 0.05, "strong {r}\nguarded {guarded}");
        assert_eq!(r["after"]["invalid"], 0, "{r}");
        assert_eq!(invalid(&mut s), 0);
        check_on_geometry(&mm(&perturbed), &mm(&s.model));
        let (perturbed_mm, now) = (mm(&perturbed), mm(&s.model));
        for (n, p) in now.nodes() {
            let old = perturbed_mm.node(*n).unwrap().xyz;
            if on_faces(old).len() >= 2 {
                assert_eq!(old, p.xyz, "node {n} on the hole's edge moved");
            }
        }
    }
}

#[test]
fn strong_mode_parameters_are_checked() {
    let mut s = plate();
    s.execute("mesh.hex.generate", &json!({"size": 4.0})).unwrap();
    for p in [
        json!({"mode": "wild", "nodes": [1]}),
        json!({"mode": "strong", "below": 0}),
        json!({"mode": "strong", "below": 1.5}),
        json!({"mode": "strong", "below": 0.5, "reach": 0}),
        json!({"mode": "strong", "below": 0.5, "reach": -1}),
        json!({"mode": "strong", "bodies": ["nope"]}),
        json!({"mode": "strong", "nodes": [1], "iterations": 0}),
        json!({"mode": "strong", "nodes": [1], "iterations": 201}),
        json!({"below": 0.5}),
        json!({"nodes": [1], "reach": 1}),
        json!({"nodes": [1], "guards": "strict"}),
        json!({"mode": "strong", "nodes": [1], "guards": "loose"}),
        json!({"nodes": [1], "mode": "strong", "layers": 1}),
    ] {
        assert!(s.execute("mesh.node.relax", &p).is_err(), "{p}");
    }
    // Nothing below the threshold: nothing moves.
    let r = s.execute("mesh.node.relax", &json!({"mode": "strong", "below": 0.01})).unwrap();
    assert_eq!(r["moved"], 0, "{r}");
}

/// The mesh and the strong relaxation are the same bits on every platform and thread count:
/// the plate's hex mesh differed between macOS and Linux (the system maths library's `atan2`,
/// `sin` … are not correctly rounded), so the distorted-plate test passed on one and failed on
/// the other (#438). The fingerprint is of the nodes after mesh.hex.generate, after the
/// perturbation and after the strong relaxation; CI checks it on Linux, it was taken on macOS
/// (again for the model stored in SI, #579, and when the plate's O-grid kept the cut that
/// rounding used to refuse, #605).
#[test]
fn mesh_and_strong_relaxation_are_the_same_bits_everywhere() {
    let fingerprint = |s: &Session, h: &mut u64| {
        // FNV-1a over every node's label and coordinate bits.
        for (n, p) in s.model.nodes() {
            for word in [u64::from(*n), p.xyz[0].to_bits(), p.xyz[1].to_bits(), p.xyz[2].to_bits()] {
                for b in word.to_le_bytes() {
                    *h = (*h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
                }
            }
        }
    };
    let run = || {
        let mut h = 0xcbf2_9ce4_8422_2325;
        let mut s = plate();
        s.execute("mesh.hex.generate", &json!({"size": 2.0})).unwrap();
        fingerprint(&s, &mut h);
        perturb(&mut s, ON_HOLE, 6.0, 0.35);
        fingerprint(&s, &mut h);
        s.execute("mesh.node.relax", &json!({"mode": "strong", "below": 0.5, "reach": 4.0})).unwrap();
        fingerprint(&s, &mut h);
        h
    };
    let one = fieldcraft_parallel::with_threads(1, run);
    let many = fieldcraft_parallel::with_threads(4, run);
    assert_eq!(one, many, "1 and 4 threads differ");
    // The meshers alone give what develop gave, stored in metres. A relaxation (guarded or
    // strong) starts from the millimetre image of the stored metres, which differs from the
    // millimetres develop kept by an ulp for 6–18% of computed coordinates, and amplifies that
    // (on develop a 1-ulp nudge moves it as much): its bits differ from develop's. The plate's
    // tessellation comes from SolidCraft, so a SolidCraft bump moves the fingerprint too (v0.25.1).
    assert_eq!(one, 0x042c_c785_a236_378a, "fingerprint {one:#018x}");
}

/// A mesher's result depends on the stored model alone: relaxing, undoing and relaxing again
/// gives the same mesh, and so does relaxing the model saved and opened again (what a replay
/// and a crash replay see).
#[test]
fn relaxing_depends_on_the_stored_model_alone() {
    let mut s = plate();
    s.execute("mesh.hex.generate", &json!({"size": 2.0})).unwrap();
    perturb(&mut s, ON_HOLE, 6.0, 0.35);
    // A command in another frame between (the display length), and a display override.
    s.execute("model.units", &json!({"system": "mm-t-s-N-MPa", "overrides": {"length": "m"}})).unwrap();
    let dir = std::env::temp_dir().join(format!("fieldcraft-relax-history-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    s.execute("mesh.export.stl", &json!({"path": dir.join("skin.stl")})).unwrap();
    let relax = json!({"mode": "strong", "below": 0.5, "reach": 4.0});
    s.execute("mesh.node.relax", &relax).unwrap();
    let live = s.model.clone();
    s.execute("edit.undo", &json!({})).unwrap();
    s.execute("mesh.node.relax", &relax).unwrap();
    assert!(s.model == live, "relax, undo, relax again");
    s.execute("edit.undo", &json!({})).unwrap();
    let path = dir.join("plate.fieldcraft");
    s.execute("file.save_as", &json!({"path": path})).unwrap();
    let mut opened = Session::default();
    opened.execute("file.open", &json!({"path": path})).unwrap();
    opened.execute("mesh.node.relax", &relax).unwrap();
    assert!(opened.model.nodes() == live.nodes(), "saved, opened and relaxed");
    let _ = std::fs::remove_dir_all(&dir);
}
