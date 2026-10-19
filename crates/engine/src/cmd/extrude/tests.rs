use super::*;
use crate::tests::{TempDir, design_file};

fn open(commands: Value) -> (TempDir, Session) {
    let bytes = fieldcraft_geometry::script_design(&json!({ "commands": commands }).to_string()).unwrap();
    let (dir, path) = design_file(&bytes);
    let mut s = Session::default();
    s.execute("geometry.open", &json!({"path": path})).unwrap();
    (dir, s)
}

/// A square post whose top is cut by a cylinder (axis X, radius 20, at z = 30), a square post
/// whose top is a sphere cap (centre z = −5, radius 25), and a block under the first post.
fn parts() -> (TempDir, Session) {
    open(json!([
        {"command": "solid.box", "params": {"corner": [-4, -4, 0], "length": 8, "width": 8, "height": 30, "body_name": "Post"}},
        {"command": "solid.cylinder", "params": {"base": [-20, 0, 30], "axis": [1, 0, 0], "radius": 20, "height": 40, "body_name": "Bar", "operation": "new"}},
        {"command": "solid.combine", "params": {"target": "Post", "tools": ["Bar"], "operation": "cut"}},
        {"command": "solid.box", "params": {"corner": [-6, -6, 0], "length": 12, "width": 12, "height": 30, "body_name": "Dome", "operation": "new"}},
        {"command": "solid.sphere", "params": {"center": [0, 0, -5], "radius": 25, "body_name": "Ball", "operation": "new"}},
        {"command": "solid.combine", "params": {"target": "Dome", "tools": ["Ball"], "operation": "intersect"}},
        {"command": "solid.box", "params": {"corner": [-4, -4, -10], "length": 8, "width": 8, "height": 10, "body_name": "Base", "operation": "new"}},
    ]))
}

/// `{"body", "face"}` of each face of `body` whose surface `pick` accepts.
fn faces(s: &Session, body: &str, pick: impl Fn(&FaceSurface) -> bool) -> Vec<Value> {
    let b = s.design.as_deref().unwrap().body(body).unwrap();
    b.faces()
        .iter()
        .zip(b.face_surfaces().unwrap())
        .filter(|(_, f)| f.as_ref().is_some_and(&pick))
        .map(|(f, _)| json!({"body": body, "face": f.name}))
        .collect()
}

fn plane(n: [f64; 3]) -> impl Fn(&FaceSurface) -> bool {
    move |f| matches!(f, FaceSurface::Plane { normal, .. } if (0..3).all(|k| (normal[k] - n[k]).abs() < 1e-9))
}

fn sides(s: &Session, body: &str) -> Vec<Value> {
    faces(s, body, |f| matches!(f, FaceSurface::Plane { normal, .. } if normal[2].abs() < 1e-9))
}

/// The coordinates of the nodes on face `face` (through a set by geometry).
fn nodes_on(s: &mut Session, name: &str, face: &Value) -> Vec<[f64; 3]> {
    s.execute("model.set.create", &json!({"name": name, "geometry": [face]})).unwrap();
    let r = s.model.resolve_set(name).unwrap();
    // In millimetres, as the design is.
    r.nodes.iter().map(|n| crate::units::to_geometry_point(s.model.nodes()[n].xyz)).collect()
}

#[test]
fn a_square_post_extrudes_up_to_a_cylinder() {
    let (_d, mut s) = parts();
    let bottom = faces(&s, "Post", plane([0.0, 0.0, -1.0]));
    let top = faces(&s, "Post", |f| matches!(f, FaceSurface::Cylinder { .. }));
    let r = s
        .execute("mesh.hex.extrude", &json!({"source": bottom, "to_surface": top, "sides": sides(&s, "Post"), "direction": [0, 0, 1], "size": 1.0}))
        .unwrap();
    assert_eq!((r["body"].as_str(), r["section"]["from"].as_str()), (Some("Post"), Some("meshed")), "{r}");
    assert_eq!(r["layers"], 10, "{r}");
    let min_sj = r["quality"]["min_scaled_jacobian"].as_f64().unwrap();
    assert!(min_sj >= 0.5, "{r}");
    assert!(r["max_target_gap"].as_f64().unwrap() < 1e-6, "{r}");
    // ∫∫ (30 − √(400 − y²)) dx dy over the 8 × 8 square.
    let exact = 8.0 * (30.0 * 8.0 - (4.0 * 384.0f64.sqrt() + 400.0 * 0.2f64.asin()));
    let volume = r["volume"].as_f64().unwrap();
    assert!((volume - exact).abs() < 1e-3 * exact, "{volume} vs {exact}");
    // Every node of the last layer is on the cylinder: the target faces are tagged.
    let on = nodes_on(&mut s, "Top", &top[0]);
    assert_eq!(on.len() as u64, r["section"]["nodes"].as_u64().unwrap(), "one node per section node");
    for p in on {
        let radius = (p[1].powi(2) + (p[2] - 30.0).powi(2)).sqrt();
        assert!((radius - 20.0).abs() < 1e-6, "{p:?}: {radius}");
    }
    let q = s.execute("mesh.quality", &json!({})).unwrap();
    assert_eq!(q["invalid"]["count"], 0, "{q}");
}

#[test]
fn a_square_post_extrudes_up_to_a_sphere_cap_along_its_normal() {
    let (_d, mut s) = parts();
    let bottom = faces(&s, "Dome", plane([0.0, 0.0, -1.0]));
    let cap = faces(&s, "Dome", |f| matches!(f, FaceSurface::Sphere { .. }));
    let r = s
        .execute("mesh.hex.extrude", &json!({"source": bottom, "to_surface": cap, "sides": sides(&s, "Dome"), "layers": 8, "bias": 0.5, "size": 1.5}))
        .unwrap();
    assert_eq!(r["layers"], 8, "{r}");
    assert!(r["quality"]["min_scaled_jacobian"].as_f64().unwrap() >= 0.5, "{r}");
    let on = nodes_on(&mut s, "Cap", &cap[0]);
    assert!(on.len() > 50, "{}", on.len());
    for p in on {
        let radius = (p[0].powi(2) + p[1].powi(2) + (p[2] + 5.0).powi(2)).sqrt();
        assert!((radius - 25.0).abs() < 1e-6, "{p:?}: {radius}");
    }
}

#[test]
fn the_extrusion_conforms_with_the_neighbour_it_starts_from() {
    let (_d, mut s) = parts();
    s.execute("mesh.hex.generate", &json!({"body": "Base", "size": 2.0})).unwrap();
    let base_nodes = s.model.nodes().len();
    let bottom = faces(&s, "Post", plane([0.0, 0.0, -1.0]));
    let top = faces(&s, "Post", |f| matches!(f, FaceSurface::Cylinder { .. }));
    let r = s
        .execute("mesh.hex.extrude", &json!({"source": bottom, "to_surface": top, "sides": sides(&s, "Post"), "direction": [0, 0, 1], "size": 2.0}))
        .unwrap();
    assert_eq!(r["section"]["from"], "neighbour", "{r}");
    assert_eq!(r["section"]["shared_nodes"], 25, "the base's 5 × 5 top nodes: {r}");
    assert_eq!(s.model.nodes().len(), base_nodes + r["new_nodes"].as_u64().unwrap() as usize);
    let d = s.execute("mesh.check.duplicates", &json!({})).unwrap();
    assert_eq!((d["counts"]["coincident"].as_u64(), d["counts"]["cracks"].as_u64()), (Some(0), Some(0)), "{d}");
    // The base keeps its mesh; extruding again replaces only the post's.
    let base = s.model.body_elements("Base").len();
    s.execute("mesh.hex.extrude", &json!({"source": bottom, "to_surface": top, "direction": [0, 0, 1], "layers": 3})).unwrap();
    assert_eq!(s.model.body_elements("Base").len(), base);
    assert_eq!(s.model.body_elements("Post").len(), 16 * 3);
}

#[test]
fn a_2d_mesh_extrudes_from_its_set() {
    let (_d, mut s) = parts();
    s.execute("mesh.planar", &json!({"polygon": [[-3, -3], [3, -3], [3, 3], [-3, 3]], "size": 1.0, "name": "Plane"})).unwrap();
    s.execute("model.set.create", &json!({"name": "Section", "geometry": [{"body": "Plane", "face": "area"}]})).unwrap();
    let top = faces(&s, "Post", |f| matches!(f, FaceSurface::Cylinder { .. }));
    let r = s.execute("mesh.hex.extrude", &json!({"source_set": "Section", "to_surface": top, "body": "Column", "layers": 5})).unwrap();
    assert_eq!(r["section"]["from"], "set", "{r}");
    let quads = r["section"]["quadrilaterals"].as_u64().unwrap() as usize;
    assert_eq!(r["section"]["shared_nodes"], r["section"]["nodes"], "{r}");
    assert_eq!(s.model.body_elements("Column").len(), quads * 5);
    assert!(r["quality"]["min_scaled_jacobian"].as_f64().unwrap() >= 0.5, "{r}");
}

#[test]
fn hostile_params_are_errors_and_change_nothing() {
    let (_d, mut s) = parts();
    let bottom = faces(&s, "Post", plane([0.0, 0.0, -1.0]));
    let top = faces(&s, "Post", |f| matches!(f, FaceSurface::Cylinder { .. }));
    let ok = json!({"source": bottom, "to_surface": top, "direction": [0, 0, 1], "layers": 2});
    let with = |k: &str, v: Value| {
        let mut p = ok.clone();
        p[k] = v;
        p
    };
    let without = |k: &str| {
        let mut p = ok.clone();
        p.as_object_mut().unwrap().remove(k);
        p
    };
    let cases = [
        without("to_surface"),
        without("source"),
        with("to_surface", json!([])),
        with("to_surface", json!("F1")),
        with("to_surface", json!([{"body": "Nope", "face": "F1"}])),
        with("to_surface", json!([{"body": "Post", "face": "nope"}])),
        with("to_surface", json!([{"body": "Post"}])),
        with("to_surface", bottom.clone().into()),
        with("source_set", json!("Nope")),
        with("direction", json!([0, 0, 0])),
        with("direction", json!([0, 0, -1])),
        with("direction", json!([1, 0])),
        with("direction", json!([f64::MAX, 0, 0])),
        with("normal", json!(false)),
        with("layers", json!(0)),
        with("layers", json!(1e9)),
        with("layers", json!(-3)),
        with("bias", json!(0)),
        with("bias", json!(1e9)),
        with("bias", json!("x")),
        with("size", json!(-1)),
        with("tolerance", json!(1e-30)),
        with("body", json!(7)),
        with("sides", json!([{"body": "Dome", "face": "F4:-z"}])),
        with("source", json!([{"body": "Post", "face": "F1:-x"}, {"body": "Post", "face": "F1:+x"}])),
        // A curved face that is not meshed cannot be meshed here.
        with("source", top.clone().into()),
    ];
    let before = s.model.clone();
    for (i, p) in cases.iter().enumerate() {
        assert!(s.execute("mesh.hex.extrude", p).is_err(), "case {i}: {p}");
        assert_eq!(s.model, before, "case {i} changed the model");
    }
    assert!(s.execute("mesh.hex.extrude", &ok).is_ok());
}
