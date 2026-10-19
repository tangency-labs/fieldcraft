//! `mesh.mirror`: the mesh of some bodies reflected across a plane, joined to them on it.
//!
//! The geometry tags of the image name what it is the image of: the body of the open design that
//! occupies the reflected place (the other half of a part cut on its symmetry plane, found by its
//! volume and centroid, or named in `as`), and on it the face each reflected boundary face lies
//! on (found by position). Without such a body the image is `<body>/mirror` with the source's
//! face names.

use crate::schema::ProblemKind;
use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_geometry::{Design, Shape};
use fieldcraft_mesh::mirror::{mirror, reflected_edges, reflected_faces};
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use fieldcraft_mesh::tet::geometry::Tessellation;
use fieldcraft_model::geometry::Entity;
use fieldcraft_model::{ElementFace, GeometryRef, Model};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const MIRROR: CommandSpec = CommandSpec {
    id: "mesh.mirror",
    label: "Mirror Mesh",
    params: r#"{"plane": "Sym" (a datum plane) or {"point": [x, y, z], "normal": [nx, ny, nz]}, "bodies": ["Part/a"] (default: every body the mesh is tagged with), "as": {"Part/a": "Part/b"} (optional: the body each image is; default the design's body at the reflected place, else "<body>/mirror")} — replaces each image body's mesh (as mesh.tet.generate) with the reflection of the bodies' elements (tet4, tet10, wedge6, pyramid5, hex8, hex20) across the plane: each copy right-handed (its nodes in reflected order), the copies of the nodes on the plane merged with them so the mesh is conforming there, geometry tags carried to the image (its body, and on it the faces the reflected faces lie on, so sets by face name select both halves). Reports the elements and nodes added, the nodes shared on the plane and each image's body"#,
    form: Some(super::forms::MIRROR),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: mirror_mesh,
};

/// Relative difference of volume and centroid (to the body's size) under which a body is the
/// image of another.
const SAME_PLACE: f64 = 1e-4;

/// Distance (relative to the mesh's bounding-box diagonal) within which nodes are on the plane.
const ON_PLANE: f64 = 1e-7;

fn not_a_plane() -> EngineError {
    EngineError::invalid("plane", r#"must be a datum plane's name or {"point": [x, y, z], "normal": [x, y, z]}"#)
}

/// The plane: a datum plane by name, or a point and a normal.
fn plane(s: &Session, p: &Value) -> Result<([f64; 3], [f64; 3])> {
    match p.get("plane") {
        Some(Value::String(name)) => super::datum::datum_plane_of(&s.model, name),
        Some(v @ Value::Object(_)) => match (super::xyz(v, "point"), super::xyz(v, "normal")) {
            (Ok(point), Ok(normal)) => Ok((point, normal)),
            _ => Err(not_a_plane()),
        },
        None | Some(Value::Null) => Err(EngineError::param("plane", ProblemKind::Required)),
        _ => Err(not_a_plane()),
    }
}

/// The distance between two points.
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt()
}

fn reflect(x: [f64; 3], point: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    let d: f64 = (0..3).map(|i| (x[i] - point[i]) * n[i]).sum();
    [x[0] - 2.0 * d * n[0], x[1] - 2.0 * d * n[1], x[2] - 2.0 * d * n[2]]
}

/// The design body at the reflected place of `source`: the same volume and the reflected
/// centroid, within [`SAME_PLACE`] (the kernel measures curved bodies on a tessellation, so the
/// two halves of a symmetric part measure alike only to its accuracy).
fn image_body<'a>(bodies: &'a [Shape], source: &Shape, point: [f64; 3], n: [f64; 3]) -> Result<Option<&'a Shape>> {
    let m = source.measure()?;
    let size = source.size()?.max(1e-9);
    let want = reflect(m.centroid, point, n);
    for b in bodies.iter().filter(|b| b.name != source.name) {
        let o = b.measure()?;
        let d = (0..3).map(|k| (o.centroid[k] - want[k]).powi(2)).sum::<f64>().sqrt();
        if (o.volume - m.volume).abs() <= SAME_PLACE * m.volume.abs().max(1e-12) && d <= SAME_PLACE * size {
            return Ok(Some(b));
        }
    }
    Ok(None)
}

/// What `mesh.mirror` reflects, across which plane, into which bodies.
struct Plan<'d> {
    point: [f64; 3],
    /// The plane's unit normal.
    n: [f64; 3],
    /// The elements reflected.
    region: BTreeSet<u32>,
    /// Each mirrored body's image: its name and, when the design has it, its shape.
    image: BTreeMap<String, (String, Option<&'d Shape>)>,
}

fn plan<'d>(s: &Session, design: Option<&'d fieldcraft_geometry::Design>, p: &Value) -> Result<Plan<'d>> {
    let (point, normal) = plane(s, p)?;
    let len = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !(len > 1e-12 && len.is_finite() && point.iter().all(|x| x.is_finite())) {
        return Err(EngineError::Params("the plane needs a finite point and a non-zero normal".into()));
    }
    let n = normal.map(|x| x / len);
    let c = s.model.classification();
    // The bodies the mesh is tagged with, by tag.
    let tagged: BTreeMap<u32, String> = c
        .entities()
        .iter()
        .enumerate()
        .filter(|(_, r)| r.entity == Entity::Body)
        .filter_map(|(t, r)| Some((u32::try_from(t).ok()?, r.body.clone())))
        .collect();
    let bodies: Vec<String> = match p.get("bodies") {
        None => {
            let mut v: Vec<String> = tagged.values().cloned().collect();
            v.dedup();
            v
        }
        Some(_) => super::req(p, "bodies")?,
    };
    if bodies.is_empty() {
        return Err(EngineError::Params("no body to mirror: the mesh is not tagged with a body".into()));
    }
    for b in &bodies {
        if !tagged.values().any(|t| t == b) {
            return Err(EngineError::Params(format!("the mesh has no elements of body `{}`", b.chars().take(64).collect::<String>())));
        }
    }
    let named: BTreeMap<String, String> = super::opt(p, "as")?.unwrap_or_default();
    let region: BTreeSet<u32> =
        c.elements().iter().filter(|(_, t)| t.region.and_then(|r| tagged.get(&r)).is_some_and(|b| bodies.contains(b))).map(|(e, _)| *e).collect();

    // The image of each body: named, found in the design, or its own `/mirror`.
    let mut image: BTreeMap<String, (String, Option<&Shape>)> = BTreeMap::new();
    for b in &bodies {
        let source = design.and_then(|d| d.body(b));
        let target = match named.get(b) {
            Some(t) => (t.clone(), design.and_then(|d| d.body(t))),
            None => match (design, source) {
                (Some(d), Some(src)) => match image_body(&d.bodies, src, point, n)? {
                    Some(t) => (t.name.clone(), Some(t)),
                    None => (format!("{b}/mirror"), None),
                },
                _ => (format!("{b}/mirror"), None),
            },
        };
        if bodies.contains(&target.0) {
            return Err(EngineError::Params(format!("the image of `{b}` would be `{}`, which is mirrored too", target.0)));
        }
        image.insert(b.clone(), target);
    }
    Ok(Plan { point, n, region, image })
}

/// The bodies `mesh.mirror` fills: the images.
pub(super) fn image_names(s: &Session, p: &Value) -> Result<Vec<String>> {
    Ok(plan(s, s.design.as_deref(), p)?.image.into_values().map(|(name, _)| name).collect())
}

fn mirror_mesh(s: &mut Session, p: &Value) -> Result<Value> {
    let design = s.design.clone();
    let Plan { point, n, region, image } = plan(s, design.as_deref(), p)?;

    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for node in s.model.nodes().values() {
        for k in 0..3 {
            lo[k] = lo[k].min(node.xyz[k]);
            hi[k] = hi[k].max(node.xyz[k]);
        }
    }
    let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt().max(1e-9);
    // An image body's earlier mesh is replaced (replace.rs).
    let model = &mut s.model;
    let cleared = super::replace::clear(model, image.values().map(|(name, _)| name.as_str()));
    let before = model.clone();
    let r = mirror(model, &region, point, n, ON_PLANE * diagonal)?;
    tag_image(model, &before, &r, &image, design.as_deref(), point, n, diagonal)?;
    let added_nodes = model.nodes().len() - before.nodes().len();
    let images: Vec<Value> = image.iter().map(|(b, (t, shape))| json!({"body": b, "image": t, "in_design": shape.is_some()})).collect();
    let reply = json!({
        "elements": r.elements.len(),
        "nodes": added_nodes,
        "on_plane": r.on_plane,
        "images": images,
        "total": {"nodes": s.model.nodes().len(), "elements": s.model.elements().len()},
    });
    Ok(cleared.report(s, reply))
}

/// Tag the image: each copy with its image body, its faces with the faces of that body they lie
/// on (or the source's names), its edges and the images of tagged nodes with the source's names.
fn tag_image(
    model: &mut Model,
    before: &Model,
    r: &fieldcraft_mesh::mirror::Mirrored,
    image: &BTreeMap<String, (String, Option<&Shape>)>,
    design: Option<&Design>,
    point: [f64; 3],
    n: [f64; 3],
    diagonal: f64,
) -> Result<()> {
    let c = before.classification();
    let entity = |t: u32| c.entity(t).cloned();
    // Each image body's faces, to find the one each reflected boundary face lies on. A source
    // face's image is the set of image faces its own reflection covers (sampled at its
    // triangles' centres: one face for an exact mirror, pieces where the image is split
    // differently); each element face goes to the one of them nearest its reflected centre. So
    // an element face the mesher gave a face narrower than the elements (its centre on a
    // neighbour) keeps that face's image.
    let tolerance = (1e-3 * diagonal).max(1e-4);
    let mut targets: BTreeMap<&str, (Tessellation, Vec<String>)> = BTreeMap::new();
    for (target, shape) in image.values() {
        let Some(shape) = shape else { continue };
        let mesh = shape.tessellate(tolerance)?;
        let options = SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() };
        let surface = BoundarySurface::new(&mesh.positions, &mesh.triangles, &mesh.tri_face, &options)
            .map_err(|e| EngineError::context(format!("body `{target}`"), e))?;
        let names = shape.faces().iter().map(|f| f.name.clone()).collect();
        targets.insert(target.as_str(), (Tessellation::new(&surface, Vec::new()), names));
    }
    // The image faces each source face's reflection covers, by (source body, face name).
    let mut covers: BTreeMap<(String, String), Vec<u32>> = BTreeMap::new();
    for (body, (target, _)) in image {
        let (Some(shape), Some((tess, _))) = (design.and_then(|d| d.body(body)), targets.get(target.as_str())) else { continue };
        let mesh = shape.tessellate(tolerance)?;
        let at = |v: u32| mesh.positions.get(v as usize).copied().unwrap_or_default();
        for (t, &f) in mesh.triangles.iter().zip(&mesh.tri_face) {
            let Some(name) = shape.faces().get(f as usize).map(|x| x.name.clone()) else { continue };
            let centre = [0, 1, 2].map(|k| (at(t[0])[k] + at(t[1])[k] + at(t[2])[k]) / 3.0);
            let Some(hit) = tess.nearest_face(reflect(centre, point, n)) else { continue };
            let list = covers.entry((body.clone(), name)).or_default();
            if !list.contains(&hit) {
                list.push(hit);
            }
        }
    }
    let image_face = |e: u32, f: u8, body: &str, name: &str| -> Option<(String, String)> {
        let (target, _) = image.get(body)?;
        let mapped = targets.get(target.as_str()).and_then(|(tess, names)| {
            let nodes = before.face_nodes(ElementFace { element: e, face: f })?;
            let pts: Vec<[f64; 3]> = nodes.iter().filter_map(|n| before.node(*n).map(|x| x.xyz)).collect();
            let k = pts.len().max(1) as f64;
            let p = reflect([0, 1, 2].map(|i| pts.iter().map(|p| p[i]).sum::<f64>() / k), point, n);
            let face = match covers.get(&(body.to_string(), name.to_string())).map(Vec::as_slice) {
                Some([only]) => Some(*only),
                Some(several @ [_, _, ..]) => several
                    .iter()
                    .filter_map(|&g| tess.project_trimmed(g, p).map(|q| (dist(p, q), g)))
                    .min_by(|a, b| a.0.total_cmp(&b.0))
                    .map(|(_, g)| g),
                _ => tess.nearest_face(p),
            };
            face.and_then(|i| names.get(i as usize)).cloned()
        });
        Some((target.clone(), mapped.unwrap_or_else(|| name.to_string())))
    };
    let mut tags: BTreeMap<GeometryRef, u32> = BTreeMap::new();
    let mut tag = |model: &mut Model, r: GeometryRef| -> Result<u32> {
        if let Some(&t) = tags.get(&r) {
            return Ok(t);
        }
        let t = model.geometry_tag(&r)?;
        tags.insert(r, t);
        Ok(t)
    };
    for (&e, &copy) in &r.elements {
        let Some(el) = before.element(e) else { continue };
        let Some(t) = c.elements().get(&e) else { continue };
        let faces = reflected_faces(el.kind).unwrap_or_default();
        let edges = reflected_edges(el.kind).unwrap_or_default();
        if let Some(GeometryRef { body, .. }) = t.region.and_then(entity)
            && let Some((target, _)) = image.get(&body)
        {
            let tg = tag(model, GeometryRef::body(target.as_str())?)?;
            model.tag_element(copy, tg)?;
        }
        for (&f, &ft) in &t.faces {
            let Some(GeometryRef { body, entity: Entity::Face(name) }) = entity(ft) else { continue };
            let Some((tb, tn)) = image_face(e, f, &body, &name) else { continue };
            let g = faces.get(usize::from(f).saturating_sub(1)).copied().unwrap_or(f);
            let tg = tag(model, GeometryRef::face(tb.as_str(), tn.as_str())?)?;
            model.tag_element_face(copy, g, tg)?;
        }
        for (&k, &et) in &t.edges {
            let Some(GeometryRef { body, entity: ent }) = entity(et) else { continue };
            let Some((target, _)) = image.get(&body) else { continue };
            let g = edges.get(usize::from(k).saturating_sub(1)).copied().unwrap_or(k);
            let tg = tag(model, GeometryRef::new(target.as_str(), ent)?)?;
            model.tag_element_edge(copy, g, tg)?;
        }
    }
    for (&node, &copy) in &r.nodes {
        if copy == node {
            continue;
        }
        let Some(GeometryRef { body, entity: ent }) = c.nodes().get(&node).and_then(|&t| entity(t)) else { continue };
        let Some((target, _)) = image.get(&body) else { continue };
        let tg = tag(model, GeometryRef::new(target.as_str(), ent)?)?;
        model.tag_node(copy, tg)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A session with a 40 × 10 × 10 bar (along x from the origin) open, partitioned at
    /// x = 20 when `half`.
    fn bar(half: bool) -> (crate::tests::TempDir, Session) {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 40, "width": 10, "height": 10, "body_name": "Bar"}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let b = s.execute("geometry.bodies", &json!({})).unwrap()[0]["bbox"].clone();
        let lo: Vec<f64> = b[0].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        let hi: Vec<f64> = b[1].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        assert_eq!(lo, [0.0; 3]);
        assert_eq!(hi, [40.0, 10.0, 10.0]);
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 20}})).unwrap();
        if half {
            s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
        }
        (dir, s)
    }

    /// The face of `body` nearest to `point`.
    fn face(s: &Session, body: &str, point: [f64; 3]) -> String {
        let shape = s.design.as_ref().unwrap().body(body).unwrap();
        let f = shape.nearest_faces(&[point], 0.01).unwrap()[0].unwrap();
        shape.faces()[f as usize].name.clone()
    }

    /// Clamped at x = 0, 1 kN down on the end at x = 40 (sets by face name on the bodies given),
    /// solved; the largest downward deflection and the mesh's volume and counts.
    fn bend(s: &mut Session, root: (&str, String), tip: (&str, String)) -> (f64, f64, usize, usize) {
        let all: Vec<Value> = s.design.as_ref().unwrap().bodies.iter().map(|b| json!({"body": b.name})).collect();
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "geometry": all}},
            {"command": "model.set.create", "params": {"name": "root", "geometry": [{"body": root.0, "face": root.1}]}},
            {"command": "model.set.create", "params": {"name": "tip", "geometry": [{"body": tip.0, "face": tip.1}]}},
            {"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}},
            {"command": "model.section.create", "params": {"name": "Bar", "set": "all", "material": "Steel", "kind": {"type": "solid"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}},
            {"command": "model.load.create", "params": {"name": "Push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -1000]}},
            {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Push"]}},
        ]);
        s.run_script(&script).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        let tip_nodes = s.model.resolve_set("tip").unwrap().nodes;
        assert!(!tip_nodes.is_empty());
        let probe = s.execute("solve.probe", &json!({"step": "Load", "nodes": tip_nodes})).unwrap();
        let w = probe["nodes"].as_array().unwrap().iter().map(|n| -n["displacement"][2].as_f64().unwrap()).fold(f64::MIN, f64::max);
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(q["valid"].as_u64().unwrap() as usize, s.model.elements().len(), "{q}");
        (w, q["volume"].as_f64().unwrap(), s.model.nodes().len(), s.model.elements().len())
    }

    /// Half the bar meshed and mirrored is the whole bar meshed: the same counts and volume, and
    /// the same bending deflection to round-off, hex8 and hex20; the image's faces carry the
    /// whole bar's face names, so the load on the far end (a face of the image only) is found.
    #[test]
    fn a_mirrored_half_is_the_whole_mesh() {
        for order in [1, 2] {
            let (_d1, mut whole) = bar(false);
            whole.execute("mesh.hex.mapped", &json!({"body": "Bar", "size": 2.5, "order": order})).unwrap();
            let (root, tip) = (face(&whole, "Bar", [0.0, 5.0, 5.0]), face(&whole, "Bar", [40.0, 5.0, 5.0]));
            let full = bend(&mut whole, ("Bar", root), ("Bar", tip));

            let (_d2, mut half) = bar(true);
            half.execute("mesh.hex.mapped", &json!({"body": "Bar/a", "size": 2.5, "order": order})).unwrap();
            let r = half.execute("mesh.mirror", &json!({"plane": "Mid"})).unwrap();
            assert_eq!(r["images"][0], json!({"body": "Bar/a", "image": "Bar/b", "in_design": true}), "{r}");
            // 5 × 5 corner nodes on the plane (and 40 mid-edge nodes on it for hex20).
            assert_eq!(r["on_plane"], json!(if order == 1 { 25 } else { 65 }), "{r}");
            let (root, tip) = (face(&half, "Bar/a", [0.0, 5.0, 5.0]), face(&half, "Bar/b", [40.0, 5.0, 5.0]));
            let mirrored = bend(&mut half, ("Bar/a", root), ("Bar/b", tip));
            assert_eq!((mirrored.2, mirrored.3), (full.2, full.3), "order {order}: node and element counts");
            assert!((mirrored.1 - full.1).abs() < 1e-9 * full.1, "order {order}: volume {} vs {}", mirrored.1, full.1);
            assert!((mirrored.0 - full.0).abs() < 1e-9 * full.0, "order {order}: deflection {} vs {}", mirrored.0, full.0);
        }
    }

    /// Tets of either order mirror into valid elements sharing the plane's nodes; the volume
    /// doubles.
    #[test]
    fn mirrored_tets_are_valid_and_conforming() {
        for order in ["linear", "quadratic"] {
            let (_d, mut s) = bar(true);
            s.execute("mesh.tet.generate", &json!({"body": "Bar/a", "size": 4, "order": order})).unwrap();
            let v = s.execute("mesh.quality", &json!({})).unwrap()["volume"].as_f64().unwrap();
            let r = s.execute("mesh.mirror", &json!({"plane": {"point": [20, 0, 0], "normal": [1, 0, 0]}, "bodies": ["Bar/a"]})).unwrap();
            assert!(r["on_plane"].as_u64().unwrap() > 0, "{r}");
            let q = s.execute("mesh.quality", &json!({})).unwrap();
            assert_eq!(q["valid"].as_u64().unwrap() as usize, s.model.elements().len(), "{q}");
            assert!((q["volume"].as_f64().unwrap() - 2.0 * v).abs() < 1e-9 * v, "{q}");
            // The plane's faces are interior now: the skin of both bodies is the bar's.
            s.execute("model.set.create", &json!({"name": "a", "geometry": [{"body": "Bar/a"}]})).unwrap();
            s.execute("model.set.create", &json!({"name": "b", "geometry": [{"body": "Bar/b"}]})).unwrap();
            let (a, b) = (s.model.resolve_set("a").unwrap(), s.model.resolve_set("b").unwrap());
            assert_eq!(a.elements.len(), b.elements.len());
            assert_eq!(a.faces.len(), b.faces.len());
        }
    }

    #[test]
    fn hostile_params_change_nothing() {
        let (_d, mut s) = bar(true);
        s.execute("mesh.hex.mapped", &json!({"body": "Bar/a", "size": 5})).unwrap();
        let model = s.model.clone();
        for p in [
            json!({}),
            json!({"plane": "Nope"}),
            json!({"plane": {"point": [0, 0, 0], "normal": [0, 0, 0]}}),
            json!({"plane": {"point": [f64::MAX, 0, 0], "normal": [1, 0, 0]}}),
            json!({"plane": {"point": [0, 0], "normal": [1, 0, 0]}}),
            json!({"plane": "Mid", "bodies": ["Nope"]}),
            json!({"plane": "Mid", "bodies": []}),
            json!({"plane": "Mid", "as": {"Bar/a": "Bar/a"}}),
            json!({"plane": "Mid", "as": "x"}),
        ] {
            assert!(s.execute("mesh.mirror", &p).is_err(), "{p}");
        }
        assert_eq!(s.model, model);
    }
}
