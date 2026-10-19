//! Virtual topology: faces merged for meshing only (`geometry.faces.merge`), and tiny or sliver
//! faces absorbed into a neighbour (`geometry.face.absorb`, a merge that keeps the neighbour's
//! name).
//!
//! A face merge is part of the model (`fieldcraft_model::virtual_faces`); the geometry layer
//! applies it to the open design's bodies (`fieldcraft_geometry::virtual_faces`), after the
//! partitions, whenever the design is opened or cut again ([`super::datum::partitioned`]), so the
//! meshers, the viewport and the checks all see the merged face. The B-rep is not changed.

use fieldcraft_geometry::virtual_faces::merged_face_name;
use fieldcraft_model::FaceMerge;
use fieldcraft_model::geometry::{Entity, GeometryRef};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const FACES_MERGE: CommandSpec = CommandSpec {
    id: "geometry.faces.merge",
    label: "Merge Faces",
    params: r#"{"name": "Fillet", "faces": [{"body": "Base", "face": "F3:blend#0"}, {"body": "Base", "face": "F3:blend#1"}]} — makes two or more adjacent faces of one body one face for meshing only (virtual topology): the meshers, the checks and the viewport take them as one face `merge:<name>`, without the edges between them (an edge they leave in pieces along a neighbour is one edge), on the surface they share (none when they lie on different surfaces: the projection mesher and node relaxation keep to each B-rep piece's own surface, the other meshers to the tessellation). The design's B-rep is not changed. A face name also selects the pieces a partition cut it into, and may be a face an earlier merge made. The model keeps the merge and applies it again, by persistent face names, whenever the design is opened (after the partitions: on a body they cut, in each cell where its faces, or their pieces, touch); one that no longer applies (a face gone, faces that no longer touch) is listed with the design's errors. The body loses its mesh (`replaced`; remesh it); sets naming the merged faces no longer resolve (`stale_sets`: name `merge:<name>` instead). Faces on several bodies, faces that do not touch along an edge, or a face the body lacks are an error. One undo step"#,
    form: Some(super::forms::FACES_MERGE),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: merge,
};

pub const FACES_MERGE_DELETE: CommandSpec = CommandSpec {
    id: "geometry.faces.merge.delete",
    label: "Delete Face Merge",
    params: r#"{"name": "Fillet"} — removes a face merge: its faces are separate faces again (a later merge that names `merge:<name>` no longer applies and is listed with the design's errors); the body loses its mesh (`replaced`)"#,
    form: Some(super::forms::FACES_MERGE_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::FaceMerge],
    run: delete,
};

pub const FACE_ABSORB: CommandSpec = CommandSpec {
    id: "geometry.face.absorb",
    label: "Absorb Faces",
    params: r#"{"name": "Sliver", "face": {"body": "Base", "face": "F4:side#2"}, "into"?: "F1:+z" (default: the neighbour it shares the longest edges with)} | {"name"?: "Absorb", "from_check": true} — absorbs a tiny or sliver face into a neighbour for meshing only (virtual topology, as geometry.faces.merge): the meshers, the checks and the viewport take the two as one face that keeps the neighbour's name, so sets on the neighbour still select it and the tiny face no longer drives the mesh. The B-rep is not changed. `from_check`: every small and sliver face the last geometry.check found (and marked) on the open design, each into its best neighbour, as `<name>-1`, `<name>-2`, … in the check's order (a face already absorbed by an earlier one is `skipped`). Kept in the model and applied again by persistent names whenever the design is opened (delete one with geometry.faces.merge.delete); the bodies lose their meshes (`replaced`); sets naming an absorbed face are `stale_sets`. Replies `absorbed` ([{name, body, face, into}]) and `skipped`. One undo step"#,
    form: Some(super::forms::FACE_ABSORB),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: absorb,
};

pub const FACE_SPLIT: CommandSpec = CommandSpec {
    id: "geometry.face.split",
    label: "Split Face",
    params: r#"{"name": "Band", "face": {"body": "Base", "face": "F1:+z"}, "datum": "Mid" (a datum plane)} — splits a face along the curve where a datum plane crosses it, for meshing only (virtual topology; partitions cut whole bodies): the meshers, the checks and the viewport take it as two faces on its surface, `split:<name>:a` on the plane's negative side and `split:<name>:b` on the other, with the cut as an edge between them (a row of nodes along it; seed it like any edge), and the face's edges cut where the plane crosses them. The B-rep is not changed. Applied after the face merges, by persistent names, whenever the design is opened (a merge cannot name a split face); moving the datum moves the cut, and the datum cannot be deleted while a split uses it. The body loses its mesh (`replaced`); sets naming the face no longer resolve (`stale_sets`). A plane that does not cross the face is an error. Delete with geometry.face.split.delete. One undo step"#,
    form: Some(super::forms::FACE_SPLIT),
    kind: CommandKind::Model,
    requires: &[Requirement::Design, Requirement::DatumPlane],
    run: split,
};

pub const FACE_SPLIT_DELETE: CommandSpec = CommandSpec {
    id: "geometry.face.split.delete",
    label: "Delete Face Split",
    params: r#"{"name": "Band"} — removes a face split: the face is one face again; the body loses its mesh (`replaced`)"#,
    form: Some(super::forms::FACE_SPLIT_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::FaceSplit],
    run: split_delete,
};

/// Most faces one merge lists.
const MAX_FACES: usize = fieldcraft_model::virtual_faces::MAX_FACE_MERGE_FACES;

/// The body and the face names of `faces` (`[{body, face}]`, all on one body).
fn faces_of(p: &Value) -> Result<(String, Vec<String>)> {
    let refs: Vec<GeometryRef> = super::req(p, "faces")?;
    if !(2..=MAX_FACES).contains(&refs.len()) {
        return Err(EngineError::invalid("faces", format!("2 to {MAX_FACES} faces of one body")));
    }
    let body = refs.first().map(|r| r.body.clone()).unwrap_or_default();
    let mut names = Vec::with_capacity(refs.len());
    for r in &refs {
        match &r.entity {
            Entity::Face(f) if r.body == body => names.push(f.clone()),
            Entity::Face(_) => return Err(EngineError::invalid("faces", "faces of one body: a merge makes one face of one body")),
            _ => {
                return Err(EngineError::invalid("faces", format!("`{}`: give {{\"body\", \"face\"}}", r.body.chars().take(80).collect::<String>())));
            }
        }
    }
    Ok((body, names))
}

/// The open design's sets that name a face of `body` it no longer has.
fn stale_sets(s: &Session, body: &str) -> Vec<String> {
    let Some(shape) = s.design.as_deref().and_then(|d| d.body(body)) else { return Vec::new() };
    let gone = |r: &GeometryRef| r.body == body && matches!(&r.entity, Entity::Face(f) if shape.faces_named(f).is_empty());
    s.model.sets().iter().filter(|(_, set)| set.geometry.iter().any(gone)).map(|(n, _)| n.clone()).collect()
}

/// Apply the model's partitions, merges and splits to the design again, failing with the error of
/// one of `features` (`merge <name>`, `split <name>`) when it has one; the meshes of `bodies` go
/// (their faces changed).
pub(super) fn reapply(s: &mut Session, features: &[String], bodies: &[String], reply: Value) -> Result<Value> {
    let mut cleared = super::datum::recut(s)?;
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let errors: Vec<String> =
        design.errors.iter().filter(|e| features.contains(&e.feature)).map(|e| format!("{}: {}", e.feature, e.message)).collect();
    if !errors.is_empty() {
        return Err(EngineError::Params(errors.join("; ")));
    }
    cleared.extend(super::replace::clear(&mut s.model, bodies.iter().map(String::as_str)));
    let stale: Vec<String> = bodies.iter().flat_map(|b| stale_sets(s, b)).collect();
    let mut reply = cleared.report_if_any(s, reply);
    if !stale.is_empty()
        && let Some(o) = reply.as_object_mut()
    {
        let mut all: Vec<String> = o.get("stale_sets").and_then(|v| serde_json::from_value(v.clone()).ok()).unwrap_or_default();
        all.extend(stale);
        all.sort();
        all.dedup();
        o.insert("stale_sets".into(), json!(all));
    }
    Ok(reply)
}

fn merge(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    let (body, faces) = faces_of(p)?;
    if s.design.as_deref().and_then(|d| d.body(&body)).is_none() {
        return Err(EngineError::Params(format!("the open design has no body `{}`", body.chars().take(200).collect::<String>())));
    }
    s.model.add_face_merge(FaceMerge { name: name.clone(), body: body.clone(), faces: faces.clone(), into: None })?;
    let face = merged_face_name(&name);
    let reply = json!({"name": name, "body": body, "face": face});
    let mut reply = reapply(s, &[format!("merge {name}")], std::slice::from_ref(&body), reply)?;
    let shape = s.design.as_deref().and_then(|d| d.body(&body));
    let brep = shape.and_then(|b| b.face_index(&face)).and_then(|f| shape?.merged_faces().into_iter().find(|(g, _)| *g == f)).map_or(0, |(_, n)| n);
    if let Some(o) = reply.as_object_mut() {
        o.insert("brep_faces".into(), json!(brep));
    }
    Ok(reply)
}

/// The faces to absorb, by body and name: the one named, or the small and sliver faces of the
/// last geometry.check.
fn to_absorb(s: &Session, p: &Value) -> Result<Vec<(String, String)>> {
    let from_check: bool = super::opt(p, "from_check")?.unwrap_or(false);
    match (p.get("face"), from_check) {
        (Some(_), false) => {
            let r: GeometryRef = super::req(p, "face")?;
            match r.entity {
                Entity::Face(f) => Ok(vec![(r.body, f)]),
                _ => Err(EngineError::invalid("face", "give {\"body\", \"face\"}")),
            }
        }
        (None, true) => {
            let marks = crate::view::geometry_marks(s)
                .ok_or_else(|| EngineError::Params("`from_check`: run geometry.check (with marks) on the open design first".into()))?;
            let mut faces: Vec<(String, String)> = Vec::new();
            for f in marks.marks.iter().filter_map(|m| m.small.clone()) {
                if !faces.contains(&f) {
                    faces.push(f);
                }
            }
            if faces.is_empty() {
                return Err(EngineError::Params("`from_check`: the last geometry.check found no small or sliver face".into()));
            }
            Ok(faces)
        }
        _ => Err(EngineError::Params("give either `face` or `from_check: true`".into())),
    }
}

fn absorb(s: &mut Session, p: &Value) -> Result<Value> {
    let faces = to_absorb(s, p)?;
    let one = faces.len() == 1 && p.get("face").is_some();
    let name = match p.get("name") {
        None if !one => "Absorb".to_string(),
        _ => super::name(p)?.to_string(),
    };
    let into: Option<String> = super::opt(p, "into")?;
    if into.is_some() && !one {
        return Err(EngineError::invalid("into", "goes with `face`"));
    }
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    // The bodies as each absorption leaves them, to find the next one's neighbour.
    let mut shapes: Vec<fieldcraft_geometry::Shape> = Vec::new();
    let (mut absorbed, mut skipped, mut merges) = (Vec::new(), Vec::new(), Vec::new());
    for (k, (body, face)) in faces.into_iter().enumerate() {
        if !shapes.iter().any(|b| b.name == body) {
            let shape = design
                .body(&body)
                .ok_or_else(|| EngineError::Params(format!("the open design has no body `{}`", body.chars().take(200).collect::<String>())))?;
            shapes.push(shape.clone());
        }
        let Some(shape) = shapes.iter_mut().find(|b| b.name == body) else { continue };
        let target = match &into {
            Some(i) => Ok(i.clone()),
            None => shape.absorbing_neighbour(&face),
        };
        let merged = target.map_err(EngineError::from).and_then(|t| Ok((shape.merge_faces(&t, &[face.clone(), t.clone()])?, t)));
        let (next, target) = match merged {
            Ok(m) => m,
            Err(e) if !one => {
                skipped.push(json!({"body": body, "face": face, "reason": e.to_string()}));
                continue;
            }
            Err(e) => return Err(e),
        };
        *shape = next;
        let merge = if one { name.clone() } else { format!("{name}-{}", k + 1) };
        s.model.add_face_merge(FaceMerge {
            name: merge.clone(),
            body: body.clone(),
            faces: vec![face.clone(), target.clone()],
            into: Some(target.clone()),
        })?;
        absorbed.push(json!({"name": merge, "body": body, "face": face, "into": target}));
        merges.push(merge);
    }
    if merges.is_empty() {
        return Err(EngineError::Params(format!("no face could be absorbed: {}", json!(skipped))));
    }
    let bodies: Vec<String> = shapes.into_iter().map(|b| b.name).collect();
    let features: Vec<String> = merges.iter().map(|m| format!("merge {m}")).collect();
    reapply(s, &features, &bodies, json!({"absorbed": absorbed, "skipped": skipped}))
}

fn split(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?.to_string();
    let r: GeometryRef = super::req(p, "face")?;
    let Entity::Face(face) = r.entity else { return Err(EngineError::invalid("face", "give {\"body\", \"face\"}")) };
    let datum: String = super::req(p, "datum")?;
    super::datum::datum_plane_of(&s.model, &datum)?;
    s.model.add_face_split(fieldcraft_model::FaceSplit { name: name.clone(), body: r.body.clone(), face: face.clone(), datum })?;
    let feature = format!("split:{name}");
    let faces = [format!("{feature}:a"), format!("{feature}:b")];
    let reply = json!({"name": name, "body": r.body, "faces": faces});
    reapply(s, &[format!("split {name}")], std::slice::from_ref(&r.body), reply)
}

fn split_delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    let removed = s.model.remove_face_split(name)?;
    let mut cleared = super::datum::recut(s)?;
    cleared.extend(super::replace::clear(&mut s.model, [removed.body.as_str()]));
    Ok(cleared.report_if_any(s, json!({"deleted": name, "body": removed.body})))
}

fn delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    let removed = s.model.remove_face_merge(name)?;
    let mut cleared = super::datum::recut(s)?;
    cleared.extend(super::replace::clear(&mut s.model, [removed.body.as_str()]));
    Ok(cleared.report_if_any(s, json!({"deleted": name, "body": removed.body})))
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use fieldcraft_geometry::FaceSurface;
    use serde_json::json;

    use crate::Session;

    /// The plate's hole halves (two half-cylinders), by name.
    fn halves(s: &Session) -> Vec<String> {
        let b = &s.design.as_ref().unwrap().bodies[0];
        b.faces()
            .iter()
            .zip(b.face_surfaces().unwrap())
            .filter(|(_, x)| matches!(x, Some(FaceSurface::Cylinder { .. })))
            .map(|(f, _)| f.name.clone())
            .collect()
    }

    fn merge_hole(s: &mut Session) -> serde_json::Value {
        let body = s.design.as_ref().unwrap().bodies[0].name.clone();
        let faces: Vec<_> = halves(s).into_iter().map(|f| json!({"body": body, "face": f})).collect();
        s.execute("geometry.faces.merge", &json!({"name": "Hole", "faces": faces})).unwrap()
    }

    fn faces(s: &Session) -> Vec<String> {
        s.design.as_ref().unwrap().bodies[0].faces().iter().map(|f| f.name.clone()).collect()
    }

    /// A merge makes the hole one face of the design; undo and redo take it out and put it back;
    /// it is saved with the model and applied again when the model and its design are opened.
    #[test]
    fn a_merge_is_one_face_undone_saved_and_replayed() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let before = faces(&s);
        assert_eq!(before.len(), 8);
        let r = merge_hole(&mut s);
        assert_eq!((r["face"].as_str(), r["brep_faces"].as_u64()), (Some("merge:Hole"), Some(2)), "{r}");
        assert_eq!(faces(&s).len(), 7);
        assert!(faces(&s).contains(&"merge:Hole".to_string()));
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(faces(&s), before);
        s.execute("edit.redo", &json!({})).unwrap();
        assert_eq!(faces(&s).len(), 7);
        // Saved, and replayed on open.
        let file = dir.join("merged.fieldcraft");
        s.execute("file.save", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(t.model.face_merges().len(), 1);
        assert_eq!(faces(&t), faces(&s));
        // Opening the design again (after an edit of it) replays the merge.
        let design = s.sources.design.as_ref().unwrap().path().to_string();
        s.execute("geometry.open", &json!({"path": design})).unwrap();
        assert_eq!(faces(&s).len(), 7);
        // Deleting it brings the halves back.
        s.execute("geometry.faces.merge.delete", &json!({"name": "Hole"})).unwrap();
        assert_eq!(faces(&s), before);
        std::fs::remove_dir_all(dir).ok();
    }

    /// The mesh nodes on the seams between the plate's hole halves (the B-rep edges between them).
    fn nodes_on_seams(s: &Session, seams: &[Vec<[f64; 3]>]) -> usize {
        let on = |p: [f64; 3]| {
            seams.iter().any(|line| {
                line.windows(2).any(|w| {
                    let (a, b) = (w[0], w[1]);
                    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
                    let l2 = ab.iter().map(|x| x * x).sum::<f64>();
                    let t = (ab.iter().zip(&ap).map(|(x, y)| x * y).sum::<f64>() / l2).clamp(0.0, 1.0);
                    (0..3).map(|k| (ap[k] - t * ab[k]).powi(2)).sum::<f64>().sqrt() < 1e-6
                })
            })
        };
        // In millimetres, as the design is.
        s.model.nodes().values().filter(|n| on(crate::units::to_geometry_point(n.xyz))).count()
    }

    /// The sample plate made square (30 × 30 × 10, the Ø10 hole in the middle): no layout of the
    /// cross-section has a reason to cut along the seams.
    fn square_plate() -> (Session, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("fieldcraft-merge-{}-{:?}", std::process::id(), std::thread::current().id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("square.solidcraft");
        let plate = fieldcraft_geometry::sample::Plate { width: 30.0, ..Default::default() };
        std::fs::write(&path, fieldcraft_geometry::sample::plate_design(&plate).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path.to_str().unwrap()})).unwrap();
        (s, dir)
    }

    /// A mesher takes the merged face as one: its boundary element faces are tagged with it, the
    /// halves are gone from the tags, and the seams between the halves (no longer edges) carry
    /// `expected(session, seams)` nodes.
    fn mesher_meshes_a_merged_face_as_one(command: &str, params: serde_json::Value, expected: impl Fn(&Session, &[Vec<[f64; 3]>]) -> usize) {
        let (mut s, dir) = square_plate();
        let hole = halves(&s);
        let b = &s.design.as_ref().unwrap().bodies[0];
        let seams: Vec<Vec<[f64; 3]>> = b
            .edges()
            .unwrap()
            .into_iter()
            .filter(|e| e.faces.len() == 2 && e.faces[0] != e.faces[1] && e.faces.iter().all(|&f| hole.contains(&b.faces()[f as usize].name)))
            .map(|e| e.points)
            .collect();
        assert_eq!(seams.len(), 2);
        s.execute(command, &params).unwrap();
        let before = nodes_on_seams(&s, &seams);
        assert!(before >= 6, "{command}: a node row on each seam before the merge: {before}");
        let r = merge_hole(&mut s);
        assert!(r["replaced"].as_array().is_some_and(|a| !a.is_empty()), "the mesh goes: {r}");
        assert!(s.model.elements().is_empty());
        s.execute(command, &params).unwrap();
        // The face, edge and vertex names the elements and nodes are tagged with.
        let c = s.model.classification();
        let used = c.elements().values().flat_map(|t| t.faces.values().chain(t.edges.values())).chain(c.nodes().values());
        let tagged: BTreeSet<String> = used.filter_map(|&t| c.entity(t)?.entity.name().map(str::to_string)).collect();
        assert!(tagged.contains("merge:Hole"), "{tagged:?}");
        assert!(!tagged.iter().any(|n| hole.iter().any(|h| n.split(['|', '~']).any(|p| p == h))), "{tagged:?}");
        let after = nodes_on_seams(&s, &seams);
        let expected = expected(&s, &seams);
        assert_eq!(after, expected, "{command}: nodes on the seams after the merge (before {before})");
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn the_tet_mesher_meshes_a_merged_face_as_one() {
        // Before: the surface mesh keeps the seams as feature lines; after: none.
        mesher_meshes_a_merged_face_as_one("mesh.tet.generate", json!({"size": 4}), |_, _| 0);
    }

    #[test]
    fn the_hex_mesher_meshes_a_merged_face_as_one() {
        // Before, the swept O-grid round the hole had a node column on each seam, and after the
        // merge one on the seam where the closed rim's chain started (#512). The rim's start is
        // no curve end now: the seams are lines like any other round the hole. The square plate's
        // layout (cuts from its corners) still puts nodes in the middle of the arcs between the
        // cuts, which on this plate are where the seams lie and also a quarter turn round: the
        // seams carry as many nodes as those turned lines, no more and no fewer.
        let turned = |s: &Session, seams: &[Vec<[f64; 3]>]| {
            let quarter: Vec<Vec<[f64; 3]>> = seams.iter().map(|l| l.iter().map(|p| [30.0 - p[1], p[0], p[2]]).collect()).collect();
            nodes_on_seams(s, &quarter)
        };
        mesher_meshes_a_merged_face_as_one("mesh.hex.generate", json!({"size": 2.5}), turned);
    }

    /// The hex-meshing plan sees the merged hole as one face: the hole's feature and its ring
    /// project onto `merge:Hole` (still a cylinder), not the halves.
    #[test]
    fn the_plan_sees_a_merged_face_as_one() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let body = s.design.as_ref().unwrap().bodies[0].name.clone();
        let hole = halves(&s);
        let before = s.execute("mesh.plan", &json!({"body": body})).unwrap();
        let ring = |r: &serde_json::Value| r["rings"][0]["faces"].clone();
        let mut halves_sorted = hole.clone();
        halves_sorted.sort();
        let mut got: Vec<String> = serde_json::from_value(ring(&before)).unwrap();
        got.sort();
        assert_eq!(got, halves_sorted, "{before}");
        merge_hole(&mut s);
        let after = s.execute("mesh.plan", &json!({"body": body})).unwrap();
        assert_eq!(ring(&after), json!(["merge:Hole"]), "{after}");
        let features = after["features"].as_array().unwrap();
        assert!(features.iter().any(|f| f["faces"] == json!(["merge:Hole"])), "{after}");
        assert!(!after.to_string().contains(&hole[0]), "{after}");
        std::fs::remove_dir_all(dir).ok();
    }

    /// The viewport draws a merged face in one colour: the hole's triangles, and only those.
    #[test]
    fn a_merged_face_is_drawn_in_one_colour() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let p = crate::view::Palette::LIGHT;
        let merged = |s: &Session| crate::view::scene(s, &p).0.tris.into_iter().filter(|t| t.color == p.merged).collect::<Vec<_>>();
        assert!(merged(&s).is_empty());
        merge_hole(&mut s);
        let tris = merged(&s);
        assert!(!tris.is_empty());
        // The hole: Ø10 about (15, 15).
        // The scene is in metres.
        let mm = crate::units::to_geometry_length;
        for t in &tris {
            assert!(t.p.iter().all(|v| ((mm(v.x) - 15.0).hypot(mm(v.y) - 15.0) - 5.0).abs() < 1e-3), "{t:?}");
        }
        std::fs::remove_dir_all(dir).ok();
    }

    /// Sets naming a merged face go stale (and are named); bad merges change nothing.
    #[test]
    fn bad_merges_change_nothing() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let body = s.design.as_ref().unwrap().bodies[0].name.clone();
        let names = faces(&s);
        let hole = halves(&s);
        s.execute("model.set.create", &json!({"name": "Bore", "geometry": [{"body": body, "face": hole[0]}]})).unwrap();
        let journal = s.journal.len();
        let f = |n: &str| json!({"body": body, "face": n});
        for p in [
            json!({}),
            json!({"name": "M"}),
            json!({"name": "M", "faces": [f(&hole[0])]}),
            json!({"name": "M", "faces": [f(&hole[0]), {"body": "Other", "face": hole[1]}]}),
            json!({"name": "M", "faces": [{"body": "Other", "face": "a"}, {"body": "Other", "face": "b"}]}),
            json!({"name": "M", "faces": [f(&hole[0]), f("nope")]}),
            json!({"name": "M", "faces": [f(&hole[0]), {"body": body}]}),
            json!({"name": "M", "faces": [f(&hole[0]), f(&hole[0])]}),
            json!({"name": "1M", "faces": [f(&hole[0]), f(&hole[1])]}),
            json!({"name": "M", "faces": "x"}),
        ] {
            assert!(s.execute("geometry.faces.merge", &p).is_err(), "{p}");
        }
        // Top and bottom do not touch.
        let b = &s.design.as_ref().unwrap().bodies[0];
        let caps: Vec<String> = b
            .faces()
            .iter()
            .zip(b.face_surfaces().unwrap())
            .filter(|(_, x)| matches!(x, Some(FaceSurface::Plane { normal, .. }) if normal[2].abs() > 0.5))
            .map(|(f, _)| f.name.clone())
            .collect();
        let e = s.execute("geometry.faces.merge", &json!({"name": "Caps", "faces": [f(&caps[0]), f(&caps[1])]})).unwrap_err();
        assert!(e.to_string().contains("does not touch"), "{e}");
        assert!(s.execute("geometry.faces.merge.delete", &json!({"name": "Caps"})).is_err());
        assert_eq!((s.journal.len(), faces(&s)), (journal, names));
        assert!(s.model.face_merges().is_empty());
        let r = merge_hole(&mut s);
        assert_eq!(r["stale_sets"], json!(["Bore"]), "{r}");
        std::fs::remove_dir_all(dir).ok();
    }

    /// A block with a 0.4 mm step on its top: the riser is a sliver at size 5.
    fn stepped() -> (crate::tests::TempDir, Session) {
        super::super::geometry::tests::design(json!([
            {"command": "solid.box", "params": {"length": 40, "width": 20, "height": 10, "body_name": "Base"}},
            {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 10.4, "operation": "join", "targets": ["Base"]}},
        ]))
    }

    fn small_faces(s: &mut Session) -> usize {
        let r = s.execute("geometry.check", &json!({"size": 5})).unwrap();
        r["findings"].as_array().unwrap().iter().filter(|f| matches!(f["kind"].as_str(), Some("small_face" | "sliver_face"))).count()
    }

    /// `from_check` absorbs every sliver the check found into its best neighbour, which keeps its
    /// name; the check then finds none, the tet mesher tags no element face with the sliver, and
    /// undo brings it back.
    #[test]
    fn absorbing_the_check_findings_clears_them() {
        let (_dir, mut s) = stepped();
        let before = small_faces(&mut s);
        assert!(before > 0);
        let start = faces(&s);
        assert!(s.execute("geometry.face.absorb", &json!({"from_check": true, "into": "x"})).is_err());
        let r = s.execute("geometry.face.absorb", &json!({"from_check": true})).unwrap();
        let absorbed = r["absorbed"].as_array().unwrap();
        assert_eq!(absorbed.len(), before, "{r}");
        assert_eq!(absorbed[0]["name"], "Absorb-1");
        let (gone, into) = (absorbed[0]["face"].as_str().unwrap().to_string(), absorbed[0]["into"].as_str().unwrap().to_string());
        assert!(!faces(&s).contains(&gone) && faces(&s).contains(&into));
        assert_eq!(faces(&s).len(), start.len() - before);
        assert_eq!(small_faces(&mut s), 0);
        s.execute("mesh.tet.generate", &json!({"size": 2})).unwrap();
        let c = s.model.classification();
        let used: BTreeSet<String> =
            c.elements().values().flat_map(|t| t.faces.values()).filter_map(|&t| c.entity(t)?.entity.name().map(str::to_string)).collect();
        assert!(used.contains(&into) && !used.contains(&gone), "{used:?}");
        s.execute("edit.undo", &json!({})).unwrap();
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(faces(&s), start);
    }

    /// One face into a named neighbour, or its best one; bad requests change nothing.
    #[test]
    fn absorbing_one_face() {
        let (_dir, mut s) = stepped();
        let r = s.execute("geometry.check", &json!({"size": 5})).unwrap();
        let sliver =
            r["findings"].as_array().unwrap().iter().find(|f| f["kind"] == "sliver_face" || f["kind"] == "small_face").unwrap()["faces"][0].clone();
        let face = json!({"body": "Base", "face": sliver});
        let journal = s.journal.len();
        for p in [
            json!({"name": "A"}),
            json!({"name": "A", "face": face, "from_check": true}),
            json!({"face": face}),
            json!({"name": "A", "face": {"body": "Base", "face": "nope"}}),
            json!({"name": "A", "face": face, "into": "nope"}),
            json!({"name": "A", "face": {"body": "Nobody", "face": "x"}}),
            json!({"name": "A", "face": face, "from_check": "yes"}),
        ] {
            assert!(s.execute("geometry.face.absorb", &p).is_err(), "{p}");
        }
        assert_eq!(s.journal.len(), journal);
        let r = s.execute("geometry.face.absorb", &json!({"name": "Riser", "face": face})).unwrap();
        assert_eq!(r["absorbed"][0]["name"], "Riser", "{r}");
        assert!(s.model.face_merges()[0].into.is_some());
    }

    /// The plate's top face (z = 10), by name.
    fn top(s: &Session) -> String {
        let b = &s.design.as_ref().unwrap().bodies[0];
        b.faces()
            .iter()
            .zip(b.face_surfaces().unwrap())
            .find(|(_, x)| matches!(x, Some(FaceSurface::Plane { normal, point }) if normal[2] > 0.5 && point[2] > 9.0))
            .unwrap()
            .0
            .name
            .clone()
    }

    /// Mesh nodes on the top along x = `x`.
    fn nodes_along(s: &Session, x: f64) -> usize {
        let mm = crate::units::to_geometry_point;
        s.model.nodes().values().filter(|n| (mm(n.xyz)[0] - x).abs() < 1e-6 && (mm(n.xyz)[2] - 10.0).abs() < 1e-6).count()
    }

    /// Splitting the top along a datum plane gives two faces the tet mesher tags, with a row of
    /// nodes along the cut; the split follows the datum, holds it, is saved and replayed, and
    /// undoes.
    #[test]
    fn a_split_face_meshes_with_a_node_row_along_the_cut() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let body = s.design.as_ref().unwrap().bodies[0].name.clone();
        let face = json!({"body": body, "face": top(&s)});
        s.execute("geometry.datum", &json!({"name": "Cut", "by": {"type": "offset", "plane": "yz", "offset": 30.3}})).unwrap();
        let mesh = json!({"size": 3});
        s.execute("mesh.tet.generate", &mesh).unwrap();
        let before = nodes_along(&s, 30.3);
        let r = s.execute("geometry.face.split", &json!({"name": "Band", "face": face, "datum": "Cut"})).unwrap();
        assert_eq!(r["faces"], json!(["split:Band:a", "split:Band:b"]), "{r}");
        assert!(r["replaced"].is_array(), "{r}");
        assert_eq!(faces(&s).len(), 9);
        s.execute("mesh.tet.generate", &mesh).unwrap();
        let after = nodes_along(&s, 30.3);
        assert!(after >= 10 && after > before, "a node row along the 30 mm cut: {before} -> {after}");
        let c = s.model.classification();
        let used: BTreeSet<String> =
            c.elements().values().flat_map(|t| t.faces.values()).filter_map(|&t| c.entity(t)?.entity.name().map(str::to_string)).collect();
        assert!(used.contains("split:Band:a") && used.contains("split:Band:b"), "{used:?}");
        // The datum holds; moving it moves the cut.
        assert!(s.execute("geometry.datum.delete", &json!({"name": "Cut"})).is_err());
        s.execute("geometry.datum.edit", &json!({"name": "Cut", "by": {"type": "offset", "plane": "yz", "offset": 35.5}})).unwrap();
        s.execute("mesh.tet.generate", &mesh).unwrap();
        assert!(nodes_along(&s, 35.5) >= 10);
        // Saved and replayed.
        let file = dir.join("split.fieldcraft");
        s.execute("file.save", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(faces(&t), faces(&s));
        // Bad splits change nothing; delete brings the face back.
        let journal = s.journal.len();
        for p in [
            json!({"name": "X", "face": face, "datum": "nope"}),
            json!({"name": "X", "face": {"body": body, "face": "nope"}, "datum": "Cut"}),
            json!({"name": "X", "face": {"body": body}, "datum": "Cut"}),
            json!({"name": "Band", "face": face, "datum": "Cut"}),
            json!({"name": "X", "face": {"body": body, "face": "split:Band:b"}, "datum": "Cut"}),
        ] {
            assert!(s.execute("geometry.face.split", &p).is_err(), "{p}");
        }
        assert_eq!(s.journal.len(), journal);
        s.execute("geometry.face.split.delete", &json!({"name": "Band"})).unwrap();
        assert_eq!(faces(&s).len(), 8);
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(faces(&s).len(), 9);
        std::fs::remove_dir_all(dir).ok();
    }
}
