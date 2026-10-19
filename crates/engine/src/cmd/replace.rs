//! Meshing a body replaces its mesh (#267). Every mesher that fills bodies goes through
//! [`clear`]: before it adds a body's new elements it removes that body's old elements and the
//! nodes only they use, so other bodies keep their meshes and two meshes never stack on one
//! body. The reply's `replaced` says what went; sets that listed the removed labels are named as
//! `stale_sets` (sets by geometry re-resolve on the new mesh), and results solved on the old
//! mesh as `results_cleared`. The command is one undo step, so undo brings the old mesh back.
//!
//! [`MESHERS`] lists those commands with the bodies each one fills, which the UI reads
//! ([`replacement`]) to ask before replacing a mesh. Programmatic calls never ask. `mesh.block`
//! is not one of them: a body is built from several blocks, each call adding one.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{BodyMesh, ElementKind, Model};
use serde_json::{Value, json};

use crate::{EngineError, Result, Session};

/// What [`clear`] removed.
#[derive(Default)]
pub(super) struct Cleared {
    /// Per body that had a mesh: its name and the element and node counts removed.
    bodies: Vec<(String, usize, usize)>,
    removed: BodyMesh,
    /// Sets that listed a removed node or element by label.
    stale_sets: Vec<String>,
}

/// Remove the meshes of `bodies` from `model` (each body once).
pub(super) fn clear<'a>(model: &mut Model, bodies: impl IntoIterator<Item = &'a str>) -> Cleared {
    let mut out = Cleared::default();
    let mut seen = BTreeSet::new();
    for body in bodies {
        if !seen.insert(body) {
            continue;
        }
        let r = model.remove_body_mesh(body);
        if !r.elements.is_empty() {
            out.bodies.push((body.to_string(), r.elements.len(), r.nodes.len()));
            out.removed.elements.extend(r.elements);
            out.removed.nodes.extend(r.nodes);
        }
    }
    out.stale_sets = model.sets_using(&out.removed);
    out
}

impl Cleared {
    /// Add what `other` removed (a body in both counts once, as `other` reports it).
    pub(super) fn extend(&mut self, other: Cleared) {
        self.bodies.retain(|(b, _, _)| !other.bodies.iter().any(|(o, _, _)| o == b));
        self.bodies.extend(other.bodies);
        self.removed.elements.extend(other.removed.elements);
        self.removed.nodes.extend(other.removed.nodes);
        self.stale_sets.extend(other.stale_sets);
        self.stale_sets.sort();
        self.stale_sets.dedup();
    }

    /// [`Cleared::report`] when anything was removed; `reply` as it is otherwise.
    pub(super) fn report_if_any(self, s: &Session, reply: Value) -> Value {
        if self.bodies.is_empty() { reply } else { self.report(s, reply) }
    }

    /// Add `replaced`, `stale_sets` and `results_cleared` to a command's reply (an object).
    pub(super) fn report(self, s: &Session, mut reply: Value) -> Value {
        if let Some(o) = reply.as_object_mut() {
            let replaced: Vec<Value> = self.bodies.iter().map(|(b, e, n)| json!({"body": b, "elements": e, "nodes": n})).collect();
            o.insert("replaced".into(), json!(replaced));
            if !self.stale_sets.is_empty() {
                o.insert("stale_sets".into(), json!(self.stale_sets));
            }
            let results: Vec<&str> = s.solved_steps().collect();
            if !results.is_empty() {
                o.insert("results_cleared".into(), json!(results));
            }
            // A seed or size control on what the design lacks shaped nothing: never in silence.
            let skipped: Vec<Value> = s
                .design
                .as_deref()
                .map(|d| super::names::stale_uses(&s.model, d))
                .unwrap_or_default()
                .into_iter()
                .filter(|u| u.owner == "edge seed" || u.owner.starts_with("size control"))
                .map(|u| json!(format!("{} was not applied: the design has no {}; pick it again", u.owner, u.reference)))
                .collect();
            if !skipped.is_empty() {
                let warnings = o.entry("warnings").or_insert_with(|| json!([]));
                if let Some(w) = warnings.as_array_mut() {
                    w.extend(skipped);
                }
            }
        }
        reply
    }
}

pub const CLEAR: crate::CommandSpec = crate::CommandSpec {
    id: "mesh.clear",
    label: "Clear Mesh",
    params: r#"{"bodies"?: ["Body1/a"]} — removes the meshes of the bodies (their elements and the nodes only they use; sets that listed them are named as `stale_sets`, results solved on them are cleared); without `bodies`, of every body the open design no longer has because its partitions changed since it was meshed (model.check lists them; none: nothing changes). Reports what went as `replaced`"#,
    form: Some(super::forms::CLEAR_MESH),
    kind: crate::CommandKind::Model,
    requires: &[crate::requirement::Requirement::Mesh],
    run: clear_meshes,
};

fn clear_meshes(s: &mut Session, p: &Value) -> Result<Value> {
    let bodies: Vec<String> = match p.get("bodies") {
        Some(_) => super::req(p, "bodies")?,
        None => s.design.as_deref().map(|d| stale(&s.model, d).into_iter().map(|(b, _)| b).collect()).unwrap_or_default(),
    };
    for b in &bodies {
        if s.model.body_elements(b).is_empty() {
            return Err(EngineError::Params(format!("body `{}` has no mesh", b.chars().take(64).collect::<String>())));
        }
    }
    let cleared = clear(&mut s.model, bodies.iter().map(String::as_str));
    let (nodes, elements) = (s.model.nodes().len(), s.model.elements().len());
    Ok(cleared.report(s, json!({"nodes": nodes, "elements": elements})))
}

/// The meshed bodies the open design no longer has because its partitions changed (#385): cells
/// of a partition that is gone (`Base/a/b`), or a body partitions now cut (`Base`), each with
/// its element count. Other bodies the design lacks (mirror images, imported surfaces, blocks)
/// are not stale.
pub(crate) fn stale(model: &Model, design: &fieldcraft_geometry::Design) -> Vec<(String, usize)> {
    let cell = |name: &str| -> bool {
        let Some(source) = design.source_of(name) else { return false };
        let rest = name.get(source.len()..).unwrap_or_default();
        rest.split('/').skip(1).all(|seg| seg.strip_prefix(['a', 'b']).is_some_and(|k| k.bytes().all(|c| c.is_ascii_digit())))
    };
    model
        .classification()
        .entities()
        .iter()
        .filter(|g| matches!(g.entity, fieldcraft_model::geometry::Entity::Body))
        .map(|g| g.body.as_str())
        .filter(|b| design.body(b).is_none() && cell(b))
        .filter_map(|b| {
            let n = model.body_elements(b).len();
            (n > 0).then(|| (b.to_string(), n))
        })
        .collect()
}

/// The bodies a meshing command fills, from its parameters, without meshing anything.
pub type Targets = fn(&Session, &Value) -> Result<Vec<String>>;

/// Every command that fills bodies with a mesh and replaces theirs, with the bodies it fills.
pub static MESHERS: &[(&str, Targets)] = &[
    ("mesh.tet.generate", body_or_all),
    ("mesh.hex.generate", body_or_all),
    ("mesh.hex.sweep", body_or_all),
    ("mesh.hex.mapped", body_or_all),
    ("mesh.hex.octree", body_or_all),
    ("mesh.hex.blocks", listed_or_all),
    ("mesh.hex.project", body_or_all),
    ("mesh.cells", listed_or_all),
    ("mesh.plan.apply", super::plan::plan_bodies),
    ("mesh.planar", super::mesh::planar_body),
    ("mesh.hex.extrude", super::extrude::extrude_body),
    ("mesh.mirror", super::mirror::image_names),
];

/// `body`, the imported `surface`, or every body of the open design.
fn body_or_all(s: &Session, p: &Value) -> Result<Vec<String>> {
    match (p.get("body"), p.get("surface")) {
        (Some(b), None) => Ok(vec![b.as_str().ok_or_else(|| EngineError::param("body", crate::schema::ProblemKind::NotAName))?.to_string()]),
        (None, Some(_)) => Ok(vec![super::hex::surface(s, p)?.0]),
        (None, None) => design_bodies(s),
        _ => Err(EngineError::Params("give `body` or `surface`, not both".into())),
    }
}

/// `bodies`, or every body of the open design.
fn listed_or_all(s: &Session, p: &Value) -> Result<Vec<String>> {
    match p.get("bodies") {
        Some(_) => super::req(p, "bodies"),
        None => design_bodies(s),
    }
}

fn design_bodies(s: &Session) -> Result<Vec<String>> {
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    Ok(design.bodies.iter().map(|b| b.name.clone()).collect())
}

/// Tetrahedra and hexahedra of either order: the kinds an order change swaps.
fn family(k: ElementKind) -> Option<u8> {
    match k {
        ElementKind::Tet4 | ElementKind::Tet10 => Some(0),
        ElementKind::Hex8 | ElementKind::Hex20 => Some(1),
        _ => None,
    }
}

/// The solid element type `command` makes with `params`, when it makes one type.
fn made_kind(command: &str, params: &Value) -> Option<ElementKind> {
    match command {
        "mesh.tet.generate" => {
            Some(if params.get("order").and_then(Value::as_str) == Some("quadratic") { ElementKind::Tet10 } else { ElementKind::Tet4 })
        }
        "mesh.hex.project" => Some(ElementKind::Hex8),
        "mesh.hex.generate" | "mesh.hex.sweep" | "mesh.hex.mapped" | "mesh.hex.blocks" | "mesh.plan.apply" | "mesh.hex.extrude" => {
            Some(if params.get("order").and_then(Value::as_i64) == Some(2) { ElementKind::Hex20 } else { ElementKind::Hex8 })
        }
        _ => None,
    }
}

/// What running a meshing command would replace, for a confirmation before it runs.
#[derive(Debug, Clone, PartialEq)]
pub struct Replacement {
    /// "CrankArm already has a mesh (7,274 hex20). Replace it?"
    pub question: String,
    /// What else goes with the old mesh: sets that list its labels, solved results.
    pub warnings: Vec<String>,
}

/// What `command` with `params` would replace: `None` when it is not a mesher, its parameters
/// do not name bodies, or none of its bodies has a mesh.
pub fn replacement(s: &Session, command: &str, params: &Value) -> Option<Replacement> {
    let targets = MESHERS.iter().find(|(id, _)| *id == command)?.1;
    let bodies = targets(s, params).ok()?;
    let mut meshed = Vec::new();
    let mut removed = BodyMesh::default();
    for body in bodies.iter().collect::<BTreeSet<_>>() {
        let mesh = s.model.body_mesh(body);
        if mesh.elements.is_empty() {
            continue;
        }
        let mut counts: BTreeMap<ElementKind, usize> = BTreeMap::new();
        for e in mesh.elements.iter().filter_map(|l| s.model.element(*l)) {
            *counts.entry(e.kind).or_default() += 1;
        }
        let kinds: Vec<String> = ElementKind::ALL.iter().filter_map(|k| counts.get(k).map(|n| format!("{} {}", thousands(*n), k.name()))).collect();
        meshed.push((body.as_str(), kinds.join(", ")));
        removed.elements.extend(mesh.elements);
        removed.nodes.extend(mesh.nodes);
    }
    let question = match meshed.as_slice() {
        [] => return None,
        [(body, counts)] => format!("{body} already has a mesh ({counts}). Replace it?"),
        many => {
            let list: Vec<String> = many.iter().map(|(body, counts)| format!("{body} ({counts})")).collect();
            format!("{} already have meshes. Replace them?", list.join(", "))
        }
    };
    let mut warnings = Vec::new();
    // The element type changing (tet10 → tet4 when the order is left out) is easy to miss.
    if let Some(new) = made_kind(command, params) {
        let old: BTreeSet<ElementKind> = removed.elements.iter().filter_map(|l| s.model.element(*l)).map(|e| e.kind).collect();
        let changed: Vec<&str> = old.iter().filter(|k| **k != new && family(**k) == family(new)).map(|k| k.name()).collect();
        if !changed.is_empty() {
            warnings.push(format!("The element type changes: {} → {}.", changed.join(", "), new.name()));
        }
    }
    let stale = s.model.sets_using(&removed);
    if !stale.is_empty() {
        warnings.push(format!("Sets that list its nodes or elements by label ({}) will name labels of the old mesh.", stale.join(", ")));
    }
    let results: Vec<&str> = s.solved_steps().collect();
    if !results.is_empty() {
        warnings.push(format!("The results of {} will be cleared.", results.join(", ")));
    }
    Some(Replacement { question, warnings })
}

/// `n` with thousands separators: 7274 → "7,274".
pub fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use fieldcraft_model::Element;

    use super::*;

    #[test]
    fn every_mesher_is_a_command() {
        for (id, _) in MESHERS {
            assert!(crate::find_command(id).is_some(), "{id}");
        }
        assert_eq!(thousands(7274), "7,274");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_000_000), "1,000,000");
    }

    /// A 20 × 10 × 10 box split at x = 10, opened as `Box.1` (x < 10) and `Box.2`.
    fn two_boxes() -> (crate::tests::TempDir, Session) {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 20, "width": 10, "height": 10, "body_name": "Box"}},
            {"command": "solid.split_body", "params": {"body": "Box", "plane": {"origin": [10, 0, 0], "normal": [1, 0, 0]}}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    /// Remeshing without the order would turn tet10 into tet4: the question says so (#343).
    #[test]
    fn an_element_type_change_is_named() {
        let (_dir, mut s) = two_boxes();
        s.execute("mesh.tet.generate", &json!({"body": "Box.1", "size": 5, "order": "quadratic"})).unwrap();
        let r = replacement(&s, "mesh.tet.generate", &json!({"body": "Box.1", "size": 4})).unwrap();
        assert!(r.warnings.iter().any(|w| w == "The element type changes: tet10 → tet4."), "{r:?}");
        let same = replacement(&s, "mesh.tet.generate", &json!({"body": "Box.1", "order": "quadratic"})).unwrap();
        assert!(same.warnings.iter().all(|w| !w.contains("element type")), "{same:?}");
    }

    /// A body's elements, by label.
    fn elements_of(s: &Session, body: &str) -> BTreeMap<u32, Element> {
        s.model.body_elements(body).into_iter().filter_map(|l| Some((l, s.model.element(l)?.clone()))).collect()
    }

    /// Every mesher in [`MESHERS`], run twice on a body: the second run replaces the first (as
    /// many elements as a fresh mesh, reported as `replaced`), leaves the other body's mesh as it
    /// was, names the set that listed the old labels, and one undo brings the first mesh back.
    #[test]
    fn every_mesher_replaces_the_mesh_of_its_body_only() {
        let tet_other = json!({"command": "mesh.tet.generate", "params": {"body": "Box.2", "size": 5}});
        let hex_box = json!({"command": "mesh.hex.generate", "params": {"body": "Box.1", "size": 2.5}});
        let square = json!([[0, 0], [10, 0], [10, 10], [0, 10]]);
        // Box.1's bottom and top faces, for the extrusion between them.
        let (_probe_dir, probe) = two_boxes();
        let face = |nz: f64| {
            let b = probe.design.as_deref().and_then(|d| d.body("Box.1")).unwrap();
            let surfaces = b.face_surfaces().unwrap();
            let i =
                surfaces.iter().position(|f| matches!(f, Some(fieldcraft_geometry::FaceSurface::Plane { normal, .. }) if normal[2] == nz)).unwrap();
            json!({"body": "Box.1", "face": b.faces()[i].name})
        };
        let extrude = json!({"source": [face(-1.0)], "to_surface": [face(1.0)], "direction": [0, 0, 1], "layers": 4, "size": 2.5});
        let cases: Vec<(&str, Vec<Value>, Value, &str, &str)> = vec![
            ("mesh.tet.generate", vec![tet_other.clone()], json!({"body": "Box.1", "size": 4}), "Box.1", "Box.2"),
            ("mesh.hex.generate", vec![tet_other.clone()], json!({"body": "Box.1", "size": 2.5}), "Box.1", "Box.2"),
            ("mesh.hex.sweep", vec![tet_other.clone()], json!({"body": "Box.1", "size": 2.5}), "Box.1", "Box.2"),
            ("mesh.hex.mapped", vec![tet_other.clone()], json!({"body": "Box.1", "size": 2.5}), "Box.1", "Box.2"),
            ("mesh.hex.octree", vec![tet_other.clone()], json!({"body": "Box.1", "size": 2.5}), "Box.1", "Box.2"),
            ("mesh.hex.blocks", vec![tet_other.clone()], json!({"bodies": ["Box.1"], "size": 2.5}), "Box.1", "Box.2"),
            ("mesh.hex.project", vec![tet_other.clone()], json!({"body": "Box.1", "size": 2.5, "merge": 0}), "Box.1", "Box.2"),
            ("mesh.cells", vec![tet_other.clone()], json!({"bodies": ["Box.1"], "size": 2.5, "default": "hex"}), "Box.1", "Box.2"),
            (
                "mesh.plan.apply",
                vec![tet_other.clone(), json!({"command": "mesh.plan", "params": {"body": "Box.1", "size": 2.5}})],
                json!({"plan": "Plan-Box_1"}),
                "Box.1",
                "Box.2",
            ),
            ("mesh.hex.extrude", vec![tet_other.clone()], extrude, "Box.1", "Box.2"),
            ("mesh.planar", vec![tet_other.clone()], json!({"polygon": square, "size": 2.5, "name": "Box.1"}), "Box.1", "Box.2"),
            // The image's tets (meshed before) are replaced by the mirrored hexahedra.
            (
                "mesh.mirror",
                vec![hex_box, tet_other],
                json!({"plane": {"point": [10, 0, 0], "normal": [1, 0, 0]}, "bodies": ["Box.1"]}),
                "Box.2",
                "Box.1",
            ),
        ];
        let covered: BTreeSet<&str> = cases.iter().map(|c| c.0).collect();
        assert_eq!(covered, MESHERS.iter().map(|m| m.0).collect::<BTreeSet<_>>(), "every mesher has a case");

        let (_dir, fresh) = two_boxes();
        for (id, prep, params, body, other) in cases {
            let mut s = fresh.clone();
            for c in &prep {
                s.execute(c["command"].as_str().unwrap(), &c["params"]).unwrap();
            }
            // What the session says it would replace matches what the command replaces.
            let first = s.execute(id, &params).unwrap();
            let fresh_count = elements_of(&s, body).len();
            assert!(fresh_count > 0, "{id}: {first}");
            let others = elements_of(&s, other);
            assert!(!others.is_empty(), "{id}");
            let some: Vec<u32> = elements_of(&s, body).keys().take(3).copied().collect();
            s.execute("model.set.create", &json!({"name": "by_label", "elements": some})).unwrap();
            let asked = replacement(&s, id, &params).unwrap_or_else(|| panic!("{id}: no question"));
            assert!(asked.question.starts_with(&format!("{body} already has a mesh (")), "{id}: {asked:?}");
            assert!(asked.warnings.iter().any(|w| w.contains("by_label")), "{id}: {asked:?}");

            let before = s.model.clone();
            let second = s.execute(id, &params).unwrap();
            assert_eq!(elements_of(&s, body).len(), fresh_count, "{id}: {second}");
            assert_eq!(elements_of(&s, other), others, "{id}: the other body's mesh is untouched");
            assert_eq!(
                second["replaced"],
                json!([{"body": body, "elements": fresh_count, "nodes": second["replaced"][0]["nodes"]}]),
                "{id}: {second}"
            );
            assert_eq!(second["stale_sets"], json!(["by_label"]), "{id}: {second}");
            // No node is left that no element uses.
            let used: BTreeSet<u32> = s.model.elements().values().flat_map(|e| e.nodes.iter().copied()).collect();
            assert!(s.model.nodes().keys().all(|n| used.contains(n)), "{id}: orphan nodes");

            // A mesher is deterministic: the same mesh again is the same model, and no undo step.
            if s.model != before {
                s.execute("edit.undo", &json!({})).unwrap();
                assert_eq!(s.model, before, "{id}: undo restores the first mesh");
            }
        }
    }

    #[test]
    fn tets_over_hexahedra_leave_only_tets() {
        let (_dir, mut s) = two_boxes();
        s.execute("mesh.hex.generate", &json!({"size": 2.5})).unwrap();
        let hexes = elements_of(&s, "Box.1").len();
        let other = elements_of(&s, "Box.2");
        let r = s.execute("mesh.tet.generate", &json!({"body": "Box.1", "size": 4})).unwrap();
        assert_eq!(r["replaced"][0]["elements"], hexes, "{r}");
        assert!(elements_of(&s, "Box.1").values().all(|e| e.kind == ElementKind::Tet4));
        assert_eq!(elements_of(&s, "Box.2"), other);
        let tets = s.model.clone();
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(elements_of(&s, "Box.1").len(), hexes, "one undo brings the hexahedra back");
        assert!(elements_of(&s, "Box.1").values().all(|e| e.kind == ElementKind::Hex8));
        s.execute("edit.redo", &json!({})).unwrap();
        assert_eq!(s.model, tets);
        // Meshing every body replaces both.
        let r = s.execute("mesh.tet.generate", &json!({"size": 4})).unwrap();
        assert_eq!(r["replaced"].as_array().map(Vec::len), Some(2), "{r}");
        assert!(s.model.elements().values().all(|e| e.kind == ElementKind::Tet4));
        // A body without a mesh asks nothing.
        let (_dir2, s2) = two_boxes();
        assert_eq!(replacement(&s2, "mesh.tet.generate", &json!({"body": "Box.1"})), None);
        assert_eq!(replacement(&s, "mesh.block", &json!({"body": "Box.1"})), None, "not a mesher");
    }

    #[test]
    fn results_solved_on_the_old_mesh_are_cleared_and_undo_brings_them_back() {
        let mut s = Session::default();
        let corners = json!([[0, 0, 0], [10, 0, 0], [10, 1, 0], [0, 1, 0], [0, 0, 1], [10, 0, 1], [10, 1, 1], [0, 1, 1]]);
        s.run_script(&json!([
            {"command": "mesh.block", "params": {"body": "Bar", "kind": "hex8", "corners": corners, "divisions": [10, 1, 1], "faces": {"6": "root", "4": "tip"}}},
            {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Bar"}]}},
            {"command": "model.set.create", "params": {"name": "root", "geometry": [{"body": "Bar", "face": "root"}]}},
            {"command": "model.set.create", "params": {"name": "tip", "geometry": [{"body": "Bar", "face": "tip"}]}},
            {"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000.0, "poisson": 0.0}}},
            {"command": "model.section.create", "params": {"name": "all", "set": "all", "material": "Mat", "kind": {"type": "solid"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "Clamp", "type": "fixed", "target": {"set": "root"}}},
            {"command": "model.load.create", "params": {"name": "Pull", "type": "force", "target": {"set": "tip"}, "vector": [1, 0, 0]}},
            {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["Clamp"], "loads": ["Pull"]}}
        ]))
        .unwrap();
        s.execute("solve.static", &json!({"step": "Load"})).unwrap();
        assert!(s.result("Load").is_some());
        let asked = replacement(&s, "mesh.planar", &json!({"polygon": [[0, 0], [10, 0], [10, 1], [0, 1]], "size": 1, "name": "Bar"})).unwrap();
        assert_eq!(asked.question, "Bar already has a mesh (10 hex8). Replace it?");
        assert!(asked.warnings.iter().any(|w| w.contains("results of Load")), "{asked:?}");
        let r = s.execute("mesh.planar", &json!({"polygon": [[0, 0], [10, 0], [10, 1], [0, 1]], "size": 1, "name": "Bar"})).unwrap();
        assert_eq!(r["results_cleared"], json!(["Load"]), "{r}");
        assert!(s.result("Load").is_none());
        assert!(s.model.elements().values().all(|e| e.kind.dimension() == 2));
        s.execute("edit.undo", &json!({})).unwrap();
        assert!(s.result("Load").is_some(), "undo brings the solved state back");
    }
}
