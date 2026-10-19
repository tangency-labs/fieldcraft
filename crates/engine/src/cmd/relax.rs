//! Node relaxation (#398): `mesh.node.relax` moves the nodes inside a sphere of influence (or
//! the nodes on listed CAD faces, or listed nodes or elements) to fix the worst elements there, nodes on the geometry
//! staying on it ([`fieldcraft_mesh::relax`]), and keeps what it did as a relaxation field of the
//! model; `mesh.node.relax.apply` does a field again on the current mesh (after a remesh), and a
//! field marked `reapply` is applied by itself after every remesh of its bodies.
//!
//! **Where a node may go** comes from the mesh's geometry tags: a node on one B-rep face slides
//! in that face's parameters, a node where two faces meet slides along their edge, a node on a
//! geometry vertex or on three or more faces stays. Faces of the same surface (the two halves a
//! kernel splits a cylinder into, or one face as cut into cells) count once. Partition cut faces
//! (`partition:<name>:a` and `:b`) are not geometry: a node on one moves freely, and the cells on
//! both sides share it, so the mesh stays conforming across the cut. A node on the skin whose
//! geometry is unknown (an imported mesh) stays.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_mesh::relax::{self, Constraints, Falloff, Freedom, RegionQuality, Relaxed};
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use fieldcraft_mesh::tet::geometry::{AnalyticSurface, Tessellation};
use fieldcraft_model::geometry::{Entity, GeometryRef};
use fieldcraft_model::{Model, Relaxation};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::schema::ProblemKind;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const RELAX: CommandSpec = CommandSpec {
    id: "mesh.node.relax",
    label: "Relax Nodes",
    params: r#"{"center": [x, y, z], "radius": mm (the sphere of influence) | "around_element": label (a sphere round that element; radius default twice its longest edge) | "faces": [{"body": "Base", "face": "F3"}] (persistent B-rep face names, as model.set.create takes them: every node on those faces and the nodes of the elements touching them; pick faces in the view) | "nodes": [labels] | "elements": [labels] (their nodes), "bodies"?: ["Body1"] (only their nodes; a partitioned body's cells count as it), "iterations"?: 5 (1..=100), "strength"?: 1 (0..1], "falloff"?: "smooth" (default: weight 1 at the centre to 0 at the radius) | "linear" | "none", "field"?: "Name" (the relaxation field to keep; default Relax1, Relax2, …), "reapply"?: false (apply the field again after each remesh of these bodies), "mode"?: "guarded" (default) | "strong", "below"?: 0.7 (strong: the elements below this scaled Jacobian are the bad ones), "reach"?: mm (strong: nodes this close to a bad node follow it; default five times the bad elements' mean edge length), "guards"?: "valid" (default) | "strict" (strong: when a cluster's morph is kept)} — moves the nodes to improve the worst elements around them (smart Laplacian, then optimization of the worst element; a move is kept only when the worst element around the node does not get worse, and never inverts one). `"mode": "strong"` morphs the region instead (`bodies` alone is the whole body, nothing listed the whole mesh; `iterations` are rounds, default 20, 1..=200): each round moves the corner nodes of the elements below `below` to better places (their own optimisation, the rest held) and carries that displacement to the nodes within `reach` by harmonic (graph Laplacian) interpolation, zero beyond; nodes on B-rep edges and vertices are held, nodes on faces put back on them, nodes outside the reach unchanged. A cluster's morph is kept whenever no element it touches inverts (`"guards": "valid"`); with `"guards": "strict"` it is kept (at full or a halved scale) only when besides no element gets longer than the longest was, and the worst element it touches improves. The result is the round with the highest smallest scaled Jacobian, never lower than the start. Its reply adds the mean scaled Jacobian and, per round, the bad elements and nodes, their clusters (bad nodes within the reach of each other, each morphed and checked on its own) and how many were kept, the region's nodes, the reach and the min and mean (`rounds`). Nodes on a B-rep face slide in its UV parameters and stay on it, nodes on a B-rep edge slide along it, nodes on vertices stay; nodes on partition cut faces move freely (shared by both cells); mid-side nodes follow, on the geometry. Reports the nodes moved, the largest displacement, the region's smallest scaled Jacobian, largest edge ratio and worst elements before and after, the field kept, and the listed faces no mesh is tagged with (`faces_without_mesh`; an error when none is). One undo step"#,
    form: Some(super::forms::NODE_RELAX),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: relax_nodes,
};

pub const APPLY: CommandSpec = CommandSpec {
    id: "mesh.node.relax.apply",
    label: "Apply Relaxation",
    params: r#"{"field": "Name", "bodies"?: ["Body1"] (default the field's), "polish"?: 0.3 (strength of one guarded smoothing sweep after; 0 for none)} — applies a relaxation field kept by mesh.node.relax to the current mesh (typically after a remesh): each node inside the relaxed region takes the displacement interpolated from the field at its parametric coordinates in the old cells, on its face or edge as mesh.node.relax moves it; a move that would invert an element is undone; then one guarded smoothing sweep. Reports the nodes moved and the largest displacement"#,
    form: Some(super::forms::NODE_RELAX_APPLY),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: apply,
};

/// The commands after which `reapply` fields are applied again.
const REMESH: &[&str] = &[
    "mesh.tet.generate",
    "mesh.hex.generate",
    "mesh.hex.sweep",
    "mesh.hex.mapped",
    "mesh.hex.blocks",
    "mesh.plan.apply",
    "mesh.hex.octree",
    "mesh.cells",
    "mesh.block",
];

const DEFAULT_ITERATIONS: usize = 5;

/// The region's quality as JSON.
pub(crate) fn quality_json(q: &RegionQuality) -> Value {
    json!({
        "elements": q.elements,
        "min_scaled_jacobian": q.min_quality,
        "max_edge_ratio": if q.max_edge_ratio.is_finite() { json!(q.max_edge_ratio) } else { Value::Null },
        "invalid": q.invalid,
        "worst": q.worst.iter().map(|(e, v)| json!({"element": e, "quality": v})).collect::<Vec<_>>(),
    })
}

/// Faces a node lies on (geometry references), from the element faces, shells and vertex tags.
struct NodeTags {
    faces: BTreeMap<u32, BTreeSet<(String, String)>>,
    vertex: BTreeSet<u32>,
    /// Nodes on a partition cut face.
    cut: BTreeSet<u32>,
}

fn node_tags(model: &Model, nodes: &BTreeSet<u32>) -> NodeTags {
    let c = model.classification();
    let mut t = NodeTags { faces: BTreeMap::new(), vertex: BTreeSet::new(), cut: BTreeSet::new() };
    let note = |n: u32, tag: u32, t: &mut NodeTags| {
        let Some(g) = c.entity(tag) else { return };
        if let Entity::Face(name) = &g.entity {
            if name.starts_with("partition:") {
                t.cut.insert(n);
            } else {
                t.faces.entry(n).or_default().insert((g.body.clone(), name.clone()));
            }
        }
    };
    for (&label, tags) in c.elements() {
        let Some(e) = model.element(label) else { continue };
        if !e.nodes.iter().any(|n| nodes.contains(n)) {
            continue;
        }
        if e.kind.dimension() == 2
            && let Some(r) = tags.region
        {
            for &n in &e.nodes {
                note(n, r, &mut t);
            }
        }
        for (&face, &tag) in &tags.faces {
            let Some(local) = e.kind.faces().get(usize::from(face).wrapping_sub(1)) else { continue };
            for &i in local.iter() {
                if let Some(&n) = e.nodes.get(i) {
                    note(n, tag, &mut t);
                }
            }
        }
    }
    for (&n, &tag) in c.nodes() {
        if nodes.contains(&n) && c.entity(tag).is_some_and(|g| matches!(g.entity, Entity::Vertex(_))) {
            t.vertex.insert(n);
        }
    }
    t
}

/// Whether two exact surfaces are the same surface.
fn same_surface(a: &AnalyticSurface, b: &AnalyticSurface, tol: f64) -> bool {
    use AnalyticSurface as S;
    let sub = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let len = |x: [f64; 3]| dot(x, x).sqrt();
    let unit = |x: [f64; 3]| {
        let l = len(x);
        if l > 0.0 { [x[0] / l, x[1] / l, x[2] / l] } else { x }
    };
    let parallel = |x: [f64; 3], y: [f64; 3]| (dot(unit(x), unit(y)).abs() - 1.0).abs() < 1e-9;
    let on_axis = |p: [f64; 3], o: [f64; 3], a: [f64; 3]| {
        let d = sub(p, o);
        let a = unit(a);
        let along = dot(d, a);
        len([d[0] - a[0] * along, d[1] - a[1] * along, d[2] - a[2] * along]) <= tol
    };
    match (*a, *b) {
        (S::Plane { normal: n1, point: p1 }, S::Plane { normal: n2, point: p2 }) => parallel(n1, n2) && dot(sub(p2, p1), unit(n1)).abs() <= tol,
        (S::Cylinder { origin: o1, axis: a1, radius: r1 }, S::Cylinder { origin: o2, axis: a2, radius: r2 }) => {
            (r1 - r2).abs() <= tol && parallel(a1, a2) && on_axis(o2, o1, a1)
        }
        (S::Sphere { center: c1, radius: r1 }, S::Sphere { center: c2, radius: r2 }) => len(sub(c1, c2)) <= tol && (r1 - r2).abs() <= tol,
        (S::Cone { apex: p1, axis: a1, half_angle: h1 }, S::Cone { apex: p2, axis: a2, half_angle: h2 }) => {
            len(sub(p1, p2)) <= tol && dot(unit(a1), unit(a2)) > 1.0 - 1e-9 && (h1 - h2).abs() <= 1e-9
        }
        (S::Torus { center: c1, axis: a1, major: m1, minor: n1 }, S::Torus { center: c2, axis: a2, major: m2, minor: n2 }) => {
            len(sub(c1, c2)) <= tol && parallel(a1, a2) && (m1 - m2).abs() <= tol && (n1 - n2).abs() <= tol
        }
        _ => false,
    }
}

/// The design's geometry as relaxation constraints.
pub(crate) struct DesignGeometry {
    /// Per body: its tessellation with exact surfaces.
    bodies: Vec<Tessellation>,
    /// Face ids: (body, B-rep face).
    faces: Vec<(usize, u32)>,
    /// Curve ids: (body, the two B-rep faces).
    curves: Vec<(usize, [u32; 2])>,
    freedom: BTreeMap<u32, Freedom>,
}

impl Constraints for DesignGeometry {
    fn freedom(&self, node: u32) -> Option<Freedom> {
        self.freedom.get(&node).copied()
    }
    fn slide_face(&self, face: u32, p: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
        let &(b, f) = self.faces.get(face as usize)?;
        self.bodies.get(b)?.slide_face(f, p, d)
    }
    fn slide_curve(&self, curve: u32, p: [f64; 3], d: [f64; 3]) -> Option<[f64; 3]> {
        let &(b, f) = self.curves.get(curve as usize)?;
        self.bodies.get(b)?.slide_curve(f, p, d)
    }
}

/// The constraints of `nodes` from the mesh's tags and the open design.
pub(crate) fn design_geometry(s: &Session, nodes: &BTreeSet<u32>) -> Result<DesignGeometry> {
    let tags = node_tags(&s.model, nodes);
    let mut g = DesignGeometry { bodies: Vec::new(), faces: Vec::new(), curves: Vec::new(), freedom: BTreeMap::new() };
    // Body name → (index into `g.bodies`, its exact surfaces), built when first needed.
    let mut loaded: BTreeMap<String, Option<(usize, Vec<Option<AnalyticSurface>>)>> = BTreeMap::new();
    let mut load = |name: &str, g: &mut DesignGeometry| -> Result<Option<(usize, Vec<Option<AnalyticSurface>>)>> {
        if let Some(l) = loaded.get(name) {
            return Ok(l.clone());
        }
        let shape = s.design.as_deref().and_then(|d| d.body(name));
        let entry = match shape {
            None => None,
            Some(shape) => {
                let tol = (shape.size()? * 1e-3).max(1e-6);
                let (mesh, pieces) = shape.tessellate_pieces(tol)?;
                let surface = BoundarySurface::new(
                    &mesh.positions,
                    &mesh.triangles,
                    &mesh.tri_face,
                    &SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() },
                )?;
                let analytic = super::mesh::analytic_surfaces(shape);
                g.bodies.push(Tessellation::new(&surface, analytic.clone()).with_pieces(&surface, super::mesh::pieces(shape, &pieces)));
                Some((g.bodies.len() - 1, analytic))
            }
        };
        loaded.insert(name.to_string(), entry.clone());
        Ok(entry)
    };
    let tol = 1e-6 * s.model.nodes().values().fold(0.0f64, |m, n| m.max(n.xyz.iter().fold(0.0f64, |a, x| a.max(x.abs())))).max(1.0);
    let mut face_ids: BTreeMap<(usize, u32), u32> = BTreeMap::new();
    let mut curve_ids: BTreeMap<(usize, [u32; 2]), u32> = BTreeMap::new();
    for &n in nodes {
        if tags.vertex.contains(&n) {
            g.freedom.insert(n, Freedom::Fixed);
            continue;
        }
        let Some(refs) = tags.faces.get(&n) else {
            if tags.cut.contains(&n) {
                g.freedom.insert(n, Freedom::Free);
            }
            continue;
        };
        // Resolve each face; group those on one surface.
        let mut groups: Vec<Vec<(usize, u32, Option<AnalyticSurface>)>> = Vec::new();
        let mut unknown = false;
        for (body, face) in refs {
            let Some((b, analytic)) = load(body, &mut g)? else {
                unknown = true;
                continue;
            };
            let Some(f) = s.design.as_deref().and_then(|d| d.body(body)).and_then(|sh| sh.face_index(face)) else {
                unknown = true;
                continue;
            };
            let surface = analytic.get(f as usize).copied().flatten();
            let group = groups.iter_mut().find(|grp| {
                grp.iter().any(|(gb, gf, gs)| (*gb == b && *gf == f) || matches!((gs, surface), (Some(x), Some(y)) if same_surface(x, &y, tol)))
            });
            match group {
                Some(grp) => grp.push((b, f, surface)),
                None => groups.push(vec![(b, f, surface)]),
            }
        }
        let freedom = match groups.as_slice() {
            _ if unknown => Freedom::Fixed,
            [one] => {
                let Some(&(b, f, _)) = one.first() else { continue };
                let next = g.faces.len() as u32;
                let id = *face_ids.entry((b, f)).or_insert(next);
                if id == next {
                    g.faces.push((b, f));
                }
                Freedom::Face(id)
            }
            [x, y] => {
                // A body with a face of each surface carries their edge.
                let pair = x.iter().find_map(|&(bx, fx, _)| y.iter().find(|(by, _, _)| *by == bx).map(|&(_, fy, _)| (bx, [fx.min(fy), fx.max(fy)])));
                match pair {
                    Some(key) => {
                        let next = g.curves.len() as u32;
                        let id = *curve_ids.entry(key).or_insert(next);
                        if id == next {
                            g.curves.push(key);
                        }
                        Freedom::Curve(id)
                    }
                    None => Freedom::Fixed,
                }
            }
            _ => Freedom::Fixed,
        };
        g.freedom.insert(n, freedom);
    }
    Ok(g)
}

/// The nodes of `bodies`' elements, when given.
fn body_nodes(model: &Model, p: &Value) -> Result<Option<(Vec<String>, BTreeSet<u32>)>> {
    let Some(bodies) = super::opt::<Vec<String>>(p, "bodies")? else { return Ok(None) };
    let mut nodes = BTreeSet::new();
    // A body cut into cells by partitions is meshed as its cells (`Base/a`, `Base/b/a`, …).
    let meshed: BTreeSet<String> = model.classification().entities().iter().map(|r| r.body.clone()).collect();
    for b in &bodies {
        let cell = format!("{b}/");
        let elements: Vec<u32> = meshed.iter().filter(|m| *m == b || m.starts_with(&cell)).flat_map(|m| model.body_elements(m)).collect();
        if elements.is_empty() {
            return Err(EngineError::Params(format!("body `{}` has no mesh", b.chars().take(80).collect::<String>())));
        }
        for e in elements {
            if let Some(el) = model.element(e) {
                nodes.extend(el.nodes.iter().copied());
            }
        }
    }
    Ok(Some((bodies, nodes)))
}

pub(crate) fn point(p: &Value, key: &str) -> Result<Option<[f64; 3]>> {
    match super::opt::<[f64; 3]>(p, key)? {
        Some(c) if c.iter().all(|x| x.is_finite()) => Ok(Some(c)),
        Some(_) => Err(EngineError::param(key, ProblemKind::NotANumber)),
        None => Ok(None),
    }
}

/// The most faces `faces` may list.
const MAX_FACES: usize = 10_000;

/// The `faces` parameter: a list of `{body, face}` B-rep face references.
pub(crate) fn face_refs(p: &Value) -> Result<Vec<GeometryRef>> {
    let refs: Vec<GeometryRef> = super::req(p, "faces")?;
    if refs.is_empty() || refs.len() > MAX_FACES {
        return Err(EngineError::invalid("faces", format!("1 to {MAX_FACES} faces")));
    }
    if let Some(r) = refs.iter().find(|r| !matches!(r.entity, Entity::Face(_))) {
        return Err(EngineError::invalid("faces", format!("`{}`: give {{\"body\", \"face\"}}", r.body.chars().take(80).collect::<String>())));
    }
    Ok(refs)
}

/// The tags of the faces `refs` select, and the faces with no tagged mesh.
fn face_tags<'a>(model: &Model, refs: &'a [GeometryRef]) -> (BTreeSet<u32>, Vec<&'a GeometryRef>) {
    let c = model.classification();
    let mut tags = BTreeSet::new();
    let mut missing = Vec::new();
    for r in refs {
        let found = c.tags_selected_by(r);
        if found.is_empty() {
            missing.push(r);
        }
        tags.extend(found);
    }
    (tags, missing)
}

/// The listed faces with no tagged mesh, as `{body, face}` (the reply names them).
pub(crate) fn unmeshed_faces(model: &Model, p: &Value) -> Vec<Value> {
    match p.get("faces").map(|_| face_refs(p)) {
        Some(Ok(refs)) => face_tags(model, &refs).1.into_iter().map(|r| json!(r)).collect(),
        _ => Vec::new(),
    }
}

/// The nodes on the B-rep faces `refs` (from the mesh's geometry tags: element faces on them and
/// shell elements meshing them). Faces with no tagged mesh are left out (the reply names them);
/// none with one is an error.
pub(crate) fn nodes_on_faces(model: &Model, refs: &[GeometryRef]) -> Result<BTreeSet<u32>> {
    let c = model.classification();
    let (tags, missing) = face_tags(model, refs);
    if tags.is_empty() {
        let names: Vec<String> =
            missing.iter().take(5).map(|r| format!("{}/{}", r.body, r.entity.name().unwrap_or_default()).chars().take(120).collect()).collect();
        return Err(EngineError::Params(format!("no mesh is tagged with the faces ({}): mesh their bodies, or mesh.tag", names.join(", "))));
    }
    let mut on = BTreeSet::new();
    for (&label, t) in c.elements() {
        let Some(e) = model.element(label) else { continue };
        if e.kind.dimension() == 2 && t.region.is_some_and(|r| tags.contains(&r)) {
            on.extend(e.nodes.iter().copied());
        }
        for (&face, tag) in &t.faces {
            if !tags.contains(tag) {
                continue;
            }
            let Some(local) = e.kind.faces().get(usize::from(face).wrapping_sub(1)) else { continue };
            on.extend(local.iter().filter_map(|&i| e.nodes.get(i).copied()));
        }
    }
    Ok(on)
}

/// The nodes on the B-rep faces `refs` ([`nodes_on_faces`]) and the nodes of every element
/// touching those.
fn face_nodes(model: &Model, refs: &[GeometryRef]) -> Result<BTreeSet<u32>> {
    let on = nodes_on_faces(model, refs)?;
    let mut out = on.clone();
    for e in model.elements().values() {
        if e.nodes.iter().any(|n| on.contains(n)) {
            out.extend(e.nodes.iter().copied());
        }
    }
    Ok(out)
}

/// The sphere around element `label`: its centroid, and twice its longest corner edge.
fn around_element(model: &Model, label: u32) -> Result<([f64; 3], f64)> {
    let e = model.element(label).ok_or_else(|| EngineError::Params(format!("no element {label}")))?;
    let x: Vec<[f64; 3]> = e.nodes.iter().take(e.kind.corner_count()).filter_map(|n| model.node(*n).map(|p| p.xyz)).collect();
    if x.is_empty() {
        return Err(EngineError::Params(format!("element {label} has no nodes")));
    }
    let c = x.iter().fold([0.0; 3], |s, p| [s[0] + p[0], s[1] + p[1], s[2] + p[2]]).map(|v| v / x.len() as f64);
    let mut longest = 0.0f64;
    for edge in e.kind.edges() {
        let at = |i: usize| edge.get(i).and_then(|&k| x.get(k)).copied();
        if let (Some(a), Some(b)) = (at(0), at(1)) {
            longest = longest.max(((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt());
        }
    }
    Ok((c, 2.0 * longest))
}

/// The nodes to move and their weights, and the sphere when there is one.
pub(crate) fn weights(model: &Model, p: &Value) -> Result<(BTreeMap<u32, f64>, Option<([f64; 3], f64)>)> {
    let strength = super::optional_number(p, "strength", 1.0)?;
    if !(strength > 0.0 && strength <= 1.0) {
        return Err(EngineError::param("strength", ProblemKind::NotInRange { min: 0.0, max: 1.0 }));
    }
    let falloff = match super::opt::<String>(p, "falloff")? {
        None => Falloff::Smooth,
        Some(f) => Falloff::from_name(&f).ok_or_else(|| EngineError::param("falloff", ProblemKind::NotAChoice))?,
    };
    let only = body_nodes(model, p)?.map(|(_, n)| n);
    let around = match super::opt::<u32>(p, "around_element")? {
        Some(e) => Some(around_element(model, e)?),
        None => None,
    };
    let center = point(p, "center")?.or(around.map(|a| a.0));
    let radius = match p.get("radius") {
        Some(_) => Some(super::optional_number(p, "radius", 0.0)?)
            .filter(|r| *r > 0.0 && r.is_finite())
            .ok_or_else(|| EngineError::param("radius", ProblemKind::NotPositive))
            .map(Some)?,
        None => around.map(|a| a.1),
    };
    let listed = ["faces", "nodes", "elements"].iter().filter(|k| p.get(**k).is_some()).count();
    if let (Some(c), Some(r)) = (center, radius) {
        if listed > 0 {
            return Err(EngineError::Params(
                "give a sphere (`center` and `radius`, or `around_element`) or `faces`, `nodes` or `elements`, not both".into(),
            ));
        }
        return Ok((relax::sphere_weights(model, c, r, strength, falloff, only.as_ref())?, Some((c, r))));
    }
    if center.is_some() != radius.is_some() {
        return Err(EngineError::Params("a sphere of influence needs both `center` and `radius`".into()));
    }
    let mut nodes: BTreeSet<u32> = BTreeSet::new();
    if p.get("nodes").is_some() {
        nodes.extend(super::labels(p, "nodes", fieldcraft_model::MAX_NODES)?);
    }
    let mut elements: Vec<u32> = Vec::new();
    if p.get("elements").is_some() {
        elements.extend(super::labels(p, "elements", fieldcraft_model::MAX_ELEMENTS)?);
    }
    if p.get("faces").is_some() {
        nodes.extend(face_nodes(model, &face_refs(p)?)?);
    }
    for e in elements {
        let el = model.element(e).ok_or_else(|| EngineError::Params(format!("no element {e}")))?;
        nodes.extend(el.nodes.iter().copied());
    }
    if listed == 0 {
        return Err(EngineError::Params(
            "say which nodes to relax: give one of `center` and `radius`, `around_element`, `faces`, `nodes` or `elements`".into(),
        ));
    }
    if let Some(only) = &only {
        nodes.retain(|n| only.contains(n));
    }
    Ok((relax::uniform_weights(model, &nodes, strength)?, None))
}

fn iterations(p: &Value) -> Result<usize> {
    match p.get("iterations") {
        None => Ok(DEFAULT_ITERATIONS),
        Some(v) => v
            .as_u64()
            .and_then(|i| usize::try_from(i).ok())
            .filter(|i| (1..=relax::MAX_ITERATIONS).contains(i))
            .ok_or_else(|| EngineError::param("iterations", ProblemKind::NotInRange { min: 1.0, max: relax::MAX_ITERATIONS as f64 })),
    }
}

/// The nodes the relaxation may touch: the weighted ones and their elements' nodes.
fn touched(model: &Model, weights: &BTreeMap<u32, f64>) -> BTreeSet<u32> {
    let mut out: BTreeSet<u32> = weights.keys().copied().collect();
    for e in model.elements().values() {
        if e.nodes.iter().any(|n| weights.contains_key(n)) {
            out.extend(e.nodes.iter().copied());
        }
    }
    out
}

/// What `mesh.node.relax` with `params` would do, without doing it: the open form's preview.
pub struct Preview {
    pub sphere: Option<([f64; 3], f64)>,
    pub relaxed: Relaxed,
}

/// The relaxation `params` (as typed, in the display system) ask for on the session's mesh,
/// without doing it: positions and the sphere in SI, as the session is stored.
pub fn preview(s: &Session, params: &Value) -> Option<Preview> {
    // The strong mode morphs a whole region: too heavy for a live preview.
    if params.get("mode").and_then(Value::as_str) == Some("strong") {
        return None;
    }
    // The meshers work in millimetres: the preview is made there and brought back.
    let framed = crate::units::geometry_session(s).ok()?;
    let params = crate::units::geometry_params(s, RELAX.id, params).ok()?;
    let (w, sphere) = weights(&framed.model, &params).ok()?;
    let g = design_geometry(&framed, &touched(&framed.model, &w)).ok()?;
    let mut relaxed = relax::relax(&framed.model, &w, iterations(&params).ok()?, &g).ok()?;
    for p in relaxed.positions.values_mut() {
        *p = crate::units::geometry_point(*p);
    }
    relaxed.max_displacement = crate::units::geometry_length(relaxed.max_displacement);
    Some(Preview { sphere: sphere.map(|(c, r)| (crate::units::geometry_point(c), crate::units::geometry_length(r))), relaxed })
}

/// The sphere `params` (as typed) describe, in SI, for the form's gizmo before the rest is valid.
pub fn sphere(s: &Session, params: &Value) -> Option<([f64; 3], f64)> {
    let length = |v: &Value| crate::units::read(v, fieldcraft_model::units::Quantity::Length, s.input_units())?.ok();
    let around = params.get("around_element").and_then(Value::as_u64).and_then(|e| around_element(&s.model, u32::try_from(e).ok()?).ok());
    let centre = params.get("center").and_then(Value::as_array).and_then(|a| match a.as_slice() {
        [x, y, z] => Some([length(x)?, length(y)?, length(z)?]),
        _ => None,
    });
    let c = centre.or(around.map(|a| a.0))?;
    let r = params.get("radius").and_then(length).filter(|r| *r > 0.0 && r.is_finite()).or(around.map(|a| a.1))?;
    Some((c, r))
}

/// Default rounds of the strong mode.
const STRONG_ROUNDS: usize = 20;

/// `mesh.node.relax` with `"mode": "strong"`: morph the region round its bad elements
/// ([`relax::morph`]).
fn relax_strong(s: &mut Session, p: &Value) -> Result<Value> {
    let only = body_nodes(&s.model, p)?;
    let below = super::optional_number(p, "below", relax::MorphOptions::default().below)?;
    if !(below > 0.0 && below <= 1.0) {
        return Err(EngineError::param("below", ProblemKind::NotInRange { min: 0.0, max: 1.0 }));
    }
    let reach = match p.get("reach") {
        None => None,
        Some(_) => Some(super::optional_number(p, "reach", 0.0)?)
            .filter(|r| *r > 0.0 && r.is_finite())
            .ok_or_else(|| EngineError::param("reach", ProblemKind::NotPositive))
            .map(Some)?,
    };
    let guards = match super::opt::<String>(p, "guards")?.as_deref() {
        None | Some("valid") => relax::Guards::Valid,
        Some("strict") => relax::Guards::Strict,
        Some(_) => return Err(EngineError::param("guards", ProblemKind::NotAChoice)),
    };
    let rounds = match p.get("iterations") {
        None => STRONG_ROUNDS,
        Some(v) => v
            .as_u64()
            .and_then(|i| usize::try_from(i).ok())
            .filter(|i| (1..=relax::MAX_ROUNDS).contains(i))
            .ok_or_else(|| EngineError::param("iterations", ProblemKind::NotInRange { min: 1.0, max: relax::MAX_ROUNDS as f64 }))?,
    };
    let listed = ["center", "radius", "around_element", "faces", "nodes", "elements"].iter().any(|k| p.get(*k).is_some());
    let nodes: BTreeSet<u32> = if listed {
        weights(&s.model, p)?.0.into_keys().collect()
    } else if let Some((_, n)) = &only {
        n.clone()
    } else {
        s.model.nodes().keys().copied().collect()
    };
    let nodes: BTreeSet<u32> = relax::uniform_weights(&s.model, &nodes, 1.0)?.into_keys().collect();
    if nodes.is_empty() {
        return Err(EngineError::Params("no mesh nodes in the region to relax".into()));
    }
    let g = design_geometry(s, &touched(&s.model, &nodes.iter().map(|&n| (n, 1.0)).collect()))?;
    let r = relax::morph(&s.model, &nodes, relax::MorphOptions { below, reach, rounds, guards }, &g)?;
    let field_name = match super::opt::<String>(p, "field")? {
        Some(n) => n,
        None => (1..).map(|i| format!("Relax{i}")).find(|n| s.model.relaxation(n).is_none()).unwrap_or_else(|| "Relax".into()),
    };
    let bodies = only.map(|(b, _)| b).unwrap_or_default();
    let reapply = super::opt::<bool>(p, "reapply")?.unwrap_or(false);
    if !r.positions.is_empty() {
        let cells = relax::field_cells(&s.model, &Relaxed { positions: r.positions.clone(), ..Relaxed::default() });
        s.model.set_relaxation(&field_name, Relaxation { bodies, cells, reapply })?;
    }
    for (n, x) in &r.positions {
        s.model.set_node(*n, *x)?;
    }
    let summary = format!(
        "strong: moved {} node{} (largest {:.3e}) in {} round{}; scaled Jacobian min {:.3} → {:.3}, mean {:.3} → {:.3}{}",
        r.moved,
        if r.moved == 1 { "" } else { "s" },
        r.max_displacement,
        r.rounds.len(),
        if r.rounds.len() == 1 { "" } else { "s" },
        r.before.min_quality,
        r.after.min_quality,
        r.mean_before,
        r.mean_after,
        if r.positions.is_empty() {
            format!("; nothing below {below} could be improved, nothing moved")
        } else {
            format!("; kept as field {field_name}")
        }
    );
    let rounds: Vec<Value> = r
        .rounds
        .iter()
        .map(|k| {
            json!({"bad_elements": k.bad_elements, "bad_nodes": k.bad_nodes, "region_nodes": k.region_nodes, "reach": k.reach,
                   "clusters": k.clusters, "kept": k.kept, "min": k.min_quality, "mean": k.mean_quality})
        })
        .collect();
    let mut before = quality_json(&r.before);
    let mut after = quality_json(&r.after);
    if let (Some(b), Some(a)) = (before.as_object_mut(), after.as_object_mut()) {
        b.insert("mean_scaled_jacobian".into(), json!(r.mean_before));
        a.insert("mean_scaled_jacobian".into(), json!(r.mean_after));
    }
    Ok(json!({
        "mode": "strong",
        "moved": r.moved,
        "nodes_changed": r.positions.len(),
        "max_displacement": r.max_displacement,
        "iterations": r.rounds.len(),
        "rounds": rounds,
        "before": before,
        "after": after,
        "field": (!r.positions.is_empty()).then_some(field_name),
        "faces_without_mesh": unmeshed_faces(&s.model, p),
        "summary": summary,
    }))
}

fn relax_nodes(s: &mut Session, p: &Value) -> Result<Value> {
    match super::opt::<String>(p, "mode")?.as_deref() {
        None | Some("guarded") => {}
        Some("strong") => return relax_strong(s, p),
        Some(_) => return Err(EngineError::param("mode", ProblemKind::NotAChoice)),
    }
    if ["below", "reach", "guards"].iter().any(|k| p.get(*k).is_some()) {
        return Err(EngineError::Params("`below`, `reach` and `guards` go with \"mode\": \"strong\"".into()));
    }
    let (w, sphere) = weights(&s.model, p)?;
    if w.is_empty() {
        return Err(EngineError::Params("no mesh nodes in the region to relax".into()));
    }
    let iterations = iterations(p)?;
    let g = design_geometry(s, &touched(&s.model, &w))?;
    let r = relax::relax(&s.model, &w, iterations, &g)?;
    let field_name = match super::opt::<String>(p, "field")? {
        Some(n) => n,
        None => (1..).map(|i| format!("Relax{i}")).find(|n| s.model.relaxation(n).is_none()).unwrap_or_else(|| "Relax".into()),
    };
    let bodies = match body_nodes(&s.model, p)? {
        Some((b, _)) => b,
        None => Vec::new(),
    };
    let reapply = super::opt::<bool>(p, "reapply")?.unwrap_or(false);
    if !r.positions.is_empty() {
        let cells = relax::field_cells(&s.model, &r);
        s.model.set_relaxation(&field_name, Relaxation { bodies, cells, reapply })?;
    }
    for (n, x) in &r.positions {
        s.model.set_node(*n, *x)?;
    }
    let summary = format!(
        "moved {} node{} (largest {:.3e}); scaled Jacobian min from {:.3} to {:.3}{}",
        r.moved,
        if r.moved == 1 { "" } else { "s" },
        r.max_displacement,
        r.before.min_quality,
        r.after.min_quality,
        if r.positions.is_empty() { String::new() } else { format!("; kept as field {field_name}") }
    );
    Ok(json!({
        "moved": r.moved,
        "nodes_changed": r.positions.len(),
        "max_displacement": r.max_displacement,
        "iterations": r.iterations,
        "sphere": sphere.map(|(c, r)| json!({"center": c, "radius": r})),
        "before": quality_json(&r.before),
        "after": quality_json(&r.after),
        "field": (!r.positions.is_empty()).then_some(field_name),
        "faces_without_mesh": unmeshed_faces(&s.model, p),
        "summary": summary,
    }))
}

/// Apply field `name` to the nodes of `bodies` (or of the whole mesh).
fn apply_field(s: &mut Session, name: &str, bodies: &[String], polish: f64) -> Result<Relaxed> {
    let field = s
        .model
        .relaxation(name)
        .cloned()
        .ok_or_else(|| EngineError::Params(format!("no relaxation field `{}`", name.chars().take(80).collect::<String>())))?;
    let mut nodes = BTreeSet::new();
    let bodies = if bodies.is_empty() { &field.bodies[..] } else { bodies };
    if bodies.is_empty() {
        nodes.extend(s.model.nodes().keys().copied());
    } else {
        for b in bodies {
            for e in s.model.body_elements(b) {
                if let Some(el) = s.model.element(e) {
                    nodes.extend(el.nodes.iter().copied());
                }
            }
        }
    }
    // Only nodes the field reaches need constraints.
    nodes.retain(|n| s.model.node(*n).is_some_and(|p| relax::field_displacement(&field.cells, p.xyz).is_some()));
    let g = design_geometry(s, &nodes)?;
    let r = relax::apply_field(&s.model, &field.cells, Some(&nodes), &g, polish)?;
    for (n, x) in &r.positions {
        s.model.set_node(*n, *x)?;
    }
    Ok(r)
}

fn apply(s: &mut Session, p: &Value) -> Result<Value> {
    let name: String = super::req(p, "field")?;
    let bodies = super::opt::<Vec<String>>(p, "bodies")?.unwrap_or_default();
    let polish = super::optional_number(p, "polish", 0.3)?;
    if !(0.0..=1.0).contains(&polish) {
        return Err(EngineError::param("polish", ProblemKind::NotInRange { min: 0.0, max: 1.0 }));
    }
    let r = apply_field(s, &name, &bodies, polish)?;
    Ok(json!({
        "moved": r.moved,
        "nodes_changed": r.positions.len(),
        "max_displacement": r.max_displacement,
        "summary": format!("applied field {name}: moved {} node{} (largest {:.3e})", r.moved, if r.moved == 1 { "" } else { "s" }, r.max_displacement),
    }))
}

/// After a remesh (`id`), apply again every field marked `reapply`; its reply lists them.
pub(crate) fn after_remesh(s: &mut Session, id: &str, mut reply: Value) -> Result<Value> {
    if !REMESH.contains(&id) {
        return Ok(reply);
    }
    let fields: Vec<String> = s.model.relaxations().iter().filter(|(_, f)| f.reapply).map(|(n, _)| n.clone()).collect();
    let mut applied = Vec::new();
    for name in fields {
        let r = apply_field(s, &name, &[], 0.3)?;
        applied.push(json!({"field": name, "moved": r.moved}));
    }
    if !applied.is_empty()
        && let Some(o) = reply.as_object_mut()
    {
        o.insert("relaxations_applied".into(), Value::Array(applied));
    }
    Ok(reply)
}

#[cfg(test)]
mod tests;
