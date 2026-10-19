//! `mesh.hex.extrude`: hexahedra from a quad section up to a target surface
//! ([`fieldcraft_mesh::hex::extrude_to_surface`]).
//!
//! The section is, in this order: the quadrilaterals of a set (`source_set`: shell elements or
//! element faces); the mesh already on the source faces (element faces tagged with them, or the
//! boundary faces of a neighbouring body's mesh lying on them), whose nodes the extrusion shares
//! so it conforms with that mesh; or a fresh quad mesh of a planar source face. Target and side
//! faces are projected onto through the trimmed B-rep projection that projection blocking uses
//! ([`fieldcraft_mesh::tet::geometry::Tessellation::project_trimmed`]): exact on analytic faces,
//! the tessellation's accuracy on free-form ones.

use std::collections::BTreeMap;

use fieldcraft_geometry::{FaceSurface, Shape};
use fieldcraft_mesh::hex::{ExtrudeLayers, ExtrudePath, ExtrudeRole, ExtrudeSpec, MAX_BIAS, MAX_LAYERS, extrude_to_surface};
use fieldcraft_mesh::planar::{face_boundary, mesh_region};
use fieldcraft_mesh::project::Geometry;
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use fieldcraft_mesh::tet::geometry::Tessellation;
use fieldcraft_model::{ElementKind, Entity, GeometryRef, Model};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::schema::ProblemKind;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const HEX_EXTRUDE: CommandSpec = CommandSpec {
    id: "mesh.hex.extrude",
    label: "Extrude to Surface",
    params: r#"{"source": [{"body": "Post", "face": "F1"}] (the section: faces already meshed, by this body or a neighbour, or one planar face meshed here) or "source_set": "Section" (a set of quad shell elements or element faces, such as a 2D mesh), "to_surface": [{"body": "Wall", "face": "F7"}] (the target faces the last layer lies on), "sides": [{"body": "Post", "face": "F3"}, …] (optional: the faces the section's sides follow), "direction": [0, 0, 1] or "normal": true (each node along the section's normal toward the target; the default), "layers": 8 or "size": 2.0 (layer thickness, a length; also the size of a section meshed here; default 1/25 of the paths' extent), "bias": 1.0 (last layer's thickness over the first's, 0.001 to 1000), "body": "Post" (the body the elements fill; default: the first source face's body), "tolerance": 0.05 (chordal tolerance of the tessellation, a length; default 1e-3 of each body's diagonal)} — replaces `body`'s mesh with hex8 layers stacked on the section, each node along its own straight path, the last layer on the target surface; with `sides`, every layer's boundary nodes lie on the side faces (which need not be straight) and the interior follows them. A section taken from an existing mesh keeps its nodes, so the extrusion conforms with the neighbour across it. Boundary element faces are tagged with `body`'s source, target and side faces. Reports the layers, counts, quality, the target gap and where the section came from; a path that misses the target, or elements that would invert, is an error that says where"#,
    form: Some(super::forms::HEX_EXTRUDE),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: hex_extrude,
};

/// Most faces in `source`, `to_surface` or `sides`.
const MAX_FACES: usize = 1000;
/// Default layer size and chordal tolerance, relative to the extent.
const DEFAULT_RELATIVE_SIZE: f64 = 1.0 / 25.0;
const DEFAULT_RELATIVE_TOLERANCE: f64 = 1e-3;
const MIN_RELATIVE_TOLERANCE: f64 = 1e-5;
/// Most quadrilaterals of a section meshed here.
const MAX_SECTION_QUADS: usize = 1_000_000;

/// The bodies the extrusion's faces lie on, and every face it names, with one projection.
struct Faces {
    bodies: Vec<Tessellation>,
    /// Per body: its name, and its faces' names.
    names: Vec<(String, Vec<String>)>,
    /// Face id → (body, B-rep face).
    faces: Vec<(usize, u32)>,
    /// The largest chordal tolerance used.
    tolerance: f64,
}

impl Geometry for Faces {
    fn project(&self, face: u32, p: [f64; 3]) -> Option<[f64; 3]> {
        let &(b, f) = self.faces.get(face as usize)?;
        self.bodies.get(b)?.project_trimmed(f, p)
    }
    fn normal(&self, face: u32, p: [f64; 3]) -> Option<[f64; 3]> {
        let &(b, f) = self.faces.get(face as usize)?;
        self.bodies.get(b)?.face_normal(f, p)
    }
}

impl Faces {
    /// The ids of the faces `refs` names (a name selects all its pieces), loading their bodies.
    fn add(&mut self, s: &Session, refs: &[GeometryRef], tolerance: Option<f64>, key: &str) -> Result<Vec<u32>> {
        let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
        let mut out = Vec::new();
        for r in refs {
            if !matches!(r.entity, Entity::Face(_)) {
                return Err(EngineError::invalid(key, format!("{r} is not a face")));
            }
            let shape = design.body(&r.body).ok_or_else(|| EngineError::invalid(key, format!("no body `{}`", r.body)))?;
            let b = match self.names.iter().position(|(n, _)| *n == shape.name) {
                Some(b) => b,
                None => self.load(shape, tolerance)?,
            };
            let names = self.names.get(b).map(|(_, f)| f.clone()).unwrap_or_default();
            let before = out.len();
            for (i, name) in names.iter().enumerate() {
                if r.selects(&GeometryRef::face(shape.name.as_str(), name.as_str())?) {
                    let id = self.id(b, i as u32);
                    if !out.contains(&id) {
                        out.push(id);
                    }
                }
            }
            if out.len() == before {
                return Err(EngineError::invalid(key, format!("body `{}` has no face `{}`", r.body, r.entity.name().unwrap_or_default())));
            }
        }
        Ok(out)
    }

    fn load(&mut self, shape: &Shape, tolerance: Option<f64>) -> Result<usize> {
        let diag = shape.size()?;
        let tol = tolerance.unwrap_or(DEFAULT_RELATIVE_TOLERANCE * diag);
        if tol < MIN_RELATIVE_TOLERANCE * diag {
            return Err(EngineError::Params(format!(
                "`tolerance` must be at least {MIN_RELATIVE_TOLERANCE:e} of body `{}`'s bounding-box diagonal",
                shape.name
            )));
        }
        let mesh = shape.tessellate(tol)?;
        let surface = BoundarySurface::new(
            &mesh.positions,
            &mesh.triangles,
            &mesh.tri_face,
            &SurfaceOptions { weld_distance: mesh.weld_distance, ..SurfaceOptions::default() },
        )
        .map_err(|e| EngineError::Params(format!("body `{}`: {e}", shape.name)))?;
        self.bodies.push(Tessellation::new(&surface, super::mesh::analytic_surfaces(shape)));
        self.names.push((shape.name.clone(), shape.faces().iter().map(|f| f.name.clone()).collect()));
        self.tolerance = self.tolerance.max(tol);
        Ok(self.bodies.len() - 1)
    }

    fn id(&mut self, b: usize, f: u32) -> u32 {
        match self.faces.iter().position(|x| *x == (b, f)) {
            Some(i) => i as u32,
            None => {
                self.faces.push((b, f));
                (self.faces.len() - 1) as u32
            }
        }
    }

    /// The body and face name of face `id`.
    fn name(&self, id: u32) -> Option<(&str, &str)> {
        let &(b, f) = self.faces.get(id as usize)?;
        let (body, faces) = self.names.get(b)?;
        Some((body.as_str(), faces.get(f as usize)?.as_str()))
    }
}

/// The section: node coordinates, the model label of each node it shares (none for a fresh
/// section), and its quadrilaterals.
struct Section {
    nodes: Vec<[f64; 3]>,
    labels: Vec<Option<u32>>,
    quads: Vec<[u32; 4]>,
    from: &'static str,
}

impl Section {
    /// A section of model nodes: each label once, in order of first use.
    fn of_labels(model: &Model, quads: &[[u32; 4]], from: &'static str) -> Result<Section> {
        let mut index: BTreeMap<u32, u32> = BTreeMap::new();
        let mut s = Section { nodes: Vec::new(), labels: Vec::new(), quads: Vec::with_capacity(quads.len()), from };
        for q in quads {
            let mut local = [0u32; 4];
            for (slot, &l) in local.iter_mut().zip(q) {
                *slot = match index.get(&l) {
                    Some(&i) => i,
                    None => {
                        let xyz = model.nodes().get(&l).map(|n| n.xyz).ok_or_else(|| EngineError::Params(format!("node {l} is not in the model")))?;
                        let i = s.nodes.len() as u32;
                        s.nodes.push(xyz);
                        s.labels.push(Some(l));
                        index.insert(l, i);
                        i
                    }
                };
            }
            s.quads.push(local);
        }
        Ok(s)
    }
}

/// The four corners of element `e`'s face `face` (1-based), when it is a quadrilateral.
fn quad_face(model: &Model, e: u32, face: u8) -> Option<[u32; 4]> {
    let el = model.elements().get(&e)?;
    let local = el.kind.faces().get(usize::from(face).checked_sub(1)?)?;
    if local.len() < 4 || (local.len() != 4 && local.len() != 8) {
        return None;
    }
    let corner = |k: usize| local.get(k).and_then(|&i| el.nodes.get(i)).copied();
    Some([corner(0)?, corner(1)?, corner(2)?, corner(3)?])
}

/// The corners of a quad shell element.
fn quad_shell(model: &Model, e: u32) -> Option<[u32; 4]> {
    let el = model.elements().get(&e)?;
    if !matches!(el.kind, ElementKind::Quad4 | ElementKind::Quad8) {
        return None;
    }
    Some([*el.nodes.first()?, *el.nodes.get(1)?, *el.nodes.get(2)?, *el.nodes.get(3)?])
}

/// The elements of `body`'s mesh.
fn body_elements(model: &Model, body: &str) -> std::collections::BTreeSet<u32> {
    model.body_elements(body).into_iter().collect()
}

/// The quadrilaterals tagged with the faces `refs` names (shell elements filling them, element
/// faces on them), leaving out `exclude`'s elements.
fn tagged_quads(model: &Model, refs: &[GeometryRef], exclude: &std::collections::BTreeSet<u32>) -> Result<Vec<[u32; 4]>> {
    let c = model.classification();
    let mut tags: Vec<u32> = refs.iter().flat_map(|r| c.tags_selected_by(r)).collect();
    tags.sort_unstable();
    tags.dedup();
    let mut out = Vec::new();
    for (&e, t) in c.elements() {
        if exclude.contains(&e) {
            continue;
        }
        if t.region.is_some_and(|r| tags.binary_search(&r).is_ok()) {
            out.push(quad_shell(model, e).ok_or_else(|| EngineError::Params(format!("element {e} on the source face is not a quadrilateral")))?);
        }
        for (&f, tag) in &t.faces {
            if tags.binary_search(tag).is_ok() {
                out.push(
                    quad_face(model, e, f)
                        .ok_or_else(|| EngineError::Params(format!("element {e}'s face {f} on the source face is not a quadrilateral")))?,
                );
            }
        }
    }
    Ok(out)
}

/// The boundary faces of other bodies' meshes (tagged faces) that lie on the source faces.
fn neighbour_quads(model: &Model, geometry: &Faces, source: &[u32], exclude: &std::collections::BTreeSet<u32>) -> Vec<[u32; 4]> {
    let tol = 10.0 * geometry.tolerance.max(1e-9);
    let on = |p: [f64; 3]| {
        source.iter().any(|&f| geometry.project(f, p).is_some_and(|q| (0..3).map(|k| (q[k] - p[k]).powi(2)).sum::<f64>().sqrt() <= tol))
    };
    let xyz = |l: u32| model.nodes().get(&l).map(|n| n.xyz);
    let mut out = Vec::new();
    for (&e, t) in model.classification().elements() {
        if exclude.contains(&e) {
            continue;
        }
        for &f in t.faces.keys() {
            let Some(q) = quad_face(model, e, f) else { continue };
            let Some(pts) = q.iter().map(|&l| xyz(l)).collect::<Option<Vec<_>>>() else { continue };
            let centre = pts.iter().fold([0.0; 3], |a, p| [a[0] + p[0] / 4.0, a[1] + p[1] / 4.0, a[2] + p[2] / 4.0]);
            if pts.iter().all(|p| on(*p)) && on(centre) {
                out.push(q);
            }
        }
    }
    out
}

/// A quad mesh of one planar face of `shape` of about `size`.
fn mesh_planar_face(shape: &Shape, face: u32, size: f64, tolerance: f64) -> Result<Section> {
    let surfaces = shape.face_surfaces()?;
    let Some(Some(FaceSurface::Plane { normal, point })) = surfaces.get(face as usize).copied() else {
        return Err(EngineError::Params(format!(
            "the source face of body `{}` is not meshed and not planar: mesh it first (a 2D mesh, or a neighbouring body's mesh)",
            shape.name
        )));
    };
    let n = unit(normal).ok_or_else(|| EngineError::Params("the source face has no normal".into()))?;
    let t = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let u = unit(cross(n, t)).ok_or_else(|| EngineError::Params("the source face has no normal".into()))?;
    let v = cross(n, u);
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let mesh = shape.tessellate(tolerance)?;
    // The tessellation in the face's frame, the face at z = 0.
    let local: Vec<[f64; 3]> = mesh
        .positions
        .iter()
        .map(|p| {
            let d = [p[0] - point[0], p[1] - point[1], p[2] - point[2]];
            [dot(d, u), dot(d, v), dot(d, n)]
        })
        .collect();
    let boundary = face_boundary(&local, &mesh.triangles, &mesh.tri_face, face, mesh.weld_distance.max(tolerance))?;
    let opts = fieldcraft_mesh::quad::QuadOptions { max_elements: MAX_SECTION_QUADS, ..Default::default() };
    let m = mesh_region(&boundary.loops, ElementKind::Quad4, size, &opts)?;
    let at = |xy: [f64; 2]| [0, 1, 2].map(|k| point[k] + xy[0] * u[k] + xy[1] * v[k]);
    let nodes: Vec<[f64; 3]> = m.nodes.iter().map(|xy| at(*xy)).collect();
    let quads = m
        .elements
        .iter()
        .map(|e| match e.as_slice() {
            [a, b, c, d, ..] => Ok([*a, *b, *c, *d]),
            _ => Err(EngineError::Params("the section mesher made an element that is not a quadrilateral".into())),
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Section { labels: vec![None; nodes.len()], nodes, quads, from: "meshed" })
}

fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    (l > 0.0 && l.is_finite()).then(|| a.map(|x| x / l))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// A list of face references.
fn face_refs(p: &Value, key: &str, required: bool) -> Result<Vec<GeometryRef>> {
    let refs: Vec<GeometryRef> = match p.get(key) {
        None if !required => return Ok(Vec::new()),
        None => return Err(EngineError::param(key, ProblemKind::Required)),
        Some(v) => serde_json::from_value(v.clone()).map_err(|e| EngineError::invalid(key, format!("a list of {{\"body\", \"face\"}}: {e}")))?,
    };
    if refs.len() > MAX_FACES || (required && refs.is_empty()) {
        return Err(EngineError::invalid(key, format!("1 to {MAX_FACES} faces")));
    }
    Ok(refs)
}

/// The body the elements fill: `body`, or the first source face's.
pub(super) fn extrude_body(_: &Session, p: &Value) -> Result<Vec<String>> {
    if let Some(b) = p.get("body") {
        let name = b.as_str().ok_or_else(|| EngineError::param("body", ProblemKind::NotAName))?;
        GeometryRef::body(name)?;
        return Ok(vec![name.to_string()]);
    }
    let source = face_refs(p, "source", false)?;
    let first = source.first().ok_or_else(|| EngineError::Params("give `body` (the body the elements fill) with `source_set`".into()))?;
    Ok(vec![first.body.clone()])
}

fn hex_extrude(s: &mut Session, p: &Value) -> Result<Value> {
    let tolerance = super::hex::number(p, "tolerance")?;
    let size = super::hex::number(p, "size")?;
    let layers = match p.get("layers") {
        None => None,
        Some(v) => Some(
            v.as_u64()
                .filter(|n| (1..=MAX_LAYERS as u64).contains(n))
                .ok_or_else(|| EngineError::invalid("layers", format!("a whole number from 1 to {MAX_LAYERS}")))? as usize,
        ),
    };
    let bias = match p.get("bias") {
        None => 1.0,
        Some(v) => v
            .as_f64()
            .filter(|b| (1.0 / MAX_BIAS..=MAX_BIAS).contains(b))
            .ok_or_else(|| EngineError::invalid("bias", format!("a number from {} to {MAX_BIAS}", 1.0 / MAX_BIAS)))?,
    };
    let path = match (p.get("direction"), p.get("normal")) {
        (Some(d), None) => {
            let d: [f64; 3] = serde_json::from_value(d.clone()).map_err(|_| EngineError::param("direction", ProblemKind::NotAVector { len: 3 }))?;
            if !d.iter().all(|x| x.is_finite()) || unit(d).is_none() {
                return Err(EngineError::invalid("direction", "a vector with a length"));
            }
            ExtrudePath::Direction(d)
        }
        (None, None) => ExtrudePath::Normal,
        (None, Some(n)) => match n.as_bool() {
            Some(true) => ExtrudePath::Normal,
            _ => return Err(EngineError::invalid("normal", "true, or give `direction`")),
        },
        (Some(_), Some(_)) => return Err(EngineError::Params("give `direction` or `normal`, not both".into())),
    };
    let source = face_refs(p, "source", false)?;
    let source_set = match p.get("source_set") {
        None => None,
        Some(v) => Some(v.as_str().ok_or_else(|| EngineError::param("source_set", ProblemKind::NotAName))?.to_string()),
    };
    if source.is_empty() == source_set.is_none() {
        return Err(EngineError::Params("give the section as `source` (faces) or `source_set` (a set), one of them".into()));
    }
    let targets = face_refs(p, "to_surface", true)?;
    let sides = face_refs(p, "sides", false)?;
    let body = extrude_body(s, p)?.into_iter().next().unwrap_or_default();
    let body_ref = GeometryRef::body(body.as_str())?;
    let mut faces = Faces { bodies: Vec::new(), names: Vec::new(), faces: Vec::new(), tolerance: 0.0 };
    let source_ids = faces.add(s, &source, tolerance, "source")?;
    let target_ids = faces.add(s, &targets, tolerance, "to_surface")?;
    let side_ids = faces.add(s, &sides, tolerance, "sides")?;
    if let Some(f) = target_ids.iter().find(|f| source_ids.contains(f) || side_ids.contains(f)) {
        let (b, n) = faces.name(*f).unwrap_or_default();
        return Err(EngineError::Params(format!("face `{n}` of body `{b}` is both a target and a source or side face")));
    }

    // The section: a set, the mesh already on the source faces, or a fresh mesh of one face.
    let own = body_elements(&s.model, &body);
    let section = match &source_set {
        Some(name) => {
            let r = s.model.resolve_set(name)?;
            let mut quads = Vec::new();
            for &e in &r.elements {
                if let Some(q) = quad_shell(&s.model, e) {
                    quads.push(q);
                }
            }
            for f in &r.faces {
                quads.push(
                    quad_face(&s.model, f.element, f.face).ok_or_else(|| {
                        EngineError::Params(format!("set `{name}`: element {}'s face {} is not a quadrilateral", f.element, f.face))
                    })?,
                );
            }
            if quads.is_empty() {
                return Err(EngineError::Params(format!("set `{name}` has no quadrilaterals (quad shell elements or hex element faces)")));
            }
            Section::of_labels(&s.model, &quads, "set")?
        }
        None => {
            let tagged = tagged_quads(&s.model, &source, &own)?;
            if !tagged.is_empty() {
                Section::of_labels(&s.model, &tagged, "mesh")?
            } else {
                let near = neighbour_quads(&s.model, &faces, &source_ids, &own);
                if !near.is_empty() {
                    Section::of_labels(&s.model, &near, "neighbour")?
                } else {
                    let [id] = source_ids.as_slice() else {
                        return Err(EngineError::Params("the source faces are not meshed: give one face to mesh here, or mesh them first".into()));
                    };
                    let (b, f) = faces.faces.get(*id as usize).copied().ok_or_else(|| EngineError::Params("no source face".into()))?;
                    let name = faces.names.get(b).map(|n| n.0.clone()).unwrap_or_default();
                    let shape = s.design.as_deref().and_then(|d| d.body(&name)).ok_or_else(|| EngineError::Params(format!("no body `{name}`")))?;
                    let extent = shape.size()?;
                    let h = size.unwrap_or(DEFAULT_RELATIVE_SIZE * extent);
                    mesh_planar_face(shape, f, h, tolerance.unwrap_or(DEFAULT_RELATIVE_TOLERANCE * extent))?
                }
            }
        }
    };

    let layers = match (layers, size) {
        (Some(n), _) => ExtrudeLayers::Count(n),
        (None, Some(h)) => ExtrudeLayers::Size(h),
        (None, None) => {
            let extent = s.design.as_deref().and_then(|d| d.body(&body)).map(Shape::size).transpose()?.unwrap_or(1.0);
            ExtrudeLayers::Size(DEFAULT_RELATIVE_SIZE * extent)
        }
    };
    let spec = ExtrudeSpec {
        nodes: &section.nodes,
        quads: &section.quads,
        path,
        targets: &target_ids,
        sides: &side_ids,
        layers,
        bias,
        tolerance: faces.tolerance.max(1e-9),
        max_elements: fieldcraft_model::MAX_ELEMENTS / 8,
    };
    let m = extrude_to_surface(&spec, &faces).map_err(|e| EngineError::of_body(&body, e))?;

    // Replace the body's mesh; the section's shared nodes stay.
    let mut model = s.model.clone();
    let cleared = super::replace::clear(&mut model, [body.as_str()]);
    if let Some(l) = section.labels.iter().flatten().find(|l| !model.nodes().contains_key(l)) {
        return Err(EngineError::Params(format!("the section's node {l} belongs to body `{body}`'s mesh, which the extrusion replaces")));
    }
    let n_src = section.nodes.len();
    let fresh = m.nodes.len() - section.labels.iter().flatten().count();
    let (first_node, first_element) = super::fresh_labels(&model, fresh, m.hexes.len())?;
    let mut labels = Vec::with_capacity(m.nodes.len());
    let mut next = first_node;
    for (i, xyz) in m.nodes.iter().enumerate() {
        let shared = if i < n_src { section.labels.get(i).copied().flatten() } else { None };
        labels.push(match shared {
            Some(l) => l,
            None => {
                model.set_node(next, *xyz)?;
                next += 1;
                next - 1
            }
        });
    }
    let body_tag = model.geometry_tag(&body_ref)?;
    for (e, h) in m.hexes.iter().enumerate() {
        let nodes: Vec<u32> = h.iter().map(|&i| labels.get(i as usize).copied().unwrap_or(0)).collect();
        let l = first_element + e as u32;
        model.set_element(l, ElementKind::Hex8, nodes)?;
        model.tag_element(l, body_tag)?;
    }
    // Boundary faces on the body's own B-rep faces get their tags.
    let source_own: Vec<u32> = source_ids.iter().copied().filter(|&f| faces.name(f).is_some_and(|(b, _)| b == body)).collect();
    let mut tags: BTreeMap<u32, u32> = BTreeMap::new();
    for b in &m.boundary {
        let face = match b.role {
            ExtrudeRole::Source => source_own.first().copied(),
            ExtrudeRole::Target(f) | ExtrudeRole::Side(Some(f)) => Some(f),
            ExtrudeRole::Side(None) => None,
        };
        let Some((owner, name)) = face.and_then(|f| faces.name(f)) else { continue };
        if owner != body {
            continue;
        }
        let key = face.unwrap_or_default();
        let tag = match tags.get(&key) {
            Some(&t) => t,
            None => {
                let t = model.geometry_tag(&GeometryRef::face(body.as_str(), name)?)?;
                tags.insert(key, t);
                t
            }
        };
        model.tag_element_face(first_element + b.element, b.face, tag)?;
    }
    let (nodes, elements) = (model.nodes().len(), model.elements().len());
    s.model = model;
    let reply = json!({
        "body": body,
        "section": {"from": section.from, "quadrilaterals": section.quads.len(), "nodes": n_src, "shared_nodes": section.labels.iter().flatten().count()},
        "layers": m.layers,
        "kind": "hex8",
        "elements": m.hexes.len(),
        "new_nodes": fresh,
        "quality": {"min_scaled_jacobian": m.min_scaled_jacobian, "mean_scaled_jacobian": m.mean_scaled_jacobian},
        "volume": m.gate.volume,
        "max_target_gap": m.max_target_gap,
        "nodes": nodes,
        "model_elements": elements,
    });
    Ok(cleared.report(s, reply))
}

#[cfg(test)]
mod tests;
