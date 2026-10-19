use crate::schema::ProblemKind;
use std::collections::BTreeSet;
use std::sync::Arc;

use fieldcraft_geometry::tagged_mesh;
use fieldcraft_geometry::{Design, MAX_DESIGN_BYTES, Shape, SurfaceMesh, open_design};
use fieldcraft_io::SourceFile;
use fieldcraft_io::roots::Access;
use fieldcraft_model::Model;
use fieldcraft_model::geometry::{Entity, GeometryRef};
use fieldcraft_model::units::Quantity;
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// Points lie on one circle when their distances from its centre agree to this fraction of its radius.
const ON_CIRCLE: f64 = 1e-3; // unit-factors: relative to the radius
/// Two unit axes are one direction when the square of the sine between them is below this.
const SAME_AXIS_SINE2: f64 = 1e-12; // unit-factors: the square of the sine of an angle
/// Points this close, relative to the body's size, are on one axis.
const SAME_AXIS: f64 = 1e-6; // unit-factors: relative to the body's size

pub const OPEN: CommandSpec = CommandSpec {
    id: "geometry.open",
    label: "Open Design",
    params: r#"{"path": "part.solidcraft", "model"?: "keep"} — rebuilds the SolidCraft design and makes its bodies the geometry. `model`: "keep" (default) keeps the mesh and every definition, so opening the same design again after an edit picks up the change and sets follow their persistent face names; "new" starts a new model with the design (no mesh, sets, materials, sections, BCs, loads or steps; the unit system stays; the model is untitled until saved, and undo restores the previous file). The reply's `related` is false when the model (before the open) refers to bodies (mesh tags, set geometry) and the design has none of them, and `stale` lists what no longer resolves: `mesh_bodies` (bodies the mesh is tagged with that the design lacks) `sets` (sets naming a body or face it lacks) and `seeds` (edge seeds on edges it lacks), and `references` every definition naming a face, edge or vertex of a body the design has but not that entity ({owner, reference, why}). Such a name is first looked for again as file.open does (`renamed`, `mismatched`). `errors` lists the features that failed to rebuild with the `bodies` each may change, and `warnings` says so in one line each. One undo step"#,
    form: Some(super::forms::OPEN_DESIGN),
    kind: CommandKind::Model,
    requires: &[],
    run: open,
};

pub const BODIES: CommandSpec = CommandSpec {
    id: "geometry.bodies",
    label: "Bodies",
    params: "{} — the open design's bodies: measures, persistent face names (by face index) and edges (persistent name `<face>|<face>`, `~k` when several separate the same faces, and length) and vertices (persistent name `<face>|<face>|<face>` after the faces that meet there, `~k` when several touch the same faces, and point)",
    form: Some(&[]),
    kind: CommandKind::Model,
    requires: &[], // No design is an answer too: no bodies.
    run: bodies,
};

pub const FEATURES: CommandSpec = CommandSpec {
    id: "geometry.features",
    label: "Feature Graph",
    params: r#"{"body"?: "Base" (default: every body)} — the features of the SolidCraft timeline that built each body of the open design, in order (read-only): each feature's `id` (`F<id>` in persistent names), `name`, `type` (extrude, box, cylinder, hole, shell, fillet, combine, pattern, mirror, split, …), `operation` (new_body, join, cut, intersect), `values` (its numeric inputs evaluated: lengths in the display units, radians, counts), `placement` (a primitive's corner, base or centre), `direction` (an extrude's along the sketch normal, a cylinder's axis, a hole's drilling direction), `two_sided`, `sketch` (plane origin, axes and normal, and the profiles used: lines, arcs and circles in sketch coordinates, holes, curve ids), `transforms` (pattern and mirror copies, 4×4 row-major), `targets` and `tools`, the stored `definition`, the persistent `faces` and `edges` of the body it made, and `error` when it failed to rebuild. A partition cell is described by the body it was cut from"#,
    form: Some(super::forms::FEATURES),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: features,
};

pub const CENTRE: CommandSpec = CommandSpec {
    id: "geometry.centre",
    label: "Centre of Selection",
    params: r#"{"geometry": [{"body": "Body1", "face": "F3:hole"} | {"body", "edge"} | {"body", "vertex"}], "nodes": [labels]} — the point a reference point or datum point goes at for what is picked: a vertex's point; the mean of nodes; the centre of a hole (cylindrical faces on one axis, or circular edges on one circle: the axis point at their mid-depth, the circle's centre); otherwise the area centroid of the faces and the length centroid of the edges. Replies `point` and `from` (`vertex`, `nodes`, `circle`, `centroid`)"#,
    form: Some(super::forms::CENTRE),
    kind: CommandKind::Model,
    requires: &[],
    run: centre,
};

pub const CHECK: CommandSpec = CommandSpec {
    id: "geometry.check",
    label: "Check Geometry",
    params: r#"{"body"?: "Body1" (default: every body of the open design), "size"?: length (the target element size; default a twentieth of the smallest checked body's bounding-box diagonal), "short_edge": 0.3, "small_face": 0.1, "sliver": 0.3, "gap": 0.3 (fractions of the size, 0..=10), "parallel_angle": 10 (degrees, 0..=60)} — before meshing, finds what is small against the size (read-only): `short_edge` edges shorter than short_edge × size, `small_face` faces with area below small_face × size², `sliver_face` faces narrower (2 × area / perimeter) than sliver × size, and `near_faces` pairs of faces within parallel_angle of parallel that are closer than gap × size without touching (a gap between bodies or cells, an overlap, or a wall thinner than the gap; faces within twice the tessellation tolerance touch). Replies `size`, `counts` per kind, `summary`, and `findings` (by kind, worst first, at most 1000; `truncated` says when more were found), each with `kind`, `body`, `edge`, `faces` (persistent names; near faces: the second on `other_body`), `value` (a length; an area for `small_face`), `limit` and `at` (a point on it). `mark` (default true): the viewport and pictures mark the findings in the marker colour (their faces, short edges as lines, and a dot at each `at`) until the design changes or a check with `mark: false` or nothing found"#,
    form: Some(super::forms::GEOMETRY_CHECK),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: check,
};

pub const TESSELLATE: CommandSpec = CommandSpec {
    id: "geometry.tessellate",
    label: "Tessellate Body",
    params: r#"{"body": "Body1", "tolerance": 0.05} — triangulates the body's boundary within the chordal tolerance; reports triangles per face"#,
    form: Some(super::forms::TESSELLATE),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: tessellate,
};

pub const IMPORT_MESH: CommandSpec = CommandSpec {
    id: "geometry.import.mesh",
    label: "Import Tagged Mesh",
    params: r#"{"path": "part.3mf"} — adds each build item of a 3MF package as a surface; triangles keep their B-rep face tags (docs/formats/tagged-mesh.md)"#,
    form: Some(super::forms::IMPORT_MESH),
    kind: CommandKind::Model,
    requires: &[],
    run: import_mesh,
};

fn open(s: &mut Session, p: &Value) -> Result<Value> {
    let new_model = match p.get("model") {
        None | Some(Value::Null) => false,
        Some(v) => match v.as_str() {
            Some("keep") => false,
            Some("new") => true,
            _ => return Err(EngineError::Params(r#"`model` must be "keep" or "new""#.into())),
        },
    };
    let path = super::file_param(s, p, "path", Access::Read)?;
    let bytes = super::file::read_file(&path, MAX_DESIGN_BYTES as u64)?;
    let source = SourceFile::of(super::file::absolute(&path), &bytes)?;
    let opened = open_design(&bytes)?;
    // The model's partitions cut the design again (a new model has none), and its names are
    // brought to the design's (#605).
    let (design, names) = if new_model {
        (opened, super::names::Report::default())
    } else {
        let cut = super::datum::partitioned(&opened, &s.model)?;
        super::names::reconcile(&mut s.model, None, &opened, cut)?
    };
    let mut fit = fit(s, &design);
    if new_model {
        let mut model = Model::default();
        model.set_display_units(s.model.display_units().cloned());
        s.model = model;
        s.surfaces.clear();
        s.sources.meshes.clear();
        s.results_view = None;
        // A new model is untitled: saving asks for a file instead of overwriting the old model's
        // (undo brings the old file back with the old model).
        s.set_document(Default::default());
        fit.mesh_bodies.clear();
        fit.sets.clear();
        fit.seeds.clear();
    }
    let out = json!({
        "name": design.name,
        "bodies": design.bodies.iter().map(|b| b.name.as_str()).collect::<Vec<_>>(),
        "errors": design.errors.iter().map(|e| json!({"feature": e.feature, "message": e.message, "bodies": e.bodies})).collect::<Vec<_>>(),
        // Shown with the open: a body that differs from the design is never silent.
        // So is a definition whose name the design no longer has, or one it now names otherwise.
        "warnings": design.errors.iter().map(|e| e.warning()).chain(names.warnings()).collect::<Vec<_>>(),
        "model": if new_model { "new" } else { "kept" },
        "related": fit.related,
        "renamed": names.json()["renamed"],
        "mismatched": names.json()["mismatched"],
        "stale": stale_json(&fit, &names),
    });
    // Another design file (or a new model) is another document; the same file opened again after
    // an edit is the same one.
    if new_model || s.sources.design.as_ref().is_none_or(|d| d.path() != source.path()) {
        s.mark_opened();
    }
    s.design = Some(Arc::new(design));
    s.sources.design = Some(source);
    Ok(out)
}

/// How the model's geometry references fit a design about to become the geometry.
struct Fit {
    /// Whether the design (or an imported surface) has a body the model refers to, or the model
    /// refers to none: `false` means the model is of another part.
    related: bool,
    /// Bodies the mesh is tagged with that neither the design nor a surface has.
    mesh_bodies: Vec<String>,
    /// Sets naming a body or a face that is not there.
    sets: Vec<String>,
    /// Edge seeds on an edge the design lacks.
    seeds: Vec<String>,
}

/// The `stale` object of geometry.open and file.open: what of the model the design does not have.
fn stale_json(fit: &Fit, names: &super::names::Report) -> Value {
    json!({"mesh_bodies": fit.mesh_bodies, "sets": fit.sets, "seeds": fit.seeds, "references": names.json()["references"]})
}

/// file.open's `stale` (see [`stale_json`]), for the session it opened.
pub(super) fn opened_stale(s: &Session, names: &super::names::Report) -> Value {
    match s.design.as_deref() {
        Some(design) => stale_json(&fit(s, design), names),
        None => stale_json(&Fit { related: true, mesh_bodies: Vec::new(), sets: Vec::new(), seeds: Vec::new() }, names),
    }
}

fn fit(s: &Session, design: &Design) -> Fit {
    let has = |body: &str| design.body(body).is_some() || s.surfaces.iter().any(|t| t.name == body);
    let resolves = |r: &GeometryRef| match (&r.entity, design.body(&r.body)) {
        (Entity::Face(face), Some(b)) => b.faces().iter().any(|f| &f.name == face),
        _ => has(&r.body),
    };
    let tagged: BTreeSet<&str> = s.model.classification().entities().iter().map(|r| r.body.as_str()).collect();
    let mut referenced = tagged.clone();
    referenced.extend(s.model.sets().values().flat_map(|set| set.geometry.iter().map(|r| r.body.as_str())));
    Fit {
        related: referenced.is_empty() || referenced.iter().any(|b| has(b)),
        mesh_bodies: tagged.into_iter().filter(|b| !has(b)).map(str::to_string).collect(),
        sets: s.model.sets().iter().filter(|(_, set)| !set.geometry.iter().all(resolves)).map(|(n, _)| n.clone()).collect(),
        seeds: s
            .model
            .seeds()
            .iter()
            .filter(|seed| super::seeds::seeded_edges(design, &seed.edge).map_or(true, |e| e.is_empty()))
            .map(|seed| seed.edge.to_string())
            .collect(),
    }
}

fn bodies(s: &mut Session, _: &Value) -> Result<Value> {
    let Some(design) = &s.design else { return Ok(json!([])) };
    let mut out = Vec::with_capacity(design.bodies.len());
    for b in &design.bodies {
        let m = b.measure()?;
        let edges: Vec<Value> = b.edges()?.iter().map(|e| json!({"name": e.name, "length": e.length()})).collect();
        let vertices: Vec<Value> = b.vertices()?.iter().map(|v| json!({"name": v.name, "point": v.point})).collect();
        out.push(json!({
            "name": b.name,
            "volume": m.volume,
            "area": m.area,
            "centroid": m.centroid,
            "bbox": [m.bbox_min, m.bbox_max],
            "faces": b.faces().iter().map(|f| json!({"name": f.name, "exact": f.exact})).collect::<Vec<_>>(),
            "edges": edges,
            "vertices": vertices,
        }));
    }
    Ok(Value::Array(out))
}

fn features(s: &mut Session, p: &Value) -> Result<Value> {
    let body = match p.get("body") {
        None | Some(Value::Null) => None,
        Some(v) => Some(v.as_str().ok_or_else(|| EngineError::param("body", ProblemKind::NotAName))?),
    };
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let mut graph = design.features(body)?;
    // Its lengths keyed by the feature's own labels, so shown here rather than by the reply
    // table: the design's millimetres in the display system.
    let (display, mm) = (s.units(), crate::units::geometry());
    let show = |q: Quantity, v: f64| mm.to_si(q, v).and_then(|si| display.express(q, si)).unwrap_or(f64::NAN);
    for b in &mut graph {
        b.convert(&|v| show(Quantity::Length, v), &|v| show(Quantity::Area, v));
    }
    Ok(json!({"bodies": graph, "units": {"length": display.symbol(Quantity::Length), "area": display.symbol(Quantity::Area)}}))
}

/// Most geometry items and nodes `geometry.centre` takes.
const MAX_CENTRE_ITEMS: usize = 10_000;

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale3(a: [f64; 3], k: f64) -> [f64; 3] {
    a.map(|v| v * k)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// The centre of the circle through three points, `None` when they lie on a line.
fn circumcentre(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<[f64; 3]> {
    let (u, v) = (sub3(b, a), sub3(c, a));
    let n = cross3(u, v);
    let nn = dot3(n, n);
    if !(nn > 1e-24 * dot3(u, u) * dot3(v, v)) {
        return None;
    }
    let k = add3(scale3(cross3(n, u), dot3(v, v)), scale3(cross3(v, n), dot3(u, u)));
    Some(add3(a, scale3(k, 0.5 / nn)))
}

/// The centre of an edge's polyline when it is an arc or circle (every point at one distance
/// from it, within a thousandth), else `None`.
fn circle_centre(points: &[[f64; 3]]) -> Option<[f64; 3]> {
    let n = points.len();
    if n < 3 {
        return None;
    }
    let closed = points.first().zip(points.last()).is_some_and(|(a, b)| dot3(sub3(*a, *b), sub3(*a, *b)) < 1e-18);
    let pick = |f: f64| points.get(((n - 1) as f64 * f).round() as usize).copied();
    let (a, b, c) = if closed { (pick(0.0)?, pick(1.0 / 3.0)?, pick(2.0 / 3.0)?) } else { (pick(0.0)?, pick(0.5)?, pick(1.0)?) };
    let centre = circumcentre(a, b, c)?;
    let r = dot3(sub3(a, centre), sub3(a, centre)).sqrt();
    let on = points.iter().all(|p| (dot3(sub3(*p, centre), sub3(*p, centre)).sqrt() - r).abs() <= ON_CIRCLE * r);
    on.then_some(centre)
}

/// The length-weighted centre of a polyline and its length.
fn polyline_centroid(points: &[[f64; 3]]) -> ([f64; 3], f64) {
    let (mut sum, mut length) = ([0.0; 3], 0.0);
    for w in points.windows(2) {
        let l = dot3(sub3(w[1], w[0]), sub3(w[1], w[0])).sqrt();
        sum = add3(sum, scale3(add3(w[0], w[1]), 0.5 * l));
        length += l;
    }
    (if length > 0.0 { scale3(sum, 1.0 / length) } else { points.first().copied().unwrap_or_default() }, length)
}

fn centre(s: &mut Session, p: &Value) -> Result<Value> {
    let geometry: Vec<GeometryRef> = super::opt(p, "geometry")?.unwrap_or_default();
    let nodes: Vec<u32> = super::opt(p, "nodes")?.unwrap_or_default();
    if geometry.len() + nodes.len() > MAX_CENTRE_ITEMS {
        return Err(EngineError::Params(format!("geometry.centre takes at most {MAX_CENTRE_ITEMS} items")));
    }
    if geometry.is_empty() && nodes.is_empty() {
        return Err(EngineError::Params("give one of `geometry` or `nodes`: what to find the centre of".into()));
    }
    if !geometry.is_empty() && !nodes.is_empty() {
        return Err(EngineError::Params("give either `geometry` or `nodes`, not both".into()));
    }
    if !nodes.is_empty() {
        let mut sum = [0.0; 3];
        for n in &nodes {
            let node = s.model.node(*n).ok_or_else(|| EngineError::Params(format!("no node {n}")))?;
            sum = add3(sum, node.xyz);
        }
        return Ok(json!({"point": scale3(sum, 1.0 / nodes.len() as f64), "from": "nodes"}));
    }
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let shape = |b: &str| design.body(b).ok_or_else(|| EngineError::Params(format!("the design has no body `{b}`")));
    if let [g] = geometry.as_slice()
        && let Entity::Vertex(v) = &g.entity
    {
        let vertex = shape(&g.body)?.vertices()?.into_iter().find(|x| &x.name == v);
        let vertex = vertex.ok_or_else(|| EngineError::Params(format!("body `{}` has no vertex `{v}`", g.body)))?;
        return Ok(json!({"point": vertex.point, "from": "vertex"}));
    }
    // Hole centres: every face a cylinder about one axis, or every edge an arc of one circle.
    let mut circles: Vec<[f64; 3]> = Vec::new();
    let mut axis: Option<([f64; 3], [f64; 3])> = None;
    let mut round = true;
    // Centroids, weighted by area (faces) or length (edges).
    let (mut faces_sum, mut faces_area, mut edges_sum, mut edges_length) = ([0.0; 3], 0.0, [0.0; 3], 0.0);
    for g in &geometry {
        let body = shape(&g.body)?;
        match &g.entity {
            Entity::Face(f) => {
                let i = body.face_index(f).ok_or_else(|| EngineError::Params(format!("body `{}` has no face `{f}`", g.body)))?;
                let c = body.face_centre(f)?;
                let mesh = body.tessellate((body.size()? * fieldcraft_geometry::EDGE_TOLERANCE).max(1e-4))?;
                let area: f64 = mesh
                    .triangles
                    .iter()
                    .zip(&mesh.tri_face)
                    .filter(|(_, t)| **t == i)
                    .filter_map(|(t, _)| {
                        let q = |k: usize| mesh.positions.get(t[k] as usize).copied();
                        let (a, b, c) = (q(0)?, q(1)?, q(2)?);
                        let n = cross3(sub3(b, a), sub3(c, a));
                        Some(0.5 * dot3(n, n).sqrt())
                    })
                    .sum();
                faces_sum = add3(faces_sum, scale3(c, area));
                faces_area += area;
                match body.face_surfaces()?.get(i as usize).copied().flatten() {
                    Some(fieldcraft_geometry::FaceSurface::Cylinder { origin, axis: d, .. }) => {
                        let len = dot3(d, d).sqrt();
                        let d = scale3(d, 1.0 / len.max(1e-300));
                        // The face centre's foot on the axis.
                        let foot = add3(origin, scale3(d, dot3(sub3(c, origin), d)));
                        let same = axis.is_none_or(|(o, a)| {
                            let off = sub3(foot, o);
                            let across = sub3(off, scale3(a, dot3(off, a)));
                            dot3(cross3(a, d), cross3(a, d)) < SAME_AXIS_SINE2
                                && dot3(across, across).sqrt() <= SAME_AXIS * body.size().unwrap_or(1.0).max(crate::units::min_length())
                        });
                        round &= same;
                        axis.get_or_insert((foot, d));
                        circles.push(foot);
                    }
                    _ => round = false,
                }
            }
            Entity::Edge(e) => {
                let edges = body.edges()?;
                let edge = edges.iter().find(|x| &x.name == e).ok_or_else(|| EngineError::Params(format!("body `{}` has no edge `{e}`", g.body)))?;
                let (c, l) = polyline_centroid(&edge.points);
                edges_sum = add3(edges_sum, scale3(c, l));
                edges_length += l;
                match circle_centre(&edge.points) {
                    Some(c) => circles.push(c),
                    None => round = false,
                }
            }
            Entity::Vertex(_) | Entity::Body => {
                return Err(EngineError::Params("geometry.centre takes one vertex, or faces and edges".into()));
            }
        }
    }
    // Arcs of one circle share a centre; faces of one axis were checked as they came.
    let spread = circles.iter().map(|c| circles.first().map_or(0.0, |f| dot3(sub3(*c, *f), sub3(*c, *f)).sqrt())).fold(0.0f64, f64::max);
    let size = geometry.first().and_then(|g| design.body(&g.body)).and_then(|b| b.size().ok()).unwrap_or(1.0).max(crate::units::min_length());
    let edges_only = geometry.iter().all(|g| matches!(g.entity, Entity::Edge(_)));
    if round && !circles.is_empty() && (!edges_only || spread <= SAME_AXIS * size) {
        let point = match axis {
            // Cylinder faces: the axis point at the mean of their centres' feet.
            Some((_, d)) => {
                let mean = scale3(circles.iter().fold([0.0; 3], |a, c| add3(a, *c)), 1.0 / circles.len() as f64);
                let o = circles.first().copied().unwrap_or(mean);
                add3(o, scale3(d, dot3(sub3(mean, o), d)))
            }
            None => circles.first().copied().unwrap_or_default(),
        };
        return Ok(json!({"point": point, "from": "circle"}));
    }
    let point = if faces_area > 0.0 {
        scale3(faces_sum, 1.0 / faces_area)
    } else if edges_length > 0.0 {
        scale3(edges_sum, 1.0 / edges_length)
    } else {
        return Err(EngineError::Params("what is picked has no area or length to find a centre of".into()));
    };
    if !point.iter().all(|v| v.is_finite()) {
        return Err(EngineError::Params("the centre is not a finite point".into()));
    }
    Ok(json!({"point": point, "from": "centroid"}))
}

fn body<'a>(s: &'a Session, p: &Value) -> Result<&'a Shape> {
    let name = p.get("body").and_then(Value::as_str).ok_or_else(|| EngineError::param("body", ProblemKind::NotAName))?;
    s.design.as_deref().and_then(|d| d.body(name)).ok_or_else(|| EngineError::Params(format!("no body `{name}` (open a design first)")))
}

/// The body named by `body`, tessellated within `tolerance`.
pub(crate) fn body_mesh<'a>(s: &'a Session, p: &Value) -> Result<(&'a Shape, SurfaceMesh)> {
    let shape = body(s, p)?;
    let tolerance = p.get("tolerance").and_then(Value::as_f64).ok_or_else(|| EngineError::param("tolerance", ProblemKind::NotANumber))?;
    Ok((shape, shape.tessellate(tolerance)?))
}

/// Most findings `geometry.check` lists.
const MAX_FINDINGS: usize = 1000;

fn check(s: &mut Session, p: &Value) -> Result<Value> {
    use fieldcraft_geometry::check::{self, CheckOptions, Kind};
    let shapes: Vec<&Shape> = match p.get("body") {
        Some(_) => vec![body(s, p)?],
        None => s.design.as_deref().map(|d| d.bodies.iter().collect()).unwrap_or_default(),
    };
    if shapes.is_empty() {
        return Err(EngineError::Params("the design has no bodies to check".into()));
    }
    let size = match p.get("size") {
        Some(v) => v.as_f64().filter(|x| *x > 0.0 && x.is_finite()).ok_or_else(|| EngineError::param("size", ProblemKind::NotPositive))?,
        None => {
            let mut smallest = f64::INFINITY;
            for b in &shapes {
                smallest = smallest.min(b.size()?);
            }
            smallest / 20.0
        }
    };
    let fraction = |key: &str, default: f64| super::mesh::number_in(p, key, default, 0.0..=10.0);
    let defaults = CheckOptions::new(size);
    let options = CheckOptions {
        size,
        short_edge: fraction("short_edge", defaults.short_edge)?,
        small_face: fraction("small_face", defaults.small_face)?,
        sliver: fraction("sliver", defaults.sliver)?,
        gap: fraction("gap", defaults.gap)?,
        parallel_angle: super::mesh::number_in(p, "parallel_angle", defaults.parallel_angle.to_degrees(), 0.0..=60.0)?.to_radians(),
    };
    let mark = match p.get("mark") {
        None => true,
        Some(v) => v.as_bool().ok_or_else(|| EngineError::Params("`mark` must be true or false".into()))?,
    };
    let findings = check::check(&shapes, &options)?;
    s.geometry_marks = s.design.as_ref().filter(|_| mark && !findings.is_empty()).map(|design| {
        let marks = findings
            .iter()
            .map(|f| {
                // A near pair's second face lies on the other body.
                let on = |k: usize| if k == 1 { f.other_body.clone().unwrap_or_else(|| f.body.clone()) } else { f.body.clone() };
                // A short edge is drawn itself: its faces can be large.
                let faces = if f.kind == Kind::ShortEdge {
                    Vec::new()
                } else if f.kind == Kind::NearFaces {
                    f.face_indices.iter().enumerate().map(|(k, &i)| (on(k), i)).collect()
                } else {
                    f.face_indices.iter().map(|&i| (f.body.clone(), i)).collect()
                };
                let small =
                    matches!(f.kind, Kind::SmallFace | Kind::SliverFace).then(|| f.faces.first().map(|n| (f.body.clone(), n.clone()))).flatten();
                crate::view::GeometryMark { faces, path: f.path.clone(), at: f.at, small }
            })
            .collect();
        Arc::new(crate::view::GeometryMarks { design: Arc::downgrade(design), marks })
    });
    let counts: serde_json::Map<String, Value> =
        Kind::ALL.iter().map(|k| (k.name().to_string(), json!(findings.iter().filter(|f| f.kind == *k).count()))).collect();
    let shown = crate::units::shown_length(s, size);
    let summary = if findings.is_empty() {
        format!("nothing small against size {shown}")
    } else {
        let parts: Vec<String> = Kind::ALL
            .iter()
            .filter_map(|k| {
                let n = findings.iter().filter(|f| f.kind == *k).count();
                (n > 0).then(|| format!("{n} {}", k.name().replace('_', " ")))
            })
            .collect();
        format!("{} against size {shown}", parts.join(", "))
    };
    let listed: Vec<Value> = findings
        .iter()
        .take(MAX_FINDINGS)
        .map(|f| {
            json!({
                "kind": f.kind.name(),
                "body": f.body,
                "edge": f.edge,
                "faces": f.faces,
                "other_body": f.other_body,
                "value": f.value,
                "limit": f.limit,
                "at": f.at,
            })
        })
        .collect();
    Ok(json!({
        "size": size,
        "counts": counts,
        "summary": summary,
        "findings": listed,
        "truncated": findings.len() > MAX_FINDINGS,
    }))
}

fn tessellate(s: &mut Session, p: &Value) -> Result<Value> {
    let (shape, mesh) = body_mesh(s, p)?;
    let mut per_face = vec![0usize; shape.faces().len()];
    for &f in &mesh.tri_face {
        if let Some(n) = per_face.get_mut(f as usize) {
            *n += 1;
        }
    }
    let faces: Vec<Value> = shape.faces().iter().zip(per_face).map(|(f, n)| json!({"name": f.name, "triangles": n})).collect();
    Ok(json!({"body": shape.name, "vertices": mesh.positions.len(), "triangles": mesh.triangles.len(), "faces": faces}))
}

fn import_mesh(s: &mut Session, p: &Value) -> Result<Value> {
    let path = super::file_param(s, p, "path", Access::Read)?;
    if s.sources.meshes.len() >= fieldcraft_io::MAX_MESH_SOURCES {
        return Err(EngineError::Params(format!("a model imports at most {} meshes", fieldcraft_io::MAX_MESH_SOURCES)));
    }
    let bytes = super::file::read_file(&path, tagged_mesh::MAX_FILE_BYTES as u64)?;
    let source = SourceFile::of(super::file::absolute(&path), &bytes)?;
    let surfaces = tagged_mesh::read_3mf(&bytes)?;
    let first = s.surfaces.len();
    let added: Vec<Value> = surfaces
        .iter()
        .enumerate()
        .map(|(k, t)| json!({"index": first + k, "name": t.name, "vertices": t.mesh.positions.len(), "triangles": t.mesh.triangles.len(), "faces": t.faces}))
        .collect();
    s.surfaces.extend(surfaces.into_iter().map(Arc::new));
    s.sources.meshes.push(source);
    // The model's stitches of surfaces by these names apply to them too.
    let (_, failed) = super::stitch::restitch(s);
    let mut reply = json!({"path": path, "surfaces": added});
    let warnings: Vec<String> =
        failed.into_iter().filter(|(_, why)| !why.starts_with("no surface")).map(|(n, why)| format!("stitch `{n}`: {why}")).collect();
    if !warnings.is_empty()
        && let Some(o) = reply.as_object_mut()
    {
        o.insert("warnings".into(), json!(warnings));
    }
    Ok(reply)
}

#[cfg(test)]
pub(crate) mod tests {
    use std::io::{Cursor, Write};

    use super::*;

    /// A session with the design the SolidCraft `commands` build open.
    pub(crate) fn design(commands: Value) -> (crate::tests::TempDir, Session) {
        let bytes = fieldcraft_geometry::script_design(&json!({"commands": commands}).to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    fn kinds(r: &Value) -> Vec<&str> {
        r["findings"].as_array().unwrap().iter().map(|f| f["kind"].as_str().unwrap()).collect()
    }

    /// A plain block is clean at a size it can be meshed with; a 0.4 mm step on its top is a
    /// sliver riser with two short edges, located on the step; the check changes nothing.
    #[test]
    fn check_finds_a_small_step() {
        let (_dir, mut s) = design(json!([{"command": "solid.box", "params": {"length": 40, "width": 20, "height": 10, "body_name": "Base"}}]));
        let clean = s.execute("geometry.check", &json!({"size": 5})).unwrap();
        assert!(kinds(&clean).is_empty(), "{clean}");
        assert_eq!(clean["summary"], "nothing small against size 5 mm");
        let (_dir, mut s) = design(json!([
            {"command": "solid.box", "params": {"length": 40, "width": 20, "height": 10, "body_name": "Base"}},
            {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 10.4, "operation": "join", "targets": ["Base"]}},
        ]));
        let journal = s.journal.len();
        let r = s.execute("geometry.check", &json!({"body": "Base", "size": 5})).unwrap();
        assert_eq!(s.journal.len(), journal, "read-only");
        assert_eq!(r["counts"]["short_edge"], 2, "{r}");
        assert_eq!(r["counts"]["sliver_face"], 1, "{r}");
        assert_eq!(r["counts"]["small_face"], 0, "{r}");
        assert_eq!(r["counts"]["near_faces"], 0, "{r}");
        let names: Vec<String> = s.execute("geometry.bodies", &json!({})).unwrap()[0]["faces"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| f["name"].as_str().unwrap().to_string())
            .collect();
        for f in r["findings"].as_array().unwrap() {
            let at: Vec<f64> = f["at"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
            assert!((at[0] - 20.0).abs() < 1e-6 && at[2] >= 10.0 - 1e-6 && at[2] <= 10.4 + 1e-6, "on the riser: {f}");
            assert!(f["faces"].as_array().unwrap().iter().all(|n| names.contains(&n.as_str().unwrap().to_string())), "{f}");
            assert!(f["value"].as_f64().unwrap() < f["limit"].as_f64().unwrap());
        }
        let short = &r["findings"][0];
        assert!((short["value"].as_f64().unwrap() - 0.4).abs() < 1e-6, "{short}");
        assert!(short["edge"].is_string());
        // The size is what makes the step small: at 1 mm it is a feature like any other.
        let fine = s.execute("geometry.check", &json!({"size": 1})).unwrap();
        assert!(kinds(&fine).is_empty(), "{fine}");
    }

    /// The findings are drawn in the picture: the sliver riser's faces in the marker colour, the
    /// short edges as lines and a dot at each finding over everything; a clean check, `mark:
    /// false` or another design takes them away.
    #[test]
    fn check_marks_its_findings_in_the_picture() {
        let (_dir, mut s) = design(json!([
            {"command": "solid.box", "params": {"length": 40, "width": 20, "height": 10, "body_name": "Base"}},
            {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 10.4, "operation": "join", "targets": ["Base"]}},
        ]));
        let p = crate::view::Palette::LIGHT;
        let marked = |s: &Session| {
            let (scene, _, _) = crate::view::scene(s, &p);
            let faces = scene.tris.iter().filter(|t| t.color == p.marked).count();
            let lines: Vec<_> = scene.on_top.iter().filter(|l| l.color == p.marked).collect();
            let dots = lines.iter().filter(|l| l.a == l.b).count();
            (faces, lines.len() - dots, dots)
        };
        assert_eq!(marked(&s), (0, 0, 0), "nothing checked yet");
        let r = s.execute("geometry.check", &json!({"size": 5})).unwrap();
        let (faces, edges, dots) = marked(&s);
        assert!(faces > 0, "the riser's faces are marked");
        assert!(edges >= 2, "the short edges are drawn: {edges}");
        assert_eq!(dots, r["findings"].as_array().unwrap().len(), "a dot at each finding");
        // Every marked triangle lies on the riser, x = 20, between z = 10 and 10.4.
        let (scene, _, _) = crate::view::scene(&s, &p);
        // The scene is in metres, the design in millimetres.
        let mm = crate::units::to_geometry_length;
        for t in scene.tris.iter().filter(|t| t.color == p.marked) {
            assert!(t.p.iter().all(|v| (mm(v.x) - 20.0).abs() < 1e-6 && mm(v.z) > 10.0 - 1e-6 && mm(v.z) < 10.4 + 1e-6), "{t:?}");
        }
        s.execute("geometry.check", &json!({"size": 5, "mark": false})).unwrap();
        assert_eq!(marked(&s), (0, 0, 0), "mark: false clears");
        s.execute("geometry.check", &json!({"size": 5})).unwrap();
        s.execute("geometry.check", &json!({"size": 1})).unwrap();
        assert_eq!(marked(&s), (0, 0, 0), "a clean check clears");
        s.execute("geometry.check", &json!({"size": 5})).unwrap();
        assert!(s.execute("geometry.check", &json!({"size": 5, "mark": "yes"})).is_err());
        assert!(marked(&s).2 > 0, "a failed check keeps the marks");
        s.design = s.design.as_deref().map(|d| Arc::new(d.clone()));
        assert_eq!(marked(&s), (0, 0, 0), "stale once another design is open");
    }

    /// Two bodies 0.5 mm apart: their facing faces are near faces; once they touch they are not.
    #[test]
    fn check_finds_a_gap_between_bodies() {
        let (_dir, mut s) = design(json!([
            {"command": "solid.box", "params": {"length": 10, "width": 10, "height": 10, "body_name": "A"}},
            {"command": "solid.box", "params": {"length": 10, "width": 10, "height": 10, "corner": [0, 0, 10.5], "body_name": "B"}},
        ]));
        let r = s.execute("geometry.check", &json!({"size": 2.5})).unwrap();
        assert_eq!(kinds(&r), ["near_faces"], "{r}");
        let f = &r["findings"][0];
        assert!((f["value"].as_f64().unwrap() - 0.5).abs() < 1e-6, "{f}");
        let bodies = [f["body"].as_str().unwrap(), f["other_body"].as_str().unwrap()];
        assert!(bodies.contains(&"A") && bodies.contains(&"B"), "{f}");
        assert!(s.execute("geometry.check", &json!({"size": 2.5, "gap": 0.1})).unwrap()["findings"].as_array().unwrap().is_empty());
        // Without a size the default is a twentieth of the smaller diagonal.
        let d = s.execute("geometry.check", &json!({})).unwrap();
        assert!((d["size"].as_f64().unwrap() - 300f64.sqrt() / 20.0).abs() < 1e-9, "{d}");
    }

    #[test]
    fn check_rejects_hostile_parameters() {
        let mut s = Session::default();
        assert!(s.execute("geometry.check", &json!({})).is_err(), "needs a design");
        let (_dir, mut s) = design(json!([{"command": "solid.box", "params": {"length": 10, "width": 10, "height": 10, "body_name": "A"}}]));
        for p in [
            json!({"size": 0}),
            json!({"size": -1}),
            json!({"size": "big"}),
            json!({"size": 1e308, "gap": 1e308}),
            json!({"short_edge": -1}),
            json!({"sliver": 1e9}),
            json!({"parallel_angle": 90}),
            json!({"body": "Nope"}),
            json!({"body": 3}),
        ] {
            assert!(s.execute("geometry.check", &p).is_err(), "{p}");
        }
    }

    pub(crate) fn package(model: &str) -> Vec<u8> {
        let mut zip = zip::ZipWriter::new(Cursor::new(Vec::new()));
        zip.start_file("3D/3dmodel.model", zip::write::SimpleFileOptions::default()).unwrap();
        zip.write_all(model.as_bytes()).unwrap();
        zip.finish().unwrap().into_inner()
    }

    /// The plate's hole (Ø10 at x = y = 15 through z 0..10): its faces give the axis point at
    /// mid-depth, a rim arc its circle's centre, a corner vertex its point, a flat face its
    /// centroid; datum points take the design's vertices by name.
    #[test]
    fn centres_of_holes_vertices_and_faces() {
        use fieldcraft_geometry::FaceSurface;
        use fieldcraft_geometry::sample::{Plate, plate_design};
        let design = Arc::new(open_design(&plate_design(&Plate::default()).unwrap()).unwrap());
        let mut s = Session { design: Some(design.clone()), ..Session::default() };
        let body = &design.bodies[0];
        let near = |v: &Value, p: [f64; 3]| {
            let q: [f64; 3] = serde_json::from_value(v["point"].clone()).unwrap();
            assert!((0..3).all(|k| (q[k] - p[k]).abs() < 1e-6), "{v} vs {p:?}");
        };
        let surfaces = body.face_surfaces().unwrap();
        let hole: Vec<Value> = body
            .faces()
            .iter()
            .zip(&surfaces)
            .filter(|(_, s)| matches!(s, Some(FaceSurface::Cylinder { .. })))
            .map(|(f, _)| json!({"body": body.name, "face": f.name}))
            .collect();
        assert_eq!(hole.len(), 2);
        let r = s.execute("geometry.centre", &json!({"geometry": hole})).unwrap();
        near(&r, [15.0, 15.0, 5.0]);
        assert_eq!(r["from"], "circle");
        let r = s.execute("geometry.centre", &json!({"geometry": [hole[0]]})).unwrap();
        near(&r, [15.0, 15.0, 5.0]);
        // A rim arc on top.
        let edges = body.edges().unwrap();
        let arc = edges.iter().find(|e| e.points.iter().all(|p| (p[2] - 10.0).abs() < 1e-9) && e.points.len() > 3).unwrap();
        let r = s.execute("geometry.centre", &json!({"geometry": [{"body": body.name, "edge": arc.name}]})).unwrap();
        near(&r, [15.0, 15.0, 10.0]);
        assert_eq!(r["from"], "circle");
        // The corner at the origin, and a datum point on it by name.
        let corner = body.vertices().unwrap().into_iter().find(|v| v.point == [0.0; 3]).unwrap();
        let r = s.execute("geometry.centre", &json!({"geometry": [{"body": body.name, "vertex": corner.name}]})).unwrap();
        near(&r, [0.0; 3]);
        s.execute("geometry.datum.point", &json!({"name": "P", "by": {"type": "vertex", "body": body.name, "vertex": corner.name}})).unwrap();
        assert_eq!(s.model.datum("P").and_then(|d| d.point()), Some([0.0; 3]));
        let listed = s.execute("geometry.bodies", &json!({})).unwrap();
        assert!(listed[0]["vertices"].as_array().unwrap().iter().any(|v| v["name"] == corner.name.as_str()), "{listed}");
        // The bottom face (z = 0): its centroid, off the hole.
        let bottom = body
            .faces()
            .iter()
            .zip(&surfaces)
            .find(|(_, s)| matches!(s, Some(FaceSurface::Plane { normal, point }) if normal[2].abs() > 0.99 && point[2].abs() < 1e-9))
            .unwrap()
            .0;
        let r = s.execute("geometry.centre", &json!({"geometry": [{"body": body.name, "face": bottom.name}]})).unwrap();
        assert_eq!(r["from"], "centroid");
        // Nodes, and what is refused.
        // The model is stored in metres: these are 1, 2, 3 mm and 3, 2, 1 mm.
        s.model.set_node(1, [0.001, 0.002, 0.003]).unwrap();
        s.model.set_node(2, [0.003, 0.002, 0.001]).unwrap();
        near(&s.execute("geometry.centre", &json!({"nodes": [1, 2]})).unwrap(), [2.0, 2.0, 2.0]);
        for bad in [
            json!({}),
            json!({"nodes": [99]}),
            json!({"geometry": [{"body": "nope"}]}),
            json!({"geometry": [{"body": body.name}]}),
            json!({"nodes": "x"}),
        ] {
            assert!(s.execute("geometry.centre", &bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn import_adds_tagged_surfaces_and_failures_change_nothing() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-import-mesh-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let model = format!(
            r#"<model unit="millimeter" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" xmlns:ft="{}"><resources>
<object id="1" name="Tet" ft:faces="2"><mesh><vertices><vertex x="0" y="0" z="0"/><vertex x="1" y="0" z="0"/><vertex x="0" y="1" z="0"/><vertex x="0" y="0" z="1"/></vertices>
<triangles><triangle v1="0" v2="2" v3="1" ft:face="0"/><triangle v1="0" v2="1" v3="3" ft:face="1"/><triangle v1="1" v2="2" v3="3" ft:face="1"/><triangle v1="0" v2="3" v3="2" ft:face="1"/></triangles></mesh></object>
</resources><build><item objectid="1"/></build></model>"#,
            tagged_mesh::FACE_TAGS_NS
        );
        let good = dir.join("tet.3mf");
        std::fs::write(&good, package(&model)).unwrap();
        let junk = dir.join("junk.3mf");
        std::fs::write(&junk, b"PK\x03\x04 not really").unwrap();

        let mut s = Session::default();
        let r = s.execute("geometry.import.mesh", &json!({"path": good.to_str().unwrap()})).unwrap();
        assert_eq!(r["surfaces"][0], json!({"index": 0, "name": "Tet", "vertices": 4, "triangles": 4, "faces": 2}));
        assert_eq!(s.surfaces[0].mesh.tri_face, vec![0, 1, 1, 1]);
        assert_eq!(s.journal.len(), 1);

        let hostile = [
            json!({"path": junk.to_str().unwrap()}),
            json!({"path": dir.join("missing.3mf").to_str().unwrap()}),
            json!({"path": dir.to_str().unwrap()}),
            json!({"path": "/dev/zero"}),
            json!({"path": ""}),
            json!({"path": "x".repeat(5000)}),
        ];
        for p in hostile {
            assert!(s.execute("geometry.import.mesh", &p).is_err(), "{p}");
        }
        assert_eq!(s.surfaces.len(), 1);
        assert_eq!(s.journal.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
