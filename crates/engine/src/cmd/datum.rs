//! Datum planes, axes and points, and partitions of the open design by datum planes.
//!
//! Datums and partitions are part of the model (`fieldcraft_model::datums`). A datum point is
//! placed by coordinates, at a mesh vertex, along a design edge, at a design face's centre or
//! where a datum axis crosses a datum plane; planes and axes can be made through datum points.
//! Each takes the coordinates of what it is made from when it is made (a later edit of a datum
//! point leaves the planes and axes made through it where they are). A partition cuts the
//! design in memory (`fieldcraft_geometry::partition`); the model's partitions are applied again,
//! in order, whenever the design is opened (`geometry.open`, `file.open`, samples), so an edit
//! of the design followed by a re-open replays them.

use crate::schema::ProblemKind;
use std::sync::Arc;

use fieldcraft_geometry::{Cut, Design, Merge};
use fieldcraft_io::{Item, ItemKind, MAX_HIDDEN};
use fieldcraft_model::{Datum, Model, Partition};
use serde::Deserialize;
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// An axis crosses a plane when the sine of their angle is at least this.
const MIN_SINE: f64 = 1e-9; // unit-factors: the sine of an angle

pub const DATUM_PLANE: CommandSpec = CommandSpec {
    id: "geometry.datum",
    label: "Datum Plane",
    params: r#"{"name": "Mid", "by": {"type": "offset", "plane": "yz" | "zx" | "xy", "offset": 15} | {"type": "three_points", "a": [x, y, z], "b": […], "c": […]} | {"type": "point_normal", "point": [x, y, z], "normal": [nx, ny, nz]} | {"type": "datum_points", "a": "P1", "b": "P2", "c": "P3"} | {"type": "normal_to_axis", "axis": "Spindle", "point"?: "P1"}} — a named datum plane (model coordinates): a principal plane moved along its normal (yz: normal +x, zx: normal +y, xy: normal +z), through three points or three datum points (normal (b − a) × (c − a)), through a point with a normal, or normal to a datum axis through a datum point (the axis's point by default). Partitions (geometry.partition) cut by it; shown in the viewport"#,
    form: Some(super::forms::DATUM_PLANE),
    kind: CommandKind::Model,
    requires: &[],
    run: datum_plane,
};

pub const DATUM_AXIS: CommandSpec = CommandSpec {
    id: "geometry.datum.axis",
    label: "Datum Axis",
    params: r#"{"name": "Spindle", "by": {"type": "principal", "axis": "x" | "y" | "z", "through"?: [x, y, z] (default the origin)} | {"type": "two_points", "a": [x, y, z], "b": […]} | {"type": "point_direction", "point": [x, y, z], "direction": [dx, dy, dz]} | {"type": "datum_points", "a": "P1", "b": "P2"}} — a named datum axis, shown in the viewport"#,
    form: Some(super::forms::DATUM_AXIS),
    kind: CommandKind::Model,
    requires: &[],
    run: datum_axis,
};

pub const DATUM_POINT: CommandSpec = CommandSpec {
    id: "geometry.datum.point",
    label: "Datum Point",
    params: r#"{"name": "P1", "by": {"type": "coordinates", "point": [x, y, z]} | {"type": "vertex", "body": "Body1", "vertex": "F1:top|F2:side:l1|F2:side:l2"} | {"type": "edge", "body": "Body1", "edge": "F1:top|F2:side", "fraction"?: 0.5} | {"type": "face_centre", "body": "Body1", "face": "F2:top"} | {"type": "axis_plane", "axis": "Spindle", "plane": "Mid"}} — a named datum point, shown in the viewport as a small cross: at coordinates, at a vertex of the open design (its persistent name, as geometry.bodies lists it and Select: Vertices picks it) or of the mesh (as its geometry tags name it), a fraction along an edge of the open design (0.5: its midpoint), at the area centre of a face of the open design, or where a datum axis crosses a datum plane. Datum planes and axes can be made through datum points"#,
    form: Some(super::forms::DATUM_POINT),
    kind: CommandKind::Model,
    requires: &[],
    run: datum_point,
};

pub const DATUM_EDIT: CommandSpec = CommandSpec {
    id: "geometry.datum.edit",
    label: "Edit Datum Plane",
    params: r#"{"name": "Mid", "by": …} — redefines a datum plane (the parameters of geometry.datum); the partitions that cut by it cut again"#,
    form: Some(super::forms::DATUM_PLANE),
    kind: CommandKind::Model,
    requires: &[Requirement::DatumPlane],
    run: |s, p| datum_edit(s, p, "plane"),
};

pub const DATUM_AXIS_EDIT: CommandSpec = CommandSpec {
    id: "geometry.datum.axis.edit",
    label: "Edit Datum Axis",
    params: r#"{"name": "Spindle", "by": …} — redefines a datum axis (the parameters of geometry.datum.axis)"#,
    form: Some(super::forms::DATUM_AXIS),
    kind: CommandKind::Model,
    requires: &[Requirement::DatumAxis],
    run: |s, p| datum_edit(s, p, "axis"),
};

pub const DATUM_POINT_EDIT: CommandSpec = CommandSpec {
    id: "geometry.datum.point.edit",
    label: "Edit Datum Point",
    params: r#"{"name": "P1", "by": …} — redefines a datum point (the parameters of geometry.datum.point); planes and axes made through it stay where they are"#,
    form: Some(super::forms::DATUM_POINT),
    kind: CommandKind::Model,
    requires: &[Requirement::DatumPoint],
    run: |s, p| datum_edit(s, p, "point"),
};

pub const DATUM_DELETE: CommandSpec = CommandSpec {
    id: "geometry.datum.delete",
    label: "Delete Datum",
    params: r#"{"name": "Mid", "cascade"?: false} — deletes a datum plane, axis or point; one a partition cuts by is refused (naming the partitions), unless `cascade` deletes those partitions too (their cells lose their meshes: `replaced`)"#,
    form: Some(super::forms::DATUM_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Datum],
    run: datum_delete,
};

pub const PARTITION: CommandSpec = CommandSpec {
    id: "geometry.partition",
    label: "Partition",
    params: r#"{"name": "Half", "datum": "Mid", "bodies"?: ["Body1"] (default: every body the plane divides; a body earlier partitions cut stands for its cells)} — cuts bodies of the open design into cells by a datum plane, in memory (the design file is not changed); the model keeps the partition and applies it again whenever the design is opened, after the partitions before it. A cut body `B` becomes the cells `B/a` (on the plane's negative side) and `B/b` (several pieces on a side: `B/a1`, `B/a2`, …). The two faces of the cut are paired by name: `partition:<name>:a` on the cells on the negative side, `partition:<name>:b` on the positive side (`…#0`, `…#1` where a cell has several); every other face keeps the persistent name of the face it is a piece of (`<face>#k` for several pieces on one cell), so sets by face name still select it. Reports the cells made and the cut faces named, and `warnings` when the plane lies on a face of a body it cuts or a cell's boundary does not tessellate into one a mesher takes (a cut through existing edges or vertices, such as a hole's axis), naming the plane, the cell and where; a body it cuts loses its mesh (`replaced`, as a mesher reports it; remesh the cells); a plane that divides no body, a listed body it does not divide, or a body the kernel fails to cut (its message) is an error"#,
    form: Some(super::forms::PARTITION),
    kind: CommandKind::Model,
    requires: &[Requirement::Design, Requirement::DatumPlane],
    run: partition,
};

pub const DATUM_SUPPRESS: CommandSpec = CommandSpec {
    id: "geometry.datum.suppress",
    label: "Suppress Datum",
    params: r#"{"name": "Mid", "suppressed"?: true} — keeps the datum but leaves it out of the model (false brings it back): it cannot be an input to another datum, a partition, a mirror or a symmetry while suppressed, and the partitions that cut by it stop cutting (the design is cut again without them)"#,
    form: Some(super::forms::DATUM_SUPPRESS),
    kind: CommandKind::Model,
    requires: &[Requirement::Datum],
    run: datum_suppress,
};

pub const PARTITION_SUPPRESS: CommandSpec = CommandSpec {
    id: "geometry.partition.suppress",
    label: "Suppress Partition",
    params: r#"{"name": "Half", "suppressed"?: true} — keeps the partition but leaves it out (false brings it back): the design is cut again without it, so its cells merge back into the bodies they came from and lose their meshes (`replaced`)"#,
    form: Some(super::forms::PARTITION_SUPPRESS),
    kind: CommandKind::Model,
    requires: &[Requirement::Partition],
    run: partition_suppress,
};

pub const PARTITION_DELETE: CommandSpec = CommandSpec {
    id: "geometry.partition.delete",
    label: "Delete Partition",
    params: r#"{"name": "Half"} — removes a partition; the design is cut again by the ones left, and the cells that are gone lose their meshes (`replaced`)"#,
    form: Some(super::forms::PARTITION_DELETE),
    kind: CommandKind::Model,
    requires: &[Requirement::Partition],
    run: partition_delete,
};

#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum Principal {
    Yz,
    Zx,
    Xy,
}

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
enum PlaneBy {
    Offset {
        plane: Principal,
        #[serde(default)]
        offset: f64,
    },
    ThreePoints {
        a: [f64; 3],
        b: [f64; 3],
        c: [f64; 3],
    },
    PointNormal {
        point: [f64; 3],
        normal: [f64; 3],
    },
    DatumPoints {
        a: String,
        b: String,
        c: String,
    },
    NormalToAxis {
        axis: String,
        #[serde(default)]
        point: Option<String>,
    },
}

#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum Axis {
    X,
    Y,
    Z,
}

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
enum AxisBy {
    Principal {
        axis: Axis,
        #[serde(default)]
        through: [f64; 3],
    },
    TwoPoints {
        a: [f64; 3],
        b: [f64; 3],
    },
    PointDirection {
        point: [f64; 3],
        direction: [f64; 3],
    },
    DatumPoints {
        a: String,
        b: String,
    },
}

fn half() -> f64 {
    0.5
}

#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
enum PointBy {
    Coordinates {
        point: [f64; 3],
    },
    Vertex {
        body: String,
        vertex: String,
    },
    Edge {
        body: String,
        edge: String,
        #[serde(default = "half")]
        fraction: f64,
    },
    FaceCentre {
        body: String,
        face: String,
    },
    AxisPlane {
        axis: String,
        plane: String,
    },
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The coordinates of datum point `name`.
fn datum_point_at(m: &Model, name: &str) -> Result<[f64; 3]> {
    let d = m.usable_datum(name)?;
    d.point().ok_or_else(|| EngineError::Params(format!("datum `{name}` is a datum {}, not a point", d.kind())))
}

/// The point and unit direction of datum axis `name`.
fn datum_axis_of(m: &Model, name: &str) -> Result<([f64; 3], [f64; 3])> {
    let d = m.usable_datum(name)?;
    d.axis().ok_or_else(|| EngineError::Params(format!("datum `{name}` is a datum {}, not an axis", d.kind())))
}

/// The point and unit normal of datum plane `name`.
pub(crate) fn datum_plane_of(m: &Model, name: &str) -> Result<([f64; 3], [f64; 3])> {
    let d = m.usable_datum(name)?;
    d.plane().ok_or_else(|| EngineError::Params(format!("datum `{name}` is a datum {}, not a plane", d.kind())))
}

/// A plane through three points; collinear points are refused relative to their spread.
fn three_points(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Result<Datum> {
    let (u, v) = (sub(b, a), sub(c, a));
    let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    let size = u.iter().chain(&v).map(|x| x * x).sum::<f64>();
    if !(n.iter().map(|x| x * x).sum::<f64>() > 1e-20 * size * size) {
        return Err(EngineError::Params("the three points must not lie on one line".into()));
    }
    Ok(Datum::Plane { point: a, normal: n })
}

fn plane_datum(m: &Model, by: PlaneBy) -> Result<Datum> {
    Ok(match by {
        PlaneBy::Offset { plane, offset } => {
            if !offset.is_finite() {
                return Err(EngineError::param("offset", ProblemKind::NotANumber));
            }
            let normal = match plane {
                Principal::Yz => [1.0, 0.0, 0.0],
                Principal::Zx => [0.0, 1.0, 0.0],
                Principal::Xy => [0.0, 0.0, 1.0],
            };
            Datum::Plane { point: normal.map(|n| n * offset), normal }
        }
        PlaneBy::ThreePoints { a, b, c } => three_points(a, b, c)?,
        PlaneBy::PointNormal { point, normal } => Datum::Plane { point, normal },
        PlaneBy::DatumPoints { a, b, c } => three_points(datum_point_at(m, &a)?, datum_point_at(m, &b)?, datum_point_at(m, &c)?)?,
        PlaneBy::NormalToAxis { axis, point } => {
            let (on_axis, direction) = datum_axis_of(m, &axis)?;
            let point = match point {
                Some(p) => datum_point_at(m, &p)?,
                None => on_axis,
            };
            Datum::Plane { point, normal: direction }
        }
    })
}

fn axis_datum(m: &Model, by: AxisBy) -> Result<Datum> {
    Ok(match by {
        AxisBy::Principal { axis, through } => {
            let direction = match axis {
                Axis::X => [1.0, 0.0, 0.0],
                Axis::Y => [0.0, 1.0, 0.0],
                Axis::Z => [0.0, 0.0, 1.0],
            };
            Datum::Axis { point: through, direction }
        }
        AxisBy::TwoPoints { a, b } => Datum::Axis { point: a, direction: sub(b, a) },
        AxisBy::PointDirection { point, direction } => Datum::Axis { point, direction },
        AxisBy::DatumPoints { a, b } => {
            let a = datum_point_at(m, &a)?;
            Datum::Axis { point: a, direction: sub(datum_point_at(m, &b)?, a) }
        }
    })
}

/// Most mesh nodes a vertex is looked for among.
const MAX_VERTEX_NODES: usize = 10_000_000;

fn point_datum(s: &Session, by: PointBy) -> Result<Datum> {
    let shape = |body: &str| -> Result<&fieldcraft_geometry::Shape> {
        let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
        design.body(body).ok_or_else(|| EngineError::Params(format!("the design has no body `{body}`")))
    };
    let point = match by {
        PointBy::Coordinates { point } => point,
        PointBy::Vertex { body, vertex } => {
            // A vertex of the open design by its persistent name, else one the mesh tags.
            let design_vertex = match s.design.as_deref().and_then(|d| d.body(&body)) {
                Some(shape) => shape.vertices()?.into_iter().find(|v| v.name == vertex).map(|v| v.point),
                None => None,
            };
            match design_vertex {
                Some(point) => point,
                None => {
                    let m = &s.model;
                    let at = fieldcraft_model::GeometryRef::new(body.clone(), fieldcraft_model::Entity::Vertex(vertex.clone()))?;
                    let tags = m.classification().tags_selected_by(&at);
                    let node = m
                        .classification()
                        .nodes()
                        .iter()
                        .take(MAX_VERTEX_NODES)
                        .find(|(_, t)| tags.contains(t))
                        .and_then(|(n, _)| m.node(*n))
                        .ok_or_else(|| EngineError::Params(format!("body `{body}` has no vertex `{vertex}` (geometry.bodies lists its vertices)")))?;
                    node.xyz
                }
            }
        }
        PointBy::Edge { body, edge, fraction } => {
            if !(0.0..=1.0).contains(&fraction) {
                return Err(EngineError::invalid("fraction", "is a number from 0 to 1".to_string()));
            }
            let edges = shape(&body)?.edges()?;
            let e = edges.iter().find(|e| e.name == edge).ok_or_else(|| EngineError::Params(format!("body `{body}` has no edge `{edge}`")))?;
            e.point_at(fraction).ok_or_else(|| EngineError::Params(format!("edge `{edge}` of body `{body}` has no points")))?
        }
        PointBy::FaceCentre { body, face } => shape(&body)?.face_centre(&face)?,
        PointBy::AxisPlane { axis, plane } => {
            let (p, d) = datum_axis_of(&s.model, &axis)?;
            let (q, n) = datum_plane_of(&s.model, &plane)?;
            let along = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            let across = along(d, n);
            // Unit vectors: the sine of the angle between the axis and the plane.
            if across.abs() < MIN_SINE {
                return Err(EngineError::Params(format!("axis `{axis}` is parallel to plane `{plane}`: they do not cross")));
            }
            let t = along(sub(q, p), n) / across;
            [p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2]]
        }
    };
    Ok(Datum::Point { point })
}

/// The datum `by` defines for a datum of `kind` (`plane`, `axis` or `point`).
fn datum_of(s: &Session, kind: &str, p: &Value) -> Result<Datum> {
    match kind {
        "plane" => plane_datum(&s.model, super::req(p, "by")?),
        "axis" => axis_datum(&s.model, super::req(p, "by")?),
        _ => point_datum(s, super::req(p, "by")?),
    }
}

/// The datum a datum form's `params` define so far, for drawing it before OK; `None` when
/// `command` makes no datum or the parameters do not define one yet. The model is not changed.
pub fn preview(s: &Session, command: &str, params: &Value) -> Option<Datum> {
    let kind = match command {
        "geometry.datum" | "geometry.datum.edit" => "plane",
        "geometry.datum.axis" | "geometry.datum.axis.edit" => "axis",
        "geometry.datum.point" | "geometry.datum.point.edit" => "point",
        _ => return None,
    };
    // Datums are made in millimetres, with the geometry; the session stores them in SI.
    let framed = crate::units::geometry_session(s).ok()?;
    let params = crate::units::geometry_params(s, command, params).ok()?;
    let mut d = datum_of(&framed, kind, &params).ok()?;
    match &mut d {
        Datum::Plane { point, .. } | Datum::Axis { point, .. } | Datum::Point { point } => *point = crate::units::geometry_point(*point),
    }
    Some(d)
}

/// The parameters of the command that edits datum `name`, describing it as it is (a plane by its
/// point and normal, an axis by its point and direction, a point by its coordinates), in SI as
/// the model stores it ([`crate::units::shown_params`] shows them in the display system); and
/// that command.
pub fn edit_params(name: &str, d: &Datum) -> (&'static str, Value) {
    match *d {
        Datum::Plane { point, normal } => (DATUM_EDIT.id, json!({"name": name, "by": {"type": "point_normal", "point": point, "normal": normal}})),
        Datum::Axis { point, direction } => {
            (DATUM_AXIS_EDIT.id, json!({"name": name, "by": {"type": "point_direction", "point": point, "direction": direction}}))
        }
        Datum::Point { point } => (DATUM_POINT_EDIT.id, json!({"name": name, "by": {"type": "coordinates", "point": point}})),
    }
}

fn datum_json(name: &str, d: &Datum) -> Value {
    let mut v = serde_json::to_value(d).unwrap_or(Value::Null);
    if let Some(o) = v.as_object_mut() {
        o.insert("name".into(), json!(name));
    }
    v
}

fn datum_create(s: &mut Session, p: &Value, kind: &str) -> Result<Value> {
    let name = super::name(p)?;
    let datum = datum_of(s, kind, p)?;
    s.model.add_datum(name, datum)?;
    Ok(datum_json(name, s.model.datum(name).unwrap_or(&datum)))
}

fn datum_plane(s: &mut Session, p: &Value) -> Result<Value> {
    datum_create(s, p, "plane")
}

fn datum_axis(s: &mut Session, p: &Value) -> Result<Value> {
    datum_create(s, p, "axis")
}

fn datum_point(s: &mut Session, p: &Value) -> Result<Value> {
    datum_create(s, p, "point")
}

/// Redefine datum `name` of `kind`; the partitions that cut by a plane cut again.
fn datum_edit(s: &mut Session, p: &Value, kind: &str) -> Result<Value> {
    let name = super::name(p)?;
    match s.model.datum(name) {
        Some(d) if d.kind() == kind => {}
        Some(d) => return Err(EngineError::Params(format!("datum `{name}` is a datum {}, not a {kind}", d.kind()))),
        None => return Err(EngineError::Params(format!("no datum {kind} `{name}`"))),
    }
    let datum = datum_of(s, kind, p)?;
    s.model.replace_datum(name, datum)?;
    let cleared = if s.model.datum_users(name).is_empty() { super::replace::Cleared::default() } else { recut(s)? };
    let reply = datum_json(name, s.model.datum(name).unwrap_or(&datum));
    Ok(cleared.report_if_any(s, reply))
}

fn datum_delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    let cascade: bool = super::opt(p, "cascade")?.unwrap_or(false);
    let partitions: Vec<String> = s.model.partitions().iter().filter(|x| x.datum == name).map(|x| x.name.clone()).collect();
    if cascade {
        for partition in &partitions {
            s.model.remove_partition(partition)?;
        }
    }
    s.model.remove_datum(name)?;
    let cleared = if cascade && !partitions.is_empty() { recut(s)? } else { super::replace::Cleared::default() };
    let reply = json!({"deleted": name, "partitions": if cascade { partitions } else { Vec::new() }});
    Ok(cleared.report_if_any(s, reply))
}

/// The model's partitions as cuts of the geometry.
fn cuts(model: &Model) -> Vec<Cut> {
    model
        .partitions()
        .iter()
        .filter(|p| model.partition_cuts(p))
        .filter_map(|p| {
            let (point, normal) = model.datum(&p.datum)?.plane()?;
            Some(Cut { name: p.name.clone(), point, normal, bodies: p.bodies.clone() })
        })
        .collect()
}

/// `design` without the features `model` defeatures, cut by its partitions, then with its face
/// merges (after those of the features defeatured virtually) and its face splits (those along a
/// suppressed datum left out; the design as opened when the model has none).
pub(crate) fn partitioned(design: &Design, model: &Model) -> Result<Design> {
    let defeatures: Vec<fieldcraft_geometry::defeature::Defeature> =
        model.defeatures().iter().map(|d| fieldcraft_geometry::defeature::Defeature { name: d.name.clone(), features: d.features.clone() }).collect();
    let (design, mut merges) = design.defeatured(&defeatures)?;
    merges.extend(model.face_merges().iter().map(|m| Merge {
        name: m.name.clone(),
        body: m.body.clone(),
        faces: m.faces.clone(),
        into: m.into.clone(),
    }));
    let splits: Vec<fieldcraft_geometry::face_split::FaceSplit> = model
        .face_splits()
        .iter()
        .filter(|s| !model.datum_suppressed(&s.datum))
        .filter_map(|s| {
            let (point, normal) = model.datum(&s.datum)?.plane()?;
            Some(fieldcraft_geometry::face_split::FaceSplit { name: s.name.clone(), body: s.body.clone(), face: s.face.clone(), point, normal })
        })
        .collect();
    Ok(design.partitioned(&cuts(model))?.merged(&merges)?.split(&splits)?)
}

/// Cut the session's design by its model's partitions again (and merge its faces again). The bodies that are gone (cells
/// whose partition went, a body a new partition cut) take their meshes with them (#385): what
/// went is in the returned [`super::replace::Cleared`], which the command reports.
pub(super) fn recut(s: &mut Session) -> Result<super::replace::Cleared> {
    let Some(d) = &s.design else { return Ok(super::replace::Cleared::default()) };
    let cut = partitioned(d, &s.model)?;
    let gone: Vec<String> = d.bodies.iter().filter(|b| cut.body(&b.name).is_none()).map(|b| b.name.clone()).collect();
    s.design = Some(Arc::new(cut));
    Ok(super::replace::clear(&mut s.model, gone.iter().map(String::as_str)))
}

fn partition(s: &mut Session, p: &Value) -> Result<Value> {
    let display = s.units().clone();
    let name = super::name(p)?.to_string();
    let datum: String = super::req(p, "datum")?;
    let bodies: Vec<String> = super::opt(p, "bodies")?.unwrap_or_default();
    datum_plane_of(&s.model, &datum)?;
    s.model.add_partition(Partition { name: name.clone(), datum: datum.clone(), bodies, suppressed: false })?;
    let cleared = recut(s)?;
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let feature = format!("partition {name}");
    let errors: Vec<&str> = design.errors.iter().filter(|e| e.feature == feature).map(|e| e.message.as_str()).collect();
    if !errors.is_empty() {
        return Err(EngineError::Params(format!("partition `{name}`: {}", errors.join("; "))));
    }
    let report = design.cuts.iter().find(|c| c.name == name);
    let split = report.map(|r| r.split.clone()).unwrap_or_default();
    if split.is_empty() {
        return Err(EngineError::Params(format!("partition `{name}`: the plane divides no body of the design")));
    }
    let cells: Vec<String> = design.bodies.iter().map(|b| b.name.clone()).filter(|b| split.iter().any(|x| b.starts_with(&format!("{x}/")))).collect();
    let (faces, bodies): (usize, Vec<String>) = (report.map_or(0, |r| r.faces), design.bodies.iter().map(|b| b.name.clone()).collect());
    let mut warnings: Vec<String> = report.map(|r| r.warnings.clone()).unwrap_or_default();
    fieldcraft_model::progress::stage("Partition: checking the cells", Some(cells.len() as f64), "cells");
    for (k, cell) in cells.iter().enumerate() {
        fieldcraft_model::progress::count(k as f64);
        fieldcraft_model::progress::check()?;
        if let Some(why) = design.body(cell).and_then(|shape| unmeshable(shape, None, &display).err()) {
            warnings.push(format!(
                "partition `{name}`: cell `{cell}` has a boundary no mesher takes ({why}): the plane passes within the kernel's tolerance of a face, edge or vertex it cut (along a face, through the axis of a hole); move the plane off it (by a fraction of the element size)"
            ));
        }
    }
    // The cut shows as edges on the cells, so the plane steps out of the way (view state, #295):
    // selecting it in the browser or its eye shows it again.
    if s.view.hidden.len() < MAX_HIDDEN {
        s.view.hidden.insert(Item { kind: ItemKind::Datum, name: datum });
    }
    let mut reply = json!({
        "name": name,
        "split": split,
        "cells": cells,
        "cut_faces": faces,
        "bodies": bodies,
    });
    if !warnings.is_empty()
        && let Some(o) = reply.as_object_mut()
    {
        o.insert("warnings".into(), json!(warnings));
    }
    Ok(cleared.report_if_any(s, reply))
}

/// Whether a cell's boundary, tessellated as the meshers take it (by default, or at `tolerance`), bounds a volume mesh:
/// the defects (open edges, edges on three triangles) when it does not.
/// Why `shape`'s boundary bounds no volume mesh, its numbers shown in `display`.
pub(super) fn unmeshable(shape: &fieldcraft_geometry::Shape, tolerance: Option<f64>, display: &fieldcraft_model::DisplaySystem) -> Result<()> {
    let input = super::hex::shape_input(shape, tolerance)?;
    let options = fieldcraft_mesh::surface::SurfaceOptions { weld_distance: input.mesh.weld_distance, ..Default::default() };
    fieldcraft_mesh::surface::BoundarySurface::new(&input.mesh.positions, &input.mesh.triangles, &input.mesh.tri_face, &options)
        .map(|_| ())
        .map_err(|e| EngineError::Mesh(e).shown(crate::units::geometry(), display))
}

fn datum_suppress(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    let suppressed = super::opt(p, "suppressed")?.unwrap_or(true);
    s.model.set_datum_suppressed(name, suppressed)?;
    let cleared = recut(s)?;
    let bodies: Vec<String> = s.design.iter().flat_map(|d| d.bodies.iter().map(|b| b.name.clone())).collect();
    Ok(cleared.report_if_any(s, json!({"name": name, "suppressed": suppressed, "bodies": bodies})))
}

fn partition_suppress(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    let suppressed = super::opt(p, "suppressed")?.unwrap_or(true);
    s.model.set_partition_suppressed(name, suppressed)?;
    let cleared = recut(s)?;
    let bodies: Vec<String> = s.design.iter().flat_map(|d| d.bodies.iter().map(|b| b.name.clone())).collect();
    Ok(cleared.report_if_any(s, json!({"name": name, "suppressed": suppressed, "bodies": bodies})))
}

fn partition_delete(s: &mut Session, p: &Value) -> Result<Value> {
    let name = super::name(p)?;
    s.model.remove_partition(name)?;
    let cleared = recut(s)?;
    let bodies: Vec<String> = s.design.iter().flat_map(|d| d.bodies.iter().map(|b| b.name.clone())).collect();
    Ok(cleared.report_if_any(s, json!({"deleted": name, "bodies": bodies})))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A session with the sample plate (40 × 30 × 10, Ø10 hole at (15, 15)) open.
    pub(crate) fn plate() -> (Session, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("fieldcraft-partition-{}-{:?}", std::process::id(), std::thread::current().id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("plate.solidcraft");
        let plate = fieldcraft_geometry::sample::Plate::default();
        std::fs::write(&path, fieldcraft_geometry::sample::plate_design(&plate).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path.to_str().unwrap()})).unwrap();
        (s, dir)
    }

    fn bodies(s: &Session) -> Vec<String> {
        s.design.as_ref().unwrap().bodies.iter().map(|b| b.name.clone()).collect()
    }

    #[test]
    fn datums_by_each_definition() {
        let mut s = Session::default();
        let r = s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 15}})).unwrap();
        assert_eq!(r, json!({"name": "Mid", "type": "plane", "point": [15.0, 0.0, 0.0], "normal": [1.0, 0.0, 0.0], "units": {"length": "mm"}}));
        let r = s
            .execute("geometry.datum", &json!({"name": "Tri", "by": {"type": "three_points", "a": [0, 0, 1], "b": [1, 0, 1], "c": [0, 1, 1]}}))
            .unwrap();
        assert_eq!(r["normal"], json!([0.0, 0.0, 1.0]));
        s.execute("geometry.datum", &json!({"name": "Tilt", "by": {"type": "point_normal", "point": [0, 0, 0], "normal": [0, 3, 4]}})).unwrap();
        assert_eq!(s.model.datum("Tilt"), Some(&Datum::Plane { point: [0.0; 3], normal: [0.0, 0.6, 0.8] }));
        let r = s.execute("geometry.datum.axis", &json!({"name": "Z", "by": {"type": "principal", "axis": "z", "through": [1, 2, 0]}})).unwrap();
        assert_eq!(r, json!({"name": "Z", "type": "axis", "point": [1.0, 2.0, 0.0], "direction": [0.0, 0.0, 1.0], "units": {"length": "mm"}}));
        s.execute("geometry.datum.axis", &json!({"name": "D", "by": {"type": "two_points", "a": [0, 0, 0], "b": [0, 2, 0]}})).unwrap();
        s.execute("geometry.datum.axis", &json!({"name": "E", "by": {"type": "point_direction", "point": [0, 0, 0], "direction": [5, 0, 0]}}))
            .unwrap();
        assert_eq!(s.model.datums().len(), 6);
        s.execute("geometry.datum.delete", &json!({"name": "D"})).unwrap();
        assert_eq!(s.model.datums().len(), 5);
    }

    #[test]
    fn hostile_datum_params_change_nothing() {
        let mut s = Session::default();
        for p in [
            json!({}),
            json!({"name": "A"}),
            json!({"name": "A", "by": {"type": "offset", "plane": "ab"}}),
            json!({"name": "A", "by": {"type": "offset", "plane": "xy", "offset": "x"}}),
            json!({"name": "A", "by": {"type": "three_points", "a": [0, 0, 0], "b": [1, 1, 1], "c": [2, 2, 2]}}),
            json!({"name": "A", "by": {"type": "point_normal", "point": [0, 0, 0], "normal": [0, 0, 0]}}),
            json!({"name": "A", "by": {"type": "point_normal", "point": [1e300, 0, 0], "normal": [0, 0, 1]}}),
            json!({"name": "A", "by": {"type": "point_normal", "point": [0, 0, 0], "normal": [0, 0, 1], "extra": 1}}),
            json!({"name": "1A", "by": {"type": "offset", "plane": "xy"}}),
            json!({"name": "A", "by": {"type": "nope"}}),
        ] {
            assert!(s.execute("geometry.datum", &p).is_err(), "{p}");
        }
        for p in [
            json!({"name": "A", "by": {"type": "two_points", "a": [1, 1, 1], "b": [1, 1, 1]}}),
            json!({"name": "A", "by": {"type": "principal", "axis": "w"}}),
        ] {
            assert!(s.execute("geometry.datum.axis", &p).is_err(), "{p}");
        }
        assert!(s.model.datums().is_empty());
        assert!(s.journal.is_empty());
    }

    /// Suppressing a partition merges its cells back; suppressing its datum does too and stops
    /// the datum being an input; unsuppressing cuts again; each is one undo step (#291).
    #[test]
    fn suppressed_partitions_and_datums_leave_the_design_uncut() {
        let (mut s, dir) = plate();
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 15}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
        let count = |s: &Session| s.design.as_ref().map_or(0, |d| d.bodies.len());
        assert_eq!(count(&s), 2);
        let r = s.execute("geometry.partition.suppress", &json!({"name": "Half"})).unwrap();
        assert_eq!((r["suppressed"].as_bool(), count(&s)), (Some(true), 1));
        s.execute("geometry.partition.suppress", &json!({"name": "Half", "suppressed": false})).unwrap();
        assert_eq!(count(&s), 2);
        s.execute("geometry.datum.suppress", &json!({"name": "Mid"})).unwrap();
        assert_eq!(count(&s), 1);
        let e = s.execute("geometry.partition", &json!({"name": "Again", "datum": "Mid"})).unwrap_err().to_string();
        assert!(e.contains("suppressed"), "{e}");
        let e = s.execute("geometry.datum.point", &json!({"name": "P", "by": {"type": "axis_plane", "axis": "Mid", "plane": "Mid"}}));
        assert!(e.is_err());
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(count(&s), 2, "undo cuts again");
        assert!(s.execute("geometry.datum.suppress", &json!({"name": "Nope"})).is_err());
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn a_partition_cuts_the_design_and_is_replayed_when_it_opens_again() {
        let (mut s, dir) = plate();
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 15}})).unwrap();
        let r = s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
        let body = r["split"][0].as_str().unwrap().to_string();
        assert_eq!(r["cells"], json!([format!("{body}/a"), format!("{body}/b")]));
        assert_eq!(r["cut_faces"], json!(4));
        // The datum is in use, and hidden: the cut shows as edges on the cells (#295). Its eye
        // shows it again.
        assert!(s.execute("geometry.datum.delete", &json!({"name": "Mid"})).is_err());
        assert!(!s.shows(ItemKind::Datum, "Mid"));
        s.execute("view.show", &json!({"items": [{"kind": "datum", "name": "Mid"}]})).unwrap();
        assert!(s.shows(ItemKind::Datum, "Mid"));
        s.execute("view.hide", &json!({"items": [{"kind": "datum", "name": "Mid"}]})).unwrap();

        // A set on a cut face selects the cut face's pieces.
        let r = s.execute("geometry.bodies", &json!({})).unwrap();
        let faces: Vec<&str> = r[0]["faces"].as_array().unwrap().iter().map(|f| f["name"].as_str().unwrap()).collect();
        assert!(faces.contains(&"partition:Half:a#0"), "{faces:?}");

        // Saved and opened again: the design is cut again by the model's partition.
        let file = dir.join("plate.fieldcraft");
        s.execute("file.save_as", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut t = Session::default();
        t.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(bodies(&t), bodies(&s));
        assert!(!t.shows(ItemKind::Datum, "Mid"), "still hidden when opened again");
        // geometry.open again keeps the model, so its partition.
        let design = dir.join("plate.solidcraft");
        let r = t.execute("geometry.open", &json!({"path": design.to_str().unwrap()})).unwrap();
        assert_eq!(r["bodies"], json!(bodies(&s)));
        // A new model drops the partitions with the rest of the model.
        let r = t.execute("geometry.open", &json!({"path": design.to_str().unwrap(), "model": "new"})).unwrap();
        assert_eq!(r["bodies"], json!([body]));

        // Undo brings the whole design back; deleting the partition does too.
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(bodies(&s), std::slice::from_ref(&body));
        s.execute("edit.redo", &json!({})).unwrap();
        assert_eq!(bodies(&s).len(), 2);
        s.execute("geometry.partition.delete", &json!({"name": "Half"})).unwrap();
        assert_eq!(bodies(&s), [body]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn datum_points_by_each_definition_and_planes_and_axes_through_them() {
        let (mut s, dir) = plate();
        // Datums are stored in metres; the plate is drawn in millimetres.
        let stored = |s: &Session, n: &str| s.model.datum(n).and_then(Datum::point).unwrap();
        let point = |s: &Session, n: &str| stored(s, n).map(crate::units::to_geometry_length);
        let close = |a: [f64; 3], b: [f64; 3]| a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-6);
        s.execute("geometry.datum.point", &json!({"name": "O", "by": {"type": "coordinates", "point": [0, 0, 0]}})).unwrap();
        s.execute("geometry.datum.point", &json!({"name": "X", "by": {"type": "coordinates", "point": [10, 0, 0]}})).unwrap();
        s.execute("geometry.datum.point", &json!({"name": "Y", "by": {"type": "coordinates", "point": [0, 10, 0]}})).unwrap();
        // The plate's edge along x at y = 0, z = 0, its midpoint and a quarter along it.
        let body = s.design.as_ref().unwrap().bodies[0].clone();
        let edge = body.edges().unwrap().into_iter().find(|e| e.points.iter().all(|p| p[1].abs() < 1e-9 && p[2].abs() < 1e-9)).unwrap();
        let by = |fraction: f64| json!({"type": "edge", "body": body.name, "edge": edge.name, "fraction": fraction});
        s.execute("geometry.datum.point", &json!({"name": "Mid", "by": {"type": "edge", "body": body.name, "edge": edge.name}})).unwrap();
        assert!(close(point(&s, "Mid"), [20.0, 0.0, 0.0]), "{:?}", point(&s, "Mid"));
        s.execute("geometry.datum.point", &json!({"name": "Q", "by": by(0.25)})).unwrap();
        let q = point(&s, "Q");
        assert!((q[0] - 10.0).abs() < 1e-6 || (q[0] - 30.0).abs() < 1e-6, "{q:?}");
        // The centre of the face at x = 0.
        let side =
            body.face_surfaces().unwrap().iter().position(
                |f| matches!(f, Some(fieldcraft_geometry::FaceSurface::Plane { normal, point }) if normal[0] < -0.5 && point[0].abs() < 1e-9),
            );
        let face = body.faces()[side.unwrap()].name.clone();
        s.execute("geometry.datum.point", &json!({"name": "C", "by": {"type": "face_centre", "body": body.name, "face": face}})).unwrap();
        assert!(close(point(&s, "C"), [0.0, 15.0, 5.0]), "{:?}", point(&s, "C"));

        // Planes and axes through datum points; a point where an axis crosses a plane.
        let r = s.execute("geometry.datum", &json!({"name": "Base", "by": {"type": "datum_points", "a": "O", "b": "X", "c": "Y"}})).unwrap();
        assert_eq!(r["normal"], json!([0.0, 0.0, 1.0]));
        s.execute("geometry.datum.axis", &json!({"name": "Ax", "by": {"type": "datum_points", "a": "O", "b": "X"}})).unwrap();
        assert_eq!(s.model.datum("Ax"), Some(&Datum::Axis { point: [0.0; 3], direction: [1.0, 0.0, 0.0] }));
        s.execute("geometry.datum", &json!({"name": "Cut", "by": {"type": "normal_to_axis", "axis": "Ax", "point": "Mid"}})).unwrap();
        assert_eq!(s.model.datum("Cut"), Some(&Datum::Plane { point: stored(&s, "Mid"), normal: [1.0, 0.0, 0.0] }));
        s.execute("geometry.datum.point", &json!({"name": "Hit", "by": {"type": "axis_plane", "axis": "Ax", "plane": "Cut"}})).unwrap();
        assert!(close(point(&s, "Hit"), [20.0, 0.0, 0.0]));

        let journal = s.journal.len();
        for p in [
            json!({"name": "Z", "by": {"type": "axis_plane", "axis": "Ax", "plane": "Base"}}),
            json!({"name": "Z", "by": {"type": "axis_plane", "axis": "Base", "plane": "Cut"}}),
            json!({"name": "Z", "by": {"type": "edge", "body": body.name, "edge": edge.name, "fraction": 2}}),
            json!({"name": "Z", "by": {"type": "edge", "body": body.name, "edge": "nope"}}),
            json!({"name": "Z", "by": {"type": "edge", "body": "Nope", "edge": edge.name}}),
            json!({"name": "Z", "by": {"type": "face_centre", "body": body.name, "face": "nope"}}),
            // No mesh: no vertex is named.
            json!({"name": "Z", "by": {"type": "vertex", "body": body.name, "vertex": "v1"}}),
            json!({"name": "Z", "by": {"type": "coordinates", "point": [f64::MAX, 0, 0]}}),
            json!({"name": "Z", "by": {"type": "coordinates", "point": [0, 0]}}),
        ] {
            assert!(s.execute("geometry.datum.point", &p).is_err(), "{p}");
        }
        for (command, p) in [
            ("geometry.datum", json!({"name": "Z", "by": {"type": "datum_points", "a": "O", "b": "X", "c": "Mid"}})),
            ("geometry.datum", json!({"name": "Z", "by": {"type": "datum_points", "a": "O", "b": "X", "c": "Ax"}})),
            ("geometry.datum", json!({"name": "Z", "by": {"type": "normal_to_axis", "axis": "Nope"}})),
            ("geometry.datum.axis", json!({"name": "Z", "by": {"type": "datum_points", "a": "O", "b": "O"}})),
        ] {
            assert!(s.execute(command, &p).is_err(), "{p}");
        }
        assert_eq!(s.journal.len(), journal);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_vertex_point_is_at_the_tagged_mesh_node() {
        let mut s = Session::default();
        s.execute(
            "mesh.block",
            &json!({"body": "Bar", "kind": "hex8", "corners": [[0,0,0],[4,0,0],[4,1,0],[0,1,0],[0,0,1],[4,0,1],[4,1,1],[0,1,1]], "divisions": [2, 1, 1], "vertices": {"7": "tip"}}),
        )
        .unwrap();
        s.execute("geometry.datum.point", &json!({"name": "Tip", "by": {"type": "vertex", "body": "Bar", "vertex": "tip"}})).unwrap();
        assert_eq!(s.model.datum("Tip").and_then(Datum::point), Some([0.004, 0.001, 0.001]), "4, 1, 1 mm, stored in metres");
    }

    #[test]
    fn editing_and_deleting_datums_follow_their_partitions() {
        let (mut s, dir) = plate();
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 15}})).unwrap();
        s.execute("geometry.datum.axis", &json!({"name": "Ax", "by": {"type": "principal", "axis": "x"}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "Half", "datum": "Mid"})).unwrap();
        let cut_at = |s: &Session| s.design.as_ref().unwrap().bodies[0].measure().unwrap().bbox_max[0];
        assert!((cut_at(&s) - 15.0).abs() < 1e-6);

        // Moving the plane cuts again; an edit of the wrong kind or to an axis is refused.
        s.execute("geometry.datum.edit", &json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 25}})).unwrap();
        assert!((cut_at(&s) - 25.0).abs() < 1e-6);
        assert!(s.execute("geometry.datum.axis.edit", &json!({"name": "Mid", "by": {"type": "principal", "axis": "y"}})).is_err());
        assert!(s.execute("geometry.datum.edit", &json!({"name": "Ax", "by": {"type": "offset", "plane": "xy"}})).is_err());
        s.execute("geometry.datum.axis.edit", &json!({"name": "Ax", "by": {"type": "principal", "axis": "y"}})).unwrap();
        assert_eq!(s.model.datum("Ax"), Some(&Datum::Axis { point: [0.0; 3], direction: [0.0, 1.0, 0.0] }));

        // The edit commands take the parameters edit_params gives back, shown in the display
        // system, and they change nothing.
        for name in ["Mid", "Ax"] {
            let before = *s.model.datum(name).unwrap();
            let (command, params) = edit_params(name, &before);
            let shown = crate::units::shown_params(command, &params, s.units()).unwrap();
            s.execute(command, &shown).unwrap();
            assert_eq!(s.model.datum(name), Some(&before));
        }

        // Deleting the plane names its partition, or deletes it too.
        let e = s.execute("geometry.datum.delete", &json!({"name": "Mid"})).unwrap_err().to_string();
        assert!(e.contains("partition `Half`"), "{e}");
        let r = s.execute("geometry.datum.delete", &json!({"name": "Mid", "cascade": true})).unwrap();
        assert_eq!(r["partitions"], json!(["Half"]));
        assert!(s.model.partitions().is_empty() && s.model.datum("Mid").is_none());
        assert_eq!(bodies(&s).len(), 1);
        // One undo step brings both back.
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!((s.model.partitions().len(), bodies(&s).len()), (1, 2));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_partition_that_cuts_nothing_changes_nothing() {
        let (mut s, dir) = plate();
        s.execute("geometry.datum", &json!({"name": "Out", "by": {"type": "offset", "plane": "yz", "offset": 100}})).unwrap();
        s.execute("geometry.datum.axis", &json!({"name": "Ax", "by": {"type": "principal", "axis": "x"}})).unwrap();
        let journal = s.journal.len();
        for p in [
            json!({"name": "P", "datum": "Out"}),
            json!({"name": "P", "datum": "Ax"}),
            json!({"name": "P", "datum": "Nope"}),
            json!({"name": "P", "datum": "Out", "bodies": ["Nope"]}),
            json!({"name": "P", "datum": "Out", "bodies": "x"}),
            json!({"name": "", "datum": "Out"}),
        ] {
            assert!(s.execute("geometry.partition", &p).is_err(), "{p}");
        }
        assert!(s.model.partitions().is_empty());
        assert_eq!(bodies(&s).len(), 1);
        assert_eq!(s.journal.len(), journal);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A session with the design `commands` (a SolidCraft script) open.
    pub(crate) fn scripted(commands: Value) -> (crate::tests::TempDir, Session) {
        let bytes = fieldcraft_geometry::script_design(&json!({"commands": commands}).to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    /// Partition the open design by the planes `x = …` / `y = …` / `z = …` (`(axis, offset)`).
    pub(crate) fn cut_at(s: &mut Session, planes: &[(usize, f64)]) {
        for (i, &(k, offset)) in planes.iter().enumerate() {
            let (name, plane) = (format!("C{i}"), ["yz", "zx", "xy"][k]);
            s.execute("geometry.datum", &json!({"name": name, "by": {"type": "offset", "plane": plane, "offset": offset}})).unwrap();
            s.execute("geometry.partition", &json!({"name": format!("P{i}"), "datum": name})).unwrap();
        }
    }

    /// A box with a plate-shaped void in its middle that touches no side (#388): the planes
    /// round the void leave a centre column whose material lies above and below the void, two
    /// pieces that become separate cells only at the last cut.
    #[test]
    fn pieces_that_do_not_touch_are_separate_cells() {
        let (_dir, mut s) = scripted(json!([
            {"command": "solid.box", "params": {"length": 106, "width": 76, "height": 20, "corner": [0, 0, 0], "body_name": "Air"}},
            {"command": "solid.box", "params": {"length": 85, "width": 56, "height": 2, "corner": [13, 13, 9], "operation": "cut"}},
        ]));
        let whole = s.design.as_ref().unwrap().bodies[0].measure().unwrap().volume;
        cut_at(&mut s, &[(0, 22.5), (0, 68.5), (1, 22.5), (1, 59.5)]);
        let design = s.design.clone().unwrap();
        let shells: Vec<(&str, usize)> = design.bodies.iter().map(|b| (b.name.as_str(), b.shells().unwrap())).collect();
        assert!(shells.iter().all(|s| s.1 == 1), "{shells:?}");
        // 3 × 3 columns, the centre one in two pieces.
        assert_eq!(shells.len(), 10, "{shells:?}");
        let volume: f64 = design.bodies.iter().map(|b| b.measure().unwrap().volume).sum();
        assert!((volume - whole).abs() < 1e-6 * whole, "{volume} vs {whole}");
        let centre: Vec<&str> = shells.iter().map(|s| s.0).filter(|n| n.ends_with('1') || n.ends_with('2')).collect();
        assert_eq!(centre.len(), 2, "{centre:?}");
        // Numbered by centroid: the piece below the void first.
        let z = |n: &str| design.body(n).unwrap().measure().unwrap().centroid[2];
        assert!(z(centre[0]) < 9.0 && z(centre[1]) > 11.0, "{centre:?}");
    }

    /// #385: the cells' meshes go with the cells (deleting, suppressing or replacing their
    /// partitions), a body's mesh goes when a partition cuts it, and meshes left over from an
    /// earlier cut are listed by model.check and removed by mesh.clear.
    #[test]
    fn meshes_go_with_the_bodies_partitions_remove() {
        let (mut s, dir) = plate();
        let whole = s.design.as_ref().unwrap().bodies[0].measure().unwrap().volume;
        let volume = |s: &mut Session| s.execute("mesh.quality", &json!({})).unwrap()["volume"].as_f64().unwrap();
        // Partitioning a meshed body takes its mesh.
        s.execute("mesh.hex.generate", &json!({"size": 2.5})).unwrap();
        s.execute("geometry.datum", &json!({"name": "X30", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
        let r = s.execute("geometry.partition", &json!({"name": "P", "datum": "X30"})).unwrap();
        assert_eq!(r["replaced"][0]["body"], "Body1", "{r}");
        assert!(s.model.elements().is_empty());
        // Deleting the datum with its partition takes the cells' meshes; the whole body's new
        // mesh is the only one.
        s.execute("mesh.cells", &json!({"default": "hex", "size": 2.5})).unwrap();
        let r = s.execute("geometry.datum.delete", &json!({"name": "X30", "cascade": true})).unwrap();
        let gone: Vec<&str> = r["replaced"].as_array().unwrap().iter().map(|b| b["body"].as_str().unwrap()).collect();
        assert_eq!(gone, ["Body1/a", "Body1/b"], "{r}");
        assert!(s.model.elements().is_empty());
        s.execute("mesh.hex.generate", &json!({"size": 2.5})).unwrap();
        let v = volume(&mut s);
        assert!((v - whole).abs() < 1e-2 * whole, "{v} vs {whole}");
        // Suppressing or deleting a partition does the same.
        s.execute("geometry.datum", &json!({"name": "X30", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "P", "datum": "X30"})).unwrap();
        s.execute("mesh.cells", &json!({"default": "hex", "size": 2.5})).unwrap();
        let r = s.execute("geometry.partition.suppress", &json!({"name": "P"})).unwrap();
        assert_eq!(r["replaced"].as_array().map(Vec::len), Some(2), "{r}");
        s.execute("edit.undo", &json!({})).unwrap();
        assert!(!s.model.elements().is_empty(), "undo brings the cells' meshes back");
        let r = s.execute("geometry.partition.delete", &json!({"name": "P"})).unwrap();
        assert_eq!(r["replaced"].as_array().map(Vec::len), Some(2), "{r}");
        // A mesh left from an earlier cut (a model saved before the cut changed): listed and
        // cleared.
        s.execute("geometry.partition", &json!({"name": "P", "datum": "X30"})).unwrap();
        s.execute("mesh.cells", &json!({"default": "hex", "size": 2.5})).unwrap();
        s.model.remove_partition("P").unwrap();
        let d = s.design.clone().unwrap();
        s.design = Some(Arc::new(partitioned(&d, &s.model).unwrap()));
        let check = s.execute("model.check", &json!({})).unwrap();
        let listed: Vec<&str> =
            check["diagnostics"].as_array().unwrap().iter().filter_map(|d| d["subject"].as_str()).filter(|d| d.starts_with("body")).collect();
        assert_eq!(listed, ["body `Body1/a`", "body `Body1/b`"], "{check}");
        let r = s.execute("mesh.clear", &json!({})).unwrap();
        assert_eq!(r["replaced"].as_array().map(Vec::len), Some(2), "{r}");
        assert!(s.model.elements().is_empty());
        assert!(s.execute("mesh.clear", &json!({"bodies": ["Body1"]})).is_err(), "nothing to clear");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The Raspberry Pi 5 enclosure analysis design (base, lid, heat sink, board, chips and the
    /// air inside), its SolidCraft recipe.
    pub(crate) fn pi5() -> (crate::tests::TempDir, Session) {
        scripted(serde_json::from_str(include_str!("../../tests/designs/pi5-enclosure.json")).unwrap())
    }

    /// #384: cuts stay on their plane (the kernel's own split put them 1.8 µm off it on the air
    /// of the enclosure, so the cut faces were not even recognized), and a plane along an
    /// existing face or through a hole's axis is a warning naming the plane, the body and where,
    /// or refused where the kernel could only cut it off the plane.
    #[test]
    fn cuts_stay_on_the_plane_and_risky_ones_are_named() {
        let (_dir, mut s) = pi5();
        let plane = |s: &mut Session, name: &str, x: f64| {
            s.execute("geometry.datum", &json!({"name": name, "by": {"type": "offset", "plane": "yz", "offset": x}})).unwrap();
            s.execute("geometry.partition", &json!({"name": format!("P{name}"), "datum": name, "bodies": ["Air"]})).unwrap()
        };
        // Beside the board's mounting holes: clean cells, cut faces on the plane.
        let r = plane(&mut s, "X", 16.5);
        assert!(r.get("warnings").is_none(), "{r}");
        let design = s.design.clone().unwrap();
        for cell in ["Air/a", "Air/b"] {
            let shape = design.body(cell).unwrap();
            let cut: Vec<f64> = shape
                .faces()
                .iter()
                .zip(shape.face_surfaces().unwrap())
                .filter(|(f, _)| f.name.starts_with("partition:PX:"))
                .filter_map(|(_, srf)| match srf {
                    Some(fieldcraft_geometry::FaceSurface::Plane { point, .. }) => Some(point[0]),
                    _ => None,
                })
                .collect();
            assert!(!cut.is_empty() && cut.iter().all(|x| (x - 16.5).abs() < 1e-9), "{cell}: {cut:?}");
        }
        s.execute("edit.undo", &json!({})).unwrap();
        // Along the side of the heat sink's pedestal: the kernel only cuts with its box nudged
        // off the plane, so the cut is refused, saying why (#427).
        s.execute("geometry.datum", &json!({"name": "Side", "by": {"type": "offset", "plane": "yz", "offset": 34.0}})).unwrap();
        let e = s.execute("geometry.partition", &json!({"name": "PSide", "datum": "Side", "bodies": ["Air"]})).unwrap_err().to_string();
        assert!(e.contains("could not cut `Air`") && e.contains("only made off the plane"), "{e}");
        // Through a standoff's hole: the cell whose boundary breaks at the hole's edge is named.
        let r = plane(&mut s, "Axis", 74.5);
        let warnings = r["warnings"].to_string();
        assert!(warnings.contains("partition `PAxis`: cell `Air/") && warnings.contains("boundary no mesher takes"), "{r}");
    }
}
