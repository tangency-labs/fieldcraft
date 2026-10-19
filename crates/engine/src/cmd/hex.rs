//! Structured hexahedral meshing of the open design's bodies or of imported surfaces.

use crate::schema::ProblemKind;
use std::borrow::Cow;

use fieldcraft_geometry::{FaceSurface, Shape, SurfaceMesh};
use fieldcraft_mesh::hex::{BoundaryFace, Cell, HexMesh, HexMethod, HexOptions, MAX_CELLS, Surface, hex_mesh, hex_mesh_cells};
use fieldcraft_mesh::surface::BoundarySurface;
use fieldcraft_model::{GeometryRef, Model};
use serde_json::{Value, json};

use super::sizes::Sizing;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const HEX_GENERATE: CommandSpec = CommandSpec {
    id: "mesh.hex.generate",
    label: "Hex Mesh",
    params: r#"{"body": "Body1" (default: every body of the open design) or "surface": 0 (an imported surface), "size": 2.0 (element size, a length; default 1/25 of the bounding-box diagonal), "order": 1 (hex8) or 2 (hex20, mid-edge nodes on the geometry), "method": "auto" (default) | "swept" (any sweep) | "extrude" | "sweep" | "revolve" | "mapped", "tolerance": 0.05 (chordal tolerance of the tessellation, a length; default 1e-3 of the diagonal)} — replaces the meshes of the bodies it meshes (as mesh.tet.generate; other bodies keep theirs) with structured hexahedra: multi-level sweeps along a direction or around an axis (their cross-section cut into mappable blocks automatically: O-grids round holes, triangle primitives, four-sided blocks; nothing is added to the model), one-to-one sweeps through linking faces, or mapped blocks. Boundary element faces are tagged with their B-rep faces (named as in geometry.bodies; by index when a face has no name; `weld`, `weld_distance` and `feature_angle` as in mesh.tet.generate), so sets by geometry resolve. Edge seeds (mesh.seed.edge) set the counts and spacing of their edges: exact seeds are held, and seeds that cannot all hold (opposite sides of a mapped face, the layers between two levels) are an error naming them. Size controls (mesh.size.set) set the division counts and spacing of the other edges through the graded size field (a structured mesh carries a count across its block, so a control refines the whole layer or row it lies in). Reports per body the method, counts and quality; a body no structured method fits is an error that says why"#,
    form: Some(super::forms::HEX_GENERATE),
    kind: CommandKind::Model,
    requires: &[Requirement::Geometry],
    run: hex_generate,
};

/// Bodies at most a command meshes.
const MAX_BODIES: usize = 1000;
/// Smallest element size and chordal tolerance, relative to the bounding-box diagonal.
const MIN_RELATIVE_SIZE: f64 = 1e-4;
const MIN_RELATIVE_TOLERANCE: f64 = 1e-5;
/// Default element size and chordal tolerance, relative to the bounding-box diagonal.
const DEFAULT_RELATIVE_SIZE: f64 = 1.0 / 25.0;
const DEFAULT_RELATIVE_TOLERANCE: f64 = 1e-3;

/// One body to mesh (also the input of `mesh.hex.octree`).
pub(super) struct Input<'a> {
    pub(super) name: String,
    pub(super) mesh: Cow<'a, SurfaceMesh>,
    pub(super) surfaces: Vec<Option<Surface>>,
    pub(super) face_names: Vec<String>,
    pub(super) diagonal: f64,
    /// The design's edge seeds ([`super::seeds::edge_seeds`]); the mesher applies those on its
    /// edges.
    pub(super) seeds: Vec<fieldcraft_mesh::seeds::EdgeSeed>,
    /// Per triangle of `mesh`, the B-rep face it lies on: its piece of a merged face (empty for
    /// a surface that is not a body's).
    pub(super) pieces: Vec<u32>,
}

fn analytic(s: FaceSurface) -> Surface {
    match s {
        FaceSurface::Plane { normal, point } => Surface::Plane { point, normal },
        FaceSurface::Cylinder { origin, axis, radius } => Surface::Cylinder { origin, axis, radius },
        FaceSurface::Cone { apex, axis, half_angle } => Surface::Cone { apex, axis, half_angle },
        FaceSurface::Sphere { center, radius } => Surface::Sphere { center, radius },
        FaceSurface::Torus { center, axis, major, minor } => Surface::Torus { center, axis, major, minor },
    }
}

pub(super) fn number(p: &Value, key: &str) -> Result<Option<f64>> {
    match p.get(key) {
        None => Ok(None),
        Some(v) => v.as_f64().filter(|x| x.is_finite() && *x > 0.0).map(Some).ok_or_else(|| EngineError::param(key, ProblemKind::NotPositive)),
    }
}

fn diagonal(lo: [f64; 3], hi: [f64; 3]) -> f64 {
    (0..3).map(|i| (hi[i] - lo[i]).powi(2)).sum::<f64>().sqrt()
}

pub(super) fn shape_input(shape: &Shape, tolerance: Option<f64>) -> Result<Input<'static>> {
    let diag = shape.size()?;
    let tol = tolerance.unwrap_or(DEFAULT_RELATIVE_TOLERANCE * diag);
    if tol < MIN_RELATIVE_TOLERANCE * diag {
        return Err(EngineError::Params(format!(
            "`tolerance` must be at least {MIN_RELATIVE_TOLERANCE:e} of body `{}`'s bounding-box diagonal",
            shape.name
        )));
    }
    let surfaces = shape.face_surfaces()?.into_iter().map(|f| f.map(analytic)).collect();
    let face_names = shape.faces().iter().map(|f| f.name.clone()).collect();
    let (mesh, pieces) = shape.tessellate_pieces(tol)?;
    Ok(Input { name: shape.name.clone(), mesh: Cow::Owned(mesh), surfaces, face_names, diagonal: diag, seeds: Vec::new(), pieces })
}

/// The imported surface `surface` names, with the body name its mesh is tagged with.
pub(super) fn surface<'a>(s: &'a Session, p: &Value) -> Result<(String, &'a fieldcraft_geometry::tagged_mesh::TaggedSurface)> {
    let n = s.surfaces.len();
    let index = p.get("surface").and_then(Value::as_u64).and_then(|i| usize::try_from(i).ok()).filter(|&i| i < n);
    let index = index.ok_or_else(|| EngineError::Params(format!("`surface` must be the index of an imported surface (there are {n})")))?;
    let t = s.surfaces.get(index).ok_or_else(|| EngineError::Params("no such surface".into()))?;
    let name = if t.name.trim().is_empty() { format!("Surface{index}") } else { t.name.clone() };
    Ok((name, t))
}

pub(super) fn inputs<'a>(s: &'a Session, p: &Value) -> Result<Vec<Input<'a>>> {
    let tolerance = number(p, "tolerance")?;
    match (p.get("body"), p.get("surface")) {
        (Some(b), None) => {
            let name = b.as_str().ok_or_else(|| EngineError::param("body", ProblemKind::NotAName))?;
            let shape = s
                .design
                .as_deref()
                .and_then(|d| d.body(name))
                .ok_or_else(|| EngineError::Params(format!("no body `{name}` (open a design first)")))?;
            Ok(vec![Input { seeds: super::seeds::edge_seeds(s)?, ..shape_input(shape, tolerance)? }])
        }
        (None, Some(_)) => {
            let (name, t) = surface(s, p)?;
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for q in &t.mesh.positions {
                for k in 0..3 {
                    lo[k] = lo[k].min(q[k]);
                    hi[k] = hi[k].max(q[k]);
                }
            }
            // Faces named by index, as mesh.tet.generate names them.
            let face_names = (0..t.faces).map(|f| f.to_string()).collect();
            Ok(vec![Input {
                name,
                mesh: Cow::Borrowed(&t.mesh),
                surfaces: Vec::new(),
                face_names,
                diagonal: diagonal(lo, hi),
                seeds: Vec::new(),
                pieces: Vec::new(),
            }])
        }
        (None, None) => {
            let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first, or give `surface`".into()))?;
            if design.bodies.is_empty() || design.bodies.len() > MAX_BODIES {
                return Err(EngineError::Params(format!("the design has {} bodies; meshing takes 1 to {MAX_BODIES}", design.bodies.len())));
            }
            let seeds = super::seeds::edge_seeds(s)?;
            design.bodies.iter().map(|b| Ok(Input { seeds: seeds.clone(), ..shape_input(b, tolerance)? })).collect()
        }
        _ => Err(EngineError::Params("give `body` or `surface`, not both".into())),
    }
}

pub const HEX_SWEEP: CommandSpec = CommandSpec {
    id: "mesh.hex.sweep",
    label: "Sweep Hex Mesh",
    params: r#"the parameters of mesh.hex.generate, with "method": "swept" (default: extrude, sweep or revolve, whichever fits first) | "extrude" | "sweep" | "revolve" — mesh.hex.generate limited to the sweeping methods"#,
    form: Some(super::forms::HEX_SWEEP),
    kind: CommandKind::Model,
    requires: &[Requirement::Geometry],
    run: hex_sweep,
};

pub const HEX_MAPPED: CommandSpec = CommandSpec {
    id: "mesh.hex.mapped",
    label: "Mapped Hex Mesh",
    params: r#"the parameters of mesh.hex.generate without "method" — mesh.hex.generate with the mapped method: a six-sided block (six faces of four edges, three edges at each of eight corners)"#,
    form: Some(super::forms::HEX_MAPPED),
    kind: CommandKind::Model,
    requires: &[Requirement::Geometry],
    run: hex_mapped,
};

pub const HEX_BLOCKS: CommandSpec = CommandSpec {
    id: "mesh.hex.blocks",
    label: "Block Hex Mesh",
    params: r#"{"bodies": ["Cell1", "Cell2"] (default: every body of the open design), "size": 2.0 (element size, a length; default 1/25 of the bounding-box diagonal of all of them), "order": 1 (hex8) or 2 (hex20, mid-edge nodes on the geometry), "tolerance": 0.05 (chordal tolerance of the tessellation, a length; default 1e-3 of that diagonal)} — replaces the cells' meshes (as mesh.tet.generate; other bodies keep theirs) with one conforming all-hex mesh of a part partitioned into cells (bodies that meet face to face, such as the pieces of solid.split_body in SolidCraft): every cell a mapped block (six four-sided faces; faces that meet smoothly, such as the two halves of a hole, count as one), division counts assigned over all cells together, each shared face meshed once and its nodes merged, radial lines graded off holes and fillets, size controls (mesh.size.set; one on a whole cell refines its surface) setting the counts and spacing through the graded size field. Elements are tagged with their cell's body and boundary faces with its B-rep faces (`weld`, `weld_distance` and `feature_angle` as in mesh.tet.generate). Reports the cells, the shared faces, counts and quality; a cell that is not a block, or cells that touch without sharing whole faces, is an error that says which"#,
    form: Some(super::forms::HEX_BLOCKS),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: hex_blocks,
};

pub(super) fn hex_blocks(s: &mut Session, p: &Value) -> Result<Value> {
    let size = number(p, "size")?;
    let quadratic = order(p)?;
    let design = s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let shapes: Vec<&Shape> = match p.get("bodies") {
        None => design.bodies.iter().collect(),
        Some(v) => {
            let names = v.as_array().ok_or_else(|| EngineError::param("bodies", ProblemKind::NotNames))?;
            // Bounded before the duplicate check, which is quadratic in the list.
            if names.len() > MAX_CELLS {
                return Err(EngineError::Params(format!("`bodies` lists {} bodies; at most {MAX_CELLS} are meshed together", names.len())));
            }
            let mut out: Vec<&Shape> = Vec::with_capacity(names.len());
            for n in names {
                let name = n.as_str().ok_or_else(|| EngineError::param("bodies", ProblemKind::NotNames))?;
                let shape = design.body(name).ok_or_else(|| EngineError::Params(format!("no body `{name}`")))?;
                if out.iter().any(|b| b.name == shape.name) {
                    return Err(EngineError::Params(format!("body `{name}` is listed twice")));
                }
                out.push(shape);
            }
            out
        }
    };
    if shapes.is_empty() || shapes.len() > MAX_CELLS {
        return Err(EngineError::Params(format!("{} bodies; mesh.hex.blocks meshes 1 to {MAX_CELLS}", shapes.len())));
    }
    // One size and one tessellation tolerance for all cells, from their joint extent, so shared
    // curves are divided alike.
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for b in &shapes {
        let m = b.measure()?;
        for k in 0..3 {
            lo[k] = lo[k].min(m.bbox_min[k]);
            hi[k] = hi[k].max(m.bbox_max[k]);
        }
    }
    let diag = diagonal(lo, hi);
    if !(diag > 0.0 && diag.is_finite()) {
        return Err(EngineError::Params("the bodies have no extent".into()));
    }
    let tolerance = number(p, "tolerance")?.unwrap_or(DEFAULT_RELATIVE_TOLERANCE * diag);
    if tolerance < MIN_RELATIVE_TOLERANCE * diag {
        return Err(EngineError::Params(format!("`tolerance` must be at least {MIN_RELATIVE_TOLERANCE:e} of the bodies' bounding-box diagonal")));
    }
    let h = size.unwrap_or(DEFAULT_RELATIVE_SIZE * diag);
    if h < MIN_RELATIVE_SIZE * diag {
        return Err(EngineError::Params(format!("`size` must be at least {MIN_RELATIVE_SIZE:e} of the bodies' bounding-box diagonal")));
    }
    let inputs: Vec<Input> = shapes.iter().map(|b| shape_input(b, Some(tolerance))).collect::<Result<_>>()?;
    let mut surfaces = Vec::with_capacity(inputs.len());
    for b in &inputs {
        let options = super::mesh::surface_options(p, &b.mesh)?;
        surfaces.push(
            BoundarySurface::new(&b.mesh.positions, &b.mesh.triangles, &b.mesh.tri_face, &options)
                .map_err(|e| EngineError::context(format!("body `{}`", b.name), e))?,
        );
    }
    let cells: Vec<Cell> = surfaces.iter().zip(&inputs).map(|(surface, b)| Cell { name: &b.name, surface, faces: &b.surfaces }).collect();
    let options = HexOptions {
        method: HexMethod::Mapped,
        quadratic,
        chordal_tolerance: tolerance,
        max_elements: fieldcraft_model::MAX_ELEMENTS / 8,
        seeds: super::seeds::edge_seeds(s)?,
        ..HexOptions::default()
    };
    let names: Vec<&str> = inputs.iter().map(|b| b.name.as_str()).collect();
    let field = super::sizes::sizing(s)?.joint(&names).field("", h)?;
    let cm = hex_mesh_cells(&cells, &field, &options).map_err(|e| match e {
        fieldcraft_mesh::Error::Model(fieldcraft_model::Error::Cancelled) => EngineError::Model(fieldcraft_model::Error::Cancelled),
        e => EngineError::Mesh(e),
    })?;
    let m = &cm.mesh;
    // Replace the cells' meshes (replace.rs); other bodies and the definitions stay.
    let model = &mut s.model;
    let cleared = super::replace::clear(model, inputs.iter().map(|b| b.name.as_str()));
    add_cells(model, &inputs, &cm)?;
    let blocks: Vec<Value> =
        inputs.iter().zip(&cm.blocks).map(|(b, d)| json!({"body": b.name, "divisions": d, "elements": d.iter().product::<usize>()})).collect();
    let reply = json!({
        "cells": inputs.len(),
        "interfaces": cm.interfaces,
        "description": crate::units::mesh_text(s, &m.report.description),
        "kind": m.kind.name(),
        "elements": m.element_count(),
        "nodes": m.nodes.len(),
        "boundary_faces": m.boundary.len(),
        "quality": {"min_scaled_jacobian": m.report.min_scaled_jacobian, "mean_scaled_jacobian": m.report.mean_scaled_jacobian},
        "volume": m.report.volume,
        "surface_volume": m.report.surface_volume,
        "size": h,
        "blocks": blocks,
    });
    Ok(cleared.report(s, reply))
}

/// The `order` parameter: false for hex8, true for hex20.
fn order(p: &Value) -> Result<bool> {
    match p.get("order") {
        None => Ok(false),
        Some(v) => match v.as_u64() {
            Some(1) => Ok(false),
            Some(2) => Ok(true),
            _ => Err(EngineError::invalid("order", "must be 1 (hex8) or 2 (hex20)")),
        },
    }
}

fn hex_generate(s: &mut Session, p: &Value) -> Result<Value> {
    let all = [HexMethod::Auto, HexMethod::Swept, HexMethod::Extrude, HexMethod::Sweep, HexMethod::Revolve, HexMethod::Mapped];
    generate(s, p, HexMethod::Auto, &all)
}

fn hex_sweep(s: &mut Session, p: &Value) -> Result<Value> {
    generate(s, p, HexMethod::Swept, &[HexMethod::Swept, HexMethod::Extrude, HexMethod::Sweep, HexMethod::Revolve])
}

fn hex_mapped(s: &mut Session, p: &Value) -> Result<Value> {
    generate(s, p, HexMethod::Mapped, &[HexMethod::Mapped])
}

/// Mesh with `default` or the `method` parameter, which must be one of `allowed`.
fn generate(s: &mut Session, p: &Value, default: HexMethod, allowed: &[HexMethod]) -> Result<Value> {
    let size = number(p, "size")?;
    let quadratic = order(p)?;
    let method = match p.get("method") {
        None => default,
        Some(v) => v.as_str().and_then(HexMethod::from_name).filter(|m| allowed.contains(m)).ok_or_else(|| {
            let names: Vec<&str> = allowed.iter().map(|m| m.name()).collect();
            EngineError::Params(format!("`method` must be one of {}", names.join(", ")))
        })?,
    };
    let bodies = inputs(s, p)?;
    let sizing = super::sizes::sizing(s)?;
    let mut meshes: Vec<(Input, HexMesh)> = Vec::with_capacity(bodies.len());
    for b in bodies {
        let mesh = hex_body(p, &b, method, quadratic, size, &sizing)?;
        meshes.push((b, mesh));
    }
    // Replace the bodies' meshes (replace.rs); other bodies and the definitions stay.
    // A copy: the inputs borrow the session.
    let mut model = s.model.clone();
    let cleared = super::replace::clear(&mut model, meshes.iter().map(|(b, _)| b.name.as_str()));
    let mut report = Vec::with_capacity(meshes.len());
    for (b, m) in &meshes {
        add_hex(&mut model, b, m)?;
        report.push(json!({
            "body": b.name,
            "method": m.report.method.name(),
            "description": crate::units::mesh_text(s, &m.report.description),
            "kind": m.kind.name(),
            "elements": m.element_count(),
            "nodes": m.nodes.len(),
            "boundary_faces": m.boundary.len(),
            "quality": {"min_scaled_jacobian": m.report.min_scaled_jacobian, "mean_scaled_jacobian": m.report.mean_scaled_jacobian},
            "volume": m.report.volume,
            "surface_volume": m.report.surface_volume,
        }));
    }
    let (nodes, elements) = (model.nodes().len(), model.elements().len());
    s.model = model;
    Ok(cleared.report(s, json!({"bodies": report, "nodes": nodes, "elements": elements})))
}

/// Mesh one body with `method` (hex20 when `quadratic`) at `size` (default: from its extent) and
/// the size controls of `sizing`; the model is not touched. Errors name the body.
pub(super) fn hex_body(p: &Value, b: &Input, method: HexMethod, quadratic: bool, size: Option<f64>, sizing: &Sizing) -> Result<HexMesh> {
    if !(b.diagonal > 0.0 && b.diagonal.is_finite()) {
        return Err(EngineError::Params(format!("body `{}` has no extent", b.name)));
    }
    let h = sizing.body_size(&b.name, size.unwrap_or(DEFAULT_RELATIVE_SIZE * b.diagonal));
    if h < MIN_RELATIVE_SIZE * b.diagonal {
        return Err(EngineError::Params(format!("`size` must be at least {MIN_RELATIVE_SIZE:e} of body `{}`'s bounding-box diagonal", b.name)));
    }
    let surface_options = super::mesh::surface_options(p, &b.mesh)?;
    let surface = BoundarySurface::new(&b.mesh.positions, &b.mesh.triangles, &b.mesh.tri_face, &surface_options)
        .map_err(|e| EngineError::context(format!("body `{}`", b.name), e))?;
    let chord = if b.surfaces.is_empty() { 0.0 } else { number(p, "tolerance")?.unwrap_or(DEFAULT_RELATIVE_TOLERANCE * b.diagonal) };
    let options = HexOptions {
        method,
        quadratic,
        chordal_tolerance: chord,
        max_elements: fieldcraft_model::MAX_ELEMENTS / 8,
        seeds: b.seeds.clone(),
        ..HexOptions::default()
    };
    let field = sizing.field(&b.name, h)?;
    hex_mesh(&surface, &b.surfaces, &field, &options).map_err(|e| EngineError::of_body(&b.name, e))
}

/// Add a mesh of several cells (`inputs[c]` is cell `c`) after the model's largest labels, its
/// nodes shared between the cells, tagging each element with its cell's body and the boundary
/// element faces with the cell's B-rep faces.
pub(super) fn add_cells(model: &mut Model, inputs: &[Input], cm: &fieldcraft_mesh::hex::CellsMesh) -> Result<()> {
    let m = &cm.mesh;
    let (first_node, first_element) = super::fresh_labels(model, m.nodes.len(), m.element_count())?;
    for (i, xyz) in m.nodes.iter().enumerate() {
        model.set_node(first_node + i as u32, *xyz)?;
    }
    let mut body_tags = Vec::with_capacity(inputs.len());
    for b in inputs {
        body_tags.push(model.geometry_tag(&GeometryRef::body(b.name.as_str())?)?);
    }
    for e in 0..m.element_count() {
        let cell = cm.element_cell.get(e).map_or(0, |&c| c as usize);
        let nodes: Vec<u32> = m.element(e).unwrap_or(&[]).iter().map(|&i| first_node + i).collect();
        let l = first_element + e as u32;
        model.set_element(l, m.kind, nodes)?;
        model.tag_element(l, body_tags.get(cell).copied().ok_or_else(|| EngineError::Params("element without a cell".into()))?)?;
    }
    tag_faces(model, inputs, |e| cm.element_cell.get(e as usize).map_or(0, |&c| c as usize), &m.boundary, first_element)
}

/// Add one body's hex mesh after the model's largest labels, tagging the elements with the body
/// and the boundary element faces with its B-rep faces; returns the first node and element
/// labels.
pub(super) fn add_hex(model: &mut Model, b: &Input, m: &HexMesh) -> Result<(u32, u32)> {
    let (first_node, first_element) = super::fresh_labels(model, m.nodes.len(), m.element_count())?;
    for (i, xyz) in m.nodes.iter().enumerate() {
        model.set_node(first_node + i as u32, *xyz)?;
    }
    let body_tag = model.geometry_tag(&GeometryRef::body(b.name.as_str())?)?;
    for e in 0..m.element_count() {
        let nodes: Vec<u32> = m.element(e).unwrap_or(&[]).iter().map(|&i| first_node + i).collect();
        let l = first_element + e as u32;
        model.set_element(l, m.kind, nodes)?;
        model.tag_element(l, body_tag)?;
    }
    tag_faces(model, std::slice::from_ref(b), |_| 0, &m.boundary, first_element)?;
    Ok((first_node, first_element))
}

/// Tag every boundary element face with the B-rep face it lies on, of the body
/// (`bodies[cell(element)]`) its element meshes.
pub(super) fn tag_faces(
    model: &mut Model,
    bodies: &[Input<'_>],
    cell: impl Fn(u32) -> usize,
    boundary: &[BoundaryFace],
    first_element: u32,
) -> Result<()> {
    let mut tags: std::collections::BTreeMap<(usize, u32), u32> = std::collections::BTreeMap::new();
    for f in boundary {
        let c = cell(f.element);
        let tag = match tags.get(&(c, f.brep_face)) {
            Some(&t) => t,
            None => {
                let b = bodies.get(c).ok_or_else(|| EngineError::Params("a boundary face without its body".into()))?;
                let i = f.brep_face as usize;
                let name = b.face_names.get(i).cloned().unwrap_or_else(|| i.to_string());
                let t = model.geometry_tag(&GeometryRef::face(b.name.as_str(), name)?)?;
                tags.insert((c, f.brep_face), t);
                t
            }
        };
        model.tag_element_face(first_element + f.element, f.face, tag)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_geometry::sample::Plate;

    fn opened() -> (crate::tests::TempDir, Session) {
        let (dir, path) = crate::tests::sample_design_file(&Plate::default());
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    #[test]
    fn plate_with_a_hole_meshes_with_hexahedra_and_face_sets_resolve() {
        let (_dir, mut s) = opened();
        let r = s.execute("mesh.hex.generate", &json!({"body": "Body1", "size": 2.5})).unwrap();
        let body = &r["bodies"][0];
        assert_eq!(body["method"], "extrude", "{r}");
        assert_eq!(body["kind"], "hex8");
        let n = body["elements"].as_u64().unwrap() as usize;
        assert!(n > 100, "{r}");
        assert!(body["quality"]["min_scaled_jacobian"].as_f64().unwrap() > 0.3, "{r}");
        assert_eq!(s.model.element_counts(), vec![(fieldcraft_model::ElementKind::Hex8, n)]);
        // Volume: the plate minus the hole, to the tessellation's accuracy.
        let v = 40.0 * 30.0 * 10.0 - std::f64::consts::PI * 25.0 * 10.0;
        assert!((body["volume"].as_f64().unwrap() - v).abs() < 0.01 * v, "{r}");
        // A set on the top face resolves to element faces.
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let faces: Vec<String> = bodies[0]["faces"].as_array().unwrap().iter().map(|f| f["name"].as_str().unwrap().to_string()).collect();
        let top = faces.iter().find(|f| f.ends_with("top") || f.contains("end")).unwrap_or(&faces[0]).clone();
        s.execute("model.set.create", &json!({"name": "top", "geometry": [{"body": "Body1", "face": top}]})).unwrap();
        let set = s.execute("model.set.resolve", &json!({"name": "top"})).unwrap();
        assert!(set["faces"]["count"].as_u64().unwrap() > 0, "{set}");
        // Every boundary face is tagged: the skin of the body set.
        s.execute("model.set.create", &json!({"name": "skin", "geometry": [{"body": "Body1"}]})).unwrap();
        let skin = s.execute("model.set.resolve", &json!({"name": "skin"})).unwrap();
        assert_eq!(skin["faces"]["count"], body["boundary_faces"], "{skin}");
        assert_eq!(skin["elements"]["count"].as_u64().unwrap() as usize, n);
        // Hex20 with mid-edge nodes on the hole.
        let r2 = s.execute("mesh.hex.generate", &json!({"order": 2, "size": 2.5})).unwrap();
        assert_eq!(r2["bodies"][0]["kind"], "hex20");
        assert_eq!(s.model.element_counts(), vec![(fieldcraft_model::ElementKind::Hex20, n)]);
        let v2 = r2["bodies"][0]["volume"].as_f64().unwrap();
        assert!((v2 - v).abs() < (body["volume"].as_f64().unwrap() - v).abs() + 1e-6 * v, "hex20 {v2} vs hex8 {} vs {v}", body["volume"]);
        // The set still resolves after the remesh.
        let set = s.execute("model.set.resolve", &json!({"name": "top"})).unwrap();
        assert!(set["faces"]["count"].as_u64().unwrap() > 0);
        // Undo returns the hex8 mesh.
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.model.element_counts(), vec![(fieldcraft_model::ElementKind::Hex8, n)]);
    }

    /// A 20 × 10 × 10 box split by SolidCraft into four cells (x = 10, then y = 5 on each half),
    /// opened as a design; with `t_junction` the second half is not split, so the first half's
    /// two cells meet it on one face.
    fn split_box(t_junction: bool) -> (crate::tests::TempDir, Session) {
        let mut commands = vec![
            json!({"command": "solid.box", "params": {"length": 20, "width": 10, "height": 10, "body_name": "Box"}}),
            json!({"command": "solid.split_body", "params": {"body": "Box", "plane": {"origin": [10, 0, 0], "normal": [1, 0, 0]}}}),
            json!({"command": "solid.split_body", "params": {"body": "Box (1)", "plane": {"origin": [0, 5, 0], "normal": [0, 1, 0]}}}),
        ];
        if !t_junction {
            commands.push(json!({"command": "solid.split_body", "params": {"body": "Box", "plane": {"origin": [0, 5, 0], "normal": [0, 1, 0]}}}));
        }
        let bytes = fieldcraft_geometry::script_design(&json!({"commands": commands}).to_string()).unwrap();
        let (dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        (dir, s)
    }

    #[test]
    fn split_cells_mesh_into_one_conforming_hex_mesh() {
        let (_dir, mut s) = split_box(false);
        let r = s.execute("mesh.hex.blocks", &json!({"size": 2.5})).unwrap();
        assert_eq!(r["cells"], 4, "{r}");
        assert_eq!(r["interfaces"], 4, "{r}");
        assert_eq!(r["elements"], 8 * 4 * 4, "{r}");
        // One 9 × 5 × 5 grid of nodes: the cells' nodes are merged on the interfaces.
        assert_eq!(r["nodes"], 9 * 5 * 5, "{r}");
        assert_eq!(s.model.element_counts(), vec![(fieldcraft_model::ElementKind::Hex8, 128)]);
        assert!((r["volume"].as_f64().unwrap() - 2000.0).abs() < 1e-6, "{r}");
        // Each cell's elements carry its body, and the outer faces their B-rep faces: a set of
        // all four bodies skins the box with no internal faces.
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let names: Vec<Value> = bodies.as_array().unwrap().iter().map(|b| json!({"body": b["name"]})).collect();
        s.execute("model.set.create", &json!({"name": "all", "geometry": names})).unwrap();
        let all = s.execute("model.set.resolve", &json!({"name": "all"})).unwrap();
        assert_eq!(all["elements"]["count"], 128, "{all}");
        assert_eq!(all["faces"]["count"], r["boundary_faces"], "{all}");
        assert_eq!(r["boundary_faces"], 2 * (8 * 4 + 8 * 4 + 4 * 4), "{r}");
        let first = bodies[0]["name"].as_str().unwrap();
        s.execute("model.set.create", &json!({"name": "one", "geometry": [{"body": first}]})).unwrap();
        assert_eq!(s.execute("model.set.resolve", &json!({"name": "one"})).unwrap()["elements"]["count"], 32);
        // Hex20 shares mid-edge nodes across the interfaces as well.
        let q = s.execute("mesh.hex.blocks", &json!({"size": 2.5, "order": 2})).unwrap();
        assert_eq!(q["kind"], "hex20");
        assert_eq!(q["nodes"], 9 * 5 * 5 + 8 * 5 * 5 + 9 * 4 * 5 + 9 * 5 * 4, "{q}");
        // A list of bodies meshes those only.
        let two = s.execute("mesh.hex.blocks", &json!({"bodies": [bodies[0]["name"], bodies[1]["name"]], "size": 2.5})).unwrap();
        assert_eq!(two["cells"], 2);
    }

    /// The Partitioned Mesh form's "All hex" option is mesh.hex.blocks.
    #[test]
    fn mesh_cells_all_blocks_is_mesh_hex_blocks() {
        let (_dir, mut s) = split_box(false);
        let blocks = s.execute("mesh.hex.blocks", &json!({"size": 2.5})).unwrap();
        let mut cells = s.execute("mesh.cells", &json!({"default": "blocks", "size": 2.5})).unwrap();
        let mut blocks = blocks;
        // The same mesh, in place of the first one.
        let replaced = cells.as_object_mut().and_then(|o| o.remove("replaced")).unwrap();
        assert_eq!(replaced.as_array().map(Vec::len), Some(4), "{replaced}");
        blocks.as_object_mut().and_then(|o| o.remove("replaced"));
        assert_eq!(cells, blocks);
        let e = s.execute("mesh.cells", &json!({"default": "blocks", "methods": {"Box": "tet"}})).unwrap_err().to_string();
        assert!(e.contains("methods"), "{e}");
    }

    #[test]
    fn cells_that_do_not_conform_and_hostile_params_are_errors() {
        let (_dir, mut s) = split_box(true);
        let before = s.model.clone();
        let e = s.execute("mesh.hex.blocks", &json!({"size": 2.5})).unwrap_err().to_string();
        assert!(e.contains("do not conform") && e.contains("`Box.1`"), "{e}");
        assert_eq!(s.model, before);
        for p in [
            json!({"size": 0}),
            json!({"size": -1}),
            json!({"size": "big"}),
            json!({"size": 1e-300}),
            json!({"order": 3}),
            json!({"tolerance": -1}),
            json!({"tolerance": 1e-300}),
            json!({"bodies": "Box"}),
            json!({"bodies": []}),
            json!({"bodies": ["Nope"]}),
            json!({"bodies": ["Box", "Box"]}),
            json!({"bodies": [7]}),
            json!({"size": 1e-3}),
        ] {
            assert!(s.execute("mesh.hex.blocks", &p).is_err(), "{p}");
            assert_eq!(s.model, before, "{p}");
        }
        assert!(Session::default().execute("mesh.hex.blocks", &json!({})).is_err());
        // One cell that is no block says which.
        let (_dir, mut plate) = opened();
        let e = plate.execute("mesh.hex.blocks", &json!({})).unwrap_err().to_string();
        assert!(e.contains("`Body1`") && e.contains("not a logical cube"), "{e}");
    }

    #[test]
    fn hostile_params_are_errors_and_change_nothing() {
        let (_dir, mut s) = opened();
        let before = s.model.clone();
        for p in [
            json!({"body": "Body1", "size": 0}),
            json!({"body": "Body1", "size": -1}),
            json!({"body": "Body1", "size": 1e-300}),
            json!({"body": "Body1", "size": "big"}),
            json!({"body": "Body1", "size": 1e300}),
            json!({"body": "Body1", "order": 3}),
            json!({"body": "Body1", "order": "2"}),
            json!({"body": "Body1", "method": "octree"}),
            json!({"body": "Body1", "method": null}),
            json!({"body": "Body1", "tolerance": 1e-300}),
            json!({"body": "Body1", "tolerance": -1}),
            json!({"body": "Nope"}),
            json!({"body": 7}),
            json!({"surface": 0}),
            json!({"surface": -1}),
            json!({"body": "Body1", "surface": 0}),
            json!({"body": "Body1", "method": "revolve"}),
            json!({"body": "Body1", "method": "mapped"}),
        ] {
            let r = s.execute("mesh.hex.generate", &p);
            if p.get("size") == Some(&json!(1e300)) {
                // One element through the whole body is still a valid request.
                if r.is_ok() {
                    s.model = before.clone();
                }
                continue;
            }
            assert!(r.is_err(), "{p}: {r:?}");
            assert_eq!(s.model, before, "{p}");
        }
        assert!(Session::default().execute("mesh.hex.generate", &json!({})).is_err());
        assert!(s.execute("mesh.hex.sweep", &json!({"body": "Body1", "method": "mapped"})).is_err());
        assert!(s.execute("mesh.hex.mapped", &json!({"body": "Body1", "method": "sweep"})).is_err());
        // The plate is not a six-sided block; it is a sweep.
        assert!(s.execute("mesh.hex.mapped", &json!({"body": "Body1"})).unwrap_err().to_string().contains("logical cube"));
        assert_eq!(s.model, before);
        let r = s.execute("mesh.hex.sweep", &json!({"body": "Body1", "size": 4})).unwrap();
        assert_eq!(r["bodies"][0]["method"], "extrude");
    }
}
