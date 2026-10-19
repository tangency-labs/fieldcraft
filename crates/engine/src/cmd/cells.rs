//! Meshing a partitioned part cell by cell (`mesh.cells`): each cell (a body of the open design,
//! typically made by `geometry.partition`) with its own method, all into one mesh, and the cells
//! joined across the partitions that divide them.
//!
//! **Sweeping together (#288).** The cells of one body of the design whose method may sweep
//! (`auto`, `hex`, `swept`, `extrude`, `revolve`) are first meshed together
//! ([`hex_mesh_swept_cells`]): one cross-section with every cell's cut faces imprinted, one
//! interval assignment and one set of layers, so cells meet node for node across every cut — a
//! cut that ends on another included — with no edge seeds, and their nodes are shared. When they
//! do not sweep together in one frame (the reply says why, per cell, as `not_swept_together`),
//! each cell is meshed alone as below.
//!
//! **Joining.** Two cells on either side of a partition `P` meet on faces named
//! `partition:P:a` and `partition:P:b`. Where the element faces of the two sides match node for
//! node (the same divisions and kind on both sides, e.g. a mesh and its mirror image) the
//! coincident nodes are merged and the mesh is conforming there. Elsewhere — a hex cell next to a
//! tet cell, or two cells meshed with different divisions — the two surfaces are tied
//! (`interaction.tie`; dual mortar constraints in the implicit analyses, so a constant stress
//! crosses the interface exactly). Each tie goes one way, by a fixed order of the cells (hex cells
//! main to tet cells), and each cell's tied cut faces make one tie, so no node is a secondary
//! node of two of them and no interface is tied both ways: either would over-constrain it (#389).
//! The tie is the transition between hex and tet cells; pyramid elements are not used, because
//! they need their own element formulation in every solver and deck, while the tie works for any
//! pair of meshes and keeps each cell's mesh as its method made it. The reply lists every
//! interface and how it was joined.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_mesh::cdt::CdtOptions;
use fieldcraft_mesh::hex::{Cell, CellsMesh, HexMesh, HexMethod, HexOptions, hex_mesh_swept_cells};
use fieldcraft_mesh::surface::BoundarySurface;
use fieldcraft_model::geometry::Entity;
use fieldcraft_model::{ElementFace, GeometryRef, Model, Set, Tie};
use serde_json::{Value, json};

use super::hex::{Input, add_cells, add_hex, hex_body};
use super::mesh::{add_tets, mesh_body, quality_options};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const CELLS: CommandSpec = CommandSpec {
    id: "mesh.cells",
    label: "Partitioned Mesh",
    params: r#"{"methods": {"Body1/a": "tet", "Body1/b": "swept"} (per cell; default `default`), "default": "auto", "bodies": ["Body1/a", "Body1/b"] (default: every body of the open design), "size": 2.0 (element size, a length; default 1/25 of the bounding-box diagonal of all of them), "order": 1 | 2 (hex8 and tet4, or hex20 and tet10), "tolerance": 0.05 (chordal tolerance, a length), "transition": "pyramid" | "tie" (how a tet cell meets hex cells; default pyramid)} — replaces the cells' meshes (as mesh.tet.generate; other bodies keep theirs) with one mesh of the cells, each with its method: "auto" (structured hex where a method fits, else tet), "hex" (any structured hex method), "swept", "extrude", "sweep", "revolve", "mapped" or "tet"; or "default": "blocks" (with no `methods`): every cell a mapped block, division counts assigned over all cells together, exactly mesh.hex.blocks with the same parameters. Cells on either side of a partition (faces `partition:<name>:a` and `:b`) are joined: where their element faces match node for node the nodes are merged (conforming); between a tet cell and the hex cells across its partitions, with "transition": "pyramid" (linear meshes), the tet cell is meshed after them with a pyramid (pyramid5) on every hex element face of the cut, its tetrahedra sharing the pyramids' sides, so the mesh is conforming there (`joined`: `pyramids`, with their count); where that cannot be built (a cut the hex faces do not cover exactly, quadratic elements, a pyramid that would reach through a thin cell) the cell says why in `no_pyramids` and is tied; otherwise, as between a hex and a tet cell with "transition": "tie", they are tied: each cell's tied cut faces make one tie, `tie_<cell>` (`/` as `_`; sets `<tie>_main`, the facing cells' cut faces, and `<tie>_secondary`, the cell's own), hex cells main to tet cells and otherwise by name, so a node is a secondary node of one such tie only; the ties an earlier run made on these cells are replaced. The cells of one design body that may sweep (auto, hex, swept, extrude, revolve) are swept together first: one cross-section and one set of layers for all of them, so they share their nodes across every cut, T-junctions included, with no edge seeds (`swept_together` names the body; `not_swept_together` says why not, and those cells are meshed one by one). Reports each cell's method (and why auto fell back to tet), its `quality` (volume, smallest and mean scaled Jacobian, largest edge ratio), and each interface: `merged` with the nodes merged or shared, or `tied` (with `why` when both cells are of one design body: the tet cell and why auto fell back, or why the body's cells were not swept together)"#,
    form: Some(super::forms::CELLS),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: cells,
};

/// Distance (relative to the cells' bounding-box diagonal) within which interface nodes coincide.
const COINCIDENT: f64 = 1e-7;

/// A cell's method.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Auto,
    Hex(HexMethod),
    Tet,
}

impl Method {
    fn parse(v: &Value) -> Result<Method> {
        let name = v.as_str().unwrap_or_default();
        Ok(match name {
            "auto" => Method::Auto,
            "hex" => Method::Hex(HexMethod::Auto),
            "tet" => Method::Tet,
            _ => match HexMethod::from_name(name).filter(|m| *m != HexMethod::Auto) {
                Some(m) => Method::Hex(m),
                None => {
                    return Err(EngineError::Params(format!(
                        "method `{}` is not one of auto, hex, swept, extrude, sweep, revolve, mapped, tet",
                        name.chars().take(40).collect::<String>()
                    )));
                }
            },
        })
    }
}

mod transition;

/// A meshed cell; tet cells are meshed after the hex cells are in the model, so they can meet
/// them through pyramids.
enum Meshed<'a> {
    Hex(Box<(Input<'a>, HexMesh)>),
    Tet,
}

/// How a tet cell joins the hex cells across its partitions.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Transition {
    /// Pyramids on the hex faces, conforming (falls back to a tie with the reason).
    Pyramid,
    /// A tie (dual mortar), each mesh as its method made it.
    Tie,
}

fn cells(s: &mut Session, p: &Value) -> Result<Value> {
    // What the reasons the meshers give are shown in.
    let display = s.units().clone();
    if p.get("default").and_then(Value::as_str) == Some("blocks") {
        if p.get("methods").is_some_and(|m| m.as_object().is_none_or(|o| !o.is_empty())) {
            return Err(EngineError::param(
                "methods",
                crate::schema::ProblemKind::Invalid("must be empty with \"default\": \"blocks\" (every cell a block)".into()),
            ));
        }
        return super::hex::hex_blocks(s, p);
    }
    let size = super::hex::number(p, "size")?;
    let quadratic = match p.get("order") {
        None => false,
        Some(v) => match v.as_u64() {
            Some(1) => false,
            Some(2) => true,
            _ => return Err(EngineError::invalid("order", "must be 1 or 2")),
        },
    };
    let transition = match p.get("transition").map(|v| v.as_str()) {
        None | Some(Some("pyramid")) => Transition::Pyramid,
        Some(Some("tie")) => Transition::Tie,
        Some(_) => return Err(EngineError::invalid("transition", "must be \"pyramid\" or \"tie\"")),
    };
    let default = p.get("default").map(Method::parse).transpose()?.unwrap_or(Method::Auto);
    let methods: BTreeMap<String, Method> = match p.get("methods") {
        None => BTreeMap::new(),
        Some(Value::Object(m)) => m.iter().map(|(k, v)| Ok((k.clone(), Method::parse(v)?))).collect::<Result<_>>()?,
        Some(_) => return Err(EngineError::invalid("methods", "must map cell names to methods")),
    };
    let design = s.design.clone().ok_or_else(|| EngineError::Params("open a design first".into()))?;
    let names: Vec<String> = match p.get("bodies") {
        None => design.bodies.iter().map(|b| b.name.clone()).collect(),
        Some(_) => super::req(p, "bodies")?,
    };
    if names.is_empty() || names.len() > fieldcraft_mesh::hex::MAX_CELLS {
        return Err(EngineError::Params(format!("{} cells; mesh.cells meshes 1 to {}", names.len(), fieldcraft_mesh::hex::MAX_CELLS)));
    }
    for n in methods.keys() {
        if !names.contains(n) {
            return Err(EngineError::Params(format!("`methods` names `{n}`, which is not a cell being meshed")));
        }
    }
    // One size and tolerance from the cells' joint extent, as mesh.hex.blocks does.
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut shapes = Vec::with_capacity(names.len());
    for (i, n) in names.iter().enumerate() {
        if names.get(..i).is_some_and(|b| b.contains(n)) {
            return Err(EngineError::Params(format!("cell `{n}` is listed twice")));
        }
        let shape = design.body(n).ok_or_else(|| EngineError::Params(format!("no body `{n}`")))?;
        let m = shape.measure()?;
        for k in 0..3 {
            lo[k] = lo[k].min(m.bbox_min[k]);
            hi[k] = hi[k].max(m.bbox_max[k]);
        }
        shapes.push(shape);
    }
    let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    if !(diagonal > 0.0 && diagonal.is_finite()) {
        return Err(EngineError::Params("the cells have no extent".into()));
    }
    let size = Some(size.unwrap_or(diagonal / 25.0));
    let tolerance = super::hex::number(p, "tolerance")?.unwrap_or(1e-3 * diagonal);
    // The tet mesher's parameters: the same size and tolerance, its own order names.
    let mut tet_params = json!({"size": size, "tolerance": tolerance, "order": if quadratic { "quadratic" } else { "linear" }});
    for key in ["weld", "weld_distance", "feature_angle"] {
        if let (Some(v), Some(o)) = (p.get(key), tet_params.as_object_mut()) {
            o.insert(key.into(), v.clone());
        }
    }
    let mut hex_params = tet_params.clone();
    if let Some(o) = hex_params.as_object_mut() {
        o.remove("order");
    }
    // The cells are meshed at one size, so the size controls apply to all of them alike.
    let sizing = super::sizes::sizing(s)?.joint(&names.iter().map(String::as_str).collect::<Vec<_>>());
    let (tet_options, quality) = quality_options(&tet_params, s.limits.max_tets as u64)?;
    let tet_options = sizing.tet_options("", &tet_params, &tet_options);
    let cdt = CdtOptions { max_steiner: s.limits.max_steiner, max_tets: tet_options.max_tets };

    let seeds = super::seeds::edge_seeds(s)?;
    let method_of = |n: &str| methods.get(n).copied().unwrap_or(default);
    // The cells of one body that may sweep are swept together first (#288): one cross-section
    // and one set of layers for all of them, so they meet node for node across their cuts. Cells
    // that do not sweep together are meshed one by one, as below.
    let mut groups: BTreeMap<(String, HexMethod), Vec<usize>> = BTreeMap::new();
    for (i, shape) in shapes.iter().enumerate() {
        let sweep = match method_of(&shape.name) {
            Method::Auto | Method::Hex(HexMethod::Auto | HexMethod::Swept) => Some(HexMethod::Swept),
            Method::Hex(m @ (HexMethod::Extrude | HexMethod::Revolve)) => Some(m),
            Method::Hex(_) | Method::Tet => None,
        };
        if let (Some(m), Some(source)) = (sweep, design.source_of(&shape.name)) {
            groups.entry((source.to_string(), m)).or_default().push(i);
        }
    }
    let mut joint: Vec<Option<usize>> = vec![None; shapes.len()];
    let mut not_together: BTreeMap<usize, String> = BTreeMap::new();
    let mut swept: Vec<(String, Vec<usize>, Vec<Input>, CellsMesh)> = Vec::new();
    for ((source, m), members) in groups {
        if members.len() < 2 {
            continue;
        }
        fieldcraft_model::progress::stage(&format!("Cells: sweeping {source}'s cells together"), None, "");
        fieldcraft_model::progress::check()?;
        let inputs: Vec<Input> = members
            .iter()
            .filter_map(|&i| shapes.get(i))
            .map(|shape| Ok(Input { seeds: seeds.clone(), ..super::hex::shape_input(shape, Some(tolerance))? }))
            .collect::<Result<_>>()?;
        match swept_together(p, &inputs, m, quadratic, size, tolerance, &seeds, &sizing, &display)? {
            Ok(cm) => {
                for &i in &members {
                    if let Some(slot) = joint.get_mut(i) {
                        *slot = Some(swept.len());
                    }
                }
                swept.push((source, members, inputs, cm));
            }
            Err(why) => {
                for &i in &members {
                    not_together.insert(i, why.clone());
                }
            }
        }
    }
    let mut meshed: Vec<(Meshed, Value)> = Vec::with_capacity(shapes.len());
    for (i, shape) in shapes.iter().enumerate() {
        if joint.get(i).copied().flatten().is_some() {
            continue;
        }
        fieldcraft_model::progress::stage(&format!("Cells: meshing {}", shape.name), Some(shapes.len() as f64), "cells");
        fieldcraft_model::progress::count(i as f64);
        fieldcraft_model::progress::check()?;
        let method = method_of(&shape.name);
        let input = || Ok::<_, EngineError>(Input { seeds: seeds.clone(), ..super::hex::shape_input(shape, Some(tolerance))? });
        let hex = |m: HexMethod| -> Result<(Input<'static>, HexMesh)> {
            let b = input()?;
            let mesh = hex_body(&hex_params, &b, m, quadratic, size, &sizing)?;
            Ok((b, mesh))
        };
        let (cell, mut info) = match method {
            Method::Hex(m) => {
                let (b, mesh) = hex(m)?;
                let info = json!({"method": mesh.report.method.name()});
                (Meshed::Hex(Box::new((b, mesh))), info)
            }
            Method::Tet => (Meshed::Tet, json!({"method": "tet"})),
            Method::Auto => match hex(HexMethod::Auto) {
                Ok((b, mesh)) => {
                    let info = json!({"method": mesh.report.method.name()});
                    (Meshed::Hex(Box::new((b, mesh))), info)
                }
                // A cancelled job stops; any other failure falls back to tet.
                Err(e) if e.is_cancelled() => return Err(e),
                Err(e) => (Meshed::Tet, json!({"method": "tet", "why": e.shown(crate::units::geometry(), &display).to_string()})),
            },
        };
        if let Some(o) = info.as_object_mut() {
            o.insert("body".into(), json!(shape.name));
            if let Some(why) = not_together.get(&i) {
                o.insert("not_swept_together".into(), json!(why));
            }
        }
        meshed.push((cell, info));
    }

    // Replace the cells' meshes (replace.rs), each after the labels before it; other bodies
    // and the definitions stay.
    let describe = |m: &fieldcraft_mesh::message::Message| m.render(&crate::units::MeshUnits { from: crate::units::geometry(), to: &display });
    let model = &mut s.model;
    let cleared = super::replace::clear(model, names.iter().map(String::as_str));
    let mut report: Vec<(usize, Value)> = Vec::with_capacity(shapes.len());
    let mut interfaces = Vec::new();
    for (source, members, inputs, cm) in &swept {
        add_cells(model, inputs, cm)?;
        let m = &cm.mesh;
        for (c, (&i, b)) in members.iter().zip(inputs).enumerate() {
            let elements: Vec<usize> = (0..m.element_count()).filter(|&e| cm.element_cell.get(e).is_some_and(|&x| x as usize == c)).collect();
            let nodes: BTreeSet<u32> = elements.iter().filter_map(|&e| m.element(e)).flatten().copied().collect();
            report.push((
                i,
                json!({"body": b.name, "method": m.report.method.name(), "swept_together": source, "description": describe(&m.report.description), "kind": m.kind.name(), "elements": elements.len(), "nodes": nodes.len()}),
            ));
        }
        interfaces.extend(swept_interfaces(&design, inputs, cm));
    }
    // The cells meshed alone, in the order of `shapes`: the hex ones first, then the tet ones,
    // which may meet the hex cells through pyramids.
    let alone: Vec<usize> = (0..shapes.len()).filter(|&i| joint.get(i).copied().flatten().is_none()).collect();
    let mut hex_cells: Vec<String> = swept.iter().flat_map(|(_, _, inputs, _)| inputs.iter().map(|b| b.name.clone())).collect();
    let mut tet_cells = Vec::new();
    for (&i, (cell, info)) in alone.iter().zip(meshed.iter_mut()) {
        let mut info = std::mem::take(info);
        match &*cell {
            Meshed::Hex(hex) => {
                let (b, m) = &**hex;
                add_hex(model, b, m)?;
                hex_cells.push(b.name.clone());
                if let Some(o) = info.as_object_mut() {
                    o.insert("kind".into(), json!(m.kind.name()));
                    o.insert("elements".into(), json!(m.element_count()));
                    o.insert("nodes".into(), json!(m.nodes.len()));
                }
                report.push((i, info));
            }
            Meshed::Tet => tet_cells.push((i, info)),
        }
    }
    for (i, mut info) in tet_cells {
        let Some(shape) = shapes.get(i) else { continue };
        fieldcraft_model::progress::stage(&format!("Cells: meshing {} (tet)", shape.name), Some(shapes.len() as f64), "cells");
        fieldcraft_model::progress::count(i as f64);
        fieldcraft_model::progress::check()?;
        let input = Input { seeds: seeds.clone(), ..super::hex::shape_input(shape, Some(tolerance))? };
        let pyramids = if transition == Transition::Tie {
            Err("`transition` is \"tie\"".to_string())
        } else if quadratic {
            Err("pyramid transitions are linear (pyramid5); order 2 is tied".to_string())
        } else {
            transition::mesh_with_pyramids(&design, model, &tet_params, &input, &hex_cells, &cdt, &tet_options, &display)?
        };
        let (kind, elements, nodes) = match pyramids {
            Ok(c) => {
                interfaces.extend(c.interfaces.iter().cloned());
                let kind = c.tets.mesh.kind;
                let (elements, nodes) = transition::add(model, c)?;
                (kind, elements, nodes)
            }
            Err(why) => {
                if let Some(o) = info.as_object_mut() {
                    o.insert("no_pyramids".into(), json!(why));
                }
                let t = mesh_body(Some(&design), &tet_params, input, &cdt, &tet_options, quality, true)?;
                add_tets(model, &t.name, &t.face_names, &t.mesh)?;
                (t.mesh.kind, t.mesh.elements.len(), t.mesh.points.len())
            }
        };
        if let Some(o) = info.as_object_mut() {
            o.insert("kind".into(), json!(kind.name()));
            o.insert("elements".into(), json!(elements));
            o.insert("nodes".into(), json!(nodes));
        }
        report.push((i, info));
    }
    report.sort_by_key(|r| r.0);
    let report: Vec<Value> = report
        .into_iter()
        .map(|(_, mut info)| {
            let quality = info.get("body").and_then(Value::as_str).map(|b| cell_quality(model, b));
            if let (Some(o), Some(q)) = (info.as_object_mut(), quality) {
                o.insert("quality".into(), q);
            }
            info
        })
        .collect();
    let tet: BTreeSet<String> = report
        .iter()
        .filter(|r| r.get("kind").and_then(Value::as_str).is_some_and(|k| k.starts_with("tet")))
        .filter_map(|r| r.get("body").and_then(Value::as_str).map(String::from))
        .collect();
    fieldcraft_model::progress::stage("Cells: joining the cells", None, "");
    fieldcraft_model::progress::check()?;
    interfaces.extend(join(model, &names, &tet, COINCIDENT * diagonal)?);
    explain_ties(&design, &report, &mut interfaces);
    let (nodes, elements) = (model.nodes().len(), model.elements().len());
    Ok(cleared.report(s, json!({"cells": report, "interfaces": interfaces, "nodes": nodes, "elements": elements, "size": size})))
}

/// Sweep the cells `inputs` together ([`hex_mesh_swept_cells`]) with `method` (extrude,
/// revolve, or swept: either). `Ok(Err(why))` when they do not sweep together; limits and a
/// cancel are errors.
fn swept_together(
    p: &Value,
    inputs: &[Input],
    method: HexMethod,
    quadratic: bool,
    size: Option<f64>,
    tolerance: f64,
    seeds: &[fieldcraft_mesh::seeds::EdgeSeed],
    sizing: &super::sizes::Sizing,
    display: &fieldcraft_model::DisplaySystem,
) -> Result<std::result::Result<CellsMesh, String>> {
    let h = size.ok_or_else(|| EngineError::Params("no element size".into()))?;
    let field = sizing.field("", h)?;
    let mut surfaces = Vec::with_capacity(inputs.len());
    for b in inputs {
        let options = super::mesh::surface_options(p, &b.mesh)?;
        match BoundarySurface::new(&b.mesh.positions, &b.mesh.triangles, &b.mesh.tri_face, &options) {
            Ok(surface) => surfaces.push(surface),
            Err(e) => return Ok(Err(format!("body `{}`: {}", b.name, crate::units::mesh_error_text(&e, display)))),
        }
    }
    let cells: Vec<Cell> = surfaces.iter().zip(inputs).map(|(surface, b)| Cell { name: &b.name, surface, faces: &b.surfaces }).collect();
    let options = HexOptions {
        method,
        quadratic,
        chordal_tolerance: tolerance,
        max_elements: fieldcraft_model::MAX_ELEMENTS / 8,
        seeds: seeds.to_vec(),
        ..HexOptions::default()
    };
    match hex_mesh_swept_cells(&cells, &field, &options) {
        Ok(cm) => Ok(Ok(cm)),
        Err(fieldcraft_mesh::Error::Model(fieldcraft_model::Error::Cancelled)) => Err(EngineError::Model(fieldcraft_model::Error::Cancelled)),
        Err(e) if e.kind() == Some(fieldcraft_mesh::Kind::Limit) => Err(e.into()),
        Err(e) => Ok(Err(crate::units::mesh_error_text(&e, display))),
    }
}

/// The interfaces of cells swept together: each pair of cells whose elements share faces, with
/// the partition between them and the nodes they share (one mesh: merged by construction).
fn swept_interfaces(design: &fieldcraft_geometry::Design, inputs: &[Input], cm: &CellsMesh) -> Vec<Value> {
    let m = &cm.mesh;
    let mut owner: BTreeMap<[u32; 4], u32> = BTreeMap::new();
    let mut pairs: BTreeMap<(u32, u32), BTreeSet<u32>> = BTreeMap::new();
    for e in 0..m.element_count() {
        let (Some(h), Some(&c)) = (m.element(e), cm.element_cell.get(e)) else { continue };
        for local in fieldcraft_model::ElementKind::Hex8.faces() {
            let mut key = [0u32; 4];
            for (slot, &i) in key.iter_mut().zip(local.iter()) {
                *slot = h.get(i).copied().unwrap_or(0);
            }
            key.sort_unstable();
            match owner.insert(key, c) {
                Some(o) if o != c => pairs.entry((o.min(c), o.max(c))).or_default().extend(key),
                _ => {}
            }
        }
    }
    // The partitions whose cut faces a cell has, with the side: two cells meet across a
    // partition they are on opposite sides of.
    let cuts = |c: u32| -> BTreeSet<(String, String)> {
        let name = inputs.get(c as usize).map_or("", |b| b.name.as_str());
        let faces = design.body(name).map(|b| b.faces().iter().map(|f| f.name.clone()).collect::<Vec<_>>()).unwrap_or_default();
        faces
            .iter()
            .filter_map(|f| {
                let (p, side) = f.strip_prefix("partition:")?.rsplit_once(':')?;
                let side = side.split_once('#').map_or(side, |(s, _)| s);
                Some((p.to_string(), side.to_string()))
            })
            .collect()
    };
    pairs
        .into_iter()
        .map(|((a, b), nodes)| {
            let (na, nb) = (inputs.get(a as usize).map_or("", |x| x.name.as_str()), inputs.get(b as usize).map_or("", |x| x.name.as_str()));
            let theirs = cuts(b);
            let partition =
                cuts(a).into_iter().find(|(p, side)| theirs.iter().any(|(q, other)| q == p && other != side)).map(|(p, _)| p).unwrap_or_default();
            json!({"partition": partition, "cells": [na, nb], "joined": "merged", "nodes": nodes.len(), "swept_together": true})
        })
        .collect()
}

/// A cell's mesh in numbers: its volume, smallest and mean scaled Jacobian and largest edge
/// ratio.
fn cell_quality(model: &Model, body: &str) -> Value {
    let (mut volume, mut min, mut sum, mut n, mut ratio) = (0.0, f64::INFINITY, 0.0, 0usize, None::<f64>);
    for label in model.body_elements(body) {
        let Some(e) = model.element(label) else { continue };
        let xyz: Vec<[f64; 3]> = e.nodes.iter().filter_map(|&l| model.node(l).map(|x| x.xyz)).collect();
        let Ok(q) = fieldcraft_mesh::quality::element_quality(e.kind, &xyz) else { continue };
        volume += q.size;
        if let Some(r) = q.get(fieldcraft_mesh::quality::Metric::EdgeRatio) {
            ratio = Some(ratio.map_or(r, |x| x.max(r)));
        }
        if let Some(sj) = q.get(fieldcraft_mesh::quality::Metric::ScaledJacobian) {
            min = min.min(sj);
            sum += sj;
            n += 1;
        }
    }
    let min = if n > 0 { Some(min) } else { None };
    json!({"volume": volume, "min_scaled_jacobian": min, "mean_scaled_jacobian": if n > 0 { Some(sum / n as f64) } else { None }, "max_edge_ratio": ratio})
}

/// One side of an interface: a cell's element faces on one side of a partition.
struct Side {
    body: String,
    /// The face name (`partition:P:a`), which selects its pieces too.
    face: String,
    faces: Vec<ElementFace>,
    nodes: BTreeSet<u32>,
    lo: [f64; 3],
    hi: [f64; 3],
}

/// The element faces of each cell tagged with a partition's cut faces, by partition and side.
fn sides(model: &Model, cells: &[String]) -> BTreeMap<String, [Vec<Side>; 2]> {
    let c = model.classification();
    // Tag → (partition, side), for cut faces of the meshed cells.
    let mut cut: BTreeMap<u32, (String, usize, String)> = BTreeMap::new();
    for (t, r) in c.entities().iter().enumerate() {
        let Entity::Face(name) = &r.entity else { continue };
        if !cells.contains(&r.body) {
            continue;
        }
        let base = name.split_once('#').map_or(name.as_str(), |(b, _)| b);
        let Some(rest) = base.strip_prefix("partition:") else { continue };
        let Some((partition, side)) = rest.rsplit_once(':') else { continue };
        let side = match side {
            "a" => 0,
            "b" => 1,
            _ => continue,
        };
        if let Ok(t) = u32::try_from(t) {
            cut.insert(t, (partition.to_string(), side, r.body.clone()));
        }
    }
    let mut grouped: BTreeMap<(String, usize, String), Vec<ElementFace>> = BTreeMap::new();
    for (&element, tags) in c.elements() {
        for (&face, t) in &tags.faces {
            if let Some(key) = cut.get(t) {
                grouped.entry(key.clone()).or_default().push(ElementFace { element, face });
            }
        }
    }
    let mut out: BTreeMap<String, [Vec<Side>; 2]> = BTreeMap::new();
    for ((partition, side, body), faces) in grouped {
        let nodes: BTreeSet<u32> = faces.iter().filter_map(|f| model.face_nodes(*f)).flatten().collect();
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for n in nodes.iter().filter_map(|&n| model.node(n)) {
            for k in 0..3 {
                lo[k] = lo[k].min(n.xyz[k]);
                hi[k] = hi[k].max(n.xyz[k]);
            }
        }
        let face = format!("partition:{partition}:{}", if side == 0 { "a" } else { "b" });
        if let Some(slot) = out.entry(partition).or_insert_with(|| [Vec::new(), Vec::new()]).get_mut(side) {
            slot.push(Side { body, face, faces, nodes, lo, hi });
        }
    }
    out
}

/// Join the cells across every partition: merge where the two sides match node for node, tie
/// where they do not. Returns the interfaces.
///
/// Ties go one way between any two cells, by a fixed order of the cells (hex before tet, then by
/// name; the earlier one is the main side), so no two ties hold one interface both ways, and
/// each cell's tied cut faces make one tie, `tie_<cell>` (its secondary surface; the main one is
/// the facing cells' cut faces), so a node where two of them meet is a secondary node once.
/// The ties of an earlier run on these cells are replaced.
fn join(model: &mut Model, cells: &[String], tet: &BTreeSet<String>, tol: f64) -> Result<Vec<Value>> {
    remove_cut_ties(model, cells)?;
    let rank = |cell: &str| (tet.contains(cell), cell.to_string());
    let mut out = Vec::new();
    // Per secondary cell: its cut faces and the facing cells' ones, as (body, face).
    let mut tied: BTreeMap<String, (BTreeSet<(String, String)>, BTreeSet<(String, String)>)> = BTreeMap::new();
    let mut tied_pairs: Vec<Value> = Vec::new();
    for (partition, [a_sides, b_sides]) in sides(model, cells) {
        for a in &a_sides {
            for b in &b_sides {
                // Faces on one plane overlap when their boxes share an area, not just an edge.
                let shared = (0..3).filter(|&i| a.lo[i].max(b.lo[i]) < a.hi[i].min(b.hi[i]) - tol).count();
                let touching = (0..3).all(|i| a.lo[i] <= b.hi[i] + tol && b.lo[i] <= a.hi[i] + tol);
                if !(touching && shared >= 2) {
                    continue;
                }
                if let Some(map) = matching(model, a, b, tol) {
                    let count = map.len();
                    model.merge_nodes(&map)?;
                    out.push(json!({"partition": partition, "cells": [a.body, b.body], "joined": "merged", "nodes": count}));
                    continue;
                }
                let (main, secondary) = if rank(&a.body) <= rank(&b.body) { (a, b) } else { (b, a) };
                let entry = tied.entry(secondary.body.clone()).or_default();
                entry.0.insert((secondary.body.clone(), secondary.face.clone()));
                entry.1.insert((main.body.clone(), main.face.clone()));
                tied_pairs.push(json!({
                    "partition": partition,
                    "cells": [a.body, b.body],
                    "joined": "tied",
                    "main": main.body,
                    "secondary": secondary.body,
                    "faces": [a.faces.len(), b.faces.len()],
                }));
            }
        }
    }
    let mut names: BTreeMap<String, String> = BTreeMap::new();
    for (cell, (own, facing)) in tied {
        let base = tie_name(&cell);
        let mut name = base.clone();
        let mut k = 1;
        while model.tie(&name).is_some() || model.set(&format!("{name}_main")).is_some() || model.set(&format!("{name}_secondary")).is_some() {
            k += 1;
            name = format!("{base}_{k}");
        }
        let set = |faces: &BTreeSet<(String, String)>| -> Result<Set> {
            Ok(Set::of_geometry(faces.iter().map(|(b, f)| GeometryRef::face(b.as_str(), f.as_str())).collect::<fieldcraft_model::Result<_>>()?))
        };
        let (main_set, secondary_set) = (format!("{name}_main"), format!("{name}_secondary"));
        model.add_set(&main_set, set(&facing)?)?;
        model.add_set(&secondary_set, set(&own)?)?;
        model.add_tie(&name, Tie::new(&main_set, &secondary_set))?;
        names.insert(cell, name);
    }
    for mut pair in tied_pairs {
        let name = pair.get("secondary").and_then(Value::as_str).and_then(|c| names.get(c)).cloned();
        if let (Some(o), Some(name)) = (pair.as_object_mut(), name) {
            o.insert("tie".into(), json!(name));
        }
        out.push(pair);
    }
    Ok(out)
}

/// Give every tied interface between cells of one design body its reason (`why`), from the
/// cells' reports (#396): a cell meshed with tet (and why auto fell back), the reason the body's
/// cells did not sweep together, or a cell whose method does not sweep with the others; so an
/// agent knows which cut or method to change for a conforming mesh.
fn explain_ties(design: &fieldcraft_geometry::Design, report: &[Value], interfaces: &mut [Value]) {
    let info = |cell: &str| report.iter().find(|r| r.get("body").and_then(Value::as_str) == Some(cell));
    let text = |r: Option<&Value>, key: &str| r.and_then(|r| r.get(key)).and_then(Value::as_str).map(String::from);
    for interface in interfaces.iter_mut() {
        if interface.get("joined").and_then(Value::as_str) != Some("tied") {
            continue;
        }
        let cells: Vec<String> = interface
            .get("cells")
            .and_then(Value::as_array)
            .map(|c| c.iter().filter_map(Value::as_str).map(String::from).collect())
            .unwrap_or_default();
        let [a, b] = cells.as_slice() else { continue };
        let Some(source) = design.source_of(a).filter(|s| design.source_of(b) == Some(*s)) else { continue };
        let mut why: Vec<String> = Vec::new();
        for cell in [a, b] {
            let r = info(cell);
            if text(r, "kind").is_some_and(|k| k.starts_with("tet")) {
                why.push(match text(r, "why") {
                    Some(w) => format!("`{cell}` is meshed with tet: {w}"),
                    None => format!("`{cell}` is meshed with tet"),
                });
            }
        }
        if why.is_empty() {
            let not_together = [a, b].into_iter().find_map(|c| text(info(c), "not_swept_together"));
            let alone = [a, b].into_iter().find(|c| text(info(c), "swept_together").is_none() && text(info(c), "not_swept_together").is_none());
            if let Some(w) = not_together {
                why.push(format!("the cells of `{source}` do not sweep together ({w}), so each was meshed alone with its own divisions"));
            } else if let Some(c) = alone {
                let method = text(info(c), "method").unwrap_or_default();
                why.push(format!("`{c}` is meshed alone (method `{method}`), not swept with the other cells of `{source}`"));
            } else {
                why.push(format!("the element faces of `{a}` and `{b}` do not match node for node"));
            }
        }
        if let Some(o) = interface.as_object_mut() {
            o.insert("why".into(), json!(why.join("; ")));
        }
    }
}

/// The name of the tie that holds `cell`'s tied cut faces.
fn tie_name(cell: &str) -> String {
    let name = format!("tie_{cell}");
    fieldcraft_model::suggest_name(&name).unwrap_or(name)
}

/// Remove the ties an earlier run of mesh.cells made on `cells` (and their two sets): ties
/// between sets `<tie>_main` and `<tie>_secondary` that hold only cut faces of these cells.
fn remove_cut_ties(model: &mut Model, cells: &[String]) -> Result<()> {
    let cut_faces_of_cells = |set: &str| {
        model.set(set).is_some_and(|s| {
            !s.geometry.is_empty()
                && s.nodes.is_empty()
                && s.elements.is_empty()
                && s.faces.is_empty()
                && s.edges.is_empty()
                && s.geometry.iter().all(|g| cells.contains(&g.body) && matches!(&g.entity, Entity::Face(f) if f.starts_with("partition:")))
        })
    };
    let stale: Vec<(String, String, String)> = model
        .ties()
        .iter()
        .filter(|(n, t)| t.main == format!("{n}_main") && t.secondary == format!("{n}_secondary"))
        .filter(|(_, t)| cut_faces_of_cells(&t.main) && cut_faces_of_cells(&t.secondary))
        .map(|(n, t)| (n.clone(), t.main.clone(), t.secondary.clone()))
        .collect();
    for (tie, main, secondary) in stale {
        model.remove_tie(&tie)?;
        for set in [main, secondary] {
            if model.set_users(&set).is_empty() {
                model.remove_set(&set)?;
            }
        }
    }
    Ok(())
}

/// When the two sides' nodes coincide one to one (within `tol`) and their faces are the same
/// in number, the merge: each node of side `b` to the node of side `a` at its place (or the other
/// way round, so the lower label survives).
fn matching(model: &Model, a: &Side, b: &Side, tol: f64) -> Option<BTreeMap<u32, u32>> {
    if a.nodes.len() != b.nodes.len() || a.faces.len() != b.faces.len() {
        return None;
    }
    let both: BTreeSet<u32> = a.nodes.union(&b.nodes).copied().collect();
    let map = fieldcraft_mesh::join::coincident_nodes(model, Some(&both), tol).ok()?;
    // Every pair joins one node of each side, and every node of either side is in a pair.
    let one_to_one = map.len() == a.nodes.len()
        && map.iter().all(|(from, to)| (a.nodes.contains(from) && b.nodes.contains(to)) || (b.nodes.contains(from) && a.nodes.contains(to)));
    one_to_one.then_some(map)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sample plate (40 × 30 × 10, a Ø10 hole at (15, 15)) open and cut at x = 30, right of
    /// the hole: `Body1/a` holds the hole, `Body1/b` is a 10 × 30 × 10 block.
    fn cut_plate() -> (crate::tests::TempDir, Session) {
        let (dir, path) = crate::tests::sample_design_file(&fieldcraft_geometry::sample::Plate::default());
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        s.execute("geometry.datum", &json!({"name": "X30", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "P", "datum": "X30"})).unwrap();
        (dir, s)
    }

    #[test]
    fn hex_and_tet_cells_are_tied_and_matching_hex_cells_merged() {
        let (_dir, mut s) = cut_plate();
        // A tet cell next to a hex cell: tied, the hex cell main.
        let r = s.execute("mesh.cells", &json!({"methods": {"Body1/a": "tet", "Body1/b": "mapped"}, "size": 2.5, "transition": "tie"})).unwrap();
        assert_eq!(r["cells"][0]["kind"], "tet4", "{r}");
        assert_eq!(r["cells"][1]["kind"], "hex8", "{r}");
        let i = &r["interfaces"][0];
        assert_eq!(i["joined"], "tied", "{r}");
        assert!(i["why"].as_str().is_some_and(|w| w.contains("`Body1/a` is meshed with tet")), "{r}");
        assert_eq!(s.model.ties().len(), 1);
        assert_eq!((i["main"].as_str(), i["secondary"].as_str(), i["tie"].as_str()), (Some("Body1/b"), Some("Body1/a"), Some("tie_Body1_a")), "{r}");
        for set in ["tie_Body1_a_main", "tie_Body1_a_secondary"] {
            assert!(s.model.resolve_set(set).unwrap().faces.len() > 4, "{set}");
        }
        // Both cells hex at one size: the plane faces match node for node and merge.
        let r = s.execute("mesh.cells", &json!({"default": "hex", "size": 2.5})).unwrap();
        assert_eq!(r["interfaces"].as_array().map(Vec::len), Some(1), "{r}");
        let i = &r["interfaces"][0];
        assert_eq!(i["joined"], "merged", "{r}");
        let separate = r["cells"].as_array().unwrap().iter().map(|c| c["nodes"].as_u64().unwrap()).sum::<u64>();
        assert_eq!(r["nodes"].as_u64().unwrap() + i["nodes"].as_u64().unwrap(), separate, "{r}");
        // Rerunning replaces the tie it made rather than adding another.
        s.execute("mesh.cells", &json!({"methods": {"Body1/a": "tet"}, "size": 2.5, "transition": "tie"})).unwrap();
        assert_eq!(s.model.ties().len(), 1);
    }

    /// A tie between cells of one design body says why they were not conforming (#396): the
    /// reason the body's cells did not sweep together, or the cell meshed alone by its method.
    #[test]
    fn ties_inside_one_body_say_why() {
        let (_dir, s) = cut_plate();
        let design = s.design.clone().unwrap();
        let tied = || vec![json!({"cells": ["Body1/a", "Body1/b"], "joined": "tied"}), json!({"cells": ["Body1/a", "Body1/b"], "joined": "merged"})];
        let mut interfaces = tied();
        let report = [
            json!({"body": "Body1/a", "kind": "hex8", "method": "extrude", "not_swept_together": "extrude: the levels differ"}),
            json!({"body": "Body1/b", "kind": "hex8", "method": "extrude", "not_swept_together": "extrude: the levels differ"}),
        ];
        explain_ties(&design, &report, &mut interfaces);
        assert_eq!(
            interfaces[0]["why"],
            "the cells of `Body1` do not sweep together (extrude: the levels differ), so each was meshed alone with its own divisions"
        );
        assert!(interfaces[1].get("why").is_none());
        let mut interfaces = tied();
        let report = [
            json!({"body": "Body1/a", "kind": "hex8", "method": "mapped"}),
            json!({"body": "Body1/b", "kind": "hex8", "method": "extrude", "swept_together": "Body1"}),
        ];
        explain_ties(&design, &report, &mut interfaces);
        assert_eq!(interfaces[0]["why"], "`Body1/a` is meshed alone (method `mapped`), not swept with the other cells of `Body1`");
    }

    /// The half crank arm (cut on its symmetry plane only) has no cell that maps or sweeps as
    /// the design gives it; the extrusion decomposes its cross-section into blocks itself, alike
    /// through mesh.hex.generate (auto) and mesh.cells, with no seeds.
    #[test]
    fn the_half_crank_extrudes_in_blocks_without_partitions() {
        let bytes = fieldcraft_geometry::script_design(crate::sample::CRANK_ARM).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        s.execute("geometry.datum", &json!({"name": "Sym", "by": {"type": "offset", "plane": "zx", "offset": 0}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "Half", "datum": "Sym"})).unwrap();
        let r = s.execute("mesh.hex.generate", &json!({"body": "CrankArm/a", "size": 3})).unwrap();
        let body = &r["bodies"][0];
        assert_eq!(body["method"], "extrude", "{r}");
        assert!(body["description"].as_str().unwrap().contains(" 0 free)"), "{r}");
        assert!(body["quality"]["min_scaled_jacobian"].as_f64().unwrap() >= 0.37, "{r}");
        let r = s.execute("mesh.cells", &json!({"bodies": ["CrankArm/a"], "default": "extrude", "size": 3})).unwrap();
        assert_eq!(r["cells"][0]["method"], "extrude", "{r}");
        assert_eq!(r["elements"], body["elements"], "{r}");
    }

    #[test]
    fn auto_falls_back_to_tet_and_says_why() {
        let (_dir, mut s) = cut_plate();
        let r = s.execute("mesh.cells", &json!({"default": "mapped", "methods": {"Body1/a": "auto"}, "size": 3})).unwrap();
        // The holed cell is not a mapped block, but auto sweeps it.
        assert!(r["cells"][0]["method"].as_str().is_some(), "{r}");
        assert_eq!(r["cells"][1]["method"], "mapped", "{r}");
    }

    #[test]
    fn hostile_params_change_nothing() {
        let (_dir, mut s) = cut_plate();
        let journal = s.journal.len();
        for p in [
            json!({"default": "pyramid"}),
            json!({"methods": {"Nope": "tet"}}),
            json!({"methods": ["tet"]}),
            json!({"bodies": []}),
            json!({"bodies": ["Body1/a", "Body1/a"]}),
            json!({"bodies": ["Nope"]}),
            json!({"order": 3}),
            json!({"size": -1}),
            json!({"default": "mapped", "bodies": ["Body1/a"]}),
        ] {
            assert!(s.execute("mesh.cells", &p).is_err(), "{p}");
        }
        assert_eq!(s.journal.len(), journal);
        assert!(s.model.elements().is_empty());
    }

    /// A 40 × 10 × 10 bar in uniaxial tension, cut in the middle: one half hex, the other tet,
    /// joined by `transition` (`joined` in the reply). With ν = 0, symmetry conditions on three
    /// faces and the end load spread by consistent weights, the stress is uniform, which the
    /// elements, the pyramids and the dual mortar tie all carry exactly: the end moves F L / (E A).
    fn hex_tet_bar_in_tension(transition: &str, joined: &str) -> Session {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 40, "width": 10, "height": 10, "body_name": "Bar"}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let bbox = s.execute("geometry.bodies", &json!({})).unwrap()[0]["bbox"].clone();
        let lo: Vec<f64> = bbox[0].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        let hi: Vec<f64> = bbox[1].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
        let (axis, length) = (0..3).map(|k| (k, hi[k] - lo[k])).max_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        assert!((length - 40.0).abs() < 1e-9);
        let mid = (lo[axis] + hi[axis]) / 2.0;
        let mut normal = [0.0; 3];
        normal[axis] = 1.0;
        let mut point = [0.0; 3];
        point[axis] = mid;
        s.execute("geometry.datum", &json!({"name": "Mid", "by": {"type": "point_normal", "point": point, "normal": normal}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "Cut", "datum": "Mid"})).unwrap();
        let r = s.execute("mesh.cells", &json!({"methods": {"Bar/a": "mapped", "Bar/b": "tet"}, "size": 2.5, "transition": transition})).unwrap();
        assert_eq!(r["interfaces"][0]["joined"], joined, "{r}");

        let m = &s.model;
        // The model is stored in metres, the geometry's bounding box in millimetres.
        let mm = crate::units::to_geometry_length;
        let on = |k: usize, v: f64| -> Vec<u32> { m.nodes().iter().filter(|(_, n)| (mm(n.xyz[k]) - v).abs() < 1e-9).map(|(l, _)| *l).collect() };
        let others: Vec<usize> = (0..3).filter(|k| *k != axis).collect();
        let end: Vec<Value> = m
            .elements()
            .iter()
            .flat_map(|(&e, el)| (1..=el.kind.faces().len() as u8).map(move |f| (e, f)))
            .filter(|&(e, f)| {
                m.face_nodes(ElementFace { element: e, face: f })
                    .is_some_and(|ns| ns.iter().all(|n| (mm(m.node(*n).unwrap().xyz[axis]) - hi[axis]).abs() < 1e-9))
            })
            .map(|(e, f)| json!({"element": e, "face": f}))
            .collect();
        let plane = |k: usize| ["x", "y", "z"][k];
        let (young, area, force) = (70_000.0, 100.0, 1000.0);
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "elements": m.elements().keys().collect::<Vec<_>>()}},
            {"command": "model.set.create", "params": {"name": "root", "nodes": on(axis, lo[axis])}},
            {"command": "model.set.create", "params": {"name": "side1", "nodes": on(others[0], lo[others[0]])}},
            {"command": "model.set.create", "params": {"name": "side2", "nodes": on(others[1], lo[others[1]])}},
            {"command": "model.set.create", "params": {"name": "end", "faces": end}},
            {"command": "model.material.create", "params": {"name": "Al", "elastic": {"type": "isotropic", "young": young, "poisson": 0}}},
            {"command": "model.section.create", "params": {"name": "Bar", "set": "all", "material": "Al", "kind": {"type": "solid"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "s0", "type": "symmetry", "target": {"set": "root"}, "plane": plane(axis)}},
            {"command": "model.boundary_condition.create", "params": {"name": "s1", "type": "symmetry", "target": {"set": "side1"}, "plane": plane(others[0])}},
            {"command": "model.boundary_condition.create", "params": {"name": "s2", "type": "symmetry", "target": {"set": "side2"}, "plane": plane(others[1])}},
            {"command": "model.load.create", "params": {"name": "Pull", "type": "force", "target": {"set": "end"}, "vector": normal.map(|n| n * force)}},
            {"command": "model.step.create", "params": {"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["s0", "s1", "s2"], "loads": ["Pull"]}},
        ]);
        let tip = on(axis, hi[axis]);
        s.run_script(&script).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        let probe = s.execute("solve.probe", &json!({"step": "Load", "nodes": tip})).unwrap();
        let expected = force * length / (young * area);
        let probed = probe["nodes"].as_array().unwrap();
        assert!(probed.len() > 4, "{probe}");
        for n in probed {
            let u = n["displacement"][axis].as_f64().unwrap();
            assert!((u - expected).abs() < 1e-6 * expected, "{u} vs {expected}: {n}");
        }
        s
    }

    #[test]
    fn a_tied_hex_tet_bar_in_tension_stretches_as_the_formula_says() {
        let s = hex_tet_bar_in_tension("tie", "tied");
        assert_eq!(s.model.ties().len(), 1);
    }

    /// Through pyramids the bar is one conforming mesh: no tie, no free faces inside, every
    /// element valid.
    #[test]
    fn a_hex_tet_bar_joined_by_pyramids_is_conforming_and_stretches_as_the_formula_says() {
        let mut s = hex_tet_bar_in_tension("pyramid", "pyramids");
        assert!(s.model.ties().is_empty());
        let kinds: BTreeSet<&str> = s.model.elements().values().map(|e| e.kind.name()).collect();
        assert_eq!(kinds.into_iter().collect::<Vec<_>>(), ["hex8", "pyramid5", "tet4"]);
        let d = s.execute("mesh.check.duplicates", &json!({})).unwrap();
        assert_eq!(d["counts"], json!({"coincident": 0, "cracks": 0, "duplicate_elements": 0, "free_edges": 0}), "{d}");
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!(q["invalid"]["count"], 0, "{q}");
    }

    /// A box cut by planes into cells, `methods` per cell, a heat source in every cell and
    /// convection on its outer faces: the heat flows of the steady solve.
    fn cut_box_heat(cuts: &[(&str, f64)], methods: Value, size: f64) -> (Value, Value) {
        let commands = json!({"commands": [
            {"command": "solid.box", "params": {"length": 30, "width": 20, "height": 10, "body_name": "Box"}},
        ]});
        let bytes = fieldcraft_geometry::script_design(&commands.to_string()).unwrap();
        let (_dir, path) = crate::tests::design_file(&bytes);
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        for (i, (plane, offset)) in cuts.iter().enumerate() {
            let name = format!("D{i}");
            s.execute("geometry.datum", &json!({"name": name, "by": {"type": "offset", "plane": plane, "offset": offset}})).unwrap();
            s.execute("geometry.partition", &json!({"name": format!("P{i}"), "datum": name})).unwrap();
        }
        let r = s.execute("mesh.cells", &json!({"methods": methods, "default": "hex", "size": size})).unwrap();
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let outer: Vec<Value> = bodies
            .as_array()
            .unwrap()
            .iter()
            .flat_map(|b| {
                b["faces"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .filter(|f| !f["name"].as_str().unwrap().starts_with("partition:"))
                    .map(|f| json!({"body": b["name"], "face": f["name"]}))
                    .collect::<Vec<_>>()
            })
            .collect();
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "elements": s.model.elements().keys().collect::<Vec<_>>()}},
            {"command": "model.set.create", "params": {"name": "skin", "geometry": outer}},
            {"command": "model.material.create", "params": {"name": "M", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0.3}, "thermal": {"conductivity": 0.5}}},
            {"command": "model.section.create", "params": {"name": "S", "set": "all", "material": "M", "kind": {"type": "solid"}}},
            {"command": "model.load.create", "params": {"name": "Q", "type": "heat_source", "set": "all", "power": 1}},
            {"command": "model.load.create", "params": {"name": "Air", "type": "convection", "set": "skin", "film": 0.01, "sink": 25}},
            {"command": "model.step.create", "params": {"name": "Steady", "procedure": {"type": "heat_transfer"}, "loads": ["Q", "Air"]}},
        ]);
        s.run_script(&script).unwrap();
        let h = s.execute("solve.heat", &json!({"step": "Steady"})).unwrap();
        (r, h)
    }

    /// A tet cell among hex cells (#389, #288): the hex cells sweep together and share their
    /// nodes, the tet cell is tied once (one tie for all its cut faces, the hex cells main), and
    /// the heat balances to round-off.
    #[test]
    fn a_tet_cell_among_swept_hex_cells_is_tied_once_and_the_heat_balances() {
        let (r, h) = cut_box_heat(&[("yz", 10.0), ("yz", 20.0), ("zx", 10.0)], json!({"Box/b/a/a": "tet"}), 2.5);
        let interfaces = r["interfaces"].as_array().unwrap();
        let tied: Vec<&Value> = interfaces.iter().filter(|i| i["joined"] == "tied").collect();
        assert_eq!(tied.len(), 3, "{r}");
        assert!(tied.iter().all(|i| i["secondary"] == "Box/b/a/a" && i["tie"] == "tie_Box_b_a_a"), "{r}");
        assert!(tied.iter().all(|i| i["why"].as_str().is_some_and(|w| w.starts_with("`Box/b/a/a` is meshed with tet"))), "{r}");
        let merged: Vec<&Value> = interfaces.iter().filter(|i| i["joined"] == "merged").collect();
        assert_eq!(merged.len(), 4, "{r}");
        // Box/a/a and Box/a/b meet across y = 10, the third cut.
        assert!(merged.iter().any(|i| i["cells"] == json!(["Box/a/a", "Box/a/b"]) && i["partition"] == "P2"), "{r}");
        let flows = &h["heat_flows"];
        let supplied = flows["loads"].as_array().unwrap().iter().find(|l| l["name"] == "Q").unwrap()["heat"].as_f64().unwrap();
        assert!(flows["imbalance"].as_f64().unwrap().abs() < 1e-9 * supplied, "{flows}");
        assert!(h["min_temperature"]["value"].as_f64().unwrap() > 25.0, "{h}");
    }

    /// Cells cut from one box by planes, one of them (`y = 10`) cutting only the cell left of
    /// `x = 10`, so a cut ends on another (a T-junction): with no seeds the cells sweep together
    /// and every interface merges node for node (#288), and the mesh fills the box.
    #[test]
    fn cells_of_a_body_sweep_together_and_merge_without_seeds() {
        let (_dir, mut s) = super::super::datum::tests::scripted(json!([
            {"command": "solid.box", "params": {"length": 30, "width": 20, "height": 10, "corner": [0, 0, 0], "body_name": "Box"}},
        ]));
        super::super::datum::tests::cut_at(&mut s, &[(0, 10.0), (0, 21.0)]);
        s.execute("geometry.datum", &json!({"name": "Y", "by": {"type": "offset", "plane": "zx", "offset": 7.5}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "PY", "datum": "Y", "bodies": ["Box/a"]})).unwrap();
        let r = s.execute("mesh.cells", &json!({"default": "auto", "size": 2.0})).unwrap();
        let cells = r["cells"].as_array().unwrap();
        assert_eq!(cells.len(), 4, "{r}");
        assert!(cells.iter().all(|c| c["swept_together"] == "Box" && c["kind"] == "hex8"), "{r}");
        let interfaces = r["interfaces"].as_array().unwrap();
        assert_eq!(interfaces.len(), 4, "{r}");
        assert!(interfaces.iter().all(|i| i["joined"] == "merged"), "{r}");
        assert!(s.model.ties().is_empty());
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        let volume = q["volume"].as_f64().unwrap();
        assert!((volume - 6000.0).abs() < 1e-6, "{q}");
        let d = s.execute("mesh.check.duplicates", &json!({})).unwrap();
        assert_eq!((d["counts"]["coincident"].as_u64(), d["counts"]["cracks"].as_u64()), (Some(0), Some(0)), "{d}");
    }
}
