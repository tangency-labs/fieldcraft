//! `mesh.block`: a structured mesh of one six-sided block (solid elements) or four-sided planar
//! patch (2D elements) from its corners and curved edges ([`fieldcraft_mesh::block`]), joined to
//! the mesh already there and tagged as a body (a patch: a face of the body) with named faces,
//! edges and vertices, so sets by geometry select it.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_mesh::block::{Block, BlockMesh, Curve, Patch, block, patch};
use fieldcraft_mesh::join::{check_conforming, weld};
use fieldcraft_model::{ElementKind, Entity, GeometryRef, Model};
use serde::Deserialize;
use serde_json::{Value, json};

use super::fresh_labels;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// The default merge distance, relative to the block's diagonal.
const DEFAULT_RELATIVE_MERGE: f64 = 1e-6; // unit-factors: relative to the block's diagonal

pub const BLOCK: CommandSpec = CommandSpec {
    id: "mesh.block",
    label: "Mesh Block",
    params: r#"{"body": "Plate", "kind": "hex8|hex20|tet4|tet10|wedge6", "corners": [[x, y, z] × 8] (hex8 order: 1-2-3-4 counter-clockwise seen from 5-6-7-8), "divisions": [n12, n14, n15] (elements along 1→2, 1→4, 1→5), "bias": [1, 1, 1] (optional: last over first element size along each direction), "edges": [{"corners": [1, 2], "arc": {"center": [x, y, z]}} | {"corners": [2, 3], "ellipse": {"center": [x, y, z], "axes": [[a], [b]] (perpendicular)}} | {"corners": [5, 6], "name": "tip_edge"}] (optional: curved edges, the shorter arc between the corners, and edge names), "faces": {"1": "bottom", "4": "tip"} (optional: names of block faces, numbered as hex8 faces: 1 = 1-2-3-4, 2 = 5-6-7-8, 3 = 1-2-6-5, 4 = 2-3-7-6, 5 = 3-4-8-7, 6 = 4-1-5-8), "vertices": {"1": "A"} (optional: names of corners), "merge": distance (optional, default 1e-6 of the block's diagonal; 0 never joins)} — adds a mapped mesh of the block (transfinite interpolation of its edges; quadratic nodes on the curves), uses existing nodes within `merge` of new ones so blocks that share a face join (the shared face must have the same elements on both sides, or nothing changes), and tags the elements with the body and element faces, edges and nodes with the named faces, edges and vertices. 2D: "kind": "quad4|quad8|tri3|tri6", four "corners" at z = 0 (counter-clockwise from +z), "divisions": [n12, n14], "bias": [r, r], "face": "area" (the face of the body the elements fill, default `area`), edges and vertices numbered among the four corners; for plane stress, plane strain and axisymmetric sections"#,
    form: Some(super::forms::BLOCK),
    kind: CommandKind::Model,
    requires: &[],
    run: block_mesh,
};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Params {
    body: String,
    kind: String,
    corners: Vec<[f64; 3]>,
    divisions: Vec<usize>,
    bias: Option<Vec<f64>>,
    #[serde(default)]
    edges: Vec<EdgeParams>,
    #[serde(default)]
    faces: BTreeMap<String, String>,
    face: Option<String>,
    #[serde(default)]
    vertices: BTreeMap<String, String>,
    merge: Option<f64>,
}

/// A list parameter of exactly `N` entries.
fn exactly<T: Copy, const N: usize>(v: &[T], what: &str, dimension: u8) -> Result<[T; N]> {
    <[T; N]>::try_from(v).map_err(|_| EngineError::Params(format!("a {dimension}D block needs {N} {what}, got {}", v.len())))
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct EdgeParams {
    corners: [usize; 2],
    name: Option<String>,
    arc: Option<ArcParams>,
    ellipse: Option<EllipseParams>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ArcParams {
    center: [f64; 3],
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct EllipseParams {
    center: [f64; 3],
    axes: [[f64; 3]; 2],
}

/// A 1-based number in `1..=n` given as a JSON object key, 0-based.
fn numbered(key: &str, n: usize, what: &str) -> Result<usize> {
    key.parse::<usize>()
        .ok()
        .filter(|k| (1..=n).contains(k))
        .map(|k| k - 1)
        .ok_or_else(|| EngineError::Params(format!("{what} `{key}` must be a number from 1 to {n}")))
}

fn block_mesh(s: &mut Session, p: &Value) -> Result<Value> {
    let p: Params = serde_json::from_value(p.clone()).map_err(|e| EngineError::Params(format!("mesh.block: {e}")))?;
    let kind = ElementKind::from_name(&p.kind).ok_or_else(|| EngineError::Params(format!("unknown element kind `{}`", p.kind)))?;
    let planar = kind.dimension() == 2;
    let (corner_count, dimension) = if planar { (4, 2) } else { (8, 3) };
    let edge_index = |pair: [usize; 2]| if planar { Patch::edge_index(pair) } else { Block::edge_index(pair) };
    let mut curves = Vec::new();
    let mut edge_names = Vec::new();
    for e in &p.edges {
        let pair = [0, 1].map(|i| e.corners[i].wrapping_sub(1));
        let index = edge_index(pair).ok_or_else(|| EngineError::Params(format!("corners {:?} are not the two ends of an edge", e.corners)))?;
        let curve = match (&e.arc, &e.ellipse) {
            (Some(a), None) => Some(Curve::Arc { center: a.center }),
            (None, Some(el)) => Some(Curve::Ellipse { center: el.center, axes: el.axes }),
            (None, None) => None,
            (Some(_), Some(_)) => return Err(EngineError::Params(format!("edge {:?} is an arc or an ellipse, not both", e.corners))),
        };
        if let Some(c) = curve {
            curves.push((pair, c));
        }
        if let Some(n) = &e.name {
            edge_names.push((index, n.clone()));
        }
        if curve.is_none() && e.name.is_none() {
            return Err(EngineError::Params(format!("edge {:?} gives neither a curve nor a name", e.corners)));
        }
    }
    // Check every name before meshing.
    let region = match (planar, &p.face) {
        (true, face) => GeometryRef::face(p.body.clone(), face.clone().unwrap_or_else(|| "area".into()))?,
        (false, None) => GeometryRef::body(p.body.clone())?,
        (false, Some(_)) => return Err(EngineError::invalid("face", "names the region of 2D elements; a solid block fills its body")),
    };
    if planar && !p.faces.is_empty() {
        return Err(EngineError::Params("a 2D patch has no faces to name: name its edges".into()));
    }
    let named = |entries: &BTreeMap<String, String>, n: usize, what: &str, entity: fn(String) -> Entity| -> Result<Vec<(usize, GeometryRef)>> {
        entries.iter().map(|(k, name)| Ok((numbered(k, n, what)?, GeometryRef::new(p.body.clone(), entity(name.clone()))?))).collect()
    };
    let face_refs = named(&p.faces, 6, "face", Entity::Face)?;
    let vertex_refs = named(&p.vertices, corner_count, "vertex", Entity::Vertex)?;
    let edge_refs: Vec<(usize, GeometryRef)> =
        edge_names.into_iter().map(|(k, n)| Ok((k, GeometryRef::new(p.body.clone(), Entity::Edge(n))?))).collect::<Result<_>>()?;

    let bias = p.bias.clone().unwrap_or_else(|| vec![1.0; dimension as usize]);
    let mesh: BlockMesh = if planar {
        let corners = exactly(&p.corners, "corners", dimension)?;
        patch(&Patch {
            kind,
            corners,
            divisions: exactly(&p.divisions, "divisions", dimension)?,
            bias: exactly(&bias, "bias values", dimension)?,
            curves,
        })?
    } else {
        let corners = exactly(&p.corners, "corners", dimension)?;
        block(&Block {
            kind,
            corners,
            divisions: exactly(&p.divisions, "divisions", dimension)?,
            bias: exactly(&bias, "bias values", dimension)?,
            curves,
        })?
    };
    let diagonal = {
        let lo = (0..3).map(|d| p.corners.iter().map(|c| c[d]).fold(f64::INFINITY, f64::min));
        let hi = (0..3).map(|d| p.corners.iter().map(|c| c[d]).fold(f64::NEG_INFINITY, f64::max));
        lo.zip(hi).map(|(l, h)| (h - l) * (h - l)).sum::<f64>().sqrt()
    };
    let merge = p.merge.unwrap_or(DEFAULT_RELATIVE_MERGE * diagonal);
    if !(merge.is_finite() && merge >= 0.0 && merge <= 1e-2 * diagonal) {
        return Err(EngineError::Params(format!("`merge` must be from 0 to 1/100 of the block's diagonal ({diagonal})")));
    }
    let (first_node, first_element) = fresh_labels(&s.model, mesh.points.len(), mesh.elements.len())?;

    // Each block point: an existing node within `merge`, or a new label.
    let welded = weld(&s.model, &mesh.points, merge)?;
    let mut labels = Vec::with_capacity(mesh.points.len());
    let mut added = 0u32;
    for (xyz, existing) in mesh.points.iter().zip(&welded) {
        let label = match existing {
            Some(label) => *label,
            None => {
                added += 1;
                s.model.set_node(first_node + added - 1, *xyz)?;
                first_node + added - 1
            }
        };
        labels.push(label);
    }
    let joined: BTreeSet<u32> = welded.iter().flatten().copied().collect();
    let label_of = |i: u32| labels.get(i as usize).copied().ok_or_else(|| EngineError::Params(format!("block point {i} has no node")));
    let elements: Vec<(ElementKind, Vec<u32>)> =
        mesh.elements.iter().map(|e| Ok((kind, e.iter().map(|&i| label_of(i)).collect::<Result<_>>()?))).collect::<Result<_>>()?;
    check_conforming(&s.model, &elements, &joined, merge)?;

    let body_tag = s.model.geometry_tag(&region)?;
    for (label, (_, nodes)) in (first_element..).zip(elements) {
        s.model.set_element(label, kind, nodes)?;
        s.model.tag_element(label, body_tag)?;
    }
    for (refs, pieces, tag) in [
        (&face_refs, &mesh.faces, Model::tag_element_face as fn(&mut Model, u32, u8, u32) -> fieldcraft_model::Result<()>),
        (&edge_refs, &mesh.edges, Model::tag_element_edge),
    ] {
        for (k, r) in refs {
            let t = s.model.geometry_tag(r)?;
            for &(e, local, _) in pieces.iter().filter(|piece| usize::from(piece.2) == *k) {
                tag(&mut s.model, first_element + e, local, t)?;
            }
        }
    }
    for (k, r) in &vertex_refs {
        let t = s.model.geometry_tag(r)?;
        let corner = mesh.corners.get(*k).copied().unwrap_or(u32::MAX);
        s.model.tag_node(label_of(corner)?, t)?;
    }
    Ok(json!({
        "body": p.body,
        "kind": kind.name(),
        "elements": mesh.elements.len(),
        "first_element": first_element,
        "nodes": added,
        "first_node": first_node,
        "joined": joined.len(),
    }))
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use crate::Session;

    fn cube(x0: f64) -> serde_json::Value {
        json!([[x0, 0, 0], [x0 + 1.0, 0, 0], [x0 + 1.0, 1, 0], [x0, 1, 0], [x0, 0, 1], [x0 + 1.0, 0, 1], [x0 + 1.0, 1, 1], [x0, 1, 1]])
    }

    #[test]
    fn blocks_join_face_to_face_and_name_their_geometry() {
        let mut s = Session::default();
        let r = s
            .execute(
                "mesh.block",
                &json!({"body": "Bar", "kind": "hex20", "corners": cube(0.0), "divisions": [2, 2, 2], "faces": {"6": "root"}, "vertices": {"7": "corner"}, "edges": [{"corners": [5, 8], "name": "top_root"}]}),
            )
            .unwrap();
        assert_eq!((r["elements"].clone(), r["nodes"].clone(), r["joined"].clone()), (json!(8), json!(81), json!(0)));
        let r = s
            .execute("mesh.block", &json!({"body": "Bar", "kind": "hex20", "corners": cube(1.0), "divisions": [2, 2, 2], "faces": {"4": "tip"}}))
            .unwrap();
        assert_eq!((r["first_element"].clone(), r["nodes"].clone(), r["joined"].clone()), (json!(9), json!(81 - 21), json!(21)));
        for (set, geometry, nodes) in [
            ("root", json!({"body": "Bar", "face": "root"}), 21),
            ("tip", json!({"body": "Bar", "face": "tip"}), 21),
            ("edge", json!({"body": "Bar", "edge": "top_root"}), 5),
            ("corner", json!({"body": "Bar", "vertex": "corner"}), 1),
        ] {
            s.execute("model.set.create", &json!({"name": set, "geometry": [geometry]})).unwrap();
            let r = s.execute("model.set.resolve", &json!({"name": set})).unwrap();
            assert_eq!(r["nodes"]["count"], nodes, "{set}: {r}");
        }
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        assert_eq!((q["valid"].clone(), q["volume"].as_f64().map(|v| (v - 2.0).abs() < 1e-12)), (json!(16), Some(true)));

        // A third block that meets the second with other divisions or another kind changes nothing.
        let model = s.model.clone();
        for bad in [
            json!({"body": "Bar", "kind": "hex20", "corners": cube(2.0), "divisions": [2, 3, 2]}),
            json!({"body": "Bar", "kind": "hex20", "corners": cube(2.0), "divisions": [1, 1, 1]}),
            json!({"body": "Bar", "kind": "hex8", "corners": cube(2.0), "divisions": [4, 4, 4]}),
            json!({"body": "Bar", "kind": "tet10", "corners": cube(2.0), "divisions": [2, 2, 2]}),
        ] {
            assert!(s.execute("mesh.block", &bad).is_err(), "{bad}");
            assert_eq!(s.model, model);
        }
        // A finer block first, then a coarser one against it: caught from the other side.
        let mut t = Session::default();
        t.execute("mesh.block", &json!({"body": "A", "kind": "hex8", "corners": cube(0.0), "divisions": [2, 2, 2]})).unwrap();
        assert!(t.execute("mesh.block", &json!({"body": "B", "kind": "hex8", "corners": cube(1.0), "divisions": [1, 1, 1]})).is_err());
        // Divisions that share no node inside the face: the hanging nodes on its edges are caught.
        assert!(t.execute("mesh.block", &json!({"body": "B", "kind": "hex8", "corners": cube(1.0), "divisions": [3, 3, 3]})).is_err());
        // A ring closed by its last block touches existing nodes at the corners of faces it does
        // not share (the ends): that is no conflict.
        let mut ring = Session::default();
        for k in 0..4 {
            let at = |r: f64, q: usize, z: f64| {
                let a = std::f64::consts::FRAC_PI_2 * q as f64;
                let c = |v: f64| if v.abs() < 1e-12 { 0.0 } else { v };
                json!([c(r * a.cos()), c(r * a.sin()), z])
            };
            let corners = json!([
                at(1.0, k, 0.0),
                at(2.0, k, 0.0),
                at(2.0, k + 1, 0.0),
                at(1.0, k + 1, 0.0),
                at(1.0, k, 1.0),
                at(2.0, k, 1.0),
                at(2.0, k + 1, 1.0),
                at(1.0, k + 1, 1.0)
            ]);
            let edges = json!([{"corners": [1, 4], "arc": {"center": [0, 0, 0]}}, {"corners": [2, 3], "arc": {"center": [0, 0, 0]}}, {"corners": [5, 8], "arc": {"center": [0, 0, 1]}}, {"corners": [6, 7], "arc": {"center": [0, 0, 1]}}]);
            ring.execute("mesh.block", &json!({"body": "Ring", "kind": "hex20", "corners": corners, "divisions": [2, 3, 2], "edges": edges}))
                .unwrap();
        }
        assert_eq!(ring.model.elements().len(), 48);
        assert!(t.execute("mesh.block", &json!({"body": "B", "kind": "hex8", "corners": cube(1.0), "divisions": [1, 1, 1], "merge": 0})).is_ok());
    }

    #[test]
    fn curved_edges_and_hostile_parameters() {
        let mut s = Session::default();
        let corners = json!([[1, 0, 0], [2, 0, 0], [0, 2, 0], [0, 1, 0], [1, 0, 1], [2, 0, 1], [0, 2, 1], [0, 1, 1]]);
        let edges = json!([
            {"corners": [2, 3], "arc": {"center": [0, 0, 0]}},
            {"corners": [6, 7], "arc": {"center": [0, 0, 1]}},
            {"corners": [1, 4], "ellipse": {"center": [0, 0, 0], "axes": [[1, 0, 0], [0, 1, 0]]}},
            {"corners": [5, 8], "arc": {"center": [0, 0, 1]}, "name": "inner_top"}
        ]);
        let r =
            s.execute("mesh.block", &json!({"body": "Ring", "kind": "tet10", "corners": corners, "divisions": [2, 4, 1], "edges": edges})).unwrap();
        assert_eq!(r["elements"], 48);
        let q = s.execute("mesh.quality", &json!({})).unwrap();
        let exact = std::f64::consts::PI * 3.0 / 4.0;
        assert!((q["volume"].as_f64().unwrap() - exact).abs() < 1e-3 * exact, "{q}");

        let journal = s.journal.len();
        for bad in [
            json!({}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1]}),
            json!({"body": "B", "kind": "quad4", "corners": cube(5.0), "divisions": [1, 1, 1]}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "faces": {"7": "x"}}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "faces": {"one": "x"}}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "vertices": {"0": "x"}}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "edges": [{"corners": [1, 7], "name": "x"}]}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "edges": [{"corners": [0, 1], "name": "x"}]}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "edges": [{"corners": [1, 2]}]}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "merge": -1}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "merge": 1}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1e9, 1e9, 1e9]}),
            json!({"body": "B", "kind": "hex8", "corners": cube(1e300), "divisions": [1, 1, 1]}),
            json!({"body": "", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1]}),
            json!({"body": "B", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "colour": "red"}),
        ] {
            assert!(s.execute("mesh.block", &bad).is_err(), "{bad}");
        }
        assert_eq!(s.journal.len(), journal);
    }

    #[test]
    fn patches_mesh_2d_regions_that_join_and_solve() {
        let mut s = Session::default();
        let o = [0, 0, 0];
        // Two quarter annuli of quad8 (a = 1, b = 2), side by side: a half annulus.
        for (corners, names, vertices) in [
            (
                json!([[1, 0, 0], [2, 0, 0], [0, 2, 0], [0, 1, 0]]),
                json!([{"corners": [1, 2], "name": "y0"}, {"corners": [1, 4], "arc": {"center": o}, "name": "bore"}]),
                json!({"3": "top"}),
            ),
            (
                json!([[0, 1, 0], [0, 2, 0], [-2, 0, 0], [-1, 0, 0]]),
                json!([{"corners": [1, 4], "arc": {"center": o}, "name": "bore"}, {"corners": [3, 4], "name": "y0"}]),
                json!({}),
            ),
        ] {
            let mut edges = names.as_array().unwrap().clone();
            edges.push(json!({"corners": [2, 3], "arc": {"center": o}}));
            s.execute(
                "mesh.block",
                &json!({"body": "Ring", "kind": "quad8", "corners": corners, "divisions": [4, 8], "edges": edges, "vertices": vertices}),
            )
            .unwrap();
        }
        assert_eq!(s.model.elements().len(), 64);
        for (cmd, params) in [
            ("model.set.create", json!({"name": "all", "geometry": [{"body": "Ring", "face": "area"}]})),
            ("model.set.create", json!({"name": "y0", "geometry": [{"body": "Ring", "edge": "y0"}]})),
            ("model.set.create", json!({"name": "bore", "geometry": [{"body": "Ring", "edge": "bore"}]})),
            ("model.set.create", json!({"name": "top", "geometry": [{"body": "Ring", "vertex": "top"}]})),
            ("model.boundary_condition.create", json!({"name": "hold", "type": "symmetry", "target": {"set": "top"}, "plane": "x"})),
            ("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 1000.0, "poisson": 0.3}})),
            ("model.section.create", json!({"name": "ring", "set": "all", "material": "Steel", "kind": {"type": "plane_strain", "thickness": 1.0}})),
            ("model.boundary_condition.create", json!({"name": "sym", "type": "symmetry", "target": {"set": "y0"}, "plane": "y"})),
            ("model.load.create", json!({"name": "p", "type": "pressure", "set": "bore", "magnitude": 10.0})),
            ("model.step.create", json!({"name": "Load", "procedure": {"type": "static"}, "boundary_conditions": ["sym", "hold"], "loads": ["p"]})),
        ] {
            s.execute(cmd, &params).unwrap_or_else(|e| panic!("{cmd}: {e}"));
        }
        s.execute("solve.static", &json!({})).unwrap();
        // Lamé, plane strain: u_r(a) = (1 + ν) a² p / (E (b² − a²)) ((1 − 2ν) a + b²/a).
        let exact = 1.3 * 10.0 / (1000.0 * 3.0) * (0.4 + 4.0);
        let r = s.execute("solve.probe", &json!({"step": "Load", "at": [[0, 1, 0]]})).unwrap();
        let u = r["nodes"][0]["displacement"][1].as_f64().unwrap();
        assert!((u - exact).abs() < 2e-3 * exact, "{u} vs {exact}");
        // 2D patches refuse solid parameters, and the other way round.
        for bad in [
            json!({"body": "P", "kind": "quad4", "corners": [[5, 0, 0], [6, 0, 0], [6, 1, 0], [5, 1, 0]], "divisions": [1, 1, 1]}),
            json!({"body": "P", "kind": "quad4", "corners": [[5, 0, 0], [6, 0, 0], [6, 1, 0], [5, 1, 0]], "divisions": [1, 1], "faces": {"1": "x"}}),
            json!({"body": "P", "kind": "quad4", "corners": [[5, 0, 1], [6, 0, 1], [6, 1, 1], [5, 1, 1]], "divisions": [1, 1]}),
            json!({"body": "P", "kind": "hex8", "corners": [[5, 0, 0], [6, 0, 0], [6, 1, 0], [5, 1, 0]], "divisions": [1, 1]}),
            json!({"body": "P", "kind": "hex8", "corners": cube(5.0), "divisions": [1, 1, 1], "face": "area"}),
        ] {
            assert!(s.execute("mesh.block", &bad).is_err(), "{bad}");
        }
    }
}
