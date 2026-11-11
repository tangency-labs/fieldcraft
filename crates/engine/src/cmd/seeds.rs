//! Seeds: `mesh.seed.edge` sets a count or a size (with bias and constraint) on design edges,
//! `mesh.seed.face` on every edge bounding design faces, `mesh.seed.clear` removes seeds and
//! `mesh.seed.list` lists them with the count each gives on its edges. Seeds live in the model
//! (`fieldcraft_model::seeds`), are saved with it, and the meshers honour them: a face seed is
//! resolved to its edges by the face's persistent name each time a mesher runs, and a seed on an
//! edge wins over a face's ([`edge_seeds`]).

use std::collections::BTreeMap;

use fieldcraft_geometry::{Design, Edge, Shape};
use fieldcraft_mesh::seeds::EdgeSeed;
use fieldcraft_model::{Bias, Divisions, Entity, GeometryRef, Seed, SeedConstraint};
use serde_json::{Value, json};

use super::opt;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// Most edges or faces one command names.
const MAX_REFS: usize = 10_000;

pub const SEED_EDGE: CommandSpec = CommandSpec {
    id: "mesh.seed.edge",
    label: "Seed Edges",
    params: r#"{"edges": [{"body": "Body1", "edge": "F2:side:hi|F2:top"}, …], "count": 8 | "size": 2.5, "bias"?: {"ratio": 4, "toward": "start" | "end" | "both"}, "constraint"?: "exact" (default) | "min" | "free"} — seeds design edges (persistent edge names, geometry.bodies): `count` divisions, or `size` (the edge's length over it, rounded up). `bias` grades the elements geometrically, the largest `ratio` times the smallest, small at the edge's start, end or both ends. `constraint`: "exact" keeps the count (a mesher whose other edges cannot match refuses, naming the conflicting edges), "min" keeps at least it, "free" makes it a goal. Replaces a seed the edge had; an edge name without its `~k` number seeds every edge of that name. A seed on an edge wins over the seeds of the faces it bounds (mesh.seed.face). The hex meshers (sweep, mapped, blocks, cells) use the counts and spacing; the tet mesher uses the spacing on the edge and grades its size field from it"#,
    form: Some(super::forms::SEED_EDGE),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: seed_edge,
};

pub const SEED_FACE: CommandSpec = CommandSpec {
    id: "mesh.seed.face",
    label: "Seed Faces",
    params: r#"{"faces": [{"body": "Body1", "face": "F2:top"}, …], "count": 8 | "size": 2.5, "bias"?: {"ratio": 4, "toward": "both"}, "constraint"?: "exact" (default) | "min" | "free"} — seeds every edge bounding design faces (persistent face names, geometry.bodies; a name without its `#k` number takes every piece) as mesh.seed.edge seeds one: `count` divisions on each edge, or `size` (each edge's length over it, rounded up); `bias` grades toward both ends of each edge only (an edge's start and end are its own curve's: seed the edge for a one-sided grading); `constraint` as mesh.seed.edge's. Resolved by the face's name each time a mesher runs, so it follows the design through edits that keep the name. A seed on an edge wins over a face's, also across bodies on an edge they share (a partition's cut); of two seeded faces sharing an edge, the one giving it more divisions wins. The reply and mesh.seed.list give each edge's count and, where another seed of its body decides it, `overridden_by` (that seed and its count). Replaces a seed the face had. A size control on a face (mesh.size.set) sets a graded element size over the face instead"#,
    form: Some(super::forms::SEED_FACE),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: seed_face,
};

pub const SEED_CLEAR: CommandSpec = CommandSpec {
    id: "mesh.seed.clear",
    label: "Clear Seeds",
    params: r#"{"edges"?: [{"body", "edge"}, …], "faces"?: [{"body", "face"}, …]} — removes the seeds on these edges and faces (by the names they were seeded with), or, with neither, every seed"#,
    form: Some(super::forms::SEED_CLEAR),
    kind: CommandKind::Model,
    requires: &[Requirement::Seed],
    run: seed_clear,
};

pub const SEED_LIST: CommandSpec = CommandSpec {
    id: "mesh.seed.list",
    label: "Seeds",
    params: "{} — the seeds on edges and faces, each with the edges of the open design it seeds, their lengths and the count it gives each (a face seed's edges with `overridden_by` where another seed decides them; `stale`: the design has no such edge or face)",
    form: Some(&[]),
    kind: CommandKind::Model,
    requires: &[],
    run: seed_list,
};

/// The `edges` or `faces` parameter (`key`, `kind` plural): a list of references to entities of
/// that kind.
fn refs(p: &Value, kind: &str, required: bool) -> Result<Option<Vec<GeometryRef>>> {
    let key = format!("{kind}s");
    let key = key.as_str();
    let refs: Option<Vec<GeometryRef>> = opt(p, key)?;
    let Some(refs) = refs else {
        return if required { Err(EngineError::Params(format!("`{key}` is required: a list of {{body, {kind}}}"))) } else { Ok(None) };
    };
    if refs.is_empty() || refs.len() > MAX_REFS {
        return Err(EngineError::Params(format!("`{key}` must list 1 to {MAX_REFS} {key}")));
    }
    if let Some(r) = refs.iter().find(|r| r.entity.kind() != kind) {
        return Err(EngineError::Params(format!("`{key}` lists the {r}; name {key} as {{body, {kind}}}")));
    }
    Ok(Some(refs))
}

/// The reference to face `index` of `body`, when it has a name.
fn face_ref(body: &Shape, index: u32) -> Option<GeometryRef> {
    GeometryRef::face(body.name.clone(), body.faces().get(index as usize)?.name.clone()).ok()
}

/// The faces `edge` of `body` bounds.
fn edge_faces(body: &Shape, edge: &Edge) -> Vec<GeometryRef> {
    edge.faces.iter().filter_map(|&f| face_ref(body, f)).collect()
}

fn edge_ref(body: &Shape, edge: &Edge) -> Option<GeometryRef> {
    GeometryRef::new(body.name.clone(), Entity::Edge(edge.name.clone())).ok()
}

/// The edges of the design a reference selects: an edge (its `~k` pieces too), or every edge
/// bounding a face (its `#k` pieces too).
pub(crate) fn seeded_edges(design: &Design, r: &GeometryRef) -> Result<Vec<Edge>> {
    let Some(body) = design.body(&r.body) else { return Ok(Vec::new()) };
    let on_face = matches!(r.entity, Entity::Face(_));
    Ok(body
        .edges()?
        .into_iter()
        .filter(|e| if on_face { edge_faces(body, e).iter().any(|f| r.selects(f)) } else { edge_ref(body, e).is_some_and(|g| r.selects(&g)) })
        .collect())
}

/// Whether the design has what a seed names: its edge, or its face.
pub(crate) fn resolves(design: &Design, r: &GeometryRef) -> Result<bool> {
    Ok(match (&r.entity, design.body(&r.body)) {
        (Entity::Face(_), Some(body)) => {
            body.faces().iter().filter_map(|f| GeometryRef::face(body.name.clone(), f.name.clone()).ok()).any(|f| r.selects(&f))
        }
        _ => !seeded_edges(design, r)?.is_empty(),
    })
}

/// Every seed on the open design's edges, as the meshers take them: each edge seed on its edges
/// (in the order of the seeds, which the meshers' results depend on), then on each edge without
/// one the face seed that decides it ([`fieldcraft_model::Model::edge_seed`]), less those an edge
/// seed or a finer face seed of another body holds on the same curve
/// ([`fieldcraft_mesh::seeds::with_face_seeds`]). A mesher matches them to its own curves by
/// position, so a seed on an edge two partition cells share holds in both.
pub(crate) fn edge_seeds(s: &Session) -> Result<Vec<EdgeSeed>> {
    let Some(design) = s.design.as_deref() else { return Ok(Vec::new()) };
    let seeds = s.model.seeds();
    // Each seeded body's edge tolerance and edges, found once.
    let mut bodies: BTreeMap<&str, (&Shape, f64, Vec<Edge>)> = BTreeMap::new();
    for body in design.bodies.iter().filter(|b| seeds.iter().any(|seed| seed.on.body == b.name)) {
        bodies.insert(body.name.as_str(), (body, (body.size()? * fieldcraft_geometry::EDGE_TOLERANCE).max(1e-4), body.edges()?));
    }
    let to_mesh = |seed: &Seed, e: &Edge, tolerance: f64, name: String| EdgeSeed {
        name,
        count: seed.count(e.length()) as usize,
        points: e.points.clone(),
        tolerance,
        bias: seed.bias,
        constraint: seed.constraint,
    };
    let mut on_edges = Vec::new();
    for seed in seeds.iter().filter(|seed| !seed.on_face()) {
        let Some((body, tolerance, edges)) = bodies.get(seed.on.body.as_str()) else { continue };
        for e in edges.iter().filter(|e| edge_ref(body, e).is_some_and(|g| seed.on.selects(&g))) {
            on_edges.push(to_mesh(seed, e, *tolerance, format!("{}:{}", body.name, e.name)));
        }
    }
    let mut on_faces = Vec::new();
    if seeds.iter().any(Seed::on_face) {
        for (body, tolerance, edges) in bodies.values() {
            for e in edges {
                let Some(edge) = edge_ref(body, e) else { continue };
                if let Some(face) = s.model.edge_seed(&edge, &edge_faces(body, e), e.length()).filter(|seed| seed.on_face()) {
                    let name = format!("{}:{} from the seed on face {}", body.name, e.name, face.on.entity.name().unwrap_or_default());
                    on_faces.push(to_mesh(face, e, *tolerance, name));
                }
            }
        }
    }
    Ok(fieldcraft_mesh::seeds::with_face_seeds(on_edges, on_faces))
}

fn design(s: &Session) -> Result<&Design> {
    s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first (geometry.open)".into()))
}

/// A seed's report: its definition and, per edge it seeds in the open design, the length and
/// the count it gives; on a face seed's edge another seed decides, `overridden_by` that seed
/// and its count.
fn report(s: &Session, seed: &Seed) -> Result<Value> {
    let mut v = json!(seed);
    let mut on = Vec::new();
    if let Some(design) = s.design.as_deref()
        && let Some(body) = design.body(&seed.on.body)
    {
        for e in seeded_edges(design, &seed.on)? {
            let length = e.length();
            let mut item = json!({"edge": e.name, "length": length, "count": seed.count(length)});
            if seed.on_face()
                && let Some(edge) = edge_ref(body, &e)
                && let Some(by) = s.model.edge_seed(&edge, &edge_faces(body, &e), length).filter(|by| by.on != seed.on)
            {
                item["overridden_by"] = json!({"seed": by.on.to_string(), "count": by.count(length)});
            }
            on.push(item);
        }
    }
    // An edge seed resolves when it has edges; a face seed when the face is there.
    v["stale"] = json!(match s.design.as_deref() {
        Some(design) if seed.on_face() => !resolves(design, &seed.on)?,
        Some(_) => on.is_empty(),
        None => true,
    });
    v["on"] = Value::Array(on);
    Ok(v)
}

/// Set seeds on the edges or faces (`kind`) `p` names, from its count or size, bias and
/// constraint.
fn seed_on(s: &mut Session, p: &Value, kind: &str) -> Result<Value> {
    let refs = refs(p, kind, true)?.unwrap_or_default();
    let divisions = match (opt::<u32>(p, "count")?, opt::<f64>(p, "size")?) {
        (Some(n), None) => Divisions::Count(n),
        (None, Some(h)) => Divisions::Size(h),
        _ => return Err(EngineError::Params("give one of `count` (divisions) or `size` (element length)".into())),
    };
    let bias: Option<Bias> = opt(p, "bias")?;
    let constraint: SeedConstraint = opt(p, "constraint")?.unwrap_or_default();
    let design = design(s)?.clone();
    let mut seeds = Vec::with_capacity(refs.len());
    for r in refs {
        if seeded_edges(&design, &r)?.is_empty() {
            return Err(EngineError::Params(if resolves(&design, &r)? {
                format!("the {r} has no edges to seed")
            } else {
                format!("the design has no {r} (geometry.bodies lists the {kind}s)")
            }));
        }
        seeds.push(Seed::new(r, divisions, bias, constraint)?);
    }
    for seed in &seeds {
        s.model.set_seed(seed.clone())?;
    }
    let out = seeds.iter().map(|seed| report(s, seed)).collect::<Result<Vec<_>>>()?;
    Ok(json!({"seeded": out.len(), "seeds": out}))
}

fn seed_edge(s: &mut Session, p: &Value) -> Result<Value> {
    seed_on(s, p, "edge")
}

fn seed_face(s: &mut Session, p: &Value) -> Result<Value> {
    seed_on(s, p, "face")
}

fn seed_clear(s: &mut Session, p: &Value) -> Result<Value> {
    let (edges, faces) = (refs(p, "edge", false)?, refs(p, "face", false)?);
    let cleared = if edges.is_none() && faces.is_none() {
        s.model.clear_seeds()
    } else {
        let named: Vec<GeometryRef> = edges.into_iter().chain(faces).flatten().collect();
        if let Some(r) = named.iter().find(|r| s.model.seed(r).is_none()) {
            return Err(EngineError::Params(format!("the {r} has no seed")));
        }
        named.iter().filter(|r| s.model.remove_seed(r)).count()
    };
    Ok(json!({"cleared": cleared, "left": s.model.seeds().len()}))
}

fn seed_list(s: &mut Session, _: &Value) -> Result<Value> {
    s.model.seeds().iter().map(|seed| report(s, seed)).collect::<Result<Vec<_>>>().map(Value::Array)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plate() -> Session {
        super::super::datum::tests::plate().0
    }

    fn thickness_edges(s: &mut Session) -> Vec<Value> {
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        bodies[0]["edges"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|e| (e["length"].as_f64().unwrap() - 10.0).abs() < 1e-6)
            .map(|e| json!({"body": "Body1", "edge": e["name"]}))
            .collect()
    }

    #[test]
    fn seeds_are_set_replaced_listed_cleared_and_undone() {
        let mut s = plate();
        let edges = thickness_edges(&mut s);
        let r = s.execute("mesh.seed.edge", &json!({"edges": edges[..2], "count": 4, "bias": {"ratio": 3, "toward": "end"}})).unwrap();
        assert_eq!(r["seeded"], 2);
        assert_eq!(r["seeds"][0]["on"][0]["count"], 4);
        assert_eq!(r["seeds"][0]["constraint"], "exact");
        // A size gives the length over it, rounded up; it replaces the count.
        let r = s.execute("mesh.seed.edge", &json!({"edges": [edges[0]], "size": 3.0, "constraint": "min"})).unwrap();
        assert_eq!(r["seeds"][0]["on"][0]["count"], 4);
        assert_eq!(s.model.seeds().len(), 2);
        let list = s.execute("mesh.seed.list", &json!({})).unwrap();
        assert_eq!(list.as_array().unwrap().len(), 2);
        assert!(list.as_array().unwrap().iter().all(|l| l["stale"] == false));

        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.model.seed(&serde_json::from_value(edges[0].clone()).unwrap()).unwrap().divisions, Divisions::Count(4));
        s.execute("mesh.seed.clear", &json!({"edges": [edges[1]]})).unwrap();
        assert_eq!(s.model.seeds().len(), 1);
        assert!(s.execute("mesh.seed.clear", &json!({"edges": [edges[1]]})).is_err());
        assert_eq!(s.execute("mesh.seed.clear", &json!({})).unwrap()["cleared"], 1);
        assert!(s.execute("mesh.seed.clear", &json!({})).is_err(), "nothing to clear");
    }

    #[test]
    fn seeds_survive_save_open_and_a_design_edit() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let edges = thickness_edges(&mut s);
        s.execute("mesh.seed.edge", &json!({"edges": edges, "count": 5, "bias": {"ratio": 2, "toward": "both"}})).unwrap();
        let file = dir.join("seeded.fieldcraft");
        s.execute("file.save", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut back = Session::default();
        back.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(back.model.seeds(), s.model.seeds());
        // The plate made thicker: the same edges, now 15 mm, keep their seeds.
        let design = dir.join("plate.solidcraft");
        let thick = fieldcraft_geometry::sample::Plate { thickness: 15.0, ..Default::default() };
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&thick).unwrap()).unwrap();
        let r = back.execute("geometry.open", &json!({"path": design.to_str().unwrap()})).unwrap();
        assert_eq!(r["stale"]["seeds"], json!([]));
        let list = back.execute("mesh.seed.list", &json!({})).unwrap();
        for l in list.as_array().unwrap() {
            assert_eq!(l["on"].as_array().unwrap().len(), 1, "{l}");
            assert!((l["on"][0]["length"].as_f64().unwrap() - 15.0).abs() < 1e-6, "{l}");
        }
    }

    /// The polyline of edge `name` of body `body`.
    fn edge_points(s: &Session, body: &str, name: &str) -> Vec<[f64; 3]> {
        s.design.as_ref().unwrap().body(body).unwrap().edges().unwrap().into_iter().find(|e| e.name == name).unwrap().points
    }

    /// The model's nodes on a straight edge, by their distance from its start, sorted.
    fn nodes_on(s: &Session, line: &[[f64; 3]]) -> Vec<f64> {
        let (a, b) = (line[0], *line.last().unwrap());
        let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let mut out: Vec<f64> = s
            .model
            .nodes()
            .values()
            .filter_map(|n| {
                // In millimetres, as the design is.
                let xyz = crate::units::to_geometry_point(n.xyz);
                let v = [xyz[0] - a[0], xyz[1] - a[1], xyz[2] - a[2]];
                let t = (v[0] * d[0] + v[1] * d[1] + v[2] * d[2]) / len;
                let off = [v[0] - t * d[0] / len, v[1] - t * d[1] / len, v[2] - t * d[2] / len];
                let off = (off[0] * off[0] + off[1] * off[1] + off[2] * off[2]).sqrt();
                (off < 1e-6 && (-1e-6..=len + 1e-6).contains(&t)).then_some(t)
            })
            .collect();
        out.sort_by(f64::total_cmp);
        out
    }

    #[test]
    fn the_meshers_honour_seeds() {
        let mut s = plate();
        let edges = thickness_edges(&mut s);
        let name = edges[0]["edge"].as_str().unwrap().to_string();
        let line = edge_points(&s, "Body1", &name);
        let seed = json!({"edges": [edges[0]], "count": 5, "bias": {"ratio": 3, "toward": "end"}});
        s.execute("mesh.seed.edge", &seed).unwrap();
        // Hex sweep: five layers through the thickness, graded geometrically toward the edge's end.
        s.execute("mesh.hex.sweep", &json!({"body": "Body1", "size": 4})).unwrap();
        let t = nodes_on(&s, &line);
        let dt: Vec<f64> = t.windows(2).map(|w| w[1] - w[0]).collect();
        assert_eq!(dt.len(), 5, "{t:?}");
        assert!((dt[0] / dt[4] - 3.0).abs() < 1e-6, "{dt:?}");
        // Tet: the edge is divided at the seed's nodes.
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 5})).unwrap();
        assert_eq!(r["bodies"][0]["seeds"][0]["mesh_edges"], 5, "{}", r["bodies"][0]["seeds"]);
        // Two exact seeds that cannot both hold are refused, naming both.
        s.execute("mesh.seed.edge", &json!({"edges": [edges[1]], "count": 3})).unwrap();
        let e = s.execute("mesh.hex.sweep", &json!({"body": "Body1", "size": 4})).unwrap_err().to_string();
        assert!(e.contains(&name) && e.contains(edges[1]["edge"].as_str().unwrap()) && e.contains("conflict"), "{e}");
    }

    /// The plate's edges whose polylines lie in the plane `x = 0` (its end face), by name.
    fn end_face_edges(s: &Session) -> Vec<(String, Vec<[f64; 3]>)> {
        let mut out: Vec<(String, Vec<[f64; 3]>)> = s
            .design
            .as_ref()
            .unwrap()
            .body("Body1")
            .unwrap()
            .edges()
            .unwrap()
            .into_iter()
            .filter(|e| e.points.iter().all(|p| p[0].abs() < 1e-9))
            .map(|e| (e.name, e.points))
            .collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }

    /// Mesh the plate with tetrahedra; its nodes on each of `edges`, as distances from its start.
    fn tet_mesh_nodes_on(s: &mut Session, edges: &[(String, Vec<[f64; 3]>)]) -> Vec<Vec<f64>> {
        s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 5})).unwrap();
        edges.iter().map(|(_, line)| nodes_on(s, line)).collect()
    }

    /// The plate's thickness edge on its x = 40 end, 10 mm.
    const SHORT_EDGE: &str = "F2:side:l2|F2:side:l3";

    /// A 30 mm edge of the plate's x = 0 end.
    fn long_end_edge(s: &Session) -> String {
        end_face_edges(s).into_iter().find(|(_, l)| (l[0][1] - l.last().unwrap()[1]).abs() > 29.0).unwrap().0
    }

    /// An exact seed of `count` alone on `edge`, the plate meshed with tetrahedra at size 5,
    /// gives `count` mesh edges, the seed's nodes at its places and no other node on the edge.
    fn assert_tet_seed_holds(edge: &str, count: u32) {
        let mut s = plate();
        let e = s.design.as_ref().unwrap().body("Body1").unwrap().edges().unwrap().into_iter().find(|e| e.name == edge).unwrap();
        s.execute("mesh.seed.edge", &json!({"edges": [{"body": "Body1", "edge": edge}], "count": count})).unwrap();
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 5})).unwrap();
        let t = nodes_on(&s, &e.points);
        assert_eq!(t.len(), count as usize + 1, "{edge} seeded {count}: {t:?}");
        assert_eq!(r["bodies"][0]["seeds"][0]["mesh_edges"], count, "{edge} seeded {count}");
        for (k, x) in t.iter().enumerate() {
            assert!((x - e.length() * k as f64 / f64::from(count)).abs() < 1e-6, "{edge} seeded {count}: {t:?}");
        }
    }

    /// An exact seed's nodes are the only mesh nodes on its edge, at the seed's places (#894:
    /// a seed of 12 on a 30 mm edge gave 11).
    #[test]
    fn an_exact_tet_seed_holds_its_count() {
        assert_tet_seed_holds(&long_end_edge(&plate()), 12);
    }

    /// A fine seed on a short edge: 8 on the 10 mm thickness edge, spacing 1.25 against size 5,
    /// gave 7 (#901: the node at 2.5 mm was refused by a flat cavity).
    #[test]
    fn a_fine_exact_tet_seed_on_a_short_edge_holds_its_count() {
        assert_tet_seed_holds(SHORT_EDGE, 8);
    }

    /// Exact seeds coarser than the size and much finer than it hold on the short edge (20 gave
    /// 19).
    #[test]
    fn exact_tet_seeds_of_any_count_hold_on_a_short_edge() {
        for count in [1, 3, 20] {
            assert_tet_seed_holds(SHORT_EDGE, count);
        }
    }

    /// And on a long edge, coarser and finer than the size.
    #[test]
    fn exact_tet_seeds_of_any_count_hold_on_a_long_edge() {
        let long = long_end_edge(&plate());
        for count in [5, 30] {
            assert_tet_seed_holds(&long, count);
        }
    }

    /// Four edges of one face, each seeded 6, give 6 each, and the mesh does not depend on the
    /// order the seeds were given in (#894: 8, 6, 5 and 9, changing with the order).
    #[test]
    fn exact_tet_seeds_on_a_face_hold_in_any_order() {
        let mut s = plate();
        let edges = end_face_edges(&s);
        assert_eq!(edges.len(), 4, "{edges:?}");
        let refs: Vec<Value> = edges.iter().map(|(n, _)| json!({"body": "Body1", "edge": n})).collect();
        s.execute("mesh.seed.edge", &json!({"edges": refs, "count": 6})).unwrap();
        let forward = tet_mesh_nodes_on(&mut s, &edges);
        for t in &forward {
            assert_eq!(t.len(), 7, "{forward:?}");
        }
        let mut back = plate();
        for r in refs.iter().rev() {
            back.execute("mesh.seed.edge", &json!({"edges": [r], "count": 6})).unwrap();
        }
        assert_eq!(tet_mesh_nodes_on(&mut back, &edges), forward);
        assert!(back.model.nodes() == s.model.nodes() && back.model.elements() == s.model.elements(), "the same mesh");
    }

    #[test]
    fn a_seed_on_a_cut_edge_holds_in_both_cells() {
        let (mut s, _dir) = super::super::datum::tests::plate();
        s.execute("geometry.datum", &json!({"name": "X30", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
        s.execute("geometry.partition", &json!({"name": "P", "datum": "X30"})).unwrap();
        // An edge of the cut face, 30 mm along y, on the cell with the hole.
        let edge = s
            .design
            .as_ref()
            .unwrap()
            .body("Body1/a")
            .unwrap()
            .edges()
            .unwrap()
            .into_iter()
            .find(|e| e.name.contains("partition:P:a") && (e.length() - 30.0).abs() < 1e-6)
            .unwrap();
        s.execute("mesh.seed.edge", &json!({"edges": [{"body": "Body1/a", "edge": edge.name}], "count": 9})).unwrap();
        let r = s.execute("mesh.cells", &json!({"default": "hex", "size": 2.5})).unwrap();
        assert_eq!(r["interfaces"][0]["joined"], "merged", "{r}");
        assert_eq!(nodes_on(&s, &edge.points).len(), 10, "{r}");
    }

    #[test]
    fn hostile_seeds_are_refused_and_change_nothing() {
        let mut s = plate();
        let edges = thickness_edges(&mut s);
        let e = edges[0].clone();
        for p in [
            json!({}),
            json!({"edges": [], "count": 2}),
            json!({"edges": [e], "count": 2, "size": 1.0}),
            json!({"edges": [e]}),
            json!({"edges": [e], "count": 0}),
            json!({"edges": [e], "count": -3}),
            json!({"edges": [e], "count": 1e12}),
            json!({"edges": [e], "size": 0.0}),
            json!({"edges": [e], "count": 2, "bias": {"ratio": 1e9, "toward": "end"}}),
            json!({"edges": [e], "count": 2, "bias": {"ratio": 2, "toward": "middle"}}),
            json!({"edges": [e], "count": 2, "constraint": "loose"}),
            json!({"edges": [{"body": "Body1", "face": "x"}], "count": 2}),
            json!({"edges": [{"body": "Body1", "edge": "no|such"}], "count": 2}),
            json!({"edges": [{"body": "Nope", "edge": "a|b"}], "count": 2}),
            json!({"edges": "all", "count": 2}),
        ] {
            assert!(s.execute("mesh.seed.edge", &p).is_err(), "{p}");
        }
        assert!(s.model.seeds().is_empty());
        assert!(s.journal.is_empty() || s.journal.iter().all(|(c, _)| c != "mesh.seed.edge"));
    }

    /// The plate's end face (x = 40): two thickness edges of 10 mm and two of 30 mm along y.
    fn end_face(s: &Session) -> Value {
        json!({"body": "Body1", "face": super::super::sizes::tests::face_at_x(s, 40.0)})
    }

    /// The (length, count) a seed report gives each of its edges, sorted.
    fn counts(report: &Value) -> Vec<(i64, i64)> {
        let mut c: Vec<(i64, i64)> =
            report["on"].as_array().unwrap().iter().map(|o| (o["length"].as_f64().unwrap().round() as i64, o["count"].as_i64().unwrap())).collect();
        c.sort();
        c
    }

    /// A thickness edge of the plate's end face, as `{body, edge}`.
    fn end_thickness_edge(s: &mut Session) -> Value {
        let design = s.design.clone().unwrap();
        thickness_edges(s)
            .into_iter()
            .find(|e| {
                let r: GeometryRef = serde_json::from_value(e.clone()).unwrap();
                seeded_edges(&design, &r).unwrap()[0].points.iter().all(|p| (p[0] - 40.0).abs() < 1e-6)
            })
            .unwrap()
    }

    #[test]
    fn face_seeds_are_set_listed_overridden_cleared_and_undone() {
        let mut s = plate();
        let end = end_face(&s);
        // A size gives each edge of the face its own count.
        let r = s.execute("mesh.seed.face", &json!({"faces": [end], "size": 2.5})).unwrap();
        assert_eq!(r["seeded"], 1);
        assert_eq!(r["seeds"][0]["face"], end, "{r}");
        assert_eq!(counts(&r["seeds"][0]), [(10, 4), (10, 4), (30, 12), (30, 12)], "{r}");
        // A count is each edge's count; it replaces the size.
        let r =
            s.execute("mesh.seed.face", &json!({"faces": [end], "count": 5, "bias": {"ratio": 2, "toward": "both"}, "constraint": "min"})).unwrap();
        assert_eq!(counts(&r["seeds"][0]), [(10, 5), (10, 5), (30, 5), (30, 5)], "{r}");
        assert_eq!(s.model.seeds().len(), 1);

        // An exact seed on one of its edges wins there, and the face seed's report names it.
        let edge = end_thickness_edge(&mut s);
        s.execute("mesh.seed.edge", &json!({"edges": [edge], "count": 3})).unwrap();
        let list = s.execute("mesh.seed.list", &json!({})).unwrap();
        let list = list.as_array().unwrap();
        assert_eq!(list.len(), 2, "{list:?}");
        let face = list.iter().find(|l| l.get("face").is_some()).unwrap();
        assert_eq!(face["stale"], false);
        let over: Vec<&Value> = face["on"].as_array().unwrap().iter().filter(|o| o.get("overridden_by").is_some()).collect();
        assert_eq!(over.len(), 1, "{face}");
        assert_eq!(over[0]["edge"], edge["edge"]);
        assert_eq!(over[0]["overridden_by"]["count"], 3);
        assert!(over[0]["overridden_by"]["seed"].as_str().unwrap().contains(edge["edge"].as_str().unwrap()), "{face}");

        // Clearing the face's seed leaves the edge's; undo brings it back.
        s.execute("mesh.seed.clear", &json!({"faces": [end]})).unwrap();
        assert_eq!(s.model.seeds().len(), 1);
        assert!(s.execute("mesh.seed.clear", &json!({"faces": [end]})).is_err(), "no seed left on the face");
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.model.seeds().len(), 2);
        assert_eq!(s.execute("mesh.seed.clear", &json!({})).unwrap()["cleared"], 2);
    }

    #[test]
    fn face_seeds_survive_save_open_and_a_design_edit() {
        let (mut s, dir) = super::super::datum::tests::plate();
        let end = end_face(&s);
        s.execute("mesh.seed.face", &json!({"faces": [end], "size": 2.5})).unwrap();
        let file = dir.join("face-seeded.fieldcraft");
        s.execute("file.save", &json!({"path": file.to_str().unwrap()})).unwrap();
        let mut back = Session::default();
        back.execute("file.open", &json!({"path": file.to_str().unwrap()})).unwrap();
        assert_eq!(back.model.seeds(), s.model.seeds());
        // The plate made thicker: the face keeps its name and its seed, now on 15 mm edges.
        let design = dir.join("plate.solidcraft");
        let thick = fieldcraft_geometry::sample::Plate { thickness: 15.0, ..Default::default() };
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&thick).unwrap()).unwrap();
        let r = back.execute("geometry.open", &json!({"path": design.to_str().unwrap()})).unwrap();
        assert_eq!(r["stale"]["seeds"], json!([]));
        let list = back.execute("mesh.seed.list", &json!({})).unwrap();
        assert_eq!(counts(&list[0]), [(15, 6), (15, 6), (30, 12), (30, 12)], "{list}");
    }

    /// A face seed sets the hex sweep's layers and the cross-section's divisions along the face's
    /// edges; an exact edge seed that disagrees with it across the sweep is refused, naming both.
    #[test]
    fn the_hex_mesher_honours_face_seeds() {
        let mut s = plate();
        let end = end_face(&s);
        s.execute("mesh.seed.face", &json!({"faces": [end], "size": 2.5})).unwrap();
        s.execute("mesh.hex.sweep", &json!({"body": "Body1", "size": 4})).unwrap();
        let on = |s: &Session, len: f64| -> Vec<usize> {
            let edges = s.design.as_ref().unwrap().body("Body1").unwrap().edges().unwrap();
            let mut n: Vec<usize> = edges
                .iter()
                .filter(|e| e.points.iter().all(|p| (p[0] - 40.0).abs() < 1e-6) && (e.length() - len).abs() < 1e-6)
                .map(|e| nodes_on(s, &e.points).len() - 1)
                .collect();
            n.sort();
            n
        };
        assert_eq!(on(&s, 10.0), [4, 4]);
        assert_eq!(on(&s, 30.0), [12, 12]);
        // Three layers on one thickness edge cannot hold beside the face seed's four.
        let edge = end_thickness_edge(&mut s);
        s.execute("mesh.seed.edge", &json!({"edges": [edge], "count": 3})).unwrap();
        let e = s.execute("mesh.hex.sweep", &json!({"body": "Body1", "size": 4})).unwrap_err().to_string();
        assert!(e.contains(edge["edge"].as_str().unwrap()) && e.contains("from the seed on face") && e.contains("conflict"), "{e}");
    }

    /// The tet mesher divides every edge of a seeded face at the seed's nodes, and an edge with its
    /// own seed at that seed's.
    #[test]
    fn the_tet_mesher_honours_face_seeds() {
        let tet = |s: &mut Session| -> Vec<(i64, i64)> {
            let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 5})).unwrap();
            let mut c: Vec<(i64, i64)> = r["bodies"][0]["seeds"]
                .as_array()
                .unwrap()
                .iter()
                .map(|e| (e["count"].as_i64().unwrap(), e["mesh_edges"].as_i64().unwrap()))
                .collect();
            c.sort();
            c
        };
        let mut s = plate();
        let end = end_face(&s);
        // An exact count of 6 holds on every edge of the face, the 10 mm ones and the 30 mm ones.
        s.execute("mesh.seed.face", &json!({"faces": [end], "count": 6})).unwrap();
        assert_eq!(tet(&mut s), [(6, 6); 4]);
        s.execute("mesh.seed.face", &json!({"faces": [end], "size": 5})).unwrap();
        assert_eq!(tet(&mut s), [(2, 2), (2, 2), (6, 6), (6, 6)]);
        // A finer seed on one thickness edge holds there in place of the face's.
        let edge = end_thickness_edge(&mut s);
        s.execute("mesh.seed.edge", &json!({"edges": [edge], "count": 8})).unwrap();
        let r = s.execute("mesh.tet.generate", &json!({"body": "Body1", "size": 5})).unwrap();
        let seeds = r["bodies"][0]["seeds"].as_array().unwrap();
        assert_eq!(seeds.len(), 4, "{seeds:?}");
        let own: Vec<&Value> = seeds.iter().filter(|e| !e["edge"].as_str().unwrap().contains("from the seed on face")).collect();
        assert_eq!(own.len(), 1, "{seeds:?}");
        assert!(own[0]["edge"].as_str().unwrap().ends_with(edge["edge"].as_str().unwrap()), "{seeds:?}");
        assert_eq!((own[0]["count"].clone(), own[0]["mesh_edges"].clone()), (json!(8), json!(8)), "{seeds:?}");
    }

    #[test]
    fn hostile_face_seeds_are_refused_and_change_nothing() {
        let mut s = plate();
        let f = end_face(&s);
        for p in [
            json!({}),
            json!({"faces": [], "count": 2}),
            json!({"faces": [f], "count": 2, "size": 1.0}),
            json!({"faces": [f]}),
            json!({"faces": [f], "count": 0}),
            json!({"faces": [f], "count": -3}),
            json!({"faces": [f], "count": 1e12}),
            json!({"faces": [f], "size": 0.0}),
            json!({"faces": [f], "size": -1.0}),
            json!({"faces": [f], "size": "big"}),
            json!({"faces": [f], "count": 2, "bias": {"ratio": 1e9, "toward": "both"}}),
            json!({"faces": [f], "count": 2, "bias": {"ratio": 2, "toward": "end"}}),
            json!({"faces": [f], "count": 2, "bias": {"ratio": 2, "toward": "start"}}),
            json!({"faces": [f], "count": 2, "constraint": "loose"}),
            json!({"faces": [{"body": "Body1", "edge": "a|b"}], "count": 2}),
            json!({"faces": [{"body": "Body1", "face": "no such"}], "count": 2}),
            json!({"faces": [{"body": "Nope", "face": "F1"}], "count": 2}),
            json!({"faces": "all", "count": 2}),
            json!({"edges": [f], "count": 2}),
        ] {
            assert!(s.execute("mesh.seed.face", &p).is_err(), "{p}");
        }
        assert!(s.model.seeds().is_empty());
        for p in [json!({"faces": [f]}), json!({"faces": []}), json!({"faces": [{"body": "Body1", "edge": "a|b"}]})] {
            assert!(s.execute("mesh.seed.clear", &p).is_err(), "{p}");
        }
        assert!(s.journal.iter().all(|(c, _)| c != "mesh.seed.face"));
    }
}
