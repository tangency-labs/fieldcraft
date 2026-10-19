//! Edge seeds: `mesh.seed.edge` sets a count or a size (with bias and constraint) on design
//! edges, `mesh.seed.clear` removes seeds and `mesh.seed.list` lists them with the count each
//! gives on its edge. Seeds live in the model (`fieldcraft_model::seeds`), are saved with it, and
//! the meshers honour them.

use fieldcraft_geometry::{Design, Edge};
use fieldcraft_model::{Bias, Divisions, Entity, GeometryRef, Seed, SeedConstraint};
use serde_json::{Value, json};

use super::opt;
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

/// Most edges one command names.
const MAX_EDGES: usize = 10_000;

pub const SEED_EDGE: CommandSpec = CommandSpec {
    id: "mesh.seed.edge",
    label: "Seed Edges",
    params: r#"{"edges": [{"body": "Body1", "edge": "F2:side:hi|F2:top"}, …], "count": 8 | "size": 2.5, "bias"?: {"ratio": 4, "toward": "start" | "end" | "both"}, "constraint"?: "exact" (default) | "min" | "free"} — seeds design edges (persistent edge names, geometry.bodies): `count` divisions, or `size` (the edge's length over it, rounded up). `bias` grades the elements geometrically, the largest `ratio` times the smallest, small at the edge's start, end or both ends. `constraint`: "exact" keeps the count (a mesher whose other edges cannot match refuses, naming the conflicting edges), "min" keeps at least it, "free" makes it a goal. Replaces a seed the edge had; an edge name without its `~k` number seeds every edge of that name. The hex meshers (sweep, mapped, blocks, cells) use the counts and spacing; the tet mesher uses the spacing on the edge and grades its size field from it"#,
    form: Some(super::forms::SEED_EDGE),
    kind: CommandKind::Model,
    requires: &[Requirement::Design],
    run: seed_edge,
};

pub const SEED_CLEAR: CommandSpec = CommandSpec {
    id: "mesh.seed.clear",
    label: "Clear Seeds",
    params: r#"{"edges"?: [{"body", "edge"}, …]} — removes the seeds on these edges (by the names they were seeded with), or every seed"#,
    form: Some(super::forms::SEED_CLEAR),
    kind: CommandKind::Model,
    requires: &[Requirement::Seed],
    run: seed_clear,
};

pub const SEED_LIST: CommandSpec = CommandSpec {
    id: "mesh.seed.list",
    label: "Seeds",
    params: "{} — the edge seeds, each with the edges of the open design it seeds, their lengths and the count it gives each (`stale`: the design has no such edge)",
    form: Some(&[]),
    kind: CommandKind::Model,
    requires: &[],
    run: seed_list,
};

/// The `edges` parameter: a list of edge references.
fn edge_refs(p: &Value, required: bool) -> Result<Option<Vec<GeometryRef>>> {
    let refs: Option<Vec<GeometryRef>> = opt(p, "edges")?;
    let Some(refs) = refs else {
        return if required { Err(EngineError::Params("`edges` is required: a list of {body, edge}".into())) } else { Ok(None) };
    };
    if refs.is_empty() || refs.len() > MAX_EDGES {
        return Err(EngineError::Params(format!("`edges` must list 1 to {MAX_EDGES} edges")));
    }
    if let Some(r) = refs.iter().find(|r| !matches!(r.entity, Entity::Edge(_))) {
        return Err(EngineError::Params(format!("`edges` lists the {r}; name edges as {{body, edge}}")));
    }
    Ok(Some(refs))
}

/// The edges of the design a reference selects (its `~k` pieces too).
pub(crate) fn seeded_edges(design: &Design, r: &GeometryRef) -> Result<Vec<Edge>> {
    let Some(body) = design.body(&r.body) else { return Ok(Vec::new()) };
    Ok(body.edges()?.into_iter().filter(|e| GeometryRef::new(r.body.clone(), Entity::Edge(e.name.clone())).is_ok_and(|g| r.selects(&g))).collect())
}

/// Every seed on the open design's edges, as the meshers take them. A mesher matches them to its
/// own curves by position, so a seed on an edge two partition cells share holds in both.
pub(crate) fn edge_seeds(s: &Session) -> Result<Vec<fieldcraft_mesh::seeds::EdgeSeed>> {
    let Some(design) = s.design.as_deref() else { return Ok(Vec::new()) };
    let mut out = Vec::new();
    // Each body's edges, found once.
    let mut edges: std::collections::BTreeMap<&str, (f64, Vec<Edge>)> = std::collections::BTreeMap::new();
    for seed in s.model.seeds() {
        let Some(body) = design.body(&seed.edge.body) else { continue };
        let (tolerance, all) = match edges.entry(body.name.as_str()) {
            std::collections::btree_map::Entry::Occupied(e) => e.into_mut(),
            std::collections::btree_map::Entry::Vacant(v) => {
                v.insert(((body.size()? * fieldcraft_geometry::EDGE_TOLERANCE).max(1e-4), body.edges()?))
            }
        };
        let selected = all.iter().filter(|e| GeometryRef::new(body.name.clone(), Entity::Edge(e.name.clone())).is_ok_and(|g| seed.edge.selects(&g)));
        for e in selected.cloned().collect::<Vec<_>>() {
            let tolerance = *tolerance;
            out.push(fieldcraft_mesh::seeds::EdgeSeed {
                name: format!("{}:{}", seed.edge.body, e.name),
                count: seed.count(e.length()) as usize,
                points: e.points,
                tolerance,
                bias: seed.bias,
                constraint: seed.constraint,
            });
        }
    }
    Ok(out)
}

fn design(s: &Session) -> Result<&Design> {
    s.design.as_deref().ok_or_else(|| EngineError::Params("open a design first (geometry.open)".into()))
}

/// A seed's report: its definition and, per edge it seeds, the length and the count it gives.
fn report(seed: &Seed, edges: &[Edge]) -> Value {
    let mut v = json!(seed);
    v["on"] = edges.iter().map(|e| json!({"edge": e.name, "length": e.length(), "count": seed.count(e.length())})).collect();
    v["stale"] = json!(edges.is_empty());
    v
}

fn seed_edge(s: &mut Session, p: &Value) -> Result<Value> {
    let refs = edge_refs(p, true)?.unwrap_or_default();
    let divisions = match (opt::<u32>(p, "count")?, opt::<f64>(p, "size")?) {
        (Some(n), None) => Divisions::Count(n),
        (None, Some(h)) => Divisions::Size(h),
        _ => return Err(EngineError::Params("give one of `count` (divisions) or `size` (element length)".into())),
    };
    let bias: Option<Bias> = opt(p, "bias")?;
    let constraint: SeedConstraint = opt(p, "constraint")?.unwrap_or_default();
    let design = design(s)?.clone();
    let mut out = Vec::with_capacity(refs.len());
    for r in refs {
        let edges = seeded_edges(&design, &r)?;
        if edges.is_empty() {
            return Err(EngineError::Params(format!("the design has no {r} (geometry.bodies lists the edges)")));
        }
        let seed = Seed::new(r, divisions, bias, constraint)?;
        out.push(report(&seed, &edges));
        s.model.set_seed(seed)?;
    }
    Ok(json!({"seeded": out.len(), "seeds": out}))
}

fn seed_clear(s: &mut Session, p: &Value) -> Result<Value> {
    let cleared = match edge_refs(p, false)? {
        None => s.model.clear_seeds(),
        Some(refs) => {
            if let Some(r) = refs.iter().find(|r| s.model.seed(r).is_none()) {
                return Err(EngineError::Params(format!("the {r} has no seed")));
            }
            refs.iter().filter(|r| s.model.remove_seed(r)).count()
        }
    };
    Ok(json!({"cleared": cleared, "left": s.model.seeds().len()}))
}

fn seed_list(s: &mut Session, _: &Value) -> Result<Value> {
    let mut out = Vec::with_capacity(s.model.seeds().len());
    for seed in s.model.seeds() {
        let edges = match s.design.as_deref() {
            Some(d) => seeded_edges(d, &seed.edge)?,
            None => Vec::new(),
        };
        out.push(report(seed, &edges));
    }
    Ok(Value::Array(out))
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
}
