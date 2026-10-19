//! Node tools (#234): `mesh.node.merge` equivalences coincident nodes, `mesh.check.duplicates`
//! finds what is not merged or does not join (coincident nodes, repeated elements, free faces and
//! edges, cracks), and `mesh.renumber` gives nodes and elements new labels. Merging and
//! renumbering move every reference with the labels (connectivity, sets, geometry tags, a step's
//! displacement limit), so boundary conditions and loads, which act on sets, stay on the same
//! entities; deck writers write the labels as they are. Both are one undo step; results solved
//! before go stale, as after any mesh edit.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_mesh::nodes::{Keep, bandwidth, duplicates, merge_plan, reverse_cuthill_mckee};
use fieldcraft_model::geometry::Entity;
use fieldcraft_model::{MAX_LABEL, Model};
use serde_json::{Value, json};

use crate::requirement::Requirement;
use crate::schema::ProblemKind;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const MERGE: CommandSpec = CommandSpec {
    id: "mesh.node.merge",
    label: "Merge Nodes",
    params: r#"{"tolerance"?: length (default a millionth of the mesh's bounding-box diagonal), "set"?: "Name" (only the set's nodes and those of its elements and faces), "keep"?: "lowest" (default) | "highest"} — equivalences coincident nodes: every node within the distance of another becomes that one (the lowest or highest label of each group survives, keeping its position) in every element, set and geometry tag, and is removed. A merge that would put two nodes of one element together (collapsing it) is left out and listed under `refused`; the rest are made. Reports `merged`, the `pairs` [node, survivor] (the first 1000) and `refused`"#,
    form: Some(super::forms::NODE_MERGE),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: merge,
};

pub const DUPLICATES: CommandSpec = CommandSpec {
    id: "mesh.check.duplicates",
    label: "Check Duplicates",
    params: r#"{"tolerance"?: length (default a millionth of the mesh's bounding-box diagonal)} — checks the mesh for what was not merged or does not join, changing nothing: `coincident` nodes within the distance of another ([node, other, distance]), `duplicate_elements` (groups on the same nodes), `free_faces` (faces of solid elements no other element shares: the skin), `cracks` (free faces lying on another free face: parts that touch without sharing nodes) and `free_edges` (edges of the skin and of 2D elements on one facet only: open borders). Lists hold the first 1000; `counts` has the totals"#,
    form: Some(super::forms::CHECK_DUPLICATES),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: check_duplicates,
};

pub const RENUMBER: CommandSpec = CommandSpec {
    id: "mesh.renumber",
    label: "Renumber",
    params: r#"{"what"?: "both" (default) | "nodes" | "elements", "by"?: "compact" (default: in label order) | "bandwidth" (nodes in reverse Cuthill–McKee order, elements by their lowest new node label) | "body" (body by body, by name; within a body in label order), "start"?: 1, "offset"?: n (instead of `by` and `start`: add n to every label)} — gives nodes and elements new labels from `start` up (or shifted by `offset`), every set, geometry tag and step following, so boundary conditions and loads stay on the same entities and decks write the new labels. Refused, changing nothing, when a label would leave 1..=99999999 or a node would take a reference point's label. Reports the counts renumbered, the label ranges and the bandwidth before and after"#,
    form: Some(super::forms::RENUMBER),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: renumber,
};

/// Most items a reply lists.
const MAX_LISTED: usize = 1000;

/// The merge distance: `tolerance`, or a millionth of the mesh's bounding-box diagonal.
pub(crate) fn tolerance(model: &Model, p: &Value) -> Result<f64> {
    match p.get("tolerance") {
        Some(v) => v.as_f64().filter(|t| *t > 0.0 && t.is_finite()).ok_or_else(|| EngineError::param("tolerance", ProblemKind::NotPositive)),
        None => {
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            for n in model.nodes().values() {
                for k in 0..3 {
                    lo[k] = lo[k].min(n.xyz[k]);
                    hi[k] = hi[k].max(n.xyz[k]);
                }
            }
            let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
            Ok(if diagonal > 0.0 && diagonal.is_finite() { 1e-6 * diagonal } else { 1e-9 })
        }
    }
}

/// The merge's candidates: the nodes of `set` (its own, and those of its elements and faces).
fn candidates(model: &Model, p: &Value) -> Result<Option<BTreeSet<u32>>> {
    let Some(name) = super::opt::<String>(p, "set")? else { return Ok(None) };
    let set = model.resolve_set(&name)?;
    let mut nodes: BTreeSet<u32> = set.nodes.into_iter().collect();
    for e in set.elements.iter().chain(set.faces.iter().map(|f| &f.element)) {
        if let Some(el) = model.element(*e) {
            nodes.extend(el.nodes.iter().copied());
        }
    }
    Ok(Some(nodes))
}

fn keep(p: &Value) -> Result<Keep> {
    match super::opt::<String>(p, "keep")?.as_deref() {
        None | Some("lowest") => Ok(Keep::Lowest),
        Some("highest") => Ok(Keep::Highest),
        Some(_) => Err(EngineError::param("keep", ProblemKind::NotAChoice)),
    }
}

/// The merges `mesh.node.merge` with `params` would make, as [node, survivor] pairs: what its
/// form previews. `None` when the parameters are not (yet) valid.
pub fn merge_preview(s: &Session, params: &Value) -> Option<Vec<[u32; 2]>> {
    // The meshers work in millimetres.
    let framed = crate::units::geometry_session(s).ok()?;
    let params = crate::units::geometry_params(s, MERGE.id, params).ok()?;
    let (s, params) = (&framed, &params);
    let tol = tolerance(&s.model, params).ok()?;
    let only = candidates(&s.model, params).ok()?;
    let plan = merge_plan(&s.model, only.as_ref(), tol, keep(params).ok()?).ok()?;
    Some(plan.map.into_iter().map(|(a, b)| [a, b]).collect())
}

fn merge(s: &mut Session, p: &Value) -> Result<Value> {
    let tol = tolerance(&s.model, p)?;
    let only = candidates(&s.model, p)?;
    let plan = merge_plan(&s.model, only.as_ref(), tol, keep(p)?)?;
    if !plan.map.is_empty() {
        s.model.merge_nodes(&plan.map)?;
    }
    let pairs: Vec<[u32; 2]> = plan.map.iter().take(MAX_LISTED).map(|(&a, &b)| [a, b]).collect();
    let refused: Vec<Value> =
        plan.refused.iter().take(MAX_LISTED).map(|r| json!({"node": r.node, "survivor": r.survivor, "element": r.element})).collect();
    let mut summary =
        format!("merged {} node{} within {}", plan.map.len(), if plan.map.len() == 1 { "" } else { "s" }, crate::units::shown_length(s, tol));
    if !plan.refused.is_empty() {
        summary.push_str(&format!("; {} left apart (they would collapse an element)", plan.refused.len()));
    }
    Ok(json!({
        "merged": plan.map.len(),
        "pairs": pairs,
        "refused": refused,
        "tolerance": tol,
        "nodes": s.model.nodes().len(),
        "summary": summary,
    }))
}

fn check_duplicates(s: &mut Session, p: &Value) -> Result<Value> {
    let tol = tolerance(&s.model, p)?;
    let d = duplicates(&s.model, tol)?;
    let coincident: Vec<Value> = d.coincident.iter().take(MAX_LISTED).map(|&(a, b, dist)| json!([a, b, dist])).collect();
    let summary = format!(
        "{} coincident node{}, {} duplicate element group{}, {} crack face{}, {} free edge{} ({} free faces)",
        d.coincident.len(),
        plural(d.coincident.len()),
        d.elements.len(),
        plural(d.elements.len()),
        d.cracks.len(),
        plural(d.cracks.len()),
        d.free_edges.len(),
        plural(d.free_edges.len()),
        d.free_faces,
    );
    Ok(json!({
        "tolerance": tol,
        "coincident": coincident,
        "duplicate_elements": d.elements.iter().take(MAX_LISTED).collect::<Vec<_>>(),
        "free_faces": d.free_faces,
        "cracks": d.cracks.iter().take(MAX_LISTED).collect::<Vec<_>>(),
        "free_edges": d.free_edges.iter().take(MAX_LISTED).collect::<Vec<_>>(),
        "counts": {
            "coincident": d.coincident.len(),
            "duplicate_elements": d.elements.len(),
            "cracks": d.cracks.len(),
            "free_edges": d.free_edges.len(),
        },
        "summary": summary,
    }))
}

fn plural(n: usize) -> &'static str {
    if n == 1 { "" } else { "s" }
}

/// The body name of each element tagged with one.
fn element_bodies(model: &Model) -> BTreeMap<u32, String> {
    let c = model.classification();
    c.elements()
        .iter()
        .filter_map(|(&e, t)| {
            let r = c.entity(t.region?)?;
            (r.entity == Entity::Body).then(|| (e, r.body.clone()))
        })
        .collect()
}

/// New labels for `order` from `start` up; `None` when they would leave the label range.
fn from_start(order: &[u32], start: u32) -> Option<BTreeMap<u32, u32>> {
    order
        .iter()
        .enumerate()
        .map(|(i, &l)| Some((l, start.checked_add(u32::try_from(i).ok()?).filter(|x| *x <= MAX_LABEL)?)))
        .filter(|x| x.as_ref().is_none_or(|(a, b)| a != b))
        .collect()
}

fn renumber(s: &mut Session, p: &Value) -> Result<Value> {
    let (nodes_too, elements_too) = match super::opt::<String>(p, "what")?.as_deref() {
        None | Some("both") => (true, true),
        Some("nodes") => (true, false),
        Some("elements") => (false, true),
        Some(_) => return Err(EngineError::param("what", ProblemKind::NotAChoice)),
    };
    let model = &s.model;
    let out_of_range = || EngineError::Params(format!("the new labels would leave 1..={MAX_LABEL}"));
    let (node_map, element_map) = if let Some(v) = p.get("offset") {
        if p.get("by").is_some() || p.get("start").is_some() {
            return Err(EngineError::Params("give `offset`, or `by` and `start`, not both".into()));
        }
        let offset = v
            .as_i64()
            .ok_or_else(|| EngineError::param("offset", ProblemKind::NotAWholeNumber { min: -i64::from(MAX_LABEL), max: i64::from(MAX_LABEL) }))?;
        let shift = |labels: &mut dyn Iterator<Item = u32>| -> Result<BTreeMap<u32, u32>> {
            labels
                .filter(|_| offset != 0)
                .map(|l| u32::try_from(i64::from(l) + offset).ok().filter(|x| (1..=MAX_LABEL).contains(x)).map(|n| (l, n)).ok_or_else(out_of_range))
                .collect()
        };
        let n = if nodes_too { shift(&mut model.nodes().keys().copied())? } else { BTreeMap::new() };
        let e = if elements_too { shift(&mut model.elements().keys().copied())? } else { BTreeMap::new() };
        (n, e)
    } else {
        let start = match p.get("start") {
            None => 1,
            Some(v) => v
                .as_u64()
                .and_then(|x| u32::try_from(x).ok())
                .filter(|x| (1..=MAX_LABEL).contains(x))
                .ok_or_else(|| EngineError::param("start", ProblemKind::NotAWholeNumber { min: 1, max: i64::from(MAX_LABEL) }))?,
        };
        let (node_order, element_order): (Vec<u32>, Vec<u32>) = match super::opt::<String>(p, "by")?.as_deref() {
            None | Some("compact") => (model.nodes().keys().copied().collect(), model.elements().keys().copied().collect()),
            Some("bandwidth") => {
                let nodes = reverse_cuthill_mckee(model);
                let rank: BTreeMap<u32, usize> = nodes.iter().enumerate().map(|(i, &l)| (l, i)).collect();
                let mut elements: Vec<(usize, u32)> = model
                    .elements()
                    .iter()
                    .map(|(&l, e)| (e.nodes.iter().filter_map(|n| rank.get(n).copied()).min().unwrap_or(usize::MAX), l))
                    .collect();
                elements.sort_unstable();
                (nodes, elements.into_iter().map(|(_, l)| l).collect())
            }
            Some("body") => {
                let bodies = element_bodies(model);
                // Elements by body name (untagged last), then label; nodes as their first
                // element in that order meets them, then the nodes no element uses.
                let mut elements: Vec<(Option<&String>, u32)> = model.elements().keys().map(|&l| (bodies.get(&l), l)).collect();
                elements.sort_by(|a, b| (a.0.is_none(), a.0, a.1).cmp(&(b.0.is_none(), b.0, b.1)));
                let elements: Vec<u32> = elements.into_iter().map(|(_, l)| l).collect();
                let mut seen = BTreeSet::new();
                let mut nodes: Vec<u32> = Vec::with_capacity(model.nodes().len());
                for e in &elements {
                    for &n in model.element(*e).map(|x| x.nodes.as_slice()).unwrap_or(&[]) {
                        if seen.insert(n) {
                            nodes.push(n);
                        }
                    }
                }
                let mut free: Vec<u32> = model.nodes().keys().copied().filter(|n| !seen.contains(n)).collect();
                free.sort_unstable();
                nodes.extend(free);
                (nodes, elements)
            }
            Some(_) => return Err(EngineError::param("by", ProblemKind::NotAChoice)),
        };
        let n = if nodes_too { from_start(&node_order, start).ok_or_else(out_of_range)? } else { BTreeMap::new() };
        let e = if elements_too { from_start(&element_order, start).ok_or_else(out_of_range)? } else { BTreeMap::new() };
        (n, e)
    };
    let before = bandwidth(model, &|l| l);
    s.model.renumber(&node_map, &element_map)?;
    let after = bandwidth(&s.model, &|l| l);
    let range = |labels: &mut dyn Iterator<Item = u32>| -> Value {
        let (lo, hi) = labels.fold((u32::MAX, 0), |(lo, hi), l| (lo.min(l), hi.max(l)));
        if hi == 0 { Value::Null } else { json!([lo, hi]) }
    };
    Ok(json!({
        "nodes": node_map.len(),
        "elements": element_map.len(),
        "node_labels": range(&mut s.model.nodes().keys().copied()),
        "element_labels": range(&mut s.model.elements().keys().copied()),
        "bandwidth": [before, after],
        "summary": format!("renumbered {} nodes and {} elements; bandwidth {before} → {after}", node_map.len(), element_map.len()),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 2 × 1 × 1 bar of hex8 blocks: whole, or as two halves that do not share their nodes
    /// on x = 1.
    fn bar(halves: bool) -> Session {
        let mut s = Session::default();
        let block = |body: &str, x0: f64, x1: f64, n: u32| {
            let c = [[x0, 0., 0.], [x1, 0., 0.], [x1, 1., 0.], [x0, 1., 0.], [x0, 0., 1.], [x1, 0., 1.], [x1, 1., 1.], [x0, 1., 1.]];
            json!({"body": body, "kind": "hex8", "corners": c, "divisions": [n, 2, 2], "merge": 0})
        };
        if halves {
            s.execute("mesh.block", &block("A", 0.0, 1.0, 2)).unwrap();
            s.execute("mesh.block", &block("B", 1.0, 2.0, 2)).unwrap();
        } else {
            s.execute("mesh.block", &block("A", 0.0, 2.0, 4)).unwrap();
        }
        s
    }

    fn positions(s: &Session) -> Vec<[u64; 3]> {
        let mut p: Vec<[u64; 3]> = s.model.nodes().values().map(|n| n.xyz.map(|x| (x + 0.0).to_bits())).collect();
        p.sort_unstable();
        p
    }

    /// Two halves meshed apart are found coincident and cracked, and merging them gives the
    /// whole bar's nodes: no crack left, and undo brings the halves back.
    #[test]
    fn merging_two_halves_gives_the_whole() {
        let whole = bar(false);
        let mut s = bar(true);
        s.execute(
            "model.set.create",
            &json!({"name": "seam", "nodes": s.model.nodes().iter().filter(|(_, n)| n.xyz[0] == 0.001).map(|(l, _)| *l).collect::<Vec<_>>()}),
        )
        .unwrap();
        let d = s.execute("mesh.check.duplicates", &json!({})).unwrap();
        assert_eq!(d["counts"]["coincident"], 9, "{d}");
        assert_eq!(d["counts"]["cracks"], 8, "{d}");
        assert_eq!(d["counts"]["free_edges"], 0, "{d}");
        assert_eq!(s.model.nodes().len(), whole.model.nodes().len() + 9);

        let r = s.execute("mesh.node.merge", &json!({})).unwrap();
        assert_eq!(r["merged"], 9, "{r}");
        assert_eq!(positions(&s), positions(&whole));
        // The seam set keeps its nine positions, on the survivors.
        assert_eq!(s.model.set("seam").unwrap().nodes.len(), 9);
        let d = s.execute("mesh.check.duplicates", &json!({})).unwrap();
        assert_eq!((d["counts"]["coincident"].as_u64(), d["counts"]["cracks"].as_u64()), (Some(0), Some(0)), "{d}");
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.model.nodes().len(), whole.model.nodes().len() + 9);

        // Keeping the highest labels; only a set's nodes; hostile parameters change nothing.
        let r = s.execute("mesh.node.merge", &json!({"keep": "highest", "set": "seam"})).unwrap();
        assert!(r["pairs"].as_array().unwrap().iter().all(|p| p[0].as_u64() < p[1].as_u64()), "{r}");
        s.execute("edit.undo", &json!({})).unwrap();
        for bad in [json!({"tolerance": -1}), json!({"tolerance": "x"}), json!({"keep": "middle"}), json!({"set": "nope"})] {
            assert!(s.execute("mesh.node.merge", &bad).is_err(), "{bad}");
        }
        assert!(merge_preview(&s, &json!({})).is_some_and(|p| p.len() == 9));
    }

    /// Renumbering keeps the analysis: the cantilever solves to the same displacements, point
    /// for point, after its nodes take the bandwidth order from 1000 and its elements are
    /// shifted; the deck writes the new labels.
    #[test]
    fn renumbering_keeps_the_solution() {
        let mut s = Session::default();
        s.run_script(&crate::sample::cantilever()).unwrap();
        let at = json!([[100.0, 10.0, 10.0], [50.0, 0.0, 20.0], [25.0, 20.0, 0.0]]);
        let probe = |s: &mut Session| {
            s.execute("solve.static", &json!({})).unwrap();
            s.execute("solve.probe", &json!({"step": "Load", "at": at})).unwrap()["nodes"].clone()
        };
        let before = probe(&mut s);
        let bw = bandwidth(&s.model, &|l| l);
        let r = s.execute("mesh.renumber", &json!({"what": "nodes", "by": "bandwidth", "start": 1000})).unwrap();
        assert_eq!(r["node_labels"][0], 1000, "{r}");
        assert!(r["bandwidth"][1].as_u64().unwrap() <= u64::from(bw), "{r}");
        s.execute("mesh.renumber", &json!({"what": "elements", "offset": 500})).unwrap();
        assert_eq!(s.model.elements().keys().next().copied(), Some(501));
        let after = probe(&mut s);
        for (a, b) in before.as_array().unwrap().iter().zip(after.as_array().unwrap()) {
            for k in 0..3 {
                let (x, y) = (a["displacement"][k].as_f64().unwrap(), b["displacement"][k].as_f64().unwrap());
                // Equal to roundoff of the tip deflection (0.24 mm).
                assert!((x - y).abs() <= 1e-9, "{a} vs {b}");
            }
        }
        let deck = fieldcraft_decks::abaqus::write(&s.model, fieldcraft_model::units::Preset::MmTonneS).unwrap();
        assert!(deck.lines().any(|l| l.starts_with("1000,")) && deck.lines().any(|l| l.starts_with("501,")), "{deck}");
        // Body by body and compact, then refusals that change nothing.
        s.execute("mesh.renumber", &json!({"by": "body"})).unwrap();
        s.execute("mesh.renumber", &json!({})).unwrap();
        assert_eq!(s.model.nodes().keys().next().copied(), Some(1));
        let model = s.model.clone();
        for bad in [
            json!({"offset": -1}),
            json!({"start": 99_999_999}),
            json!({"start": 0}),
            json!({"by": "random"}),
            json!({"what": "faces"}),
            json!({"offset": 1, "start": 3}),
        ] {
            assert!(s.execute("mesh.renumber", &bad).is_err(), "{bad}");
        }
        assert_eq!(s.model, model);
    }
}
