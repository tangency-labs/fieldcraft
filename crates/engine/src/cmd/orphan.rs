//! Atomic orphan import: a reused tag does not establish semantic identity with an old mesh.
use std::collections::BTreeMap;
use std::path::Path;

use fieldcraft_io::gmsh::{self, Unmatched};
use fieldcraft_io::roots::Access;
use fieldcraft_model::Procedure;
use fieldcraft_model::units::Quantity;
use serde_json::{Value, json};

use crate::{EngineError, Result, Session};

pub(super) fn import(s: &mut Session, p: &Value) -> Result<Value> {
    let path = super::file_param(s, p, "path", Access::Read)?;
    if !Path::new(&path).extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("msh")) {
        return Err(EngineError::invalid("path", "only .msh (Gmsh MSH 4.1) is supported"));
    }
    let unmatched = match super::opt::<&str>(p, "unmatched")?.unwrap_or("refuse") {
        "refuse" => Unmatched::Refuse,
        "shell" => Unmatched::Shell,
        _ => return Err(EngineError::invalid("unmatched", "must be refuse or shell")),
    };
    let bytes = super::file::read_file(&path, gmsh::MAX_BYTES)?;
    // Binary bodies need not be UTF-8; recognize their ASCII header before decoding text.
    let mut header = bytes.split(|b| b.is_ascii_whitespace()).filter(|b| !b.is_empty());
    if header.next() == Some(b"$MeshFormat".as_slice()) {
        let version = header.next();
        if version == Some(b"4.1".as_slice()) && header.next() == Some(b"1".as_slice()) {
            return Err(EngineError::Params("binary MSH not supported yet".into()));
        }
    }
    let text = std::str::from_utf8(&bytes).map_err(|_| EngineError::Params("ASCII MSH must be valid UTF-8 text".into()))?;
    let imported = gmsh::read_msh(
        text,
        gmsh::Limits { nodes: s.limits.max_import_nodes, elements: s.limits.max_import_elements, members: s.limits.max_import_members },
        unmatched,
    )?;
    // A still-valid label could name another point or side. Do not carry old constraints or
    // silently replace them with same-named physical groups.
    if let Some((name, set)) = s.model.sets().iter().next() {
        let kind = if set.geometry.is_empty() { "label-based" } else { "geometry-based" };
        return Err(EngineError::Params(format!(
            "orphan import would replace mesh selections of {kind} set `{name}`; import into a new model or remove the set first"
        )));
    }
    for step in s.model.steps() {
        if let Procedure::NonlinearStatic(n) = &step.procedure
            && n.arc_length.as_ref().is_some_and(|a| a.displacement_limit.is_some())
        {
            return Err(EngineError::Params(format!(
                "orphan import would retarget the node-specific displacement limit of step `{}`; remove the limit first",
                step.name
            )));
        }
    }
    for (name, point) in s.model.reference_points() {
        if imported.model.node(point.label).is_some() {
            return Err(EngineError::Params(format!("imported node {} collides with reference point `{name}`", point.label)));
        }
    }
    let replaced = json!({"nodes": s.model.nodes().len(), "elements": s.model.elements().len()});
    let results: Vec<String> = s.solved_steps().map(str::to_string).collect();
    let mut model = s.model.without_mesh();
    for (&label, node) in imported.model.nodes() {
        model.set_node(label, node.xyz)?;
    }
    let mut counts = BTreeMap::<&str, usize>::new();
    for (&label, element) in imported.model.elements() {
        model.set_element(label, element.kind, element.nodes.clone())?;
        *counts.entry(element.kind.name()).or_default() += 1;
    }
    let mut sets = Vec::new();
    for (name, set) in imported.model.sets() {
        model.add_set(name, set.clone())?;
        sets.push(json!({"name":name,"nodes":set.nodes.len(),"elements":set.elements.len(),"faces":set.faces.len(),"edges":set.edges.len()}));
    }
    let reply = json!({
        "nodes":model.nodes().len(), "elements":counts, "sets":sets,
        "skipped_sections":imported.skipped_sections, "coordinate_unit":s.input_units().symbol(Quantity::Length),
        "replaced":replaced, "results_stale":results,
    });
    s.model = model;
    Ok(reply)
}

#[cfg(test)]
mod tests;
