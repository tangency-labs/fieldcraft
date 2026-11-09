//! Constant normal element foundations, on the reference face in small-displacement steps.
//! Public syntax and configuration dependence:
//! https://abaqus.uclouvain.be/English/SIMACAEKEYRefMap/simakey-r-foundation.htm
//! https://abaqus.uclouvain.be/English/SIMACAEMODRefMap/simamod-c-foundation.htm
use super::{Limits, Read, canonical, check_work, error, history, label, lex, number, quantity, row_fields};
use crate::{Error, Result};
use fieldcraft_model::units::{Preset, Quantity};
use fieldcraft_model::{DisplaySystem, ElementFace, ElementKind, Foundation, Set};
use std::collections::BTreeSet;

pub(super) fn import(text: &str, read: &mut Read, preset: Preset, limits: Limits) -> Result<()> {
    let blocks = lex::blocks(text, super::MAX_RECORDS)?;
    if !blocks.iter().any(|b| b.keyword == "FOUNDATION") {
        return Ok(());
    }
    let mut model = read.model.clone();
    let mut consumed = BTreeSet::new();
    let mut names: BTreeSet<String> = model.sets().keys().map(|n| canonical(n)).collect();
    let mut in_step = false;
    let mut members = model.sets().values().fold(0usize, |n, s| {
        n.saturating_add(s.nodes.len()).saturating_add(s.elements.len()).saturating_add(s.faces.len()).saturating_add(s.edges.len())
    });
    let apply = (|| -> Result<()> {
        for b in &blocks {
            if b.keyword == "STEP" {
                in_step = true;
            }
            if b.keyword == "ENDSTEP" {
                in_step = false;
            }
            if b.keyword != "FOUNDATION" {
                continue;
            }
            if in_step {
                return Err(history::unsupported(b.line, "FOUNDATION is persistent model data, not step activation/replacement"));
            }
            history::allow(b, &[])?;
            if b.rows.is_empty() {
                return Err(error(b.line, "FOUNDATION needs element/set, Fn, stiffness per area"));
            }
            for row in &b.rows {
                let fields = row_fields(row.text, row.line)?;
                let [target, kind, stiffness]: [String; 3] =
                    fields.try_into().map_err(|_| error(row.line, "FOUNDATION needs exactly element/set, Fn, stiffness per area"))?;
                let id = canonical(&kind)
                    .strip_prefix('F')
                    .and_then(|v| v.parse::<u8>().ok())
                    .filter(|n| (1..=4).contains(n))
                    .ok_or_else(|| history::unsupported(row.line, format!("FOUNDATION type {kind}: only Tet4 F1..F4 normal area stiffness")))?;
                let elements = if target.bytes().next().is_some_and(|v| v.is_ascii_digit()) {
                    vec![label(&target, row.line)?]
                } else {
                    model
                        .sets()
                        .iter()
                        .find(|(n, _)| canonical(n) == canonical(&target))
                        .map(|(_, s)| s.elements.clone())
                        .ok_or_else(|| error(row.line, format!("FOUNDATION missing element set {target}")))?
                };
                if elements.is_empty() {
                    return Err(error(row.line, format!("FOUNDATION empty element target {target}")));
                }
                members = members.saturating_add(elements.len());
                if members > limits.members {
                    return Err(error(row.line, "FOUNDATION face expansion limit reached"));
                }
                let mut faces = Vec::with_capacity(elements.len());
                for (i, element) in elements.into_iter().enumerate() {
                    check_work(i)?;
                    let e = model.element(element).ok_or_else(|| error(row.line, format!("FOUNDATION missing element {element}")))?;
                    if e.kind != ElementKind::Tet4 {
                        return Err(history::unsupported(
                            row.line,
                            format!("FOUNDATION element {element} {:?}: only Tet4 normal area stiffness is mapped", e.kind),
                        ));
                    }
                    faces.push(ElementFace { element, face: id });
                }
                let name = (0..=model.sets().len())
                    .map(|i| format!("InpFoundationFaces{}_{i}", row.line))
                    .find(|n| !names.contains(&canonical(n)))
                    .ok_or_else(|| error(row.line, "FOUNDATION set name exhaustion"))?;
                names.insert(canonical(&name));
                model.add_set(&name, Set { faces, ..Set::default() }).map_err(|e| error(row.line, e.to_string()))?;
                let normal = quantity(&DisplaySystem::new(preset), Quantity::ForcePerVolume, number(&stiffness, row.line)?, row.line)?;
                model
                    .add_foundation(
                        &format!("InpFoundation{}", row.line),
                        Foundation { set: name, normal, tangential: 0., normal_damping: 0., tangential_damping: 0., loss_factor: 0. },
                    )
                    .map_err(|e| error(row.line, e.to_string()))?;
            }
            consumed.insert(b.line);
        }
        Ok(())
    })();
    match apply {
        Err(Error::Unsupported { what, .. }) => {
            read.warnings.push(what);
            return Ok(());
        }
        Err(e) => return Err(e),
        Ok(()) => (),
    }
    if !consumed.is_empty() {
        read.model = model;
        read.unknown.retain(|u| !consumed.contains(&u.line));
        read.unknown_count = read.unknown.len();
        read.warnings.push("FOUNDATION: persistent Tet4 normal area stiffness in SI N/m^3; small-displacement static/frequency/direct harmonic only, not a finite-rotation follower foundation".into());
    }
    Ok(())
}
