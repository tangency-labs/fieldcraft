//! Material frames that preserve FieldCraft's shell discretization in external formats.
use crate::{Error, Result};
use fieldcraft_model::{Elastic, ElementKind, MaterialOrientation, Model, Resolved, SectionKind, SectionOrientation, Set};
use std::collections::BTreeSet;

/// Keep the original sets for loads and constraints; only section assignments split. Each
/// default shell frame is then explicit, independent of Abaqus's global-axis default.
pub(crate) fn explicit_shell_frames(model: &Model, resolved: &Resolved, format: &'static str) -> Result<Option<Model>> {
    let selected: Vec<_> = model
        .sections()
        .iter()
        .filter(|(_, s)| {
            let shell = matches!(s.kind, SectionKind::CompositeShell { .. } | SectionKind::Shell { .. });
            let cylinder = s.orientation.is_some_and(|o| matches!(o.frame, MaterialOrientation::Cylindrical { .. }));
            shell
                && (cylinder
                    || (s.orientation.is_none()
                        && (matches!(s.kind, SectionKind::CompositeShell { .. })
                            || s.material_names()
                                .any(|n| matches!(model.material(n).and_then(|m| m.elastic), Some(Elastic::EngineeringConstants { .. }))))))
        })
        .collect();
    if selected.is_empty() {
        return Ok(None);
    }
    let additional = selected.iter().try_fold(0usize, |count, (name, s)| {
        let set = resolved
            .sets
            .get(&s.set)
            .ok_or_else(|| Error::Model(fieldcraft_model::Error::Invalid(format!("section `{name}` has no resolved set"))))?;
        Ok::<_, Error>(count.saturating_add(set.elements.len()))
    })?;
    if model.sets().len().saturating_add(additional) > fieldcraft_model::MAX_DEFINITIONS
        || model.sections().len().saturating_add(additional) > fieldcraft_model::MAX_DEFINITIONS
    {
        return Err(Error::Limit("explicit shell material frames would exceed the section/set definition cap".into()));
    }
    let mut out = model.clone();
    let mut names: BTreeSet<String> = model.sets().keys().chain(model.sections().keys()).map(|n| n.to_ascii_uppercase()).collect();
    for (name, section) in selected {
        let set = resolved
            .sets
            .get(&section.set)
            .ok_or_else(|| Error::Model(fieldcraft_model::Error::Invalid(format!("section `{name}` has no resolved set"))))?;
        out.remove_section(name)?;
        for (index, label) in set.elements.iter().enumerate() {
            if index.is_multiple_of(256) {
                fieldcraft_model::progress::check()?;
            }
            let element =
                model.element(*label).ok_or_else(|| Error::Model(fieldcraft_model::Error::Invalid(format!("element {label} is missing"))))?;
            if !matches!(element.kind, ElementKind::Tri3 | ElementKind::Quad4) {
                return Err(Error::Unsupported {
                    format,
                    what: format!(
                        "section `{name}` element {label}: an implicit material frame is defined only for tri3/quad4 shells; supply an explicit orientation"
                    ),
                });
            }
            let xyz = element
                .nodes
                .iter()
                .map(|n| model.node(*n).map(|n| n.xyz).ok_or_else(|| Error::Model(fieldcraft_model::Error::Invalid(format!("node {n} is missing")))))
                .collect::<Result<Vec<_>>>()?;
            let p = |i: usize| xyz.get(i).copied().unwrap_or_default();
            let sub = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| a[i] - b[i]);
            let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
            // Same geometry definition as FEM Shell::new, independently evaluated here to
            // keep the deck layer independent of the FEM layer. Quad4 uses opposite-edge
            // midpoints, not its first edge; triangles use the first edge.
            let (normal, along) = if element.kind == ElementKind::Tri3 {
                (cross(sub(p(1), p(0)), sub(p(2), p(0))), sub(p(1), p(0)))
            } else {
                (cross(sub(p(2), p(0)), sub(p(3), p(1))), std::array::from_fn(|i| (p(1)[i] + p(2)[i]) - (p(0)[i] + p(3)[i])))
            };
            let (frame, reference_axis) = if let Some(orientation) = section.orientation {
                let centre = std::array::from_fn(|i| xyz.iter().map(|p| p[i] / xyz.len() as f64).sum());
                let axes = axes_at(orientation.frame, centre)
                    .map_err(|e| Error::Unsupported { format, what: format!("section `{name}` element {label}: {e}") })?;
                (MaterialOrientation::Rectangular { axis1: axes[0], axis2: axes[1] }, orientation.shell_reference_axis)
            } else {
                (MaterialOrientation::Rectangular { axis1: along, axis2: cross(normal, along) }, 1)
            };
            frame
                .shell_axes_at([0.0; 3], normal, reference_axis)
                .map_err(|e| Error::Unsupported { format, what: format!("section `{name}` element {label}: {e}") })?;
            frame.validate().map_err(|e| Error::Unsupported { format, what: format!("section `{name}` element {label}: {e}") })?;
            let generated = (0..=fieldcraft_model::MAX_DEFINITIONS)
                .map(|i| format!("FC_MAT_{label}_{i}"))
                .find(|candidate| names.insert(candidate.to_ascii_uppercase()))
                .ok_or_else(|| Error::Limit("cannot name a shell material frame".into()))?;
            out.add_set(&generated, Set { elements: vec![*label], ..Set::default() })?;
            let mut single = section.clone();
            single.set = generated.clone();
            single.orientation = Some(SectionOrientation { frame, shell_reference_axis: reference_axis });
            out.add_section(&generated, single)?;
        }
    }
    Ok(Some(out))
}

/// Writer models may already be expressed in millimetres or inches. Translate and scale the
/// evaluation point for a cylinder before calling the model's direction-only calculation,
/// so the SI coordinate cap does not accidentally apply to those expressed coordinates.
fn axes_at(frame: MaterialOrientation, position: [f64; 3]) -> fieldcraft_model::Result<[[f64; 3]; 3]> {
    match frame {
        MaterialOrientation::Rectangular { .. } => frame.axes_at([0.0; 3]),
        MaterialOrientation::Cylindrical { origin, axis } => {
            let d: [f64; 3] = std::array::from_fn(|i| position[i] - origin[i]);
            let scale = d.iter().map(|v| v.abs()).fold(1.0, f64::max);
            MaterialOrientation::Cylindrical { origin: [0.0; 3], axis }.axes_at(d.map(|v| v / scale))
        }
    }
}
