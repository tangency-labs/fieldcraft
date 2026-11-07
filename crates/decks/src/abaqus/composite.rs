//! Independently mapped from Dassault Systemes' Abaqus 2026 keyword reference, *ORIENTATION
//! and *SHELL SECTION, and the model guide's "Orientations" shell projection algorithm.
//! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-orientation.htm
//! https://docs.software.vt.edu/abaqusv2026/English/SIMACAEKEYRefMap/simakey-r-shellsection.htm
use std::collections::BTreeSet;

use super::{Deck, Error, Result};
use fieldcraft_model::units::{Quantity, Unit};
use fieldcraft_model::{Elastic, MaterialOrientation, Model, Section, SectionKind, SectionOrientation};

pub(super) fn unsupported_items(model: &Model) -> BTreeSet<String> {
    let mut bad = BTreeSet::new();
    for (name, material) in model.materials() {
        if material.ply_strength.is_some() {
            bad.insert(format!("ply strengths of material `{name}`: *FAIL STRESS lacks S23; a Hashin damage-initiation law is not inferred from first-ply postprocessing strengths"));
        }
    }
    for (name, section) in model.sections() {
        let orthotropic =
            section.material_names().any(|m| matches!(model.material(m).and_then(|m| m.elastic), Some(Elastic::EngineeringConstants { .. })));
        if let SectionKind::CompositeShell { layup } = &section.kind {
            for (i, ply) in layup.plies.iter().enumerate() {
                if ply.points < 3 || ply.points.is_multiple_of(2) {
                    bad.insert(format!("section `{name}` ply {}: {} sampling points cannot preserve exact elastic integration with Abaqus Simpson quadrature (choose an odd count at least 3)", i + 1, ply.points));
                }
                if model
                    .material(&ply.material)
                    .is_some_and(|m| m.hyperelastic.is_some() || m.plastic.is_some() || m.damage.is_some() || m.erosion.is_some())
                {
                    bad.insert(format!(
                        "section `{name}` ply {}: inelastic material `{}` in a composite section is not written",
                        i + 1,
                        ply.material
                    ));
                }
            }
        } else if orthotropic && !section.kind.is_continuum() && !matches!(section.kind, SectionKind::Shell { .. }) {
            bad.insert(format!(
                "section `{name}`: orthotropic material `{}` in a {} section needs an explicit supported material orientation",
                section.material,
                section.kind.name()
            ));
        }
    }
    bad
}

/// Coordinate-vector scale has no physical meaning. Cylindrical endpoints use a scale large
/// enough that adding the axis to a remote origin does not erase the line direction.
pub(super) fn write_orientation(deck: &mut Deck, model: &Model, section: &Section, index: usize) -> Result<String> {
    let implicit_solid = matches!(section.kind, SectionKind::Solid {})
        && matches!(model.material(&section.material).and_then(|m| m.elastic), Some(Elastic::EngineeringConstants { .. }));
    let orientation = match section.orientation {
        Some(o) => o,
        None if implicit_solid => {
            SectionOrientation { frame: MaterialOrientation::Rectangular { axis1: [1.0, 0.0, 0.0], axis2: [0.0, 1.0, 0.0] }, shell_reference_axis: 1 }
        }
        None => return Ok(String::new()),
    };
    let name = format!("FC_ORI_{}", index + 1);
    let (system, coordinates) = match orientation.frame {
        MaterialOrientation::Rectangular { axis1, axis2 } => ("RECTANGULAR", [axis1, axis2].concat()),
        MaterialOrientation::Cylindrical { origin, axis } => {
            let scale = origin.iter().map(|v| v.abs()).fold(1.0, f64::max);
            let length = axis.iter().map(|v| v.abs()).fold(0.0, f64::max);
            let end = std::array::from_fn::<_, 3, _>(|i| origin[i] + (axis[i] / length) * scale);
            ("CYLINDRICAL", [origin, end].concat())
        }
    };
    deck.line(&format!("*ORIENTATION, NAME={name}, SYSTEM={system}, DEFINITION=COORDINATES"));
    deck.numbers_line("", &coordinates)?;
    if matches!(section.kind, SectionKind::Shell { .. } | SectionKind::CompositeShell { .. }) {
        // Abaqus projects the cyclic successor of its rotation axis, not that axis itself.
        let rotation_axis = if orientation.shell_reference_axis == 1 { 3 } else { orientation.shell_reference_axis - 1 };
        deck.line(&format!("{rotation_axis}, 0.0"));
    }
    Ok(format!(", ORIENTATION={name}"))
}

pub(super) fn degrees(radians: f64) -> Result<f64> {
    Unit::parse("deg")
        .and_then(|unit| Quantity::Angle.express(radians, &unit))
        .map_err(|e| Error::Model(fieldcraft_model::Error::Invalid(e.to_string())))
}
