//! MSC Nastran 2022.1 QRG: PCOMP, MAT8, CQUAD4, CORD2R and PSOLID.
//! https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf
use super::{Card, Error, FORMAT, Result};
use fieldcraft_model::units::{Quantity, Unit};
use fieldcraft_model::{MaterialOrientation, Model, SectionKind, SectionOrientation};
use std::collections::BTreeSet;

pub(super) fn check(model: &Model) -> Result<()> {
    let mut bad = BTreeSet::new();
    for (name, m) in model.materials() {
        if m.ply_strength.is_some() {
            bad.insert(format!("ply strengths of material `{name}`: MAT8 has one in-plane shear allowable, not the independent S23 used by the model's Hashin criterion"));
        }
    }
    for (name, section) in model.sections() {
        if matches!(section.kind, SectionKind::Solid {})
            && section.orientation.is_some_and(|o| matches!(o.frame, MaterialOrientation::Cylindrical { .. }))
        {
            bad.insert(format!("material orientation of solid section `{name}`: spatial cylindrical evaluation in PSOLID is not verified; it cannot be replaced by one constant frame"));
        }
        if let SectionKind::CompositeShell { layup } = &section.kind {
            for (i, ply) in layup.plies.iter().enumerate() {
                if ply.points != 1 {
                    bad.insert(format!("composite shell section `{name}` ply {}: PCOMP cannot preserve {} through-ply sampling points (only one mid-ply result is supported)",i+1,ply.points));
                }
                if model
                    .material(&ply.material)
                    .is_some_and(|m| m.hyperelastic.is_some() || m.plastic.is_some() || m.damage.is_some() || m.erosion.is_some())
                {
                    bad.insert(format!("section `{name}` ply {}: nonlinear material `{}` is not a linear PCOMP layer", i + 1, ply.material));
                }
            }
        }
    }
    crate::unsupported(FORMAT, bad)
}

pub(super) fn degrees(radians: f64) -> Result<f64> {
    Unit::parse("deg").and_then(|u| Quantity::Angle.express(radians, &u)).map_err(|e| Error::Model(fieldcraft_model::Error::Invalid(e.to_string())))
}

pub(super) fn frame(orientation: SectionOrientation, shell: bool) -> Result<MaterialOrientation> {
    let axes = orientation.frame.axes_at([0.0; 3])?;
    let i = if shell { usize::from(orientation.shell_reference_axis - 1) } else { 0 };
    Ok(MaterialOrientation::Rectangular { axis1: axes[i], axis2: axes[(i + 1) % 3] })
}

pub(super) fn coordinate_card(cid: u32, frame: MaterialOrientation) -> Result<Card> {
    let [x, _, z] = frame.axes_at([0.0; 3])?;
    Ok(Card::new("CORD2R").int(cid).blank().reals([0.0, 0.0, 0.0]).reals(z).reals(x))
}
