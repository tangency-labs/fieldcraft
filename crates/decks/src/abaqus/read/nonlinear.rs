//! Automatic STATIC increments normalized by the step period (public *STATIC reference).
//! Only the first isotropic rod or Tet4 equilibrium is mapped here; continuation paths
//! need state transfer and must not silently become independent solves.
use super::history::{unsupported, value};
use super::{Block, error, label, number, row_fields};
use crate::Result;
use fieldcraft_model::{Elastic, ElementKind, Model, NonlinearStatic, Procedure};

pub(super) fn procedure(step: &Block<'_>, block: &Block<'_>, model: &Model, nonlinear_geometry: bool) -> Result<Procedure> {
    if model.sections().values().any(|s| s.formulation.enhanced_strain13) {
        return Err(unsupported(step.line, "C3D8I with NLGEOM=YES is not mapped"));
    }
    if model.sections().values().any(|s| s.formulation.hybrid || s.formulation.b_bar || s.formulation.incompatible_modes) {
        return Err(unsupported(step.line, "NLGEOM with hybrid, B-bar or incompatible formulations is not mapped"));
    }
    if model.elements().values().any(|e| !matches!(e.kind, ElementKind::Rod2 | ElementKind::Tet4))
        || model.sections().values().any(|s| {
            s.orientation.is_some()
                || !model.material(&s.material).is_some_and(|m| {
                    matches!(m.elastic, Some(Elastic::Isotropic { .. })) && (!nonlinear_geometry || m.plastic.is_none()) && m.hyperelastic.is_none()
                })
        })
    {
        return Err(unsupported(step.line, "NLGEOM=YES import currently requires isotropic elastic T3D2/T2D2 or C3D4"));
    }
    let fields = block.rows.first().map(|r| row_fields(r.text, r.line)).transpose()?.unwrap_or_default();
    let at = |i: usize, default: f64| -> Result<f64> {
        let v = fields.get(i).filter(|s| !s.is_empty()).map(|s| number(s, block.line)).transpose()?.unwrap_or(0.0);
        if v < 0.0 {
            return Err(error(block.line, "STATIC increments must be nonnegative"));
        }
        Ok(if v == 0.0 { default } else { v })
    };
    let period = at(1, 1.0)?;
    let initial = at(0, period)?.min(period) / period;
    let min = at(2, (initial * period).min(1e-5 * period))? / period;
    let max = at(3, period)?.min(period) / period;
    let n = NonlinearStatic {
        nonlinear_geometry,
        initial_increment: initial,
        min_increment: min,
        max_increment: max,
        max_increments: label(value(step, "INC", "100")?, step.line)?,
        ..NonlinearStatic::default()
    };
    n.validate().map_err(|e| error(block.line, e.to_string()))?;
    Ok(Procedure::NonlinearStatic(n))
}

#[cfg(test)]
mod tests {
    #[test]
    fn hybrid_nonlinear_import_is_named_and_retained() {
        let text = "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n*ELEMENT,TYPE=C3D4H,ELSET=All\n1,1,2,3,4\n*MATERIAL,NAME=M\n*ELASTIC\n1000,.3\n*SOLID SECTION,ELSET=All,MATERIAL=M\n*STEP,NLGEOM=YES\n*STATIC\n*END STEP\n";
        let read = crate::abaqus::read::read_analysis_limited(text, fieldcraft_model::units::Preset::MKgS, super::super::Limits::default()).unwrap();
        assert!(!read.analysis_complete);
        assert!(read.model.steps().is_empty());
        assert!(read.warnings.iter().any(|w| w.contains("hybrid")));
        assert!(read.unknown.iter().any(|u| u.keyword == "STEP"));
    }
}
