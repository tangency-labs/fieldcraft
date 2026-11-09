//! Bounds on mappings that require more than an element label.
use super::{Elastic, Error, FORMAT, Model, Procedure, Result};

pub(super) fn enhanced(model: &Model, procedure: Option<&Procedure>) -> Result<()> {
    for (name, section) in model.sections().iter().filter(|(_, s)| s.formulation.enhanced_strain13) {
        let material = model.material(&section.material);
        if section.orientation.is_some()
            || !material.is_some_and(|m| {
                matches!(m.elastic, Some(Elastic::Isotropic { .. }))
                    && m.plastic.is_none()
                    && m.hyperelastic.is_none()
                    && m.damage.is_none()
                    && m.erosion.is_none()
                    && m.thermal.is_none()
            })
        {
            return Err(Error::Unsupported {
                format: FORMAT,
                what: format!(
                    "C3D8I section `{name}`: only constant isotropic elasticity without orientation or thermal/nonlinear material is mapped"
                ),
            });
        }
        if procedure.is_some_and(|p| !matches!(p, Procedure::Static {}))
            || model.steps().iter().any(|s| !matches!(s.procedure, Procedure::Static {}) || s.base_state.is_some())
            || !model.initial_fields().is_empty()
        {
            return Err(Error::Unsupported {
                format: FORMAT,
                what: format!("C3D8I section `{name}`: only small-strain static analysis without initial fields or base_state is mapped"),
            });
        }
    }
    Ok(())
}

/// Ordinary equilibrium followed immediately by frequency extraction; no pretension writer.
pub(super) fn ordinary_base_states(model: &Model) -> Result<()> {
    if let Some(name) = model.pretensions().keys().next() {
        return Err(Error::Unsupported { format: FORMAT, what: format!("pretension section `{name}`") });
    }
    for (i, step) in model.steps().iter().enumerate() {
        if let Some(base) = &step.base_state {
            let previous = i.checked_sub(1).and_then(|j| model.steps().get(j));
            if !matches!(step.procedure, Procedure::Frequency { .. })
                || !previous.is_some_and(|s| s.name == *base && matches!(s.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_)))
            {
                return Err(Error::Unsupported {
                    format: FORMAT,
                    what: format!("base_state of step `{}`: only immediately preceding ordinary equilibrium to FREQUENCY is mapped", step.name),
                });
            }
        }
    }
    Ok(())
}
