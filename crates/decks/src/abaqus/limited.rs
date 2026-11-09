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

/// Base states the writer maps: frequency extraction and steady-state response as perturbations
/// of an earlier general step. The pretension chain itself is checked by `validate_pretension_chain`.
pub(super) fn mapped_base_states(model: &Model) -> Result<()> {
    for (index, step) in model.steps().iter().enumerate() {
        if step.linear_base {
            let latest = model.steps()[..index].iter().rev().find(|s| matches!(s.procedure, Procedure::Static {} | Procedure::NonlinearStatic(_)));
            if latest.map(|s| s.name.as_str()) != step.base_state.as_deref()
                || model.steps().iter().any(|s| s.base_state == step.base_state && !s.linear_base)
            {
                return Err(Error::Unsupported {
                    format: FORMAT,
                    what: format!(
                        "step `{}`: linear_base must reference the latest general step and cannot share it with a stress-stiffened perturbation",
                        step.name
                    ),
                });
            }
        }
        if step.base_state.is_some() && !matches!(step.procedure, Procedure::Frequency { .. } | Procedure::Harmonic(_)) {
            return Err(Error::Unsupported {
                format: FORMAT,
                what: format!("base_state of step `{}`: only frequency and steady-state response perturbations are mapped", step.name),
            });
        }
    }
    Ok(())
}
