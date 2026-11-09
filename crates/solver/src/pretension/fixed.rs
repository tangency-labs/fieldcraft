//! Current-value control, evaluated from this solve's preceding equilibrium.
//! Rebuilding an elastic equilibrium is path independent; plastic restart history
//! is deliberately not approximated by the exported nodal stress field.
use crate::{Error, Limits, Result};
use fieldcraft_model::{Model, PretensionControl, Procedure};
use std::borrow::Cow;

pub(crate) fn resolve_fixed<'a>(model: &'a Model, name: &str, limits: &Limits) -> Result<Cow<'a, Model>> {
    let at = model.steps().iter().position(|s| s.name == name).ok_or_else(|| Error::Invalid(format!("no step `{name}`")))?;
    let mut needed = Vec::new();
    let mut i = at;
    // Strictly decreasing step indices, also bounded independently of input size.
    while model.steps().get(i).is_some_and(|s| s.pretension.values().any(|c| matches!(c, PretensionControl::Fixed))) {
        if needed.len() >= 1024 {
            return Err(Error::Limit("current-value pretension chain exceeds 1024 equilibria".into()));
        }
        let before = model
            .steps()
            .get(i)
            .and_then(|s| model.pretension_predecessor(s))
            .and_then(|previous| model.steps().iter().position(|s| s.name == previous.name))
            .ok_or_else(|| Error::Invalid(format!("step `{name}`: FIXED needs a preceding static equilibrium")))?;
        needed.push((i, before));
        i = before;
    }
    if needed.is_empty() {
        return Ok(Cow::Borrowed(model));
    }
    if model.materials().values().any(|m| m.plastic.is_some()) {
        return Err(Error::NotSupported("current-value pretension with plastic history requires a restart state".into()));
    }
    if model.contacts().values().any(|c| !c.friction.is_frictionless()) {
        return Err(Error::NotSupported("current-value pretension with frictional history requires a restart state".into()));
    }
    let mut out = model.clone();
    for (at, before) in needed.into_iter().rev() {
        let previous = out.steps().get(before).ok_or_else(|| Error::Internal("pretension predecessor".into()))?;
        let result = match previous.procedure {
            Procedure::NonlinearStatic(_) => crate::nonlinear::solve_nonlinear(&out, Some(&previous.name), limits)?,
            _ => crate::analysis::solve_static(&out, Some(&previous.name), limits)?,
        };
        if let Some(history) = &result.nonlinear
            && !history.increments.last().is_some_and(|increment| increment.converged)
        {
            return Err(Error::Invalid(format!("step `{}` has no terminal converged equilibrium for FIXED", previous.name)));
        }
        let mut step = out.steps().get(at).cloned().ok_or_else(|| Error::Internal("pretension locked step".into()))?;
        for (section, control) in &mut step.pretension {
            if matches!(control, PretensionControl::Fixed) {
                let value = result
                    .pretension
                    .iter()
                    .find(|r| r.name == *section)
                    .map(|r| r.adjustment)
                    .filter(|q| q.is_finite())
                    .ok_or_else(|| Error::Invalid(format!("step `{}` has no converged adjustment for `{section}`", previous.name)))?;
                *control = PretensionControl::Adjustment { value };
            }
        }
        out.replace_step(step)?;
    }
    Ok(Cow::Owned(out))
}
