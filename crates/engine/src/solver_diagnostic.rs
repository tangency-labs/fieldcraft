//! Local records of ordinary solver failures (#654): a solve that failed (no convergence, a
//! singular system, an unstable or diverged run, a failed device) or stopped at a cap before
//! the end of its step leaves a small record in the folder its session names
//! ([`Session::solver_diagnostics`]). A refusal of the request (a work cap such as `max_dofs`,
//! a bad parameter) and a cancellation leave none. The outcome itself is never changed: an
//! error stays the same error, and a stopped result only gains the record's path
//! (`diagnostic_report`); hosts read the path of a failure's record with
//! [`Session::take_solver_diagnostic`].

use std::sync::Arc;

use fieldcraft_model::progress::Progress;
use serde_json::{Value, json};

use crate::{EngineError, Result, Session, crash};

pub(crate) fn is_solve(id: &str) -> bool {
    matches!(
        id,
        "solve.run"
            | "solve.static"
            | "solve.nonlinear"
            | "solve.heat"
            | "solve.frequency"
            | "solve.cavity"
            | "solve.buckle"
            | "solve.explicit"
            | "solve.fluid"
    )
}

fn stopped(v: &Value) -> Option<(&Value, &'static str)> {
    match v.pointer("/stop/reason").and_then(Value::as_str) {
        Some("diverged") => Some((v, "the steady solve diverged")),
        Some("max_time") => Some((v, "the steady solve reached its time cap before settling")),
        _ if v.get("completed") == Some(&Value::Bool(false)) => Some((v, "the explicit solve stopped before completing its requested time")),
        _ => None,
    }
}

/// A failure of the solver itself, not a refusal of the request: work caps (`Limit`), bad input
/// and unsupported features are the user's to change and are not recorded.
fn solver_failure(e: &EngineError) -> bool {
    use fieldcraft_explicit::Error as E;
    use fieldcraft_lbm::Error as F;
    use fieldcraft_solver::Error as S;
    match e {
        EngineError::Context { source, .. } => solver_failure(source),
        EngineError::Solver(S::NotConverged(_) | S::Singular(_) | S::Internal(_) | S::Element { .. })
        | EngineError::Explicit(E::Unstable(_))
        | EngineError::Fluid(F::Unstable(_) | F::Backend(_) | F::NoDevice) => true,
        _ => false,
    }
}

/// Record the outcome `result` of solve command `id` (`parent`: the `solve.run` it is a step of)
/// if it is a solver failure or a stop, in the session's diagnostics folder; `result` is
/// returned as it was, a stopped result with the record's path.
pub(crate) fn finish(
    s: &mut Session,
    id: &str,
    params: &Value,
    progress: &Arc<Progress>,
    before: u64,
    parent: Option<&str>,
    result: Result<Value>,
) -> Result<Value> {
    // Each actual child of solve.run records at its own boundary. The aggregate may have
    // advanced to another step/backend by the time it returns and must not record again.
    if id == "solve.run" {
        return result;
    }
    let device = progress.backend_diagnostic();
    let (message, summary) = match &result {
        Err(e) if solver_failure(e) && !e.is_cancelled() => (e.to_string(), Value::Null),
        Ok(v) => match stopped(v) {
            Some((summary, why)) => (why.into(), summary.clone()),
            None => return result,
        },
        _ => return result,
    };
    let Some(dir) = s.solver_diagnostics.clone().and_then(|d| s.file_access.check(&d, fieldcraft_io::roots::Access::Write).ok()) else {
        return result;
    };
    let state = progress.snapshot();
    let observed = progress.observations() != before;
    let selected = params.get("step").and_then(Value::as_str);
    let fluid = (id == "solve.fluid").then(|| fieldcraft_lbm::fluid_step(&s.model, selected).ok()).flatten();
    let cells = summary.get("cells").cloned().or_else(|| fluid.and_then(|(_, domain)| domain.cells().ok()).map(|n| json!(n)));
    let actual = summary.get("backend").cloned().or_else(|| (id != "solve.fluid" && id != "solve.run").then(|| json!("cpu")));
    let brief: serde_json::Map<String, Value> = ["step", "time", "end_time", "increments", "completed", "dofs", "cells", "backend"]
        .into_iter()
        .filter_map(|k| summary.get(k).filter(|v| !v.is_array() && !v.is_object()).map(|v| (k.to_string(), v.clone())))
        .collect();
    let context = json!({
        "requested_step": params.get("step"),
        "parent_command": parent,
        "progress": observed.then_some(&state),
        "model": {"nodes": s.model.nodes().len(), "elements": s.model.elements().len(), "dofs": summary.get("dofs"), "cells": cells},
        "backend": {"requested": params.get("backend"), "actual": actual, "observed_context": device},
        "summary": brief,
        "parameter_units": s.input_units(),
        "summary_units": s.units(),
        "progress_units": "SI, with the count's unit named by progress.unit",
        "checkpoint_path": null,
        "unavailable": "Fields not observed by this solve remain null; a failed CPU fluid path has no backend identity observation, and no checkpoint is inferred from the open model path.",
    });
    let Some(path) = crash::solver_error(&dir, id, params, &context, &message) else { return result };
    log::info!("solver diagnostic of {id}: {}", path.display());
    let shown = crash::shown_path(&path);
    s.solver_diagnostic = Some(path);
    result.map(|mut v| {
        v["diagnostic_report"] = json!(shown);
        v
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backend_failures_are_classified_by_type_not_message() {
        assert!(solver_failure(&EngineError::Fluid(fieldcraft_lbm::Error::Backend("adapter unavailable".into()))));
        assert!(solver_failure(&EngineError::Fluid(fieldcraft_lbm::Error::NoDevice)));
        assert!(!solver_failure(&EngineError::Fluid(fieldcraft_lbm::Error::Invalid("GPU: invalid input".into()))));
        assert!(!solver_failure(&EngineError::Fluid(fieldcraft_lbm::Error::Invalid("no GPU backend is available input".into()))));
    }

    #[test]
    fn refusals_of_the_request_are_not_solver_failures() {
        assert!(!solver_failure(&EngineError::Solver(fieldcraft_solver::Error::Limit("max_dofs".into()))));
        assert!(!solver_failure(&EngineError::Explicit(fieldcraft_explicit::Error::Limit("max_elements".into()))));
        assert!(!solver_failure(&EngineError::Fluid(fieldcraft_lbm::Error::Limit("max_cells".into()))));
        assert!(!solver_failure(&EngineError::Params("step `Load`: no such step".into())));
        assert!(solver_failure(&EngineError::context("step `Load`", fieldcraft_solver::Error::NotConverged("increment 3".into()))));
    }
}
