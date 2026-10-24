//! Local records for ordinary solver failures (#654), without changing solve outcomes.

use std::io::Write;
use std::sync::Arc;

use fieldcraft_model::progress::Progress;
use serde_json::{Value, json};

use crate::{EngineError, Result, Session, crash};

pub(crate) fn is_solve(id: &str) -> bool {
    matches!(
        id,
        "solve.run" | "solve.static" | "solve.nonlinear" | "solve.heat" | "solve.frequency" | "solve.buckle" | "solve.explicit" | "solve.fluid"
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

fn solver_failure(e: &EngineError, device_observed: bool) -> bool {
    use fieldcraft_explicit::Error as E;
    use fieldcraft_lbm::Error as F;
    use fieldcraft_solver::Error as S;
    match e {
        EngineError::Context { source, .. } => solver_failure(source, device_observed),
        EngineError::Solver(S::Limit(_) | S::NotConverged(_) | S::Singular(_) | S::Internal(_) | S::Element { .. })
        | EngineError::Explicit(E::Limit(_) | E::Unstable(_))
        | EngineError::Fluid(F::Limit(_) | F::Unstable(_)) => true,
        // The existing GPU backend returns device errors in this documented wrapper.
        EngineError::Fluid(F::Invalid(message)) => {
            device_observed || message.starts_with("GPU:") || message.starts_with("no GPU backend is available")
        }
        _ => false,
    }
}

pub(crate) fn finish(
    s: &Session,
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
        Err(e) if solver_failure(e, device.is_some()) && !e.is_cancelled() => (e.to_string(), Value::Null),
        Ok(v) => match stopped(v) {
            Some((summary, why)) => (why.into(), summary.clone()),
            None => return result,
        },
        _ => return result,
    };
    let explicit_dir = std::env::var_os("FIELDCRAFT_CRASH_DIR").filter(|p| !p.is_empty()).is_some();
    if !s.solver_diagnostics && !crash::enabled() && !explicit_dir {
        return result;
    }
    let Some(dir) = crash::default_dir().and_then(|d| s.file_access.check(&d, fieldcraft_io::roots::Access::Write).ok()) else {
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
    let path = path.to_string_lossy().into_owned();
    let _ = writeln!(std::io::stderr().lock(), "solver diagnostic: {path}");
    match result {
        Err(source) => Err(EngineError::context(format!("solver diagnostic saved to {path}"), source)),
        Ok(mut v) => {
            v["diagnostic_report"] = json!(path);
            Ok(v)
        }
    }
}
