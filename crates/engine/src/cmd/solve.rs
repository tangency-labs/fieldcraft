//! `solve.*`: run the model's steps with the implicit, the explicit or the lattice Boltzmann
//! solver, keep the results in the session, and read them back.

use crate::schema::ProblemKind;
use fieldcraft_explicit::{ExplicitResult, solve_explicit};
use fieldcraft_lbm::{FluidResult, Obstacle};
use fieldcraft_mesh::surface::{BoundarySurface, SurfaceOptions};
use fieldcraft_model::DisplaySystem;
use fieldcraft_model::FluidDomain;
use fieldcraft_model::units::Quantity;
use fieldcraft_parallel as parallel;
use fieldcraft_solver::{
    FrequencyOptions, Limits, Location, MassKind, StaticResult, StepResult as ImplicitResult, solve_buckle, solve_frequency, solve_heat,
    solve_nonlinear, solve_static,
};
use serde_json::{Value, json};

use super::{labels, opt, points, req};
use crate::requirement::Requirement;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session, StepResult};

pub const RUN: CommandSpec = CommandSpec {
    id: "solve.run",
    label: "Solve All Steps",
    params: r#"{} — solves every step in order with its procedure's solver (solve.static, solve.frequency, solve.explicit, solve.fluid with their default caps) and keeps the results; stops at the first step that fails, naming it"#,
    form: Some(&[]),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_all,
};

pub const STATIC: CommandSpec = CommandSpec {
    id: "solve.static",
    label: "Solve Static Step",
    params: r#"{"step": "Load", "max_dofs": 6000000, "max_factor_entries": 400000000, "threads": 8} — solves a linear static or coupled temperature–displacement step (the only one when `step` is left out) and keeps the result in the session; `max_dofs` may only be lowered; `max_factor_entries` defaults to 400000000 and may be deliberately raised up to 2000000000 (8 bytes per factor value; indices and workspace are additional); `threads` (1 to 256, by default the session's: `FIELDCRAFT_THREADS` or the cores) changes only the time, never the result"#,
    form: Some(super::forms::SOLVE_STATIC),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_static,
};

pub const NONLINEAR: CommandSpec = CommandSpec {
    id: "solve.nonlinear",
    label: "Solve Nonlinear Static Step",
    params: r#"{"step": "Inflate", "max_dofs": 6000000, "max_factor_entries": 400000000, "max_newton_iterations": 100000, "threads": 8} — solves a nonlinear static step (the only one when `step` is left out) in load increments with Newton–Raphson iterations and keeps the result in the session: the final displacements, reactions, Cauchy stresses and logarithmic strains (small strains without nonlinear geometry), the equivalent plastic strain of plastic materials (the largest in the summary), every increment's load factor, iterations and residuals (force, displacement and energy, relative to their tolerances' scales), and a displacement frame per converged increment, and with contact pairs (interaction.contact.pair: mortar contact by augmented Lagrangian, Coulomb friction by a return map) each contact's force, normal force, largest pressure, closed and slipping nodes, deepest penetration and augmentations (solve.probe gives each secondary node's pressure, gap and status, and with friction its shear traction and slip); `max_dofs` and `max_newton_iterations` may only be lowered (the latter over the whole step); `max_factor_entries` and `threads` as for solve.static"#,
    form: Some(super::forms::SOLVE_NONLINEAR),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_nonlinear,
};

pub const HEAT: CommandSpec = CommandSpec {
    id: "solve.heat",
    label: "Solve Heat Transfer Step",
    params: r#"{"step": "Heat", "max_dofs": 6000000, "max_factor_entries": 400000000, "threads": 8} — steady or transient heat conduction of a heat transfer step (the only one when `step` is left out): temperatures, heat flux and the heat flow at prescribed temperatures, kept in the session; `max_dofs` may only be lowered; `max_factor_entries` and `threads` as for solve.static"#,
    form: Some(super::forms::SOLVE_STATIC),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_heat,
};

pub const FREQUENCY: CommandSpec = CommandSpec {
    id: "solve.frequency",
    label: "Solve Frequency Step",
    params: r#"{"step": "Modes", "mass": "consistent|lumped", "shift": 0.0, "max_dofs": 6000000, "max_factor_entries": 400000000, "max_iterations": 3000, "threads": 8} — natural frequencies and mode shapes of a frequency step (its modes and range); `max_iterations` caps the Lanczos vectors; the result stays in the session; `max_factor_entries` and `threads` as for solve.static"#,
    form: Some(super::forms::SOLVE_FREQUENCY),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_frequency,
};

pub const BUCKLE: CommandSpec = CommandSpec {
    id: "solve.buckle",
    label: "Solve Buckling Step",
    params: r#"{"step": "Buckle", "max_dofs": 6000000, "max_factor_entries": 400000000, "max_iterations": 3000, "threads": 8} — linear buckling of a buckling step (the only one when `step` is left out): the step's loads as a linear preload, then the lowest load factors λ with (K + λ K_σ) φ = 0 and their mode shapes (largest component 1), or those above the step's shift (whose LDLᵀ factor counts the load factors below it); negative load factors found are reported; the result stays in the session and its modes display like frequency modes; `max_factor_entries` and `threads` as for solve.static"#,
    form: Some(super::forms::SOLVE_BUCKLE),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_buckle,
};

pub const EXPLICIT: CommandSpec = CommandSpec {
    id: "solve.explicit",
    label: "Solve Explicit Dynamic Step",
    params: r#"{"step": "Impact", "max_dofs": 6000000, "max_increments": 1000000000, "max_output_values": 50000000, "max_contact_faces": 10000000, "threads": 8} — transient response of an explicit dynamic step (the only one when `step` is left out) by central differences, with the model's contacts (interaction.contact.*) and the step's rigid walls; keeps its histories and field frames in the session and reports each contact's peak force, first and last contact, deepest penetration, friction dissipated and the nodes adjusted at the start; the caps may only be lowered (the step's own increment cap applies too); `threads` as for solve.static"#,
    form: Some(super::forms::SOLVE_EXPLICIT),
    kind: CommandKind::Model,
    requires: &[Requirement::Analysable, Requirement::ElementMaterial, Requirement::ElementSection, Requirement::Step],
    run: run_explicit,
};

pub const SETTINGS: CommandSpec = CommandSpec {
    id: "solve.settings",
    label: "Run Settings",
    params: r#"{"step": "Flow", "backend": "gpu" | "cpu" | "none", "fallback"?: "cpu" | "none"} — the backend a fluid step prefers to be solved on, saved with the model (outside the definition: the step and its results are unchanged), and the one it falls back on where this machine has none (`cpu` by default for a GPU preference; `none` refuses instead); `backend: "none"` forgets the preference. solve.fluid (and solve.run) without a `backend` of its own then uses it: on a machine with no suitable GPU it falls back with a `notice`. A machine overrides every preference with FIELDCRAFT_SOLVE_BACKEND=cpu|gpu. With no step, lists every step's settings"#,
    form: Some(super::forms::SOLVE_SETTINGS),
    kind: CommandKind::Session,
    requires: &[],
    run: settings,
};

pub const HISTORY: CommandSpec = CommandSpec {
    id: "solve.history",
    label: "Result Histories",
    params: r#"{"step": "Impact", "every": 1, "sets": ["top"], "node": 21} — the time histories of a solved explicit step: time, increment, energies, momentum, constraint force, the normal force across each contact (in `contacts` order) and the history nodes' displacement and velocity; or the frame histories of a nonlinear static step (one frame per converged increment, from load factor 0): `frames` (indices for results.frame), `load_factor` and, per set of `sets` (default the sets of the step's boundary conditions), the mean displacement of its nodes and the sum of their reactions per frame (a force–travel curve), with `node` that node's displacement too; `every` keeps one sample in so many (by default enough to return at most 10000)"#,
    form: Some(super::forms::SOLVE_HISTORY),
    kind: CommandKind::Model,
    requires: &[Requirement::Result],
    run: history,
};

pub const FLUID: CommandSpec = CommandSpec {
    id: "solve.fluid",
    label: "Solve Fluid Step",
    params: r#"{"step": "Flow", "backend": "cpu", "max_cells": 16000000, "max_steps": 10000000, "max_cell_updates": 1e12, "max_frame_bytes": 2e9, "threads": 8, "estimate": false} — runs a fluid step with the lattice Boltzmann solver (the only fluid step when `step` is left out): obstacles are triangulated from the open design or taken from imported surfaces; returns the time step, relaxation time, Mach number, mass, the largest speed of all frames (`max_speed`, at `max_speed_time`) and of the last frame (`final_max_speed`), over interior fluid cells (the open-boundary layers are excluded, convective outlets too, because a wake reaching a uniform pressure outlet distorts the last layer), the final force on every wall and, with heat, the thermal relaxation time, the fluid's temperature range and the heat into every wall, and per conducting solid (`solids`) its temperature range and mean, its power, the heat its faces give the fluid, its `peak` temperature (at `peak_time`: the steady value when the run reached steady state, else the highest over the frames), its `limit`, the `margin` (limit less peak; negative over it) and `ok`, and over them `components_ok` (null when no solid has a limit); with fans and porous zones, each one's volume flow and its pressure rise (a fan's operating point: `pressure_rise` from its curve, and `measured_pressure_rise`) or `pressure_drop`; the measured values are the static-pressure difference between the fluid cell layers just outside the zone, so they include the wall friction over its thickness plus one cell (final and mean over the last quarter of the run); `backend` is `cpu` (the default: double precision, the same result on every machine) or `gpu` (the step on the GPU `FIELDCRAFT_GPU` picks, single precision, within the backend's stated tolerance of the CPU; an error when the build has none); left out, it is FIELDCRAFT_SOLVE_BACKEND's when that is set, else the step's saved preference (solve.settings), which on a machine with no suitable GPU falls back as it says, with a `notice` (the reply's `backend_choice` says what chose it), else the CPU; `max_steps` may only be lowered; `max_cells` (16000000 by default), `max_cell_updates` (1e12) and `max_frame_bytes` (2e9) may be raised up to what the machine's memory bound holds (half its physical memory), and a run with a raised `max_cells` or `max_frame_bytes` must fit in that bound; the reply's `memory` gives the run's estimate (`flow_bytes`, `frame_bytes` at about 16 bytes per cell per frame, 20 with heat, `total_bytes`) against `bound_bytes`, the `frames_that_fit` and each cap's `ceilings`, and on the GPU the population buffer against the device's limit per buffer; `"estimate": true` returns only that, without solving; a run over a cap is refused before it starts, saying what to change; `threads` as for solve.static"#,
    form: Some(super::forms::SOLVE_FLUID),
    kind: CommandKind::Model,
    requires: &[Requirement::FluidDomain, Requirement::Step],
    run: run_fluid,
};

pub const PROBE: CommandSpec = CommandSpec {
    id: "solve.probe",
    label: "Probe Results",
    params: r#"{"step": "Load", "nodes": [1, 2] (or "set": "fillet", a set's nodes), "mode": 1, "frame": 0} — the stored result at nodes (with `max_von_mises`, the largest among them, for static steps): displacement, rotation, reaction, reaction moment, stress, principal stresses (the values, largest first, without their directions) and von Mises (static, with `elements` each element's results: solid, 2D and truss nodal stresses (and, in a nonlinear step, the equivalent plastic strain at the integration points of a plastic element), shell section forces (N, M, Q per length in the element's axes) and bottom, middle and top stresses, beam end forces (N, V_y, V_z, T, M_y, M_z in its axes)), the shape of mode `mode` with its rotations (frequency, from 1), or displacement and velocity in field frame `frame` (explicit, from 0; the last by default) at `nodes` (or at the nodes nearest to points: "at": [[x, y, z], ...], lowest label on a tie, each entry then giving its `xyz` and `distance`) and, with `elements` (explicit), each element's stress, von Mises stress, strain, equivalent plastic strain, damage, temperature and whether it has eroded; a coupled reference point's label gives its six displacements (translations, rotations) and reactions (forces, moments). Fluid steps: {"step": "Flow", "points": [[x, y, z]], "frame": 10 (default the last)} gives velocity, pressure and (with heat) temperature in the cells at the points; {"step": "Flow", "forces": ["x_min", "F2:top", "Body1/F2:top"]} gives those walls' force histories and {"step": "Flow", "heat": ["y_min"]} the heat into them over time (power, positive into the wall); {"step": "Flow", "line": {"from": [x, y, z], "to": [x, y, z], "count": 65 (default one point per cell spacing), "field": "velocity_x|velocity_y|velocity_z|speed|pressure|temperature|wall_heat_flux", "face": "x_min" (wall_heat_flux: the box wall the line lies on)}, "frame": 10} gives the field at `count` points evenly along the line (`xyz`; trilinear between the cell centres, null outside the fluid; within half a cell of a box face the nearest cell centre's value, so a temperature line run onto a held wall reads the first cell, not the wall) and its `max` and `min` with their `xyz`, located between the samples by the parabola through the extreme sample and its neighbours; `wall_heat_flux` is the heat per area into a box wall (positive into the wall): at a held temperature from the second-order one-sided difference through the wall and the two cell layers beside it, at a flux face its flux, at an adiabatic face zero"#,
    form: Some(super::forms::SOLVE_PROBE),
    kind: CommandKind::Model,
    requires: &[Requirement::Result],
    run: probe,
};

/// Largest explicit factor cap accepted by solve commands; the default remains 400M entries.
pub(super) const FACTOR_ENTRY_CEILING: usize = 2_000_000_000;

/// Most nodes one probe lists.
const MAX_PROBE: usize = 10_000;

/// Most force histories one fluid probe lists.
const MAX_PROBE_FORCES: usize = 100;

/// Most history samples one `solve.history` returns.
const MAX_HISTORY_SAMPLES: usize = 10_000;

/// A whole number given as an integer or as a number with no fraction (8e9); floats only up to
/// 2⁵³, beyond which they are not exact.
fn whole(v: &Value) -> Option<u64> {
    v.as_u64().or_else(|| v.as_f64().filter(|x| x.fract() == 0.0 && (0.0..=9.007_199_254_740_992e15).contains(x)).map(|x| x as u64))
}

/// An optional cap that may only lower the default.
fn cap<T: Copy + TryInto<u64> + TryFrom<u64> + std::fmt::Display>(p: &Value, key: &str, default: T) -> Result<T> {
    let most: u64 = default.try_into().unwrap_or(u64::MAX);
    match p.get(key) {
        None => Ok(default),
        Some(v) => whole(v)
            .filter(|n| (1..=most).contains(n))
            .and_then(|n| T::try_from(n).ok())
            .ok_or_else(|| EngineError::param(key, ProblemKind::NotAWholeNumber { min: 1, max: i64::try_from(most).unwrap_or(i64::MAX) })),
    }
}

/// The `threads` parameter: 1 to [`parallel::MAX_THREADS`], by default the session's.
fn threads(s: &Session, p: &Value) -> Result<usize> {
    match p.get("threads") {
        None => Ok(s.limits.threads.clamp(1, parallel::MAX_THREADS)),
        Some(v) => v
            .as_u64()
            .and_then(|n| usize::try_from(n).ok())
            .filter(|n| (1..=parallel::MAX_THREADS).contains(n))
            .ok_or_else(|| EngineError::Params(format!("`threads` must be a whole number from 1 to {}", parallel::MAX_THREADS))),
    }
}

fn limits(p: &Value) -> Result<Limits> {
    Ok(Limits {
        max_dofs: cap(p, "max_dofs", Limits::MAX_DOFS)?,
        max_factor_entries: if p.get("max_factor_entries").is_none() {
            Limits::MAX_FACTOR_ENTRIES
        } else {
            cap(p, "max_factor_entries", FACTOR_ENTRY_CEILING)?
        },
        max_lanczos_vectors: cap(p, "max_iterations", Limits::MAX_LANCZOS_VECTORS)?,
        max_newton_iterations: cap(p, "max_newton_iterations", Limits::MAX_NEWTON_ITERATIONS)?,
    })
}

/// The node with the value that `better` prefers, as `{"node", "value"}`.
fn extreme(values: impl Iterator<Item = (u32, f64)>, better: fn(f64, f64) -> bool) -> Value {
    match values.fold(None, |best: Option<(u32, f64)>, (l, v)| if best.is_none_or(|(_, b)| better(v, b)) { Some((l, v)) } else { best }) {
        Some((node, value)) => json!({"node": node, "value": value}),
        None => Value::Null,
    }
}

/// [`parallel::with_threads`], the caller's progress reporter carried onto the pool's thread.
fn with_threads<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
    let progress = fieldcraft_model::progress::current();
    parallel::with_threads(threads, move || fieldcraft_model::progress::scoped(progress, f))
}

fn run_all(s: &mut Session, _: &Value) -> Result<Value> {
    use fieldcraft_model::Procedure;
    type Solve = fn(&mut Session, &Value) -> Result<Value>;
    // Each solve with the command whose summary it gives (its units).
    let steps: Vec<(String, &'static str, Solve, &'static str)> = s
        .model
        .steps()
        .iter()
        .map(|st| {
            let (solve, id): (Solve, &str) = match st.procedure {
                Procedure::Static {} | Procedure::CoupledTemperatureDisplacement {} => (run_static, "solve.static"),
                Procedure::HeatTransfer { .. } => (run_heat, "solve.heat"),
                Procedure::Frequency { .. } => (run_frequency, "solve.frequency"),
                Procedure::Buckle { .. } => (run_buckle, "solve.buckle"),
                Procedure::ExplicitDynamic(_) => (run_explicit, "solve.explicit"),
                Procedure::Fluid { .. } => (run_fluid, "solve.fluid"),
                Procedure::NonlinearStatic(_) => (run_nonlinear, "solve.nonlinear"),
            };
            (st.name.clone(), st.procedure.name(), solve, id)
        })
        .collect();
    if steps.is_empty() {
        return Err(EngineError::Params("the model has no steps to solve: add one (model.step.create)".into()));
    }
    let mut done = Vec::new();
    let count = steps.len();
    for (index, (name, procedure, solve, id)) in steps.into_iter().enumerate() {
        fieldcraft_model::progress::step(&name, index + 1, count);
        fieldcraft_model::progress::check()?;
        let mut summary =
            solve(s, &json!({"step": name})).map_err(|e| if e.is_cancelled() { e } else { EngineError::Params(format!("step `{name}`: {e}")) })?;
        let display = s.units().clone();
        crate::units::convert_reply(id, &mut summary, crate::units::si(), &display)?;
        done.push(json!({"step": name, "procedure": procedure, "summary": summary}));
    }
    Ok(json!({"steps": done}))
}

/// A sum's component below this fraction of its terms' magnitudes is round-off.
const CANCELLED: f64 = 1e-9; // unit-factors: relative to the terms' magnitudes

/// A sum of nodal values with its round-off set to 0: a component below a billionth of the sum
/// of the terms' magnitudes is what cancelling left (a balanced load sums to ~1e-11, not 0).
fn settled(total: [f64; 3], terms: &[[f64; 3]]) -> [f64; 3] {
    let mut out = total;
    for (k, c) in out.iter_mut().enumerate() {
        let scale: f64 = terms.iter().filter_map(|t| t.get(k)).map(|x| x.abs()).filter(|x| x.is_finite()).sum();
        if c.abs() <= CANCELLED * scale {
            *c = 0.0;
        }
    }
    out
}

fn largest(values: impl Iterator<Item = (u32, f64)>) -> Value {
    extreme(values, |v, b| v > b)
}

fn smallest(values: impl Iterator<Item = (u32, f64)>) -> Value {
    extreme(values, |v, b| v < b)
}

fn run_static(s: &mut Session, p: &Value) -> Result<Value> {
    let limits = limits(p)?;
    let step = opt::<String>(p, "step")?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || solve_static(model, step.as_deref(), &limits))?;

    let summary = json!({
        "step": r.step,
        "dofs": r.stats.dofs,
        "equations": r.stats.equations,
        "factor_entries": r.stats.factor_entries,
        "max_displacement": largest(r.nodes.iter().zip(&r.displacement).map(|(l, u)| (*l, crate::results::magnitude(*u)))),
        "max_von_mises": largest(r.nodes.iter().zip(&r.von_mises).map(|(l, v)| (*l, *v))),
        "applied_total": settled(r.applied_total, &r.reaction),
        "reaction_total": settled(r.reaction_total, &r.reaction),
        "reference_points": r.reference_points.iter().map(|p| json!({
            "name": p.name,
            "label": p.label,
            "displacement": p.displacement,
            "reaction": p.reaction,
        })).collect::<Vec<_>>(),
        "max_temperature": r.temperature.as_ref().map(|t| largest(r.nodes.iter().copied().zip(t.iter().copied()))),
        "warnings": r.warnings,
    });
    s.store_result(StepResult::Implicit(ImplicitResult::Static(r)));
    Ok(summary)
}

fn run_nonlinear(s: &mut Session, p: &Value) -> Result<Value> {
    let limits = limits(p)?;
    let step = opt::<String>(p, "step")?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || solve_nonlinear(model, step.as_deref(), &limits))?;
    let history = r.nonlinear.as_ref();
    let summary = json!({
        "step": r.step,
        "dofs": r.stats.dofs,
        "equations": r.stats.equations,
        "factor_entries": r.stats.factor_entries,
        "increments": history.map(|h| &h.increments),
        "iterations": history.map(|h| h.iterations),
        "frames": history.map(|h| h.load_factors()),
        "load_factor": history.and_then(|h| h.frames.last()).map(|f| f.load_factor),
        "critical_points": history.map(|h| &h.critical_points),
        "contacts": r.contacts.iter().map(|c| json!({
            "name": c.name,
            "normal_force": c.normal_force,
            "force": c.force,
            "max_pressure": largest(c.nodes.iter().copied().zip(c.pressure.iter().copied())),
            "closed_nodes": c.closed(),
            "slipping_nodes": c.status.iter().filter(|s| **s == fieldcraft_solver::contact::Status::Slip).count(),
            "penetration": c.penetration,
            "augmentations": c.augmentations,
        })).collect::<Vec<_>>(),
        "max_displacement": largest(r.nodes.iter().zip(&r.displacement).map(|(l, u)| (*l, crate::results::magnitude(*u)))),
        "max_von_mises": largest(r.nodes.iter().zip(&r.von_mises).map(|(l, v)| (*l, *v))),
        "max_plastic_strain": largest(r.elements.iter().flat_map(|e| e.plastic_strain.iter().map(|v| (e.label, *v)))),
        "applied_total": settled(r.applied_total, &r.reaction),
        "reaction_total": settled(r.reaction_total, &r.reaction),
        "reference_points": r.reference_points.iter().map(|p| json!({
            "name": p.name,
            "label": p.label,
            "displacement": p.displacement,
            "reaction": p.reaction,
        })).collect::<Vec<_>>(),
        "warnings": r.warnings,
    });
    s.store_result(StepResult::Implicit(ImplicitResult::Static(r)));
    Ok(summary)
}

fn run_heat(s: &mut Session, p: &Value) -> Result<Value> {
    let limits = limits(p)?;
    let step = opt::<String>(p, "step")?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || solve_heat(model, step.as_deref(), &limits))?;
    let t = r.temperature();
    let summary = json!({
        "step": r.step,
        "dofs": r.stats.dofs,
        "equations": r.stats.equations,
        "factor_entries": r.stats.factor_entries,
        "frames": r.frames.len(),
        "end_time": r.times.last(),
        "min_temperature": smallest(r.nodes.iter().copied().zip(t.iter().copied())),
        "max_temperature": largest(r.nodes.iter().copied().zip(t.iter().copied())),
        "heat_supplied": r.applied_total,
        "heat_removed": r.reaction_total,
        "heat_flows": crate::results::heat_flows(&r, crate::units::si()),
        "warnings": r.warnings,
    });
    s.store_result(StepResult::Implicit(ImplicitResult::Heat(r)));
    Ok(summary)
}

fn run_frequency(s: &mut Session, p: &Value) -> Result<Value> {
    let limits = limits(p)?;
    let step = opt::<String>(p, "step")?;
    let mass = opt::<MassKind>(p, "mass")?.unwrap_or_default();
    let shift = opt::<f64>(p, "shift")?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || solve_frequency(model, step.as_deref(), &FrequencyOptions { mass, shift }, &limits))?;
    let summary = json!({
        "step": r.step,
        "dofs": r.stats.dofs,
        "equations": r.stats.equations,
        "lanczos_vectors": r.stats.lanczos_vectors,
        "total_mass": r.total_mass,
        "modes": r.modes.iter().enumerate().map(|(i, m)| json!({
            "mode": i + 1,
            "frequency": m.frequency,
            "eigenvalue": m.eigenvalue,
            "effective_mass": m.effective_mass,
        })).collect::<Vec<_>>(),
        "warnings": r.warnings,
    });
    s.store_result(StepResult::Implicit(ImplicitResult::Frequency(r)));
    Ok(summary)
}

fn run_buckle(s: &mut Session, p: &Value) -> Result<Value> {
    let limits = limits(p)?;
    let step = opt::<String>(p, "step")?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || solve_buckle(model, step.as_deref(), &limits))?;
    let summary = json!({
        "step": r.step,
        "dofs": r.stats.dofs,
        "equations": r.stats.equations,
        "lanczos_vectors": r.stats.lanczos_vectors,
        "modes": r.modes.iter().enumerate().map(|(i, m)| json!({"mode": i + 1, "load_factor": m.eigenvalue})).collect::<Vec<_>>(),
        "buckling": r.buckling,
        "warnings": r.warnings,
    });
    s.store_result(StepResult::Implicit(ImplicitResult::Frequency(r)));
    Ok(summary)
}

fn run_explicit(s: &mut Session, p: &Value) -> Result<Value> {
    use fieldcraft_explicit::Limits as ExplicitLimits;
    let limits = ExplicitLimits {
        max_dofs: cap(p, "max_dofs", ExplicitLimits::MAX_DOFS)?,
        max_increments: cap(p, "max_increments", ExplicitLimits::MAX_INCREMENTS)?,
        max_output_values: cap(p, "max_output_values", ExplicitLimits::MAX_OUTPUT_VALUES)?,
        max_contact_faces: cap(p, "max_contact_faces", ExplicitLimits::MAX_CONTACT_FACES)?,
    };
    let step = opt::<String>(p, "step")?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || solve_explicit(model, step.as_deref(), &limits))?;
    let summary = explicit_summary(&r);
    s.store_result(StepResult::Explicit(r));
    Ok(summary)
}

fn explicit_summary(r: &ExplicitResult) -> Value {
    json!({
        "step": r.step,
        "dofs": 3 * r.nodes.len(),
        "elements": r.elements.len(),
        "increments": r.increments,
        "time": r.time,
        "completed": r.completed,
        "initial_increment": r.initial_increment,
        "controlling_elements": r.controlling,
        "total_mass": r.total_mass,
        "energy": r.energy,
        "energy_balance": r.energy.relative_error(r.initial_kinetic),
        "frames": r.frames.len(),
        "history_samples": r.history.time.len(),
        "eroded": r.eroded,
        "contacts": r.contacts,
        "mass_scaling": r.mass_scaling,
        "warnings": r.warnings,
    })
}

fn history(s: &mut Session, p: &Value) -> Result<Value> {
    let step: String = req(p, "step")?;
    let result = s.result(&step).ok_or_else(|| EngineError::Params(format!("step `{step}` has no result; solve it first")))?;
    if let crate::results::Kind::Static(r) = result.kind()
        && r.nonlinear.is_some()
    {
        return nonlinear_history(s, &step, r, p);
    }
    let StepResult::Explicit(r) = result.as_ref() else {
        return Err(EngineError::Params(format!("step `{step}` is neither an explicit dynamic nor a nonlinear static step; it has no histories")));
    };
    let h = &r.history;
    let samples = h.time.len();
    let every = match opt::<u64>(p, "every")? {
        Some(0) => return Err(EngineError::param("every", ProblemKind::Below(1.0))),
        Some(k) => usize::try_from(k).unwrap_or(usize::MAX),
        None => samples.div_ceil(MAX_HISTORY_SAMPLES).max(1),
    };
    if samples.div_ceil(every) > MAX_HISTORY_SAMPLES {
        return Err(EngineError::Params(format!("{samples} samples in steps of {every} are more than {MAX_HISTORY_SAMPLES}; raise `every`")));
    }
    let pick = |v: &dyn Fn(usize) -> Value| -> Vec<Value> { (0..samples).step_by(every).map(v).collect() };
    let nodes: Vec<Value> = h
        .nodes
        .iter()
        .enumerate()
        .map(|(k, label)| {
            json!({
                "node": label,
                "displacement": pick(&|i| json!(h.displacement.get(i).and_then(|d| d.get(k)))),
                "velocity": pick(&|i| json!(h.velocity.get(i).and_then(|d| d.get(k)))),
            })
        })
        .collect();
    Ok(json!({
        "step": step,
        "every": every,
        "time": pick(&|i| json!(h.time.get(i))),
        "increment": pick(&|i| json!(h.increment.get(i))),
        "energies": pick(&|i| json!(h.energies.get(i))),
        "momentum": pick(&|i| json!(h.momentum.get(i))),
        "constraint_force": pick(&|i| json!(h.constraint_force.get(i))),
        "contact_force": pick(&|i| json!(h.contact_force.get(i))),
        "contacts": r.contacts.iter().map(|c| c.name.as_str()).collect::<Vec<_>>(),
        "nodes": nodes,
    }))
}

/// solve.history of a nonlinear static step: per frame the load factor and, per set, the mean
/// displacement of its nodes and their reaction sum.
fn nonlinear_history(s: &Session, step: &str, r: &fieldcraft_solver::StaticResult, p: &Value) -> Result<Value> {
    let frames = r.nonlinear.as_ref().map(|h| h.frames.as_slice()).unwrap_or_default();
    let every = match opt::<u64>(p, "every")? {
        Some(0) => return Err(EngineError::param("every", ProblemKind::Below(1.0))),
        Some(k) => usize::try_from(k).unwrap_or(usize::MAX),
        None => frames.len().div_ceil(MAX_HISTORY_SAMPLES).max(1),
    };
    let node = match p.get("node") {
        None | Some(Value::Null) => None,
        Some(v) => {
            let label = v.as_u64().and_then(|n| u32::try_from(n).ok()).ok_or_else(|| EngineError::invalid("node", "must be a node label"))?;
            let i =
                r.nodes.binary_search(&label).map_err(|_| EngineError::invalid("node", format!("node {label} has no result in step `{step}`")))?;
            Some((label, i))
        }
    };
    let picked: Vec<usize> = (0..frames.len()).step_by(every).collect();
    let mut sets = Vec::new();
    for name in super::results::reaction_sets(s, step, p)? {
        let resolved = s.model.resolve_set(&name)?;
        let at: Vec<usize> = resolved.nodes.iter().filter_map(|l| r.nodes.binary_search(l).ok()).collect();
        if at.is_empty() {
            return Err(EngineError::Params(format!("set `{name}` has no nodes with results in step `{step}`")));
        }
        let (mut displacement, mut reaction) = (Vec::with_capacity(picked.len()), Vec::with_capacity(picked.len()));
        for f in picked.iter().filter_map(|k| frames.get(*k)) {
            let (mut u, mut force) = ([0.0; 3], [0.0; 3]);
            for &i in &at {
                let (d, q) = (f.displacement.get(i).copied().unwrap_or_default(), f.reaction_at(i));
                for k in 0..3 {
                    u[k] += d[k];
                    force[k] += q[k];
                }
            }
            let n = at.len().max(1) as f64;
            displacement.push(u.map(|c| c / n));
            reaction.push(force);
        }
        sets.push(json!({"set": name, "nodes": at.len(), "displacement": displacement, "reaction": reaction}));
    }
    let mut out = json!({
        "step": step,
        "every": every,
        "frames": picked,
        "load_factor": picked.iter().map(|k| frames.get(*k).map(|f| f.load_factor)).collect::<Vec<_>>(),
        "sets": sets,
    });
    if let Some((label, i)) = node {
        out["node"] =
            json!({"node": label, "displacement": picked.iter().map(|k| frames.get(*k).and_then(|f| f.displacement.get(i))).collect::<Vec<_>>()});
    }
    Ok(out)
}

/// The surfaces of a fluid domain's obstacles: bodies of the open design (triangulated at the
/// domain's tolerance) or imported surfaces, by name. The geometry is in millimetres; the
/// surfaces come out in metres, as the domain is.
fn obstacles(s: &Session, domain: &FluidDomain) -> Result<Vec<Obstacle>> {
    let tolerance = crate::units::to_geometry_length(domain.tolerance.unwrap_or(domain.spacing / 10.0));
    let mut out = Vec::with_capacity(domain.obstacles.len());
    for body in &domain.obstacles {
        let (mesh, faces) = if let Some(shape) = s.design.as_deref().and_then(|d| d.body(body)) {
            (std::borrow::Cow::Owned(shape.tessellate(tolerance)?), shape.faces().iter().map(|f| f.name.clone()).collect())
        } else if let Some(t) = s.surfaces.iter().find(|t| &t.name == body) {
            // Imported faces have no names: the solver calls them face0, face1, ….
            (std::borrow::Cow::Borrowed(&t.mesh), Vec::new())
        } else {
            return Err(EngineError::Params(format!("obstacle `{body}` is neither a body of the open design nor an imported surface")));
        };
        let positions: Vec<[f64; 3]> = mesh.positions.iter().map(|p| crate::units::geometry_point(*p)).collect();
        let options = SurfaceOptions { weld_distance: crate::units::geometry_length(mesh.weld_distance), ..SurfaceOptions::default() };
        let surface = BoundarySurface::new(&positions, &mesh.triangles, &mesh.tri_face, &options)
            .map_err(|e| EngineError::context(format!("obstacle `{body}`: its surface cannot bound a solid"), e))?;
        out.push(Obstacle { body: body.clone(), surface, faces });
    }
    Ok(out)
}

/// A cap of a fluid solve that may be raised: its default when left out; when given, a whole
/// number from 1 up to `ceiling`, what the session's memory bound (`bound`) holds.
fn raisable(p: &Value, key: &str, default: u64, ceiling: u64, bound: u64) -> Result<u64> {
    let Some(v) = p.get(key) else { return Ok(default) };
    match whole(v) {
        Some(n) if (1..=ceiling).contains(&n) => Ok(n),
        Some(n) if n > ceiling => Err(EngineError::invalid(
            key,
            if bound == 0 {
                format!("{n} is above the default {default}, and this machine's memory is unknown, so the cap cannot be raised")
            } else {
                format!(
                    "{n} is above the {ceiling} this machine allows (its bound for a run is {}, half its memory): coarsen the spacing or write fewer frames",
                    fieldcraft_lbm::memory::megabytes(bound)
                )
            },
        )),
        _ => Err(EngineError::param(key, ProblemKind::NotAWholeNumber { min: 1, max: i64::try_from(ceiling).unwrap_or(i64::MAX) })),
    }
}

/// The caps and backend of a `solve.fluid` call: `max_steps` may only be lowered;
/// `max_cells`, `max_cell_updates` and `max_frame_bytes` may be raised up to what the session's
/// memory bound holds.
fn fluid_limits(s: &Session, p: &Value) -> Result<(fieldcraft_lbm::Limits, BackendChoice)> {
    use fieldcraft_lbm::Limits as L;
    let bound = s.limits.memory;
    let ceilings = L { max_memory: bound, ..L::default() };
    let limits = L {
        max_cells: raisable(p, "max_cells", L::MAX_CELLS, ceilings.max_cells_ceiling(), bound)?,
        max_steps: cap(p, "max_steps", L::MAX_STEPS)?,
        max_cell_updates: raisable(p, "max_cell_updates", L::MAX_CELL_UPDATES, ceilings.max_cell_updates_ceiling(), bound)?,
        max_frame_bytes: raisable(p, "max_frame_bytes", L::MAX_FRAME_BYTES, ceilings.max_frame_ceiling(), bound)?,
        max_memory: bound,
    };
    Ok((limits, backend_choice(s, p, machine_backend().as_deref(), fieldcraft_lbm::device_available())?))
}

/// The machine's override of every step's preferred backend (`FIELDCRAFT_SOLVE_BACKEND`).
fn machine_backend() -> Option<String> {
    std::env::var("FIELDCRAFT_SOLVE_BACKEND").ok().filter(|v| !v.trim().is_empty())
}

/// The backend a fluid solve runs on, and why.
#[derive(Debug, Clone, PartialEq)]
pub struct BackendChoice {
    pub choice: fieldcraft_lbm::Choice,
    pub by: ChosenBy,
    /// Said when the step's preference could not be followed here.
    pub notice: Option<String>,
}

/// What chose a fluid solve's backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChosenBy {
    /// The call's own `backend`.
    Parameter,
    /// The machine's override, FIELDCRAFT_SOLVE_BACKEND.
    Machine,
    /// The step's saved preference (solve.settings).
    Step,
    Default,
}

impl ChosenBy {
    pub fn id(self) -> &'static str {
        match self {
            ChosenBy::Parameter => "parameter",
            ChosenBy::Machine => "machine",
            ChosenBy::Step => "step",
            ChosenBy::Default => "default",
        }
    }
}

impl BackendChoice {
    pub fn to_json(&self) -> Value {
        json!({"backend": choice_id(self.choice), "by": self.by.id(), "notice": self.notice})
    }
}

fn choice_id(c: fieldcraft_lbm::Choice) -> &'static str {
    match c {
        fieldcraft_lbm::Choice::Cpu => "cpu",
        fieldcraft_lbm::Choice::Device => "gpu",
    }
}

/// The backend of a `solve.fluid` call `p`: its own `backend`, else the machine's override
/// (`machine`), else the step's saved preference, falling back as it says when `device` says this
/// machine has no GPU to solve on, else the CPU.
pub fn backend_choice(s: &Session, p: &Value, machine: Option<&str>, device: std::result::Result<(), String>) -> Result<BackendChoice> {
    use fieldcraft_lbm::Choice;
    let parse = |b: &str, what: &str| {
        Choice::parse(b.trim())
            .ok_or_else(|| EngineError::Params(format!("{what} must be `cpu` or `gpu`, not `{}`", b.chars().take(40).collect::<String>())))
    };
    if let Some(b) = opt::<String>(p, "backend")? {
        return Ok(BackendChoice { choice: parse(&b, "`backend`")?, by: ChosenBy::Parameter, notice: None });
    }
    if let Some(m) = machine {
        return Ok(BackendChoice { choice: parse(m, "FIELDCRAFT_SOLVE_BACKEND")?, by: ChosenBy::Machine, notice: None });
    }
    let step = opt::<String>(p, "step")?;
    let name = fieldcraft_lbm::fluid_step(&s.model, step.as_deref()).ok().map(|(st, _)| st.name.clone());
    let Some((name, pref)) = name.and_then(|n| s.run.get(&n).copied().map(|r| (n, r))) else {
        return Ok(BackendChoice { choice: Choice::Cpu, by: ChosenBy::Default, notice: None });
    };
    match (pref.backend, device) {
        (fieldcraft_io::Backend::Cpu, _) => Ok(BackendChoice { choice: Choice::Cpu, by: ChosenBy::Step, notice: None }),
        (fieldcraft_io::Backend::Gpu, Ok(())) => Ok(BackendChoice { choice: Choice::Device, by: ChosenBy::Step, notice: None }),
        (fieldcraft_io::Backend::Gpu, Err(why)) => match pref.fallback {
            Some(fieldcraft_io::Backend::Cpu) => Ok(BackendChoice {
                choice: Choice::Cpu,
                by: ChosenBy::Step,
                notice: Some(format!(
                    "step `{name}` prefers the GPU, but {why}: it is solved on the CPU instead (slower; the same result as any CPU)"
                )),
            }),
            _ => Err(EngineError::Params(format!(
                "step `{name}` is to be solved on the GPU with no fallback, but {why} (solve.settings, or `backend: \"cpu\"`)"
            ))),
        },
    }
}

/// `solve.settings`: set or forget a fluid step's preferred backend, or list them all.
fn settings(s: &mut Session, p: &Value) -> Result<Value> {
    let listed = |s: &Session| json!({"run": s.run});
    let Some(step) = opt::<String>(p, "step")? else { return Ok(listed(s)) };
    let st =
        s.model.step(&step).ok_or_else(|| EngineError::Params(format!("step `{}` does not exist", step.chars().take(80).collect::<String>())))?;
    if !matches!(st.procedure, fieldcraft_model::Procedure::Fluid { .. }) {
        return Err(EngineError::invalid("step", format!("step `{step}` is a {} step: only fluid steps choose a backend", st.procedure.name())));
    }
    let backend: String = req(p, "backend")?;
    if backend == "none" {
        s.run.remove(&step);
        return Ok(listed(s));
    }
    let backend = fieldcraft_io::Backend::parse(&backend).ok_or_else(|| EngineError::invalid("backend", "must be gpu, cpu or none"))?;
    let fallback = match opt::<String>(p, "fallback")?.as_deref() {
        None if backend == fieldcraft_io::Backend::Gpu => Some(fieldcraft_io::Backend::Cpu),
        None | Some("none") => None,
        Some("cpu") if backend == fieldcraft_io::Backend::Gpu => Some(fieldcraft_io::Backend::Cpu),
        Some(_) => return Err(EngineError::invalid("fallback", "must be cpu (for a GPU preference) or none")),
    };
    s.run.insert(step, fieldcraft_io::RunSettings { backend, fallback });
    Ok(listed(s))
}

/// What a fluid solve would need in memory, from the model alone (no obstacle is
/// triangulated), against its caps and the machine's bound: what `solve.fluid` reports and the
/// Solve dialog shows before the run.
#[derive(Debug, Clone, PartialEq)]
pub struct FluidEstimate {
    pub memory: fieldcraft_lbm::Memory,
    pub limits: fieldcraft_lbm::Limits,
    /// The machine's physical memory, when known.
    pub physical: Option<u64>,
    /// Why the run would be refused before it starts (the cell, frame or memory caps), if it
    /// would.
    pub problem: Option<String>,
    /// The backend it would run on, and why (a step's preference, a fallback).
    pub backend: BackendChoice,
}

impl FluidEstimate {
    pub fn to_json(&self) -> Value {
        let m = &self.memory;
        let l = &self.limits;
        json!({
            "backend": match m.backend { fieldcraft_lbm::Choice::Cpu => "cpu", fieldcraft_lbm::Choice::Device => "gpu" },
            "cells": m.cells,
            "frames": m.frames,
            "flow_bytes": m.flow_bytes(),
            "frame_bytes": m.frame_bytes,
            "total_bytes": m.total(),
            "bound_bytes": l.max_memory,
            "physical_bytes": self.physical,
            "raised": l.raised(),
            "frames_that_fit": m.frames_that_fit(l),
            "ceilings": {"max_cells": l.max_cells_ceiling(), "max_cell_updates": l.max_cell_updates_ceiling(), "max_frame_bytes": l.max_frame_ceiling()},
            "gpu_buffer": m.device.map(|d| json!({"populations_bytes": d.populations, "limit_bytes": d.limit})),
            "problem": self.problem,
        })
    }

    /// The estimate in a line, for a form: what the run needs against what this machine allows.
    pub fn summary(&self) -> String {
        use fieldcraft_lbm::memory::megabytes;
        let m = &self.memory;
        let mut line = match (&self.backend.notice, self.backend.by) {
            (Some(n), _) => format!("{n}. "),
            (None, ChosenBy::Step) => format!("Backend: {}, the step's preference (solve.settings). ", choice_id(self.backend.choice).to_uppercase()),
            (None, ChosenBy::Machine) => {
                format!("Backend: {}, this machine's (FIELDCRAFT_SOLVE_BACKEND). ", choice_id(self.backend.choice).to_uppercase())
            }
            _ => String::new(),
        };
        line += &format!(
            "Memory: {} frames × {} cells: {} of frames + {} of flow ≈ {}",
            m.frames,
            crate::thousands(usize::try_from(m.cells).unwrap_or(usize::MAX)),
            megabytes(m.frame_bytes),
            megabytes(m.flow_bytes()),
            megabytes(m.total())
        );
        if let Some(d) = m.device {
            line += &format!("; GPU buffer {}", megabytes(d.populations));
            if let Some(limit) = d.limit {
                line += &format!(" of {} allowed", megabytes(limit));
            }
        }
        let bound = self.limits.max_memory;
        if bound == 0 {
            line += ". This machine's memory is unknown: the caps keep their defaults.";
        } else {
            line += &format!(". Runs with raised caps may take up to {} on this machine", megabytes(bound));
            if let Some(physical) = self.physical {
                line += &format!(" (it has {})", megabytes(physical));
            }
            line += ".";
        }
        line
    }
}

/// The memory estimate of `solve.fluid` with parameters `p`, or why there is none (no fluid
/// step, parameters out of range).
pub fn fluid_estimate(s: &Session, p: &Value) -> Result<FluidEstimate> {
    let (limits, chosen) = fluid_limits(s, p)?;
    let backend = chosen.choice;
    let step = opt::<String>(p, "step")?;
    let (_, domain) = fieldcraft_lbm::fluid_step(&s.model, step.as_deref())?;
    let memory = fieldcraft_lbm::Memory::estimate(&s.model, step.as_deref(), backend)?;
    let problem = limits.check_cells(domain.cells()?).and_then(|()| memory.check(&limits)).err().map(|e| e.to_string());
    Ok(FluidEstimate { memory, limits, physical: crate::machine::physical_memory(), problem, backend: chosen })
}

fn run_fluid(s: &mut Session, p: &Value) -> Result<Value> {
    let estimate = fluid_estimate(s, p)?;
    if opt::<bool>(p, "estimate")? == Some(true) {
        return Ok(json!({"memory": estimate.to_json()}));
    }
    let limits = estimate.limits;
    let backend = estimate.memory.backend;
    let chosen = estimate.backend.clone();
    let step = opt::<String>(p, "step")?;
    let (_, domain) = fieldcraft_lbm::fluid_step(&s.model, step.as_deref())?;
    let obstacles = obstacles(s, domain)?;
    let model = &s.model;
    let r = with_threads(threads(s, p)?, || fieldcraft_lbm::solve_fluid_on(model, step.as_deref(), &obstacles, &limits, backend))?;
    let mut summary = fluid_summary(&r, s.units());
    summary["memory"] = estimate.to_json();
    summary["backend_choice"] = chosen.to_json();
    if let Some(n) = &chosen.notice {
        summary["notice"] = json!(n);
    }
    s.store_result(StepResult::Fluid(r));
    Ok(summary)
}

/// The reply of a fluid solve, in SI but for the steady monitors' history, which is in `display`
/// (its quantities follow the monitors, so the reply table cannot name them).
fn fluid_summary(r: &FluidResult, display: &DisplaySystem) -> Value {
    let speed = |v: &[f32; 3]| (f64::from(v[0]).powi(2) + f64::from(v[1]).powi(2) + f64::from(v[2]).powi(2)).sqrt();
    let last = r.frames.last();
    // The largest speed of every frame: `max_speed` is the largest of all (the start-up often
    // passes the final state), `final_max_speed` the last frame's.
    // Over interior fluid cells: the open-boundary layers are excluded because a wake reaching a
    // uniform pressure outlet distorts the last layer (#590).
    let frame_max: Vec<f64> = r
        .frames
        .iter()
        .map(|f| f.velocity.iter().enumerate().filter(|(k, _)| !r.in_open_layer(*k)).map(|(_, v)| speed(v)).fold(0.0, f64::max))
        .collect();
    let (peak_frame, max_speed) = frame_max.iter().copied().enumerate().fold((0, 0.0), |best, (k, v)| if v > best.1 { (k, v) } else { best });
    // The mean over the last quarter of the history smooths the acoustic noise of the
    // weakly compressible scheme.
    let tail = r.history_times.len().div_ceil(4).max(1);
    let forces: Vec<Value> = r
        .forces
        .iter()
        .filter(|f| f.force.iter().any(|v| *v != [0.0; 3]))
        .map(|f| {
            let recent = f.force.get(f.force.len().saturating_sub(tail)..).unwrap_or_default();
            let n = recent.len().max(1) as f64;
            let mean = [0, 1, 2].map(|a| recent.iter().map(|v| v[a]).sum::<f64>() / n);
            json!({"name": f.name, "body": f.body, "final": f.force.last(), "mean": mean})
        })
        .collect();
    let (m0, m1) = (r.frames.first().map_or(0.0, |f| f.mass), last.map_or(0.0, |f| f.mass));
    let heat: Vec<Value> = r
        .heat
        .iter()
        .filter(|h| h.power.iter().any(|v| *v != 0.0))
        .map(|h| {
            let recent = h.power.get(h.power.len().saturating_sub(tail)..).unwrap_or_default();
            let mean = recent.iter().sum::<f64>() / recent.len().max(1) as f64;
            json!({"name": h.name, "body": h.body, "final": h.power.last(), "mean": mean})
        })
        .collect();
    let tail_mean = |v: &[f64]| {
        let recent = v.get(v.len().saturating_sub(tail)..).unwrap_or_default();
        recent.iter().sum::<f64>() / recent.len().max(1) as f64
    };
    let zones: Vec<Value> = r
        .zones
        .iter()
        .map(|z| {
            let mut v = json!({
                "name": z.name,
                "kind": if z.fan { "fan" } else { "porous" },
                "flow": {"final": z.flow.last(), "mean": tail_mean(&z.flow)},
            });
            if z.fan {
                v["pressure_rise"] = json!({"final": z.fan_rise.last(), "mean": tail_mean(&z.fan_rise)});
                v["measured_pressure_rise"] = json!({"final": z.pressure.last(), "mean": tail_mean(&z.pressure)});
            } else {
                v["pressure_drop"] = json!({"final": z.pressure.last(), "mean": tail_mean(&z.pressure)});
            }
            v
        })
        .collect();
    // Each conducting solid: its temperatures at the end, its power, the heat its faces give
    // the fluid (the sum over them of minus the heat into the wall), and its peak against its
    // limit.
    let checks = r.components();
    let solids: Vec<Value> = r
        .solids
        .iter()
        .zip(&checks)
        .map(|(sr, check)| {
            let field = last.map(|f| r.temperatures(f)).unwrap_or_default();
            let temps: Vec<f64> =
                if field.is_empty() { Vec::new() } else { sr.cells.iter().map(|c| field.get(*c as usize).copied().unwrap_or(f64::NAN)).collect() };
            let mut given = vec![0.0; r.history_times.len()];
            for h in r.heat.iter().filter(|h| h.body.as_deref() == Some(sr.body.as_str())) {
                for (g, p) in given.iter_mut().zip(&h.power) {
                    *g -= p;
                }
            }
            let mut v = json!({
                "body": sr.body,
                "cells": sr.cells.len(),
                "power": sr.power,
                "heat_to_fluid": {"final": given.last(), "mean": tail_mean(&given)},
            });
            if let Some(((_, lo), (_, hi))) = crate::results::range(&sr.cells, &temps) {
                v["temperature"] = json!({"min": lo, "max": hi, "mean": temps.iter().sum::<f64>() / temps.len() as f64});
            }
            if sr.capacity_scale != 1.0 {
                v["capacity_scale"] = json!(sr.capacity_scale);
            }
            if let Some(p) = check.peak {
                v["peak"] = json!(p.temperature);
                v["peak_time"] = json!(p.time);
            }
            v["limit"] = json!(check.limit);
            v["margin"] = json!(check.margin);
            v["ok"] = json!(check.ok);
            v
        })
        .collect();
    let temperature = last.filter(|f| !f.temperature.is_empty()).map(|f| {
        let fluid = r.temperatures(f).into_iter().zip(&r.solid).filter(|(_, solid)| !**solid).map(|(t, _)| t);
        let (lo, hi) = fluid.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| (lo.min(t), hi.max(t)));
        json!({"min": lo, "max": hi})
    });
    json!({
        "step": r.step,
        "backend": r.stats.backend,
        "cells": r.stats.cells,
        "fluid_cells": r.stats.fluid_cells,
        "steps": r.stats.steps,
        "time_step": r.stats.time_step,
        "end_time": r.stats.end_time,
        "tau": r.stats.tau_plus,
        "tau_minus": r.stats.tau_minus,
        "magic": r.stats.magic,
        "collision": r.stats.collision.name(),
        "eddy_viscosity_ratio": r.stats.eddy_viscosity_ratio,
        "velocity_scale": r.stats.velocity,
        "lattice_velocity": r.stats.lattice_velocity,
        "mach": r.stats.mach,
        "frames": r.frames.iter().map(|f| f.time).collect::<Vec<_>>(),
        "max_speed": max_speed,
        "max_speed_time": r.frames.get(peak_frame).map_or(0.0, |f| f.time),
        "final_max_speed": frame_max.last().copied().unwrap_or(0.0),
        "mass": {"first": m0, "last": m1, "relative_change": if m0 != 0.0 { (m1 - m0) / m0 } else { 0.0 }},
        "forces": forces,
        "tau_thermal": r.stats.tau_thermal,
        "temperature": temperature,
        "heat": heat,
        "solids": solids,
        "components_ok": r.components_ok(),
        "stop": r.steady.as_ref().map(|s| steady_summary(s, &r.stats, display)),
        "zones": zones,
        "warnings": r.warnings,
    })
}

/// The quantity of a steady-state monitor's values, by its name in the report (`velocity`,
/// `temperature`, `mass_flow x_max`, `force …`, `heat …`, `mean_temperature`): residuals are the
/// root mean square of the field (a temperature's of its change, so a difference).
fn monitor_quantity(name: &str) -> Option<Quantity> {
    let kind = name.split_whitespace().next().unwrap_or_default();
    Some(match kind {
        "velocity" => Quantity::Velocity,
        "temperature" => Quantity::TemperatureDifference,
        "mass_flow" => Quantity::MassFlow,
        "force" => Quantity::Force,
        "heat" => Quantity::Power,
        "mean_temperature" => Quantity::Temperature,
        _ => return None,
    })
}

/// Why a run to steady state stopped, when, and its monitor history, in `display`.
fn steady_summary(s: &fieldcraft_lbm::SteadyReport, stats: &fieldcraft_lbm::FluidStats, display: &DisplaySystem) -> Value {
    let quantities: Vec<Option<Quantity>> = s.quantities.iter().map(|q| monitor_quantity(q)).collect();
    let value = |k: usize, v: &Vec<f64>| {
        let shown: Vec<f64> = match quantities.get(k).copied().flatten() {
            Some(q) => v.iter().map(|x| crate::results::shown(q, *x, display)).collect(),
            None => v.clone(),
        };
        if let [x] = shown.as_slice() { json!(x) } else { json!(shown) }
    };
    let history: Vec<Value> = s
        .samples
        .iter()
        .map(|m| {
            // (The time converts with the reply, `units::REPLIES`.)
            json!({"step": m.step, "time": m.time,
                "values": m.values.iter().enumerate().map(|(k, v)| value(k, v)).collect::<Vec<_>>(), "changes": m.changes})
        })
        .collect();
    let units: serde_json::Map<String, Value> =
        s.quantities.iter().zip(&quantities).filter_map(|(n, q)| Some((n.clone(), json!(crate::results::unit_of((*q)?, display)?)))).collect();
    json!({
        "reason": s.reason.name(),
        "steps": stats.steps,
        "time": stats.end_time,
        "tolerance": s.tolerance,
        "window": s.window,
        "quantities": s.quantities,
        "message": s.message,
        "history": history,
        "units": units,
    })
}

fn probe_fluid(model: &fieldcraft_model::Model, r: &FluidResult, p: &Value) -> Result<Value> {
    if let Some(names) = opt::<Vec<String>>(p, "forces")? {
        if names.len() > MAX_PROBE_FORCES {
            return Err(EngineError::Params(format!("`forces` lists at most {MAX_PROBE_FORCES} walls")));
        }
        let out: Vec<Value> = names
            .iter()
            .map(|n| {
                let f =
                    r.forces.iter().find(|f| &f.name == n || f.body.as_ref().is_some_and(|b| format!("{b}/{}", f.name) == *n)).ok_or_else(|| {
                        EngineError::Params(format!("no wall `{n}` (box faces x_min … z_max, or obstacle faces by name or body/name)"))
                    })?;
                Ok(json!({"name": f.name, "body": f.body, "force": f.force}))
            })
            .collect::<Result<_>>()?;
        return Ok(json!({"step": r.step, "times": r.history_times, "forces": out}));
    }
    if let Some(names) = opt::<Vec<String>>(p, "heat")? {
        if names.len() > MAX_PROBE_FORCES {
            return Err(EngineError::Params(format!("`heat` lists at most {MAX_PROBE_FORCES} walls")));
        }
        if r.heat.is_empty() {
            return Err(EngineError::Params(format!("step `{}` carries no heat: give its domain `thermal` settings", r.step)));
        }
        let out: Vec<Value> = names
            .iter()
            .map(|n| {
                let h =
                    r.heat.iter().find(|h| &h.name == n || h.body.as_ref().is_some_and(|b| format!("{b}/{}", h.name) == *n)).ok_or_else(|| {
                        EngineError::Params(format!("no wall `{n}` (box faces x_min … z_max, or obstacle faces by name or body/name)"))
                    })?;
                Ok(json!({"name": h.name, "body": h.body, "power": h.power}))
            })
            .collect::<Result<_>>()?;
        return Ok(json!({"step": r.step, "times": r.history_times, "heat": out}));
    }
    let frame = match opt::<usize>(p, "frame")? {
        Some(k) if k < r.frames.len() => k,
        Some(_) => return Err(EngineError::Params(format!("`frame` must be from 0 to {}", r.frames.len().saturating_sub(1)))),
        None => r.frames.len().saturating_sub(1),
    };
    let time = r.frames.get(frame).map(|f| f.time);
    if let Some(line) = p.get("line") {
        return probe_line(model, r, frame, time, line);
    }
    let pts = points(p, "points", MAX_PROBE)?;
    let out: Vec<Value> = pts
        .iter()
        .map(|&pt| match r.probe(frame, pt) {
            Some((cell, velocity, pressure, temperature)) => {
                let mut v = json!({"point": pt, "cell": cell, "solid": r.solid.get(cell), "velocity": velocity, "pressure": pressure});
                if let (Some(t), Some(obj)) = (temperature, v.as_object_mut()) {
                    obj.insert("temperature".into(), json!(t));
                }
                v
            }
            None => json!({"point": pt, "cell": null}),
        })
        .collect();
    Ok(json!({"step": r.step, "frame": frame, "time": time, "points": out}))
}

/// A fluid step's field along a line (`flow::line`): `{"line": {"from", "to", "count", "field",
/// "face"}}`.
fn probe_line(model: &fieldcraft_model::Model, r: &FluidResult, frame: usize, time: Option<f64>, line: &Value) -> Result<Value> {
    use crate::flow::line::{self, LineField};
    let ends = |key: &str| -> Result<[f64; 3]> {
        line.get(key).and_then(super::xyz_value).ok_or_else(|| EngineError::invalid("line", format!("needs `{key}`: [x, y, z]")))
    };
    let (from, to) = (ends("from")?, ends("to")?);
    let count = match opt::<u64>(line, "count")? {
        // One point per cell spacing along the line, both ends included; a line too long for
        // that (or not a finite length) gets a count `line::points` refuses.
        None => {
            let n = (fieldcraft_render::Vec3::from(from).dist(fieldcraft_render::Vec3::from(to)) / r.grid.spacing).round() + 1.0;
            if n.is_finite() && n <= line::MAX_POINTS as f64 { n as usize } else { usize::MAX }
        }
        Some(n) => usize::try_from(n).unwrap_or(usize::MAX),
    };
    let points = line::points(from, to, count)?;
    let field: LineField = opt(line, "field")?.ok_or_else(|| EngineError::param("field", ProblemKind::Required))?;
    let profile = match field {
        LineField::WallHeatFlux => {
            let face = opt(line, "face")?.ok_or_else(|| EngineError::param("face", ProblemKind::Required))?;
            line::wall_heat_flux(r, frame, &line::Wall::of(model, r, face)?, points)?
        }
        _ => line::profile(r, frame, field, points)?,
    };
    let key = field.key();
    let extreme = |e: Option<line::Extreme>| e.map_or(Value::Null, |e| json!({key: e.value, "xyz": e.point}));
    Ok(json!({"step": r.step, "frame": frame, "time": time, "line": {
        "field": line.get("field"), "xyz": profile.points, key: profile.values, "max": extreme(profile.max), "min": extreme(profile.min),
    }}))
}

pub(crate) fn probe(s: &mut Session, p: &Value) -> Result<Value> {
    let step = opt::<String>(p, "step")?.ok_or_else(|| EngineError::param("step", ProblemKind::Required))?;
    let result = s.result(&step).ok_or_else(|| {
        EngineError::Params(if s.latest_result(&step).is_some() {
            format!("the result of step `{step}` is stale (the model changed since it was solved); solve it again")
        } else {
            format!("step `{step}` has no result; solve it first")
        })
    })?;
    if let StepResult::Fluid(r) = result.as_ref() {
        let mut out = probe_fluid(&s.model, r, p)?;
        let display = s.units().clone();
        show_probe(&mut out, result, &display);
        return Ok(out);
    }
    let (nodes, located) = match (p.get("nodes"), p.get("at")) {
        (Some(_), None) => (labels(p, "nodes", MAX_PROBE)?, None),
        (None, None) if p.get("set").is_some() => {
            let name = opt::<String>(p, "set")?.unwrap_or_default();
            let set = s.model.resolve_set(&name)?;
            if set.nodes.len() > MAX_PROBE {
                return Err(EngineError::Params(format!("set `{name}` has {} nodes; a probe takes at most {MAX_PROBE}", set.nodes.len())));
            }
            (set.nodes, None)
        }
        (None, Some(_)) => {
            let at = points(p, "at", MAX_PROBE)?;
            let found = nearest_nodes(&s.model, result.nodes(), &at)?;
            (found.iter().map(|(l, _, _)| *l).collect(), Some(found))
        }
        (None, None) if matches!(result.as_ref(), StepResult::Explicit(_)) || p.get("elements").is_some() => (Vec::new(), None),
        _ => return Err(EngineError::Params("give one of `nodes`, `at` or `set`".into())),
    };
    let mut out = match result.as_ref() {
        StepResult::Explicit(r) => {
            let elements = if p.get("elements").is_some() { labels(p, "elements", MAX_PROBE)? } else { Vec::new() };
            probe_frame(r, &step, &nodes, &elements, p)?
        }
        StepResult::Implicit(r) => probe_implicit(r, &step, &nodes, p)?,
        StepResult::Fluid(r) => probe_fluid(&s.model, r, p)?,
    };
    // The largest von Mises stress among the probed nodes (static results).
    let peak = out.get("nodes").and_then(Value::as_array).and_then(|entries| {
        entries
            .iter()
            .filter_map(|e| Some((e.get("node")?.as_u64()?, e.get("von_mises")?.as_f64()?)))
            .fold(None, |best: Option<(u64, f64)>, (n, v)| if best.is_none_or(|b| v > b.1) { Some((n, v)) } else { best })
    });
    if let (Some((node, value)), Some(obj)) = (peak, out.as_object_mut()) {
        obj.insert("max_von_mises".into(), json!({"node": node, "value": value}));
    }
    if let (Some(found), Some(entries)) = (located, out.get_mut("nodes").and_then(Value::as_array_mut)) {
        for (entry, (_, xyz, distance)) in entries.iter_mut().zip(found) {
            if let Some(obj) = entry.as_object_mut() {
                obj.insert("xyz".into(), json!(xyz));
                obj.insert("distance".into(), json!(distance));
            }
        }
    }
    let display = s.units().clone();
    show_probe(&mut out, result, &display);
    Ok(out)
}

/// What a value of a probe reply is, by its key: a quantity per entry of a vector (its index,
/// of `len`), or a scale (a mode shape's), or nothing (a label, a code, a strain).
fn probe_quantity(r: &StepResult, display: &fieldcraft_model::DisplaySystem, key: &str, len: usize, index: usize) -> Option<crate::units::Show> {
    use crate::units::Show;
    use fieldcraft_model::units::Quantity as Q;
    let heat = matches!(r.kind(), crate::results::Kind::Heat(_));
    let (translation, rotation) = crate::results::mode_scale(r, display);
    let q = match key {
        // A natural mode, mass-normalised: translations (and a reference point's rotations).
        "shape" => return Some(Show::Factor(if index < 3 { translation } else { rotation })),
        "rotation" if matches!(r.kind(), crate::results::Kind::Frequency(_)) => return Some(Show::Factor(rotation)),
        "displacement" if len == 6 && index >= 3 => Q::Angle,
        "displacement" | "xyz" | "distance" | "gap" | "slip" | "point" | "at" | "max_penetration" => Q::Length,
        "rotation" => Q::Angle,
        "reaction" if heat => Q::Power,
        "reaction" if len == 6 && index >= 3 => Q::Moment,
        "reaction" | "force" | "final" | "mean" if !heat => Q::Force,
        "reaction_moment" => Q::Moment,
        "stress" | "principal" | "von_mises" | "pressure" | "shear" | "nodal" | "bottom" | "middle" | "top" | "value" => Q::Pressure,
        // Shell section forces N_x, N_y, N_xy, M_x, M_y, M_xy, Q_x, Q_y, per unit length.
        "forces" if len == 8 => {
            if (3..6).contains(&index) {
                Q::Force
            } else {
                Q::ForcePerLength
            }
        }
        // Beam end forces: N, V2, V3, T, M2, M3.
        "ends" => {
            if index < 3 {
                Q::Force
            } else {
                Q::Moment
            }
        }
        "temperature" => Q::Temperature,
        "heat_flux" => Q::HeatFlux,
        "velocity" => Q::Velocity,
        "time" | "times" => Q::Time,
        "power" => Q::Power,
        _ => return None,
    };
    Some(Show::Quantity(q))
}

/// A probe reply (SI) as `display` shows it, with `units`: the unit of each key that has one.
fn show_probe(out: &mut Value, r: &StepResult, display: &fieldcraft_model::DisplaySystem) {
    let mut units = serde_json::Map::new();
    crate::units::show_by_key(out, display, &|key, len, index| probe_quantity(r, display, key, len, index), &mut units);
    if let Some(o) = out.as_object_mut() {
        o.insert("units".into(), Value::Object(units));
    }
}

/// The node with a result nearest to each point (lowest label on a tie), where it is and how far.
fn nearest_nodes(model: &fieldcraft_model::Model, with_result: &[u32], at: &[[f64; 3]]) -> Result<Vec<(u32, [f64; 3], f64)>> {
    let candidates: Vec<(u32, [f64; 3])> = with_result.iter().filter_map(|&l| Some((l, model.node(l)?.xyz))).collect();
    at.iter()
        .map(|p| {
            candidates
                .iter()
                .map(|&(l, x)| (l, x, ((x[0] - p[0]).powi(2) + (x[1] - p[1]).powi(2) + (x[2] - p[2]).powi(2)).sqrt()))
                .min_by(|a, b| a.2.total_cmp(&b.2))
                .ok_or_else(|| EngineError::Params("the result has no nodes to probe".into()))
        })
        .collect()
}

/// The result of an implicit step at nodes.
fn probe_implicit(result: &ImplicitResult, step: &str, nodes: &[u32], p: &Value) -> Result<Value> {
    let elements = if p.get("elements").is_some() { labels(p, "elements", MAX_PROBE)? } else { Vec::new() };
    let locate = |label: u32| {
        result
            .locate(label)
            .ok_or_else(|| EngineError::Params(format!("node {label} has no result (it belongs to no solved element or coupled reference point)")))
    };
    let out: Vec<Value> = match result {
        ImplicitResult::Static(r) => nodes
            .iter()
            .map(|&l| {
                Ok(match locate(l)? {
                    Location::Point(k) => {
                        let p = r.reference_points.get(k);
                        json!({"node": l, "reference_point": p.map(|p| &p.name), "displacement": p.map(|p| p.displacement), "reaction": p.map(|p| p.reaction)})
                    }
                    Location::Node(i) => {
                        let mut v = json!({
                            "node": l,
                            "displacement": r.displacement.get(i),
                            "rotation": r.rotation.get(i),
                            "reaction": r.reaction.get(i),
                            "reaction_moment": r.reaction_moment.get(i),
                            "stress": r.stress.get(i),
                            "principal": r.stress.get(i).map(fieldcraft_fem::principal),
                            "von_mises": r.von_mises.get(i),
                            "temperature": r.temperature.as_ref().and_then(|t| t.get(i)),
                        });
                        // The contacts the node is a secondary node of: pressure, gap, status.
                        let contact: Vec<Value> = r
                            .contacts
                            .iter()
                            .filter_map(|c| {
                                let k = c.nodes.binary_search(&l).ok()?;
                                let mut entry = json!({"contact": c.name, "pressure": c.pressure.get(k), "gap": c.gap.get(k), "status": c.status.get(k)});
                                if let (Some(shear), Some(slip)) = (c.shear.get(k), c.slip.get(k)) {
                                    entry["shear"] = json!(shear);
                                    entry["slip"] = json!(slip);
                                }
                                Some(entry)
                            })
                            .collect();
                        if !contact.is_empty() {
                            v["contact"] = Value::Array(contact);
                        }
                        v
                    }
                })
            })
            .collect::<Result<_>>()?,
        ImplicitResult::Frequency(r) => {
            if !elements.is_empty() {
                return Err(EngineError::Params("element results are kept for static and explicit steps".into()));
            }
            let mode = opt::<usize>(p, "mode")?.unwrap_or(1);
            let m = mode
                .checked_sub(1)
                .and_then(|k| r.modes.get(k))
                .ok_or_else(|| EngineError::Params(format!("`mode` must be from 1 to {}", r.modes.len())))?;
            nodes
                .iter()
                .map(|&l| {
                    Ok(match locate(l)? {
                        Location::Point(k) => {
                            json!({"node": l, "reference_point": r.reference_points.get(k).map(|p| &p.0), "shape": m.point_shape.get(k)})
                        }
                        Location::Node(i) => json!({"node": l, "shape": m.shape.get(i), "rotation": m.rotation.get(i)}),
                    })
                })
                .collect::<Result<_>>()?
        }
        ImplicitResult::Heat(r) => {
            let (time, t) = crate::results::heat_frame(r, crate::results::Pick { mode: 1, frame: opt::<usize>(p, "frame")? })?;
            let out = nodes
                .iter()
                .map(|&l| {
                    let i = r.nodes.binary_search(&l).map_err(|_| EngineError::Params(format!("node {l} has no temperature (it belongs to no solid element)")))?;
                    Ok(json!({"node": l, "temperature": t.get(i), "heat_flux": r.flux.get(i), "reaction": r.reaction.get(i)}))
                })
                .collect::<Result<Vec<Value>>>()?;
            return Ok(json!({"step": step, "time": time, "nodes": out}));
        }
    };
    let mut response = json!({"step": step, "nodes": out});
    if let ImplicitResult::Static(r) = result
        && !elements.is_empty()
    {
        let found: Vec<Value> = elements
            .iter()
            .map(|&l| element_result(r, l).ok_or_else(|| EngineError::Params(format!("element {l} has no result in step `{step}`"))))
            .collect::<Result<_>>()?;
        response["elements"] = Value::Array(found);
    }
    Ok(response)
}

/// The entry of a label-sorted result list.
fn by_label<T>(list: &[T], label: u32, key: impl Fn(&T) -> u32) -> Option<&T> {
    list.binary_search_by_key(&label, key).ok().and_then(|i| list.get(i))
}

/// The stored results of one element of a static step.
fn element_result(r: &StaticResult, label: u32) -> Option<Value> {
    if let Some(e) = by_label(&r.elements, label, |e| e.label) {
        let mut v = json!({"element": label, "kind": e.kind.name(), "stress": e.nodal});
        if !e.plastic_strain.is_empty() {
            v["plastic_strain"] = json!(e.plastic_strain);
        }
        return Some(v);
    }
    if let Some(e) = by_label(&r.shells, label, |e| e.label) {
        let s = &e.stress;
        return Some(
            json!({"element": label, "kind": e.kind.name(), "axes": s.axes, "forces": s.forces, "bottom": s.bottom, "middle": s.middle, "top": s.top}),
        );
    }
    let e = by_label(&r.beams, label, |e| e.label)?;
    Some(json!({"element": label, "kind": e.kind.name(), "axes": e.forces.axes, "ends": e.forces.ends}))
}

/// Displacement and velocity at nodes, and element results, in one field frame of an explicit
/// result.
fn probe_frame(r: &ExplicitResult, step: &str, nodes: &[u32], elements: &[u32], p: &Value) -> Result<Value> {
    let last = r.frames.len().saturating_sub(1);
    let k = opt::<usize>(p, "frame")?.unwrap_or(last);
    let f = r.frames.get(k).ok_or_else(|| EngineError::Params(format!("`frame` must be from 0 to {last}")))?;
    if nodes.is_empty() && elements.is_empty() {
        return Err(EngineError::Params("give `nodes`, `elements` or both".into()));
    }
    if !elements.is_empty() && f.stress.is_empty() {
        return Err(EngineError::Params(format!("step `{step}` does not output stresses")));
    }
    let element_values: Vec<Value> = elements
        .iter()
        .map(|&l| {
            let i =
                r.elements.binary_search(&l).map_err(|_| EngineError::Params(format!("element {l} has no result (it is not a solved element)")))?;
            Ok(json!({
                "element": l,
                "stress": f.stress.get(i),
                "von_mises": f.stress.get(i).map(fieldcraft_fem::von_mises),
                "strain": f.strain.get(i),
                "plastic_strain": f.plastic_strain.get(i),
                "damage": f.damage.get(i),
                "temperature": f.temperature.get(i),
                "eroded": f.failed.get(i),
            }))
        })
        .collect::<Result<_>>()?;
    if nodes.is_empty() {
        return Ok(json!({"step": step, "frame": k, "time": f.time, "elements": element_values}));
    }
    if f.displacement.is_empty() {
        return Err(EngineError::Params(format!("step `{step}` does not output displacements")));
    }
    let out: Vec<Value> = nodes
        .iter()
        .map(|&l| {
            let i =
                r.nodes.binary_search(&l).map_err(|_| EngineError::Params(format!("node {l} has no result (it belongs to no solved element)")))?;
            Ok(json!({"node": l, "displacement": f.displacement.get(i), "velocity": f.velocity.get(i)}))
        })
        .collect::<Result<_>>()?;
    Ok(json!({"step": step, "frame": k, "time": f.time, "nodes": out, "elements": element_values}))
}

#[cfg(test)]
pub(crate) mod tests {
    /// `file.sample` opens the cantilever with nothing to save; `solve.run` solves both its steps
    /// in order and the results view comes up on the first.
    #[test]
    fn samples_open_and_solve_run_solves_every_step() {
        use serde_json::json;
        let mut s = Session::default();
        s.execute("mesh.node.add", &json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let r = s.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        assert_eq!(r["sample"], "cantilever");
        assert!(!s.is_dirty() && s.path.is_none());
        assert!(s.model.node(1).is_some_and(|n| n.xyz != [0.0, 0.0, 0.0]) || s.model.nodes().len() > 1);
        assert!(s.execute("file.sample", &json!({"name": "nope"})).unwrap_err().to_string().contains("cantilever"));
        let r = s.execute("solve.run", &json!({})).unwrap();
        let steps: Vec<&str> = r["steps"].as_array().unwrap().iter().map(|x| x["step"].as_str().unwrap()).collect();
        assert_eq!(steps, ["Load", "Modes"]);
        assert!(s.result("Load").is_some() && s.result("Modes").is_some());
        assert_eq!(s.results_view.as_ref().map(|v| v.step.as_str()), Some("Load"));
        s.execute("file.sample", &json!({"name": "impact"})).unwrap();
        assert!(s.results_view.is_none(), "a new model shows no old results");
        s.execute("solve.run", &json!({})).unwrap();
        assert!(s.result("Impact").is_some());
        let mut empty = Session::default();
        assert!(empty.execute("solve.run", &json!({})).unwrap_err().to_string().contains("Needs a mesh"));
    }

    use super::*;
    use fieldcraft_model::{BoxFace, FluidBoundary, ThermalCondition};

    /// A cantilever 10 × 1 × 1 of ten hex8 elements with incompatible modes, clamped at x = 0,
    /// a total force of 1 down on the free end; a static and a frequency step.
    pub(crate) fn cantilever(s: &mut Session) {
        let mut script = Vec::new();
        let label = |i: usize, j: usize, k: usize| 1 + i + 11 * (j + 2 * k);
        for k in 0..2 {
            for j in 0..2 {
                for i in 0..11 {
                    script.push(json!({"command": "mesh.node.add", "params": {"label": label(i, j, k), "xyz": [i, j, k]}}));
                }
            }
        }
        for i in 0..10 {
            let nodes = [
                label(i, 0, 0),
                label(i + 1, 0, 0),
                label(i + 1, 1, 0),
                label(i, 1, 0),
                label(i, 0, 1),
                label(i + 1, 0, 1),
                label(i + 1, 1, 1),
                label(i, 1, 1),
            ];
            script.push(json!({"command": "mesh.element.add", "params": {"label": i + 1, "kind": "hex8", "nodes": nodes}}));
        }
        let root: Vec<usize> = (0..4).map(|n| label(0, n % 2, n / 2)).collect();
        let tip: Vec<usize> = (0..4).map(|n| label(10, n % 2, n / 2)).collect();
        script.extend([
            json!({"command": "model.set.create", "params": {"name": "all", "elements": (1..=10).collect::<Vec<_>>()}}),
            json!({"command": "model.set.create", "params": {"name": "root", "nodes": root}}),
            json!({"command": "model.set.create", "params": {"name": "tip", "nodes": tip}}),
            json!({"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0}, "density": 1}}),
            json!({"command": "model.section.create", "params": {"name": "bar", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {"incompatible_modes": true}}}),
            json!({"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}}),
            json!({"command": "model.load.create", "params": {"name": "push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -1]}}),
            json!({"command": "model.step.create", "params": {"name": "Bend", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["push"]}}),
            json!({"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 4}, "boundary_conditions": ["clamp"]}}),
        ]);
        s.run_script(&Value::Array(script)).unwrap();
    }

    #[test]
    fn solve_keep_and_probe_results() {
        let mut s = Session::default();
        cantilever(&mut s);
        let journal = s.journal.len();
        let r = s.execute("solve.static", &json!({})).unwrap();
        assert_eq!(r["step"], "Bend");
        // Euler–Bernoulli P L³ / (3 E I) = 4, plus shear; the clamp and one element through the
        // depth keep it within a few percent.
        let d = r["max_displacement"]["value"].as_f64().unwrap();
        assert!((d - 4.024).abs() / 4.024 < 0.03, "{r}");
        assert!((r["reaction_total"][2].as_f64().unwrap() - 1.0).abs() < 1e-9, "{r}");
        assert_eq!(s.journal.len(), journal, "solving does not change the model");
        let probe = s.execute("solve.probe", &json!({"step": "Bend", "nodes": [11, 1]})).unwrap();
        assert!(probe["nodes"][0]["displacement"][2].as_f64().unwrap() < -3.8, "{probe}");
        assert!(probe["nodes"][1]["reaction"][2].as_f64().unwrap() != 0.0, "{probe}");
        // The principal stresses of the nodal stress, largest first.
        let root = &probe["nodes"][1];
        // (Computed in SI and shown in MPa, as the stress: equal to round-off.)
        let stress: [f64; 6] = serde_json::from_value(root["stress"].clone()).unwrap();
        let principal: [f64; 3] = serde_json::from_value(root["principal"].clone()).unwrap();
        let scale = stress.iter().fold(0.0f64, |m, s| m.max(s.abs()));
        assert!(principal.iter().zip(fieldcraft_fem::principal(&stress)).all(|(a, b)| (a - b).abs() <= 1e-12 * scale), "{root}");
        // By location: the nearest node, with where it is.
        let at = s.execute("solve.probe", &json!({"step": "Bend", "at": [[10.1, -0.2, 0.1], [0, 0, 0]]})).unwrap();
        assert_eq!((at["nodes"][0]["node"].clone(), at["nodes"][1]["node"].clone()), (json!(11), json!(1)), "{at}");
        assert_eq!(at["nodes"][0]["displacement"], probe["nodes"][0]["displacement"]);
        assert_eq!(at["nodes"][0]["xyz"], json!([10.0, 0.0, 0.0]));
        assert!((at["nodes"][0]["distance"].as_f64().unwrap() - 0.06f64.sqrt()).abs() < 1e-12);
        for bad in [json!({"step": "Bend"}), json!({"step": "Bend", "nodes": [1], "at": [[0, 0, 0]]}), json!({"step": "Bend", "at": [[0, 0]]})] {
            assert!(s.execute("solve.probe", &bad).is_err(), "{bad}");
        }

        let f = s.execute("solve.frequency", &json!({"mass": "lumped"})).unwrap();
        assert_eq!(f["modes"].as_array().unwrap().len(), 4, "{f}");
        let shape = s.execute("solve.probe", &json!({"step": "Modes", "nodes": [44], "mode": 2})).unwrap();
        assert!(shape["nodes"][0]["shape"].is_array(), "{shape}");
        assert_eq!(s.solved_steps().count(), 2);

        // Results describe the model they were solved for: an edit clears them, a failed command
        // keeps them.
        assert!(s.execute("model.material.edit", &json!({"name": "Mat", "elastic": {"type": "isotropic", "young": -1}})).is_err());
        assert_eq!(s.solved_steps().count(), 2);
        s.execute("model.material.edit", &json!({"name": "Mat", "elastic": {"type": "isotropic", "young": 2000, "poisson": 0}, "density": 1}))
            .unwrap();
        assert_eq!(s.solved_steps().count(), 0);
        // Undo returns to the solved state, and its results are current again.
        s.execute("edit.undo", &json!({})).unwrap();
        assert_eq!(s.solved_steps().count(), 2);
        s.execute("edit.redo", &json!({})).unwrap();
        let e = s.execute("solve.probe", &json!({"step": "Bend", "nodes": [1]})).unwrap_err().to_string();
        assert!(e.contains("stale"), "{e}");
    }

    #[test]
    fn factor_cap_can_be_raised_deliberately_to_a_finite_ceiling() {
        assert_eq!(limits(&json!({})).unwrap(), Limits::default());
        for entries in [1, 400_000_000, 400_000_001, 467_000_000, 2_000_000_000] {
            assert_eq!(limits(&json!({"max_factor_entries": entries})).unwrap().max_factor_entries, entries);
        }
        for entries in
            [json!(null), json!(true), json!("467000000"), json!(-1), json!(0), json!(2.5), json!(2_000_000_001u64), json!(u64::MAX), json!(1e300)]
        {
            assert!(limits(&json!({"max_factor_entries": entries})).is_err(), "{entries}");
        }
        for params in [
            json!({"max_dofs": Limits::MAX_DOFS + 1}),
            json!({"max_iterations": Limits::MAX_LANCZOS_VECTORS + 1}),
            json!({"max_newton_iterations": Limits::MAX_NEWTON_ITERATIONS + 1}),
        ] {
            assert!(limits(&params).is_err(), "{params}");
        }
        for id in ["solve.static", "solve.nonlinear", "solve.heat", "solve.frequency", "solve.buckle"] {
            let schema = crate::schema::json_schema(crate::find_command(id).unwrap().form.unwrap());
            let factor = &schema["properties"]["max_factor_entries"];
            assert_eq!(factor["minimum"], 1, "{id}: {factor}");
            assert_eq!(factor["maximum"], 2_000_000_000u64, "{id}: {factor}");
        }
    }

    #[test]
    fn raised_factor_caps_solve_small_models_and_failures_keep_results() {
        let mut s = Session::default();
        cantilever(&mut s);
        s.run_script(&json!([
            {"command": "model.section.edit", "params": {"name": "bar", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {}}},
            {"command": "model.step.create", "params": {"name": "Nonlinear", "procedure": {"type": "nonlinear_static", "nonlinear_geometry": false}, "boundary_conditions": ["clamp"], "loads": ["push"]}},
            {"command": "model.load.create", "params": {"name": "compress", "type": "force", "target": {"set": "tip"}, "vector": [-1, 0, 0]}},
            {"command": "model.step.create", "params": {"name": "Buckle", "procedure": {"type": "buckle", "modes": 1}, "boundary_conditions": ["clamp"], "loads": ["compress"]}},
            {"command": "model.material.edit", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0}, "density": 1, "thermal": {"conductivity": 2}}},
            {"command": "model.boundary_condition.create", "params": {"name": "cold", "type": "temperature", "target": {"set": "root"}, "value": 0}},
            {"command": "model.boundary_condition.create", "params": {"name": "hot", "type": "temperature", "target": {"set": "tip"}, "value": 100}},
            {"command": "model.step.create", "params": {"name": "Heat", "procedure": {"type": "heat_transfer"}, "boundary_conditions": ["cold", "hot"]}},
        ])).unwrap();
        for (command, step) in [
            ("solve.static", "Bend"),
            ("solve.nonlinear", "Nonlinear"),
            ("solve.heat", "Heat"),
            ("solve.frequency", "Modes"),
            ("solve.buckle", "Buckle"),
        ] {
            let expected = s.execute(command, &json!({"step": step})).unwrap();
            let probe_params = json!({"step": step, "nodes": [6]});
            let probe = s.execute("solve.probe", &probe_params).unwrap();
            for entries in [400_000_001, 467_000_000, 2_000_000_000] {
                assert_eq!(s.execute(command, &json!({"step": step, "max_factor_entries": entries})).unwrap(), expected);
                assert_eq!(s.execute("solve.probe", &probe_params).unwrap(), probe);
            }
            let before = s.model.clone();
            let results = s.solved_steps().count();
            let journal = s.journal.len();
            for entries in [json!(1), json!(0), json!(-1), json!(2.5), json!(null), json!(2_000_000_001u64), json!(u64::MAX)] {
                let e = s.execute(command, &json!({"step": step, "max_factor_entries": entries})).unwrap_err().to_string();
                if entries == 1 {
                    assert!(e.contains("bytes of factor values alone") && e.contains("limit is 1 entries"), "{command}: {e}");
                }
                assert_eq!(s.model, before, "{command}: {entries}");
                assert_eq!(s.solved_steps().count(), results);
                assert_eq!(s.journal.len(), journal);
                assert_eq!(s.execute("solve.probe", &probe_params).unwrap(), probe);
            }
        }
    }

    #[test]
    fn solve_hostile_params_never_panic() {
        let mut s = Session::default();
        cantilever(&mut s);
        s.execute("solve.static", &json!({})).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        let hostile = [
            json!(null),
            json!({"step": 7}),
            json!({"step": "Nope"}),
            json!({"step": "Modes"}),
            json!({"step": "Bend", "max_dofs": 0}),
            json!({"step": "Bend", "max_dofs": 1e300}),
            json!({"step": "Bend", "max_dofs": 3}),
            json!({"step": "Bend", "max_factor_entries": -1}),
            json!({"step": "Bend", "max_factor_entries": 10}),
            json!({"step": "Modes", "max_iterations": 1}),
            json!({"step": "Modes", "max_iterations": 7}),
            json!({"step": "Modes", "mass": "heavy"}),
            json!({"step": "Modes", "shift": 1e308}),
            json!({"step": "Modes", "shift": -1e308}),
            json!({"step": "Modes", "shift": 5000.0}),
            json!({"step": "Modes", "shift": "x"}),
            json!({"step": "Modes", "nodes": [1], "mode": 0}),
            json!({"step": "Modes", "nodes": [1], "mode": 1e300}),
            json!({"step": "Bend", "nodes": [0, 99999999, 4294967296u64]}),
            json!({"step": "Bend", "nodes": vec![1; 20_000]}),
            json!({"step": "Bend", "threads": 0}),
            json!({"step": "Bend", "threads": -3}),
            json!({"step": "Bend", "threads": 257}),
            json!({"step": "Bend", "threads": 1e300}),
            json!({"step": "Bend", "threads": "all"}),
            json!({"step": "Modes", "threads": 18446744073709551615u64}),
            json!({"step": "Modes", "threads": 2.5}),
        ];
        for id in ["solve.static", "solve.frequency", "solve.probe"] {
            for p in &hostile {
                let mut t = s.clone();
                let before = t.model.clone();
                let results = t.solved_steps().count();
                if t.execute(id, p).is_err() {
                    assert_eq!(t.model, before);
                    assert_eq!(t.solved_steps().count(), results, "{id} {p}");
                }
            }
        }
    }

    /// #648: a fluid step's preferred backend is saved with the model, outside its definition;
    /// solving without a backend of its own follows it, falls back with a notice on a machine with
    /// no suitable GPU (this test binary registers none), and the machine overrides it.
    #[test]
    fn a_step_prefers_a_backend_and_falls_back_with_a_notice() {
        let mut s = Session::default();
        couette(&mut s);
        let (journal, revision) = (s.journal.len(), s.revision());
        let r = s.execute("solve.settings", &json!({"step": "Shear", "backend": "gpu"})).unwrap();
        assert_eq!(r["run"], json!({"Shear": {"backend": "gpu", "fallback": "cpu"}}));
        assert_eq!((s.journal.len(), s.revision()), (journal, revision), "not an edit");
        assert!(s.is_dirty(), "saved with the model");
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        assert_eq!(r["backend_choice"]["backend"], "cpu", "{}", r["backend_choice"]);
        assert_eq!(r["backend_choice"]["by"], "step");
        assert!(
            r["notice"].as_str().unwrap().contains("step `Shear` prefers the GPU, but this build has no GPU backend: it is solved on the CPU"),
            "{r}"
        );
        let fell_back = format!("{:?}", s.result("Shear"));
        s.execute("solve.fluid", &json!({"backend": "cpu"})).unwrap();
        assert_eq!(format!("{:?}", s.result("Shear")), fell_back, "the CPU's result");
        // The Solve dialog says so before the run.
        let note = fluid_estimate(&s, &json!({})).unwrap().summary();
        assert!(note.starts_with("step `Shear` prefers the GPU"), "{note}");
        // With a GPU here, or by the machine's word, or the call's own.
        let pick = |s: &Session, p: Value, machine: Option<&str>, gpu: bool| {
            backend_choice(s, &p, machine, if gpu { Ok(()) } else { Err("no GPU".into()) })
                .map(|c| (choice_id(c.choice), c.by.id(), c.notice.is_some()))
        };
        assert_eq!(pick(&s, json!({}), None, true).unwrap(), ("gpu", "step", false));
        assert_eq!(pick(&s, json!({}), Some("cpu"), true).unwrap(), ("cpu", "machine", false));
        assert_eq!(pick(&s, json!({"backend": "gpu"}), Some("cpu"), false).unwrap(), ("gpu", "parameter", false));
        assert!(pick(&s, json!({}), Some("tpu"), true).is_err());
        // No fallback: refused where there is no GPU.
        s.execute("solve.settings", &json!({"step": "Shear", "backend": "gpu", "fallback": "none"})).unwrap();
        assert!(pick(&s, json!({}), None, false).unwrap_err().to_string().contains("with no fallback"));
        assert!(s.execute("solve.fluid", &json!({})).unwrap_err().to_string().contains("with no fallback"));
        // Saved and opened again; the definition and its results' hash are unchanged.
        let d = std::env::temp_dir().join(format!("fieldcraft-run-settings-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let path = d.join("gap.fieldcraft");
        s.execute("solve.settings", &json!({"step": "Shear", "backend": "gpu"})).unwrap();
        s.execute("file.save_as", &json!({"path": path})).unwrap();
        let mut t = Session::default();
        let r = t.execute("file.open", &json!({"path": path})).unwrap();
        assert_eq!(r["results"], json!([{"step": "Shear", "stale": false}]), "outside the definition: {r}");
        assert_eq!(t.run, s.run);
        assert_eq!(t.execute("solve.settings", &json!({})).unwrap()["run"]["Shear"]["backend"], "gpu");
        t.execute("solve.settings", &json!({"step": "Shear", "backend": "none"})).unwrap();
        assert!(t.run.is_empty() && t.is_dirty());
        let _ = std::fs::remove_dir_all(&d);
        // Only fluid steps choose a backend; hostile parameters change nothing.
        let mut c = Session::default();
        cantilever(&mut c);
        for p in [
            json!({"step": "Bend", "backend": "gpu"}),
            json!({"step": "Nope", "backend": "gpu"}),
            json!({"step": "Shear"}),
            json!({"step": 7, "backend": "gpu"}),
        ] {
            assert!(c.execute("solve.settings", &p).is_err(), "{p}");
        }
        for p in [
            json!({"step": "Shear", "backend": "tpu"}),
            json!({"step": "Shear", "backend": "cpu", "fallback": "cpu"}),
            json!({"step": "Shear", "backend": "gpu", "fallback": "gpu"}),
        ] {
            assert!(t.execute("solve.settings", &p).is_err(), "{p}");
        }
        assert!(t.run.is_empty());
    }

    /// Couette flow between y = 0 and y = 1 (eight cells), the top wall sliding at 0.1 along x,
    /// periodic in x and z: material, domain and fluid step through commands.
    pub(crate) fn couette(s: &mut Session) {
        let periodic = json!({"type": "periodic"});
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Gap", "material": "Oil", "min": [0, 0, 0], "max": [0.125, 1, 0.125], "spacing": 0.125,
                "faces": {"x_min": periodic, "x_max": periodic, "z_min": periodic, "z_max": periodic}}},
            {"command": "model.fluid_domain.boundary", "params": {"domain": "Gap", "at": "y_max", "condition": {"type": "moving_wall", "velocity": [0.1, 0, 0]}}},
            {"command": "model.step.create", "params": {"name": "Shear", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "frames": 4}}},
        ]);
        s.run_script(&script).unwrap();
    }

    /// #590: `max_speed` is the largest speed of all frames (the start-up may pass the final
    /// state), with its time; `final_max_speed` the last frame's.
    #[test]
    fn the_fluid_summary_reports_the_largest_speed_of_all_frames() {
        let mut r = crate::flow::tests::channel();
        let mut surge = r.frames[1].clone();
        surge.time = 0.5;
        surge.velocity[0] = [30.0, 40.0, 0.0];
        r.frames.insert(1, surge);
        let summary = fluid_summary(&r, &fieldcraft_model::DisplaySystem::default());
        assert_eq!((summary["max_speed"].as_f64(), summary["max_speed_time"].as_f64()), (Some(50.0), Some(0.5)), "{summary}");
        // The last frame's fastest cell: x = 9 of ten cells, speed 10.
        assert_eq!(summary["final_max_speed"].as_f64(), Some(10.0), "{summary}");
        // With an inlet and an outlet on the x faces, their layers (x = 0 and 9) are left out.
        r.open_faces[0] = true;
        r.open_faces[1] = true;
        assert_eq!(fluid_summary(&r, &fieldcraft_model::DisplaySystem::default())["final_max_speed"].as_f64(), Some(9.0));
        let values = crate::flow::cell_values(&r, crate::results::Field::Velocity, r.frames.last().unwrap()).unwrap();
        assert_eq!(crate::flow::interior_range(&r, &values).map(|(lo, hi)| (lo.1, hi.1)), Some((2.0, 9.0)));
    }

    #[test]
    fn solve_a_fluid_step_and_probe_its_fields_and_wall_forces() {
        let mut s = Session::default();
        couette(&mut s);
        let journal = s.journal.len();
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        assert_eq!(r["step"], "Shear");
        assert_eq!(r["frames"].as_array().unwrap().len(), 5, "{r}");
        assert!((r["tau"].as_f64().unwrap() - 1.0).abs() < 1e-3, "the automatic time step keeps τ ≤ 1: {r}");
        assert!(r["mass"]["relative_change"].as_f64().unwrap().abs() < 1e-12, "{r}");
        assert_eq!(s.journal.len(), journal, "solving does not change the model");
        // The wall shear μ U / H over the wall's area 0.125².
        let bottom = r["forces"].as_array().unwrap().iter().find(|f| f["name"] == "y_min").unwrap();
        let expected = 0.1 * 0.1 / 1.0 * 0.125 * 0.125;
        assert!((bottom["mean"][0].as_f64().unwrap() / expected - 1.0).abs() < 1e-3, "{bottom}");
        assert_eq!(r["collision"], "trt", "{r}");
        // The regularised collision gives the same linear flow and wall shear, and says so.
        s.execute(
            "model.step.edit",
            &json!({"name": "Shear", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "frames": 4, "collision": "regularized"}}),
        )
        .unwrap();
        let rr = s.execute("solve.fluid", &json!({})).unwrap();
        assert_eq!(rr["collision"], "regularized", "{rr}");
        let bottom = rr["forces"].as_array().unwrap().iter().find(|f| f["name"] == "y_min").unwrap();
        assert!((bottom["mean"][0].as_f64().unwrap() / expected - 1.0).abs() < 1e-3, "{bottom}");
        assert!(
            s.execute(
                "model.step.edit",
                &json!({"name": "Shear", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "collision": "bgk"}})
            )
            .is_err()
        );

        let probe = s.execute("solve.probe", &json!({"step": "Shear", "points": [[0.06, 0.5625, 0.06], [5, 5, 5]]})).unwrap();
        let u = probe["points"][0]["velocity"][0].as_f64().unwrap();
        assert!((u - 0.05625).abs() < 1e-5, "{probe}");
        assert!(probe["points"][1]["cell"].is_null());
        let first = s.execute("solve.probe", &json!({"step": "Shear", "points": [[0.06, 0.5625, 0.06]], "frame": 0})).unwrap();
        assert_eq!(first["points"][0]["velocity"][0], 0.0, "the initial frame is at rest");
        let h = s.execute("solve.probe", &json!({"step": "Shear", "forces": ["y_min", "y_max"]})).unwrap();
        assert_eq!(h["forces"].as_array().unwrap().len(), 2);
        assert_eq!(h["times"].as_array().unwrap().len(), h["forces"][0]["force"].as_array().unwrap().len());

        // An edit makes the result stale; the step solves again.
        s.execute(
            "model.fluid_domain.boundary",
            &json!({"domain": "Gap", "at": "y_max", "condition": {"type": "moving_wall", "velocity": [0.2, 0, 0]}}),
        )
        .unwrap();
        assert!(s.execute("solve.probe", &json!({"step": "Shear", "points": [[0, 0, 0]]})).is_err());
        // The domain cannot go while a step uses it; the material cannot lose its fluid behaviour.
        assert!(s.execute("model.fluid_domain.delete", &json!({"name": "Gap"})).is_err());
        assert!(s.execute("model.material.edit", &json!({"name": "Oil", "density": 1.0})).is_err());
        // A structural step is not a fluid one, and a fluid step is not written to decks.
        assert!(s.execute("solve.static", &json!({"step": "Shear"})).is_err());
        for format in ["abaqus", "nastran", "lsdyna"] {
            let e = s.execute("deck.write", &json!({"format": format, "path": "/no/such/dir/x.deck"})).unwrap_err().to_string();
            assert!(e.contains("fluid step"), "{format}: {e}");
        }
    }

    /// A pressure outlet made convective through the boundary command solves, and the model keeps
    /// the choice; a uniform outlet does not write it.
    #[test]
    fn a_convective_outlet_is_set_by_command_and_solves() {
        let mut s = Session::default();
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Duct", "material": "Oil", "min": [0, 0, 0], "max": [2, 1, 0.125], "spacing": 0.125,
                "faces": {"x_min": {"type": "velocity", "velocity": [0.1, 0, 0], "profile": "parabolic"}, "x_max": {"type": "pressure", "pressure": 0},
                    "z_min": {"type": "periodic"}, "z_max": {"type": "periodic"}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Duct", "steps": 200, "frames": 2}}},
        ]);
        s.run_script(&script).unwrap();
        let faces = |s: &Session| serde_json::to_value(s.model.fluid_domain("Duct").unwrap().faces).unwrap();
        assert_eq!(faces(&s)["x_max"], json!({"type": "pressure", "pressure": 0.0}));
        s.execute(
            "model.fluid_domain.boundary",
            &json!({"domain": "Duct", "at": "x_max", "condition": {"type": "pressure", "pressure": 0, "outflow": "convective"}}),
        )
        .unwrap();
        assert_eq!(faces(&s)["x_max"], json!({"type": "pressure", "pressure": 0.0, "outflow": "convective"}));
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        assert!(r["max_speed"].as_f64().is_some_and(|v| v > 0.0 && v.is_finite()), "{r}");
        let e = s.execute(
            "model.fluid_domain.boundary",
            &json!({"domain": "Duct", "at": "x_max", "condition": {"type": "pressure", "pressure": 0, "outflow": "sponge"}}),
        );
        assert!(e.is_err());
    }

    /// The Couette gap with heat: the sliding top wall held at 1, the bottom at 0, the fluid's
    /// diffusivity k / (ρ c) = 0.1.
    pub(crate) fn heated_couette(s: &mut Session) {
        couette(s);
        heat_the_couette(s);
    }

    fn heat_the_couette(s: &mut Session) {
        s.execute(
            "model.material.edit",
            &json!({"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1},
            "thermal": {"conductivity": 0.1, "specific_heat": 1.0}}),
        )
        .unwrap();
        s.execute(
            "model.fluid_domain.edit",
            &json!({"name": "Gap", "material": "Oil", "min": [0, 0, 0], "max": [0.125, 1, 0.125], "spacing": 0.125,
            "faces": {"x_min": {"type": "periodic"}, "x_max": {"type": "periodic"}, "z_min": {"type": "periodic"}, "z_max": {"type": "periodic"},
                "y_max": {"type": "moving_wall", "velocity": [0.1, 0, 0]}},
            "thermal": {"initial_temperature": 0.5, "faces": {"y_min": {"type": "temperature", "value": 0}}}}),
        )
        .unwrap();
        s.execute("model.fluid_domain.thermal_boundary", &json!({"domain": "Gap", "at": "y_max", "condition": {"type": "temperature", "value": 1}}))
            .unwrap();
    }

    #[test]
    fn a_fluid_step_carries_heat_between_its_walls() {
        let mut s = Session::default();
        heated_couette(&mut s);
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        assert!(r["tau_thermal"].as_f64().unwrap() > 0.5, "{r}");
        assert!((r["temperature"]["min"].as_f64().unwrap() - 0.0625).abs() < 1e-6, "{r}");
        assert!((r["temperature"]["max"].as_f64().unwrap() - 0.9375).abs() < 1e-6, "{r}");
        // Conduction across the gap: k ΔT / H over the wall's area, out of the hot wall and
        // into the cold one.
        let expected = 0.1 * 1.0 / 1.0 * 0.125 * 0.125;
        let heat = r["heat"].as_array().unwrap();
        let cold = heat.iter().find(|h| h["name"] == "y_min").unwrap();
        let hot = heat.iter().find(|h| h["name"] == "y_max").unwrap();
        assert!((cold["mean"].as_f64().unwrap() / expected - 1.0).abs() < 1e-6, "{cold}");
        assert!((hot["mean"].as_f64().unwrap() / expected + 1.0).abs() < 1e-6, "{hot}");

        let probe = s.execute("solve.probe", &json!({"step": "Shear", "points": [[0.06, 0.5625, 0.06]]})).unwrap();
        assert!((probe["points"][0]["temperature"].as_f64().unwrap() - 0.5625).abs() < 1e-6, "{probe}");
        let h = s.execute("solve.probe", &json!({"step": "Shear", "heat": ["y_min", "y_max"]})).unwrap();
        assert_eq!(h["heat"].as_array().unwrap().len(), 2);
        assert_eq!(h["times"].as_array().unwrap().len(), h["heat"][0]["power"].as_array().unwrap().len());
        assert!(s.execute("solve.probe", &json!({"step": "Shear", "heat": ["nope"]})).is_err());
        // Along a line: the linear temperature and shear profiles, and the heat into each wall
        // from the temperature beside it, which for a linear profile is the wall's heat per area.
        let line = |field: &str, face: Option<&str>, from: [f64; 3], to: [f64; 3], count: u32| {
            let mut p = json!({"step": "Shear", "line": {"from": from, "to": to, "count": count, "field": field}});
            if let Some(f) = face {
                p["line"]["face"] = json!(f);
            }
            s.clone().execute("solve.probe", &p)
        };
        let t = line("temperature", None, [0.0625, 0.0625, 0.0625], [0.0625, 0.9375, 0.0625], 8).unwrap();
        assert!((t["line"]["max"]["temperature"].as_f64().unwrap() - 0.9375).abs() < 1e-6, "{t}");
        assert_eq!(t["line"]["max"]["xyz"], json!([0.0625, 0.9375, 0.0625]), "{t}");
        assert_eq!(t["line"]["temperature"].as_array().unwrap().len(), 8);
        let u = line("velocity_x", None, [0.0625, 0.0625, 0.0625], [0.0625, 0.9375, 0.0625], 8).unwrap();
        assert!((u["line"]["min"]["velocity"].as_f64().unwrap() - 0.1 * 0.0625).abs() < 1e-4, "{u}");
        let area = 0.125 * 0.125;
        for (face, y) in [("y_min", 0.0), ("y_max", 1.0)] {
            let q = line("wall_heat_flux", Some(face), [0.0625, y, 0.0625], [0.0625, y, 0.0625], 1).unwrap();
            let wall = heat.iter().find(|h| h["name"] == face).unwrap();
            let expected = wall["mean"].as_f64().unwrap() / area;
            let got = q["line"]["heat_flux"][0].as_f64().unwrap();
            assert!((got / expected - 1.0).abs() < 1e-5, "{face}: {got} vs {expected}");
        }
        assert!(line("wall_heat_flux", Some("x_min"), [0.0, 0.5, 0.0625], [0.0, 0.5, 0.0625], 1).is_err(), "a periodic face is no wall");
        assert!(line("wall_heat_flux", Some("y_min"), [0.0625, 0.1, 0.0625], [0.0625, 0.1, 0.0625], 1).is_err(), "off the wall");
        assert!(line("wall_heat_flux", None, [0.0625, 0.0, 0.0625], [0.0625, 0.0, 0.0625], 1).is_err(), "no face");
        assert!(line("vorticity", None, [0.0625, 0.0, 0.0625], [0.0625, 1.0, 0.0625], 2).is_err());
        assert!(line("speed", None, [0.0625, 0.0, 0.0625], [0.0625, 1.0, 0.0625], 0).is_err());
        // Hostile ends without a count: refused with the points' range, never an overflow.
        for to in [[1e300, 1e300, 1e300], [f64::MAX, 0.0, 0.0]] {
            let e = s.clone().execute("solve.probe", &json!({"step": "Shear", "line": {"from": [0, 0, 0], "to": to, "field": "speed"}}));
            assert!(e.as_ref().is_err_and(|e| e.to_string().contains("1 to 10000 points")), "{e:?}");
        }
        // A steady monitor's history is in the display units: the mean of the linear profile
        // from 0 to 1 °C is ½ °C (not 273.65 K), its heat in mW.
        let steady = json!({"tolerance": 1e-3, "monitors": [{"type": "mean_temperature"}, {"type": "heat", "faces": ["y_min"]}]});
        let mut t = s.clone();
        t.execute("model.step.create", &json!({"name": "Settle", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "steady": steady}}))
            .unwrap();
        let stop = t.execute("solve.fluid", &json!({"step": "Settle"})).unwrap()["stop"].clone();
        let last = stop["history"].as_array().unwrap().last().unwrap().clone();
        let mean = last["values"][2].as_f64().unwrap();
        assert!((mean - 0.5).abs() < 1e-3, "{stop}");
        assert!((last["values"][3].as_f64().unwrap() / expected - 1.0).abs() < 1e-2, "{stop}");
        assert_eq!((stop["units"]["mean_temperature"].as_str(), stop["units"]["heat y_min"].as_str()), (Some("°C"), Some("mW")), "{stop}");

        // The thermal boundary needs the domain's heat settings; the material must keep its
        // heat properties while a domain carries heat.
        let isothermal = json!({"name": "Gap", "material": "Oil", "min": [0, 0, 0], "max": [0.125, 1, 0.125], "spacing": 0.125});
        s.execute("model.fluid_domain.edit", &isothermal).unwrap();
        assert!(
            s.execute("model.fluid_domain.thermal_boundary", &json!({"domain": "Gap", "at": "y_max", "condition": {"type": "adiabatic"}})).is_err()
        );
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        assert!(r["temperature"].is_null() && r["heat"].as_array().unwrap().is_empty(), "{r}");
        assert!(s.execute("solve.probe", &json!({"step": "Shear", "heat": ["y_min"]})).is_err());
        heat_the_couette(&mut s);
        assert!(s.execute("model.material.edit", &json!({"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}})).is_err());
        // JSON has no NaN (`json!(f64::NAN)` is `null`): the model refuses one directly.
        let nan = ThermalCondition::Temperature { value: f64::NAN };
        assert!(s.model.set_fluid_thermal_boundary("Gap", FluidBoundary::Box(BoxFace::YMax), nan).is_err());
        let hostile = [
            json!({"domain": "Gap", "at": "y_max", "condition": {"type": "temperature", "value": "nan"}}),
            json!({"domain": "Gap", "at": "top", "condition": {"type": "temperature", "value": 1}}),
            json!({"domain": "Gap", "at": {"body": "Nope", "face": "F1"}, "condition": {"type": "temperature", "value": 1}}),
            json!({"domain": "Nope", "at": "y_max", "condition": {"type": "temperature", "value": 1}}),
            json!({"domain": "Gap", "at": "y_max", "condition": {"type": "oven"}}),
            json!({"domain": "Gap", "at": "y_max"}),
            json!(null),
        ];
        for p in &hostile {
            let before = s.model.clone();
            assert!(s.execute("model.fluid_domain.thermal_boundary", p).is_err(), "{p}");
            assert_eq!(s.model, before, "{p}");
        }
        // A finite flux the lattice cannot hold, or one on a periodic face, is refused by the
        // solve, not by the model.
        // (1e300 mW/mm² is 1e303 W/m²: finite in SI.)
        for (at, value) in [("y_max", 1e300), ("x_min", 1.0)] {
            let mut t = s.clone();
            t.execute("model.fluid_domain.thermal_boundary", &json!({"domain": "Gap", "at": at, "condition": {"type": "flux", "value": value}}))
                .unwrap();
            let r = t.execute("solve.fluid", &json!({}));
            assert!(r.is_err(), "{at}: {r:?}");
        }
    }

    #[test]
    fn a_fluid_domain_takes_its_obstacles_from_the_open_design() {
        use fieldcraft_geometry::sample::{Plate, plate_design};
        let dir = std::env::temp_dir().join(format!("fieldcraft-fluid-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("plate.solidcraft");
        std::fs::write(&path, plate_design(&Plate::default()).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path.to_str().unwrap()})).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        let bodies = s.execute("geometry.bodies", &json!({})).unwrap();
        let [lo, hi] = [0, 1].map(|k| bodies[0]["bbox"][k].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect::<Vec<_>>());
        // A box 20 mm larger than the plate on every side, 5 mm cells: a viscous liquid (ν = 100 mm²/s) flowing along x.
        let min: Vec<f64> = lo.iter().map(|v| (v / 5.0).floor() * 5.0 - 20.0).collect();
        let max: Vec<f64> = hi.iter().map(|v| (v / 5.0).ceil() * 5.0 + 20.0).collect();
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-7}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [10, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "steps": 60, "frames": 2}}},
        ]);
        s.run_script(&script).unwrap();
        let r = s.execute("solve.fluid", &json!({"step": "Flow"})).unwrap();
        assert!(r["fluid_cells"].as_u64().unwrap() < r["cells"].as_u64().unwrap(), "{r}");
        let on_plate: Vec<&Value> = r["forces"].as_array().unwrap().iter().filter(|f| f["body"] == "Body1").collect();
        assert!(!on_plate.is_empty(), "{r}");
        let faces = bodies[0]["faces"].as_array().unwrap();
        for f in &on_plate {
            assert!(faces.iter().any(|g| g["name"] == f["name"]), "forces are reported per named face: {f}");
        }
        let drag: f64 = on_plate.iter().map(|f| f["final"][0].as_f64().unwrap()).sum();
        assert!(drag > 0.0, "the flow pushes the plate downstream: {r}");
        // An obstacle that is not in the design.
        s.execute(
            "model.fluid_domain.edit",
            &json!({"name": "Tank", "material": "Water", "min": min, "max": max, "spacing": 5, "obstacles": ["Nope"]}),
        )
        .unwrap();
        assert!(s.execute("solve.fluid", &json!({})).is_err());
    }

    /// A duct 60 × 4 × 4 with walls around, periodic along x (a closed loop): a fan at x 2..4
    /// blowing +x and a filter at x 10..12 (Darcy, permeability 0.05), through commands.
    fn fan_duct(s: &mut Session) {
        let periodic = json!({"type": "periodic"});
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Air", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Duct", "material": "Air", "min": [0, 0, 0], "max": [60, 4, 4], "spacing": 1,
                "faces": {"x_min": periodic, "x_max": periodic},
                "porous": [{"name": "Grille", "min": [10, 0, 0], "max": [12, 4, 4], "axis": "x", "resistance": {"type": "darcy_forchheimer", "permeability": 0.05}}],
                "fans": [{"name": "Fan", "min": [2, 0, 0], "max": [4, 4, 4], "direction": "+x", "curve": [[0, 0.002], [1.0, 0]], "relaxation": 0.05}]}},
            {"command": "model.step.create", "params": {"name": "Blow", "procedure": {"type": "fluid", "domain": "Duct", "steps": 3000, "time_step": 1, "frames": 2}}},
        ]);
        s.run_script(&script).unwrap();
    }

    #[test]
    fn a_fan_drives_a_duct_through_a_grille_and_the_reply_reports_the_operating_point() {
        let mut s = Session::default();
        fan_duct(&mut s);
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        let zone = |name: &str| r["zones"].as_array().unwrap().iter().find(|z| z["name"] == name).cloned().unwrap();
        let (fan, grille) = (zone("Fan"), zone("Grille"));
        assert_eq!(fan["kind"], "fan");
        assert_eq!(grille["kind"], "porous");
        let q = fan["flow"]["final"].as_f64().unwrap();
        let rise = fan["pressure_rise"]["final"].as_f64().unwrap();
        assert!(q > 0.0 && rise > 0.0 && rise < 0.002, "{fan}");
        // On its curve: Δp = 0.002 (1 − Q).
        assert!((rise - 0.002 * (1.0 - q)).abs() < 0.01 * rise, "{fan}");
        // What goes through the fan goes through the filter.
        assert!((grille["flow"]["final"].as_f64().unwrap() / q - 1.0).abs() < 0.01, "{grille}");
        // The measured differences are between the fluid layers next to each zone, so they
        // include the duct's wall shear over t + Δx = 3 (G per length). Around the loop:
        // the fan loses G·3 between its push and what is measured, the filter adds G·3 to its
        // own loss, and the duct between them (60 − 2·3 = 54) loses the rest.
        let drop = grille["pressure_drop"]["final"].as_f64().unwrap();
        let measured = fan["measured_pressure_rise"]["final"].as_f64().unwrap();
        let g = (rise - measured) / 3.0;
        assert!(g > 0.0, "{fan}");
        assert!(((measured - drop) / 54.0 / g - 1.0).abs() < 0.01, "loop momentum balance: {fan} {grille}");
        // The filter's own loss is Darcy's (μ/K) u t at the mean velocity Q / 16 (linear, so
        // the duct's profile does not matter), t = 2, and dominates the filter's drop. What
        // remains above it (about 6%) is wall shear the duct's G does not cover: the filter
        // flattens the profile, which shears the walls harder inside it, and the flow redevelops
        // just after it. This is an integration bound; the precise check of the resistance law
        // is wall-free: `a_porous_plate_drops_the_pressure_by_the_darcy_forchheimer_law` in
        // crates/lbm/tests/zones.rs (within 0.5%).
        let darcy = 0.1 / 0.05 * (q / 16.0) * 2.0;
        assert!(((drop - 3.0 * g) / darcy - 1.0).abs() < 0.1, "filter loss {} vs Darcy {darcy}", drop - 3.0 * g);
        assert!(darcy > 0.5 * drop, "the filter, not the walls, dominates its drop: {darcy} of {drop}");
        assert!(r["mass"]["relative_change"].as_f64().unwrap().abs() < 1e-10, "{r}");

        // Hostile zones are refused and leave the model as it was.
        let before = s.model.clone();
        let fan_with = |patch: Value| {
            let mut f = json!({"name": "Fan", "min": [2, 0, 0], "max": [4, 4, 4], "direction": "+x", "curve": [[0, 0.002], [1.0, 0]]});
            if let (Some(f), Some(p)) = (f.as_object_mut(), patch.as_object()) {
                f.extend(p.clone());
            }
            json!({"name": "Duct", "material": "Air", "min": [0, 0, 0], "max": [20, 4, 4], "spacing": 1, "fans": [f]})
        };
        for patch in [
            json!({"curve": [[0, 1]]}),
            json!({"curve": [[0, 1], [0, 0]]}),
            json!({"curve": [[0, 1], [1, 2]]}),
            json!({"curve": [[-1, 1], [1, 0]]}),
            json!({"curve": [[0, -1], [1, -2]]}),
            json!({"curve": vec![[0.0, 1.0]; 2000]}),
            json!({"curve": "fast"}),
            json!({"direction": "up"}),
            json!({"relaxation": 0}),
            json!({"relaxation": 5}),
            json!({"min": [5, 0, 0], "max": [4, 4, 4]}),
            json!({"name": ""}),
        ] {
            assert!(s.execute("model.fluid_domain.edit", &fan_with(patch.clone())).is_err(), "{patch}");
            assert_eq!(s.model, before);
        }
        for porous in [
            json!({"name": "G", "min": [10, 0, 0], "max": [12, 4, 4], "axis": "x", "resistance": {"type": "loss_coefficient", "coefficient": -1}}),
            json!({"name": "G", "min": [10, 0, 0], "max": [12, 4, 4], "axis": "w", "resistance": {"type": "loss_coefficient", "coefficient": 1}}),
            json!({"name": "G", "min": [10, 0, 0], "max": [12, 4, 4], "axis": "x", "resistance": {"type": "darcy_forchheimer", "permeability": 0}}),
            json!({"name": "G", "min": [10, 0, 0], "max": [12, 4, 4], "axis": "x", "resistance": {"type": "darcy_forchheimer", "permeability": 1, "inertial": -3}}),
        ] {
            let p = json!({"name": "Duct", "material": "Air", "min": [0, 0, 0], "max": [20, 4, 4], "spacing": 1, "porous": [porous]});
            assert!(s.execute("model.fluid_domain.edit", &p).is_err(), "{p}");
            assert_eq!(s.model, before);
        }
        // A fan outside the box, or one too strong for any time step, is refused when the
        // step solves.
        let outside = fan_with(json!({"min": [30, 0, 0], "max": [31, 4, 4]}));
        s.execute("model.fluid_domain.edit", &outside).unwrap();
        let e = s.execute("solve.fluid", &json!({})).unwrap_err().to_string();
        assert!(e.contains("Fan"), "{e}");
        // A rise past what SI holds (f64::MAX MPa) is refused as it is given; a huge one that
        // holds, when the step solves.
        assert!(s.execute("model.fluid_domain.edit", &fan_with(json!({"curve": [[0, f64::MAX], [1, f64::MIN]]}))).is_err());
        s.execute("model.fluid_domain.edit", &fan_with(json!({"curve": [[0, 1e290], [1, -1e290]]}))).unwrap();
        assert!(s.execute("solve.fluid", &json!({})).is_err());
    }

    #[test]
    fn a_fluid_step_stops_at_steady_state_and_says_why() {
        let mut s = Session::default();
        couette(&mut s);
        let full = s.execute("solve.fluid", &json!({"step": "Shear"})).unwrap();
        assert!(full["stop"].is_null(), "{full}");
        let steady =
            json!({"tolerance": 1e-5, "window": 1, "monitors": [{"type": "force", "faces": ["y_min"]}, {"type": "mass_flow", "face": "x_max"}]});
        s.execute(
            "model.step.create",
            &json!({"name": "Settle", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "frames": 4, "steady": steady}}),
        )
        .unwrap();
        let r = s.execute("solve.fluid", &json!({"step": "Settle"})).unwrap();
        let stop = &r["stop"];
        assert_eq!(stop["reason"], "steady", "{stop}");
        assert_eq!(stop["quantities"], json!(["velocity", "force y_min", "mass_flow x_max"]));
        assert_eq!(stop["steps"], r["steps"]);
        assert!(r["steps"].as_u64().unwrap() < full["steps"].as_u64().unwrap() / 2, "{r}");
        assert_eq!(r["frames"].as_array().unwrap().last(), Some(&r["end_time"]), "a frame at the stop: {r}");
        let history = stop["history"].as_array().unwrap();
        assert!(history[0]["changes"][0].is_null() && history.last().unwrap()["changes"][0].as_f64().unwrap() < 1e-5, "{stop}");
        // The wall shear μ U / H over the wall's area, as the full run gives it.
        let shear = history.last().unwrap()["values"][1][0].as_f64().unwrap();
        assert!((shear / (0.1 * 0.1 * 0.125 * 0.125) - 1.0).abs() < 1e-3, "{shear}");
        // The history is in the display units, which it names.
        assert_eq!(stop["units"], json!({"velocity": "mm/s", "force y_min": "N", "mass_flow x_max": "t/s"}), "{stop}");

        // Hostile rules are refused and leave the model as it was.
        let hostile = [
            json!({"tolerance": -1}),
            json!({"tolerance": 2}),
            json!({"tolerance": "tight"}),
            json!({"window": 0}),
            json!({"window": -1e308}),
            json!({"min_time": -1}),
            json!({"divergence": 0}),
            json!({"monitors": [{"type": "vorticity"}]}),
            json!({"monitors": [{"type": "force", "faces": []}]}),
            json!({"monitors": [{"type": "force", "faces": vec!["x_min"; 5000]}]}),
            json!({"monitors": [{"type": "mass_flow", "face": "top"}]}),
            json!({"monitors": vec![json!({"type": "mean_temperature"}); 500]}),
            json!({"monitors": 3}),
            json!({"oops": 1}),
            json!(7),
        ];
        for steady in &hostile {
            let before = s.model.clone();
            let p = json!({"name": "Bad", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "steady": steady}});
            assert!(s.execute("model.step.create", &p).is_err(), "{steady}");
            assert_eq!(s.model, before, "{steady}");
        }
        // Rules the model accepts but the flow cannot measure fail the solve, not the session.
        for (name, monitors) in [
            ("Ghost", json!([{"type": "force", "faces": [{"body": "Ghost", "face": "F1"}]}])),
            ("Cold", json!([{"type": "mean_temperature"}])),
            ("Huge", json!([{"type": "heat", "faces": ["y_min"]}])),
        ] {
            let steady = json!({"window": 1e300, "min_time": 1e300, "monitors": monitors});
            s.execute("model.step.create", &json!({"name": name, "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "steady": steady}}))
                .unwrap();
            assert!(s.execute("solve.fluid", &json!({"step": name})).is_err(), "{name}");
        }
        // A huge window and minimum time just run the full length.
        s.execute(
            "model.step.create",
            &json!({"name": "Long", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 50, "steady": {"window": 1e300, "min_time": 1e300}}}),
        )
        .unwrap();
        let r = s.execute("solve.fluid", &json!({"step": "Long"})).unwrap();
        assert_eq!(r["stop"]["reason"], "max_time", "{r}");
        assert_eq!(r["steps"], full["steps"]);
    }

    /// #650: `max_cells`, `max_cell_updates` and `max_frame_bytes` may be raised above their
    /// defaults up to what the session's memory bound holds; above it, or with no bound, they are
    /// refused with what to change; the reply (and `estimate`) gives the run's memory.
    #[test]
    fn fluid_caps_raise_up_to_the_machine_bound_and_refusals_say_what_to_change() {
        let mut s = Session::default();
        s.limits.memory = 32_000_000_000;
        couette(&mut s);
        // Raised above the defaults, within the bound: the run solves and reports its memory.
        let r = s.execute("solve.fluid", &json!({"max_frame_bytes": 8e9, "max_cells": 20_000_000, "max_cell_updates": 2e12})).unwrap();
        let m = &r["memory"];
        assert_eq!((m["cells"].as_u64(), m["frames"].as_u64(), m["raised"].as_bool()), (Some(8), Some(5), Some(true)), "{m}");
        assert_eq!(m["frame_bytes"].as_u64(), Some(5 * 8 * 16), "{m}");
        assert_eq!(m["bound_bytes"].as_u64(), Some(32_000_000_000), "{m}");
        assert_eq!(m["ceilings"]["max_frame_bytes"].as_u64(), Some(32_000_000_000), "{m}");
        // Up to the ceilings exactly.
        let ceilings = m["ceilings"].clone();
        for key in ["max_cells", "max_cell_updates", "max_frame_bytes"] {
            let most = ceilings[key].as_u64().unwrap();
            let r = s.execute("solve.fluid", &json!({key: most, "estimate": true}));
            assert!(r.is_ok(), "{key} {most}: {r:?}");
            let e = s.execute("solve.fluid", &json!({key: most + 1, "estimate": true})).unwrap_err().to_string();
            assert!(e.contains("this machine allows") && e.contains("coarsen the spacing or write fewer frames"), "{key}: {e}");
        }
        // `estimate` alone solves nothing.
        let solved = s.solved_steps().count();
        let only = s.execute("solve.fluid", &json!({"estimate": true})).unwrap();
        assert!(only.get("frames").is_none() && only["memory"]["total_bytes"].as_u64() > Some(0), "{only}");
        assert_eq!(s.solved_steps().count(), solved);
        // Over the default cap: the error says how many frames fit and what to raise it to.
        let e = s.execute("solve.fluid", &json!({"max_frame_bytes": 300})).unwrap_err().to_string();
        assert!(e.contains("2 frames fit: set the step's `frames` to at most 1"), "{e}");
        assert!(e.contains("raise `max_frame_bytes` to 640 (this machine allows up to"), "{e}");
        let e = s.execute("solve.fluid", &json!({"max_cells": 5})).unwrap_err().to_string();
        assert!(e.contains("coarsen the spacing, or raise `max_cells` to 8"), "{e}");
        // Where the machine's memory is unknown, the caps keep their defaults.
        s.limits.memory = 0;
        let e = s.execute("solve.fluid", &json!({"max_frame_bytes": 2_000_000_001u64})).unwrap_err().to_string();
        assert!(e.contains("cannot be raised"), "{e}");
        assert!(s.execute("solve.fluid", &json!({"max_frame_bytes": 2e9, "estimate": true})).is_ok());
        // The estimate as the Solve dialog shows it.
        s.limits.memory = 32_000_000_000;
        let e = fluid_estimate(&s, &json!({"max_frame_bytes": 300})).unwrap();
        assert!(e.summary().starts_with("Memory: 5 frames × 8 cells: 0 MB of frames"), "{}", e.summary());
        assert!(e.summary().contains("Runs with raised caps may take up to 32.0 GB on this machine"), "{}", e.summary());
        assert!(e.problem.as_deref().is_some_and(|p| p.contains("frames fit")), "{e:?}");
    }

    #[test]
    fn solve_fluid_hostile_params_never_panic() {
        let mut s = Session::default();
        couette(&mut s);
        s.execute("solve.fluid", &json!({"step": "Shear"})).unwrap();
        let hostile = [
            json!(null),
            json!({"step": 7}),
            json!({"step": "Nope"}),
            json!({"step": "Shear", "max_cells": 0}),
            json!({"step": "Shear", "max_cells": 1e300}),
            json!({"step": "Shear", "max_cells": 3}),
            json!({"step": "Shear", "max_steps": 3}),
            json!({"step": "Shear", "max_steps": -1}),
            json!({"step": "Shear", "max_cell_updates": 10}),
            json!({"step": "Shear", "max_cell_updates": 0.5}),
            json!({"step": "Shear", "max_frame_bytes": 16}),
            json!({"step": "Shear", "max_frame_bytes": "lots"}),
            json!({"step": "Shear", "max_frame_bytes": 1e300}),
            json!({"step": "Shear", "max_frame_bytes": -2e9}),
            json!({"step": "Shear", "max_frame_bytes": 2.5e9 + 0.5}),
            json!({"step": "Shear", "max_frame_bytes": u64::MAX}),
            json!({"step": "Shear", "max_cells": f64::MAX}),
            json!({"step": "Shear", "max_cell_updates": u64::MAX}),
            json!({"step": "Shear", "estimate": "yes"}),
            json!({"step": "Shear", "points": [[1e308, -1e308, 0]]}),
            json!({"step": "Shear", "points": [[f64::MAX, 0, 0]], "frame": 1e300}),
            json!({"step": "Shear", "points": [[0, 0, 0]], "frame": 99}),
            json!({"step": "Shear", "points": null}),
            json!({"step": "Shear", "points": vec![[0, 0, 0]; 20_000]}),
            json!({"step": "Shear", "forces": ["nope"]}),
            json!({"step": "Shear", "forces": vec!["y_min"; 1000]}),
            json!({"step": "Shear", "forces": "y_min"}),
            json!({"step": "Shear", "heat": ["y_min"]}),
            json!({"step": "Shear", "heat": vec!["y_min"; 1000]}),
            json!({"step": "Shear", "heat": 3}),
            json!({"step": "Shear", "nodes": [1]}),
            json!({"step": "Shear", "threads": 0}),
            json!({"step": "Shear", "threads": 100000}),
            json!({"step": "Shear", "threads": null}),
            json!({"step": "Shear", "backend": "tpu"}),
            json!({"step": "Shear", "backend": 7}),
            // no device backend is registered in this test binary: an error, not a panic
            json!({"step": "Shear", "backend": "gpu"}),
        ];
        for id in ["solve.fluid", "solve.probe", "solve.static", "solve.frequency"] {
            for p in &hostile {
                let mut t = s.clone();
                let before = t.model.clone();
                let results = t.solved_steps().count();
                if t.execute(id, p).is_err() {
                    assert_eq!(t.model, before);
                    assert_eq!(t.solved_steps().count(), results, "{id} {p}");
                }
            }
        }
    }

    /// A cantilever twisted and bent through a coupled reference point in a nonlinear step,
    /// typed in mm-t-s-N-MPa and in m-kg-s-N-Pa: the same model (stored in SI), so the same
    /// increments, iterations and convergence measures.
    #[test]
    fn a_coupled_nonlinear_step_converges_the_same_in_either_display_system() {
        let script = |system: &str, mm: f64, mpa: f64| {
            let (l, w) = (0.2 * mm, 0.02 * mm);
            let corners = json!([[0, 0, 0], [l, 0, 0], [l, w, 0], [0, w, 0], [0, 0, w], [l, 0, w], [l, w, w], [0, w, w]]);
            json!([
                {"command": "model.units", "params": {"system": system}},
                {"command": "mesh.block", "params": {"body": "Beam", "kind": "hex20", "corners": corners, "divisions": [4, 1, 1], "faces": {"6": "root", "4": "tip"}}},
                {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Beam"}]}},
                {"command": "model.set.create", "params": {"name": "root", "geometry": [{"body": "Beam", "face": "root"}]}},
                {"command": "model.set.create", "params": {"name": "tip", "geometry": [{"body": "Beam", "face": "tip"}]}},
                {"command": "model.material.create", "params": {"name": "Rubber", "elastic": {"type": "isotropic", "young": 1e7 / mpa, "poisson": 0.3}}},
                {"command": "model.section.create", "params": {"name": "beam", "set": "all", "material": "Rubber", "kind": {"type": "solid"}}},
                {"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}},
                {"command": "model.reference_point.create", "params": {"name": "Tip", "xyz": [l, 0.5 * w, 0.5 * w]}},
                {"command": "model.coupling.create", "params": {"name": "end", "reference_point": "Tip", "set": "tip", "kind": "kinematic"}},
                {"command": "model.load.create", "params": {"name": "shear", "type": "force", "target": {"reference_point": "Tip"}, "vector": [0, 0, -1]}},
                {"command": "model.load.create", "params": {"name": "twist", "type": "moment", "target": {"reference_point": "Tip"}, "vector": [0.05 * mm, 0, 0]}},
                {"command": "model.step.create", "params": {"name": "Twist", "procedure": {"type": "nonlinear_static", "initial_increment": 0.5}, "boundary_conditions": ["clamp"], "loads": ["shear", "twist"]}},
                {"command": "solve.nonlinear", "params": {"step": "Twist"}}
            ])
        };
        let solve = |system: &str, mm: f64, mpa: f64| {
            let mut s = Session::default();
            let out = s.run_script(&script(system, mm, mpa)).unwrap();
            (out.last().cloned().unwrap(), s.model)
        };
        let (millimetres, a) = solve("mm-t-s-N-MPa", 1000.0, 1e6);
        let (metres, b) = solve("m-kg-s-N-Pa", 1.0, 1.0);
        assert_eq!(b.display_units().map(|u| u.preset()), Some(fieldcraft_model::units::Preset::MKgS));
        assert_eq!(a.display_units().map(|u| u.preset()), Some(fieldcraft_model::units::Preset::MmTonneS), "the model keeps its display system");
        let iterations =
            |r: &Value| r["increments"].as_array().unwrap().iter().map(|i| i["iterations"].as_array().unwrap().len()).collect::<Vec<_>>();
        assert!(iterations(&millimetres).iter().all(|n| *n > 2), "{millimetres}");
        assert_eq!(iterations(&millimetres), iterations(&metres));
        // The same convergence measures, but for round-off: the block is meshed in the geometry's
        // millimetres from corners typed in either system.
        let measures = |r: &Value| -> Vec<f64> {
            let increments = r["increments"].as_array().unwrap();
            increments
                .iter()
                .flat_map(|i| i["iterations"].as_array().unwrap().iter())
                .flat_map(|it| [&it["force"], &it["displacement"]])
                .map(|v| v.as_f64().unwrap())
                .collect()
        };
        for (x, y) in measures(&millimetres).iter().zip(measures(&metres)) {
            assert!((x - y).abs() <= 1e-6 * x.abs().max(y.abs()) + 1e-12, "{x} vs {y}");
        }
    }

    #[test]
    fn solve_and_probe_a_coupled_reference_point() {
        // The defined bar: a force on a reference point kinematically coupled to the free end.
        let mut s = Session::default();
        crate::cmd::define::tests::defined_bar(&mut s);
        let r = s.execute("solve.static", &json!({"step": "Bend"})).unwrap();
        let rp = &r["reference_points"][0];
        assert_eq!(rp["name"], "RP", "{r}");
        assert!(rp["displacement"][2].as_f64().unwrap() < 0.0, "{r}");
        assert!((r["reaction_total"][2].as_f64().unwrap() - 100.0).abs() < 1e-6, "{r}");
        let label = rp["label"].as_u64().unwrap();
        let probe = s.execute("solve.probe", &json!({"step": "Bend", "nodes": [label, 1]})).unwrap();
        assert_eq!(probe["nodes"][0]["reference_point"], "RP", "{probe}");
        assert_eq!(probe["nodes"][0]["displacement"].as_array().unwrap().len(), 6, "{probe}");
        assert!(probe["nodes"][1]["reaction"].is_array(), "{probe}");
        s.execute("solve.frequency", &json!({"step": "Modes"})).unwrap();
        let shape = s.execute("solve.probe", &json!({"step": "Modes", "nodes": [label], "mode": 1})).unwrap();
        assert_eq!(shape["nodes"][0]["shape"].as_array().unwrap().len(), 6, "{shape}");
        for p in [json!({"step": "Bend", "nodes": [label, label, 0]}), json!({"step": "Modes", "nodes": [label], "mode": 99})] {
            assert!(s.execute("solve.probe", &p).is_err(), "{p}");
        }
    }

    /// A cantilever L × 1 × 1 along x of `n[0]` × `n[1]` × `n[1]` hex8 elements, E = 1000, ν = 0,
    /// ρ = 1, clamped at x = 0; one node set per cross-section (`x0`, `x1`, …) and the tip's
    /// corner node `tip`; a frequency step `Modes` and the section `beam`.
    fn explicit_cantilever(s: &mut Session, length: f64, n: [usize; 2]) {
        use fieldcraft_model::{ElementKind, Set};
        // The mesh is built in the model as stored (SI); the commands' numbers are SI too.
        s.default_units = crate::units::si().clone();
        let [nx, ny] = n;
        let label = |i: usize, j: usize, k: usize| (1 + i + (nx + 1) * (j + (ny + 1) * k)) as u32;
        for k in 0..=ny {
            for j in 0..=ny {
                for i in 0..=nx {
                    s.model.set_node(label(i, j, k), [length * i as f64 / nx as f64, j as f64 / ny as f64, k as f64 / ny as f64]).unwrap();
                }
            }
        }
        let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        let mut e = 1;
        for k in 0..ny {
            for j in 0..ny {
                for i in 0..nx {
                    s.model.set_element(e, ElementKind::Hex8, corners.iter().map(|c| label(i + c[0], j + c[1], k + c[2])).collect()).unwrap();
                    e += 1;
                }
            }
        }
        for i in 0..=nx {
            let nodes = (0..=ny).flat_map(|k| (0..=ny).map(move |j| label(i, j, k))).collect();
            s.model.add_set(&format!("x{i}"), Set { nodes, ..Set::default() }).unwrap();
        }
        s.model.add_set("tip", Set { nodes: vec![label(nx, 0, 0)], ..Set::default() }).unwrap();
        s.model.add_set("all", Set { elements: (1..e).collect(), ..Set::default() }).unwrap();
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0}, "density": 1}},
            {"command": "model.section.create", "params": {"name": "beam", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {"incompatible_modes": true}}},
            {"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "x0"}}},
            {"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 1}, "boundary_conditions": ["clamp"]}},
        ]);
        s.run_script(&script).unwrap();
    }

    /// The stored result of a step, hashed bit for bit (Debug prints every float exactly).
    fn result_hash(s: &Session, step: &str) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::hash::DefaultHasher::new();
        format!("{:?}", s.result(step).unwrap()).hash(&mut h);
        h.finish()
    }

    /// Implicit (static, frequency), explicit and lattice Boltzmann solves large enough for every
    /// parallel path (factorization tasks and tiles, element rounds, cell blocks) give the same
    /// bits on 1, 2 and 7 threads.
    #[test]
    fn results_do_not_depend_on_the_thread_count() {
        let solve_all = |threads: usize| -> Vec<u64> {
            let mut s = Session::default();
            explicit_cantilever(&mut s, 10.0, [30, 6]);
            let script = json!([
                {"command": "model.load.create", "params": {"name": "push", "type": "force", "target": {"set": "x30"}, "vector": [0, 0.3, -1]}},
                {"command": "model.step.create", "params": {"name": "Bend", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["push"]}},
            ]);
            s.run_script(&script).unwrap();
            s.execute("solve.static", &json!({"step": "Bend", "threads": threads})).unwrap();
            s.execute("solve.frequency", &json!({"step": "Modes", "threads": threads})).unwrap();
            let mut hashes = vec![result_hash(&s, "Bend"), result_hash(&s, "Modes")];
            let step = json!({
                "type": "explicit_dynamic", "end_time": 0.15, "output_interval": 0.05,
                "initial_velocities": [{"set": "x30", "velocity": [0, 0, -2]}], "hourglass": {"type": "stiffness"}
            });
            let script = json!([
                {"command": "model.section.edit", "params": {"name": "beam", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {"integration": "reduced"}}},
                {"command": "model.step.create", "params": {"name": "Hit", "procedure": step, "boundary_conditions": ["clamp"]}},
            ]);
            s.run_script(&script).unwrap();
            s.execute("solve.explicit", &json!({"step": "Hit", "threads": threads})).unwrap();
            hashes.push(result_hash(&s, "Hit"));
            let mut f = Session::default();
            let periodic = json!({"type": "periodic"});
            let script = json!([
                {"command": "model.material.create", "params": {"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}}},
                {"command": "model.fluid_domain.create", "params": {"name": "Gap", "material": "Oil", "min": [0, 0, 0], "max": [5, 1, 2], "spacing": 0.125,
                    "faces": {"z_min": periodic, "z_max": periodic}}},
                {"command": "model.fluid_domain.boundary", "params": {"domain": "Gap", "at": "y_max", "condition": {"type": "moving_wall", "velocity": [0.1, 0, 0.02]}}},
                {"command": "model.step.create", "params": {"name": "Shear", "procedure": {"type": "fluid", "domain": "Gap", "end_time": 2, "frames": 2}}},
            ]);
            f.run_script(&script).unwrap();
            f.execute("solve.fluid", &json!({"step": "Shear", "threads": threads})).unwrap();
            hashes.push(result_hash(&f, "Shear"));
            hashes
        };
        let one = solve_all(1);
        for threads in [2, 7] {
            assert_eq!(solve_all(threads), one, "{threads} threads");
        }
        let mut s = Session::default();
        cantilever(&mut s);
        for bad in [json!(0), json!(257), json!(-1), json!(1.5), json!("8")] {
            let e = s.execute("solve.static", &json!({"threads": bad})).unwrap_err().to_string();
            assert!(e.contains("`threads` must be a whole number from 1 to 256"), "{e}");
        }
        s.limits.threads = 3;
        s.execute("solve.static", &json!({})).unwrap();
    }

    /// Upward zero crossings of a sampled signal, linearly interpolated.
    fn crossings(time: &[f64], x: &[f64]) -> Vec<f64> {
        let mut out = Vec::new();
        for k in 1..x.len() {
            let (a, b) = (x[k - 1], x[k]);
            if a < 0.0 && b >= 0.0 {
                out.push(time[k - 1] + (time[k] - time[k - 1]) * (-a) / (b - a));
            }
        }
        out
    }

    /// Free vibration of a cantilever after an initial velocity in the shape of its static
    /// deflection: the period of the explicit run (one-point hex8) against the first natural
    /// frequency of the implicit modal solver (incompatible-mode hex8, lumped mass) on the same
    /// mesh. The one-point element's mean strain misses the linear variation of the bending
    /// strain through each element: with n elements through the depth it keeps 1 − 1/n² of the
    /// bending energy, so without hourglass stiffness (viscous control) the period is longer by
    /// (1 − 1/n²)^{−½} (3.3% for n = 4); the stiffness form restores part of it.
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn explicit_cantilever_period_matches_the_modal_frequency() {
        let (length, n) = (10.0, [40usize, 4usize]);
        let mut s = Session::default();
        explicit_cantilever(&mut s, length, n);
        let modes = s.execute("solve.frequency", &json!({"mass": "lumped"})).unwrap();
        let period = 1.0 / modes["modes"][0]["frequency"].as_f64().unwrap();
        // Euler–Bernoulli: f₁ = 1.875² / (2π L²) √(E I / (ρ A)), I/A = 1/12.
        let euler = 2.0 * std::f64::consts::PI * length * length / 1.875104f64.powi(2) / (1000.0f64 / 12.0).sqrt();
        assert!((period / euler - 1.0).abs() < 0.01, "{period} vs {euler}");
        let section = json!({"name": "beam", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {"integration": "reduced"}});
        s.execute("model.section.edit", &section).unwrap();
        // v(x) ∝ x²(3L − x): the static deflection of a tip load, close to the first mode.
        let velocities: Vec<Value> = (0..=n[0])
            .map(|i| {
                let x = length * i as f64 / n[0] as f64;
                json!({"set": format!("x{i}"), "velocity": [0, 0, 0.01 * x * x * (3.0 * length - x) / (2.0 * length.powi(3))]})
            })
            .collect();
        let one_point = period / (1.0 - 1.0 / (n[1] * n[1]) as f64).sqrt();
        let mut measured = Vec::new();
        for hourglass in ["viscous", "stiffness"] {
            let step = json!({
                "type": "explicit_dynamic", "end_time": 4.2 * period, "history_interval": period / 400.0, "output_interval": period,
                "history_nodes": "tip", "initial_velocities": velocities, "hourglass": {"type": hourglass}
            });
            let name = format!("Swing-{hourglass}");
            s.execute("model.step.create", &json!({"name": name, "procedure": step, "boundary_conditions": ["clamp"]})).unwrap();
            let r = s.execute("solve.explicit", &json!({"step": name})).unwrap();
            let h = s.execute("solve.history", &json!({"step": name})).unwrap();
            let time: Vec<f64> = h["time"].as_array().unwrap().iter().map(|t| t.as_f64().unwrap()).collect();
            let uz: Vec<f64> = h["nodes"][0]["displacement"].as_array().unwrap().iter().map(|d| d[2].as_f64().unwrap()).collect();
            // The tip starts at zero moving up: it crosses zero upwards once per period after that.
            let ups = crossings(&time, &uz);
            let t = (ups[ups.len() - 1] - ups[0]) / (ups.len() - 1) as f64;
            eprintln!(
                "cantilever {n:?}, {hourglass} hourglass: explicit period {t:.5} vs modal {period:.5} ({:+.3}%), one-point estimate {one_point:.5} ({:+.3}%), Euler–Bernoulli {euler:.5}; energy balance {:.1e}, {} increments",
                100.0 * (t / period - 1.0),
                100.0 * (t / one_point - 1.0),
                r["energy_balance"].as_f64().unwrap(),
                r["increments"]
            );
            assert!(ups.len() >= 4, "{ups:?}");
            assert!(r["energy_balance"].as_f64().unwrap() < 0.01, "{r}");
            // The beam's material is deliberately unit-free (E = 1), so the only warnings are about
            // its units: none from the solve itself.
            let about_units = |w: &str| w.contains("unit system") || w.contains("check the unit it was given in");
            assert!(r["warnings"].as_array().unwrap().iter().all(|w| about_units(w.as_str().unwrap())), "{r}");
            measured.push(t);
        }
        // Viscous control adds no stiffness: the one-point estimate holds closely.
        assert!((measured[0] / one_point - 1.0).abs() < 0.005, "{measured:?} vs {one_point}");
        // The stiffness form lies between it and the modal period, and within 3% of the latter.
        assert!(measured[1] > period && measured[1] < measured[0], "{measured:?}");
        assert!((measured[1] / period - 1.0).abs() < 0.03, "{measured:?} vs {period}");
    }

    /// A small impact model through commands: a hex8 bar moving onto a rigid wall.
    pub(crate) fn explicit_impact(s: &mut Session) {
        explicit_cantilever(s, 4.0, [8, 1]);
        let step = json!({
            "type": "explicit_dynamic", "end_time": 0.2, "history_interval": 0.002, "output_interval": 0.05, "history_nodes": "tip",
            "initial_velocities": [{"set": "all_nodes", "velocity": [-0.1, 0, 0]}],
            "rigid_walls": [{"set": "x0", "point": [0, 0, 0], "normal": [1, 0, 0]}]
        });
        let nodes: Vec<u32> = s.model.nodes().keys().copied().collect();
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all_nodes", "nodes": nodes}},
            {"command": "model.section.edit", "params": {"name": "beam", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {"integration": "reduced"}}},
            {"command": "model.step.create", "params": {"name": "Impact", "procedure": step}},
        ]);
        s.run_script(&script).unwrap();
    }

    #[test]
    fn explicit_results_are_kept_probed_and_guarded() {
        let mut s = Session::default();
        explicit_impact(&mut s);
        let journal = s.journal.len();
        let r = s.execute("solve.explicit", &json!({"step": "Impact"})).unwrap();
        assert_eq!((r["completed"].as_bool(), r["frames"].as_u64()), (Some(true), Some(5)), "{r}");
        // Eight elements, so eight controlling ones; on this coarse mesh the impact excites modes
        // near the stability limit, whose energy the integer-step account measures only roughly.
        assert_eq!(r["controlling_elements"].as_array().unwrap().len(), 8, "{r}");
        assert!(r["energy_balance"].as_f64().unwrap() < 0.1, "{r}");
        assert_eq!(s.journal.len(), journal, "solving does not change the model");
        let probe = s.execute("solve.probe", &json!({"step": "Impact", "nodes": [1, 9], "frame": 0})).unwrap();
        assert_eq!(probe["nodes"][1]["velocity"][0].as_f64(), Some(-0.1), "{probe}");
        let last = s.execute("solve.probe", &json!({"step": "Impact", "nodes": [9]})).unwrap();
        assert_eq!(last["frame"], 4, "{last}");
        let el = s.execute("solve.probe", &json!({"step": "Impact", "elements": [1, 8]})).unwrap();
        let e = &el["elements"][0];
        assert!(e["stress"].as_array().unwrap().len() == 6 && e["von_mises"].as_f64().unwrap() >= 0.0 && e["eroded"] == false, "{el}");
        assert!(e["plastic_strain"].is_null(), "an elastic model has no plastic strain: {el}");
        assert_eq!(r["eroded"].as_array().map(Vec::len), Some(0), "{r}");
        let h = s.execute("solve.history", &json!({"step": "Impact", "every": 2})).unwrap();
        assert_eq!(h["time"].as_array().unwrap().len(), h["energies"].as_array().unwrap().len());
        assert!(h["nodes"][0]["velocity"].as_array().unwrap().len() > 5, "{h}");
        assert_eq!(s.solved_steps().collect::<Vec<_>>(), ["Impact"]);
        assert!(s.execute("solve.history", &json!({"step": "Modes"})).is_err(), "not solved");

        let hostile = [
            json!(null),
            json!({"step": 7}),
            json!({"step": "Nope"}),
            json!({"step": "Modes"}),
            json!({"step": "Impact", "max_dofs": 0}),
            json!({"step": "Impact", "max_dofs": 5}),
            json!({"step": "Impact", "max_increments": 0}),
            json!({"step": "Impact", "max_increments": 3}),
            json!({"step": "Impact", "max_increments": -1}),
            json!({"step": "Impact", "max_increments": 1e300}),
            json!({"step": "Impact", "max_output_values": 1}),
            json!({"step": "Impact", "every": 0}),
            json!({"step": "Impact", "every": -3}),
            json!({"step": "Impact", "every": 1e300}),
            json!({"step": "Impact", "every": u64::MAX}),
            json!({"step": "Impact", "nodes": [1], "frame": 99}),
            json!({"step": "Impact", "nodes": [1], "frame": -1}),
            json!({"step": "Impact", "nodes": [0, 4294967296u64]}),
            json!({"step": "Impact"}),
            json!({"step": "Impact", "elements": [999]}),
            json!({"step": "Impact", "elements": [0, -1, 1.5]}),
            json!({"step": "Impact", "elements": vec![1; 20_000]}),
            json!({"step": "Impact", "elements": "all", "nodes": [1]}),
        ];
        for id in ["solve.explicit", "solve.history", "solve.probe"] {
            for p in &hostile {
                let mut t = s.clone();
                let before = t.model.clone();
                let results = t.solved_steps().count();
                if t.execute(id, p).is_err() {
                    assert_eq!(t.model, before);
                    assert_eq!(t.solved_steps().count(), results, "{id} {p}");
                }
            }
        }
        // Hostile steps are refused when they are defined; solving one that a later edit makes
        // impossible (a node behind its wall) is an error, not a panic.
        for bad in [
            json!({"type": "explicit_dynamic", "end_time": -1}),
            json!({"type": "explicit_dynamic", "end_time": 1, "dt_scale": 2}),
            json!({"type": "explicit_dynamic", "end_time": 1, "max_increments": 0}),
            json!({"type": "explicit_dynamic", "end_time": 1, "output_interval": 1e-300}),
            json!({"type": "explicit_dynamic", "end_time": 1e308, "history_interval": 1e-308}),
            json!({"type": "explicit_dynamic", "end_time": 1, "rigid_walls": [{"set": "x0", "point": [0, 0, 0], "normal": [0, 0, 0]}]}),
            json!({"type": "explicit_dynamic", "end_time": 1, "initial_velocities": [{"set": "nope", "velocity": [1, 0, 0]}]}),
            json!({"type": "explicit_dynamic", "end_time": 1, "hourglass": {"type": "magic"}}),
        ] {
            assert!(s.execute("model.step.create", &json!({"name": "Bad", "procedure": bad})).is_err(), "{bad}");
        }
        let behind = json!({"type": "explicit_dynamic", "end_time": 0.1, "rigid_walls": [{"set": "x0", "point": [1, 0, 0], "normal": [1, 0, 0]}]});
        s.execute("model.step.create", &json!({"name": "Behind", "procedure": behind})).unwrap();
        assert!(s.execute("solve.explicit", &json!({"step": "Behind"})).unwrap_err().to_string().contains("behind"));
        let long = json!({"type": "explicit_dynamic", "end_time": 1e6});
        s.execute("model.step.create", &json!({"name": "Long", "procedure": long})).unwrap();
        assert!(s.execute("solve.explicit", &json!({"step": "Long"})).unwrap_err().to_string().contains("increments"));
        assert!(s.execute("model.set.delete", &json!({"name": "all_nodes"})).unwrap_err().to_string().contains("step `Impact`"));
        // Defining other steps left what the impact step reads alone: its result is current.
        assert_eq!(s.solved_steps().collect::<Vec<_>>(), ["Impact"]);
    }

    /// A mass on a spring and a dashpot, defined, solved and written through commands; editing
    /// and deleting them goes through the same checks as every definition.
    #[test]
    fn spring_mass_through_commands() {
        let mut s = Session::default();
        let script = json!([
            {"command": "mesh.node.add", "params": {"label": 1, "xyz": [0, 0, 0]}},
            {"command": "model.set.create", "params": {"name": "n", "nodes": [1]}},
            {"command": "model.mass.create", "params": {"name": "m", "target": {"set": "n"}, "mass": 4}},
            {"command": "model.spring.create", "params": {"name": "k", "type": "spring", "first": {"set": "n"}, "direction": {"dof": 1}, "stiffness": 400}},
            {"command": "model.spring.create", "params": {"name": "c", "type": "dashpot", "first": {"set": "n"}, "direction": {"vector": [1, 0, 0]}, "coefficient": 1}},
            {"command": "model.boundary_condition.create", "params": {"name": "yz", "type": "displacement", "target": {"set": "n"}, "values": [null, 0, 0, null, null, null]}},
            {"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 1}, "boundary_conditions": ["yz"]}},
        ]);
        s.run_script(&script).unwrap();
        let r = s.execute("solve.frequency", &json!({})).unwrap();
        let f = r["modes"][0]["frequency"].as_f64().unwrap();
        assert!((f - 10.0 / (2.0 * std::f64::consts::PI)).abs() < 1e-12, "{r}");
        assert!(r["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("dashpots `c`")), "{r}");
        assert!(s.execute("model.mass.edit", &json!({"name": "m", "target": {"set": "n"}, "mass": -1})).is_err());
        assert!(s.execute("model.set.delete", &json!({"name": "n"})).is_err(), "the set is in use");
        s.execute("model.mass.edit", &json!({"name": "m", "target": {"set": "n"}, "mass": 1})).unwrap();
        let r = s.execute("solve.frequency", &json!({})).unwrap();
        assert!((r["modes"][0]["frequency"].as_f64().unwrap() - 20.0 / (2.0 * std::f64::consts::PI)).abs() < 1e-12);
        s.execute("model.spring.delete", &json!({"name": "c"})).unwrap();
        assert!(s.model.springs().len() == 1 && s.model.mass("m").is_some());
    }
    #[test]
    fn solve_and_probe_shells_and_beams() {
        // A cantilever strip of two quad4 shells clamped at x = 0, stiffened along its tip edge
        // by a beam; a downward force at a tip corner.
        let mut s = Session::default();
        let mut script = Vec::new();
        for (label, xyz) in [(1, [0, 0, 0]), (2, [10, 0, 0]), (3, [20, 0, 0]), (4, [0, 5, 0]), (5, [10, 5, 0]), (6, [20, 5, 0])] {
            script.push(json!({"command": "mesh.node.add", "params": {"label": label, "xyz": xyz}}));
        }
        script.extend([
            json!({"command": "mesh.element.add", "params": {"label": 1, "kind": "quad4", "nodes": [1, 2, 5, 4]}}),
            json!({"command": "mesh.element.add", "params": {"label": 2, "kind": "quad4", "nodes": [2, 3, 6, 5]}}),
            json!({"command": "mesh.element.add", "params": {"label": 3, "kind": "beam2", "nodes": [3, 6]}}),
            json!({"command": "model.set.create", "params": {"name": "plate", "elements": [1, 2]}}),
            json!({"command": "model.set.create", "params": {"name": "stiffener", "elements": [3]}}),
            json!({"command": "model.set.create", "params": {"name": "root", "nodes": [1, 4]}}),
            json!({"command": "model.set.create", "params": {"name": "corner", "nodes": [3]}}),
            json!({"command": "model.material.create", "params": {"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}, "density": 7.85e-9}}),
            json!({"command": "model.section.create", "params": {"name": "skin", "set": "plate", "material": "Steel", "kind": {"type": "shell", "thickness": 1}}}),
            json!({"command": "model.section.create", "params": {"name": "rib", "set": "stiffener", "material": "Steel", "kind": {"type": "beam", "profile": {"type": "rectangle", "width": 1, "height": 2}, "orientation": [1, 0, 0]}}}),
            json!({"command": "model.boundary_condition.create", "params": {"name": "clamp", "type": "fixed", "target": {"set": "root"}}}),
            json!({"command": "model.load.create", "params": {"name": "push", "type": "force", "target": {"set": "corner"}, "vector": [0, 0, -1]}}),
            json!({"command": "model.step.create", "params": {"name": "Bend", "procedure": {"type": "static"}, "boundary_conditions": ["clamp"], "loads": ["push"]}}),
            json!({"command": "model.step.create", "params": {"name": "Modes", "procedure": {"type": "frequency", "modes": 3}, "boundary_conditions": ["clamp"]}}),
        ]);
        s.run_script(&Value::Array(script)).unwrap();
        let r = s.execute("solve.static", &json!({"step": "Bend"})).unwrap();
        assert!((r["reaction_total"][2].as_f64().unwrap() - 1.0).abs() < 1e-9, "{r}");
        let probe = s.execute("solve.probe", &json!({"step": "Bend", "nodes": [3, 1], "elements": [1, 3]})).unwrap();
        assert!(probe["nodes"][0]["displacement"][2].as_f64().unwrap() < 0.0, "{probe}");
        assert!(probe["nodes"][0]["rotation"][0].as_f64().unwrap() != 0.0, "{probe}");
        assert!(probe["nodes"][1]["reaction_moment"].is_array(), "{probe}");
        assert_eq!(probe["elements"][0]["forces"].as_array().unwrap().len(), 4, "{probe}");
        assert_eq!(probe["elements"][1]["ends"].as_array().unwrap().len(), 2, "{probe}");
        let only = s.execute("solve.probe", &json!({"step": "Bend", "elements": [2]})).unwrap();
        assert_eq!(only["nodes"].as_array().unwrap().len(), 0, "{only}");
        s.execute("solve.frequency", &json!({"step": "Modes"})).unwrap();
        let shape = s.execute("solve.probe", &json!({"step": "Modes", "nodes": [3], "mode": 1})).unwrap();
        assert!(shape["nodes"][0]["rotation"].is_array(), "{shape}");
        for p in [
            json!({"step": "Bend", "elements": [99]}),
            json!({"step": "Bend", "elements": [0]}),
            json!({"step": "Bend", "elements": "1"}),
            json!({"step": "Bend", "elements": vec![1; 20_000]}),
            json!({"step": "Modes", "nodes": [3], "elements": [1]}),
        ] {
            assert!(s.execute("solve.probe", &p).is_err(), "{p}");
        }
    }

    /// A heat balance per load: the convection `Air` carries away what the heater's prescribed
    /// temperature supplies (#362).
    #[test]
    fn heat_flows_balance_per_load() {
        let mut s = Session::default();
        cantilever(&mut s);
        let script = json!([
            {"command": "model.material.edit", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0}, "density": 1, "thermal": {"conductivity": 2}}},
            {"command": "model.boundary_condition.create", "params": {"name": "Heater", "type": "temperature", "target": {"set": "root"}, "value": 100}},
            {"command": "model.set.create", "params": {"name": "skin", "faces": (1..=6).map(|f| json!({"element": 10, "face": f})).collect::<Vec<_>>()}},
            {"command": "model.load.create", "params": {"name": "Air", "type": "convection", "set": "skin", "film": 0.5, "sink": 20}},
            {"command": "model.step.create", "params": {"name": "Therm", "procedure": {"type": "heat_transfer"}, "boundary_conditions": ["Heater"], "loads": ["Air"]}},
        ]);
        s.run_script(&script).unwrap();
        let r = s.execute("solve.heat", &json!({"step": "Therm"})).unwrap();
        let flows = &r["heat_flows"];
        let air = flows["loads"].as_array().unwrap().iter().find(|l| l["name"] == "Air").unwrap();
        assert_eq!(air["type"], "convection");
        let (out, supplied) = (air["heat"].as_f64().unwrap(), flows["reactions"].as_f64().unwrap());
        assert!(out < 0.0 && supplied > 0.0, "{flows}");
        assert!(flows["imbalance"].as_f64().unwrap().abs() < 1e-6 * supplied, "{flows}");
        assert_eq!(s.results_report()["heat_flows"], *flows);
    }

    #[test]
    fn solve_heat_and_probe_temperatures() {
        // The cantilever bar as a conductor: 0 at the root, 100 at the tip, linear between.
        let mut s = Session::default();
        cantilever(&mut s);
        let script = json!([
            {"command": "model.material.edit", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0}, "density": 1, "thermal": {"conductivity": 2, "specific_heat": 3, "expansion": 1e-5}}},
            {"command": "model.boundary_condition.create", "params": {"name": "cold", "type": "temperature", "target": {"set": "root"}, "value": 0}},
            {"command": "model.boundary_condition.create", "params": {"name": "hot", "type": "temperature", "target": {"set": "tip"}, "value": 100}},
            {"command": "model.step.create", "params": {"name": "Heat", "procedure": {"type": "heat_transfer"}, "boundary_conditions": ["cold", "hot"]}},
            {"command": "model.step.create", "params": {"name": "Warm", "procedure": {"type": "heat_transfer", "transient": {"end_time": 10, "increment": 1}}, "boundary_conditions": ["cold", "hot"]}},
        ]);
        s.run_script(&script).unwrap();
        let r = s.execute("solve.heat", &json!({"step": "Heat"})).unwrap();
        assert!((r["max_temperature"]["value"].as_f64().unwrap() - 100.0).abs() < 1e-9, "{r}");
        assert!(r["min_temperature"]["value"].as_f64().unwrap().abs() < 1e-9, "{r}");
        let probe = s.execute("solve.probe", &json!({"step": "Heat", "nodes": [6]})).unwrap();
        assert!((probe["nodes"][0]["temperature"].as_f64().unwrap() - 50.0).abs() < 1e-9, "{probe}");
        assert!((probe["nodes"][0]["heat_flux"][0].as_f64().unwrap() + 20.0).abs() < 1e-9, "{probe}");
        let w = s.execute("solve.heat", &json!({"step": "Warm"})).unwrap();
        assert_eq!(w["frames"], 11, "{w}");
        let early = s.execute("solve.probe", &json!({"step": "Warm", "nodes": [6], "frame": 1})).unwrap();
        assert_eq!(early["time"], 1.0);
        assert!(s.execute("solve.probe", &json!({"step": "Warm", "nodes": [6], "frame": 99})).is_err());
        assert!(s.execute("solve.probe", &json!({"step": "Warm", "nodes": [999]})).is_err());
        assert!(s.execute("solve.heat", &json!({"step": "Bend"})).is_err(), "not a heat step");
        assert!(s.execute("solve.heat", &json!({"step": "Heat", "max_dofs": 0})).is_err());
        // The results display shows temperatures, undeformed.
        let view = crate::results::ResultsView::for_result(&s.latest_result("Heat").unwrap().0.clone());
        assert_eq!((view.contour, view.deformed), (Some(crate::results::Field::Temperature), false));
        let (heat, _) = s.latest_result("Heat").unwrap();
        let fields: Vec<_> = crate::results::Field::available(heat).collect();
        assert_eq!(fields, vec![crate::results::Field::Temperature, crate::results::Field::HeatFlux]);
        assert!(crate::results::displacements(heat, crate::results::Pick::default()).is_err());
    }

    /// `solve.nonlinear`: a cantilever bent far by a dead tip load, with its increment history,
    /// a frame per converged increment shown by `results.frame`, `solve.run` dispatching to it,
    /// and the same bits on 1, 2 and 7 threads.
    #[test]
    fn nonlinear_steps_solve_with_their_history_and_frames() {
        let solve = |threads: usize| -> (Value, u64, Session) {
            let mut s = Session::default();
            cantilever(&mut s);
            let script = json!([
                {"command": "model.section.edit", "params": {"name": "bar", "set": "all", "material": "Mat", "kind": {"type": "solid"}, "formulation": {}}},
                {"command": "model.step.create", "params": {"name": "Large", "procedure": {"type": "nonlinear_static", "initial_increment": 0.25},
                    "boundary_conditions": ["clamp"], "loads": ["push"]}},
                {"command": "model.load.edit", "params": {"name": "push", "type": "force", "target": {"set": "tip"}, "vector": [0, 0, -15]}},
            ]);
            s.run_script(&script).unwrap();
            let r = s.execute("solve.nonlinear", &json!({"step": "Large", "threads": threads})).unwrap();
            let hash = result_hash(&s, "Large");
            (r, hash, s)
        };
        let (r, one, mut s) = solve(1);
        for threads in [2, 7] {
            assert_eq!(solve(threads).1, one, "{threads} threads");
        }
        let increments = r["increments"].as_array().unwrap();
        assert!(increments.iter().any(|i| i["converged"] == true) && increments.last().unwrap()["load_factor"] == 1.0, "{r}");
        assert!(increments[0]["iterations"].as_array().is_some_and(|it| !it.is_empty() && it[0]["force"].is_number()), "{r}");
        let frames = r["frames"].as_array().unwrap().len();
        assert_eq!(frames, 1 + increments.iter().filter(|i| i["converged"] == true).count());
        // Far from linear: P L² / (E I) = 15·100/(1000/12) = 18, the tip swings well down and in.
        let tip = r["max_displacement"]["value"].as_f64().unwrap();
        assert!(tip > 8.0 && tip < 12.0, "{r}");
        let view = s.execute("results.frame", &json!({"step": "Large", "frame": 1})).unwrap();
        assert_eq!(view["view"]["frame"], 1, "{view}");
        assert!(view["load_factor"].as_f64().is_some_and(|f| f > 0.0 && f < 1.0), "{view}");
        assert!(s.execute("results.frame", &json!({"step": "Large", "frame": frames})).is_err());
        let probe = s.execute("solve.probe", &json!({"step": "Large", "nodes": [11]})).unwrap();
        assert!(probe["nodes"][0]["displacement"][2].as_f64().is_some_and(|w| w < -5.0), "{probe}");
        // solve.run solves every step with its own solver; the hostile caps are refused.
        assert!(s.execute("solve.nonlinear", &json!({"step": "Large", "max_newton_iterations": 0})).is_err());
        assert!(s.execute("solve.nonlinear", &json!({"step": "Large", "max_newton_iterations": 1e300})).is_err());
        let e = s.execute("solve.nonlinear", &json!({"step": "Large", "max_newton_iterations": 2})).unwrap_err().to_string();
        assert!(e.contains("Newton iterations"), "{e}");
        assert!(s.execute("solve.nonlinear", &json!({"step": "Bend"})).is_err(), "a linear step");
        let all = s.execute("solve.run", &json!({})).unwrap();
        assert_eq!(all["steps"].as_array().map(Vec::len), Some(3), "{all}");
    }

    /// Two blocks pressed together across a contact pair, as commands.
    fn pressed_blocks(s: &mut Session) {
        let script = json!([
            {"command": "mesh.block", "params": {"body": "Lower", "kind": "hex8", "corners": [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], "divisions": [2, 2, 1], "faces": {"1": "bottom", "2": "top"}}},
            {"command": "mesh.block", "params": {"body": "Upper", "kind": "hex8", "corners": [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0, 0, 2], [1, 0, 2], [1, 1, 2], [0, 1, 2]], "divisions": [3, 3, 1], "faces": {"1": "bottom", "2": "top"}, "merge": 0}},
            {"command": "model.set.create", "params": {"name": "all", "geometry": [{"body": "Lower"}, {"body": "Upper"}]}},
            {"command": "model.set.create", "params": {"name": "lower_top", "geometry": [{"body": "Lower", "face": "top"}]}},
            {"command": "model.set.create", "params": {"name": "upper_bottom", "geometry": [{"body": "Upper", "face": "bottom"}]}},
            {"command": "model.set.create", "params": {"name": "base", "geometry": [{"body": "Lower", "face": "bottom"}]}},
            {"command": "model.set.create", "params": {"name": "top", "geometry": [{"body": "Upper", "face": "top"}]}},
            {"command": "model.material.create", "params": {"name": "Mat", "elastic": {"type": "isotropic", "young": 1000.0, "poisson": 0.0}}},
            {"command": "model.section.create", "params": {"name": "all", "set": "all", "material": "Mat", "kind": {"type": "solid"}}},
            {"command": "interaction.contact.pair", "params": {"name": "interface", "main": "lower_top", "secondary": "upper_bottom"}},
            {"command": "model.boundary_condition.create", "params": {"name": "base", "type": "fixed", "target": {"set": "base"}}},
            {"command": "model.boundary_condition.create", "params": {"name": "push", "type": "displacement", "target": {"set": "top"}, "values": [0, 0, -0.02, null, null, null]}},
            {"command": "model.step.create", "params": {"name": "Squeeze", "procedure": {"type": "nonlinear_static", "nonlinear_geometry": false}, "boundary_conditions": ["base", "push"], "loads": []}}
        ]);
        s.run_script(&script).unwrap();
    }

    #[test]
    fn solve_and_probe_a_contact() {
        let mut s = Session::default();
        pressed_blocks(&mut s);
        let r = s.execute("solve.nonlinear", &json!({"step": "Squeeze"})).unwrap();
        let c = &r["contacts"][0];
        assert_eq!(c["name"], "interface");
        // ν = 0 and held ends: uniaxial strain 1% of E = 1000, over the unit area.
        let force = c["normal_force"].as_f64().unwrap();
        assert!((force - 10.0).abs() < 2e-3, "{r}");
        assert_eq!(c["closed_nodes"], 16, "{r}");
        // A secondary node gives its contact pressure, gap and status; a main one gives none.
        let third = 1.0 / 3.0;
        let probe = s.execute("solve.probe", &json!({"step": "Squeeze", "at": [[third, third, 1.0], [0.5, 0.5, 1.0]]})).unwrap();
        let contact = &probe["nodes"][0]["contact"][0];
        assert_eq!(contact["status"], "closed", "{probe}");
        assert!((contact["pressure"].as_f64().unwrap() - 10.0).abs() < 2e-3, "{probe}");
        assert!(contact["gap"].as_f64().unwrap() < 0.0, "{probe}");
        assert!(probe["nodes"][1].get("contact").is_none(), "{probe}");
        // The contact fields on the results display.
        let v = s.execute("results.contour", &json!({"step": "Squeeze", "field": "contact_pressure"})).unwrap();
        let range = &v["contour"]["range"];
        assert!((range["max"]["value"].as_f64().unwrap() - 10.0).abs() < 2e-3 && range["min"]["value"] == 0.0, "{v}");
        let status = s.execute("results.contour", &json!({"step": "Squeeze", "field": "contact_status"})).unwrap();
        assert_eq!(status["contour"]["range"]["max"]["value"], 1.0, "{status}");
        // The linear analyses refuse the contact by name.
        s.execute(
            "model.step.create",
            &json!({"name": "Linear", "procedure": {"type": "static"}, "boundary_conditions": ["base", "push"], "loads": []}),
        )
        .unwrap();
        let e = s.execute("solve.static", &json!({"step": "Linear"})).unwrap_err().to_string();
        assert!(e.contains("interface"), "{e}");
    }

    #[test]
    fn solve_and_probe_a_frictional_contact() {
        let mut s = Session::default();
        pressed_blocks(&mut s);
        s.execute("interaction.contact.edit", &json!({"name": "interface", "surfaces": {"type": "pair", "main": "lower_top", "secondary": "upper_bottom"}, "friction": {"static": 0.2}})).unwrap();
        s.execute(
            "model.boundary_condition.edit",
            &json!({"name": "push", "type": "displacement", "target": {"set": "top"}, "values": [0.05, 0, -0.02, null, null, null]}),
        )
        .unwrap();
        let r = s.execute("solve.nonlinear", &json!({"step": "Squeeze"})).unwrap();
        let c = &r["contacts"][0];
        assert_eq!(c["slipping_nodes"], 16, "{r}");
        // Slipping, the interface carries μ N sideways.
        let (n, t) = (c["force"][2].as_f64().unwrap(), c["force"][0].as_f64().unwrap());
        assert!((t + 0.2 * n).abs() < 1e-3 * n, "{r}");
        let third = 1.0 / 3.0;
        let probe = s.execute("solve.probe", &json!({"step": "Squeeze", "at": [[third, third, 1.0]]})).unwrap();
        let contact = &probe["nodes"][0]["contact"][0];
        assert_eq!(contact["status"], "slip", "{probe}");
        assert!(contact["shear"][0].as_f64().unwrap() < 0.0 && contact["slip"][0].as_f64().unwrap() > 0.0, "{probe}");
        let v = s.execute("results.contour", &json!({"step": "Squeeze", "field": "contact_status"})).unwrap();
        assert_eq!(v["contour"]["range"]["max"]["value"], 2.0, "{v}");
        let v = s.execute("results.contour", &json!({"step": "Squeeze", "field": "contact_shear"})).unwrap();
        assert!(v["contour"]["range"]["max"]["value"].as_f64().unwrap() > 1.0, "{v}");
    }

    /// A chip with a pad on it, both in a pocket of air (tet4), the chip and pad hex8 at other
    /// sizes: tie `air` holds the chip's and the pad's sides on the air, tie `chip` the pad's
    /// underside on the chip. The pad's bottom edge is on both ties' secondary surfaces, and the
    /// chip's top edge, the main nodes there, follows the air: two constraints per edge node that
    /// nearly repeat each other.
    fn chip_pad_air(s: &mut Session) {
        let (_dir, opened) = super::super::datum::tests::scripted(json!([
            {"command": "solid.box", "params": {"length": 20, "width": 20, "height": 6, "corner": [0, 0, 0], "body_name": "Air"}},
            {"command": "solid.box", "params": {"length": 8, "width": 8, "height": 3, "corner": [6, 6, 0], "operation": "cut"}},
            {"command": "solid.box", "params": {"length": 8, "width": 8, "height": 2, "corner": [6, 6, 0], "operation": "new_body", "body_name": "Chip"}},
            {"command": "solid.box", "params": {"length": 8, "width": 8, "height": 1, "corner": [6, 6, 2], "operation": "new_body", "body_name": "Pad"}},
        ]));
        *s = opened;
        std::mem::forget(_dir);
        s.execute("mesh.tet.generate", &json!({"body": "Air", "size": 2.0})).unwrap();
        s.execute("mesh.hex.generate", &json!({"body": "Chip", "size": 0.8})).unwrap();
        s.execute("mesh.hex.generate", &json!({"body": "Pad", "size": 0.45})).unwrap();
        let m = &s.model;
        // Element faces of a body all of whose nodes satisfy `on`.
        let faces = |body: &str, on: &dyn Fn([f64; 3]) -> bool| -> Vec<Value> {
            m.body_elements(body)
                .into_iter()
                .flat_map(|e| (1..=m.element(e).unwrap().kind.faces().len() as u8).map(move |f| (e, f)))
                .filter(|&(e, f)| {
                    let ns = m.face_nodes(fieldcraft_model::ElementFace { element: e, face: f }).unwrap();
                    // In millimetres, as the design is.
                    ns.iter().all(|n| on(crate::units::to_geometry_point(m.node(*n).unwrap().xyz)))
                })
                .map(|(e, f)| json!({"element": e, "face": f}))
                .collect()
        };
        let near = |a: f64, b: f64| (a - b).abs() < 1e-6;
        let side = |p: [f64; 3]| near(p[0], 6.0) || near(p[0], 14.0) || near(p[1], 6.0) || near(p[1], 14.0);
        let pocket = |p: [f64; 3]| {
            (side(p) && p[2] < 3.0 + 1e-6 && (6.0 - 1e-6..=14.0 + 1e-6).contains(&p[0]) && (6.0 - 1e-6..=14.0 + 1e-6).contains(&p[1]))
                || (near(p[2], 3.0) && (6.0 - 1e-6..=14.0 + 1e-6).contains(&p[0]) && (6.0 - 1e-6..=14.0 + 1e-6).contains(&p[1]))
        };
        let outer = |p: [f64; 3]| near(p[0], 0.0) || near(p[0], 20.0) || near(p[1], 0.0) || near(p[1], 20.0) || near(p[2], 6.0);
        let mut sides = faces("Chip", &|p| side(p));
        sides.extend(faces("Pad", &|p| side(p)));
        let script = json!([
            {"command": "model.set.create", "params": {"name": "all", "elements": m.elements().keys().collect::<Vec<_>>()}},
            {"command": "model.set.create", "params": {"name": "chip", "elements": m.body_elements("Chip")}},
            {"command": "model.set.create", "params": {"name": "pocket", "faces": faces("Air", &pocket)}},
            {"command": "model.set.create", "params": {"name": "sides", "faces": sides}},
            {"command": "model.set.create", "params": {"name": "chip_top", "faces": faces("Chip", &|p| near(p[2], 2.0))}},
            {"command": "model.set.create", "params": {"name": "pad_bottom", "faces": faces("Pad", &|p| near(p[2], 2.0))}},
            {"command": "model.set.create", "params": {"name": "outside", "faces": faces("Air", &outer)}},
            {"command": "model.material.create", "params": {"name": "M", "elastic": {"type": "isotropic", "young": 1000, "poisson": 0.3}, "thermal": {"conductivity": 0.5}}},
            {"command": "model.section.create", "params": {"name": "S", "set": "all", "material": "M", "kind": {"type": "solid"}}},
            {"command": "interaction.tie", "params": {"name": "air", "main": "pocket", "secondary": "sides", "position_tolerance": 0.1}},
            {"command": "interaction.tie", "params": {"name": "chip", "main": "chip_top", "secondary": "pad_bottom", "position_tolerance": 0.1}},
            {"command": "model.load.create", "params": {"name": "Q", "type": "heat_source", "set": "chip", "power": 1}},
            {"command": "model.load.create", "params": {"name": "Room", "type": "convection", "set": "outside", "film": 0.01, "sink": 25}},
            {"command": "model.step.create", "params": {"name": "Steady", "procedure": {"type": "heat_transfer"}, "loads": ["Q", "Room"]}},
        ]);
        s.run_script(&script).unwrap();
    }

    /// #389: a node on the secondary surfaces of two ties follows their main surfaces together,
    /// so the heat balances to round-off and no temperature drops below the surroundings'. (With
    /// a constraint per tie, this model put 335 W more in than the 128 W supplied, and some
    /// temperatures below 0.)
    #[test]
    fn nodes_in_two_ties_keep_the_heat_balance() {
        let mut s = Session::default();
        chip_pad_air(&mut s);
        let r = s.execute("solve.heat", &json!({"step": "Steady"})).unwrap();
        let flows = &r["heat_flows"];
        let supplied = flows["loads"].as_array().unwrap().iter().find(|l| l["name"] == "Q").unwrap()["heat"].as_f64().unwrap();
        assert!((supplied - 128.0).abs() < 1e-9, "{flows}");
        assert!(flows["imbalance"].as_f64().unwrap().abs() < 1e-9 * supplied, "{flows}");
        assert!(r["min_temperature"]["value"].as_f64().unwrap() > 25.0, "{r}");
        let warnings = r["warnings"].to_string();
        assert!(warnings.contains("ties `air`, `chip`"), "{warnings}");
    }

    /// Two ties on one interface in opposite directions hold every node there twice: an error
    /// naming both, not a field.
    #[test]
    fn ties_both_ways_on_one_interface_are_an_error_naming_them() {
        let mut s = Session::default();
        chip_pad_air(&mut s);
        s.execute("interaction.tie", &json!({"name": "back", "main": "pad_bottom", "secondary": "chip_top", "position_tolerance": 0.1})).unwrap();
        let e = s.execute("solve.heat", &json!({"step": "Steady"})).unwrap_err().to_string();
        assert!(e.contains("tie `back`") && e.contains("tie `chip`"), "{e}");
    }
}
