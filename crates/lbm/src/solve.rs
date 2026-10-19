//! Running a fluid step: the time loop, output frames and force histories, in model units.

use fieldcraft_model::Model;

use crate::backend::{Backend, Choice, Fields};
use crate::flow::Flow;
use crate::setup::{Obstacle, Setup};
use crate::steady::{Monitor, SteadyReport, StopReason};
use crate::voxel::Grid;
use crate::{Error, Limits, Result};

/// Most samples a force history keeps; longer runs record every k-th step.
pub const MAX_HISTORY: u64 = 10_000;

/// A relative density variation above this at the end of a run is reported.
pub const COMPRESSIBILITY_WARNING: f64 = 0.01;

/// Bytes in a megabyte, for messages.
const MEGABYTE: f64 = 1e6; // unit-factors: bytes of memory, not a physical quantity

/// The velocity and pressure field at one time, per cell (index `i + nx (j + ny k)`; zero in
/// solid cells). Stored in single precision: output for display and probing.
#[derive(Debug, Clone, PartialEq)]
pub struct Frame {
    pub time: f64,
    pub step: u64,
    pub velocity: Vec<[f32; 3]>,
    /// Gauge pressure.
    pub pressure: Vec<f32>,
    /// Temperature above [`FluidResult::temperature_reference`] (the lattice's scale, see
    /// [`crate::setup::ThermalSetup`]; [`FluidResult::temperatures`] gives kelvin), when the flow
    /// carries heat (empty otherwise).
    pub temperature: Vec<f32>,
    /// The total fluid mass (density times volume).
    pub mass: f64,
}

/// The heat the fluid gives one wall over time (power; positive into the wall).
#[derive(Debug, Clone, PartialEq)]
pub struct HeatHistory {
    pub name: String,
    pub body: Option<String>,
    /// One value per [`FluidResult::history_times`].
    pub power: Vec<f64>,
}

/// A fan's or porous zone's flow and pressure over time (model units).
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneHistory {
    pub name: String,
    /// Whether it is a fan (else a porous zone: a vent, grille, filter or plate).
    pub fan: bool,
    /// One value per [`FluidResult::history_times`]: the volume flow through it (along its axis;
    /// a fan's along its direction).
    pub flow: Vec<f64>,
    /// The static-pressure difference between the fluid cell layers just outside it (so it
    /// includes the wall friction over its thickness plus one cell): the drop across a porous
    /// zone, the rise across a fan.
    pub pressure: Vec<f64>,
    /// A fan's pressure rise from its curve at that time (empty for a porous zone).
    pub fan_rise: Vec<f64>,
}

/// The force the fluid exerts on one wall over time.
#[derive(Debug, Clone, PartialEq)]
pub struct ForceHistory {
    /// `x_min` … for box faces, the face's persistent name for obstacle faces.
    pub name: String,
    /// The obstacle body, for obstacle faces.
    pub body: Option<String>,
    /// One value per [`FluidResult::history_times`].
    pub force: Vec<[f64; 3]>,
}

/// Run statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct FluidStats {
    pub cells: usize,
    pub fluid_cells: usize,
    pub steps: u64,
    pub time_step: f64,
    pub end_time: f64,
    pub tau_plus: f64,
    pub tau_minus: f64,
    pub magic: f64,
    pub collision: fieldcraft_model::Collision,
    /// The characteristic velocity (model units), its lattice value and its lattice Mach number.
    pub velocity: f64,
    pub lattice_velocity: f64,
    pub mach: f64,
    /// Fluid cells times steps.
    pub cell_updates: u64,
    /// What stepped the flow: `cpu`, or `gpu: <adapter>`.
    pub backend: String,
    /// The scalar lattice's relaxation time, with heat.
    pub tau_thermal: Option<f64>,
}

/// The result of a fluid step.
#[derive(Debug, Clone, PartialEq)]
pub struct FluidResult {
    pub step: String,
    pub grid: Grid,
    /// Per cell: not fluid (inside an obstacle).
    pub solid: Vec<bool>,
    /// The initial field and `frames` more, evenly spaced to the end.
    pub frames: Vec<Frame>,
    pub history_times: Vec<f64>,
    /// Every box wall and every obstacle face.
    pub forces: Vec<ForceHistory>,
    /// The heat every wall receives, with heat (empty otherwise), the same walls and times.
    pub heat: Vec<HeatHistory>,
    /// Every fan and porous zone, the same times.
    pub zones: Vec<ZoneHistory>,
    pub stats: FluidStats,
    pub warnings: Vec<String>,
    /// How the run to steady state went, when the step has a stop rule.
    pub steady: Option<SteadyReport>,
    /// What the frames' temperatures are above, in kelvin: the lattice's zero
    /// ([`crate::setup::ThermalSetup::zero`]) with heat, 0 without.
    pub temperature_reference: f64,
}

impl FluidResult {
    /// The cell containing a point and the frame's velocity, pressure and temperature (`None`
    /// without heat) there (`None` outside the grid; zeros in a solid cell).
    pub fn probe(&self, frame: usize, p: [f64; 3]) -> Option<(usize, [f32; 3], f32, Option<f64>)> {
        let cell = self.grid.cell_at(p)?;
        let f = self.frames.get(frame)?;
        let temperature = f.temperature.get(cell).map(|t| fieldcraft_model::units::kelvin_from(f64::from(*t), self.temperature_reference));
        Some((cell, *f.velocity.get(cell)?, *f.pressure.get(cell)?, temperature))
    }

    /// The temperature of every cell in `frame` (empty without heat).
    pub fn temperatures(&self, frame: &Frame) -> Vec<f64> {
        frame.temperature.iter().map(|t| fieldcraft_model::units::kelvin_from(f64::from(*t), self.temperature_reference)).collect()
    }
}

/// Solve the fluid step named `step` (or the model's only fluid step) with the obstacles' surfaces
/// prepared by the caller, on the CPU.
pub fn solve_fluid(model: &Model, step: Option<&str>, obstacles: &[Obstacle], limits: &Limits) -> Result<FluidResult> {
    solve_fluid_on(model, step, obstacles, limits, Choice::Cpu)
}

/// [`solve_fluid`] on the backend `choice` names (the registered device backend, or the CPU).
pub fn solve_fluid_on(model: &Model, step: Option<&str>, obstacles: &[Obstacle], limits: &Limits, choice: Choice) -> Result<FluidResult> {
    let setup = Setup::new(model, step, obstacles, limits)?;
    let flow = Flow::new(&setup, obstacles)?;
    let mut backend = crate::backend::backend(choice, flow, &setup)?;
    run(&setup, backend.as_mut(), limits)
}

fn frame(setup: &Setup, flow: &mut dyn Backend, step: u64) -> Result<Frame> {
    let n = setup.grid.len();
    let mut velocity = Vec::new();
    let mut pressure = Vec::new();
    velocity.try_reserve_exact(n).map_err(|_| Error::Limit("not enough memory for an output frame".into()))?;
    pressure.try_reserve_exact(n).map_err(|_| Error::Limit("not enough memory for an output frame".into()))?;
    let vs = setup.velocity_scale();
    let Fields { density, velocity: u, temperature: t } = flow.fields()?;
    let mut temperature = Vec::new();
    if !t.is_empty() {
        temperature.try_reserve_exact(n).map_err(|_| Error::Limit("not enough memory for an output frame".into()))?;
    }
    for cell in 0..n {
        if flow.is_fluid(cell) {
            velocity.push(u.get(cell).copied().unwrap_or([0.0; 3]).map(|c| (c * vs) as f32));
            pressure.push(setup.gauge_pressure(density.get(cell).copied().unwrap_or(1.0)) as f32);
        } else {
            velocity.push([0.0; 3]);
            pressure.push(0.0);
        }
        if !t.is_empty() {
            let v = t.get(cell).copied().unwrap_or(0.0);
            // Single precision holds about ±3.4e38: a temperature beyond it cannot be stored.
            let stored = v as f32;
            if v.is_finite() && !stored.is_finite() {
                return Err(Error::Limit(format!(
                    "a temperature of {v:.3e} °C at time step {step} is beyond what the results can hold; check the heat fluxes"
                )));
            }
            temperature.push(stored);
        }
    }
    let mass = flow.mass()? * setup.density * setup.grid.spacing.powi(3);
    Ok(Frame { time: step as f64 * setup.time_step, step, velocity, pressure, temperature, mass })
}

/// Run a prepared flow (fresh from [`Flow::new`], possibly with its own initial state, on any
/// backend) through the setup's steps, recording frames and force histories. The steps go to
/// the backend in chunks, up to the next frame, history sample, monitor sample or progress
/// check. With a stop rule ([`crate::steady`]) the run may end early: steady, or diverged
/// (an unstable flow then ends the run with the frames saved so far instead of an error); a
/// frame is added at the stop when the flow is still sound.
pub fn run(setup: &Setup, flow: &mut dyn Backend, limits: &Limits) -> Result<FluidResult> {
    let n = setup.grid.len();
    let total = setup.steps;
    let frames = u64::from(setup.frames).min(total).max(1);
    let frame_bytes = (frames + 1) as f64 * n as f64 * if setup.thermal.is_some() { 20.0 } else { 16.0 };
    if !(frame_bytes <= limits.max_frame_bytes as f64) {
        return Err(Error::Limit(format!(
            "{} frames of {n} cells need {:.0} MB, above the cap of {:.0} MB; write fewer frames or raise `max_frame_bytes`",
            frames + 1,
            frame_bytes / MEGABYTE,
            limits.max_frame_bytes as f64 / MEGABYTE
        )));
    }
    let every = total.div_ceil(MAX_HISTORY).max(1);
    let fs = setup.force_scale();
    let hs = setup.thermal.as_ref().map_or(0.0, |t| t.heat_scale);

    let mut out_frames = vec![frame(setup, flow, 0)?];
    let mut history_times = Vec::new();
    let mut forces: Vec<Vec<[f64; 3]>> = vec![Vec::new(); flow.groups().len()];
    let mut heats: Vec<Vec<f64>> = vec![Vec::new(); if setup.thermal.is_some() { flow.groups().len() } else { 0 }];
    let mut zones: Vec<ZoneHistory> = setup
        .zones
        .iter()
        .map(|z| ZoneHistory { name: z.name.clone(), fan: z.is_fan(), flow: Vec::new(), pressure: Vec::new(), fan_rise: Vec::new() })
        .collect();
    let (qs, ps) = (setup.volume_flow_scale(), setup.pressure_scale());
    let mut next_frame = 1u64;
    // Frame j (1..=frames) at step round(j · total / frames).
    let frame_step = |j: u64| (j * total + frames / 2) / frames;
    fieldcraft_model::progress::stage("Fluid steps", Some(total as f64), "");
    let mut monitor = Monitor::new(setup, flow)?;
    if let Some(m) = monitor.as_mut() {
        m.sample(setup, flow, 0)?;
    }
    let mut sound = true;
    let mut k = 0u64;
    while k < total {
        fieldcraft_model::progress::check()?;
        // Up to the next event: a history sample, a frame, a progress count.
        let mut next = total.min(k + 64 - k % 64);
        next = next.min((k / every + 1) * every);
        if next_frame <= frames {
            next = next.min(frame_step(next_frame).max(k + 1));
        }
        if let Some(m) = &monitor {
            next = next.min((k / m.setup.sample + 1) * m.setup.sample);
        }
        match (flow.advance(next - k), monitor.as_mut()) {
            (Ok(()), _) => {}
            (Err(Error::Unstable(e)), Some(m)) => {
                m.unstable(e);
                sound = false;
                break;
            }
            (Err(e), _) => return Err(e),
        }
        k = next;
        if k.is_multiple_of(64) || k == total {
            fieldcraft_model::progress::count(k as f64);
        }
        if k.is_multiple_of(every) || k == total {
            history_times.push(k as f64 * setup.time_step);
            for (h, f) in forces.iter_mut().zip(flow.forces()?) {
                h.push(f.map(|c| c * fs));
            }
            if !zones.is_empty() {
                for (h, z) in zones.iter_mut().zip(flow.zones()?) {
                    h.flow.push(z.flow * qs);
                    h.pressure.push(z.pressure * ps);
                    if h.fan {
                        h.fan_rise.push(z.fan_rise * ps);
                    }
                }
            }
            if !heats.is_empty() {
                for (h, q) in heats.iter_mut().zip(flow.heat()?) {
                    h.push(q * hs);
                }
            }
        }
        while next_frame <= frames && k == frame_step(next_frame) {
            out_frames.push(frame(setup, flow, k)?);
            next_frame += 1;
        }
        if let Some(m) = monitor.as_mut()
            && m.samples_at(k)
            && let Some(reason) = m.sample(setup, flow, k)?
        {
            sound = reason != StopReason::Diverged;
            if sound && out_frames.last().is_some_and(|f| f.step != k) {
                out_frames.push(frame(setup, flow, k)?);
            }
            break;
        }
    }

    let mut warnings = setup.warnings.clone();
    let steady = monitor.map(Monitor::report);
    if let Some(r) = &steady {
        match r.reason {
            StopReason::Diverged => warnings.push(format!(
                "the run diverged at time {}: {}; the frames end there",
                k as f64 * setup.time_step,
                r.message.as_deref().unwrap_or("a monitored quantity blew up")
            )),
            StopReason::MaxTime => warnings.push(format!(
                "the flow did not reach steady state (relative changes below {:e}) within the step's length; run longer or loosen the tolerance",
                r.tolerance
            )),
            StopReason::Steady => {}
        }
    }
    // The weakly compressible scheme turns pressure differences into density differences of
    // Δp / (ρ c_s²) in lattice units; large ones show the time step is too long for the flow.
    let spread = if sound {
        let end = flow.fields()?;
        (0..n).filter(|c| flow.is_fluid(*c)).filter_map(|c| end.density.get(c)).map(|rho| (rho - 1.0).abs()).fold(0.0, f64::max)
    } else {
        0.0
    };
    if spread > COMPRESSIBILITY_WARNING {
        warnings.push(format!(
            "the density varies by {:.1}% at the end: the weakly compressible error is visible; a shorter time step reduces it (quadratically)",
            100.0 * spread
        ));
    }
    for body in flow.empty_obstacles() {
        warnings.push(format!("obstacle `{body}` covers no cell centre: it is smaller than the spacing or outside the box"));
    }
    if flow.name() != "cpu" {
        warnings
            .push(format!("solved on {}: single precision, deterministic on this device, within the backend's tolerance of the CPU", flow.name()));
    }
    let solid: Vec<bool> = (0..n).map(|c| !flow.is_fluid(c)).collect();
    Ok(FluidResult {
        step: setup.step.clone(),
        grid: setup.grid,
        solid,
        frames: out_frames,
        history_times,
        forces: flow.groups().iter().zip(forces).map(|(g, force)| ForceHistory { name: g.name.clone(), body: g.body.clone(), force }).collect(),
        zones,
        heat: flow.groups().iter().zip(heats).map(|(g, power)| HeatHistory { name: g.name.clone(), body: g.body.clone(), power }).collect(),
        stats: FluidStats {
            cells: n,
            fluid_cells: flow.fluid_cells(),
            steps: k,
            time_step: setup.time_step,
            end_time: k as f64 * setup.time_step,
            tau_plus: setup.tau_plus,
            tau_minus: setup.tau_minus,
            magic: setup.magic,
            collision: setup.collision,
            velocity: setup.velocity,
            lattice_velocity: setup.lattice_velocity,
            mach: setup.mach,
            cell_updates: (flow.fluid_cells() as u64).saturating_mul(k),
            backend: flow.name(),
            tau_thermal: setup.thermal.as_ref().map(|t| t.tau),
        },
        warnings,
        steady,
        temperature_reference: setup.thermal.as_ref().map_or(0.0, |t| t.zero),
    })
}
