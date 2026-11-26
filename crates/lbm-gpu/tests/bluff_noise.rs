//! The density noise of a bluff body's LES against its resolution and its relaxation time (#1110): a 2 × 2 design on the
//! GPU, one cylinder case at Re 3900 at 16 and 32 cells per diameter and at two lattice speeds.
//!
//! The protocol is #1106's: the start at rest on the body (a divergence-free field that stops at the body's surface,
//! `psi = U (y − y_c)(1 − e^(−d/0.5D))`), and the inlet raised from rest over 2.5 fundamental periods of the box,
//! T = m T₀ with T₀ = 4 L / c_s (the manual's `inlet_ramp`, `docs/manual/fluid/boundary-conditions.md`, "Raising the inlet").
//! The walls are halfway, as in the #785 probes. The box is 24 D × 12 D in one cell of span, the collision Hrr with
//! Smagorinsky (C = 0.1), at Re 3900; the run is 120 D/U.
//!
//! The relaxation time is τ⁺ = ½ + 3 ν_lat with ν_lat = ν Δt / Δx² = U_lat · cells per D / Re (D = U = 1), so
//! τ⁺ = ½ + 3 U_lat · cells per D / Re. The two τ⁺ levels of the design are two lattice speeds at each resolution.
//!
//! Every D/U of the history is recorded: the density departure of the fluid cells, its peak |ρ − 1| over the box and the
//! rms of ρ − 1 over the wake (from 1 D behind the body to 2 D before the outlet). The statistics are over the window
//! from `WINDOW_START_D_PER_U` on, and from the end of the ramp on, so the ramp never enters them. A run that goes
//! unstable reports the step it stopped at and the history before it.
//!
//! Heavy: one GPU job per design, run with `--ignored` (`docs/verification.md`, "Bluff-body LES noise (#1110)").

#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use std::sync::Mutex;
use std::time::Instant;

use common::*;
use fieldcraft_gpu::select::Policy;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile, Turbulence, WallTreatment};

/// The box and the flow in D = 1 units: 24 D long, 12 D high, the cylinder's centre 6 D from the inlet.
const LENGTH_D: f64 = 24.0;
const HEIGHT_D: f64 = 12.0;
const UP_D: f64 = 6.0;
const RE: f64 = 3900.0;
/// The Smagorinsky constant of the LES study's cylinder case.
const SMAGORINSKY: f64 = 0.1;
/// The lattice sound speed, 1/√3.
const SOUND_SPEED: f64 = 0.577_350_269_189_625_8;
/// The inlet's ramp in fundamental periods of the box (#1106's default, the manual's lowest-residual row).
const RAMP_PERIODS: f64 = 2.5;
/// The run in D/U, and the first D/U of the statistics window (after the ramp, which ends near 21 D/U at 0.05).
const END_D_PER_U: u64 = 120;
const WINDOW_START_D_PER_U: f64 = 60.0;
/// The wake is sampled from this many D behind the body's centre to 2 D before the outlet.
const WAKE_FROM_D: f64 = 1.0;
const WAKE_TO_OUTLET_D: f64 = 2.0;
/// The design: (cells per diameter, lattice speed).
const DESIGN: [(f64, f64); 4] = [(16.0, 0.05), (16.0, 0.10), (32.0, 0.05), (32.0, 0.10)];

/// The GPU is one device at a time in a process (#614).
static TURN: Mutex<()> = Mutex::new(());

/// The divergence-free start at rest on the body, as #1106 sets it: the velocity in D/U from the stream function, with
/// its derivatives taken over a quarter cell.
fn start_velocity(x: f64, y: f64, cells_per_d: f64) -> [f64; 2] {
    let (xc, yc) = (UP_D, HEIGHT_D / 2.0);
    let dist = |x: f64, y: f64| (((x - xc).powi(2) + (y - yc).powi(2)).sqrt() - 0.5).max(0.0);
    let psi = |x: f64, y: f64| (y - yc) * (1.0 - (-dist(x, y) / 0.5).exp());
    let h = 0.25 / cells_per_d;
    let u = (psi(x, y + h) - psi(x, y - h)) / (2.0 * h);
    let v = -(psi(x + h, y) - psi(x - h, y)) / (2.0 * h);
    [u, v]
}

/// One D/U of the history: the density departure of the fluid cells.
struct Sample {
    d_per_u: f64,
    peak: f64,
    wake_rms: f64,
}

/// What one design point produced.
struct Outcome {
    mach: f64,
    tau_plus: f64,
    dt: f64,
    ramp_end: f64,
    samples: Vec<Sample>,
    diverged: Option<String>,
}

/// Run one design point on the GPU.
fn run(cells_per_d: f64, speed: f64) -> Outcome {
    let dx = 1.0 / cells_per_d;
    let span = dx;
    let (xc, yc) = (UP_D, HEIGHT_D / 2.0);
    let dt = speed * dx;
    let steps_per_d_per_u = (cells_per_d / speed).round() as u64;
    // T₀ = 4 L / c_s with L the box's length in cells; the ramp is m T₀ steps, and its end in D/U is its time (D = U = 1).
    let length_cells = (LENGTH_D * cells_per_d).round();
    let ramp_time = RAMP_PERIODS * 4.0 * length_cells / SOUND_SPEED * dt;
    let mut d = domain([0.0; 3], [LENGTH_D, HEIGHT_D, span], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective, sigma: None };
    d.faces.y_min = FaceCondition::Slip {};
    d.faces.y_max = FaceCondition::Slip {};
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    d.obstacles = vec!["Body".into()];
    let body = obstacle("Body", cylinder([xc, yc], 0.5, -dx, span + dx, 1024), &["side", "bottom", "top"]);
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(END_D_PER_U as f64),
        steps: None,
        time_step: Some(dt),
        frames: 1,
        magic: None,
        sigma: None,
        collision: Some(Collision::Hrr),
        turbulence: Some(Turbulence::Smagorinsky { constant: SMAGORINSKY, prandtl: 0.9, van_driest: None, wall_model: None }),
        steady: None,
        walls: WallTreatment::Halfway,
        inlet_ramp: Some(ramp_time),
    };
    let m = model(1.0, 1.0 / RE, d, procedure);
    let obstacles = [body];
    let (setup, mut flow) = prepare(&m, &obstacles);
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| {
        let [x, y, _] = grid.centre(cell);
        let [u, v] = start_velocity(x, y, cells_per_d);
        (1.0, [speed * u, speed * v, 0.0])
    });
    let fluid: Vec<usize> = (0..grid.len()).filter(|&c| flow.is_fluid(c)).collect();
    let centres: Vec<[f64; 3]> = fluid.iter().map(|&c| grid.centre(c)).collect();
    let wake: Vec<bool> = centres.iter().map(|p| p[0] >= xc + WAKE_FROM_D && p[0] <= LENGTH_D - WAKE_TO_OUTLET_D).collect();
    let tau_plus = 0.5 + 3.0 * speed * cells_per_d / RE;

    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let _turn = TURN.lock().unwrap_or_else(|e| e.into_inner());
    let mut gpu = GpuFlow::with_options(&flow, policy, Options { streaming: Streaming::InPlace, profile: false })
        .expect("a GPU is required for the bluff-body study");
    println!(
        "CASE cells_per_d {cells_per_d} U_lat {speed}: GPU {} | τ⁺ {tau_plus:.5} | Mach {:.3} | cells {} | steps {} | start: stream-function at rest on the body | inlet ramp {RAMP_PERIODS} T₀ = {ramp_time:.2} D/U | halfway walls",
        gpu.adapter_name(),
        speed * 3f64.sqrt(),
        grid.len(),
        setup.steps
    );

    let mut samples = Vec::new();
    for d_per_u in 1..=END_D_PER_U {
        if let Err(e) = gpu.advance(steps_per_d_per_u) {
            println!("DIVERGED in D/U block {d_per_u}: {e}");
            return Outcome { mach: speed * 3f64.sqrt(), tau_plus, dt, ramp_end: ramp_time, samples, diverged: Some(e.to_string()) };
        }
        let fields = gpu.fields().unwrap();
        let mut peak: f64 = 0.0;
        let (mut sum2, mut n) = (0.0, 0usize);
        for (i, &c) in fluid.iter().enumerate() {
            let x = fields.density[c] - 1.0;
            peak = peak.max(x.abs());
            if wake[i] {
                sum2 += x * x;
                n += 1;
            }
        }
        let wake_rms = (sum2 / n as f64).sqrt();
        println!("  t = {d_per_u} D/U: peak |ρ−1| {peak:.4e}, wake rms {wake_rms:.4e}");
        samples.push(Sample { d_per_u: d_per_u as f64, peak, wake_rms });
    }
    Outcome { mach: speed * 3f64.sqrt(), tau_plus, dt, ramp_end: ramp_time, samples, diverged: None }
}

/// The mean and the largest value of a series (NaN for an empty one).
fn summary(values: impl Iterator<Item = f64>) -> (f64, f64) {
    let (mut sum, mut n, mut max) = (0.0, 0usize, 0.0f64);
    for v in values {
        sum += v;
        n += 1;
        max = max.max(v);
    }
    (sum / n as f64, max)
}

#[test]
#[ignore = "heavy: GPU study, cargo test -- --ignored (docs/verification.md, \"Bluff-body LES noise (#1110)\")"]
fn bluff_body_density_noise_two_by_two() {
    let start = Instant::now();
    for (cells_per_d, speed) in DESIGN {
        let t0 = Instant::now();
        let o = run(cells_per_d, speed);
        let window_start = WINDOW_START_D_PER_U.max(o.ramp_end);
        let in_window = o.samples.iter().filter(|s| s.d_per_u >= window_start);
        let (peak_mean, peak_max) = summary(in_window.clone().map(|s| s.peak));
        let (rms_mean, rms_max) = summary(in_window.map(|s| s.wake_rms));
        let last = o.samples.last().map_or(f64::NAN, |s| s.d_per_u);
        println!(
            "ROW cells_per_d {cells_per_d} U_lat {speed} Mach {:.3} τ⁺ {:.5} dt {:.4e} ramp_end {:.2} D/U | window {window_start:.1}–{END_D_PER_U} D/U | peak |ρ−1| mean {:.2} % max {:.2} % | wake rms mean {:.2} % max {:.2} % | last D/U recorded {last} | diverged {} | wall {:.1} s",
            o.mach,
            o.tau_plus,
            o.dt,
            o.ramp_end,
            100.0 * peak_mean,
            100.0 * peak_max,
            100.0 * rms_mean,
            100.0 * rms_max,
            o.diverged.is_some(),
            t0.elapsed().as_secs_f64()
        );
    }
    println!("design wall time {:.1} s", start.elapsed().as_secs_f64());
}
