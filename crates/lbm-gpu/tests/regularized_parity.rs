//! The regularized and hybrid collisions rebuild every population from (δρ, j, P); the GPU's
//! rebuild keeps mass and momentum only to its f32 rounding, so the kernels project the rebuilt
//! populations back onto the exact budget (mass unchanged, momentum + F).
//!
//! Gates, all on the device against the CPU:
//! - One collision step of a plane channel (16 cells, walls on y, periodic elsewhere, starting
//!   at lattice speed 0.05 so the velocity is large against its rounding) for the regularized
//!   collision and the hybrid one (σ 0.99), at two relaxation times. The density and velocity of
//!   every cell match the CPU within a bar derived from f32 rounding (see `ROUNDINGS` and
//!   `density_bar`).
//! - Conservation of mass and momentum by one step of a fully periodic box (64 cells, uniform
//!   flow, no force), where streaming is an exact permutation: the only change in Σf and Σf c is
//!   the collision's, and it must stay within the f32 summation bound.
//! - The 10 000-step channel is printed, not asserted: its numbers are the drift the projection
//!   removes, the figures the kernel's comment quotes.
//!
//! Skipped without a GPU adapter unless `FIELDCRAFT_REQUIRE_GPU=1`.
#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::lattice::Q;
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure};

/// The unit roundoff of f32, 2^-24: half the spacing of the f32 numbers at 1.
const UNIT_ROUNDOFF: f64 = f32::EPSILON as f64 / 2.0;

/// The roundings one population takes in one collision step on the device, each relative to
/// its own magnitude: the f32 storage of the input population, the rebuild, the projection's
/// add and the storage of the output.
const ROUNDINGS: f64 = 4.0;

/// The lattice speed the single-step channel starts at (a Mach number of 0.05·√3, about 0.087).
const CHANNEL_SPEED: f64 = 0.05;

/// The lattice velocity of the periodic box's uniform start.
const BOX_VELOCITY: [f64; 3] = [0.05, 0.02, 0.01];

/// A plane channel of 16 cells (walls on y, periodic elsewhere), driven by a body force, with
/// relaxation time minus one half `tau_half`, run for `steps` steps from lattice speed
/// `start_speed` along x.
fn channel(collision: Collision, tau_half: f64, steps: u64, start_speed: f64) -> (fieldcraft_lbm::Setup, Flow) {
    let (h, n) = (1.0, 16usize);
    let dx = h / n as f64;
    let dt = dx * dx;
    // Lattice ν = (τ − ½)/3 in cells and steps, so ν = tau_half/3 · dx²/dt = tau_half/3.
    let nu = tau_half / 3.0;
    let umax = 0.05 * dx / dt;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.acceleration = [8.0 * nu * umax / (h * h), 0.0, 0.0];
    d.initial_velocity = [start_speed * dx / dt, 0.0, 0.0];
    prepare(&model(1.0, nu, d, procedure(collision, steps as f64 * dt, dt)), &[])
}

/// A fully periodic 4 × 4 × 4 box at lattice velocity [`BOX_VELOCITY`], no force.
fn periodic_box(collision: Collision, tau_half: f64) -> (fieldcraft_lbm::Setup, Flow) {
    let dx = 1.0 / 16.0;
    let dt = dx * dx;
    let mut d = domain([0.0; 3], [4.0 * dx; 3], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.initial_velocity = BOX_VELOCITY.map(|c| c * dx / dt);
    prepare(&model(1.0, tau_half / 3.0, d, procedure(collision, dt, dt)), &[])
}

/// The fluid step: the collision, with σ 0.99 for the hybrid one.
fn procedure(collision: Collision, end_time: f64, dt: f64) -> Procedure {
    let mut p = collide(end_time, dt, None, collision);
    if let (Procedure::Fluid { sigma, .. }, Collision::Hrr) = (&mut p, collision) {
        *sigma = Some(HrrSigma::Value(0.99));
    }
    p
}

/// The device flow of `cpu`, or `None` when no GPU is present and none is required.
fn open_gpu(cpu: &Flow, policy: Policy, streaming: Streaming) -> Option<GpuFlow> {
    match GpuFlow::with_options(cpu, policy, Options { streaming, profile: false }) {
        Ok(g) => Some(g),
        Err(e) if gpu_required() => panic!("a GPU is required: {e}"),
        Err(e) if e.to_string().contains("no GPU device") => {
            println!("skipping: {e}");
            None
        }
        Err(e) => panic!("{e}"),
    }
}

/// Σ over the populations of every cell of |f|, one value per cell: the magnitudes the rounding
/// bars scale with.
fn magnitudes(cpu: &Flow, cells: usize) -> Vec<f64> {
    (0..cells).map(|c| (0..Q).map(|q| cpu.population(q, c).abs()).sum()).collect()
}

/// The bar on a cell's density after one collision step, in lattice units. The density is a sum
/// of Q populations, which carries (Q − 1) roundings of the sum, and each population carries
/// [`ROUNDINGS`] roundings of its own magnitude: |Δρ| ≤ (Q + ROUNDINGS) · u · Σ|f|.
fn density_bar(magnitude: f64) -> f64 {
    (Q as f64 + ROUNDINGS) * UNIT_ROUNDOFF * magnitude
}

/// The bar on a cell's velocity after one collision step: the momentum error Σ δf c is bounded
/// by the same sum (|c_a| ≤ 1), and u = (j + F/2)/ρ carries δu ≤ (δj + |u| δρ)/ρ.
fn velocity_bar(magnitude: f64, density: f64, velocity: [f64; 3]) -> f64 {
    let speed = velocity.iter().fold(0.0, |m: f64, v| m.max(v.abs()));
    density_bar(magnitude) * (1.0 + speed) / density
}

/// The bar on a sum Σf or Σf c over `terms` values of total magnitude `magnitude`, taken before
/// and after a step. Each sum carries at most (terms − 1) roundings of its own. The momentum is
/// recovered by the device as ρ·u per cell, two more roundings of that cell's moment, which over
/// all cells sum to at most 2 · u · Σ|f|. The two sums' errors add: 2 · (terms + 2) · u · Σ|f|.
fn conservation_bar(terms: usize, magnitude: f64) -> f64 {
    2.0 * (terms + 2) as f64 * UNIT_ROUNDOFF * magnitude
}

#[test]
fn one_collision_step_matches_the_cpu_within_the_f32_bar() {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let mut failed = Vec::new();
    for collision in [Collision::Regularized, Collision::Hrr] {
        for tau_half in [0.1, 0.005] {
            for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
                let (setup, mut cpu) = channel(collision, tau_half, 1, CHANNEL_SPEED);
                let Some(mut gpu) = open_gpu(&cpu, policy, streaming) else { return };
                let cells = cpu.fields().unwrap().density.len();
                let magnitude = magnitudes(&cpu, cells);
                cpu.advance(setup.steps).unwrap();
                gpu.advance(setup.steps).unwrap();
                let (a, b) = (cpu.fields().unwrap(), gpu.fields().unwrap());
                let what = format!("{collision:?} τ − ½ {tau_half} {streaming:?}");
                let (mut worst_rho, mut worst_u, mut bar_rho, mut bar_u) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                let mut finite = true;
                for (c, &mag) in magnitude.iter().enumerate() {
                    finite &= b.density[c].is_finite() && b.velocity[c].iter().all(|v| v.is_finite());
                    if !finite {
                        break;
                    }
                    let (rb, ub) = (density_bar(mag), velocity_bar(mag, a.density[c], a.velocity[c]));
                    let drho = (a.density[c] - b.density[c]).abs();
                    let du = (0..3).map(|k| (a.velocity[c][k] - b.velocity[c][k]).abs()).fold(0.0, f64::max);
                    worst_rho = worst_rho.max(drho / rb);
                    worst_u = worst_u.max(du / ub);
                    bar_rho = bar_rho.max(rb);
                    bar_u = bar_u.max(ub);
                }
                let (drho, du) = (
                    a.density.iter().zip(&b.density).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max),
                    a.velocity.iter().zip(&b.velocity).map(|(x, y)| (0..3).map(|k| (x[k] - y[k]).abs()).fold(0.0, f64::max)).fold(0.0, f64::max),
                );
                println!("{what}: density {drho:.2e} (bar {bar_rho:.2e}), velocity {du:.2e} (bar {bar_u:.2e}) ({})", gpu.adapter_name());
                if !finite || worst_rho > 1.0 || worst_u > 1.0 {
                    failed.push(format!("{what}: density {drho:.2e} vs {bar_rho:.2e}, velocity {du:.2e} vs {bar_u:.2e}, finite {finite}"));
                }
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}

#[test]
fn the_projection_conserves_mass_and_momentum_to_f32_rounding() {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let mut failed = Vec::new();
    for collision in [Collision::Regularized, Collision::Hrr] {
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let (setup, mut cpu) = periodic_box(collision, 0.1);
            let Some(mut gpu) = open_gpu(&cpu, policy, streaming) else { return };
            let cells = cpu.fields().unwrap().density.len();
            let magnitude: f64 = magnitudes(&cpu, cells).iter().sum();
            let terms = Q * cells;
            let bar = conservation_bar(terms, magnitude);
            let (m0, p0) = (gpu.mass().unwrap(), gpu.population_momentum().unwrap());
            gpu.advance(setup.steps).unwrap();
            let (m1, p1) = (gpu.mass().unwrap(), gpu.population_momentum().unwrap());
            let dm = (m1 - m0).abs();
            let dp = (0..3).map(|k| (p1[k] - p0[k]).abs()).fold(0.0, f64::max);
            let what = format!("{collision:?} {streaming:?}");
            println!("{what}: mass {dm:.2e}, momentum {dp:.2e} (bar {bar:.2e}) ({})", gpu.adapter_name());
            if !(dm <= bar && dp <= bar && m0.is_finite() && m1.is_finite()) {
                failed.push(format!("{what}: mass {dm:.2e}, momentum {dp:.2e} vs bar {bar:.2e}"));
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}

/// Mean momentum of a flow's fields, Σ ρ u / cells.
fn mean_momentum(f: &fieldcraft_lbm::backend::Fields) -> [f64; 3] {
    let n = f.density.len().max(1) as f64;
    let mut m = [0.0; 3];
    for (r, v) in f.density.iter().zip(&f.velocity) {
        for a in 0..3 {
            m[a] += r * v[a] / n;
        }
    }
    m
}

/// The 10 000-step channel, printed: the drift the projection keeps small over a long run.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_ten_thousand_step_channel_drift_is_printed() {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    for collision in [Collision::Regularized, Collision::Hrr] {
        for tau_half in [0.1, 0.005] {
            for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
                let (setup, mut cpu) = channel(collision, tau_half, 10_000, 0.0);
                let Some(mut gpu) = open_gpu(&cpu, policy, streaming) else { return };
                cpu.advance(setup.steps).unwrap();
                gpu.advance(setup.steps).unwrap();
                let (a, b) = (cpu.fields().unwrap(), gpu.fields().unwrap());
                let scale = a.velocity.iter().map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).fold(0.0, f64::max).max(1e-30);
                let du = a
                    .velocity
                    .iter()
                    .zip(&b.velocity)
                    .map(|(x, y)| ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt())
                    .fold(0.0, f64::max)
                    / scale;
                let drho = a.density.iter().zip(&b.density).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
                let (ma, mb) = (mean_momentum(&a), mean_momentum(&b));
                let dm = (ma[0] - mb[0]).abs() / ma[0].abs().max(1e-30);
                println!(
                    "{collision:?} τ − ½ {tau_half} {streaming:?}: velocity {du:.2e}, density {drho:.2e}, mean momentum {dm:.2e} ({})",
                    gpu.adapter_name()
                );
            }
        }
    }
}
