//! The regularized and hybrid collisions rebuild every population from (δρ, j, P); the GPU's
//! rebuild keeps mass and momentum only to its f32 rounding, which beside walls in a nearly
//! steady flow has one sign every step, so the backends' mean momentum drifts apart. The kernels
//! project the rebuilt populations back onto the exact budget (mass unchanged, momentum + F).
//! Gate: a forced plane channel of 16 cells between no-slip walls, 10 000 steps, both layouts,
//! the regularized collision and the hybrid one (σ 0.99) at the smallest relaxation time a bare
//! flow takes and at a moderate one: the GPU's velocity within 1e-5 of the largest CPU speed and
//! its density within 3e-7 of the CPU's, mean momentum within 1e-5 of the CPU's (measured on an
//! RTX 4070 Laptop GPU: 3e-7 to 6.4e-6; without the projection 2.1e-5 to 8.9e-5). Skipped without a
//! GPU adapter unless `FIELDCRAFT_REQUIRE_GPU=1`.
#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure};

/// A plane channel of 16 cells (walls on y, periodic elsewhere), driven by a body force, with
/// relaxation time minus one half `tau_half`.
fn channel(collision: Collision, tau_half: f64, steps: u64) -> (fieldcraft_lbm::Setup, Flow) {
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
    let mut p = collide(steps as f64 * dt, dt, None, collision);
    if let (Procedure::Fluid { sigma, .. }, Collision::Hrr) = (&mut p, collision) {
        *sigma = Some(HrrSigma::Value(0.99));
    }
    prepare(&model(1.0, nu, d, p), &[])
}

fn mean_momentum(f: &fieldcraft_lbm::Fields) -> [f64; 3] {
    let n = f.density.len().max(1) as f64;
    let mut m = [0.0; 3];
    for (r, v) in f.density.iter().zip(&f.velocity) {
        for a in 0..3 {
            m[a] += r * v[a] / n;
        }
    }
    m
}

#[test]
fn the_regularized_and_hybrid_collisions_keep_mass_and_momentum_on_the_gpu() {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let mut failed = Vec::new();
    for collision in [Collision::Regularized, Collision::Hrr] {
        for tau_half in [0.1, 0.005] {
            for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
                let (setup, mut cpu) = channel(collision, tau_half, 10_000);
                let mut gpu = match GpuFlow::with_options(&cpu, policy, Options { streaming, profile: false }) {
                    Ok(g) => g,
                    Err(e) if gpu_required() => panic!("a GPU is required: {e}"),
                    Err(e) if e.to_string().contains("no GPU device") => {
                        println!("skipping: {e}");
                        return;
                    }
                    Err(e) => panic!("{e}"),
                };
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
                let what = format!("{collision:?} τ − ½ {tau_half} {streaming:?}");
                println!("{what}: velocity {du:.2e}, density {drho:.2e}, mean momentum {dm:.2e} ({})", gpu.adapter_name());
                if du > 1e-5 || drho > 3e-7 || dm > 1e-5 {
                    failed.push(format!("{what}: {du:.2e} {drho:.2e} {dm:.2e}"));
                }
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}
