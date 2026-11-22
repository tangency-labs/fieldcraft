//! CPU/device parity of the hybrid collision's opt-in wall treatment: σ = 1 in the walls' layer (`Flow::set_wall_layer_sigma_one`) and the
//! trace relaxation at ω_b = 1 (`Flow::set_trace_bulk_one`), alone and together, on a bed whose layer is a box face (a square duct) and on one
//! with a wall of obstacle cells and open faces (a ball in a box with an inlet and a pressure outlet). Both backends advance 500 steps from
//! rest, 100 at a time; the fields must agree within the bar of the backend (|Δu| < 2e-5 lattice, |Δρ| < 1e-6) and the treatment must act: the
//! CPU flow with it on differs from the one with it off by more than that. Skips without a GPU.
#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use common::*;
use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_lbm_gpu::{GpuFlow, Options};
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Outflow, Procedure, Profile, Turbulence};

const STEPS: u64 = 500;
const CHUNK: u64 = 100;

fn smagorinsky() -> Turbulence {
    Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None }
}

fn hybrid(dt: f64) -> Procedure {
    let mut p = collide(STEPS as f64 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { turbulence, sigma, .. } = &mut p {
        *turbulence = Some(smagorinsky());
        *sigma = Some(HrrSigma::Value(0.99));
    }
    p
}

fn duct() -> (Setup, Vec<fieldcraft_lbm::Obstacle>) {
    let dx = 1.0 / 8.0;
    let mut d = domain([0.0; 3], [4.0 * dx, 2.0, 2.0], dx);
    d.faces.x_min = FaceCondition::Periodic {};
    d.faces.x_max = FaceCondition::Periodic {};
    d.acceleration = [2.0, 0.0, 0.0];
    let dt = 0.1 * dx / 14.6;
    let m = model(1.0, 1.0 / 151.0, d, hybrid(dt));
    (Setup::new(&m, None, &[], &Limits::default()).expect("the duct's setup"), Vec::new())
}

/// A ball in a box with a velocity inlet and a pressure outlet. `x` puts the ball: 0.9 clear of the margin at the inlet, 0.5 with part of
/// its layer inside it, where it stays.
fn ball(x: f64) -> (Setup, Vec<fieldcraft_lbm::Obstacle>) {
    let dx = 0.1;
    let dt = 0.03 * dx / 0.1;
    let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_max = FaceCondition::Slip {};
    d.obstacles = vec!["Ball".into()];
    let bodies = vec![obstacle("Ball", sphere([x, 0.5, 0.5], 0.2, 3), &["skin"])];
    let m = model(1.0, 6e-3 / 3.0 * dx * dx / dt, d, hybrid(dt));
    (Setup::new(&m, None, &bodies, &Limits::default()).expect("the ball bed's setup"), bodies)
}

fn parity(name: &str, setup: &Setup, obstacles: &[fieldcraft_lbm::Obstacle], layer: bool, trace: bool) -> Result<(), String> {
    let mut cpu = Flow::new(setup, obstacles).map_err(|e| e.to_string())?;
    cpu.set_wall_layer_sigma_one(layer);
    cpu.set_trace_bulk_one(trace);
    let mut off = Flow::new(setup, obstacles).map_err(|e| e.to_string())?;
    let mut gpu = match GpuFlow::with_options(&cpu, Policy::from_env().map_err(|e| e.to_string())?, Options::default()) {
        Ok(g) => g,
        Err(e) if !gpu_required() && e.to_string().contains("no GPU device") => {
            println!("skipping: {e}");
            return Ok(());
        }
        Err(e) => return Err(e.to_string()),
    };
    let in_layer = cpu.wall_sigma_layer_mask().iter().filter(|&&m| m).count();
    println!("{name}: device {}, {} cells, {in_layer} in the layer, layer {layer}, trace {trace}", gpu.adapter_name(), cpu.grid().len());
    assert!(in_layer > 0, "{name}: the layer is empty");
    let (mut worst_u, mut worst_r) = (0.0f64, 0.0f64);
    let mut done = 0;
    while done < STEPS {
        for b in [&mut cpu as &mut dyn Backend, &mut gpu, &mut off] {
            b.advance(CHUNK).map_err(|e| e.to_string())?;
        }
        done += CHUNK;
        let (a, b) = (Backend::fields(&mut cpu).map_err(|e| e.to_string())?, Backend::fields(&mut gpu).map_err(|e| e.to_string())?);
        for c in (0..a.velocity.len()).filter(|&c| Backend::is_fluid(&cpu, c)) {
            let e = [0, 1, 2].map(|k| a.velocity[c][k] - b.velocity[c][k]);
            worst_u = worst_u.max(e.iter().map(|v| v * v).sum::<f64>().sqrt());
            worst_r = worst_r.max((a.density[c] - b.density[c]).abs());
        }
    }
    let (a, o) = (Backend::fields(&mut cpu).map_err(|e| e.to_string())?, Backend::fields(&mut off).map_err(|e| e.to_string())?);
    let rule = (0..a.velocity.len())
        .filter(|&c| Backend::is_fluid(&cpu, c))
        .map(|c| [0, 1, 2].map(|k| (a.velocity[c][k] - o.velocity[c][k]).abs()).into_iter().fold(0.0, f64::max))
        .fold(0.0, f64::max);
    println!("{name}: CPU vs device max |Δu| {worst_u:.3e}, max |Δρ| {worst_r:.3e}; treatment on vs off (CPU) max |Δu| {rule:.3e}");
    if worst_u >= 2e-5 || worst_r >= 1e-6 {
        return Err(format!("{name}: |Δu| {worst_u:.3e} (bar 2e-5), |Δρ| {worst_r:.3e} (bar 1e-6)"));
    }
    if rule <= 1e-9 {
        return Err(format!("{name}: the treatment did not act (on vs off max |Δu| {rule:.3e})"));
    }
    Ok(())
}

#[test]
fn the_cpu_and_the_device_agree_on_the_wall_treatment() {
    let mut failures = Vec::new();
    for (layer, trace) in [(true, false), (false, true), (true, true)] {
        for (name, (setup, obstacles)) in
            [("square duct, 8 cells", duct()), ("ball in a box, open faces", ball(0.9)), ("ball near the inlet", ball(0.5))]
        {
            if let Err(e) = parity(&format!("{name}, layer {layer}, trace {trace}"), &setup, &obstacles, layer, trace) {
                failures.push(e);
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("; "));
}
