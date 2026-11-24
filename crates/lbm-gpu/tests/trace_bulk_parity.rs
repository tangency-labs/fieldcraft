//! CPU/device parity of the hybrid collision's opt-in trace relaxation (the step's `trace_bulk_one`), on a square duct and on a ball in a box
//! with open faces (a velocity inlet and a pressure outlet), the ball clear of the inlet and near it. Both backends advance 500 steps from
//! rest, 100 at a time; the fields must agree within the suite's stated LES tolerances (|Δu| < 1e-4 lattice, |Δρ| < 1e-5: `check_les` of
//! verification.rs) and, on the beds that carry a trace mode, the relaxation must act: the CPU flow with it on differs from the one with it
//! off by more than 1e-9 (the laminar duct has no trace mode to remove, so there the difference is printed only). Skips without a GPU.
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

fn hybrid(dt: f64, sigma: f64, on: bool) -> Procedure {
    let mut p = collide(STEPS as f64 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { turbulence, sigma: s, trace_bulk_one, .. } = &mut p {
        *turbulence = Some(smagorinsky());
        *s = Some(HrrSigma::Value(sigma));
        *trace_bulk_one = on;
    }
    p
}

fn duct(sigma: f64, on: bool) -> (Setup, Vec<fieldcraft_lbm::Obstacle>) {
    let dx = 1.0 / 8.0;
    let mut d = domain([0.0; 3], [4.0 * dx, 2.0, 2.0], dx);
    d.faces.x_min = FaceCondition::Periodic {};
    d.faces.x_max = FaceCondition::Periodic {};
    d.acceleration = [2.0, 0.0, 0.0];
    let dt = 0.1 * dx / 14.6;
    let m = model(1.0, 1.0 / 151.0, d, hybrid(dt, sigma, on));
    (Setup::new(&m, None, &[], &Limits::default()).expect("the duct's setup"), Vec::new())
}

/// A ball in a box with a velocity inlet and a pressure outlet. `x` puts the ball's centre: 0.9 well clear of the inlet, 0.5 close to it.
fn ball(x: f64, sigma: f64, on: bool) -> (Setup, Vec<fieldcraft_lbm::Obstacle>) {
    let dx = 0.1;
    let dt = 0.03 * dx / 0.1;
    let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.faces.z_max = FaceCondition::Slip {};
    d.obstacles = vec!["Ball".into()];
    let bodies = vec![obstacle("Ball", sphere([x, 0.5, 0.5], 0.2, 3), &["skin"])];
    let m = model(1.0, 6e-3 / 3.0 * dx * dx / dt, d, hybrid(dt, sigma, on));
    (Setup::new(&m, None, &bodies, &Limits::default()).expect("the ball bed's setup"), bodies)
}

type Bed = (Setup, Vec<fieldcraft_lbm::Obstacle>);

fn open_gpu(cpu: &Flow) -> Result<Option<GpuFlow>, String> {
    match GpuFlow::with_options(cpu, Policy::from_env().map_err(|e| e.to_string())?, Options::default()) {
        Ok(g) => Ok(Some(g)),
        Err(e) if !gpu_required() && e.to_string().contains("no GPU device") => {
            println!("skipping: {e}");
            Ok(None)
        }
        Err(e) => Err(e.to_string()),
    }
}

/// The largest |Δu| over the fluid cells between two backends' fields.
fn max_du(fluid: &[bool], a: &mut dyn Backend, b: &mut dyn Backend) -> Result<f64, String> {
    let (fa, fb) = (a.fields().map_err(|e| e.to_string())?, b.fields().map_err(|e| e.to_string())?);
    Ok((0..fa.velocity.len())
        .filter(|&c| fluid[c])
        .map(|c| [0, 1, 2].map(|k| (fa.velocity[c][k] - fb.velocity[c][k]).abs()).into_iter().fold(0.0, f64::max))
        .fold(0.0, f64::max))
}

fn parity(name: &str, make: &dyn Fn(bool) -> Bed, on: bool, must_act: bool) -> Result<(), String> {
    let (setup, obstacles) = make(on);
    let (off_setup, off_obstacles) = make(false);
    let mut cpu = Flow::new(&setup, &obstacles).map_err(|e| e.to_string())?;
    let mut off = Flow::new(&off_setup, &off_obstacles).map_err(|e| e.to_string())?;
    let Some(mut gpu) = open_gpu(&cpu)? else {
        return Ok(());
    };
    let Some(mut gpu_off) = open_gpu(&off)? else {
        return Ok(());
    };
    println!("{name}: device {}, {} cells", gpu.adapter_name(), cpu.grid().len());
    let (mut worst_u, mut worst_r) = (0.0f64, 0.0f64);
    let mut done = 0;
    while done < STEPS {
        for b in [&mut cpu as &mut dyn Backend, &mut gpu, &mut off, &mut gpu_off] {
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
    let fluid: Vec<bool> = (0..cpu.grid().len()).map(|c| Backend::is_fluid(&cpu, c)).collect();
    let rule = max_du(&fluid, &mut cpu, &mut off)?;
    let device_rule = max_du(&fluid, &mut gpu, &mut gpu_off)?;
    println!(
        "{name}: CPU vs device max |Δu| {worst_u:.3e}, max |Δρ| {worst_r:.3e}; treatment on vs off max |Δu| {rule:.3e} (CPU), {device_rule:.3e} (device)"
    );
    if worst_u >= 1e-4 || worst_r >= 1e-5 {
        return Err(format!("{name}: |Δu| {worst_u:.3e} (bar 1e-4), |Δρ| {worst_r:.3e} (bar 1e-5)"));
    }
    if must_act && rule <= 1e-9 {
        return Err(format!("{name}: the treatment did not act on the CPU (on vs off max |Δu| {rule:.3e})"));
    }
    if must_act && device_rule <= 1e-9 {
        return Err(format!("{name}: the treatment did not act on the device (on vs off max |Δu| {device_rule:.3e})"));
    }
    Ok(())
}

#[test]
fn the_cpu_and_the_device_agree_on_the_trace_relaxation() {
    let mut failures = Vec::new();
    let beds: [(&str, &dyn Fn(bool) -> Bed, bool, bool); 6] = [
        ("square duct, 8 cells, sigma 0.99", &|on| duct(0.99, on), true, false),
        ("square duct, 8 cells, sigma 1", &|on| duct(1.0, on), true, false),
        ("ball in a box, open faces, sigma 0.99", &|on| ball(0.9, 0.99, on), true, true),
        ("ball in a box, open faces, sigma 1", &|on| ball(0.9, 1.0, on), true, true),
        ("ball in a box, open faces, sigma 1, control (relaxation off on both)", &|on| ball(0.9, 1.0, on), false, false),
        ("ball near the inlet, sigma 0.99", &|on| ball(0.5, 0.99, on), true, true),
    ];
    for (name, make, on, must_act) in beds {
        if let Err(e) = parity(&format!("{name}, trace relaxation"), make, on, must_act) {
            failures.push(e);
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("; "));
}
