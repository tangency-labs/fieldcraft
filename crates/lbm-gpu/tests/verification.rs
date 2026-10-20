//! The GPU backend against the CPU scheme and the analytic cases of `fieldcraft-lbm`'s
//! verification suite. The stated tolerances of the backend, with populations stored as
//! deviations from their weights in single precision (measured on lavapipe at 2–5e-6, with
//! a margin): velocity fields within 2e-5 of the largest speed, densities within 1e-6, wall
//! forces within 5e-5 of their magnitude (a wall's shear within 5e-5 of itself), mass drift below
//! 1e-6; with heat, temperatures within 2e-5 of the field's range and wall heats within 1e-4 of
//! the largest (a heat is a difference of populations, so it keeps fewer digits). Skipped without
//! a GPU adapter unless `FIELDCRAFT_REQUIRE_GPU=1`, which makes a missing adapter a failure;
//! `FIELDCRAFT_GPU` picks the adapter (`hardware` refuses a software one). The tests take the
//! device one at a time ([`Gpu`]): several devices open at once in one process hung NVIDIA's
//! driver on a GB10 (#614).

#[path = "../../lbm/tests/common/mod.rs"]
mod common;

use std::ops::{Deref, DerefMut};
use std::sync::{Mutex, MutexGuard};

use common::*;
use fieldcraft_gpu::select::{Policy, gpu_required};
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::{Flow, Setup};
use fieldcraft_lbm_gpu::{GpuFlow, Options, Streaming};
use fieldcraft_model::{
    BoxFace, Directional, FaceCondition, Fluid, FluidDomain, FluidThermal, Material, Model, Procedure, Profile, Step, Thermal, ThermalCondition,
};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// Whose turn it is on the GPU.
static TURN: Mutex<()> = Mutex::new(());

/// The next turn on the GPU (a test that failed while holding it leaves it to the next).
fn turn() -> MutexGuard<'static, ()> {
    TURN.lock().unwrap_or_else(|e| e.into_inner())
}

/// A flow on the GPU, holding the turn until it is dropped (the flow first, then the turn).
struct Gpu {
    flow: GpuFlow,
    _turn: MutexGuard<'static, ()>,
}

impl Deref for Gpu {
    type Target = GpuFlow;
    fn deref(&self) -> &GpuFlow {
        &self.flow
    }
}

impl DerefMut for Gpu {
    fn deref_mut(&mut self) -> &mut GpuFlow {
        &mut self.flow
    }
}

/// The flow on the GPU `FIELDCRAFT_GPU` picks, or `None` to skip (a failure when a GPU is
/// required).
fn gpu(flow: &Flow, streaming: Streaming) -> Option<Gpu> {
    let policy = Policy::from_env().expect("FIELDCRAFT_GPU names a policy");
    let turn = turn();
    match GpuFlow::with_options(flow, policy, Options { streaming, profile: false }) {
        Ok(g) => {
            println!("GPU: {} ({streaming:?})", g.adapter_name());
            Some(Gpu { flow: g, _turn: turn })
        }
        Err(e) if gpu_required() => panic!("a GPU is required: {e}"),
        // Only a missing device skips: a shader that does not compile is a failure anywhere.
        Err(e) if e.to_string().contains("no GPU device") => {
            println!("skipping: {e}");
            None
        }
        Err(e) => panic!("{e}"),
    }
}

/// Both backends through every step; the CPU's fields and the GPU's.
fn both(setup: &Setup, mut cpu: Flow, streaming: Streaming) -> Option<(Flow, Gpu)> {
    let mut gpu = gpu(&cpu, streaming)?;
    cpu.advance(setup.steps).unwrap();
    gpu.advance(setup.steps).unwrap();
    Some((cpu, gpu))
}

/// The largest difference of the velocity fields over the largest CPU speed, and of the density
/// fields.
fn field_difference(cpu: &mut Flow, gpu: &mut GpuFlow) -> (f64, f64) {
    let (a, b) = (cpu.fields().unwrap(), gpu.fields().unwrap());
    let scale = a.velocity.iter().map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).fold(0.0, f64::max).max(1e-30);
    let du = a
        .velocity
        .iter()
        .zip(&b.velocity)
        .map(|(x, y)| ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt())
        .fold(0.0, f64::max);
    let drho = a.density.iter().zip(&b.density).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
    (du / scale, drho)
}

#[test]
fn plane_poiseuille_matches_the_analytic_profile_and_the_cpu() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (h, nu, umax, n) = (1.0, 1.0 / 6.0, 0.1, 16usize);
        let dx = h / n as f64;
        let a = 8.0 * nu * umax / (h * h);
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        d.faces.x_min = PERIODIC;
        d.faces.x_max = PERIODIC;
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        d.acceleration = [a, 0.0, 0.0];
        let m = model(1.0, nu, d, step(4.0 * h * h / nu, dx * dx, Some(3.0 / 16.0)));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let f = gpu.fields().unwrap();
        let mut err: f64 = 0.0;
        for j in 0..n {
            let y = (j as f64 + 0.5) * dx;
            let exact = a / (2.0 * nu) * y * (h - y);
            err = err.max((f.velocity[j][0] * setup.velocity_scale() - exact).abs() / umax);
        }
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!("Poiseuille on the GPU: max error {err:.2e} of the peak; CPU/GPU velocity {du:.2e}, density {drho:.2e}");
        assert!(err < 2e-5, "{err}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

#[test]
fn couette_wall_forces_match_the_cpu() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (h, nu, uw, n) = (1.0, 0.1, 0.05, 16usize);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        d.faces.x_min = PERIODIC;
        d.faces.x_max = PERIODIC;
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        d.faces.y_max = FaceCondition::MovingWall { velocity: [uw, 0.0, 0.0] };
        let m = model(1.0, nu, d, step(3.0 * h * h / nu, 0.5 * dx * dx / (3.0 * nu), None));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let bottom = BoxFace::ALL.iter().position(|f| f.name() == "y_min").unwrap();
        let top = BoxFace::ALL.iter().position(|f| f.name() == "y_max").unwrap();
        let tau = nu * uw / h / setup.force_scale() * dx * dx;
        println!("Couette: wall force CPU {:?} / GPU {:?} (expected ±{tau:.3e} per step)", fc[bottom], fg[bottom]);
        // The tangential force is the shear: a cancellation of the links' exchanges to 1e-6 per link,
        // summed from the stored deviations alone (adding the weights back in single precision left
        // it 1e-4 of itself on an RTX 4090; 3e-6 on lavapipe this way). The normal force is the
        // pressure.
        for g in [bottom, top] {
            println!("Couette group {g}: CPU {:.9e}, GPU {:.9e}, shear parity {:.9e}", fc[g][0], fg[g][0], (fc[g][0] - fg[g][0]).abs() / tau.abs());
            assert!((fc[g][0] - fg[g][0]).abs() <= 5e-5 * tau.abs(), "group {g} shear: {} vs {}", fc[g][0], fg[g][0]);
            assert!((fc[g][1] - fg[g][1]).abs() <= 1e-6 * fc[g][1].abs(), "group {g} normal: {} vs {}", fc[g][1], fg[g][1]);
        }
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

#[test]
fn taylor_green_decays_like_the_cpu_and_conserves_mass() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (nu, u0, l, n) = (0.1, 0.1, 2.0 * std::f64::consts::PI, 32usize);
        let dx = l / n as f64;
        let mut d = domain([0.0; 3], [l, l, dx], dx);
        for f in BoxFace::ALL {
            *d.faces.get_mut(f) = PERIODIC;
        }
        let dt = 0.3 * dx * dx / (3.0 * nu);
        let half_life = std::f64::consts::LN_2 / (2.0 * nu);
        let m = model(1.0, nu, d, step(half_life, dt, None));
        let (setup, mut flow) = prepare(&m, &[]);
        let grid = *flow.grid();
        let to_lattice = 1.0 / setup.velocity_scale();
        flow.set_equilibrium(|c| {
            let [x, y, _] = grid.centre(c);
            let p = -u0 * u0 / 4.0 * ((2.0 * x).cos() + (2.0 * y).cos());
            (setup.lattice_density(p), [-x.cos() * y.sin() * u0 * to_lattice, x.sin() * y.cos() * u0 * to_lattice, 0.0])
        });
        let Some(mut gpu) = gpu(&flow, streaming) else { continue };
        let mass0 = gpu.mass().unwrap();
        let mut cpu = flow;
        cpu.advance(setup.steps).unwrap();
        gpu.advance(setup.steps).unwrap();
        let drift = (gpu.mass().unwrap() - mass0).abs() / mass0;
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let f = gpu.fields().unwrap();
        let (mut num, mut den) = (0.0, 0.0);
        for c in 0..grid.len() {
            let [x, y, _] = grid.centre(c);
            let mode = [-x.cos() * y.sin(), x.sin() * y.cos()];
            num += f.velocity[c][0] * mode[0] + f.velocity[c][1] * mode[1];
            den += mode[0] * mode[0] + mode[1] * mode[1];
        }
        let amplitude = num / den * setup.velocity_scale();
        let exact = u0 * (-2.0 * nu * setup.steps as f64 * setup.time_step).exp();
        let err = (amplitude / exact - 1.0).abs();
        println!("Taylor–Green on the GPU: amplitude error {err:.2e}; CPU/GPU velocity {du:.2e}, density {drho:.2e}; mass drift {drift:.2e}");
        // The amplitude error is the scheme's at 32 cells (the CPU test shows second order).
        assert!(err < 5e-3, "{err}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(drift < 1e-6, "{drift}");
    }
}

#[test]
fn flow_past_a_sphere_matches_the_cpu_fields_and_forces() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let l = 24usize;
        let side = l as f64;
        let mut d = domain([0.0; 3], [side; 3], 1.0);
        for f in BoxFace::ALL {
            *d.faces.get_mut(f) = PERIODIC;
        }
        d.acceleration = [1e-5, 0.0, 0.0];
        d.obstacles = vec!["Ball".into()];
        let ball = obstacle("Ball", sphere([side / 2.0; 3], 6.0, 3), &["skin"]);
        let m = model(
            1.0,
            1.0 / 6.0,
            d,
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: Some(300),
                time_step: Some(1.0),
                frames: 1,
                magic: None,
                collision: None,
                steady: None,
            },
        );
        let (setup, flow) = prepare(&m, std::slice::from_ref(&ball));
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let skin = cpu.groups().iter().position(|g| g.name == "skin").unwrap();
        let scale = fc[skin].iter().map(|v| v.abs()).fold(0.0, f64::max);
        let dforce = (0..3).map(|a| (fc[skin][a] - fg[skin][a]).abs()).fold(0.0, f64::max) / scale;
        println!("Sphere: CPU/GPU velocity {du:.2e}, density {drho:.2e}; drag CPU {:?} GPU {:?} ({dforce:.2e})", fc[skin], fg[skin]);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(dforce < 5e-5, "{dforce}");
    }
}

#[test]
fn zou_he_channel_matches_the_cpu() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (h, l, umax, nu) = (1.0, 4.0, 0.1, 0.05);
        let dx = h / 10.0;
        let mut d = domain([0.0; 3], [l, h, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let m = model(1.0, nu, d, step(2.0 * l / umax, 0.05 * dx / umax, None));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!("Zou–He channel: CPU/GPU velocity {du:.2e}, density {drho:.2e}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

// ---- Heat ----

/// A fluid with heat: density ρ, kinematic viscosity ν, thermal diffusivity α (conductivity
/// ρ α, unit specific heat).
fn hot_model(density: f64, nu: f64, alpha: f64, domain: FluidDomain, procedure: Procedure) -> Model {
    let mut m = Model::default();
    m.add_material(
        "Fluid",
        Material {
            density: Some(density),
            fluid: Some(Fluid::Newtonian { viscosity: nu * density }),
            thermal: Some(Thermal { conductivity: Some(Directional::Isotropic(density * alpha)), specific_heat: Some(1.0), ..Default::default() }),
            ..Material::default()
        },
    )
    .unwrap();
    m.add_fluid_domain("Flow", domain).unwrap();
    m.add_step(Step::new("Run", procedure), None).unwrap();
    m
}

/// The largest difference of the temperature fields over the CPU's temperature range, and of
/// the groups' heats over the largest CPU heat.
fn heat_difference(cpu: &mut Flow, gpu: &mut GpuFlow) -> (f64, f64) {
    let (a, b) = (cpu.fields().unwrap(), gpu.fields().unwrap());
    assert_eq!(a.temperature.len(), b.temperature.len());
    let (lo, hi) = a.temperature.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| (lo.min(*t), hi.max(*t)));
    let dt = a.temperature.iter().zip(&b.temperature).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max) / (hi - lo).max(1e-30);
    let (ha, hb) = (cpu.heat().unwrap(), gpu.heat().unwrap());
    assert_eq!(ha.len(), hb.len());
    let scale = ha.iter().map(|h| h.abs()).fold(0.0, f64::max).max(1e-30);
    let dh = ha.iter().zip(&hb).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max) / scale;
    (dt, dh)
}

/// Plates held at 20 and 21 with the fluid at rest: the steady temperature is linear between
/// them (anti-bounce-back puts the walls halfway between cells), the heat through both walls
/// is k ΔT A / H, and the GPU agrees with the CPU although it stores the small differences
/// about T₀ = 20 in single precision.
#[test]
fn conduction_between_plates_matches_the_analytic_profile_and_the_cpu() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (h, alpha, n) = (1.0, 0.05, 16usize);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        d.faces.x_min = PERIODIC;
        d.faces.x_max = PERIODIC;
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal { initial_temperature: kelvin(20.0), ..FluidThermal::default() };
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(20.0) };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(21.0) };
        d.thermal = Some(t);
        let m = hot_model(1.0, 0.1, alpha, d, step(8.0 * h * h / alpha, dx * dx / (8.0 * alpha), None));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let f = gpu.fields().unwrap();
        let grid = *cpu.grid();
        let mut err: f64 = 0.0;
        for j in 0..n {
            let expected = 20.0 + (j as f64 + 0.5) / n as f64;
            err = err.max((f.temperature[grid.index([0, j, 0])] - expected).abs());
        }
        let q = alpha * 1.0 / h * dx * dx;
        let heat = gpu.heat().unwrap();
        let cold = cpu.groups().iter().position(|g| g.name == "y_min").unwrap();
        let hot = cpu.groups().iter().position(|g| g.name == "y_max").unwrap();
        // lattice heat (temperature × faces per step) to model units: ρ c h³ / dt
        let scale = setup.thermal.as_ref().unwrap().heat_scale;
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        println!(
            "conduction: profile error {err:.2e}; heat {:.6e} / {:.6e} vs {q:.6e}; GPU vs CPU: temperature {dt:.2e} heat {dh:.2e}",
            heat[cold] * scale,
            heat[hot] * scale
        );
        assert!(err < 5e-5, "{err}");
        assert!((heat[cold] * scale / q - 1.0).abs() < 1e-4 && (heat[hot] * scale / q + 1.0).abs() < 1e-4, "{heat:?}");
        assert!(dt < 2e-5 && dh < 1e-4, "{dt} {dh}");
    }
}

/// A heated channel: a velocity inlet carrying 1 in, a pressure outlet, a flux into the fluid
/// at the bottom wall, the top wall held at 0, periodic across. Every rule of the scalar
/// lattice (pulls, the three wall conditions, the inlet and the outlet) agrees with the CPU.
#[test]
fn a_heated_channel_matches_the_cpu() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (h, l, umax, nu, alpha) = (1.0, 2.0, 0.1, 0.05, 0.02);
        let dx = h / 8.0;
        let mut d = domain([0.0; 3], [l, h, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal::default();
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.y_min = ThermalCondition::Flux { value: 0.05 };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(0.0) };
        d.thermal = Some(t);
        let m = hot_model(1.0, nu, alpha, d, step(2.0 * l / umax, 0.05 * dx / umax, None));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        println!("heated channel: GPU vs CPU: velocity {du:.2e} density {drho:.2e} temperature {dt:.2e} heat {dh:.2e}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(dt < 2e-5 && dh < 1e-4, "{dt} {dh}");
    }
}

/// A differentially heated cavity with buoyancy (the de Vahl Davis geometry at Ra 10³ on 16
/// cells, over the first diffusive time): the temperature drives the flow on both backends
/// alike.
#[test]
fn natural_convection_matches_the_cpu() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (h, pr, alpha, n, u_max) = (1.0, 0.71, 1.0, 16usize, 3.65);
        let nu = pr * alpha;
        let g_beta = 1e3 * nu * alpha / (h * h * h);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [h, h, dx], dx);
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -g_beta, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
        d.thermal = Some(t);
        let dt = (0.05 * dx / (u_max * alpha / h)).min(dx * dx / (4.0 * alpha));
        let m = hot_model(1.0, nu, alpha, d, step(h * h / alpha, dt, None));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        // the fields are in lattice units: the time step puts the cavity's fastest flow at 0.05
        let speed = cpu.fields().unwrap().velocity.iter().map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt()).fold(0.0, f64::max);
        println!("cavity: fastest lattice speed {speed:.3e}; GPU vs CPU: velocity {du:.2e} density {drho:.2e} temperature {dt:.2e} heat {dh:.2e}");
        assert!(speed > 0.01, "the buoyancy drives a flow: {speed}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(dt < 2e-5 && dh < 1e-4, "{dt} {dh}");
    }
}

/// Runs to steady state take their decisions at the same steps on both backends: the plane
/// Poiseuille channel and the de Vahl Davis cavity (Ra 10³, 16 cells) stop as steady at the same
/// sample, with the same monitored values within the stated tolerances.
#[test]
fn runs_to_steady_state_stop_where_the_cpu_does() {
    let steady = |tolerance: f64, window: f64, monitors: Vec<fieldcraft_model::FlowMonitor>| {
        Some(fieldcraft_model::Steady { monitors, tolerance, window: Some(window), min_time: None, divergence: None })
    };
    let fluid = |end_time: f64, dt: f64, steady: Option<fieldcraft_model::Steady>| Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(dt),
        frames: 2,
        magic: None,
        collision: None,
        steady,
    };
    let (h, nu, n) = (1.0, 1.0 / 6.0, 16usize);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [BoxFace::XMin, BoxFace::XMax, BoxFace::ZMin, BoxFace::ZMax] {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [8.0 * nu * 0.1 / (h * h), 0.0, 0.0];
    let walls = vec![fieldcraft_model::FlowMonitor::Force {
        faces: vec![fieldcraft_model::FluidBoundary::Box(BoxFace::YMin), fieldcraft_model::FluidBoundary::Box(BoxFace::YMax)],
    }];
    let channel = model(1.0, nu, d, fluid(24.0, dx * dx, steady(1e-4, 0.5, walls)));

    let (pr, alpha, u_max) = (0.71, 1.0, 3.65);
    let nu = pr * alpha;
    let mut d = domain([0.0; 3], [h, h, dx], dx);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let mut t =
        FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -1e3 * nu * alpha, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    d.thermal = Some(t);
    let dt = (0.05 * dx / u_max).min(dx * dx / 4.0);
    let heat = vec![fieldcraft_model::FlowMonitor::Heat { faces: vec![fieldcraft_model::FluidBoundary::Box(BoxFace::XMin)] }];
    let cavity = hot_model(1.0, nu, alpha, d, fluid(4.0, dt, steady(1e-3, 0.05, heat)));

    for (name, m) in [("channel", channel), ("cavity", cavity)] {
        let limits = fieldcraft_lbm::Limits::default();
        let (setup, flow) = prepare(&m, &[]);
        let Some(mut g) = gpu(&flow, Streaming::default()) else { return };
        let mut cpu = flow;
        let a = fieldcraft_lbm::run(&setup, &mut cpu, &limits).unwrap();
        let b = fieldcraft_lbm::run(&setup, &mut *g, &limits).unwrap();
        let (ra, rb) = (a.steady.unwrap(), b.steady.unwrap());
        println!("{name}: CPU {:?} at step {}, GPU {:?} at step {}", ra.reason, a.stats.steps, rb.reason, b.stats.steps);
        assert_eq!(ra.reason, fieldcraft_lbm::StopReason::Steady, "{ra:?}");
        assert_eq!((rb.reason, b.stats.steps), (ra.reason, a.stats.steps), "{rb:?}");
        let (sa, sb) = (ra.samples.last().unwrap(), rb.samples.last().unwrap());
        for (va, vb) in sa.values.iter().zip(&sb.values) {
            let scale = va.iter().map(|v| v.abs()).fold(0.0, f64::max);
            for (x, y) in va.iter().zip(vb) {
                assert!((x - y).abs() <= 1e-4 * scale, "{name}: {sa:?} vs {sb:?}");
            }
        }
    }
}

/// Readbacks and force caches must describe the same state on both AA orientations,
/// including calls that cross a 64-step submission boundary. The CPU defines the result.
fn compare_incremental(mut cpu: Flow, heated: bool, streaming: Streaming) {
    let Some(mut gpu) = gpu(&cpu, streaming) else { return };
    let mut total = 0;
    for steps in [1, 1, 61, 1, 1, 64, 65] {
        cpu.advance(steps).unwrap();
        gpu.advance(steps).unwrap();
        total += steps;
        assert_eq!(gpu.steps_done(), total);
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!("incremental step {total}: velocity {du:.9e}, density {drho:.9e}");
        assert!(du < 2e-5 && drho < 1e-6, "step {total}: {du} {drho}");
        let mass = Backend::mass(&mut cpu).unwrap();
        assert!((gpu.mass().unwrap() - mass).abs() < 1e-6 * mass.abs());
        let a = Backend::forces(&mut cpu).unwrap();
        let b = gpu.forces().unwrap();
        let scale = a.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max).max(1e-30);
        let error = a.iter().flatten().zip(b.iter().flatten()).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        assert!(error < 5e-5 * scale, "step {total}: force {error} / {scale}");
        gpu.advance(0).unwrap();
        assert_eq!(gpu.forces().unwrap(), b, "a zero-step advance preserves the force cache");
        if heated {
            let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
            assert!(dt < 2e-5 && dh < 1e-4, "step {total}: heat {dt} {dh}");
        }
    }
}

#[test]
fn odd_even_and_chunk_readbacks_match_the_cpu_with_heat_open_and_slip_faces() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let mut d = domain([0.0; 3], [16.0, 8.0, 1.0], 1.0);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.01, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
        d.faces.y_max = FaceCondition::Slip {};
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal { initial_temperature: kelvin(0.5), ..FluidThermal::default() };
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.y_min = ThermalCondition::Flux { value: 0.001 };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(0.0) };
        d.thermal = Some(t);
        let m = hot_model(
            1.0,
            0.1,
            0.05,
            d,
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: Some(194),
                time_step: Some(1.0),
                frames: 1,
                magic: None,
                collision: None,
                steady: None,
            },
        );
        let (_, flow) = prepare(&m, &[]);
        compare_incremental(flow, true, streaming);
    }
}

#[test]
fn odd_even_and_chunk_readbacks_match_the_cpu_around_an_obstacle() {
    // Both layouts: the default depends on the device (see `Streaming`).
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let mut d = domain([0.0; 3], [12.0; 3], 1.0);
        for face in BoxFace::ALL {
            *d.faces.get_mut(face) = PERIODIC;
        }
        d.acceleration = [1e-5, 0.0, 0.0];
        d.obstacles = vec!["Ball".into()];
        let ball = obstacle("Ball", sphere([6.0; 3], 3.0, 2), &["skin"]);
        let m = model(
            1.0,
            1.0 / 6.0,
            d,
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: Some(194),
                time_step: Some(1.0),
                frames: 1,
                magic: None,
                collision: None,
                steady: None,
            },
        );
        let (_, flow) = prepare(&m, &[ball]);
        compare_incremental(flow, false, streaming);
    }
}

/// Fans and porous zones run on the CPU only for now: the backend refuses them by name, before it
/// looks for a device.
#[test]
fn fans_and_porous_zones_are_refused_clearly() {
    let mut d = domain([0.0; 3], [10.0, 2.0, 2.0], 1.0);
    d.porous = vec![fieldcraft_model::PorousZone {
        name: "Vent".into(),
        min: [4.0, 0.0, 0.0],
        max: [6.0, 2.0, 2.0],
        axis: fieldcraft_model::Axis::X,
        resistance: fieldcraft_model::Resistance::LossCoefficient { coefficient: 1.0 },
    }];
    let m = model(1.0, 0.1, d, step(10.0, 1.0, None));
    let (_, flow) = prepare(&m, &[]);
    let _turn = turn();
    let e = GpuFlow::with_options(&flow, Policy::from_env().expect("FIELDCRAFT_GPU names a policy"), Options::default()).err().map(|e| e.to_string());
    assert!(e.as_deref().is_some_and(|e| e.contains("Vent") && e.contains("cpu")), "{e:?}");
}
