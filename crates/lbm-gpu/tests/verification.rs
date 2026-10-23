//! The GPU backend against the CPU scheme and the analytic cases of `fieldcraft-lbm`'s
//! verification suite, with both collision operators. The stated tolerances of the backend, with populations stored as
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
    BoxFace, Collision, Directional, FaceCondition, Fluid, FluidDomain, FluidThermal, Material, Model, Outflow, Procedure, Profile, Step, Thermal,
    ThermalCondition,
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

/// Every collision with every layout.
fn cases() -> impl Iterator<Item = (Collision, Streaming)> {
    COLLISIONS.into_iter().flat_map(|c| [Streaming::InPlace, Streaming::TwoBuffers].map(move |s| (c, s)))
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
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
        let (h, nu, umax, n) = (1.0, 1.0 / 6.0, 0.1, 16usize);
        let dx = h / n as f64;
        let a = 8.0 * nu * umax / (h * h);
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        d.faces.x_min = PERIODIC;
        d.faces.x_max = PERIODIC;
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        d.acceleration = [a, 0.0, 0.0];
        let m = model(1.0, nu, d, collide(4.0 * h * h / nu, dx * dx, Some(3.0 / 16.0), collision));
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
        println!("{collision:?} Poiseuille on the GPU: max error {err:.2e} of the peak; CPU/GPU velocity {du:.2e}, density {drho:.2e}");
        // Exact with Λ = 3/16; the regularised collision's wall is second-order accurate (its
        // CPU test), 1.3e-3 of the peak on 16 cells.
        assert!(err < if collision == Collision::Trt { 2e-5 } else { 2e-3 }, "{err}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

#[test]
fn couette_wall_forces_match_the_cpu() {
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
        let (h, nu, uw, n) = (1.0, 0.1, 0.05, 16usize);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        d.faces.x_min = PERIODIC;
        d.faces.x_max = PERIODIC;
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        d.faces.y_max = FaceCondition::MovingWall { velocity: [uw, 0.0, 0.0] };
        let m = model(1.0, nu, d, collide(3.0 * h * h / nu, 0.5 * dx * dx / (3.0 * nu), None, collision));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let bottom = BoxFace::ALL.iter().position(|f| f.name() == "y_min").unwrap();
        let top = BoxFace::ALL.iter().position(|f| f.name() == "y_max").unwrap();
        let tau = nu * uw / h / setup.force_scale() * dx * dx;
        println!("{collision:?} Couette: wall force CPU {:?} / GPU {:?} (expected ±{tau:.3e} per step)", fc[bottom], fg[bottom]);
        // The tangential force is the shear: a cancellation of the links' exchanges to 1e-6 per link,
        // summed from the stored deviations alone (adding the weights back in single precision left
        // it 1e-4 of itself on an RTX 4090; 3e-6 on lavapipe this way). The normal force is the
        // pressure.
        for g in [bottom, top] {
            println!(
                "{collision:?} Couette group {g}: CPU {:.9e}, GPU {:.9e}, shear parity {:.9e}",
                fc[g][0],
                fg[g][0],
                (fc[g][0] - fg[g][0]).abs() / tau.abs()
            );
            assert!((fc[g][0] - fg[g][0]).abs() <= 5e-5 * tau.abs(), "group {g} shear: {} vs {}", fc[g][0], fg[g][0]);
            assert!((fc[g][1] - fg[g][1]).abs() <= 1e-6 * fc[g][1].abs(), "group {g} normal: {} vs {}", fc[g][1], fg[g][1]);
        }
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

#[test]
fn taylor_green_decays_like_the_cpu_and_conserves_mass() {
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
        let (nu, u0, l, n) = (0.1, 0.1, 2.0 * std::f64::consts::PI, 32usize);
        let dx = l / n as f64;
        let mut d = domain([0.0; 3], [l, l, dx], dx);
        for f in BoxFace::ALL {
            *d.faces.get_mut(f) = PERIODIC;
        }
        let dt = 0.3 * dx * dx / (3.0 * nu);
        let half_life = std::f64::consts::LN_2 / (2.0 * nu);
        let m = model(1.0, nu, d, collide(half_life, dt, None, collision));
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
        println!(
            "{collision:?} Taylor–Green on the GPU: amplitude error {err:.2e}; CPU/GPU velocity {du:.2e}, density {drho:.2e}; mass drift {drift:.2e}"
        );
        // The amplitude error is the scheme's at 32 cells (the CPU test shows second order).
        assert!(err < 5e-3, "{err}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(drift < 1e-6, "{drift}");
    }
}

#[test]
fn flow_past_a_sphere_matches_the_cpu_fields_and_forces() {
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
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
                collision: Some(collision),
                turbulence: None,
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
        println!("{collision:?} Sphere: CPU/GPU velocity {du:.2e}, density {drho:.2e}; drag CPU {:?} GPU {:?} ({dforce:.2e})", fc[skin], fg[skin]);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(dforce < 5e-5, "{dforce}");
    }
}

#[test]
fn zou_he_channel_matches_the_cpu() {
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
        let (h, l, umax, nu) = (1.0, 4.0, 0.1, 0.05);
        let dx = h / 10.0;
        let mut d = domain([0.0; 3], [l, h, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let m = model(1.0, nu, d, collide(2.0 * l / umax, 0.05 * dx / umax, None, collision));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!("{collision:?} Zou–He channel: CPU/GPU velocity {du:.2e}, density {drho:.2e}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

/// A convective outlet (fieldcraft_lbm::outlet) that a vortex crosses with reversed flow: the
/// cells where the flow leaves convect, the reversed ones hold the pressure (Zou–He). Uniform
/// stream U = 0.05 (lattice), periodic sides, a vortex of peak swirl 1.8 U placed at the outlet.
/// Both collisions: under the regularised one the reversed cells hold the equilibrium, as the
/// uniform outlet does, on both backends.
#[test]
fn convective_outlet_matches_the_cpu_with_a_vortex_and_backflow() {
    for (collision, streaming) in
        [Collision::Trt, Collision::Regularized].into_iter().flat_map(|c| [(c, Streaming::InPlace), (c, Streaming::TwoBuffers)])
    {
        let n = 40usize;
        let (u0, r, x0) = (0.05, 6.0, 30.0);
        let side = n as f64;
        let mut d = domain([0.0; 3], [side, side, 1.0], 1.0);
        d.faces.x_min = FaceCondition::Velocity { velocity: [u0, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
        for f in [BoxFace::YMin, BoxFace::YMax, BoxFace::ZMin, BoxFace::ZMax] {
            *d.faces.get_mut(f) = PERIODIC;
        }
        let m = model(1.0, 0.01, d, collide(400.0, 1.0, None, collision));
        let (setup, mut flow) = prepare(&m, &[]);
        let grid = *flow.grid();
        let a = 1.8 * u0 * r / (-0.5f64).exp();
        flow.set_equilibrium(|c| {
            let [x, y, _] = grid.centre(c);
            let (dx, dy) = (x - x0, y - 0.5 * side);
            let e = (-(dx * dx + dy * dy) / (2.0 * r * r)).exp();
            let p = -a * a / (2.0 * r * r) * e * e;
            (setup.lattice_density(p), [u0 - a * dy / (r * r) * e, a * dx / (r * r) * e, 0.0])
        });
        let Some(mut gpu) = gpu(&flow, streaming) else { continue };
        let mut cpu = flow;
        let mut reversed = 0;
        for chunk in [1, 1, 61, 1, 64, 72, 200] {
            cpu.advance(chunk).unwrap();
            gpu.advance(chunk).unwrap();
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            reversed += (0..n).filter(|j| cpu.moments(grid.index([n - 1, *j, 0])).is_some_and(|(_, u)| u[0] < 0.0)).count();
            println!("convective outlet with a vortex ({collision:?}), step {}: CPU/GPU velocity {du:.2e}, density {drho:.2e}", cpu.steps_done());
            assert!(du < 2e-5 && drho < 1e-6, "step {}: {du} {drho}", cpu.steps_done());
        }
        assert!(reversed > 0, "the vortex reverses the flow at the outlet");
    }
}

/// Pressure faces meeting in an edge (x_max and y_max of a 16 × 8 box with an inlet at x_min):
/// each edge cell belongs to one face, so the GPU's one thread per open cell never writes a cell
/// twice. Two entries for one cell raced: a uniform edge diverged on the GPU at step 128 and a
/// convective one differed from the CPU by the whole velocity scale.
#[test]
fn edges_of_two_pressure_faces_match_the_cpu() {
    for outflow in [Outflow::Uniform, Outflow::Convective] {
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let mut d = domain([0.0; 3], [16.0, 8.0, 1.0], 1.0);
            d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
            d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow };
            d.faces.y_max = FaceCondition::Pressure { pressure: 0.0, outflow };
            d.faces.z_min = PERIODIC;
            d.faces.z_max = PERIODIC;
            let m = model(1.0, 0.05, d, step(800.0, 1.0, None));
            let (setup, flow) = prepare(&m, &[]);
            let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            println!("edges of two {outflow:?} faces ({streaming:?}): CPU/GPU velocity {du:.2e}, density {drho:.2e}");
            assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        }
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
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
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
        let m = hot_model(1.0, 0.1, alpha, d, collide(8.0 * h * h / alpha, dx * dx / (8.0 * alpha), None, collision));
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
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
        let (h, l, umax, nu, alpha) = (1.0, 2.0, 0.1, 0.05, 0.02);
        let dx = h / 8.0;
        let mut d = domain([0.0; 3], [l, h, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal::default();
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.y_min = ThermalCondition::Flux { value: 0.05 };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(0.0) };
        d.thermal = Some(t);
        let m = hot_model(1.0, nu, alpha, d, collide(2.0 * l / umax, 0.05 * dx / umax, None, collision));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        println!("{collision:?} heated channel: GPU vs CPU: velocity {du:.2e} density {drho:.2e} temperature {dt:.2e} heat {dh:.2e}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        assert!(dt < 2e-5 && dh < 1e-4, "{dt} {dh}");
    }
}

/// Plates held at 1 that start at a held inlet carrying 0 (the Graetz channel of `fieldcraft-lbm`'s
/// thermal tests): the wall links into the inlet and outlet cells give heat the open boundary
/// resets, which neither backend counts (#675), so the GPU's wall heat is the CPU's.
#[test]
fn wall_heat_at_a_held_inlet_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (h, l, umax, n) = (1.0, 4.0, 0.1, 8usize);
        let dx = h / n as f64;
        let alpha = 2.0 / 3.0 * umax * h / 200.0;
        let mut d = domain([0.0; 3], [l, h, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal { initial_temperature: kelvin(0.0), ..FluidThermal::default() };
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(0.0) };
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        d.thermal = Some(t);
        let m = hot_model(1.0, 0.05, alpha, d, collide(l / umax, 0.05 * dx / umax, None, collision));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        let (hc, hg) = (cpu.heat().to_vec(), gpu.heat().unwrap());
        let walls =
            |heat: &[f64]| ["y_min", "y_max"].map(|w| cpu.groups().iter().position(|g| g.name == w).map_or(0.0, |k| heat[k])).iter().sum::<f64>();
        let (hc, hg) = (walls(&hc), walls(&hg));
        println!("{collision:?} {streaming:?} held walls at a held inlet: wall heat CPU {hc:.6e} GPU {hg:.6e}; temperature {dt:.2e} heat {dh:.2e}");
        assert!(dt < 2e-5 && dh < 1e-4, "{dt} {dh}");
    }
}

/// A second solve on the same device compiles nothing: the kernels stay with the device for the
/// process (#588: DX12 recompiled them for 9 s on every solve).
#[test]
fn a_second_solve_reuses_the_compiled_kernels() {
    let mut d = domain([0.0; 3], [1.0, 1.0, 0.125], 0.125);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let m = model(1.0, 0.1, d, step(1.0, 0.01, None));
    let (_, flow) = prepare(&m, &[]);
    let Some(first) = gpu(&flow, Streaming::InPlace) else { return };
    let compiled = first.pipeline_seconds();
    drop(first);
    let Some(mut second) = gpu(&flow, Streaming::InPlace) else { return };
    println!(
        "kernels: {compiled:.3} s for the first solve (zero if an earlier test compiled them), {:.3} s for the second",
        second.pipeline_seconds()
    );
    assert_eq!(second.pipeline_seconds(), 0.0);
    second.advance(3).unwrap();
    assert!(second.mass().unwrap().is_finite());
}

/// A differentially heated cavity with buoyancy (the de Vahl Davis geometry at Ra 10³ on 16
/// cells, over the first diffusive time): the temperature drives the flow on both backends
/// alike.
#[test]
fn natural_convection_matches_the_cpu() {
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
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
        let m = hot_model(1.0, nu, alpha, d, collide(h * h / alpha, dt, None, collision));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        // the fields are in lattice units: the time step puts the cavity's fastest flow at 0.05
        let speed = cpu.fields().unwrap().velocity.iter().map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt()).fold(0.0, f64::max);
        println!(
            "{collision:?} cavity: fastest lattice speed {speed:.3e}; GPU vs CPU: velocity {du:.2e} density {drho:.2e} temperature {dt:.2e} heat {dh:.2e}"
        );
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
        Some(fieldcraft_model::Steady { monitors, tolerance, window: Some(window), min_time: None, divergence: None, averaging: None })
    };
    let fluid = |end_time: f64, dt: f64, steady: Option<fieldcraft_model::Steady>| Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(dt),
        frames: 2,
        magic: None,
        collision: None,
        turbulence: None,
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
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
        let mut d = domain([0.0; 3], [16.0, 8.0, 1.0], 1.0);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.01, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
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
                collision: Some(collision),
                turbulence: None,
                steady: None,
            },
        );
        let (_, flow) = prepare(&m, &[]);
        compare_incremental(flow, true, streaming);
    }
}

#[test]
fn odd_even_and_chunk_readbacks_match_the_cpu_around_an_obstacle() {
    // Both collisions and both layouts (the default layout depends on the device, see `Streaming`).
    for (collision, streaming) in cases() {
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
                collision: Some(collision),
                turbulence: None,
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

// ---- Smagorinsky's eddy viscosity (#665) ----

/// `p` with Smagorinsky's model, constant `c` and turbulent Prandtl number 0.9.
fn smagorinsky(mut p: Procedure, c: f64) -> Procedure {
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(fieldcraft_model::Turbulence::Smagorinsky { constant: c, prandtl: 0.9 });
    }
    p
}

/// Couette flow with the eddy viscosity a few per cent of the fluid's (`fieldcraft-lbm`'s
/// turbulence tests): the fields, the wall shear and the reply's largest eddy viscosity ratio
/// agree with the CPU.
///
/// This thin fluid (τ = 0.55, 80 000 steps) keeps fewer single-precision digits than the
/// suite's other cases with or without the model. Measured on an RTX 4070 Laptop GPU, GPU
/// against CPU, model off / on: velocity 5.1e-6 / 5.3e-6 (TRT) and 4.2e-5 / 2.3e-5
/// (regularised), wall shear 2.6e-5 / 2.1e-5 and 9.6e-5 / 5.7e-5. So the tolerances here are
/// the model-off ones of this case, with a margin: velocity 1e-4, shear 2e-4.
#[test]
fn couette_flow_with_smagorinsky_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (h, nu, rho, u_wall, n, c) = (1.0, 1e-3, 2.0, 0.1, 16usize, 0.4);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
        let m = model(rho, nu, d, smagorinsky(collide(5.0 * h * h / nu, dx, None, collision), c));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let bottom = BoxFace::ALL.iter().position(|f| f.name() == "y_min").unwrap();
        let shear = (fc[bottom][0] - fg[bottom][0]).abs() / fc[bottom][0].abs();
        let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
        println!(
            "{collision:?} Smagorinsky Couette, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}, wall shear {shear:.2e}; ν_t/ν CPU {rc:.6} GPU {rg:.6}"
        );
        assert!(rc > 0.05, "the eddy viscosity shows: {rc}");
        assert!(du < 1e-4 && drho < 1e-6, "{du} {drho}");
        assert!(shear < 2e-4, "{shear}");
        // The GPU's ratio is the chunk's largest, the CPU's the last step's: at steady state the
        // same, to the f32 cancellation in 1/ω − 1/ω₀ (1.4e-4 of itself measured, regularised).
        assert!((rg - rc).abs() < 5e-4 * rc, "{rg} vs {rc}");
    }
}

/// Plane Poiseuille driven by a body force with the model on (the forced path, and Λ kept by
/// the two relaxation times).
#[test]
fn forced_poiseuille_with_smagorinsky_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (h, nu, umax, n) = (1.0, 1.0 / 6.0, 0.1, 16usize);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d.acceleration = [8.0 * nu * umax / (h * h), 0.0, 0.0];
        let m = model(1.0, nu, d, smagorinsky(collide(4.0 * h * h / nu, dx * dx, Some(3.0 / 16.0), collision), 0.4));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!("{collision:?} Smagorinsky Poiseuille, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}");
        assert!(Backend::eddy_viscosity_ratio(&cpu) > 0.0);
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
    }
}

/// Heated Couette flow with the model on: the turbulent diffusivity ν_t / Pr_t reaches the
/// scalar lattice through `moments`' w, so the temperature and the wall heat agree with the CPU.
/// The thin fluid of the Couette case above, 96 000 steps; measured model off / on: temperature
/// 1.2e-4 / 8.7e-5 of the range, wall heat 1.3e-4 / 1.1e-4 of the largest (TRT); the
/// tolerances are the model-off ones with a margin.
#[test]
fn heated_couette_flow_with_smagorinsky_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (h, nu, alpha, rho, u_wall, n, c) = (1.0, 1e-3, 1e-3, 2.0, 0.1, 16usize, 0.4);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
        let mut t = FluidThermal::default();
        t.faces.y_min = ThermalCondition::Flux { value: 0.3 };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        d.thermal = Some(t);
        let m = hot_model(rho, nu, alpha, d, smagorinsky(collide(6.0 * h * h / alpha, dx, None, collision), c));
        let (setup, cpu) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        println!("{collision:?} Smagorinsky heated Couette, GPU vs CPU: velocity {du:.2e} density {drho:.2e} temperature {dt:.2e} heat {dh:.2e}");
        assert!(du < 1e-4 && drho < 1e-6, "{du} {drho}");
        assert!(dt < 3e-4 && dh < 3e-4, "{dt} {dh}");
    }
}

/// An unsteady case the lattice cannot hold without the closure (`fieldcraft-lbm`'s 48² cavity at
/// Re 5000, cell Reynolds number 104): over its first 1000 steps the GPU follows the CPU.
#[test]
fn a_cavity_held_by_the_closure_follows_the_cpu() {
    for (collision, streaming) in cases() {
        let (l, re, steps) = (48.0, 5000.0, 1000.0);
        let mut d = domain([0.0; 3], [l, l, 1.0], 1.0);
        d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let m = model(1.0, 0.1 * l / re, d, smagorinsky(collide(steps, 1.0, None, collision), 0.17));
        let (setup, flow) = prepare(&m, &[]);
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!(
            "{collision:?} Smagorinsky cavity Re {re}, {steps} steps, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}, ν_t/ν CPU {:.3} GPU {:.3}",
            Backend::eddy_viscosity_ratio(&cpu),
            gpu.eddy_viscosity_ratio()
        );
        assert!(du < 1e-4 && drho < 1e-5, "{du} {drho}");
    }
}

// ---- Conducting solids ----

/// A column across y with a powered conducting slab against the cold plate and another against
/// the hot side's flux, the fluid at rest between them: the solid–solid faces, a held and a flux
/// box face on the solids, the power and the interface links, against the CPU part way to steady
/// state (both solids sub-step: the slab's conductivity is 4 times the fluid's). lavapipe:
/// temperature 2.7e-6 of the range, heat 2.7e-6 of the largest.
#[test]
fn a_composite_wall_with_power_matches_the_cpu() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (n, alpha) = (16usize, 0.05);
        let dx = 1.0 / n as f64;
        let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
        for f in [BoxFace::XMin, BoxFace::XMax, BoxFace::ZMin, BoxFace::ZMax] {
            *d.faces.get_mut(f) = PERIODIC;
        }
        d.obstacles = vec!["Low".into(), "Mid".into()];
        let mut t = FluidThermal { initial_temperature: kelvin(0.5), ..FluidThermal::default() };
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(0.0) };
        t.faces.y_max = ThermalCondition::Flux { value: 0.02 };
        t.solids = vec![conducting("Low", 0.0), conducting("Mid", 0.5 * dx * dx)];
        d.thermal = Some(t);
        let dt = dx * dx / (8.0 * alpha);
        let m = cht_model([1.0, 0.1, alpha], (1.0, [4.0 * alpha; 3], 1.0), d, step(0.5 / alpha, dt, None));
        let obstacles = [
            obstacle("Low", cuboid([-1.0, -1.0, -1.0], [2.0, 0.25, 2.0]), &CUBOID_FACES),
            obstacle("Mid", cuboid([-1.0, 0.5, -1.0], [2.0, 2.0, 2.0]), &CUBOID_FACES),
        ];
        let (setup, cpu) = prepare(&m, &obstacles);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        println!("composite wall with power ({streaming:?}): GPU vs CPU: temperature {dt:.2e} heat {dh:.2e}");
        assert!(dt < 2e-5 && dh < 2e-5, "{dt} {dh}");
    }
}

/// A powered die and its lid across y, joined by a thermal contact, the fluid at rest on both
/// sides and the plates held at 0 and 1: the contact's lower face conductance, which the GPU
/// reads with the others, against the CPU part way to steady state, with the profile's jump
/// across the contact already clear (Apple M4 Max: temperature 4.8e-6, heat 4.9e-6).
#[test]
fn a_thermal_contact_matches_the_cpu() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (n, alpha) = (16usize, 0.05);
        let dx = 1.0 / n as f64;
        let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
        for f in [BoxFace::XMin, BoxFace::XMax, BoxFace::ZMin, BoxFace::ZMax] {
            *d.faces.get_mut(f) = PERIODIC;
        }
        d.obstacles = vec!["Die".into(), "Lid".into()];
        let mut t = FluidThermal { initial_temperature: kelvin(0.5), ..FluidThermal::default() };
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(0.0) };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.solids = vec![
            conducting("Die", 0.5 * dx * dx),
            fieldcraft_model::ConductingSolid { body: "Lid".into(), material: "Other".into(), power: 0.0, limit: None },
        ];
        t.contacts = vec![fieldcraft_model::ThermalContact {
            bodies: ["Die".into(), "Lid".into()],
            resistance: Some(1.25),
            thickness: None,
            conductivity: None,
        }];
        d.thermal = Some(t);
        let dt = dx * dx / (8.0 * alpha);
        let mut m = Model::default();
        m.add_material("Fluid", hot_fluid(1.0, 0.1, alpha)).unwrap();
        add_solid_material(&mut m, "Solid", 1.0, [4.0 * alpha; 3], 1.0);
        add_solid_material(&mut m, "Other", 1.0, [10.0 * alpha; 3], 1.0);
        m.add_fluid_domain("Flow", d).unwrap();
        m.add_step(Step::new("Run", step(0.5 / alpha, dt, None)), None).unwrap();
        let obstacles = [
            obstacle("Die", cuboid([-1.0, 0.25, -1.0], [2.0, 0.5, 2.0]), &CUBOID_FACES),
            obstacle("Lid", cuboid([-1.0, 0.5, -1.0], [2.0, 0.75, 2.0]), &CUBOID_FACES),
        ];
        let (setup, cpu) = prepare(&m, &obstacles);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        // the jump across the contact (cells 7 and 8 sit on either side of y = 0.5)
        let jump = cpu.temperature(8).unwrap() - cpu.temperature(7).unwrap();
        println!("thermal contact ({streaming:?}): GPU vs CPU: temperature {dt:.2e} heat {dh:.2e}; jump across the contact {jump:.4}");
        assert!(jump.abs() > 1e-2, "{jump}");
        assert!(dt < 2e-5 && dh < 2e-5, "{dt} {dh}");
    }
}

/// A powered block conducting in a channel flow (parabolic inlet, outlet): the interface links
/// under advection, against the CPU in velocity, temperature and every face's heat. lavapipe:
/// temperature 2.5–2.7e-6, heat 3.2–3.3e-6, velocity 1.6–2.0e-6.
#[test]
fn a_heated_block_in_a_channel_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (n, nu, alpha) = (8usize, 0.02, 0.02);
        let dx = 1.0 / n as f64;
        let mut d = domain([0.0; 3], [3.0, 1.0, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        d.obstacles = vec!["Block".into()];
        let mut t = FluidThermal { initial_temperature: kelvin(1.0), ..FluidThermal::default() };
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.solids = vec![conducting("Block", 0.01)];
        d.thermal = Some(t);
        let dt = 0.25 * dx * dx / nu;
        let m = cht_model([1.0, nu, alpha], (2.0, [0.5; 3], 1.0), d, collide(40.0, dt, None, collision));
        let block = [obstacle("Block", cuboid([0.75, 0.375, -1.0], [1.25, 0.625, 1.0]), &CUBOID_FACES)];
        let (setup, cpu) = prepare(&m, &block);
        let Some((mut cpu, mut gpu)) = both(&setup, cpu, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
        println!("{collision:?} heated block in a channel: GPU vs CPU: velocity {du:.2e} density {drho:.2e} temperature {dt:.2e} heat {dh:.2e}");
        // (the regularised collision round the block: density 1.4e-6 on lavapipe)
        assert!(du < 2e-5 && drho < 2e-6, "{du} {drho}");
        assert!(dt < 2e-5 && dh < 2e-5, "{dt} {dh}");
    }
}

/// A powered plate along the whole channel floor, through the inlet and the outlet: the open
/// faces' cells above it have links into it, so the inlet must hold its temperature on them and
/// the outlet copy its neighbour's populations after the interface, as on the CPU.
#[test]
fn a_conductor_through_the_inlet_and_the_outlet_matches_the_cpu() {
    let (n, nu, alpha) = (8usize, 0.02, 0.02);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [2.0, 1.0, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Plate".into()];
    let mut t = FluidThermal { initial_temperature: kelvin(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.solids = vec![conducting("Plate", 0.02)];
    d.thermal = Some(t);
    let dt = 0.25 * dx * dx / nu;
    let m = cht_model([1.0, nu, alpha], (2.0, [0.5; 3], 1.0), d, step(20.0, dt, None));
    let plate = [obstacle("Plate", cuboid([-1.0, -1.0, -1.0], [3.0, 0.25, 1.0]), &CUBOID_FACES)];
    let (setup, cpu) = prepare(&m, &plate);
    let Some((mut cpu, mut gpu)) = both(&setup, cpu, Streaming::InPlace) else { return };
    let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
    println!("plate through the open faces: GPU vs CPU: temperature {dt:.2e} heat {dh:.2e}");
    assert!(dt < 2e-5 && dh < 2e-5, "{dt} {dh}");
}

/// A powered conductor filling the whole box, a held and a flux face on it: no fluid and so no
/// fluid links, and the box faces' heat still comes from the solid as on the CPU.
#[test]
fn a_conductor_filling_the_box_reports_its_box_face_heat() {
    let n = 8usize;
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [BoxFace::XMin, BoxFace::XMax, BoxFace::ZMin, BoxFace::ZMax] {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.obstacles = vec!["Fill".into()];
    let mut t = FluidThermal { initial_temperature: kelvin(0.5), ..FluidThermal::default() };
    t.faces.y_min = ThermalCondition::Temperature { value: kelvin(0.0) };
    t.faces.y_max = ThermalCondition::Flux { value: 0.02 };
    t.solids = vec![conducting("Fill", 0.5 * dx * dx)];
    d.thermal = Some(t);
    let dt = dx * dx / 0.4;
    let m = cht_model([1.0, 0.1, 0.05], (1.0, [0.2; 3], 1.0), d, step(0.5, dt, None));
    let fill = [obstacle("Fill", cuboid([-1.0; 3], [2.0; 3]), &CUBOID_FACES)];
    let (setup, cpu) = prepare(&m, &fill);
    let Some((mut cpu, mut gpu)) = both(&setup, cpu, Streaming::InPlace) else { return };
    let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
    println!("conductor filling the box: GPU vs CPU: temperature {dt:.2e} heat {dh:.2e}");
    assert!(dt < 2e-5 && dh < 2e-5, "{dt} {dh}");
}

/// A heated cube of small heat capacity with strong links (τ_T 0.52, conducting 8 times the
/// fluid, 0.03 of its heat capacity): the implicit links on the GPU stay stable and agree with
/// the CPU (lavapipe: temperature 2.9e-6, heat 5.3e-6).
#[test]
fn a_low_capacity_solid_matches_the_cpu() {
    let (n, alpha) = (12usize, 1.0);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [1.0; 3], dx);
    d.obstacles = vec!["Cube".into()];
    let mut t = FluidThermal::default();
    for f in BoxFace::ALL {
        *t.faces.get_mut(f) = ThermalCondition::Temperature { value: kelvin(0.0) };
    }
    t.solids = vec![conducting("Cube", 1.0)];
    d.thermal = Some(t);
    let dt = 0.005 * dx * dx / alpha;
    let m = cht_model([1.0, alpha, alpha], (0.03, [8.0 * alpha; 3], 1.0), d, step(300.0 * dt, dt, None));
    let cube = [obstacle("Cube", cuboid([0.25; 3], [0.75; 3]), &CUBOID_FACES)];
    let (setup, cpu) = prepare(&m, &cube);
    let Some((mut cpu, mut gpu)) = both(&setup, cpu, Streaming::InPlace) else { return };
    let (dt, dh) = heat_difference(&mut cpu, &mut gpu);
    println!(
        "low-capacity cube: GPU vs CPU: temperature {dt:.2e} heat {dh:.2e}; conducting {:?}",
        gpu.conducting().iter().map(Vec::len).collect::<Vec<_>>()
    );
    assert_eq!(gpu.conducting(), cpu.conducting());
    assert!(dt < 2e-5 && dh < 5e-5, "{dt} {dh}");
}
