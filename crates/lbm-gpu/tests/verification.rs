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
    BoxFace, Collision, Directional, FaceCondition, Fluid, FluidDomain, FluidThermal, HrrSigma, Material, Model, Outflow, Procedure, Profile, Step,
    Thermal, ThermalCondition, Turbulence,
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

/// Every collision (the multiple-relaxation-time one too) with every layout.
fn cases() -> impl Iterator<Item = (Collision, Streaming)> {
    COLLISIONS.into_iter().chain([Collision::Mrt]).flat_map(|c| [Streaming::InPlace, Streaming::TwoBuffers].map(move |s| (c, s)))
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
        // Exact with Λ = 3/16; the regularised and multiple-relaxation-time collisions' walls are
        // second-order accurate (their CPU tests), 1.3e-3 and 3.0e-3 of the peak on 16 cells.
        let bound = match collision {
            Collision::Trt => 2e-5,
            Collision::Regularized => 2e-3,
            Collision::Mrt => 4e-3,
            // Not among `cases()`: the GPU refuses the hybrid collision.
            Collision::Hrr => continue,
        };
        assert!(err < bound, "{err}");
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
        // The amplitude error is the scheme's at 32 cells (the CPU tests show second order; the
        // multiple-relaxation-time collision's is 5.5e-3 there).
        assert!(err < if collision == Collision::Mrt { 6e-3 } else { 5e-3 }, "{err}");
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
                sigma: None,
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
        sigma: None,
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
                sigma: None,
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
                sigma: None,
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
        *turbulence = Some(fieldcraft_model::Turbulence::Smagorinsky { constant: c, prandtl: 0.9, van_driest: None, wall_model: None });
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
        // Both backends report the last step's ratio: the same, to the f32 cancellation in 1/ω − 1/ω₀ (1.4e-4 of itself measured, regularised).
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
/// The regularized collision's single-precision gap peaks in the lid's start-up and then decays
/// (measured on an NVIDIA GB10, of the largest speed: 1.2e-6 at step 100, 2.3e-4 at 500,
/// 1.9e-4 at 1000, 3.5e-5 at 1500; llvmpipe 7.8e-5 at 1000), while the eddy viscosity agrees to
/// four digits throughout, so its bound is the measured peak with a margin; TRT stays near 3e-6.
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
        let (du_max, drho_max) = match collision {
            Collision::Regularized => (5e-4, 1e-5),
            _ => (1e-4, 1e-5),
        };
        assert!(du < du_max && drho < drho_max, "{du} {drho}");
        // The reply's ratio is the last step's on both backends, not the chunk's largest.
        let (ratio_cpu, ratio_gpu) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
        assert!((ratio_gpu - ratio_cpu).abs() < 1e-3 * ratio_cpu, "{ratio_cpu} {ratio_gpu}");
    }
}

/// A body force changed while the flow runs (`Backend::set_acceleration`, for a flow held at a
/// constant flow rate): plane Poiseuille flow driven by a for half its steps and 2a for the rest
/// on both backends, which agree, as do their mean velocities (`Backend::mean_velocity`).
#[test]
fn a_body_force_changed_mid_run_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (h, nu, umax, n) = (1.0, 1.0 / 6.0, 0.05, 16usize);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d.acceleration = [8.0 * nu * umax / (h * h), 0.0, 0.0];
        let m = model(1.0, nu, d, collide(2.0 * h * h / nu, dx * dx, None, collision));
        let (setup, mut cpu) = prepare(&m, &[]);
        let Some(mut gpu) = gpu(&cpu, streaming) else { continue };
        let half = setup.steps / 2;
        let doubled = setup.acceleration.map(|v| 2.0 * v);
        cpu.advance(half).unwrap();
        gpu.advance(half).unwrap();
        Backend::set_acceleration(&mut cpu, doubled).unwrap();
        gpu.set_acceleration(doubled).unwrap();
        cpu.advance(setup.steps - half).unwrap();
        gpu.advance(setup.steps - half).unwrap();
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let peak = gpu.fields().unwrap().velocity.iter().map(|u| u[0]).fold(0.0, f64::max) * setup.velocity_scale();
        println!("{collision:?} force doubled mid-run: peak {peak:.4} (from {umax}), GPU vs CPU velocity {du:.2e}, density {drho:.2e}");
        assert!(peak > 1.5 * umax, "the doubled force shows: {peak}");
        assert!(du < 2e-5 && drho < 1e-6, "{du} {drho}");
        // The mean velocity a constant flow rate controls, on both backends.
        let (mc, mg) = (Backend::mean_velocity(&mut cpu, 0).unwrap(), gpu.mean_velocity(0).unwrap());
        assert!((mg / mc - 1.0).abs() < 3e-5 && mc > 0.0, "{mg} vs {mc}");
        assert!(gpu.mean_velocity(1).unwrap().abs() < 1e-6 * mc && gpu.mean_velocity(3).is_err());
        assert!(gpu.set_acceleration([f64::NAN, 0.0, 0.0]).is_err());
    }
}

// ---- The wall model ----

/// `p` with Smagorinsky's model, constant `c`, and the Spalding wall model.
fn wall_modelled(mut p: Procedure, c: f64) -> Procedure {
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(fieldcraft_model::Turbulence::Smagorinsky {
            constant: c,
            prandtl: 0.9,
            van_driest: None,
            wall_model: Some(fieldcraft_model::WallModel::Musker),
        });
    }
    p
}

/// Couette flow on 8 cells at ν = 5e-6 under the wall model (the first cells above the
/// viscous sublayer), both walls sampled at their second cells and rebuilt, the moving one
/// relative to its velocity, 16 000 steps into the start-up: the fields and the largest eddy
/// viscosity ratio (the wall cells' mixing length) agree with the CPU. Measured on the GB10:
/// velocity 1.8e-6 (TRT) and 3.1e-5 (regularised), ν_t/ν to 1e-5. The wall shear is not compared
/// here: at ν = 5e-6 it is about 1e-6 of the wall's pressure force (lattice units), below what
/// the f32 momentum exchange resolves; the channel below compares the wall drag.
#[test]
fn couette_flow_with_the_wall_model_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (h, nu, rho, u_wall, n, c) = (1.0, 5e-6, 1.2, 0.1, 8usize, 0.4);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
        let m = model(rho, nu, d, wall_modelled(collide(16000.0 * dx, dx, None, collision), c));
        let (setup, flow) = prepare(&m, &[]);
        assert!(flow.wall_model());
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
        println!("{collision:?} wall-modelled Couette, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}; ν_t/ν CPU {rc:.4} GPU {rg:.4}");
        assert!(rc > 1.0, "the wall model's viscosity shows: {rc}");
        // The multiple-relaxation-time collision's density keeps fewer f32 digits here (1.2e-5
        // measured on the GB10; its velocity 6.3e-6).
        let drho_max = if collision == Collision::Mrt { 3e-5 } else { 1e-6 };
        assert!(du < 1e-4 && drho < drho_max, "{du} {drho}");
        assert!((rg - rc).abs() < 1e-3 * rc, "{rg} vs {rc}");
    }
}

/// A coarse turbulent-like channel (walls at y = 0 and 2, periodic in x and z, 6 cells per
/// half-height, first cell at y⁺ ≈ 15 at Re_τ 180) started from a perturbed blunt profile and run
/// 150 steps under the wall model: a fully three-dimensional flow with every wall cell's u_τ
/// moving, against the CPU before f32 round-off grows chaotically. The wall drag is the
/// bounce-back's exchange plus what the reconstruction takes, two parts about 450 times the body
/// force in this start-up that nearly cancel (the net drag is 0.1 to 9 times the body force), so
/// the states' f32 differences show in it magnified. Measured on the GB10 after 150 steps: velocity
/// 2.3e-6 (TRT), 2.7e-6 (regularised) and 1.0e-6 (MRT), the wall cells at their imposed velocity
/// within as much; wall drag 2.8e-4, 3.2e-4 and 9.8e-6 of itself, about 5e-6 of either part. The
/// steady channel below compares the drag to the force tests' tolerance.
#[test]
fn a_coarse_wall_modelled_channel_follows_the_cpu() {
    for (collision, streaming) in cases() {
        let (setup, flow) = coarse_channel(|p| wall_modelled(p, 0.1), collision);
        let grid = *flow.grid();
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        // The GPU's fields report each modelled wall cell at the velocity its reconstruction
        // imposes, as the CPU's do (not the streamed state the next reconstruction overwrites).
        let fields = gpu.fields().unwrap();
        let walls = cpu.walls().unwrap();
        let scale = fields.velocity.iter().map(|v| v[0].abs()).fold(0.0, f64::max);
        let modelled: Vec<[f64; 3]> =
            (0..grid.len()).filter_map(|c| walls.modelled(c).map(|(u, _)| [0, 1, 2].map(|k| fields.velocity[c][k] - u[k]))).collect();
        let imposed = modelled.iter().flatten().map(|d| d.abs()).fold(0.0, f64::max) / scale;
        assert!(!modelled.is_empty() && imposed < 1e-4, "GPU wall cells vs the imposed velocity: {imposed}");
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let wall = |f: &[[f64; 3]]| f[BoxFace::YMin.index()][0] + f[BoxFace::YMax.index()][0];
        let drag = (wall(&fc) - wall(&fg)).abs() / wall(&fc).abs();
        println!(
            "{collision:?} wall-modelled channel, {} cells, 150 steps, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}, wall drag {drag:.2e}; wall cells vs the imposed velocity {imposed:.2e}",
            grid.len()
        );
        assert!(du < 1e-4 && drho < 1e-5, "{du} {drho}");
        assert!(drag < 1e-3, "{drag}");
    }
}

/// The coarse channel above (walls at y = 0 and 2, periodic in x and z, 6 cells per
/// half-height, Re_τ 180), started from a perturbed blunt profile, 150 steps, with the
/// turbulence `turbulence` sets.
fn coarse_channel(turbulence: impl Fn(Procedure) -> Procedure, collision: Collision) -> (Setup, Flow) {
    channel_at(180.0, turbulence, collision)
}

/// [`coarse_channel`] at friction Reynolds number `re_tau`, the bulk speed scaled with it (30 runs
/// without a turbulence model: τ 0.52).
fn channel_at(re_tau: f64, turbulence: impl Fn(Procedure) -> Procedure, collision: Collision) -> (Setup, Flow) {
    use std::f64::consts::PI;
    let (n, ub) = (6usize, 15.7 * re_tau / 180.0);
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let m = model(1.0, 1.0 / re_tau, d, turbulence(collide(150.0 * dt, dt, None, collision)));
    let (setup, mut flow) = prepare(&m, &[]);
    let grid = *flow.grid();
    let to_lattice = 1.0 / setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0;
        let wall = 1.0 - eta * eta;
        let u = 1.15 * ub * (1.0 - eta.abs().powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u * to_lattice, v * to_lattice, w * to_lattice])
    });
    (setup, flow)
}

// ---- Van Driest's damping ----

/// `p` with Smagorinsky's model, constant `c`, Van Driest's damping at A⁺ 26 and, with
/// `wall_model`, the wall model too.
fn van_driest(mut p: Procedure, c: f64, wall_model: bool) -> Procedure {
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(fieldcraft_model::Turbulence::Smagorinsky {
            constant: c,
            prandtl: 0.9,
            van_driest: Some(26.0),
            wall_model: wall_model.then_some(fieldcraft_model::WallModel::Musker),
        });
    }
    p
}

/// The coarse channel above with Van Driest's damping, on the bounce-back walls alone and under
/// the wall model: every wall cell's u_τ from Spalding's law and every cell's damped constant
/// change each step in a three-dimensional flow. Against the CPU after 150 steps, at the wall
/// model's tolerances: the fields, the largest eddy viscosity ratio and the wall cells' largest
/// y⁺ (`Backend::wall_y_plus`). Measured on the GB10, TRT / regularised / MRT: on the
/// bounce-back walls, velocity 3.3e-6 / 5.6e-6 / 5.5e-7, density 5.3e-7 / 2.6e-7 / 7.1e-8, ν_t/ν
/// 1.3e-5 / 1.4e-5 / 1.2e-4 of itself, y⁺ (11.8 to 12.6) 3.9e-6 / 6.3e-6 / 1.3e-6; under the wall
/// model, velocity 2.1e-6 / 2.8e-6 / 9.8e-7, density 4.0e-7 / 2.0e-7 / 1.0e-7, ν_t/ν 1.1e-5 /
/// 1.3e-5 / 2.3e-5, y⁺ (13.6 to 14.6) within 2e-7. With the damping left out of the GPU's
/// collision the velocity differs by 0.12 of the largest speed.
#[test]
fn a_coarse_channel_with_van_driest_damping_follows_the_cpu() {
    for wall_model in [false, true] {
        for (collision, streaming) in cases() {
            let (setup, flow) = coarse_channel(|p| van_driest(p, 0.1, wall_model), collision);
            assert!(flow.van_driest() && flow.wall_model() == wall_model);
            let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
            let (yc, yg) = (Backend::wall_y_plus(&cpu).unwrap(), gpu.wall_y_plus().unwrap());
            println!(
                "{collision:?} {streaming:?} channel with Van Driest damping{}, 150 steps, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}; ν_t/ν CPU {rc:.5} GPU {rg:.5} ({:.1e}); wall y⁺ CPU {yc:.5} GPU {yg:.5} ({:.1e})",
                if wall_model { " and the wall model" } else { "" },
                (rg - rc).abs() / rc,
                (yg - yc).abs() / yc
            );
            assert!(rc > 0.0 && yc > 1.0, "the damping acts: ν_t/ν {rc}, y⁺ {yc}");
            assert!(du < 1e-4 && drho < 1e-5, "{du} {drho}");
            assert!((rg - rc).abs() < 1e-3 * rc, "{rg} vs {rc}");
            assert!((yg - yc).abs() < 1e-3 * yc, "{yg} vs {yc}");
        }
    }
}

// ---- Han, Ooka and Kikumoto's on-node wall-function bounce (`spalding_bounce`, their eq. 14) ----

/// The subgrid models the on-node bounce is checked with: Smagorinsky's (C_s 0.1) with Van
/// Driest's damping (A⁺ 26), as the paper runs it, and Vreman's (C_s 0.18).
fn on_node_turbulence(vreman: bool) -> Turbulence {
    let wall_model = Some(fieldcraft_model::WallModel::SpaldingBounce);
    if vreman {
        Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model }
    } else {
        Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: Some(26.0), wall_model }
    }
}

/// The collisions the on-node bounce is checked with, and their σ: regularised, hybrid at σ 0.99,
/// MRT and TRT.
const ON_NODE_COLLISIONS: [(Collision, Option<HrrSigma>); 4] =
    [(Collision::Regularized, None), (Collision::Hrr, Some(HrrSigma::Value(0.99))), (Collision::Mrt, None), (Collision::Trt, None)];

/// The CPU's coarse on-node channel (`crates/lbm/tests/wall_bounce.rs`): 6 cells per half-height,
/// Re_τ 180, periodic in x and z, the box 2δ + Δx high (13 cells) so that the walls through the
/// outermost cells' centres are 2δ apart, started from a perturbed blunt profile, driven by a body
/// force, `steps` steps.
fn on_node_channel(turbulence: Turbulence, collision: Collision, sigma: Option<HrrSigma>, steps: f64) -> (Setup, Flow) {
    use std::f64::consts::PI;
    let (n, re_tau, ub) = (6usize, 180.0, 15.7);
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0 + dx, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let p = les(sigma, Some(turbulence))(collide(steps * dt, dt, None, collision));
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / re_tau, d, p), &[]);
    let grid = *flow.grid();
    assert_eq!(grid.cells[1], 2 * n + 1);
    let to_lattice = 1.0 / setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0 - 0.5 * dx;
        let wall = (1.0 - eta * eta).max(0.0);
        let u = 1.15 * ub * (1.0 - eta.abs().min(1.0).powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u * to_lattice, v * to_lattice, w * to_lattice])
    });
    (setup, flow)
}

/// What one GPU-against-CPU run of the on-node bounce measured: the velocity (over the largest
/// CPU speed) and density differences; the walls' drag (x, both walls) on the GPU against the
/// CPU's; the GPU walls' carried stress against the stress the GPU sampled (x and z, over the x
/// stress); and the GPU's sampled stress against the CPU's.
struct OnNode {
    du: f64,
    drho: f64,
    drag: f64,
    carried: f64,
    sampled: f64,
}

/// The stated tolerances of the on-node bounce on the GPU, those of the coarse wall-modelled
/// channel: velocity within 1e-4 of the largest speed, density within 1e-5, the walls' drag within
/// 1e-3 of the CPU's, the carried stress within 1e-3 of the sampled one (its own and the CPU's).
fn check_on_node(what: &str, r: &OnNode) {
    println!(
        "{what}, GPU vs CPU: velocity {:.2e}, density {:.2e}, wall drag {:.2e}; GPU carried vs its sampled stress {:.2e}, GPU vs CPU sampled stress {:.2e}",
        r.du, r.drho, r.drag, r.carried, r.sampled
    );
    assert!(r.du < 1e-4 && r.drho < 1e-5, "{what}: {} {}", r.du, r.drho);
    assert!(r.drag < 1e-3 && r.carried < 1e-3 && r.sampled < 1e-3, "{what}: {} {} {}", r.drag, r.carried, r.sampled);
}

/// `steps` steps of `flow` on both from step 0, and what [`OnNode`] measures.
fn on_node_difference(flow: Flow, steps: u64, streaming: Streaming) -> Option<OnNode> {
    let mut cpu = flow;
    assert!(cpu.on_node_bounce());
    let mut gpu = gpu(&cpu, streaming)?;
    cpu.advance(steps).unwrap();
    gpu.advance(steps).unwrap();
    let (du, drho) = field_difference(&mut cpu, &mut gpu);
    let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
    let walls = |f: &[[f64; 3]], k: usize| f[BoxFace::YMin.index()][k] + f[BoxFace::YMax.index()][k];
    let drag = (walls(&fc, 0) - walls(&fg, 0)).abs() / walls(&fc, 0).abs();
    // The GPU's sampled stresses (its wall rows) and the CPU's (its boundary nodes).
    let rows = gpu.wall_targets().unwrap();
    let sg = |k: usize| rows.iter().map(|r| r.0[k]).sum::<f64>();
    let sc = |k: usize| cpu.wall_nodes().unwrap().nodes().iter().map(|w| w.stress[k]).sum::<f64>();
    let carried = [0, 2].map(|k| (walls(&fg, k) - sg(k)).abs() / sg(0).abs()).into_iter().fold(0.0, f64::max);
    let sampled = [0, 2].map(|k| (sg(k) - sc(k)).abs() / sc(0).abs()).into_iter().fold(0.0, f64::max);
    Some(OnNode { du, drho, drag, carried, sampled })
}

/// Han, Ooka and Kikumoto's on-node bounce (their eq. 14) in the perturbed coarse channel, 150
/// steps from the start (no warm-up), on the regularised, hybrid (σ 0.99), MRT and TRT
/// collisions, with Smagorinsky's model and Van Driest's damping and with Vreman's, both layouts:
/// the fields, the walls' drag and stresses at [`check_on_node`]'s tolerances.
#[test]
fn the_on_node_wall_function_bounce_follows_the_cpu() {
    for vreman in [false, true] {
        for (collision, sigma) in ON_NODE_COLLISIONS {
            for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
                let (setup, flow) = on_node_channel(on_node_turbulence(vreman), collision, sigma, 150.0);
                let Some(r) = on_node_difference(flow, setup.steps, streaming) else { return };
                let model = if vreman { "Vreman" } else { "Smagorinsky + Van Driest" };
                check_on_node(&format!("on-node bounce, {} + {model}, {streaming:?}, coarse channel, 150 steps", collision.name()), &r);
            }
        }
    }
}

/// The same in Han et al.'s A20 channel (Re_τ 640, 20 cells per δ, 126 × 41 × 63 cells, their
/// Mach number), 300 steps from its perturbed start, on every collision and both subgrid models,
/// both layouts.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_on_node_wall_function_bounce_follows_the_cpu_in_han_s_a20_channel() {
    let case = &HAN_CASES[0];
    for vreman in [false, true] {
        for (collision, sigma) in ON_NODE_COLLISIONS {
            for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
                let HanFlow { flow, .. } = han_flow(case, collision, sigma, on_node_turbulence(vreman), true, 1.0).unwrap();
                let Some(r) = on_node_difference(flow, 300, streaming) else { return };
                let model = if vreman { "Vreman" } else { "Smagorinsky + Van Driest" };
                check_on_node(&format!("on-node bounce, {} + {model}, {streaming:?}, Han A20, 300 steps", collision.name()), &r);
            }
        }
    }
}
/// The on-node bounce twice on one device gives the same fields and wall forces bit for bit (the
/// paper's MRT with Smagorinsky's model and Van Driest's damping, and the hybrid collision with
/// Vreman's, in the coarse channel), both layouts.
#[test]
fn the_on_node_wall_function_bounce_is_deterministic_on_its_device() {
    for (vreman, collision, sigma) in [(false, Collision::Mrt, None), (true, Collision::Hrr, Some(HrrSigma::Value(0.99)))] {
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let (setup, flow) = on_node_channel(on_node_turbulence(vreman), collision, sigma, 150.0);
            let run = || {
                let mut g = gpu(&flow, streaming)?;
                g.advance(setup.steps).unwrap();
                let mut bits = field_bits(&g.fields().unwrap());
                bits.extend(g.forces().unwrap().iter().flatten().map(|v| v.to_bits()));
                Some(bits)
            };
            let (Some(a), Some(b)) = (run(), run()) else { return };
            println!("on-node bounce, {}, {streaming:?}: two runs bit for bit: {}", collision.name(), a == b);
            assert!(a == b);
        }
    }
}

/// Spalding's y⁺(u⁺) with Han et al.'s constants (κ 0.41, B 5.5).
fn spalding_y_plus(up: f64) -> f64 {
    let k = 0.41 * up;
    up + (-0.41f64 * 5.5).exp() * (k.exp() - 1.0 - k - k * k / 2.0 - k * k * k / 6.0)
}

/// The k_x = π mode's energy `Σ_rows Σ_components |û|²`, `û = Σ_{x,z} (−1)^x u / (n_x n_z)`, over
/// the rows `rows` of a velocity field on `grid`.
fn odd_even_energy(grid: &fieldcraft_lbm::Grid, velocity: &[[f64; 3]], rows: std::ops::RangeInclusive<usize>) -> f64 {
    let [nx, _, nz] = grid.cells;
    let mut e = 0.0;
    for j in rows {
        let mut a = [0.0; 3];
        for k in 0..nz {
            for i in 0..nx {
                let u = velocity[grid.index([i, j, k])];
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                for c in 0..3 {
                    a[c] += sign * u[c];
                }
            }
        }
        e += a.iter().map(|x| (x / (nx * nz) as f64).powi(2)).sum::<f64>();
    }
    e
}

/// The least-squares slope per step of `y` over its second half.
fn tail_slope(y: &[f64]) -> f64 {
    let tail = &y[y.len() / 2..];
    let n = tail.len() as f64;
    let (mean_t, mean_y) = ((n - 1.0) / 2.0, tail.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (t, v) in tail.iter().enumerate() {
        sxy += (t as f64 - mean_t) * (v - mean_y);
        sxx += (t as f64 - mean_t).powi(2);
    }
    sxy / sxx
}

/// The CPU's stationary k_x = π toy (`crates/lbm/tests/wall_stability.rs`, whose notes derive it)
/// on the GPU for its trunk case: Han et al.'s on-node eq. 14 on the regularised collision at
/// τ̄ 0.5015 with the first node at u⁺ 5, no subgrid model, in a periodic box 16 × 9 × 2 (N = 4),
/// the laminar base held steady by the body force that balances both walls' Spalding stress. A
/// checkerboard in x of u_x is seeded in the first fluid row at each wall, and the growth rate per
/// step of its amplitude (the least-squares slope of ½ ln E over the second half of 4000 steps)
/// must not be positive, as on the CPU, and must be the CPU's (run beside it) within 10 %. The
/// seed is 1e-3 of a lattice speed (the base is about 0.03): with the CPU's 1e-6, or 1e-5, single
/// precision's rounding of the base (ε |u| ≈ 2e-9 per step, against a decay of 1.2e-4 per step)
/// holds the mode at about 4e-6 on the GB10 and on llvmpipe, so the rate measured is the
/// rounding's. At 1e-3 the CPU's mode is linear (its energy is 10⁶ times that at 1e-6 to 4 digits).
/// Measured on the GB10: −1.195e-4 per step on both layouts, the CPU −1.169e-4.
#[test]
fn the_on_node_bounce_damps_the_odd_even_mode_of_its_wall_rows_on_the_gpu() {
    const N: usize = 4;
    let (nx, nz) = (16usize, 2usize);
    let (tau, steps, seed) = (0.5015, 4000u64, 1e-3);
    let nu = (tau - 0.5) / 3.0;
    let (ny, rows, d) = (2 * N + 1, 2 * N - 1, 1.0);
    let h = N as f64;
    let up = 5.0;
    let ut = spalding_y_plus(up) * nu / d;
    let g = 2.0 * ut * ut / rows as f64;
    let u1 = ut * up;
    let profile = |y: f64| u1 + g / (2.0 * nu) * ((h - d).powi(2) - (h - y).powi(2));
    assert!(profile(h) <= 0.1, "the centreline under lattice speed 0.1");
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let mut dom = domain([0.0; 3], [nx as f64, ny as f64, nz as f64], 1.0);
        for f in [&mut dom.faces.x_min, &mut dom.faces.x_max, &mut dom.faces.z_min, &mut dom.faces.z_max] {
            *f = PERIODIC;
        }
        dom.acceleration = [g, 0.0, 0.0];
        let t =
            Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(fieldcraft_model::WallModel::SpaldingBounce) };
        let p = les(None, Some(t))(collide(steps as f64, 1.0, None, Collision::Regularized));
        let mut setup = Setup::new(&model(1.0, nu, dom, p), None, &[], &fieldcraft_lbm::Limits::default()).unwrap();
        assert!((setup.tau_plus - tau).abs() < 1e-12);
        // The wall treatment alone: no subgrid model.
        setup.smagorinsky = 0.0;
        let mut flow = Flow::new(&setup, &[]).unwrap();
        assert!(!flow.turbulent() && flow.on_node_bounce());
        let grid = *flow.grid();
        let (first, last) = (1, ny - 2);
        flow.set_equilibrium(|cell| {
            let [i, j, _] = grid.coords(cell);
            let y = (j as f64).min(2.0 * h - j as f64);
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            let s = if j == first || j == last { seed * sign } else { 0.0 };
            (1.0, [profile(y) + s, 0.0, 0.0])
        });
        let Some(mut gpu) = gpu(&flow, streaming) else { return };
        let mut cpu = flow;
        let (mut on_gpu, mut on_cpu) = (Vec::with_capacity(steps as usize), Vec::with_capacity(steps as usize));
        for _ in 0..steps {
            gpu.advance(1).unwrap();
            cpu.advance(1).unwrap();
            for (log_e, e) in [
                (&mut on_gpu, odd_even_energy(&grid, &gpu.fields().unwrap().velocity, first..=last)),
                (&mut on_cpu, odd_even_energy(&grid, &cpu.fields().unwrap().velocity, first..=last)),
            ] {
                assert!(e.is_finite() && e > 0.0, "{e}");
                log_e.push(0.5 * e.ln());
            }
        }
        let (rate, cpu_rate) = (tail_slope(&on_gpu), tail_slope(&on_cpu));
        println!(
            "on-node, regularised, τ̄ {tau}, u⁺ {up}, seed {seed:e}, {streaming:?}: growth {rate:+.4e} per step on the GPU, {cpu_rate:+.4e} on the CPU"
        );
        assert!(rate <= 0.0 && cpu_rate <= 0.0, "{rate} {cpu_rate}");
        assert!((rate - cpu_rate).abs() <= 0.1 * cpu_rate.abs(), "{rate} vs {cpu_rate}");
    }
}

/// Han, Ooka and Kikumoto's wall-function bounce (`spalding_bounce_halfway`) with Smagorinsky's model and
/// Van Driest's damping, as their channel runs, in the coarse channel for 150 steps on both, every
/// collision and both layouts: the fields within the wall-model test's tolerances, and the walls'
/// drag (their sampled stresses; the reflection carries none along the walls) within 1e-3 of the
/// CPU's.
#[test]
fn a_coarse_channel_with_the_wall_function_bounce_follows_the_cpu() {
    let bounce = |p: Procedure| {
        let mut p = van_driest(p, 0.1, false);
        if let Procedure::Fluid { turbulence: Some(fieldcraft_model::Turbulence::Smagorinsky { wall_model, .. }), .. } = &mut p {
            *wall_model = Some(fieldcraft_model::WallModel::SpaldingBounceHalfway);
        }
        p
    };
    for (collision, streaming) in cases() {
        let (setup, flow) = coarse_channel(bounce, collision);
        assert!(flow.wall_bounce() && flow.van_driest());
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let wall = |f: &[[f64; 3]]| f[BoxFace::YMin.index()][0] + f[BoxFace::YMax.index()][0];
        let stress: f64 = cpu.walls().unwrap().stresses().map(|s| s.3[0]).sum();
        let drag = (wall(&fc) - wall(&fg)).abs() / wall(&fc).abs();
        println!(
            "{collision:?} {streaming:?} channel with the wall-function bounce, 150 steps, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}, wall drag {drag:.2e} (CPU drag over its sampled stress {:.15})",
            wall(&fc) / stress
        );
        assert!((wall(&fc) / stress - 1.0).abs() < 1e-12, "the CPU's walls carry their stress");
        assert!(du < 1e-4 && drho < 1e-5, "{du} {drho}");
        assert!(drag < 1e-3, "{drag}");
    }
}

/// A laminar channel under the wall model (8 cells across, body force, periodic along it) run to
/// steady state on both: the wall drag, the bounce-back's exchange plus what the reconstruction
/// takes, carries the body force on the GPU as on the CPU, within the force tests' 5e-5.
#[test]
fn a_steady_wall_modelled_channel_s_drag_matches_the_cpu() {
    for (collision, streaming) in cases() {
        let (n, h, nu) = (8usize, 1.0, 0.02);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        d.acceleration = [1e-3, 0.0, 0.0];
        let dt = 0.25 * dx;
        let m = model(1.0, nu, d, wall_modelled(collide(5.0 * h * h / nu, dt, None, collision), 0.1));
        let (setup, flow) = prepare(&m, &[]);
        assert!(flow.wall_model());
        let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { continue };
        // The reconstruction's density does not conserve mass: the body force on the final mass.
        let body = cpu.acceleration()[0] * cpu.mass();
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        let wall = |f: &[[f64; 3]]| f[BoxFace::YMin.index()][0] + f[BoxFace::YMax.index()][0];
        let (dc, dg) = (wall(&fc), wall(&fg));
        println!(
            "{collision:?} steady wall-modelled channel ({streaming:?}): drag / body force CPU {:.9} GPU {:.9}, GPU vs CPU {:.2e}",
            dc / body,
            dg / body,
            (dg - dc).abs() / dc.abs()
        );
        assert!((dc / body - 1.0).abs() < 1e-9, "{dc} vs {body}");
        assert!((dg - dc).abs() <= 5e-5 * dc.abs(), "{dg} vs {dc}");
    }
}

// ---- The hybrid recursive regularised collision, Vreman's model and dynamic σ ----

/// `p` with the hybrid collision's σ `sigma` and the subgrid model `turbulence`.
fn les(sigma: Option<HrrSigma>, turbulence: Option<Turbulence>) -> impl Fn(Procedure) -> Procedure + Copy {
    move |mut p| {
        if let Procedure::Fluid { sigma: s, turbulence: t, .. } = &mut p {
            *s = sigma;
            *t = turbulence;
        }
        p
    }
}

const VREMAN: Turbulence = Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: None };
const SMAGORINSKY: Turbulence = Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None };

/// A periodic box of 24³ cells with a ball of radius 6 at its centre (an obstacle: the gradient's
/// one-sided differences at its surface), a body force along x, ν 0.005 (τ 0.515), started from a
/// three-dimensional field of lattice speed 0.03, 200 steps.
fn ball_box(turbulence: impl Fn(Procedure) -> Procedure, collision: Collision) -> (Setup, Flow) {
    use std::f64::consts::PI;
    let l = 24usize;
    let side = l as f64;
    let mut d = domain([0.0; 3], [side; 3], 1.0);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [1e-5, 0.0, 0.0];
    d.obstacles = vec!["Ball".into()];
    let ball = obstacle("Ball", sphere([side / 2.0; 3], 6.0, 3), &["skin"]);
    let mut p = collide(200.0, 1.0, None, collision);
    if let Procedure::Fluid { end_time, steps, frames, .. } = &mut p {
        (*end_time, *steps, *frames) = (None, Some(200), 1);
    }
    let m = model(1.0, 0.005, d, turbulence(p));
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&ball));
    let grid = *flow.grid();
    let k = 2.0 * PI / side;
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        (1.0, [0.03 * (k * y).sin() * (k * z).cos(), 0.03 * (k * z).sin() * (k * x).cos(), 0.03 * (k * x).sin() * (k * y).cos()])
    });
    (setup, flow)
}

/// The CPU and the GPU through `setup`'s steps from `flow`: the velocity (over the largest speed)
/// and density differences, and the largest eddy viscosity ratio on each.
fn les_difference(setup: &Setup, flow: Flow, streaming: Streaming) -> Option<(f64, f64, f64, f64)> {
    let (mut cpu, mut gpu) = both(setup, flow, streaming)?;
    let (du, drho) = field_difference(&mut cpu, &mut gpu);
    Some((du, drho, Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio()))
}

/// The stated tolerances of the LES cases below, those of the coarse wall-modelled channel:
/// velocity within 1e-4 of the largest speed, density within 1e-5, the largest ν_t/ν within 1e-3
/// of itself.
fn check_les(what: &str, (du, drho, rc, rg): (f64, f64, f64, f64)) {
    let dr = if rc > 0.0 { (rg - rc).abs() / rc } else { rg.abs() };
    println!("{what}, GPU vs CPU: velocity {du:.2e}, density {drho:.2e}; ν_t/ν CPU {rc:.5} GPU {rg:.5} ({dr:.1e})");
    assert!(du < 1e-4 && drho < 1e-5, "{what}: {du} {drho}");
    assert!(dr < 1e-3, "{what}: ν_t/ν {rg} vs {rc}");
}

/// Every field's bits, velocity then density.
fn field_bits(f: &fieldcraft_lbm::backend::Fields) -> Vec<u64> {
    f.velocity.iter().flatten().chain(&f.density).map(|v| v.to_bits()).collect()
}

/// Both cases (the forced channel at friction Reynolds number `re_tau`, and the box with a ball)
/// with `turbulence` and `collision`, both layouts, against the CPU.
fn les_cases(name: &str, re_tau: f64, turbulence: impl Fn(Procedure) -> Procedure + Copy, collision: Collision) {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (setup, flow) = channel_at(re_tau, turbulence, collision);
        let Some(r) = les_difference(&setup, flow, streaming) else { return };
        check_les(&format!("{name}, {collision:?}, {streaming:?}, channel at Re_τ {re_tau}, 150 steps"), r);
        let (setup, flow) = ball_box(turbulence, collision);
        let Some(r) = les_difference(&setup, flow, streaming) else { return };
        check_les(&format!("{name}, {collision:?}, {streaming:?}, box with a ball, 200 steps"), r);
    }
}

/// The hybrid collision (σ 0.9, alone and with Smagorinsky's model) in the coarse channel and the
/// box with a ball, both layouts, against the CPU at [`check_les`]'s tolerances.
#[test]
fn the_hybrid_collision_follows_the_cpu() {
    // Without a subgrid model the channel runs at Re_τ 30 (180 needs one).
    les_cases("hrr σ 0.9", 30.0, les(Some(HrrSigma::Value(0.9)), None), Collision::Hrr);
    les_cases("hrr σ 0.9 + Smagorinsky 0.1", 180.0, les(Some(HrrSigma::Value(0.9)), Some(SMAGORINSKY)), Collision::Hrr);
}

/// σ = 1 on the GPU is its recursive regularised collision bit for bit: a fixed σ of 1 compiles to
/// the regularised kernels (with the hybrid kernels, which skip the blend there, the driver rounds
/// the rest differently: 1.5e-8 of a lattice speed of 0.03 on the GB10).
#[test]
fn sigma_one_is_the_regularised_collision_on_the_gpu() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let run = |collision: Collision, sigma: Option<HrrSigma>| {
            let (setup, flow) = ball_box(les(sigma, None), collision);
            let mut g = gpu(&flow, streaming)?;
            g.advance(setup.steps).unwrap();
            Some(g.fields().unwrap())
        };
        let Some(rr) = run(Collision::Regularized, None) else { return };
        let Some(hrr) = run(Collision::Hrr, Some(HrrSigma::Value(1.0))) else { return };
        let same = field_bits(&rr) == field_bits(&hrr);
        println!("{streaming:?}: hrr σ 1 and regularized on the GPU bit for bit: {same}");
        assert!(same);
    }
}

/// Vreman's model with every collision that takes it, in both cases and layouts, against the CPU.
#[test]
fn vreman_s_model_follows_the_cpu() {
    for collision in [Collision::Trt, Collision::Regularized, Collision::Mrt] {
        les_cases("Vreman 0.18", 180.0, les(None, Some(VREMAN)), collision);
    }
    les_cases("Vreman 0.18", 180.0, les(Some(HrrSigma::Value(0.9)), Some(VREMAN)), Collision::Hrr);
}

/// The hybrid collision's dynamic σ, from Vreman's and from Smagorinsky's eddy viscosity, in both
/// cases and layouts, against the CPU.
#[test]
fn dynamic_sigma_follows_the_cpu() {
    les_cases("dynamic σ, Vreman 0.18", 180.0, les(Some(HrrSigma::Dynamic), Some(VREMAN)), Collision::Hrr);
    les_cases("dynamic σ, Smagorinsky 0.1", 180.0, les(Some(HrrSigma::Dynamic), Some(SMAGORINSKY)), Collision::Hrr);
}

/// The hybrid collision (σ 0.9) with Smagorinsky's model and Van Driest's damping, and with the
/// wall-function bounce (`spalding_bounce_halfway`, damping too), in the coarse channel, both layouts:
/// the fields and ν_t/ν at [`check_les`]'s tolerances, the bounce's wall drag within 1e-3.
#[test]
fn the_hybrid_collision_with_wall_damping_and_the_wall_function_bounce_follows_the_cpu() {
    for bounce in [false, true] {
        let turbulence = move |p: Procedure| {
            let mut p = van_driest(p, 0.1, false);
            if let Procedure::Fluid { sigma, turbulence: Some(Turbulence::Smagorinsky { wall_model, .. }), .. } = &mut p {
                *sigma = Some(HrrSigma::Value(0.9));
                *wall_model = bounce.then_some(fieldcraft_model::WallModel::SpaldingBounceHalfway);
            }
            p
        };
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let (setup, flow) = coarse_channel(turbulence, Collision::Hrr);
            assert!(flow.van_driest() && flow.wall_bounce() == bounce);
            let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { return };
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
            let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
            let wall = |f: &[[f64; 3]]| f[BoxFace::YMin.index()][0] + f[BoxFace::YMax.index()][0];
            let drag = (wall(&fc) - wall(&fg)).abs() / wall(&fc).abs();
            let what = if bounce { "the wall-function bounce" } else { "Van Driest's damping" };
            check_les(&format!("hrr σ 0.9 + Smagorinsky with {what}, {streaming:?}, channel, 150 steps"), (du, drho, rc, rg));
            println!("    wall drag {drag:.2e}");
            if bounce {
                assert!(drag < 1e-3, "{drag}");
            }
        }
    }
}

/// Vreman's model under the wall-function bounce on the hybrid collision, at σ 0.99 and with
/// dynamic σ (the Han et al. channel's `hrr-vreman` and `hrr-dynamic` set-ups), in the coarse
/// channel, both layouts: the fields and ν_t/ν at [`check_les`]'s tolerances, the bounce's wall
/// drag within 1e-3.
#[test]
fn vreman_s_model_with_the_wall_function_bounce_follows_the_cpu() {
    for han in [HanSetup::HrrVreman, HanSetup::HrrDynamic] {
        let (sigma, t) = (han.sigma(), han.turbulence(true));
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let (setup, flow) = coarse_channel(les(sigma, Some(t)), han.collision());
            assert!(flow.wall_bounce() && flow.vreman() > 0.0);
            let Some((mut cpu, mut gpu)) = both(&setup, flow, streaming) else { return };
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
            let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
            let wall = |f: &[[f64; 3]]| f[BoxFace::YMin.index()][0] + f[BoxFace::YMax.index()][0];
            let drag = (wall(&fc) - wall(&fg)).abs() / wall(&fc).abs();
            check_les(&format!("{} with the wall-function bounce, {streaming:?}, channel, 150 steps", han.name()), (du, drho, rc, rg));
            println!("    wall drag {drag:.2e}");
            assert!(rc > 0.0 && drag < 1e-3, "{rc} {drag}");
        }
    }
}

/// The same run twice on one device gives the same fields bit for bit (dynamic σ with Vreman's
/// model in the box with a ball, and the hybrid collision with Smagorinsky's in the channel).
#[test]
fn the_les_kernels_are_deterministic_on_their_device() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        for (name, case) in [
            ("dynamic σ + Vreman, box with a ball", ball_box(les(Some(HrrSigma::Dynamic), Some(VREMAN)), Collision::Hrr)),
            ("hrr σ 0.9 + Smagorinsky, channel", coarse_channel(les(Some(HrrSigma::Value(0.9)), Some(SMAGORINSKY)), Collision::Hrr)),
        ] {
            let (setup, flow) = case;
            let run = || {
                let mut g = gpu(&flow, streaming)?;
                g.advance(setup.steps).unwrap();
                Some(field_bits(&g.fields().unwrap()))
            };
            let (Some(a), Some(b)) = (run(), run()) else { return };
            println!("{name}, {streaming:?}: two runs bit for bit: {}", a == b);
            assert!(a == b, "{name}");
        }
    }
}

/// The CPU's shear-wave gate on the GPU (`crates/lbm/tests/hrr.rs`): N 16, τ̄ 0.8, σ 0.9, the wave
/// u_x = 1e-4 sin(k y) in a periodic box of 2 × 16 × 2 cells decays at the exact rate of the
/// linear analysis within 1e-3 (f32), both layouts (`exact_rate`, the 2 × 2 map of the momentum and
/// the xy flux derived in `common`).
#[test]
fn a_shear_wave_decays_at_the_exact_rate_on_the_gpu() {
    let (n, tau, sigma, amplitude) = (16usize, 0.8, 0.9, 1e-4);
    assert!((exact_rate(std::f64::consts::PI / 4.0, 0.8, 0.5) / 6.4112e-2 - 1.0).abs() < 1e-4, "the CPU test's independent value");
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let mut d = domain([0.0; 3], [2.0, n as f64, 2.0], 1.0);
        for f in BoxFace::ALL {
            *d.faces.get_mut(f) = PERIODIC;
        }
        let p = les(Some(HrrSigma::Value(sigma)), None)(collide(1.0, 1.0, None, Collision::Hrr));
        let (setup, mut flow) = prepare(&model(1.0, (tau - 0.5) / 3.0, d, p), &[]);
        assert!((setup.tau_plus - tau).abs() < 1e-12);
        let k = 2.0 * std::f64::consts::PI / n as f64;
        let grid = *flow.grid();
        let mode = |cell: usize| (k * grid.coords(cell)[1] as f64).sin();
        flow.set_equilibrium(|cell| (1.0, [amplitude * mode(cell), 0.0, 0.0]));
        let Some(mut g) = gpu(&flow, streaming) else { return };
        let measure = |g: &mut Gpu| {
            let f = g.fields().unwrap();
            let (num, den) = (0..grid.len()).fold((0.0, 0.0), |(a, b), cell| (a + f.velocity[cell][0] * mode(cell), b + mode(cell) * mode(cell)));
            num / den
        };
        let exact = exact_rate(k, tau, sigma);
        let (t0, t1) = (200u64, 200 + (2.0 / exact).round() as u64);
        g.advance(t0).unwrap();
        let a0 = measure(&mut g);
        g.advance(t1 - t0).unwrap();
        let a1 = measure(&mut g);
        let measured = -(a1 / a0).ln() / (t1 - t0) as f64;
        println!("{streaming:?} shear wave N {n}, τ̄ {tau}, σ {sigma}: measured {measured:.6e}, exact {exact:.6e} ({:+.2e})", measured / exact - 1.0);
        assert!((measured / exact - 1.0).abs() < 1e-3, "{measured} vs {exact}");
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

// ---- The wall-modelled channel against Lee & Moser (#664) ----

/// Lee & Moser's mean velocity U⁺ at y/δ (M. Lee, R. D. Moser, *Direct numerical simulation of
/// turbulent channel flow up to Re_τ ≈ 5200*, J. Fluid Mech. 774:395–415, 2015,
/// doi:10.1017/jfm.2015.268; sampled from the mean profiles of turbulence.oden.utexas.edu,
/// `LM_Channel_0550_mean_prof.dat`, `LM_Channel_1000_mean_prof.dat` and
/// `LM_Channel_2000_mean_prof.dat`), their Re_τ and bulk
/// U_b⁺ = 1 / u_τ (the files' U_mean = 1).
struct LeeMoser {
    re_tau: f64,
    bulk: f64,
    profile: [(f64, f64); 16],
}

const LM_550: LeeMoser = LeeMoser {
    re_tau: 543.496,
    bulk: 1.0 / 5.43496e-2,
    profile: [
        (0.01, 5.203),
        (0.02, 8.917),
        (0.03, 11.062),
        (0.05, 13.163),
        (0.07, 14.178),
        (0.1, 15.068),
        (0.15, 16.015),
        (0.2, 16.714),
        (0.3, 17.781),
        (0.4, 18.608),
        (0.5, 19.288),
        (0.6, 19.856),
        (0.7, 20.326),
        (0.8, 20.688),
        (0.9, 20.919),
        (1.0, 21.000),
    ],
};

const LM_1000: LeeMoser = LeeMoser {
    re_tau: 1000.512,
    bulk: 1.0 / 5.00256e-2,
    profile: [
        (0.01, 8.437),
        (0.02, 11.968),
        (0.03, 13.455),
        (0.05, 14.839),
        (0.07, 15.624),
        (0.1, 16.456),
        (0.15, 17.450),
        (0.2, 18.200),
        (0.3, 19.320),
        (0.4, 20.168),
        (0.5, 20.858),
        (0.6, 21.434),
        (0.7, 21.909),
        (0.8, 22.272),
        (0.9, 22.509),
        (1.0, 22.593),
    ],
};

const LM_2000: LeeMoser = LeeMoser {
    re_tau: 1994.756,
    bulk: 1.0 / 4.58794e-2,
    profile: [
        (0.01, 11.924),
        (0.02, 14.241),
        (0.03, 15.233),
        (0.05, 16.417),
        (0.07, 17.228),
        (0.1, 18.129),
        (0.15, 19.195),
        (0.2, 19.976),
        (0.3, 21.121),
        (0.4, 21.971),
        (0.5, 22.659),
        (0.6, 23.237),
        (0.7, 23.716),
        (0.8, 24.081),
        (0.9, 24.315),
        (1.0, 24.396),
    ],
};

impl LeeMoser {
    /// U⁺ at `y` (over δ, in (0, 1]), interpolated linearly in ln y.
    fn at(&self, y: f64) -> f64 {
        let p = &self.profile;
        let i = p.windows(2).position(|w| y <= w[1].0).unwrap_or(p.len() - 2);
        let ((y0, u0), (y1, u1)) = (p[i], p[i + 1]);
        let t = ((y.max(y0)).ln() - y0.ln()) / (y1.ln() - y0.ln());
        u0 + t.clamp(0.0, 1.0) * (u1 - u0)
    }
}

/// Lee & Moser interpolated to Malaspinas and Sagaut's Re_τ 950: linearly in Re_τ between the
/// 550 and 1000 channels, point by point in y/δ and for the bulk speed; shown beside their own
/// reference ([`HJ_950`]).
fn lm_950() -> LeeMoser {
    let w = (950.0 - LM_550.re_tau) / (LM_1000.re_tau - LM_550.re_tau);
    let mut profile = LM_1000.profile;
    for (p, (a, b)) in profile.iter_mut().zip(LM_550.profile.iter().zip(&LM_1000.profile)) {
        p.1 = a.1 + w * (b.1 - a.1);
    }
    LeeMoser { re_tau: 950.0, bulk: LM_550.bulk + w * (LM_1000.bulk - LM_550.bulk), profile }
}

/// The mean profile Malaspinas and Sagaut compare with at Re_τ 950: Hoyas and Jiménez's
/// statistics (S. Hoyas, J. Jiménez, *Reynolds number effects on the Reynolds-stress budgets in
/// turbulent channels*, Phys. Fluids 20:101511, 2008, doi:10.1063/1.3005862) of del Álamo and
/// Jiménez's channel at Re_τ 934 (J. C. del Álamo, J. Jiménez, *Spectra of the very large
/// anisotropic scales in turbulent channels*, Phys. Fluids 15:L41–L44, 2003,
/// doi:10.1063/1.1570830), from the file `Re950.prof` of the Torroja DNS database
/// (https://torroja.dmt.upm.es/channels/data/): (y/δ, U⁺) at every row from y/δ 0.01 to the
/// centre. Re_τ 934 against the nominal 950 is the paper's own mismatch, kept.
#[rustfmt::skip]
const HJ_950: [(f64, f64); 175] = [
    (0.010823, 8.4979), (0.012057, 9.1019), (0.013357, 9.6639), (0.014722, 10.1824),
    (0.016154, 10.6580), (0.017651, 11.0922), (0.019215, 11.4877), (0.020844, 11.8476),
    (0.022538, 12.1752), (0.024298, 12.4737), (0.026123, 12.7462), (0.028013, 12.9956),
    (0.029969, 13.2245), (0.031989, 13.4353), (0.034074, 13.6302), (0.036224, 13.8112),
    (0.038438, 13.9799), (0.040717, 14.1379), (0.043060, 14.2866), (0.045467, 14.4270),
    (0.047937, 14.5604), (0.050472, 14.6874), (0.053070, 14.8090), (0.055731, 14.9258),
    (0.058456, 15.0385), (0.061244, 15.1474), (0.064094, 15.2530), (0.067007, 15.3557),
    (0.069983, 15.4558), (0.073021, 15.5536), (0.076120, 15.6492), (0.079282, 15.7429),
    (0.082505, 15.8349), (0.085790, 15.9253), (0.089136, 16.0142), (0.092543, 16.1016),
    (0.096011, 16.1878), (0.099539, 16.2727), (0.103127, 16.3565), (0.106776, 16.4392),
    (0.110484, 16.5209), (0.114252, 16.6016), (0.118079, 16.6812), (0.121965, 16.7599),
    (0.125910, 16.8377), (0.129913, 16.9146), (0.133975, 16.9906), (0.138094, 17.0659),
    (0.142271, 17.1404), (0.146506, 17.2141), (0.150798, 17.2871), (0.155146, 17.3593),
    (0.159552, 17.4309), (0.164013, 17.5017), (0.168530, 17.5719), (0.173103, 17.6414),
    (0.177732, 17.7103), (0.182415, 17.7786), (0.187153, 17.8463), (0.191946, 17.9133),
    (0.196792, 17.9797), (0.201693, 18.0455), (0.206647, 18.1108), (0.211654, 18.1754),
    (0.216713, 18.2396), (0.221825, 18.3033), (0.226990, 18.3665), (0.232206, 18.4292),
    (0.237473, 18.4914), (0.242791, 18.5531), (0.248160, 18.6143), (0.253580, 18.6750),
    (0.259049, 18.7352), (0.264568, 18.7949), (0.270136, 18.8541), (0.275753, 18.9128),
    (0.281418, 18.9711), (0.287132, 19.0289), (0.292893, 19.0864), (0.298702, 19.1434),
    (0.304557, 19.2001), (0.310460, 19.2563), (0.316408, 19.3122), (0.322402, 19.3676),
    (0.328441, 19.4226), (0.334525, 19.4772), (0.340654, 19.5312), (0.346827, 19.5848),
    (0.353044, 19.6378), (0.359304, 19.6905), (0.365607, 19.7427), (0.371952, 19.7945),
    (0.378339, 19.8460), (0.384768, 19.8970), (0.391239, 19.9478), (0.397750, 19.9981),
    (0.404301, 20.0480), (0.410892, 20.0975), (0.417522, 20.1466), (0.424192, 20.1954),
    (0.430900, 20.2437), (0.437646, 20.2917), (0.444430, 20.3393), (0.451251, 20.3866),
    (0.458108, 20.4335), (0.465002, 20.4800), (0.471932, 20.5261), (0.478897, 20.5718),
    (0.485897, 20.6171), (0.492932, 20.6620), (0.500000, 20.7066), (0.507102, 20.7507),
    (0.514237, 20.7945), (0.521404, 20.8379), (0.528603, 20.8809), (0.535834, 20.9236),
    (0.543096, 20.9659), (0.550389, 21.0078), (0.557711, 21.0494), (0.565064, 21.0907),
    (0.572445, 21.1317), (0.579855, 21.1723), (0.587293, 21.2125), (0.594759, 21.2523),
    (0.602252, 21.2918), (0.609771, 21.3309), (0.617317, 21.3695), (0.624888, 21.4077),
    (0.632484, 21.4454), (0.640105, 21.4826), (0.647750, 21.5194), (0.655419, 21.5557),
    (0.663110, 21.5916), (0.670824, 21.6269), (0.678561, 21.6617), (0.686318, 21.6960),
    (0.694097, 21.7299), (0.701896, 21.7632), (0.709715, 21.7959), (0.717554, 21.8281),
    (0.725411, 21.8596), (0.733287, 21.8904), (0.741181, 21.9206), (0.749092, 21.9501),
    (0.757020, 21.9790), (0.764964, 22.0072), (0.772924, 22.0347), (0.780899, 22.0614),
    (0.788888, 22.0875), (0.796892, 22.1128), (0.804910, 22.1373), (0.812940, 22.1610),
    (0.820983, 22.1839), (0.829038, 22.2060), (0.837105, 22.2273), (0.845182, 22.2476),
    (0.853270, 22.2671), (0.861367, 22.2856), (0.869474, 22.3032), (0.877589, 22.3199),
    (0.885713, 22.3356), (0.893844, 22.3503), (0.901983, 22.3640), (0.910128, 22.3767),
    (0.918279, 22.3884), (0.926436, 22.3990), (0.934597, 22.4085), (0.942763, 22.4169),
    (0.950932, 22.4243), (0.959105, 22.4305), (0.967281, 22.4356), (0.975459, 22.4396),
    (0.983638, 22.4424), (0.991819, 22.4441), (1.000000, 22.4447),
];

/// The bulk U⁺ of [`HJ_950`]'s file, by the trapezoid rule over its rows: its C_f = 2 / U_b⁺².
const HJ_950_BULK: f64 = 19.827;

/// [`HJ_950`] at `y` (δ = 1), linear between its rows, held beyond its ends.
fn hj_950(y: f64) -> f64 {
    let i = HJ_950.partition_point(|p| p.0 < y).clamp(1, HJ_950.len() - 1);
    let ((y0, u0), (y1, u1)) = (HJ_950[i - 1], HJ_950[i]);
    u0 + (u1 - u0) * ((y - y0) / (y1 - y0)).clamp(0.0, 1.0)
}

/// Malaspinas and Sagaut's average error of the mean profile (their eq. 69), in per cent:
/// E = 100 / M · √(Σ (ũ⁺ − ũ⁺_ref)² / ũ⁺_ref²). The paper does not say which points it sums or
/// which u_τ scales ũ⁺; the reading that reproduces its Table 1 is taken: `half` is the profile
/// folded about the centreline, scaled by the carried u_τ, counted on both halves (M = 2N).
fn profile_error(half: &[(f64, f64)], reference: impl Fn(f64) -> f64) -> f64 {
    let m = 2.0 * half.len().max(1) as f64;
    100.0 / m * (2.0 * half.iter().map(|&(y, u)| ((u - reference(y)) / reference(y)).powi(2)).sum::<f64>()).sqrt()
}

/// One grid of the study: what it measured against Lee & Moser.
struct ChannelResult {
    bulk: f64,
    centre: f64,
    /// The largest relative error of U⁺ from the second cell to y/δ 0.3, and the first cell's.
    log_error: f64,
    first_error: f64,
    cells: usize,
    steps: u64,
    seconds: f64,
    stop: String,
    /// The wall model's requested stress over the carried one (NaN where not measured).
    requested: f64,
    /// C_f from the wall model's time-averaged u_τ, 2 (⟨u_τ⟩ / U_b)², and Malaspinas and
    /// Sagaut's profile error E (NaN where not measured).
    cf_model: f64,
    e_paper: f64,
}

/// The wall-modelled channel at `lm`'s Re_τ on `n` cells per δ, on the GPU: 2πδ × 2δ × πδ,
/// periodic, driven by u_τ²/δ (δ = u_τ = 1), regularised collision, C = 0.1, the Spalding wall
/// model; started from a blunt profile at the reference's bulk speed with streaks to trip it,
/// run under the running-mean stop rule on the mass flow (10 δ/u_τ discarded, windows of
/// 2 δ/u_τ, 2e-3, at most 60 δ/u_τ), the plane-averaged profile summed from the start-up on.
fn wall_modelled_channel(lm: &LeeMoser, n: usize) -> Option<ChannelResult> {
    use fieldcraft_lbm::steady::Monitor;
    use std::f64::consts::PI;
    let (start, window, end, tolerance) = (10.0, 2.0, 60.0, 2e-3);
    let ub = lm.bulk;
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    // The protocol's collision is the regularised one; FIELDCRAFT_WALL_MODEL_COLLISION=trt runs
    // the two-relaxation-time one instead, as a diagnostic.
    let collision = match std::env::var("FIELDCRAFT_WALL_MODEL_COLLISION").as_deref() {
        Ok("trt") => Collision::Trt,
        Ok("mrt") => Collision::Mrt,
        _ => Collision::Regularized,
    };
    let mut p = wall_modelled(collide(end, dt, None, collision), 0.1);
    if let Procedure::Fluid { steady, .. } = &mut p {
        *steady = Some(fieldcraft_model::Steady {
            monitors: vec![fieldcraft_model::FlowMonitor::MassFlow { face: BoxFace::XMax }],
            tolerance,
            window: Some(window),
            min_time: None,
            divergence: None,
            averaging: Some(fieldcraft_model::Averaging { start }),
        });
    }
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / lm.re_tau, d, p), &[]);
    let grid = *flow.grid();
    let vs = setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0;
        let wall = 1.0 - eta * eta;
        let u = 1.15 * ub * (1.0 - eta.abs().powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u / vs, v / vs, w / vs])
    });
    let mut gpu = gpu(&flow, Streaming::InPlace)?;
    let mut monitor = Monitor::new(&setup, &*gpu).unwrap().unwrap();
    let ny = grid.cells[1];
    let from = (start / dt).ceil() as u64;
    let every = ((0.25 / dt) as u64).max(1);
    let (mut profile, mut samples) = (vec![0.0; ny], 0usize);
    let clock = std::time::Instant::now();
    let mut k = 0u64;
    let mut stop = "the step's end".to_string();
    while k < setup.steps {
        let mut next = setup.steps.min((k / every + 1) * every);
        next = next.min((k / monitor.setup.sample + 1) * monitor.setup.sample);
        if let Err(e) = gpu.advance(next - k) {
            println!("  Re_τ {} {n}/δ: unstable before {:.2} δ/u_τ: {e}", lm.re_tau.round(), next as f64 * dt);
            let seconds = clock.elapsed().as_secs_f64();
            return Some(ChannelResult {
                bulk: f64::NAN,
                centre: f64::NAN,
                log_error: f64::NAN,
                first_error: f64::NAN,
                cells: grid.len(),
                steps: next,
                seconds,
                stop: format!("unstable before {:.2} δ/u_τ", next as f64 * dt),
                requested: f64::NAN,
                cf_model: f64::NAN,
                e_paper: f64::NAN,
            });
        }
        k = next;
        if k >= from && k.is_multiple_of(every) {
            let fields = gpu.fields().unwrap();
            for (cell, u) in fields.velocity.iter().enumerate() {
                profile[grid.coords(cell)[1]] += u[0] * vs;
            }
            samples += 1;
        }
        if monitor.samples_at(k)
            && let Some(reason) = monitor.sample(&setup, &mut *gpu, k).unwrap()
        {
            stop = format!("{} at {:.1} δ/u_τ", reason.name(), k as f64 * dt);
            break;
        }
        if k % (setup.steps / 20).max(1) < every {
            println!("  Re_τ {} {n}/δ: t = {:.1} δ/u_τ, {:.0} s", lm.re_tau.round(), k as f64 * dt, clock.elapsed().as_secs_f64());
        }
    }
    let seconds = clock.elapsed().as_secs_f64();
    let report = monitor.report();
    let mass = report.averages.iter().zip(&report.quantities).find(|(_, q)| q.contains("mass")).map(|(a, _)| a.mean[0]);
    let area = 2.0 * lz;
    let bulk = mass.map_or(f64::NAN, |m| m.abs() / area);
    let (nx, nz) = (grid.cells[0], grid.cells[2]);
    let mean: Vec<f64> = profile.iter().map(|s| s / (samples * nx * nz) as f64).collect();
    // Symmetric about the centre: both halves at each distance from the nearer wall.
    let half: Vec<(f64, f64)> = (0..ny / 2).map(|j| ((j as f64 + 0.5) * dx, 0.5 * (mean[j] + mean[ny - 1 - j]))).collect();
    let error = |(y, u): (f64, f64)| (u / lm.at(y) - 1.0).abs();
    let first_error = half.first().map_or(f64::NAN, |c| error(*c));
    let log_error = half.iter().skip(1).filter(|(y, _)| *y <= 0.3).map(|c| error(*c)).fold(0.0, f64::max);
    let centre = half.last().map_or(f64::NAN, |c| c.1);
    for (y, u) in &half {
        println!("    y/δ {y:.4} y⁺ {:7.1}: U⁺ {u:7.3}, Lee & Moser {:7.3} ({:+.1} %)", y * lm.re_tau, lm.at(*y), 100.0 * (u / lm.at(*y) - 1.0));
    }
    Some(ChannelResult {
        bulk,
        centre,
        log_error,
        first_error,
        cells: grid.len(),
        steps: k,
        seconds,
        stop,
        requested: f64::NAN,
        cf_model: f64::NAN,
        e_paper: f64::NAN,
    })
}

/// The environment variable that runs the wall-modelled channel study on the GPU.
const WALL_MODEL_STUDY: &str = "FIELDCRAFT_WALL_MODEL_STUDY";

/// The wall-modelled channel at Re_τ 550 (10, 14, 18 cells per δ) and 1000 (16, 24, 32), first
/// cell sizes Δy⁺ 30 to 63, against Lee & Moser, on the GPU only (hours): docs/verification.md,
/// "Long studies". Acceptance, set before the runs: C_f within 5 % on every grid and the log
/// layer's mean velocity within 5 % (from the second cell to y/δ 0.3). Set
/// `FIELDCRAFT_WALL_MODEL_STUDY` to `1` for every grid, or to `550` or `1000`, or to
/// `1000:16` for one grid.
#[test]
#[ignore = "hours, GPU: set FIELDCRAFT_WALL_MODEL_STUDY=1 (docs/verification.md, \"Long studies\")"]
fn study_a_wall_modelled_channel_against_lee_and_moser() {
    let Some(which) = std::env::var_os(WALL_MODEL_STUDY).map(|v| v.to_string_lossy().into_owned()) else {
        println!("the wall-modelled channel study runs only with {WALL_MODEL_STUDY} set (hours, GPU): skipped");
        return;
    };
    let grids: [(&LeeMoser, &str, usize); 6] =
        [(&LM_550, "550", 10), (&LM_550, "550", 14), (&LM_550, "550", 18), (&LM_1000, "1000", 16), (&LM_1000, "1000", 24), (&LM_1000, "1000", 32)];
    let mut failed = Vec::new();
    for (lm, name, n) in grids {
        let chosen = which == "1" || which == name || which == format!("{name}:{n}");
        if !chosen {
            continue;
        }
        let Some(r) = wall_modelled_channel(lm, n) else { return };
        let cf_error = (lm.bulk / r.bulk).powi(2) - 1.0;
        let dy_plus = lm.re_tau / n as f64;
        println!(
            "wall-modelled channel Re_τ {name}, {n}/δ (Δy⁺ {dy_plus:.0}, first cell y⁺ {:.0}): U_b⁺ {:.3} (Lee & Moser {:.2}), C_f {:.4e} ({:+.1} %), U_c⁺ {:.2} (LM {:.2}), log layer within {:.1} %, first cell {:.1} %; {} cells, {} steps, {:.0} s, {:.0} MLUPS; stopped: {}",
            dy_plus / 2.0,
            r.bulk,
            lm.bulk,
            2.0 / (r.bulk * r.bulk),
            100.0 * cf_error,
            r.centre,
            lm.at(1.0),
            100.0 * r.log_error,
            100.0 * r.first_error,
            r.cells,
            r.steps,
            r.seconds,
            r.cells as f64 * r.steps as f64 / r.seconds / 1e6,
            r.stop
        );
        if !(cf_error.abs() <= 0.05 && r.log_error <= 0.05) {
            failed.push(format!("Re_τ {name} {n}/δ"));
        }
    }
    assert!(failed.is_empty(), "outside the acceptance: {failed:?}");
}

/// The wall model's own consistency in the turbulent channel (Re_τ 1000, 16 cells per δ, the
/// study's set-up), on demand: after 20 δ/u_τ, 20 samples a quarter δ/u_τ apart of the plane
/// averages at the second cells (U₂), the first cells as the fields read them after streaming
/// (U₁), the u_τ and target the wall law gives from the second cells' velocities on the host, and
/// the targets the GPU's `wall_sample` computed from the same state (read back one step later).
#[test]
#[ignore = "a diagnostic: set FIELDCRAFT_WALL_MODEL_DIAGNOSTIC=1"]
fn diagnose_the_wall_model_in_the_channel() {
    use std::f64::consts::PI;
    if std::env::var_os("FIELDCRAFT_WALL_MODEL_DIAGNOSTIC").is_none() {
        return;
    }
    let (lm, n) = (&LM_1000, 16usize);
    let ub = lm.bulk;
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let collision = match std::env::var("FIELDCRAFT_WALL_MODEL_DIAGNOSTIC").as_deref() {
        Ok("trt") => Collision::Trt,
        Ok("mrt") => Collision::Mrt,
        _ => Collision::Regularized,
    };
    let p = wall_modelled(collide(30.0, dt, None, collision), 0.1);
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / lm.re_tau, d, p), &[]);
    let grid = *flow.grid();
    let vs = setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    println!("{collision:?}");
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0;
        let wall = 1.0 - eta * eta;
        let u = 1.15 * ub * (1.0 - eta.abs().powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u / vs, v / vs, w / vs])
    });
    let walls: Vec<(usize, usize)> = flow.walls().unwrap().wall_cells().filter_map(|(cell, _, _, _, second)| second.map(|s| (cell, s))).collect();
    let mut gpu = gpu(&flow, Streaming::InPlace).unwrap();
    let nu = (setup.tau_plus - 0.5) / 3.0;
    gpu.advance((20.0 / dt) as u64).unwrap();
    let every = (0.25 / dt) as u64;
    let [nx, ny, nz] = grid.cells;
    let plane = |fields: &fieldcraft_lbm::backend::Fields, j: usize| {
        let mut s = 0.0;
        for k in 0..nz {
            for i in 0..nx {
                s += fields.velocity[grid.index([i, j, k])][0];
            }
        }
        s * vs / (nx * nz) as f64
    };
    let mut sums = [0.0f64; 6];
    let samples = 20;
    // Per plane from the wall (both halves): Σ u, Σ v (towards the centre), Σ u v, Σ u², samples.
    let half = ny / 2;
    let mut stats = vec![[0.0f64; 4]; half];
    let mut cells_per_plane = 0.0;
    for _ in 0..samples {
        gpu.advance(every - 1).unwrap();
        let fields = gpu.fields().unwrap();
        let (u2, u1) = (0.5 * (plane(&fields, 1) + plane(&fields, ny - 2)), 0.5 * (plane(&fields, 0) + plane(&fields, ny - 1)));
        // The host's wall law at each wall cell from this state, averaged over the wall cells.
        let (mut ut, mut target) = (0.0, 0.0);
        for &(_, second) in &walls {
            let u = fields.velocity[second];
            let speed = (u[0] * u[0] + u[2] * u[2]).sqrt();
            let (u1_law, _) = fieldcraft_lbm::wall::wall_law(speed, 0.5, 1.5, nu);
            ut += fieldcraft_lbm::wall::friction_velocity(speed, 1.5, nu) * vs;
            target += u1_law * if speed > 0.0 { u[0] / speed } else { 0.0 } * vs;
        }
        let m = walls.len() as f64;
        for (j, st) in stats.iter_mut().enumerate() {
            for (row, sign) in [(j, 1.0), (ny - 1 - j, -1.0)] {
                for k in 0..nz {
                    for i in 0..nx {
                        let u = fields.velocity[grid.index([i, row, k])];
                        let (a, b) = (u[0] * vs, sign * u[1] * vs);
                        st[0] += a;
                        st[1] += b;
                        st[2] += a * b;
                        st[3] += a * a;
                    }
                }
            }
        }
        cells_per_plane += 2.0 * (nx * nz) as f64;
        gpu.advance(1).unwrap();
        let gpu_target = gpu.wall_targets().unwrap().iter().map(|t| t.0[0]).sum::<f64>() / m * vs;
        let row = [u2, u1, ut / m, target / m, gpu_target, 0.0];
        for (s, v) in sums.iter_mut().zip(row) {
            *s += v;
        }
        println!(
            "  U₂⁺ {:.3}  U₁⁺ (fields) {:.3}  u_τ (law, host) {:.3}  target U₁⁺ (host) {:.3}  target U₁⁺ (GPU) {:.3}",
            row[0], row[1], row[2], row[3], row[4]
        );
    }
    let s = samples as f64;
    println!(
        "mean over {samples}: U₂⁺ {:.3} (LM at y⁺ {:.0}: {:.3}), U₁⁺ (fields) {:.3}, u_τ (law) {:.3}, target U₁⁺ host {:.3} / GPU {:.3} (LM at y⁺ {:.0}: {:.3})",
        sums[0] / s,
        1.5 * dx * lm.re_tau,
        lm.at(1.5 * dx),
        sums[1] / s,
        sums[2] / s,
        sums[3] / s,
        sums[4] / s,
        0.5 * dx * lm.re_tau,
        lm.at(0.5 * dx)
    );
    // The resolved shear stress −u'v' against the total the channel carries, 1 − y/δ (u_τ = 1):
    // what the resolved motion leaves to the viscous and subgrid stresses.
    for (j, st) in stats.iter().enumerate().take(6) {
        let y = (j as f64 + 0.5) * dx;
        let (mu, mv) = (st[0] / cells_per_plane, st[1] / cells_per_plane);
        let uv = st[2] / cells_per_plane - mu * mv;
        let urms = (st[3] / cells_per_plane - mu * mu).max(0.0).sqrt();
        println!("    y⁺ {:6.1}: U⁺ {mu:7.3}, −u'v'⁺ {:.3} of the total {:.3}, u'⁺ rms {urms:.3}", y * lm.re_tau, -uv, 1.0 - y);
    }
}

/// Non-overlapping block means of a statistically steady signal sampled at a fixed interval: the
/// mean over the completed blocks and its 95 % interval, 2 s / √B with s the block means' standard
/// deviation and B their number, which holds while the blocks are longer than the signal's
/// correlation time (the block means' lag-1 correlation, reported, stays small then). Each sample
/// counts once, unlike a mean of running means, which weighs the early samples more and can make
/// a drift look settled (bluecmd's review on #664).
#[derive(Debug, Clone, Default)]
struct Blocks {
    len: usize,
    sum: f64,
    filled: usize,
    means: Vec<f64>,
}

impl Blocks {
    fn new(len: usize) -> Blocks {
        Blocks { len: len.max(1), ..Blocks::default() }
    }

    fn push(&mut self, v: f64) {
        self.sum += v;
        self.filled += 1;
        if self.filled == self.len {
            self.means.push(self.sum / self.len as f64);
            (self.sum, self.filled) = (0.0, 0);
        }
    }

    /// The mean over the completed blocks (NaN before the first).
    fn mean(&self) -> f64 {
        self.means.iter().sum::<f64>() / self.means.len() as f64
    }

    /// The half-width of the mean's 95 % interval (infinite with fewer than two blocks).
    fn half_width(&self) -> f64 {
        let b = self.means.len();
        if b < 2 {
            return f64::INFINITY;
        }
        let m = self.mean();
        let var = self.means.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (b - 1) as f64;
        2.0 * (var / b as f64).sqrt()
    }

    /// The lag-1 correlation of the block means (0 with fewer than three blocks).
    fn lag1(&self) -> f64 {
        let b = self.means.len();
        if b < 3 {
            return 0.0;
        }
        let m = self.mean();
        let var: f64 = self.means.iter().map(|x| (x - m) * (x - m)).sum();
        let cov: f64 = self.means.windows(2).map(|w| (w[0] - m) * (w[1] - m)).sum();
        if var > 0.0 { cov / var } else { 0.0 }
    }
}

/// The block means of a known sequence: each sample once, the partial block left out, and the
/// interval of independent blocks.
#[test]
fn block_means_count_each_sample_once() {
    let mut b = Blocks::new(4);
    for v in [1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 100.0] {
        b.push(v);
    }
    assert_eq!(b.means, vec![1.0, 3.0]);
    assert_eq!(b.mean(), 2.0);
    // s = √2, B = 2: 2 √2 / √2 = 2.
    assert!((b.half_width() - 2.0).abs() < 1e-12);
    assert_eq!(Blocks::new(4).half_width(), f64::INFINITY);
    // A steady signal plus a slow drift: the mean of running means lags the plain mean.
    let (mut blocks, mut running, mut sum) = (Blocks::new(10), 0.0, 0.0);
    for k in 1..=200 {
        let v = 1.0 + 1e-3 * k as f64;
        blocks.push(v);
        sum += v;
        running += sum / k as f64;
    }
    let plain = sum / 200.0;
    assert!((blocks.mean() - plain).abs() < 1e-12);
    assert!(running / 200.0 < plain - 0.03, "the old statistic lags: {} vs {plain}", running / 200.0);
    // Alternating block means are anti-correlated.
    let mut alt = Blocks::new(1);
    for k in 0..20 {
        alt.push(if k % 2 == 0 { 1.0 } else { -1.0 });
    }
    assert!(alt.lag1() < -0.9);
}

/// A channel's stress budget and the wall model's requested stress, from the fields sampled
/// (lattice units in, model units out: δ = 1, ρ = 1). Per wall-normal row: the mean u, v, uv and
/// an estimate of the subgrid shear ν_t (∂u/∂y + ∂v/∂x) with Smagorinsky's ν_t = (C Δ)² |S| from
/// the resolved strain by central differences (the solver takes |S| from the non-equilibrium
/// moments instead, so this is an estimate); and the stress the wall model requests,
/// u_τ² with u_τ from Musker's profile at each second cell's speed, as the model computes it.
struct Budget {
    /// Per row: Σu, Σv, Σuv, Σ subgrid shear.
    rows: Vec<[f64; 4]>,
    /// Σ requested u_τ² and Σ u_τ over the wall cells, and how many.
    requested: f64,
    friction: f64,
    walls: usize,
    samples: usize,
}

impl Budget {
    fn new(ny: usize) -> Budget {
        Budget { rows: vec![[0.0; 4]; ny], requested: 0.0, friction: 0.0, walls: 0, samples: 0 }
    }

    /// Add one sample of `velocity` (lattice units) on `grid`, periodic in x and z, walls at
    /// both y ends; `nu` the lattice viscosity, `c` Smagorinsky's constant.
    fn sample(&mut self, grid: &fieldcraft_lbm::voxel::Grid, velocity: &[[f64; 3]], nu: f64, c: f64) {
        let [nx, ny, nz] = grid.cells;
        let at = |i: usize, j: usize, k: usize| velocity.get(grid.index([i % nx, j, k % nz])).copied().unwrap_or([0.0; 3]);
        for k in 0..nz {
            for i in 0..nx {
                for j in 0..ny {
                    let u = at(i, j, k);
                    let row = &mut self.rows[j];
                    row[0] += u[0];
                    row[1] += u[1];
                    row[2] += u[0] * u[1];
                    if j == 0 || j + 1 == ny {
                        continue;
                    }
                    // The velocity gradient by central differences (lattice spacing 1).
                    let (xp, xm, yp, ym, zp, zm) =
                        (at(i + 1, j, k), at(i + nx - 1, j, k), at(i, j + 1, k), at(i, j - 1, k), at(i, j, k + 1), at(i, j, k + nz - 1));
                    let d = |a: [f64; 3], b: [f64; 3]| [0, 1, 2].map(|q| 0.5 * (a[q] - b[q]));
                    let g = [d(xp, xm), d(yp, ym), d(zp, zm)]; // g[a][q] = ∂u_q/∂x_a
                    let ss: f64 = (0..9)
                        .map(|ab| {
                            let sab = 0.5 * (g[ab / 3][ab % 3] + g[ab % 3][ab / 3]);
                            sab * sab
                        })
                        .sum();
                    let nu_t = c * c * (2.0 * ss).sqrt();
                    row[3] += nu_t * (g[1][0] + g[0][1]);
                }
            }
            // The requested stress on both walls: the second cell's speed parallel to the wall.
            for i in 0..nx {
                for (second, y1) in [(1usize, 0.5), (ny - 2, 0.5)] {
                    let u = at(i, second, k);
                    let speed = (u[0] * u[0] + u[2] * u[2]).sqrt();
                    let ut = fieldcraft_lbm::wall::musker_friction(speed, 3.0 * y1, nu);
                    self.requested += ut * ut;
                    self.friction += ut;
                    self.walls += 1;
                }
            }
        }
        self.samples += 1;
    }

    /// Print the budget's first rows off the walls (both halves folded) against the total the
    /// channel carries, g (1 − y) with g the mean force, and the requested over the carried stress.
    /// `vs` turns lattice speeds into model ones, `nu` is the lattice viscosity, `dx` the spacing.
    fn report(&self, grid: &fieldcraft_lbm::voxel::Grid, vs: f64, nu: f64, dx: f64, g: f64, re_tau: f64) -> f64 {
        let [nx, ny, nz] = grid.cells;
        let per = (self.samples * nx * nz).max(1) as f64;
        let mean = |j: usize, q: usize| self.rows[j][q] / per;
        let requested = self.requested / self.walls.max(1) as f64 * vs * vs;
        println!("  stress at the wall: requested by the model ⟨u_τ²⟩ {requested:.5}, carried (mean force) {g:.5}, ratio {:.4}", requested / g);
        // Lattice stresses to model ones: a lattice stress (per unit density) times vs².
        for j in 1..(ny / 2).min(8) {
            let y = (j as f64 + 0.5) * dx;
            let fold = |q: usize| 0.5 * (mean(j, q) + if q == 0 || q == 2 { mean(ny - 1 - j, q) } else { -mean(ny - 1 - j, q) });
            let (u, v, uv) = (fold(0), fold(1), 0.5 * (mean(j, 2) - mean(ny - 1 - j, 2)));
            let dudy = 0.5 * (0.5 * (mean(j + 1, 0) - mean(j - 1, 0)) + 0.5 * (mean(ny - 2 - j, 0) - mean(ny - j, 0)));
            let viscous = nu * dudy * vs * vs;
            let resolved = -(uv - u * v) * vs * vs;
            let sgs = 0.5 * (mean(j, 3) - mean(ny - 1 - j, 3)) * vs * vs;
            let total = g * (1.0 - y);
            println!(
                "    y⁺ {:6.1}: total {total:.4} = viscous {viscous:.4} + subgrid (estimate) {sgs:.4} + resolved {resolved:.4} + residual {:.4}",
                y * re_tau,
                total - viscous - sgs - resolved
            );
        }
        requested / g
    }
}

/// The bulk speed a channel study holds: Dean's correlation at the nominal Re_τ, as Malaspinas
/// and Sagaut do (the paper reproduction), or Lee & Moser's, the DNS the acceptance judges against
/// (the control: no difference of reference conditions inside the 5 %).
#[derive(Debug, Clone, Copy, PartialEq)]
enum Held {
    Dean,
    LeeMoser,
}

/// Malaspinas and Sagaut's channel set-up (J. Comput. Phys. 275:25–40, 2014, §5) at `lm`'s
/// bulk Reynolds number on `n` cells per δ, on the GPU: 2πδ × 2δ × 2πδ, periodic, the D3Q19
/// multiple-relaxation-time collision, Smagorinsky's model at C = 0.08 and the wall model; the
/// flow rate held by the driving force g = u_τ²/δ + (U_m − ⟨u⟩) U_m / δ, updated every step (their
/// eq. 67), U_m = Lee & Moser's bulk speed (δ = 1, the
/// nominal u_τ = 1). u_τ and C_f come from the mean force (the walls carry it all: ρ u_τ² = ρ ḡ δ).
/// The stop rule is on the force's block means ([`Blocks`]): 20 δ/u_τ discarded, blocks of 2 δ/u_τ,
/// at least 50 δ/u_τ averaged, then stopped once the mean force's 95 % interval is within 1 % of it
/// and the block means' lag-1 correlation is below 0.5; at most 150 δ/u_τ. `held` picks the bulk
/// speed held ([`Held`]).
fn paper_channel(lm: &LeeMoser, n: usize, c: f64, cpu: bool, held: Held) -> Option<(ChannelResult, f64)> {
    use std::f64::consts::PI;
    // The paper does not state its averaging.
    let (start, min_average, block, end, tolerance) = (20.0, 50.0, 2.0, 150.0, 1e-2);
    // The target bulk speed, with δ = 1 and the nominal u_τ = 1: from Dean's correlation at the
    // nominal Re_τ (their eq. 68: Re_m = 2 U_m / ν, ν = 1 / Re_τ), or Lee & Moser's.
    let um = match held {
        Held::Dean => dean_bulk(lm.re_tau),
        Held::LeeMoser => lm.bulk,
    };
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (2.0 * PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / um;
    let p = wall_modelled(collide(end, dt, None, Collision::Mrt), c);
    let m = model(1.0, 1.0 / lm.re_tau, d, p);
    // A long study: past the default work cap, which guards interactive runs.
    let limits = fieldcraft_lbm::Limits { max_cell_updates: 10_000_000_000_000, ..fieldcraft_lbm::Limits::default() };
    let setup = match Setup::new(&m, None, &[], &limits) {
        Ok(s) => s,
        Err(e) => {
            println!("  M&S set-up Re_τ {} {n}/δ: refused: {e}", lm.re_tau.round());
            return Some((
                ChannelResult {
                    bulk: f64::NAN,
                    centre: f64::NAN,
                    log_error: f64::NAN,
                    first_error: f64::NAN,
                    cells: 0,
                    steps: 0,
                    seconds: 0.0,
                    stop: format!("refused: {e}"),
                    requested: f64::NAN,
                    cf_model: f64::NAN,
                    e_paper: f64::NAN,
                },
                f64::NAN,
            ));
        }
    };
    let mut flow = Flow::new(&setup, &[]).unwrap();
    let grid = *flow.grid();
    let vs = setup.velocity_scale();
    let lattice_acceleration = setup.acceleration[0];
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0;
        let wall = 1.0 - eta * eta;
        let u = 1.15 * um * (1.0 - eta.abs().powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u / vs, v / vs, w / vs])
    });
    // On the GPU, or (`cpu`) on the CPU in double precision, the GPU's reference.
    let mut device = if cpu { None } else { Some(gpu(&flow, Streaming::InPlace)?) };
    let gpu: &mut dyn Backend = match device.as_mut() {
        Some(g) => &mut **g,
        None => &mut flow,
    };
    // The profile and the stop rule's samples every 0.05 δ/U_m; the force every step.
    let control = ((0.05 / dt) as u64).max(1);
    let from = (start / dt) as u64;
    // The force every step after the start-up, in blocks.
    let mut forces = Blocks::new(((block / dt) as usize).max(1));
    let [nx, ny, nz] = grid.cells;
    let (mut profile, mut samples) = (vec![0.0; ny], 0usize);
    let (mut ub_sum, mut count) = (0.0, 0usize);
    let mut budget = Budget::new(ny);
    let mut skipped = 0usize;
    let lattice_nu = (setup.tau_plus - 0.5) / 3.0;
    let clock = std::time::Instant::now();
    let (mut k, mut g) = (0u64, 1.0);
    let mut stop = "the step's end".to_string();
    while k < setup.steps {
        let step = control.min(setup.steps - k);
        let mut failed = None;
        for _ in 0..step {
            if let Err(e) = gpu.advance(1) {
                failed = Some(e);
                break;
            }
            // Their eq. 67 (δ = 1): the nominal u_τ² plus the flow-rate error times U_m / δ.
            let mean = gpu.mean_velocity(0).unwrap() * vs;
            g = 1.0 + (um - mean) * um;
            gpu.set_acceleration([g * lattice_acceleration, 0.0, 0.0]).unwrap();
            k += 1;
            if k >= from {
                forces.push(g);
            }
        }
        k -= step;
        if let Some(e) = failed {
            println!("  M&S set-up Re_τ {} {n}/δ: unstable before {:.2} δ/u_τ: {e}", lm.re_tau.round(), (k + step) as f64 * dt);
            return Some((
                ChannelResult {
                    bulk: f64::NAN,
                    centre: f64::NAN,
                    log_error: f64::NAN,
                    first_error: f64::NAN,
                    cells: grid.len(),
                    steps: k + step,
                    seconds: clock.elapsed().as_secs_f64(),
                    stop: "unstable".into(),
                    requested: f64::NAN,
                    cf_model: f64::NAN,
                    e_paper: f64::NAN,
                },
                f64::NAN,
            ));
        }
        k += step;
        // The bulk speed from the device's own sum, and the fields only when two reads agree:
        // a read of the GPU's fields during a long run has come back corrupt (a few cells with
        // impossible speeds) a few times in a thousand, while the device's sum stayed right
        // (work item on solver-LES). The skipped samples are counted.
        let mean = gpu.mean_velocity(0).unwrap() * vs;
        let fields = (0..3).find_map(|_| {
            let (a, b) = (gpu.fields().ok()?, gpu.fields().ok()?);
            (a.velocity == b.velocity && a.velocity.iter().all(|u| u.iter().all(|v| v.abs() < 1.0))).then_some(a)
        });
        if fields.is_none() {
            skipped += 1;
        }
        if k >= from {
            ub_sum += mean;
            count += 1;
            if let Some(fields) = &fields {
                for (cell, u) in fields.velocity.iter().enumerate() {
                    profile[grid.coords(cell)[1]] += u[0] * vs;
                }
                samples += 1;
                budget.sample(&grid, &fields.velocity, lattice_nu, c);
            }
            // Steady once the mean force is known to the tolerance from uncorrelated blocks.
            if k as f64 * dt >= start + min_average && forces.half_width() < tolerance * forces.mean() && forces.lag1() < 0.5 {
                stop = format!("steady at {:.1} δ/u_τ", k as f64 * dt);
                break;
            }
        }
        if k % (setup.steps / 20).max(1) < control {
            println!(
                "  M&S set-up Re_τ {} {n}/δ: t = {:.1}, U_b {mean:.3}, g {g:.4}, {:.0} s",
                lm.re_tau.round(),
                k as f64 * dt,
                clock.elapsed().as_secs_f64()
            );
        }
    }
    // The force averaged over the completed blocks after the start-up, and its interval.
    let g_mean = forces.mean();
    let ut = g_mean.sqrt();
    println!(
        "  force: {:.5} ± {:.5} (95 %, {} blocks of {block} δ/u_τ, lag-1 correlation {:.2}); C_f ± {:.2} %",
        g_mean,
        forces.half_width(),
        forces.means.len(),
        forces.lag1(),
        100.0 * forces.half_width() / g_mean
    );
    let requested = budget.report(&grid, vs, lattice_nu, dx, g_mean, g_mean.sqrt() * lm.re_tau);
    println!("  field samples: {samples} averaged, {skipped} skipped (two reads disagreed or a speed was impossible)");
    let ub = ub_sum / count.max(1) as f64 / ut;
    let mean: Vec<f64> = profile.iter().map(|s| s / (samples * nx * nz) as f64 / ut).collect();
    let half: Vec<(f64, f64)> = (0..ny / 2).map(|j| ((j as f64 + 0.5) * dx, 0.5 * (mean[j] + mean[ny - 1 - j]))).collect();
    let error = |(y, u): (f64, f64)| (u / lm.at(y) - 1.0).abs();
    let log_error = half.iter().skip(1).filter(|(y, _)| *y <= 0.3).map(|c| error(*c)).fold(0.0, f64::max);
    for (y, u) in &half {
        println!("    y/δ {y:.4} y⁺ {:7.1}: U⁺ {u:7.3}, Lee & Moser {:7.3} ({:+.1} %)", y * ut * lm.re_tau, lm.at(*y), 100.0 * (u / lm.at(*y) - 1.0));
    }
    // Malaspinas and Sagaut's comparison (their eqs. 69 and 70): C_f = 2 u_τ² / U_b², with u_τ the
    // wall model's own time average (the second cell's reading, both walls) or the carried force's,
    // against Dean's correlation at the held Re_m and the DNS's bulk C_f; E against the DNS (Hoyas
    // and Jiménez at Re_τ 950, else Lee & Moser) and Lee & Moser.
    let ub_phys = ub_sum / count.max(1) as f64;
    let ut_model = budget.friction / budget.walls.max(1) as f64 * vs;
    let cf_model = 2.0 * (ut_model / ub_phys).powi(2);
    let cf_force = 2.0 * g_mean / (ub_phys * ub_phys);
    let cf_dean = 0.073 * (2.0 * ub_phys * lm.re_tau).powf(-0.25);
    let hoyas = lm.re_tau == 950.0;
    let cf_dns = 2.0 / if hoyas { HJ_950_BULK } else { lm.bulk }.powi(2);
    let e_lm = profile_error(&half, |y| lm.at(y));
    let e_paper = if hoyas { profile_error(&half, hj_950) } else { e_lm };
    println!(
        "  M&S comparison: ⟨u_τ⟩ model {ut_model:.4}, carried {ut:.4} (gap {:+.2} %); C_f model {cf_model:.5} ({:+.2} % of Dean {cf_dean:.5}, {:+.2} % of the DNS {cf_dns:.5}), carried {cf_force:.5} ({:+.2} %, {:+.2} %); E {e_paper:.2} % against {} ({e_lm:.2} % against Lee & Moser)",
        100.0 * (ut_model / ut - 1.0),
        100.0 * (cf_model / cf_dean - 1.0),
        100.0 * (cf_model / cf_dns - 1.0),
        100.0 * (cf_force / cf_dean - 1.0),
        100.0 * (cf_force / cf_dns - 1.0),
        if hoyas { "Hoyas & Jiménez" } else { "Lee & Moser" }
    );
    let result = ChannelResult {
        bulk: ub,
        centre: half.last().map_or(f64::NAN, |c| c.1),
        log_error,
        first_error: half.first().map_or(f64::NAN, |c| error(*c)),
        cells: grid.len(),
        steps: k,
        seconds: clock.elapsed().as_secs_f64(),
        stop,
        requested,
        cf_model,
        e_paper,
    };
    Some((result, ut * lm.re_tau))
}

/// Dean's bulk speed U_m / u_τ at friction Reynolds number `re_tau` (Malaspinas and Sagaut's eq.
/// 68: Re_m = (8 / 0.073)^(4/7) Re_τ^(8/7), Re_m = 2 U_m δ / ν).
fn dean_bulk(re_tau: f64) -> f64 {
    (8.0f64 / 0.073).powf(4.0 / 7.0) * re_tau.powf(8.0 / 7.0) / (2.0 * re_tau)
}

/// Malaspinas and Sagaut's set-up (see [`paper_channel`]) on their grids, 10, 20 and 30 cells per
/// δ, at Re_τ 1000 (their 950; Lee & Moser's nearest, C = 0.08 as theirs) and 2000 (C = 0.1 as
/// theirs), with the study's acceptance against Lee & Moser (C_f and the log layer within 5 %);
/// C_f against Dean's correlation is printed too, the paper's own comparison. Re_τ 550 (`550`)
/// is run only when asked, outside the acceptance (the paper tests Re_τ ≥ 950). Set
/// `FIELDCRAFT_PAPER_CHANNEL` to `1` (1000 and 2000), `1000`, `2000`, `550`, or `1000:10` for one
/// grid. `950` (or `950:10`) reproduces the paper itself at its Re_τ 950: E (their eq. 69, read as
/// [`profile_error`] does) against Hoyas and Jiménez ([`HJ_950`]) within 15 % of their Table 1
/// (2.11, 1.07 and 0.83 on 10, 20 and 30 cells), and C_f from the wall model's u_τ (their eq. 70)
/// within the 4.22 % of Dean's correlation they report; measured on the GB10: E 2.08, 0.99 and
/// 0.87, C_f +1.8, +0.8 and +3.5 %, while the carried drag's C_f is 14 to 17 % below Dean's (the
/// wall model's u_τ is 8.1 to 10.6 % above the carried one). The bulk speed held is Dean's, the
/// paper's; `FIELDCRAFT_PAPER_CHANNEL_HELD=lm` holds Lee & Moser's instead, the control against
/// the DNS the acceptance uses. Each grid reports its first
/// cell's centre and the wall model's matching point (the second cell, 3 y₁) in wall units.
#[test]
#[ignore = "hours, GPU: set FIELDCRAFT_PAPER_CHANNEL=1 (docs/verification.md, \"Long studies\")"]
fn study_malaspinas_and_sagaut_s_channel_against_lee_and_moser() {
    let Some(which) = std::env::var_os("FIELDCRAFT_PAPER_CHANNEL").map(|v| v.to_string_lossy().into_owned()) else {
        println!("the paper's channel set-up runs only with FIELDCRAFT_PAPER_CHANNEL set (hours, GPU): skipped");
        return;
    };
    let lm950 = lm_950();
    let grids: [(&LeeMoser, &str, usize, f64, bool); 12] = [
        (&lm950, "950", 10, 0.08, false),
        (&lm950, "950", 20, 0.08, false),
        (&lm950, "950", 30, 0.08, false),
        (&LM_1000, "1000", 10, 0.08, true),
        (&LM_1000, "1000", 20, 0.08, true),
        (&LM_1000, "1000", 30, 0.08, true),
        (&LM_2000, "2000", 10, 0.1, true),
        (&LM_2000, "2000", 20, 0.1, true),
        (&LM_2000, "2000", 30, 0.1, true),
        (&LM_550, "550", 10, 0.08, false),
        (&LM_550, "550", 14, 0.08, false),
        (&LM_550, "550", 18, 0.08, false),
    ];
    let mut failed = Vec::new();
    for (lm, name, n, c, counted) in grids {
        if !((which == "1" && counted) || which == name || which == format!("{name}:{n}")) {
            continue;
        }
        // FIELDCRAFT_PAPER_CHANNEL_CPU=1 runs the grid on the CPU in double precision instead.
        let cpu = std::env::var_os("FIELDCRAFT_PAPER_CHANNEL_CPU").is_some();
        let held = if std::env::var("FIELDCRAFT_PAPER_CHANNEL_HELD").is_ok_and(|v| v == "lm") { Held::LeeMoser } else { Held::Dean };
        let Some((r, re_tau)) = paper_channel(lm, n, c, cpu, held) else { return };
        let cf = 2.0 / (r.bulk * r.bulk);
        let cf_lm = 2.0 / (lm.bulk * lm.bulk);
        let cf_error = cf / cf_lm - 1.0;
        // Dean's C_f at the bulk Reynolds number held (nominal u_τ = δ = 1): C_f = 0.073 Re_m^(−1/4).
        let um = if held == Held::Dean { dean_bulk(lm.re_tau) } else { lm.bulk };
        let cf_dean = 0.073 * (2.0 * um * lm.re_tau).powf(-0.25);
        println!(
            "M&S set-up ({}, {held:?} bulk held), Lee & Moser Re_τ {name}, {n}/δ, C {c} (Δy⁺ {:.0}, first cell centre y⁺ {:.0}, matching point y⁺ {:.0}): measured Re_τ {re_tau:.0}, U_b⁺ {:.3} (LM {:.2}), C_f {cf:.4e} ({:+.1} % on LM, {:+.1} % on Dean), U_c⁺ {:.2} (LM {:.2}), log layer within {:.1} %, requested / carried wall stress {:.3}; {} cells, {} steps, {:.0} s, {:.0} MLUPS; stopped: {}",
            if cpu { "CPU, f64" } else { "GPU, f32" },
            re_tau / n as f64,
            0.5 * re_tau / n as f64,
            1.5 * re_tau / n as f64,
            r.bulk,
            lm.bulk,
            100.0 * cf_error,
            100.0 * (cf / cf_dean - 1.0),
            r.centre,
            lm.at(1.0),
            100.0 * r.log_error,
            r.requested,
            r.cells,
            r.steps,
            r.seconds,
            r.cells as f64 * r.steps as f64 / r.seconds / 1e6,
            r.stop
        );
        println!("M&S eqs. 69/70, Re_τ {name} {n}/δ: C_f by the model's ⟨u_τ⟩ {:.4e}, E {:.2} %", r.cf_model, r.e_paper);
        // The reproduction of the paper's own results at Re_τ 950 (Dean's bulk held, as theirs).
        if name == "950" && held == Held::Dean {
            let table = match n {
                10 => 2.11,
                20 => 1.07,
                _ => 0.83,
            };
            if !((r.e_paper / table - 1.0).abs() <= 0.15 && (r.cf_model / cf_dean - 1.0).abs() <= 0.0422) {
                failed.push(format!(
                    "Re_τ 950 {n}/δ against the paper: E {:.2} (theirs {table}), C_f {:+.1} % of Dean",
                    r.e_paper,
                    100.0 * (r.cf_model / cf_dean - 1.0)
                ));
            }
        }
        // Counted where the method is meant to work: the first cell's centre at y⁺ ≥ 30 for the
        // nominal Re_τ (a grid the τ floor refuses is reported, not counted).
        let in_range = 0.5 * lm.re_tau / n as f64 >= 30.0;
        if counted && in_range && !r.stop.starts_with("refused") && !(cf_error.abs() <= 0.05 && r.log_error <= 0.05) {
            failed.push(format!("Re_τ {name} {n}/δ"));
        }
    }
    assert!(failed.is_empty(), "outside the acceptance (first cell at y⁺ ≥ 30, or the paper's own at Re_τ 950): {failed:?}");
}

/// One of Han, Ooka and Kikumoto's channel cases (Fluid Dyn. Res. 53:045506, 2021, §3, Table 1):
/// Re_τ, cells per half-height, the paper's Mach number and its centreline Reynolds number U_c D / ν
/// (which sets the lattice speed), and its results for the wall-function bounce (WFB) and the
/// bounce-back (BB): the L2 error ε_u⁺ against Spalding's law in the buffer layer, the log layer
/// and the whole channel (Table 2), and the log law's intercept b (Table 3).
struct HanCase {
    name: &'static str,
    re_tau: f64,
    n: usize,
    mach: f64,
    re_centre: f64,
    /// Whether the case may run at the smallest Mach number that keeps τ at the floor when the
    /// paper's falls under it (a labelled deviation; B40 only, by the owner's decision).
    raise_to_floor: bool,
    wfb: ([f64; 3], f64),
    bb: ([f64; 3], f64),
}

const HAN_CASES: [HanCase; 4] = [
    HanCase {
        name: "A20",
        raise_to_floor: false,
        re_tau: 640.0,
        n: 20,
        mach: 0.060,
        re_centre: 13800.0,
        wfb: ([0.069, 0.055, 0.057], 6.4),
        bb: ([0.183, 0.009, 0.056], 5.6),
    },
    HanCase {
        name: "A40",
        raise_to_floor: false,
        re_tau: 640.0,
        n: 40,
        mach: 0.060,
        re_centre: 13800.0,
        wfb: ([0.071, 0.067, 0.067], 6.6),
        bb: ([0.103, 0.039, 0.048], 6.2),
    },
    HanCase {
        name: "B40",
        raise_to_floor: true,
        re_tau: 2003.0,
        n: 40,
        mach: 0.069,
        re_centre: 48000.0,
        wfb: ([0.049, 0.080, 0.079], 6.6),
        bb: ([0.911, 0.420, 0.429], -0.1),
    },
    HanCase {
        name: "B80",
        raise_to_floor: false,
        re_tau: 2003.0,
        n: 80,
        mach: 0.069,
        re_centre: 48000.0,
        wfb: ([0.057, 0.066, 0.065], 6.5),
        bb: ([0.617, 0.224, 0.235], -4.1),
    },
];

/// The bulk U⁺ of Hoyas and Jiménez's channel at Re_τ 2003 (`Re2000.prof` of the Torroja
/// database, the trapezoid rule over its rows), the reference of Han et al. at that Re_τ.
const HJ_2000_BULK: f64 = 21.748;

/// Spalding's u⁺ at `y_plus` with Han et al.'s constants (κ 0.41, B 5.5): the root of their eq.
/// 16 by bisection (y⁺(u⁺) rises with u⁺; u⁺ ≤ y⁺).
fn spalding_u_plus(y_plus: f64) -> f64 {
    let (kappa, b) = (0.41f64, 5.5f64);
    let e = (-kappa * b).exp();
    let y = |u: f64| {
        let k = kappa * u;
        u + e * (k.exp() - 1.0 - k - k * k / 2.0 - k * k * k / 6.0)
    };
    let (mut lo, mut hi) = (0.0, y_plus.max(0.0));
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if y(mid) > y_plus { hi = mid } else { lo = mid }
    }
    0.5 * (lo + hi)
}

/// What one Han et al. case measured.
struct HanResult {
    /// U⁺ at the half-channel's cell centres (y/δ, U⁺), folded about the centreline, by the
    /// nominal u_τ = 1; u′rms⁺ likewise.
    profile: Vec<(f64, f64)>,
    rms: Vec<f64>,
    bulk: f64,
    /// ε_u⁺ in the buffer layer (y⁺ < 100), the log layer (y⁺ ≥ 100) and the whole channel, and b.
    eps: [f64; 3],
    b: f64,
    /// The walls' mean force over the body force on the fluid (1 when the walls carry it).
    carried: f64,
    window: (f64, f64),
    stop: String,
    steps: u64,
    seconds: f64,
    /// The Mach number run and τ.
    mach: f64,
    tau: f64,
    /// The stress budget, eddy viscosity, dissipation and σ per row, and the spectra
    /// ([`ChannelStats`]).
    rows: Vec<HanRow>,
    spectra: Vec<HanSpectrum>,
    /// How many field snapshots the window averaged, and the streamwise cells and length.
    snapshots: usize,
    nx: usize,
    lx: f64,
}

/// The configuration a Han et al. channel runs with (`FIELDCRAFT_HAN_SETUP`). The wall treatment
/// (the wall-function bounce or bounce-back) is chosen per run on top of it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HanSetup {
    /// The paper's: MRT, Smagorinsky's model at C_s 0.1 with Van Driest's damping (A⁺ 26).
    Paper,
    /// J. Jacob, O. Malaspinas, P. Sagaut, J. Turbul. 19:1051, 2018, their explicit LES (case 3):
    /// the hybrid recursive regularised collision at σ 0.99 with Vreman's model at C_s 0.18.
    HrrVreman,
    /// Their case 4: the hybrid collision's dynamic σ from Vreman's ν_t (C_s 0.18), which then
    /// only sets σ (it stays out of τ).
    HrrDynamic,
    /// A labelled extra row: Han et al.'s subgrid model (Smagorinsky 0.1, Van Driest 26) on the
    /// hybrid collision at σ 0.99.
    HrrSmagorinsky,
}

const HAN_SETUPS: [HanSetup; 4] = [HanSetup::Paper, HanSetup::HrrVreman, HanSetup::HrrDynamic, HanSetup::HrrSmagorinsky];

impl HanSetup {
    fn name(self) -> &'static str {
        match self {
            HanSetup::Paper => "paper",
            HanSetup::HrrVreman => "hrr-vreman",
            HanSetup::HrrDynamic => "hrr-dynamic",
            HanSetup::HrrSmagorinsky => "hrr-smagorinsky",
        }
    }

    fn collision(self) -> Collision {
        if self.sigma().is_some() { Collision::Hrr } else { Collision::Mrt }
    }

    fn sigma(self) -> Option<HrrSigma> {
        match self {
            HanSetup::Paper => None,
            HanSetup::HrrDynamic => Some(HrrSigma::Dynamic),
            HanSetup::HrrVreman | HanSetup::HrrSmagorinsky => Some(HrrSigma::Value(0.99)),
        }
    }

    /// The subgrid model, with the wall-function bounce (`wfb`) as its wall model or none.
    fn turbulence(self, wfb: bool) -> Turbulence {
        let wall_model = wfb.then_some(fieldcraft_model::WallModel::SpaldingBounceHalfway);
        match self {
            HanSetup::Paper | HanSetup::HrrSmagorinsky => Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: Some(26.0), wall_model },
            HanSetup::HrrVreman | HanSetup::HrrDynamic => Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model },
        }
    }

    /// Whether the subgrid ν_t enters the cells' relaxation time (not with dynamic σ, where it
    /// only sets σ).
    fn subgrid_applied(self) -> bool {
        self.sigma() != Some(HrrSigma::Dynamic)
    }

    /// Whether a run counts for the gate: the paper's set-up on both walls (as recorded), the
    /// hybrid ones on the wall-function bounce only.
    fn gated(self, wfb: bool) -> bool {
        wfb || self == HanSetup::Paper
    }

    /// What the record says of ν_t's source, and the note printed with the rows.
    fn nu_t_source(self) -> (&'static str, &'static str) {
        match self.turbulence(false) {
            Turbulence::Vreman { .. } => (
                "vreman: fieldcraft_lbm::vreman::eddy_viscosity of fieldcraft_lbm::hrr::gradient on the readback velocity (the model's own formula)",
                "",
            ),
            Turbulence::Smagorinsky { .. } => (
                "smagorinsky: estimate (C D)^2 sqrt(2 S:S) of the finite-difference strain, D Van Driest's at the nominal y+ (the device's is Hou's, from the flux)",
                " (ν_t a finite-difference estimate)",
            ),
        }
    }
}

/// One row of a Han et al. channel, the two halves folded about the centreline (the wall-normal
/// velocity and the shear stresses change sign across it), in wall units of the nominal u_τ = 1.
struct HanRow {
    y: f64,
    u: f64,
    /// u′, v′, w′ rms⁺ (each half's own about its mean, then averaged).
    rms: [f64; 3],
    /// −⟨u′v′⟩⁺, the viscous ⟨∂u⁺/∂y⁺⟩ and the subgrid ⟨2 ν_t S_xy⟩⁺ (from the finite-difference
    /// gradient: [`ChannelStats`]).
    resolved: f64,
    viscous: f64,
    subgrid: f64,
    /// ⟨ν_t⟩/ν, the molecular ⟨2 ν S:S⟩⁺ and the subgrid ⟨2 ν_t S:S⟩⁺ (ε⁺ = ε ν / u_τ⁴).
    nu_t: f64,
    eps_molecular: f64,
    eps_subgrid: f64,
    /// The hybrid collision's σ: the row's mean and smallest (the fixed σ, or 1 without it).
    sigma: (f64, f64),
}

impl HanRow {
    /// The subgrid shear the budget counts (none where ν_t stays out of τ: `applied` false) and
    /// what the total 1 − y/δ leaves after it, the resolved and the viscous stress.
    fn closure(&self, applied: bool) -> (f64, f64) {
        let subgrid = if applied { self.subgrid } else { 0.0 };
        (subgrid, 1.0 - self.y - self.resolved - self.viscous - subgrid)
    }
}

/// u′'s streamwise spectrum at one row (both halves), E_uu⁺ per discrete wavenumber
/// k_x = 2π m / L_x, m = 0 … n_x / 2, one-sided so that Σ_m E = ⟨u′²⟩⁺ about the plane's mean.
struct HanSpectrum {
    row: usize,
    y_plus: f64,
    energy: Vec<f64>,
}

impl HanSpectrum {
    /// The wavenumber index m ≥ 1 of the premultiplied spectrum's peak, k_x E_uu.
    fn peak(&self) -> usize {
        (1..self.energy.len()).fold(1, |best, m| if m as f64 * self.energy[m] > best as f64 * self.energy[best] { m } else { best })
    }

    /// The share of ⟨u′²⟩ (m ≥ 1) in the upper half of the wavenumbers, k_x Δx > π/2.
    fn upper_share(&self, nx: usize) -> f64 {
        let total: f64 = self.energy.iter().skip(1).sum();
        let upper: f64 = self.energy.iter().enumerate().filter(|(m, _)| 4 * m > nx).map(|(_, e)| e).sum();
        upper / total
    }
}

/// Per-row sums of a channel's snapshots ([`ChannelStats`]).
#[derive(Clone, Copy, Default)]
struct RowSums {
    u: [f64; 3],
    uu: [f64; 3],
    uv: f64,
    nu_t: f64,
    eps_molecular: f64,
    eps_subgrid: f64,
    viscous: f64,
    subgrid: f64,
    sigma: f64,
    sigma_min: f64,
}

/// The statistics of a Han et al. channel over its averaging window, from field readbacks of the
/// velocity (lattice units; δ = 1, u_τ = 1 nominal; periodic in x and z, walls at both y ends).
///
/// **Gradient.** The library's own finite differences, [`fieldcraft_lbm::hrr::gradient`]: centred,
/// first-order one-sided at the wall rows (the wall's side is not a fluid cell), periodic images
/// across x and z, the stencil the hybrid collision and Vreman's model use on the device.
///
/// **ν_t.** With Vreman's model, exactly the model's: [`fieldcraft_lbm::vreman::eddy_viscosity`] of
/// that gradient with `c = 2.5 C_s²` ([`fieldcraft_lbm::vreman::model_constant`]), from the readback
/// velocity (the device computes it in f32 from the same velocities). With Smagorinsky's model an
/// **estimate**: the device takes Hou's ν_t from the cell's own non-equilibrium flux, which no
/// readback carries, so here `(C D Δx)² √(2 S:S)` of the finite-difference strain, D Van Driest's
/// `1 − e^{−y⁺/A⁺}` at the nominal y⁺ of the cell centre (the device damps with each wall cell's
/// u_τ from Spalding's law).
///
/// **σ.** Fixed σ as given; dynamic σ by [`fieldcraft_lbm::hrr::dynamic_sigma`] of the same
/// gradient, [`fieldcraft_lbm::hrr::laplacian`], Vreman's ν_t and the molecular ν, with the cap
/// L_max the grid's smallest cell count (as [`fieldcraft_lbm::Flow`] does): the library's formula
/// on the readback field, since the device keeps no per-cell σ.
///
/// **Units.** u⁺ = u v_s, ∂u⁺/∂y⁺ = g v_s / (Δx Re_τ) with g the lattice gradient; ε⁺ = 2 (ν_t/ν)
/// S⁺:S⁺ and 2 S⁺:S⁺; the subgrid shear ⟨2 ν_t S_xy⟩⁺ = (ν_t/ν)(∂u⁺/∂y⁺ + ∂v⁺/∂x⁺).
///
/// **Spectra.** At each snapshot, for each line along x of the rows chosen (and their mirror
/// rows), the discrete Fourier transform of u⁺ less the row's plane mean, |û_m|² / n_x² doubled for
/// 0 < m < n_x/2; averaged over the lines (z and both halves) and the snapshots.
struct ChannelStats {
    grid: fieldcraft_lbm::voxel::Grid,
    vs: f64,
    re_tau: f64,
    dx: f64,
    /// The molecular lattice viscosity (τ − ½)/3 and dynamic σ's cap.
    nu: f64,
    l_max: f64,
    turbulence: Turbulence,
    sigma: Option<HrrSigma>,
    rows: Vec<RowSums>,
    samples: usize,
    /// The rows (lower half) whose spectra are taken, their sums, and sin/cos of 2π q / n_x.
    spectrum_rows: Vec<usize>,
    spectra: Vec<Vec<f64>>,
    trig: Vec<(f64, f64)>,
    /// The cell layers at each wall that are not fluid (1 under the on-node bounce, whose
    /// boundary nodes are the wall, 0 otherwise), and the offset of the fluid rows' distance:
    /// row j lies at (j + half) Δx from its wall.
    wall_rows: usize,
    half: f64,
}

impl ChannelStats {
    /// The statistics of a channel 2δ wide whose walls are halfway outside the outermost cells
    /// (`wall_rows` 0: 2N cells across) or through their centres (the on-node bounce,
    /// `wall_rows` 1: 2N + 1 cells across, the outermost layers the boundary nodes).
    fn new(
        grid: fieldcraft_lbm::voxel::Grid,
        vs: f64,
        re_tau: f64,
        nu: f64,
        turbulence: Turbulence,
        sigma: Option<HrrSigma>,
        wall_rows: usize,
    ) -> ChannelStats {
        use std::f64::consts::PI;
        let [nx, ny, _] = grid.cells;
        let (dx, half) = if wall_rows > 0 { (2.0 / (ny - 1) as f64, 0.0) } else { (2.0 / ny as f64, 0.5) };
        // The fluid rows of the lower half: from the first off the wall to the centreline's.
        let lower = wall_rows..ny / 2 + wall_rows;
        // The row nearest y⁺ 48 (both, when two are as near) and the centreline row.
        let y_plus = |j: usize| (j as f64 + half) * dx * re_tau;
        let nearest = lower.clone().map(|j| (y_plus(j) - 48.0).abs()).fold(f64::INFINITY, f64::min);
        let mut spectrum_rows: Vec<usize> = lower.clone().filter(|&j| (y_plus(j) - 48.0).abs() <= nearest + 1e-9).collect();
        if !spectrum_rows.contains(&(lower.end - 1)) {
            spectrum_rows.push(lower.end - 1);
        }
        let trig = (0..nx).map(|q| (2.0 * PI * q as f64 / nx as f64).sin_cos()).collect();
        ChannelStats {
            grid,
            vs,
            re_tau,
            dx,
            nu,
            l_max: grid.cells.iter().copied().min().unwrap_or(1) as f64,
            turbulence,
            sigma,
            rows: vec![RowSums { sigma_min: f64::INFINITY, ..RowSums::default() }; ny],
            samples: 0,
            spectra: vec![vec![0.0; nx / 2 + 1]; spectrum_rows.len()],
            spectrum_rows,
            trig,
            wall_rows,
            half,
        }
    }

    /// Add one snapshot of `velocity` (lattice units, every cell in grid order).
    fn add(&mut self, velocity: &[[f64; 3]]) {
        use fieldcraft_lbm::{hrr, vreman};
        let grid = self.grid;
        let [nx, ny, nz] = grid.cells;
        let at = |i: usize, j: usize, k: usize| velocity.get(grid.index([i, j, k])).copied().unwrap_or([0.0; 3]);
        // Lattice gradients to wall units.
        let to_plus = self.vs / (self.dx * self.re_tau);
        // Per row: Vreman's c, or Smagorinsky's (C D)² with Van Driest's D at the row's nominal y⁺.
        let constant: Vec<f64> = (0..ny)
            .map(|j| match self.turbulence {
                Turbulence::Vreman { constant, .. } => vreman::model_constant(constant),
                Turbulence::Smagorinsky { constant, van_driest, .. } => {
                    let y_plus = (j as f64 + self.half).min((ny - 1 - j) as f64 + self.half) * self.dx * self.re_tau;
                    (constant * van_driest.map_or(1.0, |a| 1.0 - (-y_plus / a).exp())).powi(2)
                }
            })
            .collect();
        let (first, last) = (self.wall_rows, ny - 1 - self.wall_rows);
        for (cell, &u) in velocity.iter().enumerate() {
            let [i, j, k] = grid.coords(cell);
            if j < first || j > last {
                continue;
            }
            let nb: hrr::Neighbours = [
                [Some(at((i + nx - 1) % nx, j, k)), Some(at((i + 1) % nx, j, k))],
                [(j > first).then(|| at(i, j - 1, k)), (j < last).then(|| at(i, j + 1, k))],
                [Some(at(i, j, (k + nz - 1) % nz)), Some(at(i, j, (k + 1) % nz))],
            ];
            let g = hrr::gradient(u, &nb);
            let s = hrr::strain(&g);
            let ss = s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
            let c = constant.get(j).copied().unwrap_or(0.0);
            let nu_t = match self.turbulence {
                Turbulence::Vreman { .. } => vreman::eddy_viscosity(&g, c),
                Turbulence::Smagorinsky { .. } => c * (2.0 * ss).sqrt(),
            };
            let ratio = nu_t / self.nu;
            let sigma = match self.sigma {
                Some(HrrSigma::Dynamic) => hrr::dynamic_sigma(&g, &hrr::laplacian(u, &nb), nu_t, self.nu, self.l_max),
                Some(HrrSigma::Value(s)) => s,
                None => 1.0,
            };
            let ss_plus = ss * to_plus * to_plus;
            let Some(row) = self.rows.get_mut(j) else { continue };
            for ((m, mm), x) in row.u.iter_mut().zip(&mut row.uu).zip(u) {
                *m += x;
                *mm += x * x;
            }
            row.uv += u[0] * u[1];
            row.nu_t += ratio;
            row.eps_molecular += 2.0 * ss_plus;
            row.eps_subgrid += 2.0 * ratio * ss_plus;
            row.viscous += g[1][0] * to_plus;
            row.subgrid += ratio * (g[1][0] + g[0][1]) * to_plus;
            row.sigma += sigma;
            row.sigma_min = row.sigma_min.min(sigma);
        }
        // The spectra: a direct transform of each line (n_x is not a power of two).
        let mut line = vec![0.0; nx];
        for (spectrum, &j) in self.spectra.iter_mut().zip(&self.spectrum_rows) {
            for row in [j, ny - 1 - j] {
                let plane = (0..nz).flat_map(|k| (0..nx).map(move |i| (i, k))).map(|(i, k)| at(i, row, k)[0]).sum::<f64>() / (nx * nz) as f64;
                for k in 0..nz {
                    for (i, v) in line.iter_mut().enumerate() {
                        *v = (at(i, row, k)[0] - plane) * self.vs;
                    }
                    for (m, e) in spectrum.iter_mut().enumerate() {
                        let (mut re, mut im) = (0.0, 0.0);
                        for (i, v) in line.iter().enumerate() {
                            let (s, c) = self.trig[(m * i) % nx];
                            re += v * c;
                            im -= v * s;
                        }
                        let both = if m == 0 || 2 * m == nx { 1.0 } else { 2.0 };
                        *e += both * (re * re + im * im) / (nx * nx) as f64;
                    }
                }
            }
        }
        self.samples += 1;
    }

    /// The folded rows and the spectra.
    fn finish(&self) -> (Vec<HanRow>, Vec<HanSpectrum>) {
        let [nx, ny, nz] = self.grid.cells;
        let per = (self.samples * nx * nz).max(1) as f64;
        let vs = self.vs;
        let rows = (self.wall_rows..ny / 2 + self.wall_rows)
            .map(|j| {
                let (b, t) = (self.rows[j], self.rows[ny - 1 - j]);
                let mean = |r: &RowSums| r.u.map(|x| x / per);
                let (mb, mt) = (mean(&b), mean(&t));
                let rms = |r: &RowSums, m: [f64; 3], a: usize| (r.uu[a] / per - m[a] * m[a]).max(0.0).sqrt() * vs;
                let uv = |r: &RowSums, m: [f64; 3]| (r.uv / per - m[0] * m[1]) * vs * vs;
                HanRow {
                    y: (j as f64 + self.half) * self.dx,
                    u: 0.5 * (mb[0] + mt[0]) * vs,
                    rms: [0, 1, 2].map(|a| 0.5 * (rms(&b, mb, a) + rms(&t, mt, a))),
                    // The upper half's wall-normal velocity points the other way.
                    resolved: 0.5 * (-uv(&b, mb) + uv(&t, mt)),
                    viscous: 0.5 * (b.viscous - t.viscous) / per,
                    subgrid: 0.5 * (b.subgrid - t.subgrid) / per,
                    nu_t: 0.5 * (b.nu_t + t.nu_t) / per,
                    eps_molecular: 0.5 * (b.eps_molecular + t.eps_molecular) / per,
                    eps_subgrid: 0.5 * (b.eps_subgrid + t.eps_subgrid) / per,
                    sigma: (0.5 * (b.sigma + t.sigma) / per, b.sigma_min.min(t.sigma_min)),
                }
            })
            .collect();
        let lines = (self.samples * 2 * nz).max(1) as f64;
        let spectra = self
            .spectrum_rows
            .iter()
            .zip(&self.spectra)
            .map(|(&row, e)| HanSpectrum {
                row,
                y_plus: (row as f64 + self.half) * self.dx * self.re_tau,
                energy: e.iter().map(|x| x / lines).collect(),
            })
            .collect();
        (rows, spectra)
    }
}

/// The hybrid collision's own viscosity beyond the one it is run with, at wavenumber `k` (radians
/// per cell): `(r(k)/k² − (τ̄ − ½)/3) / ν`, with r the exact decay rate per step of a transverse
/// shear wave under the collision at relaxation time `tau_bar` and weight `sigma`
/// ([`common::exact_rate`], the 2 × 2 map of crates/lbm/tests/common/mod.rs, checked against the
/// measured decay in crates/lbm/tests/hrr.rs) and ν the molecular lattice viscosity: the scheme's
/// term beside the subgrid ⟨ν_t⟩/ν. σ = 1 gives the recursive regularised collision's own.
fn scheme_viscosity(k: f64, tau_bar: f64, sigma: f64, nu: f64) -> f64 {
    (exact_rate(k, tau_bar, sigma) / (k * k) - (tau_bar - 0.5) / 3.0) / nu
}

/// The wavenumbers k_x Δx of the band the MRT channel piled energy into (docs/verification.md,
/// "Wall-function bounce, method two"), where the scheme's term is reported.
const HAN_BAND: [f64; 3] = [1.5, 1.8, 2.1];

/// The scheme's term at a spectrum's row ([`scheme_viscosity`]) for a hybrid-collision run: τ̄ is
/// the row's ½ + 3 (ν + ⟨ν_t⟩) where the subgrid viscosity enters τ, else the molecular τ; σ the
/// fixed one or the row's mean dynamic σ. Per wavenumber of [`HAN_BAND`] then the premultiplied
/// peak's: (k Δx, ν_s/ν at the run's σ, ν_s/ν at σ = 1). `None` for the MRT (the map is the
/// regularised family's).
fn han_scheme(setup: HanSetup, row: &HanRow, peak_k: f64, nu: f64) -> Option<(f64, f64, Vec<[f64; 3]>)> {
    let sigma = match setup.sigma()? {
        HrrSigma::Value(s) => s,
        HrrSigma::Dynamic => row.sigma.0,
    };
    let applied = if setup.subgrid_applied() { row.nu_t } else { 0.0 };
    let tau_bar = 0.5 + 3.0 * nu * (1.0 + applied);
    let at = HAN_BAND
        .iter()
        .copied()
        .chain([peak_k])
        .map(|k| [k, scheme_viscosity(k, tau_bar, sigma, nu), scheme_viscosity(k, tau_bar, 1.0, nu)])
        .collect();
    Some((tau_bar, sigma, at))
}

/// ε_u⁺ (Han et al.'s eq. 19, √(Σ (u⁺ − u⁺_S)² / Σ u⁺_S²) against Spalding's law u⁺_S) over the
/// half-channel's cell centres (their u⁺ "at y⁺ = i", read as the grid points: interpolating
/// between them would add an error of its own, 2 % in the buffer layer on Spalding's own profile
/// at 40 cells), in the buffer layer (y⁺ < 100), the log layer (y⁺ ≥ 100) and the whole half
/// channel; and b, the mean of U⁺ − ln(y⁺) / 0.41 over the cells with 100 ≤ y⁺ ≤ 0.3 Re_τ (the
/// paper states neither the points of its sums nor the range of its fit: these are our readings).
fn han_measures(profile: &[(f64, f64)], re_tau: f64) -> ([f64; 3], f64) {
    let mut sums = [[0.0f64; 2]; 3];
    for &(y, u) in profile {
        let yp = y * re_tau;
        let s = spalding_u_plus(yp);
        for (k, inside) in [yp < 100.0, yp >= 100.0, true].into_iter().enumerate() {
            if inside {
                sums[k][0] += (u - s) * (u - s);
                sums[k][1] += s * s;
            }
        }
    }
    let eps = sums.map(|[d, s]| (d / s).sqrt());
    let log: Vec<f64> =
        profile.iter().filter(|(y, _)| (100.0..=0.3 * re_tau).contains(&(y * re_tau))).map(|(y, u)| u - (y * re_tau).ln() / 0.41).collect();
    (eps, log.iter().sum::<f64>() / log.len() as f64)
}

/// Whether the study's wall-function bounce is Han et al.'s on-node eq. 14 (`spalding_bounce`,
/// CPU only): `FIELDCRAFT_HAN_WALL=onnode`; otherwise our halfway variant
/// (`spalding_bounce_halfway`, CPU or GPU).
fn han_on_node() -> bool {
    std::env::var("FIELDCRAFT_HAN_WALL").is_ok_and(|v| v == "onnode")
}

/// The subgrid model of a Han et al. run: `setup`'s, with the wall-function bounce of
/// [`han_on_node`] when `wfb`.
fn han_turbulence(setup: HanSetup, wfb: bool) -> Turbulence {
    match setup.turbulence(wfb) {
        Turbulence::Smagorinsky { constant, prandtl, van_driest, wall_model: Some(_) } if han_on_node() => {
            Turbulence::Smagorinsky { constant, prandtl, van_driest, wall_model: Some(fieldcraft_model::WallModel::SpaldingBounce) }
        }
        Turbulence::Vreman { constant, prandtl, wall_model: Some(_) } if han_on_node() => {
            Turbulence::Vreman { constant, prandtl, wall_model: Some(fieldcraft_model::WallModel::SpaldingBounce) }
        }
        t => t,
    }
}

/// A Han et al. channel's flow, started (see [`han_channel`]).
struct HanFlow {
    setup: Setup,
    flow: Flow,
    /// The Mach number run, the time step (δ/u_τ) and the box's streamwise length (δ).
    mach: f64,
    dt: f64,
    lx: f64,
}

/// [`han_channel`]'s box, set-up and perturbed start for `case` on `collision` (σ `sigma`) with
/// `turbulence_of`, run to `end` δ/u_τ; with `on_node` the box is 2δ + Δx high, the walls through
/// the outermost cells' centres (Han et al.'s on-node bounce).
fn han_flow(
    case: &HanCase,
    collision: Collision,
    sigma: Option<HrrSigma>,
    turbulence_of: Turbulence,
    on_node: bool,
    end: f64,
) -> Result<HanFlow, String> {
    use std::f64::consts::PI;
    let n = case.n;
    let dx = 1.0 / n as f64;
    let wall = if on_node { 0.5 * dx } else { 0.0 };
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0 + 2.0 * wall, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let uc = case.re_centre / case.re_tau;
    // τ − ½ = 3 ν Δt / Δx² grows with the lattice speed: the Mach number that puts τ on the floor.
    let floor = fieldcraft_lbm::setup::MIN_TAU_TURBULENT;
    let tau_at = |mach: f64| 0.5 + 3.0 / case.re_tau * (mach / 3f64.sqrt() / uc) / dx;
    let mach = if case.raise_to_floor && tau_at(case.mach) < floor {
        case.mach * (floor - 0.5) / (tau_at(case.mach) - 0.5) * (1.0 + 1e-9)
    } else {
        case.mach
    };
    let dt = mach / 3f64.sqrt() * dx / uc;
    let p = les(sigma, Some(turbulence_of))(collide(end, dt, None, collision));
    let m = model(1.0, 1.0 / case.re_tau, d, p);
    let limits = fieldcraft_lbm::Limits { max_cell_updates: 10_000_000_000_000, ..fieldcraft_lbm::Limits::default() };
    let setup = Setup::new(&m, None, &[], &limits).map_err(|e| format!("refused: {e}"))?;
    let mut flow = Flow::new(&setup, &[]).map_err(|e| format!("refused: {e}"))?;
    let grid = *flow.grid();
    let vs = setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    let ub0 = 0.87 * uc;
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0 - wall;
        let wall = (1.0 - eta * eta).max(0.0);
        let u = 1.15 * ub0 * (1.0 - eta.abs().powi(8)) + 2.0 * (2.0 * kz * z).sin() * (kx * x).cos() * wall;
        let v = 1.5 * wall * (2.0 * kz * z).cos() * (3.0 * kx * x).sin();
        let w = 1.5 * wall * eta * (kx * x).sin() * (kz * z).cos();
        (1.0, [u / vs, v / vs, w / vs])
    });
    Ok(HanFlow { setup, flow, mach, dt, lx })
}

/// Han et al.'s channel (their §3): 2πD × 2D × πD on `case.n` cells per half-height D, periodic
/// in x and z, driven by the constant pressure gradient of the nominal u_τ = 1 (δ = 1; the paper
/// normalises by the nominal u_τ), with the wall-function bounce (`wfb`) or plain bounce-back. The
/// collision and subgrid model are `setup`'s ([`HanSetup`]); the paper's is D3Q19 MRT with
/// Smagorinsky's model at C_s 0.1 and Van Driest's damping (A⁺ 26: the paper names the damping,
/// not its constant). The lattice speed is the paper's Mach number at the centreline speed its Re
/// implies (U_c⁺ = Re / Re_τ). `times.0` δ/u_τ (20) discarded, then U_b's block means of 2 δ/u_τ
/// until their 95 % interval is under 1 % with a lag-1 correlation under 0.5 and at least 50 δ/u_τ
/// averaged, at most `times.1` (150). Every 0.05 δ/u_τ in the window the velocity field is read
/// back for the profile and [`ChannelStats`]. A case whose Mach number puts τ under the floor is
/// refused, except B40 (Re_τ 2003, 40 cells: τ 0.500099 at their Ma 0.069), which runs at the
/// smallest Mach number that reaches the floor, 0.0693, a deviation printed and recorded with
/// both Mach numbers.
fn han_channel(case: &HanCase, setup_of: HanSetup, wfb: bool, cpu: bool, times: (f64, f64)) -> Option<Result<HanResult, String>> {
    let ((start, end), min_average, block, tolerance) = (times, 50.0, 2.0, 1e-2);
    // Under the on-node bounce the walls pass through the outermost cells' centres: the box is
    // 2δ + Δx high (2N + 1 cells) for 2δ between the walls; it runs on the CPU only.
    let on_node = wfb && han_on_node();
    let cpu = cpu || on_node;
    let turbulence_of = han_turbulence(setup_of, wfb);
    let HanFlow { setup, mut flow, mach, dt, lx } = match han_flow(case, setup_of.collision(), setup_of.sigma(), turbulence_of, on_node, end) {
        Ok(h) => h,
        Err(e) => return Some(Err(e)),
    };
    let kind = if wfb { "WFB" } else { "BB" };
    let label = setup_of.name();
    println!(
        "  Han {} {kind} [{label}]: τ {:.6}, Mach {mach:.5} at the centreline (theirs {}{})",
        case.name,
        setup.tau_plus,
        case.mach,
        if mach != case.mach { "; raised to the τ floor, a deviation" } else { "" }
    );
    let grid = *flow.grid();
    let vs = setup.velocity_scale();
    let mut device = if cpu { None } else { Some(gpu(&flow, Streaming::InPlace)?) };
    let control = ((0.05 / dt) as u64).max(1);
    let from = (start / dt) as u64;
    let mut bulk = Blocks::new(((block / 0.05) as usize).max(1));
    let nu = (setup.tau_plus - 0.5) / 3.0;
    let mut stats = ChannelStats::new(grid, vs, case.re_tau, nu, turbulence_of, setup_of.sigma(), usize::from(on_node));
    let (mut force_sum, mut body_sum) = (0.0, 0.0);
    let clock = std::time::Instant::now();
    let mut k = 0u64;
    let mut stop = "the step's end".to_string();
    let mut window = (f64::NAN, f64::NAN);
    while k < setup.steps {
        let step = control.min(setup.steps - k);
        if let Err(e) = backend(&mut device, &mut flow).advance(step) {
            return Some(Err(format!("unstable before {:.2} δ/u_τ: {e}", (k + step) as f64 * dt)));
        }
        k += step;
        let b = backend(&mut device, &mut flow);
        let mean = b.mean_velocity(0).unwrap() * vs;
        if k >= from {
            if window.0.is_nan() {
                window.0 = k as f64 * dt;
            }
            bulk.push(mean);
            // The walls' force against the body force on the fluid, this step.
            let forces = b.forces().unwrap();
            force_sum += forces.iter().map(|f| f[0]).sum::<f64>();
            body_sum += setup.acceleration[0] * b.mass().unwrap();
            let fields = (0..3).find_map(|_| {
                let (a, c) = (b.fields().ok()?, b.fields().ok()?);
                (a.velocity == c.velocity && a.velocity.iter().all(|u| u.iter().all(|v| v.abs() < 1.0))).then_some(a)
            });
            if let Some(fields) = fields {
                stats.add(&fields.velocity);
            }
            if k as f64 * dt >= start + min_average && bulk.half_width() < tolerance * bulk.mean() && bulk.lag1() < 0.5 {
                stop = format!("steady at {:.1} δ/u_τ", k as f64 * dt);
                break;
            }
        }
        if k % (setup.steps / 20).max(1) < control {
            println!("  Han {} {kind} [{label}]: t = {:.1}, U_b⁺ {mean:.3}, {:.0} s", case.name, k as f64 * dt, clock.elapsed().as_secs_f64());
        }
    }
    window.1 = k as f64 * dt;
    let (rows, spectra) = stats.finish();
    let profile: Vec<(f64, f64)> = rows.iter().map(|r| (r.y, r.u)).collect();
    let rms = rows.iter().map(|r| r.rms[0]).collect();
    let (eps, b) = han_measures(&profile, case.re_tau);
    // Lattice momentum per step against the body force per step: the same scale.
    let carried = force_sum / body_sum;
    Some(Ok(HanResult {
        profile,
        rms,
        bulk: bulk.mean(),
        eps,
        b,
        carried,
        window,
        stop,
        steps: k,
        seconds: clock.elapsed().as_secs_f64(),
        mach,
        tau: setup.tau_plus,
        rows,
        spectra,
        snapshots: stats.samples,
        nx: grid.cells[0],
        lx,
    }))
}

/// Method two: Han, Ooka and Kikumoto's channel (see [`han_channel`]), each case with the
/// wall-function bounce (WFB) and its bounce-back control (BB), on each configuration of
/// `FIELDCRAFT_HAN_SETUP` ([`HanSetup`], comma-separated; `paper` when unset): `paper` (MRT,
/// Smagorinsky 0.1, Van Driest 26), `hrr-vreman` (the hybrid collision at σ 0.99, Vreman 0.18:
/// Jacob et al. 2018's case 3), `hrr-dynamic` (dynamic σ from Vreman 0.18: their case 4) and
/// `hrr-smagorinsky` (σ 0.99 with the paper's subgrid model), or `all`. The gate is Han et al.'s:
/// ε_u⁺ over the whole channel and b within 15 % of their Tables 2 and 3; on `paper` for both
/// walls (as recorded), on the hybrid set-ups for the WFB runs, the BB runs reported. C_f =
/// 2 / U_b⁺² is reported against Dean's correlation and, at Re_τ 2003, Hoyas and Jiménez.
///
/// Every run prints, and records, per row of the folded half channel U⁺, u′, v′, w′ rms⁺, the
/// shear stress budget (−⟨u′v′⟩⁺, viscous, subgrid against the total 1 − y/δ), ⟨ν_t⟩/ν, the
/// molecular and subgrid dissipation and, on the hybrid collision, σ's row mean and minimum
/// ([`ChannelStats`]: which are exact and which estimates); u′'s streamwise spectrum at the rows
/// nearest y⁺ 48 and the centreline, averaged over z and the window's snapshots; and on the hybrid
/// collision the scheme's own viscosity at k_x Δx 1.5, 1.8, 2.1 and the spectrum's premultiplied
/// peak ([`han_scheme`]), beside the subgrid ⟨ν_t⟩/ν, so a pass is not credited to the subgrid
/// model when the collision's filter did the work. Each run writes a JSON record
/// `han-<setup>-<case>-<wall>.json` into `FIELDCRAFT_STUDY_RECORDS` (a directory) when set. Set
/// `FIELDCRAFT_HAN` to a case (`A20`, `A40`, `B40`, `B80`), several separated by commas, or `1`
/// for all; `_BB=0` skips the controls; `FIELDCRAFT_HAN_CPU=1` runs on the CPU in double
/// precision; `FIELDCRAFT_HAN_WALL=onnode` runs the WFB cases with Han et al.'s on-node eq. 14
/// (`spalding_bounce`, on the CPU, the box 2δ + Δx high so the walls through the outermost
/// cells' centres stay 2δ apart; [`han_on_node`]) instead of our halfway variant;
/// `FIELDCRAFT_HAN_START` and `FIELDCRAFT_HAN_END` (δ/u_τ, 20 and 150) shorten the
/// discarded start and the run for a smoke test (a run shorter than the protocol is not a result).
#[test]
#[ignore = "hours, GPU: set FIELDCRAFT_HAN=A20,A40,B40 (docs/verification.md, \"Long studies\")"]
fn study_han_ooka_and_kikumoto_s_wall_function_bounce() {
    let Ok(which) = std::env::var("FIELDCRAFT_HAN") else {
        println!("Han et al.'s channel runs only with FIELDCRAFT_HAN set (hours, GPU): skipped");
        return;
    };
    let setups: Vec<HanSetup> = match std::env::var("FIELDCRAFT_HAN_SETUP").as_deref() {
        Err(_) => vec![HanSetup::Paper],
        Ok("all") => HAN_SETUPS.to_vec(),
        Ok(list) => list
            .split(',')
            .map(|s| {
                HAN_SETUPS
                    .into_iter()
                    .find(|h| h.name() == s.trim())
                    .unwrap_or_else(|| panic!("FIELDCRAFT_HAN_SETUP: no set-up `{s}` ({} or all)", HAN_SETUPS.map(HanSetup::name).join(", ")))
            })
            .collect(),
    };
    let time = |var: &str, default: f64| {
        std::env::var(var).ok().and_then(|v| v.parse::<f64>().ok()).filter(|t| t.is_finite() && *t >= 0.0).unwrap_or(default)
    };
    let times = (time("FIELDCRAFT_HAN_START", 20.0), time("FIELDCRAFT_HAN_END", 150.0));
    let cpu = std::env::var_os("FIELDCRAFT_HAN_CPU").is_some();
    let controls = std::env::var("FIELDCRAFT_HAN_BB").map_or(true, |v| v != "0");
    let records = std::env::var_os("FIELDCRAFT_STUDY_RECORDS").map(std::path::PathBuf::from);
    let commit = std::process::Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string());
    let mut failed = Vec::new();
    for case in HAN_CASES.iter().filter(|c| which == "1" || which.split(',').any(|w| w == c.name)) {
        for &setup in &setups {
            let label = setup.name();
            for wfb in [true, false].into_iter().filter(|w| *w || controls) {
                let kind = if wfb { "WFB" } else { "BB" };
                let (paper_eps, paper_b) = if wfb { case.wfb } else { case.bb };
                let Some(result) = han_channel(case, setup, wfb, cpu, times) else { return };
                let r = match result {
                    Ok(r) => r,
                    Err(e) => {
                        println!("Han {} {kind} [{label}]: {e}", case.name);
                        failed.push(format!("{} {kind} [{label}]: {e}", case.name));
                        continue;
                    }
                };
                han_report(case, setup, wfb, cpu || (wfb && han_on_node()), &r, records.as_deref(), commit.as_deref());
                // The paper's set-up is gated on both walls (as recorded); the hybrid ones on WFB.
                let gated = setup.gated(wfb);
                let eps_ok = (r.eps[2] / paper_eps[2] - 1.0).abs() <= 0.15;
                let b_ok = (r.b / paper_b - 1.0).abs() <= 0.15;
                if gated && !(eps_ok && b_ok) {
                    failed.push(format!(
                        "{} {kind} [{label}]: ε {:.3} (theirs {}), b {:.2} (theirs {})",
                        case.name, r.eps[2], paper_eps[2], r.b, paper_b
                    ));
                }
            }
        }
    }
    assert!(failed.is_empty(), "outside Han et al.'s results by more than 15 %: {failed:?}");
}

/// Print one Han et al. run and write its record (see [`study_han_ooka_and_kikumoto_s_wall_function_bounce`]).
fn han_report(case: &HanCase, setup: HanSetup, wfb: bool, cpu: bool, r: &HanResult, records: Option<&std::path::Path>, commit: Option<&str>) {
    let kind = if wfb { "WFB" } else { "BB" };
    let label = setup.name();
    let (paper_eps, paper_b) = if wfb { case.wfb } else { case.bb };
    let cf = 2.0 / (r.bulk * r.bulk);
    let cf_dean = 0.073 * (2.0 * r.bulk * case.re_tau).powf(-0.25);
    let cf_dns = (case.re_tau == 2003.0).then(|| 2.0 / (HJ_2000_BULK * HJ_2000_BULK));
    let nu = (r.tau - 0.5) / 3.0;
    let hybrid = setup.sigma().is_some();
    println!(
        "Han {} {kind} [{label}] ({}, Re_τ {}, {}/δ): ε_u⁺ buffer {:.3} (theirs {}), log {:.3} ({}), whole {:.3} ({}, {:+.0} %); b {:.2} ({}, {:+.0} %); U_b⁺ {:.3}, C_f {cf:.4e} ({:+.1} % on Dean{}); walls carry {:.4} of the body force; Mach {:.4} (theirs {}), τ {:.6}; averaged {:.0}–{:.0} δ/u_τ ({} snapshots), {}; {} steps, {:.0} s",
        case.name,
        if cpu { "CPU, f64" } else { "GPU, f32" },
        case.re_tau,
        case.n,
        r.eps[0],
        paper_eps[0],
        r.eps[1],
        paper_eps[1],
        r.eps[2],
        paper_eps[2],
        100.0 * (r.eps[2] / paper_eps[2] - 1.0),
        r.b,
        paper_b,
        100.0 * (r.b / paper_b - 1.0),
        r.bulk,
        100.0 * (cf / cf_dean - 1.0),
        cf_dns.map_or(String::new(), |d| format!(", {:+.1} % on Hoyas & Jiménez", 100.0 * (cf / d - 1.0))),
        r.carried,
        r.mach,
        case.mach,
        r.tau,
        r.window.0,
        r.window.1,
        r.snapshots,
        r.stop,
        r.steps,
        r.seconds
    );
    let (nu_t_source, estimate) = setup.nu_t_source();
    let applied = if setup.subgrid_applied() { "" } else { " (ν_t sets σ only, not in τ)" };
    println!("  Han {} {kind} [{label}] per row{estimate}{applied}: total 1 − y/δ = −⟨u′v′⟩⁺ + viscous + subgrid + residual", case.name);
    for row in &r.rows {
        let yp = row.y * case.re_tau;
        let (subgrid, residual) = row.closure(setup.subgrid_applied());
        let total = 1.0 - row.y;
        let sigma = if hybrid { format!(", σ mean {:.4} min {:.4}", row.sigma.0, row.sigma.1) } else { String::new() };
        println!(
            "    y/δ {:.4} y⁺ {yp:7.1}: U⁺ {:7.3} (Spalding {:7.3}), rms⁺ u′ {:.3} v′ {:.3} w′ {:.3}; {total:.3} = {:.3} + {:.3} + {subgrid:.3} + {:.3}; ν_t/ν {:.3}; ε⁺ molecular {:.3e} subgrid {:.3e}{sigma}",
            row.y,
            row.u,
            spalding_u_plus(yp),
            row.rms[0],
            row.rms[1],
            row.rms[2],
            row.resolved,
            row.viscous,
            residual,
            row.nu_t,
            row.eps_molecular,
            row.eps_subgrid,
        );
    }
    let k_delta = |m: usize| 2.0 * std::f64::consts::PI * m as f64 / r.lx;
    let k_dx = |m: usize| 2.0 * std::f64::consts::PI * m as f64 / r.nx as f64;
    let mut scheme_records = Vec::new();
    for s in &r.spectra {
        let peak = s.peak();
        println!(
            "  Han {} {kind} [{label}] E_uu⁺(k_x) at y⁺ {:.0} (row {}): premultiplied peak at k_x δ {:.1} (k_x Δx {:.2}), {:.0} % of ⟨u′²⟩ above k_x Δx π/2; k_x δ = {:.4} m: [{}]",
            case.name,
            s.y_plus,
            s.row,
            k_delta(peak),
            k_dx(peak),
            100.0 * s.upper_share(r.nx),
            k_delta(1),
            s.energy.iter().map(|e| format!("{e:.3e}")).collect::<Vec<_>>().join(", ")
        );
        let Some(row) = r.rows.get(s.row) else { continue };
        if let Some((tau_bar, sigma, at)) = han_scheme(setup, row, k_dx(peak), nu) {
            let list: Vec<String> = at.iter().map(|[k, v, v1]| format!("k_x Δx {k:.2}: {v:.3} (σ 1: {v1:.3})")).collect();
            println!("    scheme's ν_s/ν at τ̄ {tau_bar:.6}, σ {sigma:.4}: {}; subgrid ⟨ν_t⟩/ν {:.3}{applied}", list.join(", "), row.nu_t);
            scheme_records.push(serde_json::json!({
                "row": s.row, "y_plus": s.y_plus, "tau_bar": tau_bar, "sigma": sigma, "nu_t_over_nu": row.nu_t, "subgrid_applied": setup.subgrid_applied(),
                "at": at.iter().map(|[k, v, v1]| serde_json::json!({"k_x_dx": k, "nu_scheme_over_nu": v, "nu_scheme_over_nu_sigma_one": v1})).collect::<Vec<_>>(),
            }));
        }
    }
    let Some(dir) = records else { return };
    let rows: Vec<_> = r
        .rows
        .iter()
        .map(|row| {
            let (_, residual) = row.closure(setup.subgrid_applied());
            serde_json::json!({
                "y": row.y, "y_plus": row.y * case.re_tau, "u_plus": row.u, "spalding": spalding_u_plus(row.y * case.re_tau),
                "u_rms_plus": row.rms[0], "v_rms_plus": row.rms[1], "w_rms_plus": row.rms[2],
                "resolved_shear": row.resolved, "viscous_shear": row.viscous, "subgrid_shear": row.subgrid, "total_shear": 1.0 - row.y,
                "residual": residual,
                "nu_t_over_nu": row.nu_t, "eps_molecular_plus": row.eps_molecular, "eps_subgrid_plus": row.eps_subgrid,
                "sigma_mean": hybrid.then_some(row.sigma.0), "sigma_min": hybrid.then_some(row.sigma.1),
            })
        })
        .collect();
    let spectra: Vec<_> = r
        .spectra
        .iter()
        .map(|s| {
            serde_json::json!({
                "row": s.row, "y_plus": s.y_plus, "k_x_delta": (0..s.energy.len()).map(k_delta).collect::<Vec<_>>(),
                "k_x_dx": (0..s.energy.len()).map(k_dx).collect::<Vec<_>>(), "e_uu_plus": s.energy,
                "premultiplied_peak_k_x_delta": k_delta(s.peak()), "upper_half_share": s.upper_share(r.nx),
            })
        })
        .collect();
    let turbulence = han_turbulence(setup, wfb);
    let record = serde_json::json!({
        "study": "han_ooka_kikumoto_2021", "setup": label, "case": case.name, "wall": kind,
        "wall_form": if !wfb { "bounce-back" } else if han_on_node() { "on-node eq. 14 (spalding_bounce)" } else { "halfway (spalding_bounce_halfway)" }, "re_tau": case.re_tau, "cells_per_delta": case.n,
        "collision": format!("{:?}", setup.collision()), "sigma": setup.sigma().map(|s| format!("{s:?}")),
        "turbulence": serde_json::to_value(turbulence).ok(), "subgrid_applied": setup.subgrid_applied(), "nu_t_source": nu_t_source,
        "sigma_source": (!setup.subgrid_applied()).then_some("fieldcraft_lbm::hrr::dynamic_sigma on the readback velocity (no per-cell readback on the device)"),
        "backend": if cpu { "cpu-f64" } else { "gpu-f32" }, "commit": commit,
        "profile": r.profile.iter().zip(&r.rms).map(|((y, u), rms)| serde_json::json!({"y": y, "y_plus": y * case.re_tau, "u_plus": u, "u_rms_plus": rms})).collect::<Vec<_>>(),
        "rows": rows, "spectra": spectra, "scheme": scheme_records, "snapshots": r.snapshots,
        "u_b_plus": r.bulk, "cf_carried": cf, "cf_dean": cf_dean, "cf_dns": cf_dns, "carried_over_body_force": r.carried,
        "eps": {"buffer": r.eps[0], "log": r.eps[1], "whole": r.eps[2]}, "eps_paper": paper_eps, "b": r.b, "b_paper": paper_b,
        "gated": setup.gated(wfb),
        "mach": r.mach, "mach_paper": case.mach, "tau": r.tau, "window": [r.window.0, r.window.1], "stop": r.stop, "steps": r.steps, "seconds": r.seconds,
    });
    let _ = std::fs::create_dir_all(dir);
    let _ = std::fs::write(dir.join(format!("han-{label}-{}-{kind}.json", case.name)), serde_json::to_string_pretty(&record).unwrap_or_default());
}

/// The stepper of a study: the GPU flow when there is one, else the CPU flow.
fn backend<'a>(device: &'a mut Option<Gpu>, flow: &'a mut Flow) -> &'a mut dyn Backend {
    match device.as_mut() {
        Some(g) => &mut **g,
        None => flow,
    }
}

/// The Han et al. channel's statistics ([`ChannelStats`]) on a field they can be worked out for:
/// u = U(y) + a sin(2π · 3 x / L_x) with U symmetric about the centreline, v = w = 0, two
/// snapshots. u′'s spectrum holds a²/2 at m = 3 and nothing else; the row means and rms are
/// the field's, the folded viscous stress the lower half's ∂U⁺/∂y⁺; −⟨u′v′⟩ vanishes; only u
/// varies, so the gradient has rank one and Vreman's ν_t vanishes, and dynamic σ is 1.
/// Smagorinsky's estimate is positive, and the scheme's term behaves as its map says.
#[test]
fn han_s_channel_statistics_read_a_known_field() {
    use std::f64::consts::PI;
    let (nx, ny, nz) = (12usize, 8usize, 3usize);
    let grid = fieldcraft_lbm::voxel::Grid { origin: [0.0; 3], spacing: 0.25, cells: [nx, ny, nz] };
    let (re_tau, vs, nu, a) = (100.0, 2.0, 1e-3, 0.01);
    let profile = |j: usize| 0.05 * (j as f64 + 0.5) * (ny as f64 - j as f64 - 0.5);
    let velocity: Vec<[f64; 3]> = (0..grid.len())
        .map(|c| {
            let [i, j, _] = grid.coords(c);
            [profile(j) + a * (2.0 * PI * 3.0 * i as f64 / nx as f64).sin(), 0.0, 0.0]
        })
        .collect();
    let vreman = Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: None };
    let mut stats = ChannelStats::new(grid, vs, re_tau, nu, vreman, Some(HrrSigma::Dynamic), 0);
    stats.add(&velocity);
    stats.add(&velocity);
    let (rows, spectra) = stats.finish();
    assert_eq!(rows.len(), ny / 2);
    // dx = 2 / ny = 0.25: y⁺ of rows 0 … 3 is 12.5, 37.5, 62.5, 87.5; 37.5 and 62.5 are both
    // 10.5 and 14.5 from 48, so row 1, and the centreline row 3.
    assert_eq!(spectra.iter().map(|s| s.row).collect::<Vec<_>>(), vec![1, 3]);
    for s in &spectra {
        assert_eq!(s.energy.len(), nx / 2 + 1);
        let want = 0.5 * (a * vs) * (a * vs);
        for (m, e) in s.energy.iter().enumerate() {
            let expected = if m == 3 { want } else { 0.0 };
            assert!((e - expected).abs() < 1e-12 * want, "m {m}: {e} vs {expected}");
        }
        assert_eq!(s.peak(), 3);
        assert!(s.upper_share(nx) < 1e-9, "m = 3 is k Δx π/2 exactly, not above it");
    }
    for (j, row) in rows.iter().enumerate() {
        assert!((row.u - profile(j) * vs).abs() < 1e-12, "{j}: {}", row.u);
        assert!((row.rms[0] - a * vs / 2f64.sqrt()).abs() < 1e-12 && row.rms[1] == 0.0 && row.rms[2] == 0.0);
        assert!(row.resolved.abs() < 1e-15);
        // Only u varies: α has one nonzero column, so Vreman's ν_t is 0 and dynamic σ 1.
        assert!(row.nu_t == 0.0 && row.eps_subgrid == 0.0 && row.subgrid == 0.0, "{j}: {}", row.nu_t);
        assert_eq!(row.sigma, (1.0, 1.0));
        // The viscous stress: the lower half's ∂U⁺/∂y⁺, the upper half's mirrored (centred
        // differences, one-sided at the wall row).
        let dudy = if j == 0 { profile(1) - profile(0) } else { 0.5 * (profile(j + 1) - profile(j - 1)) };
        let to_plus = vs / (0.25 * re_tau);
        assert!((row.viscous - dudy * to_plus).abs() < 1e-12, "{j}: {} vs {}", row.viscous, dudy * to_plus);
        assert!(row.eps_molecular > 0.0);
    }
    // Smagorinsky's estimate damps towards the wall and is positive in a shear.
    let smagorinsky = Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: Some(26.0), wall_model: None };
    let mut stats = ChannelStats::new(grid, vs, re_tau, nu, smagorinsky, None, 0);
    stats.add(&velocity);
    let (rows, _) = stats.finish();
    assert!(rows[0].nu_t > 0.0 && rows.iter().all(|r| r.sigma == (1.0, 1.0)));
    // The scheme's term: σ = 1 at small k is the regularised collision's own (near 0), σ < 1 adds.
    let tau = 0.50015;
    let nu = (tau - 0.5) / 3.0;
    assert!(scheme_viscosity(0.05, tau, 1.0, nu).abs() < 1e-3);
    assert!(scheme_viscosity(1.8, tau, 0.99, nu) > scheme_viscosity(1.8, tau, 1.0, nu));
}

/// Han et al.'s measures on Spalding's own profile: ε_u⁺ vanishes, and b is 5.44 at Re_τ 2003 on
/// 40 cells, Spalding's intercept 5.5 less its approach to the log law from below over 100 ≤ y⁺ ≤
/// 600.
#[test]
fn han_s_measures_read_spalding_s_profile_as_exact() {
    let re_tau = 2003.0;
    let profile: Vec<(f64, f64)> = (0..40).map(|j| ((j as f64 + 0.5) / 40.0, spalding_u_plus((j as f64 + 0.5) / 40.0 * re_tau))).collect();
    let (eps, b) = han_measures(&profile, re_tau);
    assert!(eps.iter().all(|e| *e < 1e-12), "{eps:?}");
    assert!((b - 5.44).abs() < 0.01, "{b}");
    // The bisection inverts Spalding's y⁺(u⁺).
    let u = spalding_u_plus(300.0);
    assert!((u - (300f64.ln() / 0.41 + 5.5)).abs() < 0.06, "{u}");
}

/// The collision's dissipation with Smagorinsky's model, without walls: the Taylor–Green vortex
/// (box 2πL, V₀ = L = 1, ν = 1 / Re), periodic, from the classical initial field with its pressure
/// (M. E. Brachet, D. I. Meiron, S. A. Orszag, B. G. Nickel, R. H. Morf, U. Frisch, J. Fluid
/// Mech. 130:411–452, 1983), stepped to `FIELDCRAFT_TGV_END` (20) with each variant in
/// `FIELDCRAFT_TGV` (`mrt:0.1,regularized:0.1,hrr@0.97:0.1,hrr@0.99:v0.18`: collision, for
/// `hrr` its σ after `@` (0.99 when left out, or `dynamic`), and the subgrid model's constant: Smagorinsky's C,
/// or with a `v` before it Vreman's model with that Smagorinsky-equivalent C_s), on
/// `FIELDCRAFT_TGV_N` cells per side (64) at lattice speed `FIELDCRAFT_TGV_SPEED` (0.05) and
/// Reynolds number `FIELDCRAFT_TGV_RE` (1600; 12000 puts τ below the Han et al. channel's).
/// Prints k(t) and −dk/dt (with dynamic σ on the CPU, its smallest and mean value in the step
/// printed); at Re 1600 the spectral DNS peaks at about 0.0127 near t ≈ 9.
/// `FIELDCRAFT_TGV_CPU=1` runs on the CPU.
#[test]
#[ignore = "minutes, GPU: set FIELDCRAFT_TGV=mrt:0.1 (docs/verification.md, \"Long studies\")"]
fn study_taylor_green_vortex_dissipation() {
    use std::f64::consts::PI;
    let Ok(which) = std::env::var("FIELDCRAFT_TGV") else {
        println!("the Taylor–Green study runs only with FIELDCRAFT_TGV set: skipped");
        return;
    };
    let n: usize = std::env::var("FIELDCRAFT_TGV_N").ok().and_then(|v| v.parse().ok()).unwrap_or(64);
    let speed: f64 = std::env::var("FIELDCRAFT_TGV_SPEED").ok().and_then(|v| v.parse().ok()).unwrap_or(0.05);
    let cpu = std::env::var_os("FIELDCRAFT_TGV_CPU").is_some();
    let re: f64 = std::env::var("FIELDCRAFT_TGV_RE").ok().and_then(|v| v.parse().ok()).unwrap_or(1600.0);
    let end: f64 = std::env::var("FIELDCRAFT_TGV_END").ok().and_then(|v| v.parse().ok()).unwrap_or(20.0);
    for variant in which.split(',') {
        let (name, c) = variant.split_once(':').unwrap_or((variant, "0.1"));
        // A `v` before the constant: Vreman's model with that Smagorinsky-equivalent C_s.
        let (vreman, c) = match c.strip_prefix('v') {
            Some(c) => (true, c),
            None => (false, c),
        };
        let c: f64 = c.parse().unwrap();
        let (name, sigma) = match name.split_once('@') {
            Some((n, "dynamic")) => (n, Some(fieldcraft_model::HrrSigma::Dynamic)),
            Some((n, s)) => (n, Some(fieldcraft_model::HrrSigma::Value(s.parse().unwrap()))),
            None => (name, None),
        };
        let collision = match name {
            "regularized" => Collision::Regularized,
            "trt" => Collision::Trt,
            "hrr" => Collision::Hrr,
            _ => Collision::Mrt,
        };
        let side = 2.0 * PI;
        let dx = side / n as f64;
        let dt = speed * dx;
        let mut d = domain([0.0; 3], [side; 3], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = PERIODIC;
        }
        let mut p = collide(end, dt, None, collision);
        if let Procedure::Fluid { turbulence, sigma: s, .. } = &mut p {
            *turbulence = Some(if vreman {
                fieldcraft_model::Turbulence::Vreman { constant: c, prandtl: 0.9, wall_model: None }
            } else {
                fieldcraft_model::Turbulence::Smagorinsky { constant: c, prandtl: 0.9, van_driest: None, wall_model: None }
            });
            *s = sigma;
        }
        let m = model(1.0, 1.0 / re, d, p);
        let limits = fieldcraft_lbm::Limits { max_cell_updates: 10_000_000_000_000, ..fieldcraft_lbm::Limits::default() };
        let setup = match Setup::new(&m, None, &[], &limits) {
            Ok(s) => s,
            Err(e) => {
                println!("TGV {variant}: refused: {e}");
                continue;
            }
        };
        let mut flow = Flow::new(&setup, &[]).unwrap();
        let grid = *flow.grid();
        let vs = setup.velocity_scale();
        flow.set_equilibrium(|cell| {
            let [x, y, z] = grid.centre(cell);
            let u = x.sin() * y.cos() * z.cos();
            let v = -x.cos() * y.sin() * z.cos();
            // p − p₀ = (ρ₀V₀²/16)(cos 2x + cos 2y)(cos 2z + 2), as a lattice density (c_s² = ⅓).
            let pressure = ((2.0 * x).cos() + (2.0 * y).cos()) * ((2.0 * z).cos() + 2.0) / 16.0;
            (1.0 + 3.0 * pressure / (vs * vs), [u / vs, v / vs, 0.0])
        });
        let mut device = if cpu { None } else { Some(gpu(&flow, Streaming::InPlace).unwrap()) };
        println!("TGV {variant} ({}, {n}³, lattice speed {speed}, Re {re}): τ {:.6}", if cpu { "CPU" } else { "GPU" }, setup.tau_plus);
        let every = ((0.25 / dt).round() as u64).max(1);
        let energy = |b: &mut dyn Backend| {
            let f = b.fields().unwrap();
            0.5 * f.velocity.iter().map(|u| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * vs * vs).sum::<f64>() / f.velocity.len() as f64
        };
        let mut k = 0u64;
        let mut last = (0.0, energy(backend(&mut device, &mut flow)));
        let mut peak = (0.0f64, 0.0f64, 0.0f64);
        while k < setup.steps {
            let step = every.min(setup.steps - k);
            if let Err(e) = backend(&mut device, &mut flow).advance(step) {
                println!("TGV {variant}: unstable at t {:.2}: {e}", (k + step) as f64 * dt);
                break;
            }
            k += step;
            let t = k as f64 * dt;
            let e = energy(backend(&mut device, &mut flow));
            let rate = -(e - last.1) / (t - last.0);
            if rate > peak.1 {
                peak = (0.5 * (t + last.0), rate, backend(&mut device, &mut flow).eddy_viscosity_ratio());
            }
            if k % (4 * every) < every {
                let sigma = flow.dynamic_sigma().filter(|_| cpu).map_or(String::new(), |(min, mean)| format!(", σ min {min:.4} mean {mean:.4}"));
                println!("    t {:5.2}: k {e:.5}, −dk/dt {rate:.5}{sigma}", 0.5 * (t + last.0));
            }
            last = (t, e);
        }
        println!(
            "TGV {variant} Re {re}: peak −dk/dt {:.5} at t {:.2} (spectral DNS ≈ 0.0127 near t ≈ 9), largest ν_t/ν there {:.3}",
            peak.1, peak.0, peak.2
        );
    }
}
