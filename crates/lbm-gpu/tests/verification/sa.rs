//! The Spalart–Allmaras model on the GPU against the CPU, with tolerances frozen here at their first
//! measurement: ν̃'s transport alone within 1e-5 of the largest CPU ν̃ after 200 steps, the
//! coupled RANS channel within 5e-5 for ν̃, 1e-4 / 1e-5 for velocity / density (the LES bounds),
//! 1e-3 for the largest ν_t/ν and the wall shear, SA-neg's smallest ν̃ within 1e-5 of ν; the wall
//! distance on the device is the CPU's rounded to f32, bit for bit; two runs agree bit for bit.
//! Both run the same fixed sub-step count (`Flow::set_sa_substeps`, `GpuFlow::fix_sa_substeps`),
//! so the sub-step splitting is not part of the difference. Each test prints what it measured.

use super::*;
use fieldcraft_model::{FaceCondition, WallModel};

/// The model with inflow ν̃/ν `ratio` and wall model `wall_model`.
fn sa(ratio: f64, wall_model: Option<WallModel>) -> Turbulence {
    Turbulence::SpalartAllmaras { prandtl: 0.9, wall_model, inflow_ratio: ratio }
}

/// `p` with turbulence model `t`.
fn with(mut p: Procedure, t: Turbulence) -> Procedure {
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(t);
    }
    p
}

/// The largest difference of two ν̃ fields over the CPU's largest |ν̃|, and that largest.
fn nu_tilde_difference(cpu: &mut Flow, gpu: &mut GpuFlow) -> (f64, f64) {
    let (a, b) = (cpu.nu_tilde().unwrap().unwrap(), gpu.nu_tilde().unwrap().unwrap());
    let scale = a.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
    (a.iter().zip(&b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max) / scale, scale)
}

/// A channel 24 × 16 × 8 cells: a uniform inflow on −x, a pressure outlet on +x, walls on y,
/// periodic in z (every face kind of the model's neighbour rule), lattice speed 0.05, ν 0.05 in
/// lattice units, inflow ν̃/ν `ratio`.
fn inflow_channel(collision: Collision, steps: f64, ratio: f64) -> Model {
    let dx = 1.0 / 16.0;
    let mut d = domain([0.0; 3], [1.5, 1.0, 0.5], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    model(1.0, 0.05 * dx, d, with(collide(steps * dx, dx, None, collision), sa(ratio, None)))
}

/// A forced channel two high at Re_τ about `re`, six cells per half-height, periodic in x and z
/// (the CPU tests' channel at Re_τ 180), `steps` time steps.
fn forced_channel(collision: Collision, steps: f64, re: f64) -> Model {
    let dx = 1.0 / 6.0;
    let mut d = domain([0.0; 3], [4.0 * dx, 2.0, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [15.7, 0.0, 0.0];
    let dt = 0.1 * dx / 15.7;
    model(1.0, 1.0 / re, d, with(collide(steps * dt, dt, None, collision), sa(3.0, None)))
}

/// ν̃'s transport alone (the model's sources off on both): a ν̃ front entering a channel from its
/// inflow (ν̃ 0 inside, ν at the inflow), carried by the flow, diffused, with the c_b2 term,
/// past walls, an outlet and periodic faces; within 1e-5 of the largest CPU ν̃ after 200 steps,
/// both layouts, the two-relaxation-time and multiple-relaxation-time collisions. The inflow is ν,
/// not the usual 3 ν: without the model's wall destruction, the c_b2 term at the corner where the
/// inflow meets a wall feeds on its own jump, and at 2 ν or more the CPU's transport alone runs
/// past 64 sub-steps (at 3 ν within 8 steps; measured, so not a device matter).
#[test]
fn sa_transport_follows_the_cpu() {
    for (collision, streaming) in
        [Collision::Trt, Collision::Mrt].into_iter().flat_map(|c| [Streaming::InPlace, Streaming::TwoBuffers].map(move |s| (c, s)))
    {
        let (setup, mut cpu) = prepare(&inflow_channel(collision, 200.0, 1.0), &[]);
        cpu.set_sa_model(false);
        cpu.set_sa_substeps(Some(2));
        cpu.set_sa_nu_tilde(&vec![0.0; cpu.grid().len()]).unwrap();
        let Some(mut gpu) = gpu(&cpu, streaming) else { return };
        gpu.set_sa_model(false).unwrap();
        gpu.fix_sa_substeps(Some(2)).unwrap();
        cpu.advance(setup.steps).unwrap();
        gpu.advance(setup.steps).unwrap();
        let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        println!("SA transport ({collision:?}, {streaming:?}), 200 steps: |Δν̃| {dv:.2e} of {scale:.3e}; velocity {du:.2e}, density {drho:.2e}");
        assert!(dv < 1e-5, "{dv}");
    }
}

/// SA-neg: a negative ν̃ patch in a periodic box at rest decays towards 0 on both backends and
/// never crosses it; the smallest ν̃ agrees within 1e-5 ν after 200 steps.
#[test]
fn sa_neg_recovers_like_the_cpu() {
    let n = 8usize;
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [1.0; 3], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let nu = 0.05 * dx;
    let (setup, mut cpu) = prepare(&model(1.0, nu, d, with(collide(200.0 * dx, dx, None, Collision::Trt), sa(3.0, None))), &[]);
    let nu_lattice = (setup.tau_plus - 0.5) / 3.0;
    let grid = *cpu.grid();
    let field: Vec<f64> = (0..grid.len())
        .map(|c| {
            let p = grid.coords(c);
            if p.iter().all(|&i| (2..6).contains(&i)) { -2.0 * nu_lattice } else { 0.0 }
        })
        .collect();
    cpu.set_sa_nu_tilde(&field).unwrap();
    cpu.set_sa_substeps(Some(2));
    let Some(mut gpu) = gpu(&cpu, Streaming::InPlace) else { return };
    gpu.fix_sa_substeps(Some(2)).unwrap();
    cpu.advance(setup.steps).unwrap();
    gpu.advance(setup.steps).unwrap();
    let low = |v: &[f64]| v.iter().fold(f64::INFINITY, |m, x| m.min(*x));
    let high = |v: &[f64]| v.iter().fold(f64::NEG_INFINITY, |m, x| m.max(*x));
    let (fa, fb) = (cpu.nu_tilde().unwrap().unwrap(), gpu.nu_tilde().unwrap().unwrap());
    let (a, b) = (low(&fa), low(&fb));
    assert!(high(&fa) <= 0.0 && high(&fb) <= 0.0, "never above 0: {} {}", high(&fa), high(&fb));
    println!("SA-neg after 200 steps: min ν̃/ν CPU {:.6e}, GPU {:.6e}, |Δ|/ν {:.2e}", a / nu_lattice, b / nu_lattice, (a - b).abs() / nu_lattice);
    assert!(a < 0.0 && a > -2.0 * nu_lattice && b <= 0.0 && b > -2.0 * nu_lattice, "decays towards 0 without crossing it");
    assert!((a - b).abs() < 1e-5 * nu_lattice, "{a} vs {b}");
}

/// The coupled RANS channel (ν_t in the collision) at Re_τ 20 (τ 0.5057): velocity, density,
/// ν̃, the largest ν_t/ν and each wall's shear (its mean over the last 200 steps, as a skin
/// friction is read) follow the CPU over 2000 steps, both layouts, the two-relaxation-time and
/// multiple-relaxation-time collisions, and the hybrid collision (σ 0.99) in place; and the hybrid
/// collision again at Re_τ 1000 (τ − ½ 1.15e-4, the gated cases' τ). A single step's shear is
/// not the measure: the two-relaxation-time collision's acoustic noise moves it 1e-3 between
/// backends whose fields agree to 5e-6 (the mean, 1e-5).
///
/// The two-relaxation-time collision is not run below τ − ½ ≈ 2e-3: at Λ 3/16 its antisymmetric
/// modes relax at ω⁻ = (τ⁺ − ½)/Λ, almost not at all, and their f32 rounding builds up (velocity
/// 2.0e-4 at Re_τ 180, 1.5e-3 at Re_τ 1000, on an RTX 4070 Laptop), where the multiple-relaxation-time and hybrid collisions meet the
/// bounds down to τ − ½ 1.15e-4.
#[test]
fn the_sa_channel_follows_the_cpu() {
    let cases = [Collision::Trt, Collision::Mrt]
        .into_iter()
        .flat_map(|c| [Streaming::InPlace, Streaming::TwoBuffers].map(move |s| (c, s)))
        .map(|(c, s)| (c, s, 20.0))
        .chain([(Collision::Hrr, Streaming::InPlace, 20.0), (Collision::Hrr, Streaming::InPlace, 1000.0)]);
    let mut failed = Vec::new();
    for (collision, streaming, re) in cases {
        let (setup, mut cpu) = prepare(&forced_channel(collision, 2000.0, re), &[]);
        cpu.set_sa_substeps(Some(2));
        let Some(mut gpu) = gpu(&cpu, streaming) else { return };
        gpu.fix_sa_substeps(Some(2)).unwrap();
        let tail = 200;
        cpu.advance(setup.steps - tail).unwrap();
        gpu.advance(setup.steps - tail).unwrap();
        // Each wall's shear, the instantaneous and its mean over the last 200 steps.
        let walls = [BoxFace::YMin.index(), BoxFace::YMax.index()];
        let (mut mc, mut mg, mut last) = ([0.0f64; 2], [0.0f64; 2], 0.0f64);
        for _ in 0..tail {
            cpu.advance(1).unwrap();
            gpu.advance(1).unwrap();
            let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
            last = 0.0;
            for (k, &w) in walls.iter().enumerate() {
                mc[k] += fc[w][0] / tail as f64;
                mg[k] += fg[w][0] / tail as f64;
                last = last.max((fc[w][0] - fg[w][0]).abs() / fc[w][0].abs().max(1e-300));
            }
        }
        let dshear = (0..2).map(|k| (mc[k] - mg[k]).abs() / mc[k].abs().max(1e-300)).fold(0.0, f64::max);
        let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
        println!(
            "SA channel ({collision:?}, {streaming:?}), Re_τ {re}, 2000 steps: ν̃ {dv:.2e} of {scale:.3e}; velocity {du:.2e}, density {drho:.2e}; ν_t/ν {rc:.4} vs {rg:.4}; wall shear, mean of 200 steps {dshear:.2e} (last step {last:.2e}); CPU {mc:?} GPU {mg:?}"
        );
        if !(rc > 0.0 && ((rg - rc) / rc).abs() < 1e-3 && dv < 5e-5 && du < 1e-4 && drho < 1e-5 && dshear < 1e-3) {
            failed.push(format!("{collision:?} {streaming:?} Re_τ {re}"));
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}

/// The RANS channel at Re_τ 1000 (τ − ½ 1.15e-4, the gated cases' τ) over 20 000 steps, the hybrid
/// collision (σ 0.99) and the multiple-relaxation-time one in place: whether their parity stays
/// bounded past the 2000 steps of [`the_sa_channel_follows_the_cpu`], with that test's bounds at
/// every 5000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_sa_channel_stays_with_the_cpu_for_20000_steps() {
    let mut failed = Vec::new();
    for collision in [Collision::Hrr, Collision::Mrt] {
        let (_, mut cpu) = prepare(&forced_channel(collision, 20_000.0, 1000.0), &[]);
        cpu.set_sa_substeps(Some(2));
        let Some(mut gpu) = gpu(&cpu, Streaming::InPlace) else { return };
        gpu.fix_sa_substeps(Some(2)).unwrap();
        let walls = [BoxFace::YMin.index(), BoxFace::YMax.index()];
        let (tail, mut done) = (200u64, 0u64);
        for mark in [5_000u64, 10_000, 15_000, 20_000] {
            cpu.advance(mark - tail - done).unwrap();
            gpu.advance(mark - tail - done).unwrap();
            let (mut mc, mut mg) = ([0.0f64; 2], [0.0f64; 2]);
            for _ in 0..tail {
                cpu.advance(1).unwrap();
                gpu.advance(1).unwrap();
                let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
                for (k, &w) in walls.iter().enumerate() {
                    mc[k] += fc[w][0] / tail as f64;
                    mg[k] += fg[w][0] / tail as f64;
                }
            }
            done = mark;
            let dshear = (0..2).map(|k| (mc[k] - mg[k]).abs() / mc[k].abs().max(1e-300)).fold(0.0, f64::max);
            let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
            let (du, drho) = field_difference(&mut cpu, &mut gpu);
            let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
            println!(
                "SA channel ({collision:?}), Re_τ 1000, {mark} steps: ν̃ {dv:.2e} of {scale:.3e}; velocity {du:.2e}, density {drho:.2e}; ν_t/ν {rc:.4} vs {rg:.4}; wall shear, mean of 200 steps {dshear:.2e}"
            );
            if !(rc > 0.0 && ((rg - rc) / rc).abs() < 1e-3 && dv < 5e-5 && du < 1e-4 && drho < 1e-5 && dshear < 1e-3) {
                failed.push(format!("{collision:?} at {mark}"));
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}

/// The wall distance reaches the device as the CPU's rounded to f32, bit for bit (a channel with
/// a ball, an inlet, an outlet and a slip face); and two runs of the model agree bit for bit.
#[test]
fn the_wall_distance_reaches_the_device_and_runs_repeat() {
    let (m, ball) = ball_channel(Collision::Trt, None, Some(sa(3.0, None)));
    let (setup, flow) = prepare(&m, std::slice::from_ref(&ball));
    let cpu = flow.sa_wall_distance().unwrap().to_vec();
    let run = || {
        let mut g = gpu(&flow, Streaming::InPlace)?;
        let d = g.sa_wall_distance().unwrap().unwrap();
        g.advance(setup.steps.min(200)).unwrap();
        let bits: Vec<u64> =
            g.nu_tilde().unwrap().unwrap().iter().chain(g.fields().unwrap().velocity.iter().flatten()).map(|v| v.to_bits()).collect();
        Some((d, bits, g.sa_substeps()))
    };
    let Some((d, first, n)) = run() else { return };
    assert!(cpu.iter().zip(&d).all(|(c, g)| f64::from(*c as f32).to_bits() == g.to_bits()), "the wall distance on the device");
    let Some((_, second, _)) = run() else { return };
    println!("SA ball channel on the GPU: {} cells, {n:?} sub-steps; two runs equal bit for bit", d.len());
    assert!(first == second, "two runs differ");
}

/// The planned sub-step count covers the chunk's monotone bound: the first chunk's is probed from
/// the initial field, a hot inflow on a coarse grid takes more than one, and a field past the
/// cap is refused with the CPU's wording, never a NaN.
#[test]
fn the_sub_step_count_is_planned_or_refused() {
    let (setup, flow) = prepare(&inflow_channel(Collision::Trt, 64.0, 3.0), &[]);
    let Some(mut g) = gpu(&flow, Streaming::InPlace) else { return };
    assert_eq!(g.sa_substeps(), None);
    g.advance(setup.steps).unwrap();
    let planned = g.sa_substeps().unwrap();
    println!("planned sub-steps after one chunk: {planned}");
    assert!((1..=64).contains(&planned));
    drop(g);
    let (setup, mut hot) = prepare(&inflow_channel(Collision::Trt, 64.0, 3.0), &[]);
    let nu = (setup.tau_plus - 0.5) / 3.0;
    hot.set_sa_nu_tilde(&vec![30.0 * nu; hot.grid().len()]).unwrap();
    let Some(mut g) = gpu(&hot, Streaming::InPlace) else { return };
    let cpu_hot = hot.advance(setup.steps).map_err(|e| e.to_string());
    println!("a hot field on the CPU: {cpu_hot:?}");
    match g.advance(setup.steps) {
        Ok(()) => {
            let n = g.sa_substeps().unwrap();
            println!("a hot field takes {n} sub-steps");
            assert!(n > planned, "{n}");
            assert!(g.nu_tilde().unwrap().unwrap().iter().all(|v| v.is_finite()));
        }
        Err(e) => {
            let e = e.to_string();
            println!("a hot field: {e}");
            assert!(e.contains("Spalart–Allmaras"), "{e}");
        }
    }
    drop(g);
    let (_, mut runaway) = prepare(&inflow_channel(Collision::Trt, 64.0, 3.0), &[]);
    runaway.set_sa_nu_tilde(&vec![1e4 * nu; runaway.grid().len()]).unwrap();
    let Some(mut g) = gpu(&runaway, Streaming::InPlace) else { return };
    let e = g.advance(1).unwrap_err().to_string();
    println!("a runaway field: {e}");
    assert!(e.contains("more than 64 sub-steps"), "{e}");
}

/// What the device does not run with the model yet is refused before a device is opened, naming
/// the model, the wall treatment and the CPU: Malaspinas and Sagaut's wall model (`musker`) and
/// heat. Han's on-node wall runs with it (the halfway wall-function bounce is refused by the
/// setup on either backend).
#[test]
fn what_the_gpu_does_not_run_with_spalart_allmaras_is_refused_clearly() {
    let refused = |flow: &Flow| {
        GpuFlow::with_options(flow, Policy::from_env().expect("FIELDCRAFT_GPU names a policy"), Options::default()).err().map(|e| e.to_string())
    };
    let mut channel = forced_channel(Collision::Mrt, 10.0, 180.0);
    let (model, name) = (WallModel::Musker, "`musker`");
    let mut step = channel.step("Run").cloned().unwrap();
    if let Procedure::Fluid { turbulence, .. } = &mut step.procedure {
        *turbulence = Some(sa(3.0, Some(model)));
    }
    channel.replace_step(step).unwrap();
    let (_, flow) = prepare(&channel, &[]);
    let e = refused(&flow);
    println!("{model:?}: {e:?}");
    assert!(e.as_deref().is_some_and(|e| e.contains("Spalart–Allmaras") && e.contains(name) && e.contains("cpu")), "{e:?}");
    let mut d = domain([0.0; 3], [4.0, 2.0, 2.0], 1.0);
    d.thermal = Some(FluidThermal::default());
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.01, 0.0, 0.0] };
    let m = hot_model(1.0, 0.1, 0.05, d, with(collide(10.0, 1.0, None, Collision::Trt), sa(3.0, None)));
    let (_, flow) = prepare(&m, &[]);
    let e = refused(&flow);
    println!("heat: {e:?}");
    assert!(e.as_deref().is_some_and(|e| e.contains("Spalart–Allmaras") && e.contains("heat") && e.contains("cpu")), "{e:?}");
}

/// The model on Han, Ooka and Kikumoto's on-node wall (`spalding_bounce`), inflow ν̃/ν 3.
fn on_han() -> Turbulence {
    sa(3.0, Some(WallModel::SpaldingBounce))
}

/// The same wall with the model's own constant-stress wall law (`spalart_allmaras_bounce`).
fn on_han_sa_law() -> Turbulence {
    sa(3.0, Some(WallModel::SpalartAllmarasBounce))
}

thread_local! {
    /// Whether the cases below run the model's own wall law instead of Spalding's.
    static OWN_LAW: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Runs `body` with the cases on the model's own wall law.
fn with_own_law<T>(body: impl FnOnce() -> T) -> T {
    OWN_LAW.with(|c| c.set(true));
    let out = body();
    OWN_LAW.with(|c| c.set(false));
    out
}

/// `p` with the model on Han's wall, and the hybrid collision's σ 0.99 when it is the hybrid one.
fn han_procedure(p: Procedure, collision: Collision) -> Procedure {
    let mut p = with(p, if OWN_LAW.with(|c| c.get()) { on_han_sa_law() } else { on_han() });
    if let Procedure::Fluid { sigma, .. } = &mut p
        && collision == Collision::Hrr
    {
        *sigma = Some(HrrSigma::Value(0.99));
    }
    p
}

/// The forced channel on Han's wall, 4 cells per half-height, 12 × 9 × 6 cells, periodic in x and
/// z, started from a blunt profile at the bulk speed 15.7 (lattice 0.1), ν set for τ − ½ =
/// `tau_half`: the walls are the box's y faces, or (`slabs`) obstacle slabs one cell thick in a
/// box periodic in y too. Every node is flat: the plane path, `Side::Node` and the wall function.
fn han_channel(slabs: bool, collision: Collision, tau_half: f64, steps: f64) -> (Setup, Flow) {
    let (n, ub) = (4usize, 15.7);
    let dx = 1.0 / n as f64;
    let (lx, lz) = (12.0 * dx, 6.0 * dx);
    let mut d = if slabs { domain([0.0, -dx, 0.0], [lx, 2.0 + 2.0 * dx, lz], dx) } else { domain([0.0; 3], [lx, 2.0 + dx, lz], dx) };
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = if slabs {
        d.faces.y_min = PERIODIC;
        d.faces.y_max = PERIODIC;
        vec![
            obstacle("Floor", cuboid([-1.0, -dx, -1.0], [lx + 1.0, 0.0, lz + 1.0]), &CUBOID_FACES),
            obstacle("Roof", cuboid([-1.0, 2.0 + dx, -1.0], [lx + 1.0, 2.0 + 2.0 * dx, lz + 1.0]), &CUBOID_FACES),
        ]
    } else {
        Vec::new()
    };
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let nu = tau_half / 3.0 * dx * dx / dt;
    let (setup, mut flow) = prepare(&model(1.0, nu, d, han_procedure(collide(steps * dt, dt, None, collision), collision)), &obstacles);
    let grid = *flow.grid();
    let to_lattice = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|cell| {
        let [_, y, _] = grid.centre(cell);
        let eta = y - 1.0 - 0.5 * dx;
        (1.0, [1.15 * ub * (1.0 - eta.abs().min(1.0).powi(8)) * to_lattice, 0.0, 0.0])
    });
    (setup, flow)
}

/// The junction duct (64 × 16 × 1 cells): a uniform inflow on −x, a convective outlet on +x, Han's
/// walls on y, periodic in z, lattice speed 0.05, τ − ½ = `tau_half`. The four cells where the
/// walls meet the inflow and the outlet keep the open rule, and the four nodes beside them are
/// not flat (their link along the wall to the open cell is not Han's): the general path.
fn han_junction_duct(collision: Collision, tau_half: f64, steps: f64) -> (Setup, Flow) {
    let (n, u_lat) = (16usize, 0.05);
    let dx = 1.0 / n as f64;
    let dt = u_lat * dx;
    let nu = tau_half / 3.0 * dx * dx / dt;
    let mut d = domain([0.0; 3], [4.0, 1.0, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.initial_velocity = [1.0, 0.0, 0.0];
    prepare(&model(1.0, nu, d, han_procedure(collide(steps * dt, dt, None, collision), collision)), &[])
}

/// A block on a forced channel's floor: 32 × 16 × 16 cells at Δ 1/8, periodic in x and z, Han's
/// walls at the y faces, the block 8 × 5 × 4 cells, started at the bulk speed 15.7 (lattice 0.1),
/// τ − ½ = `tau_half`. Its nodes: flat ones on the faces, nodes beside the block's convex edges
/// with one wall direction (not flat), and corner nodes in its concave edges with two (no
/// normal, never sampled); interior cells shared by two sampled nodes, so the wall function
/// reduces rows of several nodes.
fn han_floor_block(collision: Collision, tau_half: f64, steps: f64) -> (Setup, Flow) {
    let (dx, ub) = (1.0 / 8.0, 15.7);
    let mut d = domain([0.0; 3], [4.0, 2.0, 2.0], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = vec![obstacle("Block", cuboid([1.5, -1.0, 0.75], [2.5, 5.0 * dx, 1.25]), &CUBOID_FACES)];
    d.obstacles = vec!["Block".into()];
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [ub, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let nu = tau_half / 3.0 * dx * dx / dt;
    prepare(&model(1.0, nu, d, han_procedure(collide(steps * dt, dt, None, collision), collision)), &obstacles)
}

/// A flow's boundary nodes: how many, how many flat, how many sampled.
fn node_counts(flow: &Flow) -> [usize; 3] {
    let nodes = flow.wall_nodes().map(|w| w.nodes()).unwrap_or_default();
    [nodes.len(), nodes.iter().filter(|w| w.flat).count(), nodes.iter().filter(|w| w.sampled()).count()]
}

/// What one run of the model on Han's wall measured, GPU against CPU: velocity (over the largest
/// CPU speed), density, ν̃ (over the largest CPU ν̃); the wall function's friction velocities at
/// the last step (over the largest CPU one); and every group's force averaged over the last 200
/// steps (over the largest CPU group force).
struct HanParity {
    du: f64,
    drho: f64,
    dv: f64,
    dut: f64,
    dforce: f64,
    rows: usize,
}

/// `steps` steps of `flow` on the CPU and the GPU (`streaming`), both with 2 fixed sub-steps of ν̃,
/// and what [`HanParity`] measures; the wall function's rows are the same cells on both.
fn han_parity(flow: Flow, steps: u64, streaming: Streaming) -> Option<HanParity> {
    let mut cpu = flow;
    cpu.set_sa_substeps(Some(2));
    let mut gpu = gpu(&cpu, streaming)?;
    gpu.fix_sa_substeps(Some(2)).unwrap();
    let tail = 200.min(steps);
    cpu.advance(steps - tail).unwrap();
    gpu.advance(steps - tail).unwrap();
    let groups = cpu.groups().len();
    let (mut mc, mut mg) = (vec![[0.0f64; 3]; groups], vec![[0.0f64; 3]; groups]);
    for _ in 0..tail {
        cpu.advance(1).unwrap();
        gpu.advance(1).unwrap();
        let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
        assert_eq!(fg.len(), groups);
        for g in 0..groups {
            for k in 0..3 {
                mc[g][k] += fc[g][k] / tail as f64;
                mg[g][k] += fg[g][k] / tail as f64;
            }
        }
    }
    let norm = |v: &[f64; 3]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let scale = mc.iter().map(norm).fold(0.0, f64::max).max(1e-300);
    let dforce = mc.iter().zip(&mg).map(|(a, b)| norm(&[a[0] - b[0], a[1] - b[1], a[2] - b[2]])).fold(0.0, f64::max) / scale;
    let (wc, wg) = (cpu.sa_wall_function().1.to_vec(), gpu.sa_wall_function().unwrap());
    assert_eq!(wc.iter().map(|r| r.0).collect::<Vec<_>>(), wg.iter().map(|r| r.0).collect::<Vec<_>>(), "the wall function's rows");
    let ut = wc.iter().map(|r| r.2).fold(0.0, f64::max).max(1e-300);
    let dut = wc.iter().zip(&wg).map(|(a, b)| (a.2 - b.2).abs()).fold(0.0, f64::max) / ut;
    let (dv, _) = nu_tilde_difference(&mut cpu, &mut gpu);
    let (du, drho) = field_difference(&mut cpu, &mut gpu);
    Some(HanParity { du, drho, dv, dut, dforce, rows: wc.len() })
}

/// The bounds: velocity 1e-4, density 1e-5, ν̃ 5e-5, the friction velocities the wall function
/// asked for and the groups' 200-step mean forces 1e-3. Prints what it measured; false past any.
fn han_check(what: &str, r: &HanParity) -> bool {
    println!(
        "{what}: velocity {:.2e}, density {:.2e}, ν̃ {:.2e}; wall function u_τ {:.2e} over {} rows; group forces, mean of 200 steps {:.2e}",
        r.du, r.drho, r.dv, r.dut, r.rows, r.dforce
    );
    r.du < 1e-4 && r.drho < 1e-5 && r.dv < 5e-5 && r.dut < 1e-3 && r.dforce < 1e-3
}

/// The model on Han's on-node wall, GPU against CPU, at τ − ½ 1.2e-4 (the gated cases' τ), both
/// layouts: the forced channel with box-face walls and with obstacle slabs (144 flat nodes and
/// 144 wall-function rows each), on the hybrid collision (σ 0.99) and the multiple-relaxation-time
/// one, 2000 steps; and the junction duct (120 flat and 4 edge nodes, 124 rows) on the hybrid
/// collision, 2000 steps (with the multiple-relaxation-time collision the CPU itself goes unstable
/// at this τ, at step 255), and on both at τ 0.6.
fn the_model_on_han_s_wall_follows_the_cpu_cases(steps: f64) {
    let mut failed = Vec::new();
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        for slabs in [false, true] {
            for collision in [Collision::Hrr, Collision::Mrt] {
                let (setup, flow) = han_channel(slabs, collision, 1.2e-4, steps);
                assert_eq!(node_counts(&flow), [144, 144, 144]);
                let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
                let walls = if slabs { "obstacle slabs" } else { "box faces" };
                let what = format!("SA on Han, channel ({walls}), {collision:?}, {streaming:?}, τ − ½ 1.2e-4, {steps} steps");
                assert_eq!(r.rows, 144);
                if !han_check(&what, &r) {
                    failed.push(what);
                }
            }
        }
        for (collision, tau_half) in [(Collision::Hrr, 1.2e-4), (Collision::Hrr, 0.1), (Collision::Mrt, 0.1)] {
            let (setup, flow) = han_junction_duct(collision, tau_half, steps);
            assert_eq!(node_counts(&flow), [124, 120, 124]);
            let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
            let what = format!("SA on Han, junction duct, {collision:?}, {streaming:?}, τ − ½ {tau_half}, {steps} steps");
            assert_eq!(r.rows, 124);
            if !han_check(&what, &r) {
                failed.push(what);
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}

#[test]
fn the_model_on_han_s_wall_follows_the_cpu() {
    the_model_on_han_s_wall_follows_the_cpu_cases(300.0);
}

/// [`the_model_on_han_s_wall_follows_the_cpu`] over 2000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_model_on_han_s_wall_follows_the_cpu_for_2000_steps() {
    the_model_on_han_s_wall_follows_the_cpu_cases(2000.0);
}

/// The forced channel on the half-way wall: [`han_channel`]'s box (4 cells per half-height, 12 × 8 × 6
/// open cells, periodic in x and z) between two obstacle slabs two cells thick, periodic in y too,
/// with `wall` (`spalding_bounce_surface` or `spalart_allmaras_bounce_surface`). Every node is flat
/// (72 per wall); the law reads the second open cell off each wall (y = 1½) and ν̃ is imposed at
/// the first two.
fn half_way_channel(wall: WallModel, collision: Collision, tau_half: f64, steps: f64) -> (Setup, Flow) {
    let (n, ub) = (4usize, 15.7);
    let dx = 1.0 / n as f64;
    let (lx, lz) = (12.0 * dx, 6.0 * dx);
    let mut d = domain([0.0, -2.0 * dx, 0.0], [lx, 2.0 + 2.0 * dx, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = vec![
        obstacle("Floor", cuboid([-1.0, -2.0 * dx, -1.0], [lx + 1.0, 0.0, lz + 1.0]), &CUBOID_FACES),
        obstacle("Roof", cuboid([-1.0, 2.0, -1.0], [lx + 1.0, 2.0 + 2.0 * dx, lz + 1.0]), &CUBOID_FACES),
    ];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let nu = tau_half / 3.0 * dx * dx / dt;
    let mut p = with(collide(steps * dt, dt, None, collision), sa(3.0, Some(wall)));
    if let Procedure::Fluid { sigma, .. } = &mut p
        && collision == Collision::Hrr
    {
        *sigma = Some(HrrSigma::Value(0.99));
    }
    let (setup, mut flow) = prepare(&model(1.0, nu, d, p), &obstacles);
    let grid = *flow.grid();
    let to_lattice = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|cell| {
        let [_, y, _] = grid.centre(cell);
        let eta = y - 1.0;
        (1.0, [1.15 * ub * (1.0 - eta.abs().min(1.0).powi(8)) * to_lattice, 0.0, 0.0])
    });
    (setup, flow)
}

/// The half-way wall (law at the second open cell, ν̃ at the first two), GPU against CPU: the forced
/// channel between slabs with the model's own law and with Spalding's, on both collisions and both
/// streamings, 2000 steps at τ − ½ 1.2e-4, and at τ 0.6 with the multiple-relaxation-time collision.
fn the_half_way_wall_follows_the_cpu_cases(steps: f64) {
    let mut failed = Vec::new();
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        for wall in [WallModel::SpalartAllmarasBounceSurface, WallModel::SpaldingBounceSurface] {
            for (collision, tau_half) in [(Collision::Hrr, 1.2e-4), (Collision::Mrt, 1.2e-4), (Collision::Mrt, 0.1)] {
                let (setup, flow) = half_way_channel(wall, collision, tau_half, steps);
                assert_eq!(node_counts(&flow), [144, 144, 144]);
                assert!(flow.wall_nodes().is_some_and(|n| n.second_cell() && n.nodes().iter().all(|w| w.sample_cell != w.interior)));
                // Half the steps in, the fields agree to 2e-6 of the largest speed (f32 at τ − ½ 1.2e-4).
                let early = han_parity(half_way_channel(wall, collision, tau_half, steps / 2.0).1, (steps / 2.0) as u64, streaming);
                let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
                let what = format!("half-way wall {wall:?}, {collision:?}, {streaming:?}, τ − ½ {tau_half}, {steps} steps");
                if let Some(e) = early {
                    println!("{what}: first {} steps: velocity {:.2e}, density {:.2e}", steps / 2.0, e.du, e.drho);
                    if e.du > 2e-6 || e.drho > 1e-6 {
                        failed.push(format!("{what}: pointwise at half the steps {:.2e} / {:.2e}", e.du, e.drho));
                    }
                }
                assert_eq!(r.rows, 144);
                if !han_check(&what, &r) {
                    failed.push(what);
                }
            }
        }
    }
    assert!(failed.is_empty(), "{failed:?}");
}

#[test]
fn the_half_way_wall_follows_the_cpu() {
    the_half_way_wall_follows_the_cpu_cases(300.0);
}

/// [`the_half_way_wall_follows_the_cpu`] over 2000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_half_way_wall_follows_the_cpu_for_2000_steps() {
    the_half_way_wall_follows_the_cpu_cases(2000.0);
}

/// A passage two cells wide between two slabs: each wall node's second open cell is the opposite
/// wall's first, and the device sets every cell's ν̃ once from the law at y = ½ as the CPU does
/// (lattice viscosity 0.05, a mean speed of 0.02).
#[test]
fn a_two_cell_passage_follows_the_cpu() {
    let dx = 0.5;
    let span = 4.0 * dx;
    let mut d = domain([0.0, -2.0 * dx, 0.0], [span, 1.0 + 2.0 * dx, span], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let obstacles = vec![
        obstacle("Floor", cuboid([-1.0, -2.0 * dx, -1.0], [span + 1.0, 0.0, span + 1.0]), &CUBOID_FACES),
        obstacle("Roof", cuboid([-1.0, 1.0, -1.0], [span + 1.0, 1.0 + 2.0 * dx, span + 1.0]), &CUBOID_FACES),
    ];
    d.obstacles = obstacles.iter().map(|o| o.body.clone()).collect();
    let dt = 0.05 * dx * dx;
    d.acceleration = [12.0 * 0.02 * dx / dt, 0.0, 0.0];
    let p = with(collide(100.0 * dt, dt, None, Collision::Mrt), sa(3.0, Some(WallModel::SpalartAllmarasBounceSurface)));
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (_, flow) = prepare(&model(1.0, 1.0, d.clone(), p.clone()), &obstacles);
        let Some(r) = han_parity(flow, 100, streaming) else { return };
        assert!(han_check(&format!("two-cell passage, {streaming:?}"), &r));
    }
}

/// The device reports what the half-way wall made of the obstacles' walls as the CPU does, so the
/// solve's warnings (narrow passages, unsealed surface) are the same on both backends.
#[test]
fn the_gpu_reports_the_half_way_walls_layer() {
    let (_, flow) = half_way_channel(WallModel::SpaldingBounceSurface, Collision::Mrt, 0.1, 10.0);
    let Some(gpu) = gpu(&flow, Streaming::InPlace) else { return };
    let want = Backend::wall_layer(&flow);
    assert!(want.is_some_and(|l| l.nodes == 144));
    assert_eq!(Backend::wall_layer(&gpu.flow), want);
}

/// The model on Han's wall with its own constant-stress law, GPU against CPU (the same cases as
/// above, and the table in f32 on the device): the forced channel with box-face walls and with
/// obstacle slabs, on both collisions, 2000 steps at τ − ½ 1.2e-4 where the friction velocity
/// asked for runs from the sublayer through the buffer layer; the junction duct on the hybrid
/// collision at both τ; the block on a floor (rows of several nodes) for 400 steps.
fn the_model_s_own_wall_law_follows_the_cpu_cases(steps: f64) {
    let mut failed = Vec::new();
    with_own_law(|| {
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            for slabs in [false, true] {
                for collision in [Collision::Hrr, Collision::Mrt] {
                    let (setup, flow) = han_channel(slabs, collision, 1.2e-4, steps);
                    assert_eq!(flow.wall_nodes().map(|n| n.law()), Some(fieldcraft_lbm::sa::WallLaw::SpalartAllmaras));
                    let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
                    let walls = if slabs { "obstacle slabs" } else { "box faces" };
                    let what = format!("SA own wall law, channel ({walls}), {collision:?}, {streaming:?}, τ − ½ 1.2e-4, {steps} steps");
                    assert_eq!(r.rows, 144);
                    if !han_check(&what, &r) {
                        failed.push(what);
                    }
                }
            }
            for (collision, tau_half) in [(Collision::Hrr, 1.2e-4), (Collision::Hrr, 0.1)] {
                let (setup, flow) = han_junction_duct(collision, tau_half, steps);
                let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
                let what = format!("SA own wall law, junction duct, {collision:?}, {streaming:?}, τ − ½ {tau_half}, {steps} steps");
                assert_eq!(r.rows, 124);
                if !han_check(&what, &r) {
                    failed.push(what);
                }
            }
            let (setup, flow) = han_floor_block(Collision::Hrr, 1.2e-4, 400.0);
            let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
            let what = format!("SA own wall law, block on a floor, Hrr, {streaming:?}, τ − ½ 1.2e-4, 400 steps");
            assert_eq!(r.rows, 1072);
            if !han_check(&what, &r) {
                failed.push(what);
            }
        }
    });
    assert!(failed.is_empty(), "{failed:?}");
}

#[test]
fn the_model_s_own_wall_law_follows_the_cpu() {
    the_model_s_own_wall_law_follows_the_cpu_cases(300.0);
}

/// [`the_model_s_own_wall_law_follows_the_cpu`] over 2000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_model_s_own_wall_law_follows_the_cpu_for_2000_steps() {
    the_model_s_own_wall_law_follows_the_cpu_cases(2000.0);
}

/// The block on a channel floor `steps` steps, GPU against CPU, at τ − ½ 1.2e-4 on the hybrid
/// (σ 0.99) and multiple-relaxation-time collisions, both layouts: 984 flat, 112 edge and 24
/// corner nodes, 1096 sampled, 1072 wall-function rows. The CPU's stresses have no component
/// along a node's normal axis (the sample takes the velocity parallel to the wall), which the
/// device's general path relies on.
fn floor_block_cases(steps: f64) {
    let mut failed = Vec::new();
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        for collision in [Collision::Hrr, Collision::Mrt] {
            let (setup, flow) = han_floor_block(collision, 1.2e-4, steps);
            assert_eq!(node_counts(&flow), [1120, 984, 1096]);
            let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
            let what = format!("SA on Han, block on a floor, {collision:?}, {streaming:?}, τ − ½ 1.2e-4, {steps} steps");
            assert_eq!(r.rows, 1072);
            if !han_check(&what, &r) {
                failed.push(what);
            }
        }
    }
    let (setup, mut cpu) = han_floor_block(Collision::Hrr, 1.2e-4, steps);
    cpu.set_sa_substeps(Some(2));
    cpu.advance(setup.steps).unwrap();
    for w in cpu.wall_nodes().unwrap().nodes() {
        assert_eq!(w.stress[usize::from(w.axis)], 0.0, "node {}", w.cell);
    }
    assert!(failed.is_empty(), "{failed:?}");
}

/// [`floor_block_cases`] over 400 steps, the trunk gate's row.
#[test]
fn the_model_on_han_s_wall_follows_the_cpu_around_a_block() {
    floor_block_cases(400.0);
}

/// [`floor_block_cases`] over 2000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_model_on_han_s_wall_follows_the_cpu_around_a_block_for_2000_steps() {
    floor_block_cases(2000.0);
}

/// The model on Han's wall twice on one device gives the same fields, ν̃, forces and wall-function
/// rows bit for bit: the block on a floor (the general path, rows of several nodes), 100 steps,
/// hybrid collision, both layouts.
#[test]
fn the_model_on_han_s_wall_is_deterministic_on_its_device() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (_, mut flow) = han_floor_block(Collision::Hrr, 1.2e-4, 100.0);
        flow.set_sa_substeps(Some(2));
        let run = || {
            let mut g = gpu(&flow, streaming)?;
            g.fix_sa_substeps(Some(2)).unwrap();
            g.advance(100).unwrap();
            let mut bits = field_bits(&g.fields().unwrap());
            bits.extend(g.nu_tilde().unwrap().unwrap().iter().map(|v| v.to_bits()));
            bits.extend(g.forces().unwrap().iter().flatten().map(|v| v.to_bits()));
            bits.extend(g.sa_wall_function().unwrap().iter().map(|r| r.2.to_bits()));
            Some(bits)
        };
        let (Some(a), Some(b)) = (run(), run()) else { return };
        println!("SA on Han, block on a floor, {streaming:?}: two runs bit for bit: {}", a == b);
        assert!(a == b);
    }
}
