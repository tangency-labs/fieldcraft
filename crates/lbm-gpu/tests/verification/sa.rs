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
