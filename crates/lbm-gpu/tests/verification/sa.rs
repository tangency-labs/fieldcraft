//! The Spalart–Allmaras model on the GPU against the CPU. Each bar is derived from f32 arithmetic (Higham's γ_k and
//! the constants below: the errors of N sub-steps add, to first order, since the perturbation obeys the same stencil,
//! whose coefficients sum to at most 1, so a sub-step does not amplify it) or is a stated tolerance of this suite,
//! named at its use; none is set from a measurement. The wall distance on the device is the CPU's rounded to f32, bit
//! for bit; two runs agree bit for bit. The parity tests run the same fixed sub-step count
//! (`Flow::set_sa_substeps`, `GpuFlow::fix_sa_substeps`), so the sub-step splitting is not part of the difference;
//! `the_device_sub_steps_within_the_bound_as_the_cpu_does` runs the device's own sub-stepping against the CPU's.
//! Each test prints what it measured.

use super::*;
use fieldcraft_model::{FaceCondition, WallModel};

/// The unit round-off of f32, 2⁻²⁴ (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM 2002, §2.1).
const UNIT_ROUND_OFF: f64 = f32::EPSILON as f64 / 2.0;

/// Higham's γ_k = k u / (1 − k u) (2002, §3.1, Lemma 3.1): the relative error bound of a result that `k` rounded
/// operations lead to.
fn gamma(k: f64) -> f64 {
    k * UNIT_ROUND_OFF / (1.0 - k * UNIT_ROUND_OFF)
}

/// One sub-step is ν̃′ = ν̃ + dt Σ_f c_f (ν̃_f − ν̃) + dt s (`SaField::rate`, `update`). Its f32 error is at most γ_k
/// times the sum of the terms' sizes, |ν̃| + dt Σ c_f |ν̃_f − ν̃| + dt s. With dt Σ c_f ≤ 1 (the monotone bound) and a
/// field of one sign no larger than S, the first two are at most S each, and dt s = ν̃′ − ν̃ − dt Σ c_f (ν̃_f − ν̃) at
/// most 3 S: 5 S in all.
const SUBSTEP_TERMS: f64 = 5.0;
/// The rounded operations on the longest path of a sub-step with a fixed count of 2 (dt = ½, exact), ν̃ ≥ 0: a
/// face's diffusivity ν + ν̃ (1), the two cells' sum (1), over σ (1) and σ's own rounding (1), plus the upwind speed
/// (1), the difference ν̃_f − ν̃ (1), their product (1), the upwind speed's own mean (1); the six faces' sum (5);
/// the c_b2 term added (1); times dt (1); added to ν̃ (1).
const OPS_TRANSPORT: f64 = 16.0;
/// With ν̃ < 0 the diffusivity is ν + ν̃ f_n(χ): χ (1), χ³ (2), the numerator, the denominator and the quotient (3),
/// times ν̃ (1) in place of nothing, 7 more; and the SA-neg decay ν̃ / (1 + dt (k − p)): the sum and the quotient (2)
/// and the rate's own error, below one more (1).
const OPS_NEGATIVE: f64 = OPS_TRANSPORT + 7.0 + 3.0;
/// With the device's own dt the rounded operations also count the product most · remaining's (1) and dt's own:
/// `remaining` is what k ≤ parts − 1 subtractions left, each within 2⁻²⁵ = u/2 of a number below 1, and dt is at
/// least 1/parts of the step, so dt is within (k parts / 2) u of the CPU's: parts (parts − 1) / 2 at most.
fn ops_adaptive(parts: f64) -> f64 {
    OPS_TRANSPORT + 1.0 + parts * (parts - 1.0) / 2.0
}

/// The largest difference of the two velocity fields. The time step of this file's models is their spacing, so it
/// is a lattice speed; through the upwind coefficients (one per face, six faces) it moves a sub-step's rate by at
/// most six times itself times the largest |ν̃_f − ν̃|.
fn velocity_gap(cpu: &mut Flow, gpu: &mut GpuFlow) -> f64 {
    let (a, b) = (cpu.fields().unwrap(), gpu.fields().unwrap());
    a.velocity
        .iter()
        .zip(&b.velocity)
        .map(|(x, y)| ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt())
        .fold(0.0, f64::max)
}

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
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    model(1.0, 0.05 * dx, d, with(collide(steps * dx, dx, None, collision), sa(ratio, None)))
}

/// A forced channel two high at Re_τ about `re`, six cells per half-height, periodic in x and z
/// (τ⁺ above the laminar floor 0.505 for `re` up to 20), `steps` time steps.
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

/// A duct 24 × 4 × 4 cells: a uniform flow along x from the start (lattice speed 0.05), an inflow on −x, a pressure
/// outlet on +x, periodic in y and z, ν 0.05 in lattice units, inflow ν̃/ν `ratio`. The time step is the spacing, so
/// a velocity is its lattice value.
fn inflow_duct(steps: f64, ratio: f64) -> Model {
    let dx = 1.0 / 4.0;
    let mut d = domain([0.0; 3], [6.0, 1.0, 1.0], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.initial_velocity = [0.05, 0.0, 0.0];
    model(1.0, 0.05 * dx, d, with(collide(steps * dx, dx, None, Collision::Trt), sa(ratio, None)))
}

/// ν̃'s transport alone (the model's sources off on both): a ν̃ front entering a channel from its
/// inflow (ν̃ 0 inside, ν at the inflow), carried by the flow, diffused, with the c_b2 term,
/// past walls, an outlet and periodic faces; within the f32 bar of one step
/// (2·5·γ₁₆ plus six times the velocity gap) and of 200 steps (200 times that), relative to the larger of the largest ν̃ and the inflow's,
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
        let inflow = (setup.tau_plus - 0.5) / 3.0;
        let one_step = 2.0 * SUBSTEP_TERMS * gamma(OPS_TRANSPORT);
        cpu.advance(1).unwrap();
        gpu.advance(1).unwrap();
        let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
        let gap = velocity_gap(&mut cpu, &mut gpu);
        let size = scale.max(inflow);
        let bar = one_step + 6.0 * gap;
        println!("SA transport ({collision:?}, {streaming:?}), 1 step: |Δν̃| {:.2e} of {size:.3e} (bar {bar:.2e})", dv * scale / size);
        assert!(dv * scale <= bar * size, "one step: {} vs {}", dv * scale, bar * size);
        cpu.advance(setup.steps - 1).unwrap();
        gpu.advance(setup.steps - 1).unwrap();
        let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
        let gap = velocity_gap(&mut cpu, &mut gpu);
        let (du, drho) = field_difference(&mut cpu, &mut gpu);
        let size = scale.max(inflow);
        let bar = setup.steps as f64 * (one_step + 6.0 * gap);
        println!(
            "SA transport ({collision:?}, {streaming:?}), 200 steps: |Δν̃| {:.2e} of {size:.3e} (bar {bar:.2e}); velocity {du:.2e}, density {drho:.2e}",
            dv * scale / size
        );
        assert!(dv * scale <= bar * size, "{} vs {}", dv * scale, bar * size);
    }
}

/// The device takes every sub-step within the monotone bound, as the CPU does, where the bound moves inside one
/// submission. A ν̃ front enters an empty duct (the transport alone): the largest coefficient sum goes from
/// 10.5 ν + 3 ν r + u to about 10.5 ν (1 + r) + u (r the inflow ν̃/ν; six faces of (ν + ν̃)/σ, the inflow face twice),
/// so at r = 5 from 1.3 to 3.2 and at r = 20 from 3.6 to 11 within the first 64 steps. Every sub-step's
/// dt × coefficient sum stays at or below 1, a run the CPU completes is completed, and ν̃ follows the CPU's.
#[test]
fn the_device_sub_steps_within_the_bound_as_the_cpu_does() {
    for ratio in [5.0, 20.0] {
        let (setup, mut cpu) = prepare(&inflow_duct(128.0, ratio), &[]);
        cpu.set_sa_model(false);
        cpu.set_sa_nu_tilde(&vec![0.0; cpu.grid().len()]).unwrap();
        let Some(mut gpu) = gpu(&cpu, Streaming::InPlace) else { return };
        gpu.set_sa_model(false).unwrap();
        assert_eq!((gpu.sa_substeps(), gpu.sa_bound()), (None, None));
        let inflow = ratio * (setup.tau_plus - 0.5) / 3.0;
        let (mut rounding, mut seconds) = (0.0, 0.0);
        // Two submissions of 64 steps (native.rs, MAX_CHUNK).
        for chunk in 0..2 {
            cpu.advance(64).unwrap();
            let started = std::time::Instant::now();
            gpu.advance(64).unwrap();
            seconds += started.elapsed().as_secs_f64();
            let (bound, most) = (gpu.sa_bound().unwrap(), gpu.sa_substeps().unwrap());
            println!("SA duct, inflow ν̃/ν {ratio}, chunk {chunk}: largest dt × coefficient sum {bound:.7}, at most {most} sub-steps per step");
            assert!(bound <= 1.0, "inflow ν̃/ν {ratio}, chunk {chunk}: a sub-step at {bound} times the monotone bound");
            assert!((2..=64).contains(&most), "{most}");
            rounding += 64.0 * f64::from(most) * SUBSTEP_TERMS * gamma(ops_adaptive(f64::from(most)));
        }
        let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
        let gap = velocity_gap(&mut cpu, &mut gpu);
        let size = scale.max(inflow);
        let bar = (rounding + 128.0 * 6.0 * gap) * size;
        println!(
            "SA duct, inflow ν̃/ν {ratio}, 128 steps: |Δν̃| {:.2e} of {size:.3e} (bar {:.2e}); velocity gap {gap:.2e}; {:.2} ms per GPU step",
            dv * scale / size,
            bar / size,
            1e3 * seconds / 128.0
        );
        assert!(dv * scale <= bar, "{} vs {bar}", dv * scale);
    }
}

/// SA-neg: a negative ν̃ patch in a periodic box at rest decays towards 0 on both backends and
/// never crosses it; the smallest ν̃ agrees within the f32 bar (400 · 5 · γ₂₆ · 2 ν) after 200 steps.
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
    // 200 steps of 2 sub-steps, the patch's 2 ν the largest |ν̃|.
    let bar = 400.0 * SUBSTEP_TERMS * gamma(OPS_NEGATIVE) * 2.0 * nu_lattice;
    assert!((a - b).abs() <= bar, "{a} vs {b} (bar {bar})");
}

/// The coupled RANS channel (ν_t in the collision) at Re_τ 20 (τ 0.5057): velocity, density,
/// ν̃, the largest ν_t/ν and each wall's shear (its mean over the last 200 steps, as a skin
/// friction is read) follow the CPU over 2000 steps, both layouts, the two-relaxation-time and
/// multiple-relaxation-time collisions, and the hybrid collision in place. A single step's shear
/// is not the measure: the two-relaxation-time collision's acoustic noise moves it 1e-3 between
/// backends whose fields agree to 5e-6 (the mean, 1e-5).
///
/// Not at Re_τ 180 (τ 0.50064): there the backends' own f32 drift, with Vreman's model as much
/// as with this one, is above the suite's stated LES tolerances (`check_les`) (measured on an RTX 4070 Laptop: velocity 2.0e-4
/// with Vreman's, 1.4e-4 with this model; each wall's shear 6 to 40 %, the two walls no longer
/// equal), where at Re_τ 20 both models meet them; the bounds are to be
/// checked at the case's τ.
#[test]
fn the_sa_channel_follows_the_cpu() {
    let cases = [Collision::Trt, Collision::Mrt]
        .into_iter()
        .flat_map(|c| [Streaming::InPlace, Streaming::TwoBuffers].map(move |s| (c, s)))
        .chain([(Collision::Hrr, Streaming::InPlace)]);
    for (collision, streaming) in cases {
        let (setup, mut cpu) = prepare(&forced_channel(collision, 2000.0, 20.0), &[]);
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
            "SA channel ({collision:?}, {streaming:?}), Re_τ 20, 2000 steps: ν̃ {dv:.2e} of {scale:.3e}; velocity {du:.2e}, density {drho:.2e}; ν_t/ν {rc:.4} vs {rg:.4}; wall shear, mean of 200 steps {dshear:.2e} (last step {last:.2e}); CPU {mc:?} GPU {mg:?}"
        );
        check_les(&format!("SA channel ({collision:?}, {streaming:?})"), (du, drho, rc, rg));
        // ν̃ at the eddy viscosity's tolerance (ν_t = ν̃ f_v1 moves at least as fast as ν̃); each wall's mean shear
        // at the wall drag's of the coarse wall-modelled channel (`check_on_node`).
        assert!(rc > 0.0 && dv < 1e-3 && dshear < 1e-3, "{collision:?} {streaming:?}: ν̃ {dv}, shear {dshear}");
    }
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

/// The device takes as many sub-steps as each lattice step needs: none reported before the first chunk, more than one
/// for a hot inflow on a coarse grid, a hot field completed or refused as the CPU does, and a field past the cap
/// refused with the CPU's wording, never a NaN.
#[test]
fn the_sub_steps_follow_the_field_or_are_refused() {
    let (setup, flow) = prepare(&inflow_channel(Collision::Trt, 64.0, 3.0), &[]);
    let Some(mut g) = gpu(&flow, Streaming::InPlace) else { return };
    assert_eq!(g.sa_substeps(), None);
    g.advance(setup.steps).unwrap();
    let planned = g.sa_substeps().unwrap();
    println!("most sub-steps after one chunk: {planned}");
    assert!((1..=64).contains(&planned));
    assert!(g.sa_bound().unwrap() <= 1.0);
    drop(g);
    let (setup, mut hot) = prepare(&inflow_channel(Collision::Trt, 64.0, 3.0), &[]);
    let nu = (setup.tau_plus - 0.5) / 3.0;
    hot.set_sa_nu_tilde(&vec![30.0 * nu; hot.grid().len()]).unwrap();
    let Some(mut g) = gpu(&hot, Streaming::InPlace) else { return };
    let cpu_hot = hot.advance(setup.steps).map_err(|e| e.to_string());
    println!("a hot field on the CPU: {cpu_hot:?}");
    let gpu_hot = g.advance(setup.steps).map_err(|e| e.to_string());
    println!("a hot field on the GPU: {gpu_hot:?}");
    assert_eq!(cpu_hot.is_ok(), gpu_hot.is_ok(), "CPU {cpu_hot:?}, GPU {gpu_hot:?}");
    match gpu_hot {
        Ok(()) => {
            assert!(g.sa_substeps().unwrap() > planned && g.sa_bound().unwrap() <= 1.0);
            assert!(g.nu_tilde().unwrap().unwrap().iter().all(|v| v.is_finite()));
        }
        Err(e) => assert!(e.contains("more than 64 sub-steps"), "{e}"),
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
    let mut channel = forced_channel(Collision::Hrr, 10.0, 180.0);
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
/// `tau_half`: the walls are the box's y faces. Every node is flat: the plane path, `Side::Node`
/// and the wall function.
fn han_channel(collision: Collision, tau_half: f64, steps: f64) -> (Setup, Flow) {
    let (n, ub) = (4usize, 15.7);
    let dx = 1.0 / n as f64;
    let (lx, lz) = (12.0 * dx, 6.0 * dx);
    let mut d = domain([0.0; 3], [lx, 2.0 + dx, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let nu = tau_half / 3.0 * dx * dx / dt;
    let (setup, mut flow) = prepare(&model(1.0, nu, d, han_procedure(collide(steps * dt, dt, None, collision), collision)), &[]);
    let grid = *flow.grid();
    let to_lattice = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|cell| {
        let [_, y, _] = grid.centre(cell);
        let eta = y - 1.0 - 0.5 * dx;
        (1.0, [1.15 * ub * (1.0 - eta.abs().min(1.0).powi(8)) * to_lattice, 0.0, 0.0])
    });
    (setup, flow)
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
/// layouts: the forced channel with box-face walls (144 nodes and 144 wall-function rows) on the
/// hybrid collision (σ 0.99), `steps` steps. With `own_law` the wall is the model's own
/// constant-stress law (the table in f32 on the device), where the friction velocity asked for
/// runs from the sublayer through the buffer layer.
fn the_model_on_han_s_wall_follows_the_cpu_cases(steps: f64, own_law: bool) {
    let mut failed = Vec::new();
    let mut run = || {
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let (setup, flow) = han_channel(Collision::Hrr, 1.2e-4, steps);
            assert_eq!(flow.wall_nodes().map(|n| n.nodes().len()), Some(144));
            let law = if own_law { fieldcraft_lbm::sa::WallLaw::SpalartAllmaras } else { fieldcraft_lbm::sa::WallLaw::Spalding };
            assert_eq!(flow.wall_nodes().map(|n| n.law()), Some(law));
            let Some(r) = han_parity(flow, setup.steps, streaming) else { return };
            let what = format!(
                "SA {}, channel (box faces), Hrr, {streaming:?}, τ − ½ 1.2e-4, {steps} steps",
                if own_law { "own wall law" } else { "on Han" }
            );
            assert_eq!(r.rows, 144);
            if !han_check(&what, &r) {
                failed.push(what);
            }
        }
    };
    if own_law {
        with_own_law(run)
    } else {
        run()
    }
    assert!(failed.is_empty(), "{failed:?}");
}

#[test]
fn the_model_on_han_s_wall_follows_the_cpu() {
    the_model_on_han_s_wall_follows_the_cpu_cases(300.0, false);
}

/// [`the_model_on_han_s_wall_follows_the_cpu`] over 2000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_model_on_han_s_wall_follows_the_cpu_for_2000_steps() {
    the_model_on_han_s_wall_follows_the_cpu_cases(2000.0, false);
}

#[test]
fn the_model_s_own_wall_law_follows_the_cpu() {
    the_model_on_han_s_wall_follows_the_cpu_cases(300.0, true);
}

/// [`the_model_s_own_wall_law_follows_the_cpu`] over 2000 steps.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_model_s_own_wall_law_follows_the_cpu_for_2000_steps() {
    the_model_on_han_s_wall_follows_the_cpu_cases(2000.0, true);
}

/// The model on Han's wall twice on one device gives the same fields, ν̃, forces and wall-function
/// rows bit for bit: the forced channel, 100 steps, hybrid collision, both layouts.
#[test]
fn the_model_on_han_s_wall_is_deterministic_on_its_device() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (_, mut flow) = han_channel(Collision::Hrr, 1.2e-4, 100.0);
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
        println!("SA on Han, channel, {streaming:?}: two runs bit for bit: {}", a == b);
        assert!(a == b);
    }
}
