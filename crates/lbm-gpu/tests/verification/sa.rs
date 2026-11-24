//! The Spalart–Allmaras model on the GPU against the CPU. Each bar is derived from f32 arithmetic (Higham's γ_k and
//! the constants below: the errors of N sub-steps add, to first order, since the perturbation obeys the same stencil,
//! whose coefficients sum to at most 1, so a sub-step does not amplify it) or is a stated tolerance of this suite,
//! named at its use; none is set from a measurement. The wall law's bars also take WGSL's stated accuracy of `atan2` and `log` (W3C WGSL §15.7.4.1), which IEEE
//! arithmetic does not cover. The wall distance on the device is the CPU's rounded to f32, bit
//! for bit; two runs agree bit for bit. The parity tests run the same fixed sub-step count
//! (`Flow::set_sa_substeps`, `GpuFlow::fix_sa_substeps`), so the sub-step splitting is not part of the difference;
//! `the_device_sub_steps_within_the_bound_as_the_cpu_does` runs the device's own sub-stepping against the CPU's.
//! Each test prints what it measured.

use super::*;
use fieldcraft_lbm::sa::KAPPA;
use fieldcraft_lbm::sa_law::{factors, u_plus};
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

/// What the device does not run yet is refused before a device is opened, naming the model
/// and the CPU: heat. (The setup refuses every wall model but the model's own law on both backends.)
#[test]
fn what_the_gpu_does_not_run_with_spalart_allmaras_is_refused_clearly() {
    let refused = |flow: &Flow| {
        GpuFlow::with_options(flow, Policy::from_env().expect("FIELDCRAFT_GPU names a policy"), Options::default()).err().map(|e| e.to_string())
    };
    let mut d = domain([0.0; 3], [4.0, 2.0, 2.0], 1.0);
    d.thermal = Some(FluidThermal::default());
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.01, 0.0, 0.0] };
    let m = hot_model(1.0, 0.1, 0.05, d, with(collide(10.0, 1.0, None, Collision::Trt), sa(3.0, None)));
    let (_, flow) = prepare(&m, &[]);
    let e = refused(&flow);
    println!("heat: {e:?}");
    assert!(e.as_deref().is_some_and(|e| e.contains("Spalart–Allmaras") && e.contains("heat") && e.contains("cpu")), "{e:?}");
}

// ---- The model's own wall law (`wall_model: "spalart_allmaras_bounce"`) ----

/// WGSL's accuracy of `atan2` in f32, in units in the last place (W3C, *WebGPU Shading Language*, §15.7.4.1,
/// "Accuracy of concrete floating point built-in functions": `atan2(y, x)`: "4096 ULP"). An ULP is at most 2 u of
/// the value. The law's two arctangents are at most the whole of F = κ u⁺ (its logarithms' part is positive).
const ATAN2_ULPS: f64 = 4096.0;
/// The law's builtin logarithms (κ y⁺ ≥ 1), in u of max(1, |log|) each: WGSL's `log` is within 3 ULP (6 u) outside
/// [0.5, 2] and within 2⁻²¹ (8 u) inside (the same table): 8; the argument 1 + (p + q/x)/x of the second is within
/// 8 u of itself: 8; their two additions: 2.
const LOG_UNITS: f64 = 18.0;
/// The rounded operations on the longest of the law's four terms (lbm.wgsl, `sa_law_factor`), the arctangent's
/// c θ. d = √(4q − p²) is within 4.4 u (q and p rounded once, p², the difference with the amplification
/// (4q + 3p²)/(4q − p²) ≤ 5.8, a rounding, halved by the root, a rounding). c = (2B − A p)/d is within 14 u (the
/// numerator within 8.2 u, its amplification (2|B| + 3|A p|)/|2B − A p| ≤ 7.2 and a rounding; d; the quotient). The
/// angle θ = atan2(Y, X), Y = 2x/d within 5.4 u and X = 1 + p(2x + p)/d² within (15 + 14.8|X|) u, moves by at most
/// u|Y|(20.2|X| + 15)/(X² + Y²), which is at most 36 u of θ (θ ≥ XY/(X² + Y²) where Y < X, and there X ≥ 0.965;
/// θ ≥ π/4 elsewhere, and there Y ≥ 0.965). Their product: 1.
const OPS_LAW_TERM: f64 = 51.0;
/// The operations common to the four terms: x = κ y⁺ (κ's rounding and the product), the terms' sums (2), F / κ (2).
const OPS_LAW_COMMON: f64 = 6.0;
/// From u⁺ to u_τ: the residual y⁺u⁺ − Re's product and the stop within an ulp of y⁺ (3); Re = U y / ν (ν's
/// rounding, a product, a quotient: 3); u_τ = ν y⁺ / y (3).
const OPS_INVERSION: f64 = 9.0;
/// The interior node's speed is evaluated twice on the device, in `bounce_sample` and in the readback the gap is
/// taken from. Each sums 19 populations (18), adds the half force (2) and divides (1) on terms whose sizes sum to
/// at most twice the momentum (the populations are deviations from the weights), 42, and takes the tangential
/// length (3): 2 · 45.
const OPS_SAMPLE: f64 = 90.0;

/// The relative error bound of the device's u⁺ (lbm.wgsl, `sa_law_u_plus`), and so of the root y⁺ of y⁺u⁺ = Re
/// (the product's log-derivative is at least 1) and of u_τ: the arctangents' allowance on a share of at most 1;
/// the logarithms' on Σ ½|A| (2 m(ln x) + m(ln(1 + (p + q/x)/x)) + m(ln q)) / F, m = max(1, |·|), which is largest
/// at x = κ y⁺ = 1; and the rounded operations of a term, the terms' sizes summing to at most (|A₁| + |A₂|) F.
fn law_bar() -> f64 {
    let f = factors();
    let size = |v: f64| v.abs().max(1.0);
    // At x = 1: ln x = 0, the second logarithm's argument is 1 + p + q, and F = κ u⁺(1/κ).
    let logs =
        f.p.iter().zip(&f.q).zip(&f.a).map(|((p, q), a)| 0.5 * a.abs() * (2.0 * size(0.0) + size((1.0 + p + q).ln()) + size(q.ln()))).sum::<f64>()
            / (KAPPA * u_plus(1.0 / KAPPA));
    let terms = f.a.iter().map(|a| a.abs()).sum::<f64>();
    2.0 * ATAN2_ULPS * UNIT_ROUND_OFF + LOG_UNITS * UNIT_ROUND_OFF * logs + gamma(OPS_LAW_TERM * terms + OPS_LAW_COMMON)
}

/// The CPU's wall-law column (`crates/lbm/tests/wall_bounce.rs`, `sa_law_channel`): 1 × 21 × 1 cells, 10 per δ at
/// Re_τ 1000 (the first interior node at y⁺ 100, τ⁺ 0.50015), on-node walls on y, periodic in x and z, a body
/// force, started from the law's own profile, `steps` steps.
fn wall_law_column(collision: Collision, sigma: Option<HrrSigma>, steps: f64) -> (Setup, Flow) {
    let (re_tau, n, ub) = (1000.0, 10usize, 20.0);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [dx, 2.0 + dx, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = 0.1 * dx / ub;
    let p = les(sigma, Some(sa(3.0, Some(WallModel::SpalartAllmarasBounce))))(collide(steps * dt, dt, None, collision));
    let (setup, mut flow) = prepare(&model(1.0, 1.0 / re_tau, d, p), &[]);
    let grid = *flow.grid();
    assert_eq!(grid.cells, [1, 2 * n + 1, 1]);
    let to_lattice = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|cell| {
        // Wall units: u_τ = δ = 1, the walls through the outermost centres.
        let y = grid.centre(cell)[1] - 0.5 * dx;
        (1.0, [u_plus(re_tau * y.min(2.0 - y).max(0.0)) * to_lattice, 0.0, 0.0])
    });
    (setup, flow)
}

/// A slot 1 × 3 × 1 cells: the two walls' nodes share the one interior node (the max face's node sets its ν̃, as
/// the CPU's list in cell order does). ν 0.02 and speed 0.05 in lattice units: y⁺ about 1.6, κ y⁺ below 1, the
/// law's series branch. 10 steps, the two-relaxation-time collision (τ⁺ 0.56).
fn wall_law_slot() -> (Setup, Flow) {
    let mut d = domain([0.0; 3], [1.0, 3.0, 1.0], 1.0);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1e-5, 0.0, 0.0];
    d.initial_velocity = [0.05, 0.0, 0.0];
    let t = sa(3.0, Some(WallModel::SpalartAllmarasBounce));
    prepare(&model(1.0, 0.02, d, with(collide(10.0, 1.0, None, Collision::Trt), t)), &[])
}

/// `flow` through `setup`'s steps on both backends, 2 sub-steps of ν̃ per step on each, and the wall law's bars.
/// After all but the last step the two velocity fields are read: they are what the last step's sample reads. Then,
/// per boundary node, with gap the two backends' difference of the interior node's tangential velocity over the
/// CPU's and bar = [`law_bar`] + γ(OPS_INVERSION + OPS_SAMPLE) + gap (d ln u_τ / d ln U is between ½ and 1):
/// - u_τ within bar of the CPU's;
/// - the stress ρ u_τ² along that velocity within 2 bar + gap + the density's difference + γ₄;
/// - ν̃ at the interior node within γ₄ of κ times the device's own u_τ times the wall distance (κ's rounding, two
///   products, the distance's rounding), which no builtin enters, and within bar + γ₄ of the CPU's;
/// - ν̃ at the node itself exactly 0.
///
/// The fields and the largest ν_t/ν at the suite's stated tolerances (`check_les`), ν̃ everywhere at the eddy
/// viscosity's (1e-3, as `the_sa_channel_follows_the_cpu`), the device's wall force within 1e-3 of its own sampled
/// stress (`check_on_node`'s stated tolerance), and so the walls' drag within 1e-3 plus the largest stress bar of the
/// CPU's. `None` without a GPU.
fn check_wall_law(what: &str, setup: &Setup, flow: Flow, streaming: Streaming) -> Option<()> {
    let mut cpu = flow;
    assert!(cpu.on_node_bounce() && cpu.spalart_allmaras().is_some() && setup.steps > 1, "{what}");
    cpu.set_sa_substeps(Some(2));
    let mut gpu = gpu(&cpu, streaming)?;
    gpu.fix_sa_substeps(Some(2)).unwrap();
    cpu.advance(setup.steps - 1).unwrap();
    gpu.advance(setup.steps - 1).unwrap();
    let (before_c, before_g) = (cpu.fields().unwrap(), gpu.fields().unwrap());
    cpu.advance(1).unwrap();
    gpu.advance(1).unwrap();
    let nodes = cpu.wall_nodes().unwrap().nodes().to_vec();
    let rows = gpu.wall_targets().unwrap();
    assert_eq!(rows.len(), nodes.len(), "{what}");
    let distance = cpu.sa_wall_distance().unwrap().to_vec();
    let (vc, vg) = (cpu.nu_tilde().unwrap().unwrap(), gpu.nu_tilde().unwrap().unwrap());
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let less = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let law = law_bar();
    // The largest of each measure, with its bar where the bar is the node's own.
    let (mut friction, mut stress, mut own, mut imposed, mut stress_bar) = ((0.0f64, 0.0f64), (0.0f64, 0.0f64), 0.0f64, 0.0f64, 0.0f64);
    for (w, row) in nodes.iter().zip(&rows) {
        let (i, axis) = (w.interior as usize, usize::from(w.axis));
        let tangent = |mut v: [f64; 3]| {
            v[axis] = 0.0;
            v
        };
        let (tc, tg) = (tangent(before_c.velocity[i]), tangent(before_g.velocity[i]));
        assert!(w.friction > 0.0 && norm(tc) > 0.0, "{what}: node {} at rest", w.cell);
        let gap = norm(less(tg, tc)) / norm(tc);
        let density = (before_g.density[i] - before_c.density[i]).abs() / before_c.density[i];
        let bar = law + gamma(OPS_INVERSION + OPS_SAMPLE) + gap;
        let bar_s = 2.0 * bar + gap + density + gamma(4.0);
        let (ut, du) = (row.1, (row.1 - w.friction).abs() / w.friction);
        let ds = norm(less(row.0, w.stress)) / norm(w.stress);
        assert!(du <= bar, "{what}: node {} u_τ {ut} vs {} ({du:.2e}, bar {bar:.2e})", w.cell, w.friction);
        assert!(ds <= bar_s, "{what}: node {} stress {:?} vs {:?} ({ds:.2e}, bar {bar_s:.2e})", w.cell, row.0, w.stress);
        let value = KAPPA * ut * distance[i];
        let (d_own, d_cpu) = ((vg[i] - value).abs() / value, (vg[i] - vc[i]).abs() / vc[i]);
        assert!(d_own <= gamma(4.0), "{what}: node {} ν̃ {} vs κ u_τ d {value} ({d_own:.2e})", w.cell, vg[i]);
        assert!(d_cpu <= bar + gamma(4.0), "{what}: node {} ν̃ {} vs {} ({d_cpu:.2e}, bar {bar:.2e})", w.cell, vg[i], vc[i]);
        assert_eq!((vc[w.cell as usize], vg[w.cell as usize]), (0.0, 0.0), "{what}: a boundary node has no ν̃");
        if du > friction.0 {
            friction = (du, bar);
        }
        if ds > stress.0 {
            stress = (ds, bar_s);
        }
        (own, imposed, stress_bar) = (own.max(d_own), imposed.max(d_cpu), stress_bar.max(bar_s));
    }
    let (du, drho) = field_difference(&mut cpu, &mut gpu);
    let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
    let (rc, rg) = (Backend::eddy_viscosity_ratio(&cpu), gpu.eddy_viscosity_ratio());
    let (fc, fg) = (Backend::forces(&mut cpu).unwrap(), gpu.forces().unwrap());
    let walls = |f: &[[f64; 3]]| f[BoxFace::YMin.index()][0] + f[BoxFace::YMax.index()][0];
    let sampled = rows.iter().map(|r| r.0[0]).sum::<f64>();
    let carried = (walls(&fg) - sampled).abs() / sampled.abs();
    let drag = (walls(&fc) - walls(&fg)).abs() / walls(&fc).abs();
    println!(
        "{what}, {} steps, {} nodes: u_τ {:.2e} (bar {:.2e}, the law's {law:.2e}), stress {:.2e} (bar {:.2e}), ν̃ at the interior nodes {own:.2e} of κ u_τ d (bar {:.2e}) and {imposed:.2e} of the CPU's; ν̃ {dv:.2e} of {scale:.3e}; wall force {carried:.2e} of the sampled stress, drag {drag:.2e} of the CPU's",
        setup.steps,
        nodes.len(),
        friction.0,
        friction.1,
        stress.0,
        stress.1,
        gamma(4.0)
    );
    check_les(what, (du, drho, rc, rg));
    assert!(rc > 0.0 && dv < 1e-3, "{what}: ν_t/ν {rc}, ν̃ {dv}");
    assert!(carried < 1e-3 && drag <= 1e-3 + stress_bar, "{what}: carried {carried}, drag {drag} (bar {})", 1e-3 + stress_bar);
    Some(())
}

/// The model's own wall law in the CPU's column at Re_τ 1000 (the first node at y⁺ 100, the law's logarithmic
/// branch), 40 steps from the law's profile, on the regularised and hybrid (σ 0.99) collisions, both layouts, at
/// [`check_wall_law`]'s bars.
#[test]
fn the_wall_law_follows_the_cpu_in_its_column() {
    for (collision, sigma) in [(Collision::Regularized, None), (Collision::Hrr, Some(HrrSigma::Value(0.99)))] {
        for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
            let (setup, flow) = wall_law_column(collision, sigma, 40.0);
            let what = format!("SA wall law, column at Re_τ 1000, {}, {streaming:?}", collision.name());
            if check_wall_law(&what, &setup, flow, streaming).is_none() {
                return;
            }
        }
    }
}

/// The same in three dimensions: the on-node bounce's perturbed coarse channel (38 × 13 × 19 cells, Re_τ 180, the
/// first node at y⁺ 30, periodic in x and z), 150 steps from the start, the regularised collision, both layouts.
#[test]
fn the_wall_law_follows_the_cpu_in_a_box() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (setup, flow) = on_node_channel(sa(3.0, Some(WallModel::SpalartAllmarasBounce)), Collision::Regularized, None, 150.0);
        if check_wall_law(&format!("SA wall law, coarse channel at Re_τ 180, regularized, {streaming:?}"), &setup, flow, streaming).is_none() {
            return;
        }
    }
}

/// The same in a slot 3 cells across, where the two walls share the one interior node and the law is in its
/// sublayer branch, the two-relaxation-time collision, both layouts.
#[test]
fn the_wall_law_follows_the_cpu_where_two_nodes_share_an_interior_node() {
    for streaming in [Streaming::InPlace, Streaming::TwoBuffers] {
        let (setup, flow) = wall_law_slot();
        if check_wall_law(&format!("SA wall law, slot 3 cells across, trt, {streaming:?}"), &setup, flow, streaming).is_none() {
            return;
        }
    }
}

/// The wall law runs on the device, and at rest it gives exactly nothing on both backends: no friction velocity,
/// no stress, ν̃ = 0 at the interior nodes (a column 9 cells high, no force, 10 steps, the two-relaxation-time
/// collision); ν̃ elsewhere at the eddy viscosity's stated tolerance.
#[test]
fn the_wall_law_at_rest_gives_no_stress_on_the_gpu() {
    let mut d = domain([0.0; 3], [1.0, 9.0, 1.0], 1.0);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let m = model(1.0, 0.1, d, with(collide(10.0, 1.0, None, Collision::Trt), sa(3.0, Some(WallModel::SpalartAllmarasBounce))));
    let (setup, mut cpu) = prepare(&m, &[]);
    let Some(mut gpu) = gpu(&cpu, Streaming::InPlace) else { return };
    cpu.advance(setup.steps).unwrap();
    gpu.advance(setup.steps).unwrap();
    let nodes = cpu.wall_nodes().unwrap().nodes().to_vec();
    let rows = gpu.wall_targets().unwrap();
    assert_eq!((nodes.len(), rows.len()), (2, 2));
    let (vc, vg) = (cpu.nu_tilde().unwrap().unwrap(), gpu.nu_tilde().unwrap().unwrap());
    for (w, row) in nodes.iter().zip(&rows) {
        assert_eq!((w.friction, w.stress), (0.0, [0.0; 3]), "{w:?}");
        assert_eq!((row.1, row.0), (0.0, [0.0; 3]), "{row:?}");
        assert_eq!((vc[w.interior as usize], vg[w.interior as usize]), (0.0, 0.0), "{w:?}");
    }
    let (dv, scale) = nu_tilde_difference(&mut cpu, &mut gpu);
    println!("SA wall law at rest, {} steps: u_τ, stress and the interior nodes' ν̃ are 0 on both; ν̃ {dv:.2e} of {scale:.3e}", setup.steps);
    assert!(dv < 1e-3, "{dv}");
}

/// The wall law twice on one device gives the same ν̃, velocities, stresses and friction velocities bit for bit
/// (the column, the regularised collision, in place).
#[test]
fn the_wall_law_repeats_bit_for_bit_on_its_device() {
    let (setup, flow) = wall_law_column(Collision::Regularized, None, 40.0);
    let run = || {
        let mut g = gpu(&flow, Streaming::InPlace)?;
        g.advance(setup.steps).unwrap();
        let mut bits: Vec<u64> =
            g.nu_tilde().unwrap().unwrap().iter().chain(g.fields().unwrap().velocity.iter().flatten()).map(|v| v.to_bits()).collect();
        bits.extend(g.wall_targets().unwrap().iter().flat_map(|r| r.0.into_iter().chain([r.1])).map(|v| v.to_bits()));
        Some(bits)
    };
    let (Some(a), Some(b)) = (run(), run()) else { return };
    println!("SA wall law, column: two runs bit for bit: {}", a == b);
    assert!(a == b);
}

/// A field smaller than 16 bytes reads back: a column of 3 cells has 3 ν̃ (12 bytes), and the readback used to pad
/// every copy to 16 bytes, past the buffer's end, where wgpu's validation panicked. Before any step the device's ν̃
/// is the CPU's rounded to f32, exactly.
#[test]
fn a_field_of_three_cells_reads_back() {
    let dx = 1.0 / 3.0;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    let m = model(1.0, 0.05 * dx, d, with(collide(10.0 * dx, dx, None, Collision::Trt), sa(3.0, None)));
    let (_, mut cpu) = prepare(&m, &[]);
    assert_eq!(cpu.grid().len(), 3);
    let Some(mut gpu) = gpu(&cpu, Streaming::InPlace) else { return };
    let (a, b) = (cpu.nu_tilde().unwrap().unwrap(), gpu.nu_tilde().unwrap().unwrap());
    println!("3-cell column: CPU ν̃ {a:?}, GPU ν̃ {b:?}");
    assert_eq!(b.len(), 3);
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(*y, f64::from(*x as f32), "the device holds the CPU's ν̃ rounded to f32");
    }
}
