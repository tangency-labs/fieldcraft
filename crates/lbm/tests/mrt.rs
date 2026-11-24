//! The multiple-relaxation-time collision (`collision: "mrt"`) against the analytic laminar cases:
//! plane Poiseuille and Couette flow and the decaying Taylor–Green vortex, its determinism, and
//! Smagorinsky's eddy viscosity through its stresses' rate.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test mrt -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Limits, Setup};
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const MRT: Collision = Collision::Mrt;

/// Plane Poiseuille flow (walls at y = 0 and 1, driven by a body force, peak 0.1), `n` cells
/// across, diffusive time steps (τ fixed at 0.8): the largest error of the velocity against the
/// parabola, over its peak.
fn poiseuille(n: usize) -> f64 {
    let (h, nu, umax) = (1.0, 0.1, 0.1);
    let dx = h / n as f64;
    let a = 8.0 * nu * umax / (h * h);
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    let (setup, mut flow) = prepare(&model(1.0, nu, d, collide(6.0 * h * h / nu, dx * dx, None, MRT)), &[]);
    run_all(&setup, &mut flow);
    (0..n)
        .map(|j| {
            let y = (j as f64 + 0.5) * dx;
            (flow.moments(j).unwrap().1[0] * setup.velocity_scale() - a / (2.0 * nu) * y * (h - y)).abs() / umax
        })
        .fold(0.0, f64::max)
}

/// With halfway bounce-back the MRT wall lies half a cell out up to an error of order Δ² that
/// depends on τ (the TRT scheme's Λ = 3/16 makes it exact; here the rates are fixed), so the
/// error falls like Δ² at a fixed τ.
#[test]
fn poiseuille_flow_converges_to_the_parabola_like_dx_squared() {
    let errors: Vec<f64> = [8, 16, 32].iter().map(|&n| poiseuille(n)).collect();
    let p = orders(&errors, 2.0);
    println!("MRT Poiseuille, 8 → 32 cells: errors {errors:?}, observed orders {p:.3?}");
    assert!(errors[2] < 2e-3, "{errors:?}");
    assert!(p.iter().all(|o| (o - 2.0).abs() < 0.2), "second order: {p:.3?}");
}

/// Plane Couette flow (a moving wall at y = 1, 16 cells): the linear profile, to the scheme's
/// compressibility error at lattice speed 0.1.
#[test]
fn couette_flow_is_linear() {
    let (h, nu, u_wall, n) = (1.0, 0.05, 0.1, 16);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
    let (setup, mut flow) = prepare(&model(1.0, nu, d, collide(6.0 * h * h / nu, dx, None, MRT)), &[]);
    run_all(&setup, &mut flow);
    let err = (0..n)
        .map(|j| (flow.moments(j).unwrap().1[0] * setup.velocity_scale() - u_wall * (j as f64 + 0.5) * dx / h).abs() / u_wall)
        .fold(0.0, f64::max);
    println!("MRT Couette, {n} cells: largest error {err:.2e} of the wall speed");
    assert!(err < 2e-3, "{err}");
}

/// The decaying Taylor–Green vortex on [0, 2π]² (fully periodic): the amplitude's relative error
/// after its half-life against e^(−2νt), at a fixed τ, falls like Δ².
fn taylor_green(n: usize) -> f64 {
    let (nu, u0, l) = (0.1, 0.1, 2.0 * std::f64::consts::PI);
    let dx = l / n as f64;
    let mut d = domain([0.0; 3], [l, l, dx], dx);
    for f in fieldcraft_model::BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let dt = 0.3 * dx * dx / (3.0 * nu);
    let half_life = std::f64::consts::LN_2 / (2.0 * nu);
    let (setup, mut flow) = prepare(&model(1.0, nu, d, collide(half_life, dt, None, MRT)), &[]);
    let grid = *flow.grid();
    let mode = |cell: usize| {
        let [x, y, _] = grid.centre(cell);
        [-x.cos() * y.sin(), x.sin() * y.cos()]
    };
    let to_lattice = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|cell| {
        let [x, y, _] = grid.centre(cell);
        let pressure = -u0 * u0 / 4.0 * ((2.0 * x).cos() + (2.0 * y).cos());
        let m = mode(cell);
        (setup.lattice_density(pressure), [u0 * m[0] * to_lattice, u0 * m[1] * to_lattice, 0.0])
    });
    run_all(&setup, &mut flow);
    let (mut num, mut den) = (0.0, 0.0);
    for cell in 0..grid.len() {
        let (_, u) = flow.moments(cell).unwrap();
        let m = mode(cell);
        num += u[0] * m[0] + u[1] * m[1];
        den += m[0] * m[0] + m[1] * m[1];
    }
    let amplitude = num / den * setup.velocity_scale();
    (amplitude / (u0 * (-2.0 * nu * setup.steps as f64 * setup.time_step).exp()) - 1.0).abs()
}

#[test]
fn taylor_green_decays_at_the_viscous_rate_to_second_order() {
    let errors: Vec<f64> = [16, 32, 64].iter().map(|&n| taylor_green(n)).collect();
    let p = orders(&errors, 2.0);
    println!("MRT Taylor–Green, 16 → 64 cells: amplitude errors {errors:?}, observed orders {p:.3?}");
    // Measured 2.2e-2, 5.5e-3, 1.4e-3 (orders 2.015, 2.004).
    assert!(errors[2] < 2e-3, "{errors:?}");
    assert!(p.iter().all(|o| *o > 1.8), "second order: {p:.3?}");
}

/// Couette flow with Smagorinsky's model through the stresses' rate: the strain rate is uniform,
/// so the eddy viscosity is the closed form (C Δx)² U / h and the wall shear ρ (ν + ν_t) U / h
/// (as `turbulence.rs` checks for the other collisions).
#[test]
fn smagorinsky_couette_flow_takes_the_closed_form_eddy_viscosity() {
    let (h, nu, rho, u_wall, n, c) = (1.0, 1e-3, 2.0, 0.1, 16, 0.4);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
    let mut p = collide(5.0 * h * h / nu, dx, None, MRT);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: c, prandtl: 0.9, van_driest: None, wall_model: None });
    }
    let (setup, mut flow) = prepare(&model(rho, nu, d, p), &[]);
    run_all(&setup, &mut flow);
    let gamma = u_wall / h;
    let nu_t = (c * dx) * (c * dx) * gamma;
    let tau = rho * (nu + nu_t) * gamma * dx * dx;
    let g = |name: &str| flow.groups().iter().position(|g| g.name == name).unwrap();
    let bottom = flow.forces()[g("y_min")][0] * setup.force_scale();
    let ratio = flow.eddy_viscosity_ratio();
    println!("MRT Smagorinsky Couette: wall shear {bottom:.6e} vs {tau:.6e}; ν_t/ν {ratio:.5} vs {:.5}", nu_t / nu);
    assert!((bottom - tau).abs() < 2e-3 * tau, "{bottom} vs {tau}");
    assert!((ratio - nu_t / nu).abs() < 2e-2 * nu_t / nu, "{ratio}");
}

/// The same result bit for bit on one and on several threads.
#[test]
fn results_are_bit_for_bit_the_same_on_one_and_many_threads() {
    let l = 48.0;
    let mut d = domain([0.0; 3], [l, l, 1.0], 1.0);
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let m = model(1.0, 0.1 * l / 1000.0, d, collide(300.0, 1.0, None, MRT));
    let solve =
        |threads: usize| fieldcraft_parallel::with_threads(threads, || fieldcraft_lbm::solve_fluid(&m, None, &[], &Limits::default()).unwrap());
    let one = solve(1);
    assert_eq!(one, solve(7));
    assert_eq!(one, solve(3));
}

/// A periodic 8³ box at rest at the lattice viscosity `nu` (Δx = Δt = 1) with Smagorinsky's model, which is what
/// lets a fluid's own τ⁺ go below 0.505.
fn box_at_rest(nu: f64, collision: Collision) -> Result<Setup, fieldcraft_lbm::Error> {
    let mut d = domain([0.0; 3], [8.0; 3], 1.0);
    for f in fieldcraft_model::BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let mut p = collide(10.0, 1.0, None, collision);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None });
    }
    Setup::new(&model(1.0, nu, d, p), None, &[], &Limits::default())
}

/// Below `mrt::MIN_VISCOSITY` the collision is linearly unstable at rest (`mrt_stability.rs`), so the step is
/// refused whatever the subgrid model: at Han, Ooka & Kikumoto's A20 viscosity (5e-5) and one per cent below the
/// bound. From the bound up it is set up, with a warning until the bound at lattice Mach 0.1. The other
/// collisions keep the floor of 0.5001.
#[test]
fn below_its_stability_bound_the_collision_is_refused() {
    use fieldcraft_lbm::mrt::{MIN_VISCOSITY, MIN_VISCOSITY_AT_MACH_0_1};
    for nu in [5.0e-5, 0.99 * MIN_VISCOSITY] {
        let err = match box_at_rest(nu, MRT) {
            Ok(_) => panic!("ν {nu}: set up"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("stability bound of the MRT collision") && err.contains("0.50354") && err.contains("regularized"), "{err}");
    }
    let warned = |s: &Setup| s.warnings.iter().any(|w| w.contains("the MRT collision is stable at rest"));
    assert!(warned(&box_at_rest(MIN_VISCOSITY, MRT).unwrap()));
    assert!(!warned(&box_at_rest(MIN_VISCOSITY_AT_MACH_0_1, MRT).unwrap()));
    for collision in [Collision::Trt, Collision::Regularized, Collision::Hrr] {
        let setup = box_at_rest(5.0e-5, collision).unwrap();
        assert!(!warned(&setup), "{collision:?}");
    }
}
