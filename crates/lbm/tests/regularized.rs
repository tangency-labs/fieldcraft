//! The recursive regularised collision (`collision: "regularized"`) against the analytic cases of
//! the verification suite, its conservation and determinism, and the stability it adds over the
//! two-relaxation-time scheme as τ⁺ approaches ½.
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test regularized -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Error, Flow, Limits, Setup, solve_fluid};
use fieldcraft_model::{BoxFace, Collision, FaceCondition, Outflow, Procedure, Profile};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const REGULARIZED: Collision = Collision::Regularized;

fn sci(v: &[f64]) -> String {
    v.iter().map(|x| format!("{x:.3e}")).collect::<Vec<_>>().join(", ")
}

/// A fluid step of `steps` time steps of 1 (lattice units).
fn steps(steps: u64, collision: Collision) -> Procedure {
    Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(steps),
        time_step: Some(1.0),
        frames: 1,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    }
}

/// Plane Poiseuille flow between walls at y = 0 and 1 driven by a body force, `n` cells across,
/// relaxation time `tau`; the largest error of the velocity over its peak.
fn plane_poiseuille(n: usize, tau: f64) -> f64 {
    let (h, umax) = (1.0, 0.05);
    let dx = h / n as f64;
    // ν Δt / Δx² = (τ − ½)/3 with Δt = Δx²: ν = (τ − ½)/3.
    let nu = (tau - 0.5) / 3.0;
    let a = 8.0 * nu * umax / (h * h);
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    let m = model(1.0, nu, d, collide(6.0 * h * h / nu, dx * dx, None, REGULARIZED));
    let (setup, mut flow) = prepare(&m, &[]);
    assert!((setup.tau_plus - tau).abs() < 1e-12);
    run_all(&setup, &mut flow);
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let exact = a / (2.0 * nu) * y * (h - y);
        let (_, u) = flow.moments(j).unwrap();
        err = err.max((u[0] * setup.velocity_scale() - exact).abs() / umax);
        // round-off only, over up to 10⁵ steps
        assert!(u[1].abs() < 1e-12 && u[2].abs() < 1e-12, "{:e} {:e}", u[1], u[2]);
    }
    err
}

#[test]
fn plane_poiseuille_converges_at_second_order_with_the_body_force() {
    // The forcing is Guo's, re-derived for the regularised update: a wrong force term shows as a
    // first-order error or a profile that does not converge.
    for tau in [1.0, 0.6] {
        let errors: Vec<f64> = [8, 16, 32].iter().map(|&n| plane_poiseuille(n, tau)).collect();
        let p = orders(&errors, 2.0);
        println!("regularized Poiseuille, τ = {tau}, 8 → 32 cells: errors {}, orders {p:.3?}", sci(&errors));
        assert!(errors[2] < 2e-3, "{errors:?}");
        for o in p {
            assert!((o - 2.0).abs() < 0.1, "second order: {o}");
        }
    }
}

#[test]
fn couette_flow_is_linear_and_the_wall_shear_is_mu_u_over_h() {
    let (h, nu, rho, u_wall, n) = (1.0, 1.0 / 6.0, 2.0, 0.1, 16);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
    let m = model(rho, nu, d, collide(5.0 * h * h / nu, dx * dx, None, REGULARIZED));
    let (setup, mut flow) = prepare(&m, &[]);
    run_all(&setup, &mut flow);
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let (_, u) = flow.moments(j).unwrap();
        err = err.max((u[0] * setup.velocity_scale() - u_wall * y / h).abs() / u_wall);
    }
    let tau = rho * nu * u_wall / h * dx * dx;
    let g = |name: &str| flow.groups().iter().position(|g| g.name == name).unwrap();
    let fs = setup.force_scale();
    let (bottom, top) = (flow.forces()[g("y_min")][0] * fs, flow.forces()[g("y_max")][0] * fs);
    println!("regularized Couette, {n} cells: max error {err:.2e}; wall shear {bottom:.6e} and {top:.6e}, expected ±{tau:.6e}");
    assert!(err < 1e-10, "{err}");
    assert!((bottom - tau).abs() < 1e-6 * tau, "{bottom} vs {tau}");
    assert!((top + tau).abs() < 1e-6 * tau, "{top} vs {tau}");
}

/// The decaying Taylor–Green vortex (as in the verification suite) with the regularised
/// collision: the amplitude's relative error after its half-life, and the relative mass change.
fn taylor_green(n: usize) -> (f64, f64) {
    let (nu, u0, l) = (0.1, 0.1, 2.0 * std::f64::consts::PI);
    let dx = l / n as f64;
    let mut d = domain([0.0; 3], [l, l, dx], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let dt = 0.3 * dx * dx / (3.0 * nu);
    let half_life = std::f64::consts::LN_2 / (2.0 * nu);
    let m = model(1.0, nu, d, collide(half_life, dt, None, REGULARIZED));
    let (setup, mut flow) = prepare(&m, &[]);
    let grid = *flow.grid();
    let mode = |c: usize| {
        let [x, y, _] = grid.centre(c);
        [-x.cos() * y.sin(), x.sin() * y.cos()]
    };
    let to_lattice = 1.0 / setup.velocity_scale();
    flow.set_equilibrium(|c| {
        let [x, y, _] = grid.centre(c);
        let p = -u0 * u0 / 4.0 * ((2.0 * x).cos() + (2.0 * y).cos());
        let m = mode(c);
        (setup.lattice_density(p), [u0 * m[0] * to_lattice, u0 * m[1] * to_lattice, 0.0])
    });
    let mass0 = flow.mass();
    run_all(&setup, &mut flow);
    let (mut num, mut den) = (0.0, 0.0);
    for c in 0..grid.len() {
        let (_, u) = flow.moments(c).unwrap();
        let m = mode(c);
        num += u[0] * m[0] + u[1] * m[1];
        den += m[0] * m[0] + m[1] * m[1];
    }
    let amplitude = num / den * setup.velocity_scale();
    let exact = u0 * (-2.0 * nu * setup.steps as f64 * setup.time_step).exp();
    ((amplitude / exact - 1.0).abs(), (flow.mass() - mass0).abs() / mass0)
}

#[test]
fn taylor_green_vortex_decays_at_the_viscous_rate_to_second_order_and_conserves_mass() {
    let runs: Vec<(f64, f64)> = [16, 32, 64].iter().map(|&n| taylor_green(n)).collect();
    let errors: Vec<f64> = runs.iter().map(|r| r.0).collect();
    let p = orders(&errors, 2.0);
    println!("regularized Taylor–Green, 16 → 64 cells: errors {}, orders {p:.3?}", sci(&errors));
    assert!(errors[2] < 2e-3, "{errors:?}");
    for o in p {
        assert!(o > 1.8, "second order: {o}");
    }
    for (_, drift) in runs {
        assert!(drift < 1e-13, "mass is conserved to round-off: {drift}");
    }
}

/// A sphere in a periodic box driven by a body force: bounce-back walls and the regularised
/// collision conserve the mass to round-off, and at steady state the wall takes the driving force.
#[test]
fn bounce_back_walls_conserve_mass_and_balance_the_driving_force() {
    let (l, a, steps_run) = (12usize, 1e-6, 1500);
    let side = l as f64;
    let mut d = domain([0.0; 3], [side; 3], 1.0);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    d.obstacles = vec!["Ball".into()];
    let ball = obstacle("Ball", sphere([side / 2.0; 3], 2.5, 3), &["skin"]);
    let m = model(1.0, 5.0 / 6.0, d, steps(steps_run, REGULARIZED));
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&ball));
    let mass0 = flow.mass();
    run_all(&setup, &mut flow);
    let drift = (flow.mass() - mass0).abs() / mass0;
    let skin = flow.groups().iter().position(|g| g.name == "skin").unwrap();
    let (force, drive) = (flow.forces()[skin][0], a * flow.fluid_cells() as f64);
    println!("regularized sphere array: mass drift {drift:.2e}; wall force {force:.6e}, driving force {drive:.6e}");
    assert!(drift < steps_run as f64 * 1e-14, "{drift}");
    assert!((force - drive).abs() < 1e-3 * drive, "{force} vs {drive}");
}

/// A channel with a sphere, a parabolic inlet, a pressure outlet and a slip face, large enough
/// for several parallel blocks.
fn ball_channel(collision: Collision) -> (fieldcraft_model::Model, fieldcraft_lbm::Obstacle) {
    let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], 0.05);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.02, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_max = FaceCondition::Slip {};
    d.obstacles = vec!["Ball".into()];
    d.acceleration = [0.0, 0.0, -0.01];
    let ball = obstacle("Ball", sphere([0.7, 0.5, 0.45], 0.2, 3), &["skin"]);
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(0.5),
        steps: None,
        time_step: None,
        frames: 2,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    };
    (model(1.0, 0.01, d, p), ball)
}

#[test]
fn results_are_bit_for_bit_the_same_on_one_and_many_threads() {
    let (m, ball) = ball_channel(REGULARIZED);
    let solve = |threads: usize| {
        fieldcraft_parallel::with_threads(threads, || solve_fluid(&m, None, std::slice::from_ref(&ball), &Limits::default()).unwrap())
    };
    let one = solve(1);
    assert!(one.grid.len() > 3 * 4096, "several parallel blocks: {}", one.grid.len());
    assert_eq!(one, solve(7));
    assert_eq!(one, solve(3));
    let skin = one.forces.iter().find(|f| f.name == "skin").unwrap();
    assert!(skin.force.last().unwrap()[0] > 0.0, "the flow drags the ball downstream");
    // The two collisions differ (the step really uses the one it names).
    let (trt, _) = ball_channel(Collision::Trt);
    assert_ne!(one.frames, solve_fluid(&trt, None, std::slice::from_ref(&ball), &Limits::default()).unwrap().frames);
}

/// The lid-driven cavity of `n` cells at lid speed 0.1 (lattice) and Reynolds number `re`, run
/// `steps`: `None` when it survives, else the time step it went unstable at.
fn cavity_failure(n: usize, re: f64, steps_run: u64, collision: Collision) -> Option<u64> {
    let l = n as f64;
    let mut d = domain([0.0; 3], [l, l, 1.0], 1.0);
    d.faces.y_max = FaceCondition::MovingWall { velocity: [0.1, 0.0, 0.0] };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let m = model(1.0, 0.1 * l / re, d, steps(steps_run, collision));
    let (_, mut flow) = prepare(&m, &[]);
    (0..steps_run).find(|_| match flow.step() {
        Ok(()) => false,
        Err(Error::Unstable(_)) => true,
        Err(e) => panic!("{e}"),
    })
}

#[test]
fn a_lid_driven_cavity_at_tau_0_51_stays_stable_where_the_two_relaxation_time_scheme_fails() {
    // 48² cells, Re 1500: τ⁺ = 0.5096, a cell Reynolds number of 31 at the lid.
    let (n, re, run) = (48, 1500.0, 1500);
    let trt = cavity_failure(n, re, run, Collision::Trt);
    let regularized = cavity_failure(n, re, run, REGULARIZED);
    println!("cavity {n}², Re {re}: two relaxation times fail at {trt:?}, regularized at {regularized:?} (of {run} steps)");
    assert!(trt.is_some(), "the two-relaxation-time scheme fails here");
    assert_eq!(regularized, None);
}

/// A cylinder of `cells` diameter in the Schäfer–Turek channel (2.2 × 0.41, parabolic inlet with
/// mean 1, pressure outlet, quasi-2D), starting at the mean inlet velocity, at Reynolds number
/// `re` and inlet peak lattice speed 0.1, run `steps`: the time step it failed at, if it did.
fn cylinder_failure(cells: f64, re: f64, steps_run: u64, collision: Collision) -> Option<u64> {
    let d = 0.1;
    let dx = d / cells;
    let mut dom = domain([0.0; 3], [2.2, 0.41, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [1.5, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    dom.faces.z_min = PERIODIC;
    dom.faces.z_max = PERIODIC;
    dom.initial_velocity = [1.0, 0.0, 0.0];
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, 512), &["side", "bottom", "top"]);
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(steps_run),
        time_step: Some(0.1 * dx / 1.5),
        frames: 1,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    };
    let m = model(1.0, d / re, dom, p);
    let setup = Setup::new(&m, None, std::slice::from_ref(&cyl), &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, std::slice::from_ref(&cyl)).unwrap();
    (0..steps_run).find(|_| match flow.step() {
        Ok(()) => false,
        Err(Error::Unstable(_)) => true,
        Err(e) => panic!("{e}"),
    })
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_coarse_cylinder_wake_stays_stable_where_the_two_relaxation_time_scheme_fails() {
    // Ten cells across the cylinder: Re 200 is τ⁺ = 0.510, cell Reynolds number 15 at the peak
    // inlet speed; the inlet and outlet reflect the start-up's pressure waves.
    let run = 6000;
    let trt = cylinder_failure(10.0, 200.0, run, Collision::Trt);
    let regularized = cylinder_failure(10.0, 200.0, run, REGULARIZED);
    println!("cylinder, D = 10 cells, Re 200: two relaxation times fail at {trt:?}, regularized at {regularized:?} (of {run} steps)");
    assert!(trt.is_some(), "the two-relaxation-time scheme fails here");
    assert_eq!(regularized, None);
}

/// The inlet/outlet channel of the verification suite: a parabolic inlet, a pressure outlet,
/// `n` cells across. Returns the largest velocity error over the inlet peak at mid-length and in
/// the cell next to the outlet, the pressure gradient over Poiseuille's, and the relative mass
/// change over the last quarter of the run (inflow and outflow balance once it is steady).
fn channel(n: usize, collision: Collision) -> (f64, f64, f64, f64) {
    let (h, l, nu, umax) = (1.0, 3.0, 0.05, 0.1);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let dt = 0.1 / 3.0 * dx * dx / nu;
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(4.0 * h * h / nu),
        steps: None,
        time_step: Some(dt),
        frames: 4,
        magic: None,
        collision: Some(collision),
        turbulence: None,
        steady: None,
    };
    let r = solve_fluid(&model(1.0, nu, d, p), None, &[], &Limits::default()).unwrap();
    let last = r.frames.last().unwrap();
    let profile_error = |i: usize| {
        (0..n)
            .map(|j| {
                let y = (j as f64 + 0.5) * dx;
                (f64::from(last.velocity[r.grid.index([i, j, 0])][0]) - 4.0 * umax * y * (h - y) / (h * h)).abs() / umax
            })
            .fold(0.0, f64::max)
    };
    let p = |i: usize| f64::from(last.pressure[r.grid.index([i, n / 2, 0])]);
    let (i0, i1) = (r.grid.cells[0] / 4, 3 * r.grid.cells[0] / 4);
    let gradient = (p(i1) - p(i0)) / ((i1 - i0) as f64 * dx);
    let expected = -12.0 * nu * (2.0 / 3.0 * umax) / (h * h);
    let quarter = &r.frames[r.frames.len() - 2..];
    let mass = (quarter[1].mass - quarter[0].mass).abs() / quarter[0].mass;
    (profile_error(r.grid.cells[0] / 2), profile_error(r.grid.cells[0] - 2), gradient / expected, mass)
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn channel_with_a_parabolic_inlet_and_a_pressure_outlet_develops_poiseuille_flow() {
    for collision in COLLISIONS {
        let (mid, outlet, gradient, mass) = channel(16, collision);
        println!(
            "{collision:?} inlet/outlet channel: profile error {mid:.2e} at mid-length, {outlet:.2e} next to the outlet; pressure gradient / Poiseuille {gradient:.5}; mass change over the last quarter {mass:.1e}"
        );
        assert!(mid < 5e-3, "{mid}");
        assert!((gradient - 1.0).abs() < 0.01, "{gradient}");
        assert!(mass < 1e-4, "{mass}");
        // The regularised outlet holds the equilibrium (see fieldcraft_lbm::flow), so the
        // stress is off near it: measured 3.2e-2 at this τ⁺ = 0.6, 7e-2 at 0.52.
        assert!(outlet < if collision == REGULARIZED { 0.04 } else { 5e-3 }, "{outlet}");
    }
}

/// A walled channel, 16 cells across and 48 long (quasi-2D), a parabolic inlet with peak lattice
/// speed 0.1, a pressure outlet, starting at that speed, τ⁺ = 0.52: whether it runs `steps` steps.
/// Under the regularised collision the outlet's populations rebuilt from Zou–He's bounce-back
/// non-equilibrium part grow an oscillation along the outlet that flips sign every step (without
/// the equilibrium outlet this run went unstable at step 301).
fn outlet_channel_failure(steps_run: u64, collision: Collision) -> Option<u64> {
    let mut d = domain([0.0; 3], [48.0, 16.0, 1.0], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.initial_velocity = [0.1, 0.0, 0.0];
    let (_, mut flow) = prepare(&model(1.0, 0.02 / 3.0, d, steps(steps_run, collision)), &[]);
    (0..steps_run).find(|_| match flow.step() {
        Ok(()) => false,
        Err(Error::Unstable(_)) => true,
        Err(e) => panic!("{e}"),
    })
}

#[test]
fn a_regularised_channel_with_a_pressure_outlet_stays_stable_at_tau_0_52() {
    assert_eq!(outlet_channel_failure(2000, REGULARIZED), None);
}

#[test]
fn a_step_names_its_collision_in_the_setup_and_the_result() {
    let mut d = domain([0.0; 3], [1.0, 1.0, 0.25], 0.25);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.01, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    let m = model(1.0, 0.1, d.clone(), collide(0.1, 0.01, None, REGULARIZED));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    assert_eq!(setup.collision, REGULARIZED);
    assert!(setup.warnings.is_empty(), "{:?}", setup.warnings);
    assert_eq!(solve_fluid(&m, None, &[], &Limits::default()).unwrap().stats.collision, REGULARIZED);
    // Λ means nothing to the regularised collision: setting it is pointed out.
    let m = model(1.0, 0.1, d.clone(), collide(0.1, 0.01, Some(0.25), REGULARIZED));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    assert!(setup.warnings.iter().any(|w| w.contains("Λ")), "{:?}", setup.warnings);
    // Below its measured limit the regularised collision is warned about too.
    let m = model(1.0, 0.0125, d.clone(), collide(0.1, 0.01, None, REGULARIZED)); // τ⁺ 0.506
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    assert!(setup.tau_plus < fieldcraft_lbm::setup::REGULARIZED_TAU_LIMIT, "{}", setup.tau_plus);
    assert!(setup.warnings.iter().any(|w| w.contains("inlet and outlet")), "{:?}", setup.warnings);
    // The default is the two-relaxation-time scheme.
    let mut default = collide(0.1, 0.01, None, Collision::Trt);
    if let Procedure::Fluid { collision, .. } = &mut default {
        *collision = None;
    }
    assert_eq!(Setup::new(&model(1.0, 0.1, d, default), None, &[], &Limits::default()).unwrap().collision, Collision::Trt);
}
