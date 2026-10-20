//! Verification of the lattice Boltzmann solver against analytic solutions: plane Poiseuille
//! (exact with Λ = 3/16, second order otherwise), Couette (exact, with the wall shear force),
//! Taylor–Green vortex decay (second order), Poiseuille pipe flow on a staircase wall, Stokes drag
//! of a periodic array of spheres (Hasimoto 1959), inlet/outlet channel flow, mass conservation
//! to round-off, and the stability limits.
//!
//! Each test prints its errors and measured orders (`cargo test -p fieldcraft-lbm -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Error, Limits, Setup, solve_fluid};
use fieldcraft_model::{FaceCondition, Profile};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

fn sci(v: &[f64]) -> String {
    v.iter().map(|x| format!("{x:.3e}")).collect::<Vec<_>>().join(", ")
}

/// Plane Poiseuille flow between walls at y = 0 and y = 1, driven by a body force; periodic in x
/// and z (one cell each). τ = 1; returns the largest error of the velocity over its peak.
fn plane_poiseuille(n: usize, magic: f64) -> f64 {
    let (h, nu, umax) = (1.0, 1.0 / 6.0, 0.1);
    let dx = h / n as f64;
    let a = 8.0 * nu * umax / (h * h);
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.acceleration = [a, 0.0, 0.0];
    let m = model(1.0, nu, d, step(4.0 * h * h / nu, dx * dx, Some(magic)));
    let (setup, mut flow) = prepare(&m, &[]);
    assert!((setup.tau_plus - 1.0).abs() < 1e-12);
    run_all(&setup, &mut flow);
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let exact = a / (2.0 * nu) * y * (h - y);
        let (_, u) = flow.moments(j).unwrap();
        err = err.max((u[0] * setup.velocity_scale() - exact).abs() / umax);
        assert!(u[1].abs() < 1e-14 && u[2].abs() < 1e-14);
    }
    err
}

#[test]
fn plane_poiseuille_is_exact_with_the_magic_parameter_and_second_order_otherwise() {
    // Λ = 3/16: halfway bounce-back is exact for the parabola (Ginzburg et al. 2008).
    for n in [4, 8, 16] {
        let e = plane_poiseuille(n, 3.0 / 16.0);
        println!("Poiseuille channel, Λ = 3/16, {n} cells: max error {e:.2e}");
        assert!(e < 1e-10, "{n}: {e}");
    }
    // Λ = 1/4 (BGK at τ = 1): the wall moves off the halfway point by O(Δx²).
    let errors: Vec<f64> = [8, 16, 32, 64].iter().map(|&n| plane_poiseuille(n, 0.25)).collect();
    let p = orders(&errors, 2.0);
    println!("Poiseuille channel, Λ = 1/4, 8 → 64 cells: errors {}, orders {p:.3?}", sci(&errors));
    assert!(errors[0] > 1e-6, "the error is real at Λ = 1/4");
    for o in p {
        assert!((o - 2.0).abs() < 0.05, "second order: {o}");
    }
}

#[test]
fn couette_flow_is_linear_and_the_wall_shear_is_mu_u_over_h() {
    let (h, nu, rho, u_wall, n) = (1.0, 1.0 / 6.0, 2.0, 0.1, 16);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
    let m = model(rho, nu, d, step(5.0 * h * h / nu, dx * dx, None));
    let (setup, mut flow) = prepare(&m, &[]);
    run_all(&setup, &mut flow);
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let (_, u) = flow.moments(j).unwrap();
        err = err.max((u[0] * setup.velocity_scale() - u_wall * y / h).abs() / u_wall);
    }
    // The shear stress μ U / H on each wall, over its area Δx · Δx.
    let tau = rho * nu * u_wall / h * dx * dx;
    let g = |name: &str| flow.groups().iter().position(|g| g.name == name).unwrap();
    let fs = setup.force_scale();
    let bottom = flow.forces()[g("y_min")][0] * fs;
    let top = flow.forces()[g("y_max")][0] * fs;
    println!("Couette, {n} cells: max error {err:.2e}; wall shear {bottom:.6e} and {top:.6e}, expected ±{tau:.6e}");
    assert!(err < 1e-10, "{err}");
    assert!((bottom - tau).abs() < 1e-6 * tau, "{bottom} vs {tau}");
    assert!((top + tau).abs() < 1e-6 * tau, "{top} vs {tau}");
}

/// The decaying Taylor–Green vortex on [0, 2π]², periodic: returns the relative error of the
/// velocity amplitude after its half-life, and the relative change of the total mass.
fn taylor_green(n: usize) -> (f64, f64) {
    let (nu, u0, l) = (0.1, 0.1, 2.0 * std::f64::consts::PI);
    let dx = l / n as f64;
    let mut d = domain([0.0; 3], [l, l, dx], dx);
    for f in fieldcraft_model::BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    // τ = 0.8: Δt = 0.3 Δx² / (3 ν).
    let dt = 0.3 * dx * dx / (3.0 * nu);
    let half_life = std::f64::consts::LN_2 / (2.0 * nu);
    let m = model(1.0, nu, d, step(half_life, dt, None));
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
    let t = setup.steps as f64 * setup.time_step;
    let exact = u0 * (-2.0 * nu * t).exp();
    ((amplitude / exact - 1.0).abs(), (flow.mass() - mass0).abs() / mass0)
}

#[test]
fn taylor_green_vortex_decays_at_the_viscous_rate_to_second_order_and_conserves_mass() {
    let runs: Vec<(f64, f64)> = [16, 32, 64].iter().map(|&n| taylor_green(n)).collect();
    let errors: Vec<f64> = runs.iter().map(|r| r.0).collect();
    let p = orders(&errors, 2.0);
    println!(
        "Taylor–Green, 16 → 64 cells: amplitude errors {}, orders {p:.3?}; mass drift {}",
        sci(&errors),
        sci(&runs.iter().map(|r| r.1).collect::<Vec<_>>())
    );
    assert!(errors[2] < 2e-3, "{errors:?}");
    for o in p {
        assert!(o > 1.8, "second order: {o}");
    }
    for (_, drift) in runs {
        assert!(drift < 1e-13, "mass is conserved to round-off: {drift}");
    }
}

/// Poiseuille flow in a pipe of radius 1 along x (a tube obstacle, staircase wall), driven by a
/// body force; `k` cells per radius. Returns the relative L2 error of the velocity over the pipe's
/// cross-section and the relative error of the flow rate.
fn pipe(k: usize) -> (f64, f64) {
    let (r, nu, umax) = (1.0, 1.0 / 6.0, 0.1);
    let dx = r / k as f64;
    let half = r + 2.0 * dx;
    let mut d = domain([0.0, -half, -half], [dx, half, half], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    let a = 4.0 * nu * umax / (r * r);
    d.acceleration = [a, 0.0, 0.0];
    d.obstacles = vec!["Tube".into()];
    let tube = obstacle("Tube", tube(-1.0, dx + 1.0, [0.0, 0.0], r, 3.0 * half, 512), &["inner", "outer", "start", "end"]);
    let m = model(1.0, nu, d, step(5.0 * r * r / nu, dx * dx, None));
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&tube));
    run_all(&setup, &mut flow);
    let grid = *flow.grid();
    let (mut e2, mut u2, mut q) = (0.0, 0.0, 0.0);
    for c in 0..grid.len() {
        let [_, y, z] = grid.centre(c);
        let rr = y * y + z * z;
        let exact = if rr < r * r { a / (4.0 * nu) * (r * r - rr) } else { 0.0 };
        let u = flow.moments(c).map_or(0.0, |(_, u)| u[0] * setup.velocity_scale());
        e2 += (u - exact) * (u - exact);
        u2 += exact * exact;
        q += u * dx * dx;
    }
    let q_exact = std::f64::consts::PI * a * r.powi(4) / (8.0 * nu);
    // Momentum balance: the wall takes the whole driving force.
    let wall: f64 = flow.groups().iter().zip(flow.forces()).filter(|(g, _)| g.body.is_some()).map(|(_, f)| f[0]).sum::<f64>() * setup.force_scale();
    let drive = a * flow.fluid_cells() as f64 * dx.powi(3);
    assert!((wall - drive).abs() < 1e-6 * drive, "{wall} vs {drive}");
    ((e2 / u2).sqrt(), (q / q_exact - 1.0).abs())
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn pipe_poiseuille_on_a_staircase_wall_converges() {
    let runs: Vec<(f64, f64)> = [4, 8, 16].iter().map(|&k| pipe(k)).collect();
    let l2: Vec<f64> = runs.iter().map(|r| r.0).collect();
    let flow: Vec<f64> = runs.iter().map(|r| r.1).collect();
    println!(
        "Poiseuille pipe, 4 → 16 cells per radius: L2 errors {} (orders {:.2?}), flow-rate errors {} (orders {:.2?})",
        sci(&l2),
        orders(&l2, 2.0),
        sci(&flow),
        orders(&flow, 2.0)
    );
    assert!(l2[2] < l2[0] && l2[2] < 0.05, "{l2:?}");
}

/// The Stokes drag coefficient K = F / (6 π μ R U) of a simple cubic array of spheres at volume
/// fraction c (H. Hasimoto, JFM 5(2):317–328, 1959, with the higher terms of A. S. Sangani,
/// A. Acrivos, Int. J. Multiphase Flow 8(4):343–360, 1982), U the mean velocity over the cell.
fn hasimoto(c: f64) -> f64 {
    1.0 / (1.0 - 1.7601 * c.cbrt() + c - 1.5593 * c * c + 3.9799 * c.powf(8.0 / 3.0) - 3.0734 * c.powf(10.0 / 3.0))
}

/// A sphere of radius `r` cells at the centre of a periodic box of `l` cells, the fluid driven by
/// a small body force (Re ≪ 1); returns the measured K over Hasimoto's, and the relative mass
/// drift.
fn sphere_drag(r: f64, l: usize, steps: u64) -> (f64, f64) {
    // τ = 3: with Λ = 3/16 the steady solution does not depend on τ, and a high viscosity
    // reaches it sooner (the mean flow relaxes over ρ V / (6 π μ R K) steps).
    let (nu, a) = (5.0 / 6.0, 1e-6);
    let side = l as f64;
    let mut d = domain([0.0; 3], [side; 3], 1.0);
    for f in fieldcraft_model::BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    d.obstacles = vec!["Ball".into()];
    let ball = obstacle("Ball", sphere([side / 2.0; 3], r, 4), &["skin"]);
    let m = model(
        1.0,
        nu,
        d,
        fieldcraft_model::Procedure::Fluid {
            domain: "Flow".into(),
            end_time: None,
            steps: Some(steps),
            time_step: Some(1.0),
            frames: 1,
            magic: None,
            collision: None,
            steady: None,
        },
    );
    let (setup, mut flow) = prepare(&m, std::slice::from_ref(&ball));
    let mass0 = flow.mass();
    run_all(&setup, &mut flow);
    let cells = (l * l * l) as f64;
    let u: f64 = (0..l * l * l).filter_map(|c| flow.moments(c)).map(|(_, u)| u[0]).sum::<f64>() / cells;
    let skin = flow.groups().iter().position(|g| g.name == "skin").unwrap();
    let force = flow.forces()[skin][0];
    let fluid = flow.fluid_cells() as f64;
    // Momentum balance at steady state.
    assert!((force - a * fluid).abs() < 1e-3 * a * fluid, "{force} vs {}", a * fluid);
    // The array's drag per sphere is the driving force over the whole cell (the body force acts
    // on the fluid only, so the measured force is that times the fluid fraction).
    let k = force * cells / fluid / (6.0 * std::f64::consts::PI * nu * r * u);
    let c = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3) / cells;
    let re = u * 2.0 * r / nu;
    println!("sphere array, R = {r}, L = {l}: c = {c:.4}, Re = {re:.1e}, K = {k:.4}, Hasimoto {:.4}", hasimoto(c));
    assert!(re < 0.1);
    (k / hasimoto(c), (flow.mass() - mass0).abs() / mass0)
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn stokes_drag_of_a_sphere_array_matches_hasimoto() {
    let coarse = sphere_drag(3.0, 16, 1500);
    let fine = sphere_drag(6.0, 32, 3000);
    println!("sphere drag over Hasimoto: R = 3: {:.4}, R = 6: {:.4}", coarse.0, fine.0);
    assert!((fine.0 - 1.0).abs() < (coarse.0 - 1.0).abs(), "refining approaches the analytic drag");
    assert!((fine.0 - 1.0).abs() < 0.05, "{}", fine.0);
    // Round-off: each step's rounding moves the sum by a few ulps at most.
    assert!(coarse.1 < 1500.0 * 1e-14 && fine.1 < 3000.0 * 1e-14, "mass is conserved to round-off with bounce-back walls: {} {}", coarse.1, fine.1);
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn channel_with_a_parabolic_inlet_and_a_pressure_outlet_develops_poiseuille_flow() {
    // Plates 1 apart, 6 long; inlet peak 0.1; outlet at gauge pressure 0.
    let (h, l, nu, umax, n) = (1.0, 3.0, 0.05, 0.1, 16);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    // τ = 0.6. The pressure drop changes the density along the channel by Δp / (ρ c_s²) in lattice
    // units (the weakly compressible error), here 0.07%; the automatic time step would give 1.6%.
    let dt = 0.1 / 3.0 * dx * dx / nu;
    let procedure = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(4.0 * h * h / nu),
        steps: None,
        time_step: Some(dt),
        frames: 4,
        magic: None,
        collision: None,
        steady: None,
    };
    let m = model(1.0, nu, d, procedure);
    let r = solve_fluid(&m, Some("Run"), &[], &Limits::default()).unwrap();
    println!("inlet/outlet channel: {:?}", r.stats);
    assert_eq!(r.frames.len(), 5);
    assert!((r.stats.tau_plus - 0.6).abs() < 1e-12, "{:?}", r.stats);
    let last = r.frames.last().unwrap();
    // The profile at mid-length, against the parabola with the inlet's peak.
    let i = r.grid.cells[0] / 2;
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let exact = 4.0 * umax * y * (h - y) / (h * h);
        let c = r.grid.index([i, j, 0]);
        err = err.max((f64::from(last.velocity[c][0]) - exact).abs() / umax);
    }
    // The pressure gradient of plane Poiseuille flow, −12 μ Ū / H² with Ū = 2/3 of the peak.
    let p = |i: usize| f64::from(last.pressure[r.grid.index([i, n / 2, 0])]);
    let (i0, i1) = (r.grid.cells[0] / 4, 3 * r.grid.cells[0] / 4);
    let gradient = (p(i1) - p(i0)) / ((i1 - i0) as f64 * dx);
    let expected = -12.0 * nu * (2.0 / 3.0 * umax) / (h * h);
    println!("inlet/outlet channel: profile error {err:.2e}, pressure gradient {gradient:.5} vs {expected:.5}");
    assert!(err < 5e-3, "{err}");
    assert!((gradient / expected - 1.0).abs() < 0.01, "{gradient} vs {expected}");
    // The walls take the pressure drop: their drag balances it.
    let wall: f64 = r.forces.iter().filter(|f| f.name == "y_min" || f.name == "y_max").map(|f| f.force.last().unwrap()[0]).sum();
    assert!(wall > 0.0);
}

#[test]
fn stability_limits_are_enforced_with_clear_errors() {
    let channel = |u: f64, dt: f64, nu: f64| {
        let mut d = domain([0.0; 3], [4.0, 1.0, 0.25], 0.25);
        d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
        model(1.0, nu, d, step(10.0 * dt, dt, None))
    };
    let limits = Limits::default();
    // Mach: u Δt / Δx = 0.4 → Ma = 0.69.
    let e = Setup::new(&channel(1.0, 0.1, 1.0), None, &[], &limits).unwrap_err();
    assert!(e.to_string().contains("Mach") && e.to_string().contains("lower the time step"), "{e}");
    // τ: ν Δt / Δx² = 1e-4 → τ = 0.5003.
    let e = Setup::new(&channel(0.01, 0.1, 6.25e-5), None, &[], &limits).unwrap_err();
    assert!(e.to_string().contains("τ") && e.to_string().contains("refine the spacing"), "{e}");
    // Caps.
    let tight = Limits { max_cells: 10, ..limits };
    assert!(matches!(Setup::new(&channel(0.01, 0.01, 0.1), None, &[], &tight), Err(Error::Limit(_))));
    let tight = Limits { max_steps: 5, ..limits };
    assert!(matches!(Setup::new(&channel(0.01, 0.01, 0.1), None, &[], &tight), Err(Error::Limit(_))));
    let tight = Limits { max_cell_updates: 100, ..limits };
    assert!(matches!(Setup::new(&channel(0.01, 0.01, 0.1), None, &[], &tight), Err(Error::Limit(_))));
    let tight = Limits { max_frame_bytes: 100, ..limits };
    assert!(matches!(solve_fluid(&channel(0.01, 0.01, 0.1), None, &[], &tight), Err(Error::Limit(_))));
    // A flow pushed beyond the stable speed stops with an error, never a panic or NaNs.
    let (setup, mut flow) = prepare(&channel(0.01, 0.01, 0.1), &[]);
    flow.set_equilibrium(|_| (1.0, [0.6, 0.0, 0.0]));
    assert!(matches!(flow.step(), Err(Error::Unstable(_))));
    let _ = setup;
    // An obstacle the domain lists but the caller did not prepare.
    let mut m = channel(0.01, 0.01, 0.1);
    let mut d = m.fluid_domain("Flow").unwrap().clone();
    d.obstacles = vec!["Ghost".into()];
    m.replace_fluid_domain("Flow", d).unwrap();
    assert!(matches!(Setup::new(&m, None, &[], &limits), Err(Error::NotFound(_))));
}

#[test]
fn results_are_bit_for_bit_reproducible() {
    let mut d = domain([0.0; 3], [2.0, 1.0, 1.0], 0.125);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.1, 0.02, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0 };
    d.faces.z_max = FaceCondition::Slip {};
    d.obstacles = vec!["Ball".into()];
    let ball = obstacle("Ball", sphere([0.7, 0.5, 0.45], 0.2, 3), &["skin"]);
    let m = model(
        1.0,
        0.01,
        d,
        fieldcraft_model::Procedure::Fluid {
            domain: "Flow".into(),
            end_time: Some(1.0),
            steps: None,
            time_step: None,
            frames: 3,
            magic: None,
            collision: None,
            steady: None,
        },
    );
    let a = solve_fluid(&m, None, std::slice::from_ref(&ball), &Limits::default()).unwrap();
    let b = solve_fluid(&m, None, std::slice::from_ref(&ball), &Limits::default()).unwrap();
    assert_eq!(a, b);
    let skin = a.forces.iter().find(|f| f.name == "skin").unwrap();
    assert_eq!(skin.body.as_deref(), Some("Ball"));
    assert!(skin.force.last().unwrap()[0] > 0.0, "the flow drags the ball downstream");
    assert!(a.solid.iter().any(|s| *s));
    let probe = a.probe(3, [1.8, 0.5, 0.5]).unwrap();
    assert!(probe.1[0] > 0.0);
}
