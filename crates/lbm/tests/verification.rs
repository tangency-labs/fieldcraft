//! Verification of the lattice Boltzmann solver against analytic solutions: plane Poiseuille
//! (exact with Λ = 3/16, second order otherwise), Couette (exact, with the wall shear force),
//! Taylor–Green vortex decay (second order), Poiseuille pipe flow on a staircase wall, Stokes drag
//! of a periodic array of spheres (Hasimoto 1959), inlet/outlet and pressure-driven channel flow, backflow at a short outlet, the wake-at-outlet limitation, mass conservation
//! to round-off, and the stability limits.
//!
//! Each test prints its errors and measured orders (`cargo test -p fieldcraft-lbm -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Error, Limits, Setup, solve_fluid};
use fieldcraft_model::{FaceCondition, Outflow, Profile};

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
            turbulence: None,
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
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
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
        turbulence: None,
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

/// Plane Poiseuille flow driven by the pressures alone: Zou–He pressure faces at both ends
/// (Δp = 12 μ Ū L / H²), plates 1 apart, 3 long, 16 cells across. The profile is the parabola
/// within 3 % and each end holds its pressure: a guard for any later change of the open
/// boundaries (an anti-bounce-back outlet tried in #598 doubled the error and shifted both ends).
#[test]
fn pressure_driven_poiseuille_holds_its_end_pressures() {
    let (h, l, nu, umax, n) = (1.0, 3.0, 0.05, 0.1, 16);
    let dx = h / n as f64;
    let dp = 12.0 * nu * (2.0 / 3.0 * umax) * l / (h * h);
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Pressure { pressure: dp, outflow: Outflow::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let dt = 0.1 / 3.0 * dx * dx / nu;
    let procedure = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(4.0 * h * h / nu),
        steps: None,
        time_step: Some(dt),
        frames: 1,
        magic: None,
        steady: None,
        collision: None,
        turbulence: None,
    };
    let r = solve_fluid(&model(1.0, nu, d, procedure), Some("Run"), &[], &Limits::default()).unwrap();
    let last = r.frames.last().unwrap();
    let nx = r.grid.cells[0];
    let mut err: f64 = 0.0;
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let exact = 4.0 * umax * y * (h - y) / (h * h);
        err = err.max((f64::from(last.velocity[r.grid.index([nx / 2, j, 0])][0]) - exact).abs() / umax);
    }
    let p = |i: usize| f64::from(last.pressure[r.grid.index([i, n / 2, 0])]);
    println!("pressure-driven Poiseuille: profile error {err:.3e}, end pressures {:.5} / {:.5} of {dp:.5} / 0", p(0), p(nx - 1));
    assert!(err < 0.03, "{err}");
    assert!((p(0) / dp - 1.0).abs() < 1e-3 && (p(nx - 1) / dp).abs() < 1e-3, "{} {}", p(0), p(nx - 1));
}

/// Schäfer–Turek's channel at Re = 100, cut 1 diameter behind the cylinder (l = 0.3): the
/// shedding wake reaches the pressure outlet with reversed flow in some of its cells, and the run
/// stays stable and conserves mass (an anti-bounce-back outlet tried in #598 diverged here).
#[test]
fn a_short_channel_at_re_100_survives_backflow_at_the_outlet() {
    let (h, l, d, um, nu, dx) = (0.41, 0.3, 0.1, 1.5, 1e-3, 0.01);
    let mut dom = domain([0.0; 3], [l, h, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [um, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    dom.faces.z_min = PERIODIC;
    dom.faces.z_max = PERIODIC;
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, 256), &["side", "bottom", "top"]);
    let procedure = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(4.0),
        steps: None,
        time_step: Some(0.05 * dx / um),
        frames: 8,
        magic: None,
        steady: None,
        collision: None,
        turbulence: None,
    };
    let r = solve_fluid(&model(1.0, nu, dom, procedure), None, std::slice::from_ref(&cyl), &Limits::default()).unwrap();
    let [nx, ny, _] = r.grid.cells;
    let reversed = r.frames.iter().map(|f| (0..ny).filter(|j| f.velocity[r.grid.index([nx - 1, *j, 0])][0] < 0.0).count()).max().unwrap_or(0);
    let (m0, m1) = (r.frames[0].mass, r.frames.last().unwrap().mass);
    println!("Re 100 short channel: up to {reversed} reversed outlet cells, mass change {:.2e}", m1 / m0 - 1.0);
    // The case must actually reverse flow at the outlet, or the test passes vacuously.
    assert!(reversed > 0, "no reversed outlet cells: the wake no longer reaches the outlet");
    assert!((m1 / m0 - 1.0).abs() < 0.01);
}

/// Known limitation (#590): a Zou–He pressure outlet forces its layer to one pressure, so where a
/// wake has not recovered (Schäfer–Turek 2D-1 cut 1.5 diameters behind the cylinder) the last
/// layer departs from the profile of the layers inside it: the speed rises through that one-cell
/// pressure step. The reported extremes leave the open-boundary layers out for this reason. The
/// interior must stay smooth; the outlet layer's departure is printed, not checked.
#[test]
fn a_wake_at_a_pressure_outlet_distorts_only_the_last_layer() {
    let (h, l, d, um, nu, dx) = (0.41, 0.4, 0.1, 0.3, 1e-3, 0.01);
    let mut dom = domain([0.0; 3], [l, h, dx], dx);
    dom.faces.x_min = FaceCondition::Velocity { velocity: [um, 0.0, 0.0], profile: Profile::Parabolic };
    dom.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    dom.faces.z_min = PERIODIC;
    dom.faces.z_max = PERIODIC;
    dom.obstacles = vec!["Cylinder".into()];
    let cyl = obstacle("Cylinder", cylinder([0.2, 0.2], d / 2.0, -0.1, 0.1, 256), &["side", "bottom", "top"]);
    let procedure = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(4.0),
        steps: None,
        time_step: Some(0.025 * dx / um),
        frames: 1,
        magic: None,
        steady: None,
        collision: None,
        turbulence: None,
    };
    let r = solve_fluid(&model(1.0, nu, dom, procedure), None, std::slice::from_ref(&cyl), &Limits::default()).unwrap();
    let last = r.frames.last().unwrap();
    let [nx, ny, _] = r.grid.cells;
    let ux = |i: usize, j: usize| f64::from(last.velocity[r.grid.index([i, j, 0])][0]);
    let departure = |i: usize| (0..ny).map(|j| (ux(i, j) - (2.0 * ux(i - 1, j) - ux(i - 2, j))).abs()).fold(0.0, f64::max);
    println!("wake at the outlet: interior departure {:.2e}, outlet layer {:.2e} (of {um})", departure(nx - 2), departure(nx - 1));
    assert!(ux(nx - 1, ny / 2) < 0.5 * um, "the wake reaches the outlet");
    assert!(departure(nx - 2) < 0.03 * um, "{}", departure(nx - 2));
    assert!(r.in_open_layer(r.grid.index([nx - 1, 0, 0])) && !r.in_open_layer(r.grid.index([nx - 2, 0, 0])));
}

/// #650: the memory a run is refused or admitted on is the memory it allocates: the arrays the
/// CPU flow holds once it steps, one frame's fields and the frames, within 1% (measured +0.5%
/// for the channel, +0.2% with heat and buoyancy, −0.1% with a sphere, whose links the
/// estimate cannot know before voxelising). The frames and the fields are exact.
#[test]
fn the_memory_estimate_matches_what_a_run_allocates() {
    use fieldcraft_lbm::{Backend, Choice, Memory};
    let channel = |heat: bool, sphere: bool| {
        let mut d = domain([0.0; 3], [3.0, 1.5, 1.5], 0.0625);
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
        let mut obstacles = Vec::new();
        if sphere {
            d.obstacles = vec!["Ball".into()];
            obstacles.push(obstacle("Ball", common::sphere([1.0, 0.75, 0.75], 0.4, 3), &["surface"]));
        }
        let procedure = step(0.4, 0.01, None);
        let m = if heat {
            d.thermal = Some(fieldcraft_model::FluidThermal {
                initial_temperature: kelvin(20.0),
                gravity: [0.0, 0.0, -9.81],
                expansion: Some(1e-3),
                ..Default::default()
            });
            hot_model(1.0, 0.01, 0.01, d, procedure)
        } else {
            model(1.0, 0.01, d, procedure)
        };
        (m, obstacles)
    };
    for (heat, sphere) in [(false, false), (true, false), (false, true)] {
        let (m, obstacles) = channel(heat, sphere);
        let setup = Setup::new(&m, None, &obstacles, &Limits::default()).unwrap();
        let memory = Memory::of(&setup, Choice::Cpu).unwrap();
        // From the model alone (what a form shows) it is the same.
        assert_eq!(Memory::estimate(&m, None, Choice::Cpu).unwrap(), memory);
        let mut flow = fieldcraft_lbm::Flow::new(&setup, &obstacles).unwrap();
        flow.advance(1).unwrap();
        let n = setup.grid.len() as u64;
        // The flow and the result's solid flags.
        let held = flow.heap_bytes() + n;
        let fields = flow.fields().unwrap();
        let readback = (fields.density.capacity() * 8 + fields.velocity.capacity() * 24 + fields.temperature.capacity() * 8) as u64;
        let r = solve_fluid(&m, None, &obstacles, &Limits::default()).unwrap();
        let frames: u64 =
            r.frames.iter().map(|f| (f.velocity.capacity() * 12 + f.pressure.capacity() * 4 + f.temperature.capacity() * 4) as u64).sum();
        assert_eq!((memory.frames, memory.frame_bytes), (r.frames.len() as u64, frames), "heat {heat}");
        assert_eq!(memory.readback, readback, "heat {heat}");
        let measured = held + readback + frames;
        let over = memory.total() as f64 / measured as f64 - 1.0;
        assert!(memory.state >= held || sphere, "heat {heat}: {} below {held}", memory.state);
        assert!(over.abs() < 0.01, "heat {heat} sphere {sphere}: estimate {} vs {measured} allocated ({over:+.3})", memory.total());
    }
}

/// #650: the caps refuse with what to change: the frames that fit, the value to raise a cap
/// to and how far this machine allows it, or coarsening when even the machine's bound is
/// exceeded; a raised cap is checked against the bound before anything is allocated.
#[test]
fn caps_refuse_with_what_to_change_and_raised_caps_fit_the_machine() {
    use fieldcraft_lbm::{Choice, Memory};
    let mut d = domain([0.0; 3], [4.0, 1.0, 0.25], 0.25);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.01, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    let m = model(1.0, 0.1, d, step(0.1, 0.01, None));
    let machine = Limits { max_memory: 64_000_000_000, ..Limits::default() };
    let setup = Setup::new(&m, None, &[], &machine).unwrap();
    let memory = Memory::of(&setup, Choice::Cpu).unwrap();
    // 16 × 4 × 1 cells, 3 frames of 16 bytes per cell.
    assert_eq!((memory.cells, memory.frames, memory.frame_bytes), (64, 3, 3 * 64 * 16));

    let e = memory.check(&Limits { max_frame_bytes: 2100, ..machine }).unwrap_err().to_string();
    assert!(e.contains("2 frames fit: set the step's `frames` to at most 1, or raise `max_frame_bytes` to 3072 (this machine allows up to"), "{e}");
    let e = memory.check(&Limits { max_frame_bytes: 100, ..machine }).unwrap_err().to_string();
    assert!(e.contains("not even two frames fit: coarsen the spacing"), "{e}");
    // Without a bound (no machine known to the solver), no ceiling is named.
    let e = memory.check(&Limits { max_frame_bytes: 2100, ..Limits::default() }).unwrap_err().to_string();
    assert!(e.ends_with("or raise `max_frame_bytes` to 3072"), "{e}");

    // The cells: raise to what is needed, or coarsen when the machine cannot hold them.
    let e = Setup::new(&m, None, &[], &Limits { max_cells: 10, ..machine }).unwrap_err().to_string();
    assert!(e.contains("coarsen the spacing, or raise `max_cells` to 64 (this machine allows up to"), "{e}");
    let small = Limits { max_memory: 1_000_000_000, ..Limits::default() };
    assert_eq!(small.max_cells_ceiling(), Limits::MAX_CELLS, "never below the default");
    let e = small.check_cells([1000, 1000, 1000]).unwrap_err().to_string();
    assert!(e.ends_with("coarsen the spacing"), "{e}");
    assert!(machine.check_cells([1000, 400, 200]).unwrap_err().to_string().contains("raise `max_cells` to 80000000"));
    assert!(machine.max_cells_ceiling() > 80_000_000);
    let e = Setup::new(&m, None, &[], &Limits { max_cell_updates: 100, ..machine }).unwrap_err().to_string();
    assert!(e.contains("coarsen the spacing or shorten the step, or raise `max_cell_updates` to 640"), "{e}");

    // A raised cap admits the run while it fits the machine's bound, and refuses it, before
    // anything is allocated, once it does not.
    let raised = Limits { max_frame_bytes: Limits::MAX_FRAME_BYTES + 1, ..machine };
    assert!(raised.raised() && memory.check(&raised).is_ok());
    let tight = Limits { max_memory: memory.total() - 1, ..raised };
    let e = memory.check(&tight).unwrap_err().to_string();
    assert!(e.contains("above this machine's bound") && e.contains("frames fit: set the step's `frames` to at most"), "{e}");
    let e = solve_fluid(&m, None, &[], &tight).unwrap_err().to_string();
    assert!(e.contains("above this machine's bound"), "{e}");
    // Within the default caps the bound does not apply: such runs solve as they always have.
    assert!(memory.check(&Limits { max_memory: 1, ..Limits::default() }).is_ok());
}

#[test]
fn stability_limits_are_enforced_with_clear_errors() {
    let channel = |u: f64, dt: f64, nu: f64| {
        let mut d = domain([0.0; 3], [4.0, 1.0, 0.25], 0.25);
        d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
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
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
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
            turbulence: None,
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
