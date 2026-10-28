//! The hybrid recursive regularised collision (`collision: "hrr"`, Jacob, Malaspinas and Sagaut
//! 2018) against the analytic laminar cases (plane Poiseuille and Couette flow, Smagorinsky's
//! Couette flow, the decaying Taylor–Green vortex), its determinism, σ = 1 against the recursive
//! regularised collision bit for bit, the refusals of σ, and the decay of a shear wave against the
//! exact linear analysis of the scheme and the paper's hyperviscosity (its eq. 47).
//!
//! Each test prints its numbers (`cargo test -p fieldcraft-lbm --test hrr -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_model::{BoxFace, Collision, FaceCondition, HrrSigma, Model, Procedure, Step, Turbulence};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};
const HRR: Collision = Collision::Hrr;

/// [`collide`] with the hybrid collision and σ (`None`: the default).
fn hrr(end_time: f64, time_step: f64, sigma: Option<f64>) -> Procedure {
    let mut p = collide(end_time, time_step, None, HRR);
    if let Procedure::Fluid { sigma: s, .. } = &mut p {
        *s = sigma.map(HrrSigma::Value);
    }
    p
}

/// Plane Poiseuille flow (walls at y = 0 and 1, driven by a body force, peak 0.1), `n` cells
/// across, diffusive time steps (τ fixed at 0.8), σ 0.99: the largest error of the velocity
/// against the parabola, over its peak.
fn poiseuille(n: usize) -> f64 {
    let (h, nu, umax) = (1.0, 0.1, 0.1);
    let dx = h / n as f64;
    let a = 8.0 * nu * umax / (h * h);
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [a, 0.0, 0.0];
    let (setup, mut flow) = prepare(&model(1.0, nu, d, hrr(6.0 * h * h / nu, dx * dx, None)), &[]);
    assert_eq!(setup.sigma, 0.99, "the default σ");
    run_all(&setup, &mut flow);
    (0..n)
        .map(|j| {
            let y = (j as f64 + 0.5) * dx;
            (flow.moments(j).unwrap().1[0] * setup.velocity_scale() - a / (2.0 * nu) * y * (h - y)).abs() / umax
        })
        .fold(0.0, f64::max)
}

/// The parabola's second derivative is constant, so the centred differences are exact inside and
/// only the wall cells' one-sided differences (weighted 1 − σ = 1 %) add to the regularised
/// collision's second-order wall error: the error still falls like Δ².
#[test]
fn poiseuille_flow_converges_to_the_parabola_like_dx_squared() {
    let errors: Vec<f64> = [8, 16, 32].iter().map(|&n| poiseuille(n)).collect();
    let p = orders(&errors, 2.0);
    println!("HRR Poiseuille, 8 → 32 cells: errors {errors:?}, observed orders {p:.3?}");
    assert!(errors[2] < 2e-3, "{errors:?}");
    assert!(p.iter().all(|o| (o - 2.0).abs() < 0.2), "second order: {p:.3?}");
}

/// Plane Couette flow (a moving wall at y = 1, 16 cells): the linear profile, whose strain rate
/// the differences give exactly.
#[test]
fn couette_flow_is_linear() {
    let (h, nu, u_wall, n) = (1.0, 0.05, 0.1, 16);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
    let (setup, mut flow) = prepare(&model(1.0, nu, d, hrr(6.0 * h * h / nu, dx, Some(0.9))), &[]);
    run_all(&setup, &mut flow);
    let err = (0..n)
        .map(|j| (flow.moments(j).unwrap().1[0] * setup.velocity_scale() - u_wall * (j as f64 + 0.5) * dx / h).abs() / u_wall)
        .fold(0.0, f64::max);
    println!("HRR Couette (σ 0.9), {n} cells: largest error {err:.2e} of the wall speed");
    assert!(err < 2e-3, "{err}");
}

/// Couette flow with Smagorinsky's model: the strain rate is uniform, so the eddy viscosity is the
/// closed form (C Δx)² U / h and the wall shear ρ (ν + ν_t) U / h (as `mrt.rs` checks).
#[test]
fn smagorinsky_couette_flow_takes_the_closed_form_eddy_viscosity() {
    let (h, nu, rho, u_wall, n, c) = (1.0, 1e-3, 2.0, 0.1, 16, 0.4);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.faces.y_max = FaceCondition::MovingWall { velocity: [u_wall, 0.0, 0.0] };
    let mut p = hrr(5.0 * h * h / nu, dx, Some(0.9));
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
    println!("HRR Smagorinsky Couette: wall shear {bottom:.6e} vs {tau:.6e}; ν_t/ν {ratio:.5} vs {:.5}", nu_t / nu);
    assert!((bottom - tau).abs() < 2e-3 * tau, "{bottom} vs {tau}");
    assert!((ratio - nu_t / nu).abs() < 2e-2 * nu_t / nu, "{ratio}");
}

/// The decaying Taylor–Green vortex on [0, 2π]² (fully periodic): the amplitude's relative error
/// after its half-life against e^(−2νt), at a fixed τ.
fn taylor_green(n: usize) -> f64 {
    let (nu, u0, l) = (0.1, 0.1, 2.0 * std::f64::consts::PI);
    let dx = l / n as f64;
    let mut d = domain([0.0; 3], [l, l, dx], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let dt = 0.3 * dx * dx / (3.0 * nu);
    let half_life = std::f64::consts::LN_2 / (2.0 * nu);
    let (setup, mut flow) = prepare(&model(1.0, nu, d, hrr(half_life, dt, None)), &[]);
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
    let mass = flow.mass();
    run_all(&setup, &mut flow);
    assert!((flow.mass() / mass - 1.0).abs() < 1e-12, "mass");
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

/// σ 0.99 adds a hyperviscosity of order (1 − σ)/σ (k Δx)² ν, which falls like Δ² too.
#[test]
fn taylor_green_decays_at_the_viscous_rate_to_second_order() {
    let errors: Vec<f64> = [16, 32, 64].iter().map(|&n| taylor_green(n)).collect();
    let p = orders(&errors, 2.0);
    println!("HRR Taylor–Green, 16 → 64 cells: amplitude errors {errors:?}, observed orders {p:.3?}");
    assert!(errors[2] < 2e-3, "{errors:?}");
    assert!(p.iter().all(|o| *o > 1.8), "second order: {p:.3?}");
}

/// σ = 1 leaves the populations' flux untouched: the recursive regularised collision, bit for
/// bit, in every cell's density and velocity and every wall force after many steps (with and
/// without Smagorinsky's model).
#[test]
fn sigma_one_is_the_regularised_collision_bit_for_bit() {
    let smagorinsky = Turbulence::Smagorinsky { constant: 0.17, prandtl: 0.9, van_driest: None, wall_model: None };
    for turbulence in [None, Some(smagorinsky)] {
        let run = |collision: Collision, sigma: Option<f64>| {
            let (m, ball) = ball_channel(collision, sigma.map(HrrSigma::Value), turbulence);
            let (setup, mut flow) = prepare(&m, std::slice::from_ref(&ball));
            for _ in 0..200 {
                flow.step().unwrap();
            }
            (setup, flow)
        };
        let (_, rr) = run(Collision::Regularized, None);
        let (setup, one) = run(HRR, Some(1.0));
        assert_eq!(setup.sigma, 1.0);
        let mut cells = 0;
        for cell in 0..rr.grid().len() {
            match (rr.moments(cell), one.moments(cell)) {
                (Some((ra, ua)), Some((rb, ub))) => {
                    assert_eq!(ra.to_bits(), rb.to_bits(), "density of cell {cell}");
                    for k in 0..3 {
                        assert_eq!(ua[k].to_bits(), ub[k].to_bits(), "velocity of cell {cell}");
                    }
                    cells += 1;
                }
                (None, None) => {}
                _ => panic!("cell {cell} is fluid in one run only"),
            }
        }
        for (a, b) in rr.forces().iter().zip(one.forces()) {
            assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
        }
        assert!(rr.mass().to_bits() == one.mass().to_bits());
        // And σ < 1 does change the flow.
        let (_, blended) = run(HRR, Some(0.9));
        let moved = (0..rr.grid().len()).any(|c| rr.moments(c).map(|m| m.1) != blended.moments(c).map(|m| m.1));
        println!("σ = 1 against regularized ({turbulence:?}): {cells} fluid cells equal bit for bit after 200 steps; σ 0.9 differs: {moved}");
        assert!(moved);
    }
}

/// The same result bit for bit on one and on several threads.
#[test]
fn results_are_bit_for_bit_the_same_on_one_and_many_threads() {
    let (m, ball) = ball_channel(HRR, Some(HrrSigma::Value(0.95)), None);
    let solve = |threads: usize| {
        fieldcraft_parallel::with_threads(threads, || fieldcraft_lbm::solve_fluid(&m, None, std::slice::from_ref(&ball), &Limits::default()).unwrap())
    };
    let one = solve(1);
    assert!(one.grid.len() > 3 * 4096, "several parallel blocks: {}", one.grid.len());
    assert_eq!(one, solve(7));
    assert_eq!(one, solve(3));
    assert_eq!(one.stats.collision, HRR);
}

/// σ out of (0, 1], or not a number, is refused by the model with a message that names σ; σ = 0
/// by name. "dynamic" without a subgrid model is refused by the solver, naming `turbulence`. σ with
/// another collision is pointed out. Nothing panics.
#[test]
fn bad_sigmas_are_refused_and_dynamic_needs_a_subgrid_model() {
    let d = || {
        let mut d = domain([0.0; 3], [1.0, 1.0, 0.25], 0.25);
        d.faces.y_max = FaceCondition::MovingWall { velocity: [0.01, 0.0, 0.0] };
        d
    };
    let base = || {
        let mut m = Model::default();
        m.add_material(
            "Fluid",
            fieldcraft_model::Material {
                density: Some(1.0),
                fluid: Some(fieldcraft_model::Fluid::Newtonian { viscosity: 0.1 }),
                ..fieldcraft_model::Material::default()
            },
        )
        .unwrap();
        m.add_fluid_domain("Flow", d()).unwrap();
        m
    };
    let with = |sigma: HrrSigma| {
        let mut p = hrr(0.1, 0.01, None);
        if let Procedure::Fluid { sigma: s, .. } = &mut p {
            *s = Some(sigma);
        }
        p
    };
    for bad in [0.0, -0.0, 1.5, -0.2, f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1.0 + f64::EPSILON] {
        let e = base().add_step(Step::new("Run", with(HrrSigma::Value(bad))), None).unwrap_err().to_string();
        println!("σ = {bad}: {e}");
        assert!(e.contains("σ"), "{e}");
        if bad == 0.0 {
            assert!(e.contains("cannot be 0") && e.contains("diverges"), "{e}");
        }
    }
    let mut m = base();
    m.add_step(Step::new("Run", with(HrrSigma::Dynamic)), None).unwrap();
    let e = Setup::new(&m, None, &[], &Limits::default()).unwrap_err().to_string();
    println!("σ = dynamic: {e}");
    assert!(e.contains("dynamic σ needs a subgrid model") && e.contains("`turbulence`"), "{e}");
    assert!(fieldcraft_lbm::solve_fluid(&m, None, &[], &Limits::default()).is_err());
    // The edges of the range are taken.
    for good in [1.0, 1e-3, 0.5] {
        let mut m = base();
        m.add_step(Step::new("Run", with(HrrSigma::Value(good))), None).unwrap();
        assert_eq!(Setup::new(&m, None, &[], &Limits::default()).unwrap().sigma, good);
    }
    // σ with another collision means nothing, and is pointed out.
    let mut p = collide(0.1, 0.01, None, Collision::Regularized);
    if let Procedure::Fluid { sigma, .. } = &mut p {
        *sigma = Some(HrrSigma::Value(0.9));
    }
    let setup = Setup::new(&model(1.0, 0.1, d(), p), None, &[], &Limits::default()).unwrap();
    assert!(setup.warnings.iter().any(|w| w.contains("σ")), "{:?}", setup.warnings);
    let setup = Setup::new(&model(1.0, 0.1, d(), hrr(0.1, 0.01, None)), None, &[], &Limits::default()).unwrap();
    assert!(setup.warnings.is_empty(), "{:?}", setup.warnings);
}

/// The leading-order rate of Jacob et al.'s eq. 47: ν k² [1 + ((1 − σ)/σ)(1 − sin k / k)],
/// ν = (τ̄ − ½)/3.
fn eq47_rate(k: f64, tau: f64, sigma: f64) -> f64 {
    let nu = (tau - 0.5) / 3.0;
    nu * k * k * (1.0 + (1.0 - sigma) / sigma * (1.0 - k.sin() / k))
}

struct Wave {
    n: usize,
    tau: f64,
    sigma: f64,
    measured: f64,
    exact: f64,
    eq47: f64,
}

/// The shear wave u_x = U sin(k y), k = 2π/N, U = 1e-4 lattice speeds, in a periodic box of
/// 2 × N × 2 cells (ρ = 1, no force), set up in model units so the time step makes τ̄ the given
/// value (spacing 1, time step 1, ν = (τ̄ − ½)/3). The amplitude is the projection of u_x on
/// sin(k y) over every cell; the rate is −ln(A(t₁)/A(t₀))/(t₁ − t₀) per step, t₀ = 200 and
/// t₁ = t₀ + round(2 / predicted rate), at most 60 000 steps.
fn shear_wave(n: usize, tau: f64, sigma: f64) -> Wave {
    let (dx, dt, amplitude) = (1.0, 1.0, 1e-4);
    let nu = (tau - 0.5) / 3.0 * dx * dx / dt;
    let mut d = domain([0.0; 3], [2.0 * dx, n as f64 * dx, 2.0 * dx], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    let (setup, mut flow) = prepare(&model(1.0, nu, d, hrr(1.0, dt, Some(sigma))), &[]);
    assert!((setup.tau_plus - tau).abs() < 1e-12, "{}", setup.tau_plus);
    assert_eq!(flow.grid().cells, [2, n, 2]);
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let grid = *flow.grid();
    let mode = |cell: usize| (k * grid.coords(cell)[1] as f64).sin();
    flow.set_equilibrium(|cell| (1.0, [amplitude * mode(cell), 0.0, 0.0]));
    let measure = |flow: &Flow| {
        let (mut num, mut den) = (0.0, 0.0);
        for cell in 0..grid.len() {
            let m = mode(cell);
            num += flow.moments(cell).unwrap().1[0] * m;
            den += m * m;
        }
        num / den
    };
    let exact = exact_rate(k, tau, sigma);
    let t0 = 200u64;
    let t1 = t0 + ((2.0 / exact).round() as u64).clamp(1, 60_000 - t0);
    for _ in 0..t0 {
        flow.step().unwrap();
    }
    let a0 = measure(&flow);
    for _ in t0..t1 {
        flow.step().unwrap();
    }
    let a1 = measure(&flow);
    let measured = -(a1 / a0).ln() / (t1 - t0) as f64;
    Wave { n, tau, sigma, measured, exact, eq47: eq47_rate(k, tau, sigma) }
}

/// The gate's tolerances: within 1 % of the exact rate always, and of eq. 47 at k Δx = 2π/64.
fn check(waves: &[Wave]) {
    println!("   N   kΔx     τ̄     σ      measured     exact −ln λ  eq. 47       measured/exact − 1  measured/eq47 − 1  pass");
    let mut failed = Vec::new();
    for w in waves {
        let k = 2.0 * std::f64::consts::PI / w.n as f64;
        let (r_exact, r_47) = (w.measured / w.exact, w.measured / w.eq47);
        let ok_exact = (r_exact - 1.0).abs() < 0.01;
        let ok_47 = w.n != 64 || (r_47 - 1.0).abs() < 0.01;
        let pass = ok_exact && ok_47;
        println!(
            "{:>4}  {k:.4}  {:.2}  {:.3}  {:.5e}  {:.5e}  {:.5e}  {:+.2e}           {:+.2e}{}          {}",
            w.n,
            w.tau,
            w.sigma,
            w.measured,
            w.exact,
            w.eq47,
            r_exact - 1.0,
            r_47 - 1.0,
            if w.n == 64 { " " } else { "*" },
            if pass { "ok" } else { "FAIL" }
        );
        if !pass {
            failed.push((w.n, w.tau, w.sigma));
        }
    }
    println!("(* eq. 47 is leading order in k Δx and not part of the gate there)");
    assert!(failed.is_empty(), "outside 1 %: {failed:?}");
}

/// The orientation values computed independently of this code: at τ̄ 0.8, k Δx π/4, σ 0.5 the
/// exact rate is 6.4112e-2 per step and eq. 47 gives 6.7834e-2; at τ̄ 0.55, k Δx 2π/64, σ 0.9
/// the exact rate is 1.6093e-4.
#[test]
fn the_shear_wave_analysis_gives_the_independent_values() {
    let pi = std::f64::consts::PI;
    assert!((exact_rate(pi / 4.0, 0.8, 0.5) / 6.4112e-2 - 1.0).abs() < 1e-4, "{}", exact_rate(pi / 4.0, 0.8, 0.5));
    assert!((eq47_rate(pi / 4.0, 0.8, 0.5) / 6.7834e-2 - 1.0).abs() < 1e-4, "{}", eq47_rate(pi / 4.0, 0.8, 0.5));
    assert!((exact_rate(2.0 * pi / 64.0, 0.55, 0.9) / 1.6093e-4 - 1.0).abs() < 1e-4, "{}", exact_rate(2.0 * pi / 64.0, 0.55, 0.9));
    // σ = 1: the regularised collision, whose small-k rate is ν k².
    let k = 2.0 * pi / 256.0;
    assert!((exact_rate(k, 0.8, 1.0) / (0.1 * k * k) - 1.0).abs() < 1e-3);
}

/// The gate's fast representative case: N = 16, τ̄ 0.8, σ 0.9.
#[test]
fn a_shear_wave_decays_at_the_exact_rate_of_the_linear_analysis() {
    check(&[shear_wave(16, 0.8, 0.9)]);
}

/// THE GATE: every N ∈ {64, 16, 8}, σ ∈ {0.5, 0.9, 0.97, 0.985}, τ̄ ∈ {0.8, 0.55}.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn shear_waves_decay_at_the_exact_rate_and_at_eq_47_when_resolved() {
    let mut waves = Vec::new();
    for n in [64, 16, 8] {
        for tau in [0.8, 0.55] {
            for sigma in [0.5, 0.9, 0.97, 0.985] {
                waves.push(shear_wave(n, tau, sigma));
            }
        }
    }
    check(&waves);
}
