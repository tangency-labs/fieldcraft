//! Heat in the fluid: the temperature lattice against closed-form solutions and a published
//! benchmark (de Vahl Davis, heavy).

mod common;

use common::*;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::{Limits, run};
use fieldcraft_model::{BoxFace, FaceCondition, FluidThermal, Outflow, ThermalCondition};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// Plates at y = 0 and y = 1, periodic in x and z, no flow: held at 0 and 1, the steady
/// temperature is linear, and anti-bounce-back puts the walls exactly halfway between cells.
#[test]
fn conduction_between_plates_is_linear() {
    for n in [8usize, 16] {
        let (h, alpha) = (1.0, 0.05);
        let dx = h / n as f64;
        let mut d = domain([0.0; 3], [dx, h, dx], dx);
        d.faces.x_min = PERIODIC;
        d.faces.x_max = PERIODIC;
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal::default();
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(0.0) };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        d.thermal = Some(t);
        // τ_T = 1: Δt = Δx² / (8 α)
        let dt = dx * dx / (8.0 * alpha);
        let m = hot_model(1.0, 0.1, alpha, d, step(5.0 * h * h / alpha, dt, None));
        let (setup, mut flow) = prepare(&m, &[]);
        assert!((setup.thermal.as_ref().unwrap().tau - 1.0).abs() < 1e-12);
        run_all(&setup, &mut flow);
        let mut err: f64 = 0.0;
        for j in 0..n {
            let y = (j as f64 + 0.5) * dx;
            let t = flow.temperature(j).unwrap();
            err = err.max((t - y / h).abs());
        }
        println!("conduction between plates, {n} cells: max error {err:.2e}");
        assert!(err < 1e-9, "{n}: {err}");
    }
}

/// A flux q through the cold plate and a held temperature on the other: the steady drop over
/// the gap is q H / k, and the heat the hot plate receives is the flux times its area.
#[test]
fn a_heat_flux_sets_the_gradient_and_the_wall_heat_balances() {
    let (n, h, alpha, rho) = (16usize, 1.0, 0.05, 2.0);
    let dx = h / n as f64;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let q = 0.3;
    let mut t = FluidThermal::default();
    t.faces.y_min = ThermalCondition::Flux { value: q };
    t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
    d.thermal = Some(t);
    let dt = dx * dx / (8.0 * alpha);
    let m = hot_model(rho, 0.1, alpha, d, step(8.0 * h * h / alpha, dt, None));
    let setup = fieldcraft_lbm::Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = fieldcraft_lbm::Flow::new(&setup, &[]).unwrap();
    let r = run(&setup, &mut flow, &Limits::default()).unwrap();
    let k = rho * alpha;
    // the temperature halfway below the first cell (the wall) against the hot wall at 1
    let t0 = f64::from(r.frames.last().unwrap().temperature[0]) + 0.5 * dx * q / k;
    let drop = t0 - 1.0;
    let expected = q * h / k;
    println!("flux plate: drop {drop:.6} vs q H / k = {expected:.6}");
    assert!((drop - expected).abs() < 1e-6 * expected, "{drop} vs {expected}");
    let hot = r.heat.iter().find(|w| w.name == "y_max").unwrap();
    let last = *hot.power.last().unwrap();
    let area = dx * dx;
    println!("heat into the hot wall: {last:.6e} vs q A = {:.6e}", q * area);
    assert!((last - q * area).abs() < 1e-6 * q * area, "{last} vs {}", q * area);
    let cold = r.heat.iter().find(|w| w.name == "y_min").unwrap();
    assert!((cold.power.last().unwrap() + q * area).abs() < 1e-6 * q * area, "the flux face gives q A");
}

/// A channel with a velocity inlet at 1, slip (adiabatic) walls and a pressure outlet: the
/// plug flow carries the inlet temperature through the whole channel. The flow starts at
/// its steady state (an impulsive start in a frictionless channel rings for a long time).
#[test]
fn the_flow_carries_the_inlet_temperature() {
    let (h, l, umax, nu, alpha) = (1.0, 4.0, 0.1, 0.05, 0.01);
    let dx = h / 8.0;
    let mut d = domain([0.0; 3], [l, h, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: fieldcraft_model::Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.y_min = FaceCondition::Slip {};
    d.faces.y_max = FaceCondition::Slip {};
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let mut t = FluidThermal::default();
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    d.thermal = Some(t);
    let lattice_u = 0.05;
    let m = hot_model(1.0, nu, alpha, d, step(4.0 * l / umax, lattice_u * dx / umax, None));
    let (setup, mut flow) = prepare(&m, &[]);
    flow.set_equilibrium(|_| (1.0, [lattice_u, 0.0, 0.0]));
    run_all(&setup, &mut flow);
    let grid = *flow.grid();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for c in 0..grid.len() {
        let t = flow.temperature(c).unwrap();
        lo = lo.min(t);
        hi = hi.max(t);
    }
    println!("inlet temperature carried through: {lo:.6}..{hi:.6}");
    assert!(lo > 1.0 - 1e-3 && hi < 1.0 + 1e-3, "{lo}..{hi}");
    let f = flow.fields().unwrap();
    assert_eq!(f.temperature.len(), grid.len());
}

/// Natural convection in a square cavity (de Vahl Davis, Int. J. Numer. Meth. Fluids 3 (1983)
/// 249–264): hot left wall at 1, cold right wall at 0, adiabatic top and bottom, Pr = 0.71.
/// The mean Nusselt number on the hot wall at Ra = 10³ is 1.118 (10⁴: 2.243); `u_max` is the
/// fastest flow in the cavity in units of α / H (3.65 and 19.6). Returns the Nusselt number
/// and the steady power into the hot wall and into the cold wall (the mean over the last
/// tenth of the run).
fn cavity(ra: f64, n: usize, u_max: f64, collision: fieldcraft_model::Collision) -> (f64, f64, f64) {
    let (h, pr) = (1.0, 0.71);
    let dx = h / n as f64;
    let alpha = 1.0;
    let nu = pr * alpha;
    // g β ΔT H³ / (ν α) = Ra, with ΔT = 1, H = 1
    let g_beta = ra * nu * alpha / (h * h * h);
    let mut d = domain([0.0; 3], [h, h, dx], dx);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let mut t = FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -g_beta, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    d.thermal = Some(t);
    // the fastest flow at lattice speed 0.05, and a lattice diffusivity of at most 1/4 (τ_T ≤ 1.5)
    let dt = (0.05 * dx / (u_max * alpha / h)).min(dx * dx / (4.0 * alpha));
    let m = hot_model(1.0, nu, alpha, d, collide(4.0 * h * h / alpha, dt, None, collision));
    let setup = fieldcraft_lbm::Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = fieldcraft_lbm::Flow::new(&setup, &[]).unwrap();
    let r = run(&setup, &mut flow, &Limits::default()).unwrap();
    // Nu = (heat into the fluid through the hot wall) H / (k ΔT A): the hot wall's history
    // holds the power the wall receives (negative when it heats the fluid).
    let steady = |name: &str| {
        let w = r.heat.iter().find(|w| w.name == name).unwrap();
        let m = w.power.len();
        let tail = &w.power[m - m / 10..];
        tail.iter().sum::<f64>() / tail.len() as f64
    };
    let (hot, cold) = (steady("x_min"), steady("x_max"));
    let k = alpha;
    let area = h * dx;
    let nu_mean = -hot * h / (k * 1.0 * area);
    println!("{collision:?} cavity Ra = {ra:.0e}, {n} cells: Nu = {nu_mean:.4}, hot wall {hot:.6e}, cold wall {cold:.6e}; {:?}", r.warnings);
    (nu_mean, hot, cold)
}

/// At steady state the heat the hot wall gives the fluid is the heat the cold wall takes.
fn heat_balances(hot: f64, cold: f64, tolerance: f64) {
    assert!(hot < 0.0 && cold > 0.0, "the hot wall heats the fluid ({hot}) and the cold wall cools it ({cold})");
    assert!((hot + cold).abs() < tolerance * cold, "heat in {hot:e} vs heat out {cold:e}");
}

/// The published cavity on a coarse grid, every run: Ra = 10³ on 16 cells is within 2 % of
/// de Vahl Davis (the heavy test refines it), and the heat balances between the walls.
#[test]
fn natural_convection_cavity_is_close_on_a_coarse_grid() {
    // Both collisions: buoyancy enters through Guo's force term, which each re-derives.
    for collision in COLLISIONS {
        let (nu3, hot, cold) = cavity(1e3, 16, 3.65, collision);
        assert!((nu3 / 1.118 - 1.0).abs() < 0.02, "{collision:?} Ra 1e3: {nu3}");
        heat_balances(hot, cold, 1e-3);
    }
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn natural_convection_cavity_matches_de_vahl_davis() {
    for collision in COLLISIONS {
        let (nu3, hot, cold) = cavity(1e3, 32, 3.65, collision);
        assert!((nu3 / 1.118 - 1.0).abs() < 0.01, "{collision:?} Ra 1e3: {nu3}");
        heat_balances(hot, cold, 1e-3);
        let (nu4, hot, cold) = cavity(1e4, 48, 19.6, collision);
        assert!((nu4 / 2.243 - 1.0).abs() < 0.01, "{collision:?} Ra 1e4: {nu4}");
        heat_balances(hot, cold, 1e-3);
    }
}

/// A periodic box with a uniform flow and a blob of heat: the total temperature is conserved
/// to round-off and the blob moves with the flow.
#[test]
fn a_periodic_flow_conserves_and_carries_heat() {
    let (l, u0, nu, alpha) = (4.0, 0.1, 0.05, 0.01);
    let dx = l / 32.0;
    let mut d = domain([0.0; 3], [l, dx, dx], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.initial_velocity = [u0, 0.0, 0.0];
    d.thermal = Some(FluidThermal::default());
    let m = hot_model(1.0, nu, alpha, d, step(l / u0, 0.05 * dx / u0, None));
    let (setup, mut flow) = prepare(&m, &[]);
    let grid = *flow.grid();
    flow.set_temperature(|c| if grid.coords(c)[0] < 8 { 1.0 } else { 0.0 });
    let total = |f: &fieldcraft_lbm::Flow| (0..grid.len()).map(|c| f.temperature(c).unwrap()).sum::<f64>();
    let t0 = total(&flow);
    run_all(&setup, &mut flow);
    let t1 = total(&flow);
    let row: Vec<String> = (0..32).map(|i| format!("{:.3}", flow.temperature(grid.index([i, 0, 0])).unwrap())).collect();
    println!("periodic blob after one lap: total {t0:.6} -> {t1:.6}; {}", row.join(" "));
    assert!((t1 - t0).abs() < 1e-10 * t0, "{t0} vs {t1}");
}

/// A sine of temperature in a periodic box drifting with a uniform flow: it keeps its shape,
/// moves with the flow and decays as e^(−α k² t). At rest the decay converges to the
/// diffusivity at second order. With the flow the D3Q7 lattice's linear equilibrium diffuses
/// less along the flow by (τ_T − ½) u² in lattice units, 4 u² of the diffusivity (Chopard,
/// Falcone and Latt, Eur. Phys. J. Special Topics 171 (2009) 245–249): the measured decay
/// shows that deficit, which bounds the lattice speed a heat problem should run at.
#[test]
fn a_drifting_sine_decays_at_the_diffusivity() {
    let (l, nu, alpha, amplitude) = (4.0, 0.05, 0.01, 0.25);
    let k = 2.0 * std::f64::consts::PI / l;
    let decay_time = 1.0 / (alpha * k * k);
    let mut at_rest = Vec::new();
    for (n, lattice_u) in [(32usize, 0.0), (64, 0.0), (64, 0.1), (64, 0.15)] {
        let dx = l / n as f64;
        // At rest τ_T = 1; with the flow at `lattice_u` the time step follows the speed.
        let u0 = if lattice_u > 0.0 { 0.1 } else { 0.0 };
        let dt = if lattice_u > 0.0 { lattice_u * dx / u0 } else { dx * dx / (8.0 * alpha) };
        let mut d = domain([0.0; 3], [l, dx, dx], dx);
        for f in BoxFace::ALL {
            *d.faces.get_mut(f) = PERIODIC;
        }
        d.initial_velocity = [u0, 0.0, 0.0];
        d.thermal = Some(FluidThermal { initial_temperature: kelvin(0.5), ..FluidThermal::default() });
        let m = hot_model(1.0, nu, alpha, d, step(decay_time, dt, None));
        let (setup, mut flow) = prepare(&m, &[]);
        let tau = setup.thermal.as_ref().unwrap().tau;
        let grid = *flow.grid();
        let x = |c: usize| (grid.coords(c)[0] as f64 + 0.5) * dx;
        flow.set_temperature(|c| 0.5 + amplitude * (k * x(c)).sin());
        run_all(&setup, &mut flow);
        let t = setup.steps as f64 * setup.time_step;
        // The sine's amplitude and phase from the Fourier coefficient at k, then the shape.
        let (mut a, mut b) = (0.0, 0.0);
        for i in 0..n {
            let c = grid.index([i, 0, 0]);
            let temp = flow.temperature(c).unwrap();
            a += temp * (k * x(c)).sin();
            b += temp * (k * x(c)).cos();
        }
        let measured = (a * a + b * b).sqrt() * 2.0 / n as f64;
        let drift = -b.atan2(a) / k;
        let mut shape: f64 = 0.0;
        for i in 0..n {
            let c = grid.index([i, 0, 0]);
            shape = shape.max((flow.temperature(c).unwrap() - 0.5 - measured * (k * (x(c) - drift)).sin()).abs());
        }
        let alpha_measured = -(measured / amplitude).ln() / (k * k * t);
        let error = alpha_measured / alpha - 1.0;
        let lag = (drift - u0 * t).rem_euclid(l);
        let lag = lag.min(l - lag);
        println!(
            "drifting sine, {n} cells, lattice speed {lattice_u}: τ_T = {tau:.4}, α measured {alpha_measured:.5} vs {alpha} ({error:+.2e}), lag {lag:.2e}, shape error {shape:.2e}"
        );
        // A Fourier mode of a linear scheme on a periodic grid stays one, to round-off.
        assert!(shape < 1e-9 * amplitude, "{n} cells at {lattice_u}: the sine lost its shape by {shape}");
        assert!(lag < dx / 5.0, "{n} cells at {lattice_u}: drift {drift} vs {}", u0 * t);
        if lattice_u > 0.0 {
            let deficit = -4.0 * lattice_u * lattice_u;
            assert!((error - deficit).abs() < 0.01, "{n} cells at {lattice_u}: α error {error} vs the deficit {deficit}");
        } else {
            assert!(error.abs() < 1e-3, "{n} cells at rest: α error {error}");
            at_rest.push(error.abs());
        }
    }
    let order = orders(&at_rest, 2.0);
    println!("decay at rest converges at order {order:?}");
    assert!(order.iter().all(|o| *o > 1.8), "{order:?}");
}
