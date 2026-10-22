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
    for n in [8usize, 16, 32] {
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

/// A heated parallel-plate channel (the Graetz problem, Pe 200): a developed parabolic flow
/// enters at T_in = 0 between plates held at T_w = 1 that start at the inlet. The open boundary
/// sets every scalar population of the inlet and outlet layers after streaming, so the heat a
/// wall link gives a cell there never stays in the fluid; counted, a held wall meeting a held
/// inlet reported several times its heat (#675). The monitor must report exactly what the inner
/// layers received, which in the steady state they pass on: out toward the outlet as the flow's
/// enthalpy ρ c_p Q (T_bulk − T_in), and (a tenth of it at Pe 200, from the corners where the hot
/// walls meet the cold inlet) by conduction back into the inlet layer.
#[test]
fn wall_heat_from_plates_that_start_at_the_inlet_is_what_the_fluid_keeps() {
    for n in [8usize, 16] {
        let (h, l, umax, peclet) = (1.0, 4.0, 0.1, 200.0);
        let dx = h / n as f64;
        let alpha = 2.0 / 3.0 * umax * h / peclet;
        let mut d = domain([0.0; 3], [l, h, dx], dx);
        d.faces.x_min = FaceCondition::Velocity { velocity: [umax, 0.0, 0.0], profile: fieldcraft_model::Profile::Parabolic };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
        d.faces.z_min = PERIODIC;
        d.faces.z_max = PERIODIC;
        let mut t = FluidThermal { initial_temperature: kelvin(0.0), ..FluidThermal::default() };
        t.faces.x_min = ThermalCondition::Temperature { value: kelvin(0.0) };
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(1.0) };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        d.thermal = Some(t);
        let lattice_u = 0.05;
        // three flow-through times
        let m = hot_model(1.0, 0.05, alpha, d, step(3.0 * l / (2.0 / 3.0 * umax), lattice_u * dx / umax, None));
        let (setup, mut flow) = prepare(&m, &[]);
        let grid = *flow.grid();
        flow.set_equilibrium(|c| {
            let s = (grid.coords(c)[1] as f64 + 0.5) / n as f64;
            (1.0, [lattice_u * 4.0 * s * (1.0 - s), 0.0, 0.0])
        });
        run_all(&setup, &mut flow);
        // Lattice units: temperature times cell faces per step; the lattice runs in degrees
        // Celsius (T_w = 1, T_in = 0). The heat each wall face received is what it took.
        let heat = flow.heat();
        let wall = |name: &str| flow.groups().iter().position(|g| g.name == name).and_then(|k| heat.get(k).copied()).unwrap();
        let reported = -(wall("y_min") + wall("y_max"));
        // The heat the walls gave layer i in the last step, from its streamed populations: a
        // held wall's link brings back g = −g* + 2 w T_w for the g* that left, so the fluid
        // gained g − g* = 2 g − 2 w T_w (w = ⅛) through each.
        let layer = |i: usize| -> f64 {
            let from_wall = |j: usize, q: usize| 2.0 * flow.thermal_population(q, grid.index([i, j, 0])) - 0.25;
            from_wall(0, 3) + from_wall(n - 1, 4)
        };
        let last = grid.cells[0] - 2;
        let inner: f64 = (1..=last).map(layer).sum();
        // What streamed across the face between layers i and i + 1: +x into i + 1 minus −x into i
        // (inside: the open boundary overwrites what streamed into the open layers).
        let across = |i: usize| -> f64 {
            (0..n).map(|j| flow.thermal_population(1, grid.index([i + 1, j, 0])) - flow.thermal_population(2, grid.index([i, j, 0]))).sum()
        };
        // Steady: what the walls give layers 2..last is what crosses out of them.
        let kept: f64 = (2..last).map(layer).sum();
        let passed = across(last - 1) - across(1);
        // The enthalpy the flow carries toward the outlet, ρ c_p Q (T_bulk − T_in).
        let enthalpy: f64 = (0..n)
            .map(|j| {
                let c = grid.index([last - 1, j, 0]);
                flow.moments(c).unwrap().1[0] * flow.temperature(c).unwrap()
            })
            .sum();
        println!(
            "heated channel, {n} cells across: reported wall heat {reported:.6e}, given to the inner layers {inner:.6e}; \
             layers 2..: given {kept:.6e}, passed on {passed:.6e}; out {:.6e}, enthalpy rise {enthalpy:.6e}; reported / enthalpy {:.4}",
            across(last - 1),
            reported / enthalpy
        );
        assert!((reported / inner - 1.0).abs() < 1e-9, "{n}: the monitor reports {reported}, the inner layers got {inner}");
        // (three flow-through times leave 2e-3 of transient; twelve 1e-4)
        assert!((kept / passed - 1.0).abs() < 5e-3, "{n}: {kept} given, {passed} passed on");
        assert!((across(last - 1) / enthalpy - 1.0).abs() < 0.02, "{n}: {} out vs {enthalpy}", across(last - 1));
    }
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
/// diffusivity at second order: with Λ_T = 1/12 at τ_T = 1 the lattice decays a mode 0.105 (k Δx)²
/// faster than α k² (the exact amplification factor of the D3Q7 two-relaxation-time lattice:
/// 4.0e-3 on 32 cells per wavelength, 1.0e-3 on 64; a single relaxation time at τ_T = 1 had
/// −8e-4 and −2e-4, at τ_T = 1.5 −1.5e-2 and −3.8e-3). With the flow the D3Q7 lattice's linear equilibrium diffuses
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
        // The lattice's own phase speed is 1.5e-3 and 3.7e-3 fast at lattice speeds 0.1 and 0.15
        // (Λ_T = 1/12; the exact amplification factor), 0.10 and 0.24 cells over the run (a single
        // relaxation time at these τ_T: −0.9e-3 and −2.4e-3).
        assert!(lag < dx / 3.0, "{n} cells at {lattice_u}: drift {drift} vs {}", u0 * t);
        if lattice_u > 0.0 {
            let deficit = -4.0 * lattice_u * lattice_u;
            assert!((error - deficit).abs() < 0.01, "{n} cells at {lattice_u}: α error {error} vs the deficit {deficit}");
        } else {
            let kdx = k * dx;
            assert!(error.abs() < 0.12 * kdx * kdx, "{n} cells at rest: α error {error} vs 0.105 (k Δx)² = {:.2e}", 0.105 * kdx * kdx);
            at_rest.push(error.abs());
        }
    }
    let order = orders(&at_rest, 2.0);
    println!("decay at rest converges at order {order:?}");
    assert!(order.iter().all(|o| *o > 1.8), "{order:?}");
}

/// The de Vahl Davis cavity at Rayleigh number `ra` on `n` cells, run to 1.5 H²/α (steady to
/// well below the differences reported) at thermal relaxation time `tau` (`None`: the automatic
/// time step) with the thermal lattice's parameter Λ_T `magic(τ_T)`. Returns τ_T, the mean
/// Nusselt number of the hot wall (the heat its links take) and two profiles of the local one
/// along the wall: as `solve.probe`'s `wall_heat_flux` reads it, from the held wall temperature
/// and the two cell layers beside it, `(8 T_w − 9 T₁ + T₂) / (3 Δx)`, and the lattice's own heat
/// through each wall link, `2 w T_w − 2 g` per step (`g` the population the wall sent back),
/// which is what a conducting neighbour would receive.
fn cavity_wall_flux(ra: f64, n: usize, tau: Option<f64>, magic: &dyn Fn(f64) -> f64) -> (f64, f64, [Vec<f64>; 2]) {
    let (h, pr) = (1.0, 0.71);
    let dx = h / n as f64;
    let (alpha, nu) = (1.0, pr);
    let mut d = domain([0.0; 3], [h, h, dx], dx);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let mut t =
        FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -ra * nu * alpha, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    d.thermal = Some(t);
    let mut procedure = step(1.5 * h * h / alpha, 1.0, None);
    if let fieldcraft_model::Procedure::Fluid { time_step, .. } = &mut procedure {
        *time_step = tau.map(|tau| (tau - 0.5) * dx * dx / (4.0 * alpha));
    }
    let m = hot_model(1.0, nu, alpha, d, procedure);
    let mut setup = fieldcraft_lbm::Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let ts = setup.thermal.as_mut().unwrap();
    let tau = ts.tau;
    ts.magic = magic(tau);
    let mut flow = fieldcraft_lbm::Flow::new(&setup, &[]).unwrap();
    let r = run(&setup, &mut flow, &Limits::default()).unwrap();
    let hot = *r.heat.iter().find(|w| w.name == "x_min").unwrap().power.last().unwrap();
    let nu_mean = -hot * h / (alpha * h * dx);
    let grid = *flow.grid();
    let temp = |i: usize, j: usize| flow.temperature(grid.index([i, j, 0])).unwrap();
    let probe = (0..n).map(|j| (8.0 - 9.0 * temp(0, j) + temp(1, j)) / (3.0 * dx)).collect();
    // the heat per step through the link along +x of the wall cell, over the lattice
    // diffusivity: Nu = −heat / (α_lattice Δx) with H = ΔT = 1
    let alpha_lattice = (tau - 0.5) / 4.0;
    let link = (0..n).map(|j| -(2.0 * 0.125 - 2.0 * flow.thermal_population(1, grid.index([0, j, 0]))) / (alpha_lattice * dx)).collect();
    (tau, nu_mean, [probe, link])
}

/// The extreme of a profile (largest for `sign` 1, smallest for −1) on the parabola through the
/// extreme sample and its neighbours (`solve.probe`'s rule); at an end, the sample itself.
fn extreme(local: &[f64], sign: f64) -> f64 {
    let n = local.len();
    let k = (0..n).max_by(|a, b| (sign * local[*a]).total_cmp(&(sign * local[*b]))).unwrap();
    if k == 0 || k == n - 1 {
        return local[k];
    }
    let (a, b, c) = (local[k - 1], local[k], local[k + 1]);
    let curvature = a - 2.0 * b + c;
    if sign * curvature >= 0.0 {
        return b;
    }
    let d = (0.5 * (a - c) / curvature).clamp(-0.5, 0.5);
    b - 0.25 * (a - c) * d
}

/// The local wall heat flux of the de Vahl Davis cavity (Ra = 10³, 24 cells) at τ_T 0.8 and 1.4
/// (#629): with the thermal lattice's Λ_T fixed its largest and smallest values agree to the
/// numerical diffusion of the faster lattice flow at τ_T 1.4, while a single relaxation time
/// (Λ_T = (τ_T − ½)²) moves the wall and them by percents.
#[test]
fn local_wall_heat_flux_does_not_depend_on_the_thermal_relaxation_time() {
    let bgk = |tau: f64| (tau - 0.5) * (tau - 0.5);
    let change = |magic: &dyn Fn(f64) -> f64| {
        let [(_, _, a), (_, _, b)] = [Some(0.8), Some(1.4)].map(|tau| cavity_wall_flux(1e3, 24, tau, magic));
        [(&a[0], &b[0], 1.0), (&a[0], &b[0], -1.0), (&a[1], &b[1], 1.0), (&a[1], &b[1], -1.0)]
            .map(|(x, y, sign)| (extreme(y, sign) / extreme(x, sign) - 1.0).abs())
    };
    let (fixed, single) = (change(&fieldcraft_lbm::setup::thermal_magic), change(&bgk));
    println!(
        "cavity, 24 cells: relative change of Nu_max, Nu_min (probe, links) from τ_T 0.8 to 1.4: Λ_T fixed {fixed:.5?}, single relaxation time {single:.5?}"
    );
    assert!(fixed.iter().all(|d| *d < 3e-3), "{fixed:?}");
    assert!(single.iter().all(|d| *d > 2e-2), "{single:?}");
}

/// The wall-flux study behind `fieldcraft_lbm::setup::THERMAL_MAGIC` (#629): the local wall heat
/// flux of the de Vahl Davis cavity (de Vahl Davis 1983, Table V) against the thermal relaxation
/// time with a single relaxation time (Λ_T = (τ_T − ½)², the scheme before #629), with fixed
/// values of Λ_T and with the scheme as shipped (`thermal_magic`), at Ra = 10³ on 32 and 64 cells
/// and at Ra = 10⁵ on 64. Each row gives the errors of the mean Nusselt number and of the largest
/// and smallest local ones, read as the probe reads them and from the lattice's own link heat.
/// The table and the choice are in `docs/lbm-collision.md`.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn local_wall_heat_flux_with_two_relaxation_times() {
    let bgk = |tau: f64| (tau - 0.5) * (tau - 0.5);
    let shipped = fieldcraft_lbm::setup::thermal_magic;
    let fixed: [(&str, f64); 4] = [("1/24", 1.0 / 24.0), ("1/12", 1.0 / 12.0), ("1/6", 1.0 / 6.0), ("1/4", 0.25)];
    let all = [Some(0.6), Some(0.8), Some(1.0), Some(1.5), None];
    println!("| Ra | cells | Λ_T | τ_T | Nu | Nu_max | Nu_min | Nu_max (links) | Nu_min (links) |");
    println!("|---|---|---|---|---|---|---|---|---|");
    // (Ra, cells, the single relaxation time's and the shipped scheme's τ_T, the fixed Λ_T's,
    // reference Nu, Nu_max, Nu_min)
    let plan: [(f64, usize, &[Option<f64>], &[Option<f64>], [f64; 3]); 3] = [
        (1e3, 32, &all, &[Some(0.6), Some(1.5)], [1.118, 1.505, 0.692]),
        (1e3, 64, &[Some(0.8), Some(1.5)], &[Some(1.5)], [1.118, 1.505, 0.692]),
        (1e5, 64, &[Some(0.8)], &[Some(0.8)], [4.519, 7.717, 0.729]),
    ];
    // the largest change of each error over the relaxation times
    let spread = |errors: &[[f64; 5]]| errors.iter().fold([0.0f64; 5], |d, e| std::array::from_fn(|k| d[k].max((e[k] - errors[0][k]).abs())));
    let mut spreads = Vec::new();
    for (ra, n, taus, fixed_taus, reference) in plan {
        let row = |label: &str, tau: Option<f64>, magic: &dyn Fn(f64) -> f64| {
            let (tau, nu_mean, [probe, link]) = cavity_wall_flux(ra, n, tau, magic);
            let values = [nu_mean, extreme(&probe, 1.0), extreme(&probe, -1.0), extreme(&link, 1.0), extreme(&link, -1.0)];
            let refs = [reference[0], reference[1], reference[2], reference[1], reference[2]];
            let e: [f64; 5] = std::array::from_fn(|k| 100.0 * (values[k] / refs[k] - 1.0));
            let cells: Vec<String> = e.iter().map(|e| format!("{e:+.2} %")).collect();
            println!("| {ra:.0e} | {n}² | {label} | {tau:.3} | {} |", cells.join(" | "));
            // the links' heat sums to the wall's
            let link_mean = link.iter().sum::<f64>() / n as f64;
            assert!((link_mean / nu_mean - 1.0).abs() < 1e-9, "{link_mean} vs {nu_mean}");
            e
        };
        let single: Vec<[f64; 5]> = taus.iter().map(|&tau| row("(τ_T − ½)²", tau, &bgk)).collect();
        for &(label, magic) in &fixed {
            for &tau in fixed_taus {
                row(label, tau, &|_| magic);
            }
        }
        // the shipped scheme over the τ_T at which Λ_T is fixed (from 0.79)
        let ours: Vec<[f64; 5]> = taus.iter().filter(|t| t.is_none_or(|t| t > 0.79)).map(|&tau| row("shipped", tau, &shipped)).collect();
        if ours.len() > 1 {
            spreads.push((ra, n, spread(&ours), spread(&single)));
        }
    }
    for (ra, n, ours, single) in spreads {
        println!(
            "Ra {ra:.0e}, {n} cells: change over τ_T of each error (percentage points), shipped {ours:.2?}, single relaxation time {single:.2?}"
        );
        // With Λ_T fixed the wall closure does not depend on τ_T. What is left is the linear
        // equilibrium's numerical diffusion along the flow, 4 u² of the diffusivity in lattice
        // units (`a_drifting_sine_decays_at_the_diffusivity`), which grows with the time step:
        // 0.33 % at the fastest flow at τ_T = 1.5 on 32 cells.
        assert!(ours.iter().all(|d| *d < 0.25), "Ra {ra:.0e}, {n} cells: {ours:?}");
        let largest = |d: &[f64; 5]| d.iter().skip(1).copied().fold(0.0, f64::max);
        assert!(largest(&single) > 10.0 * largest(&ours), "Ra {ra:.0e}, {n} cells: {single:?} vs {ours:?}");
    }
}
