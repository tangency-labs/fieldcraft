//! Conjugate heat transfer: heat conducting through solid obstacles, coupled to the temperature
//! lattice of the fluid, against closed-form solutions and the heat balance.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{
    BoxFace, ConductingSolid, FaceCondition, FluidDomain, FluidThermal, Model, Outflow, Procedure, Profile, Step, ThermalCondition,
};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// A model with fluid `Fluid` (density ρ, kinematic viscosity ν, diffusivity α: conductivity
/// ρ α, unit specific heat), conducting material `Solid` (density ρ_s, conductivity k_s per axis,
/// specific heat c_s), domain `Flow` and step `Run`.
fn cht_model(fluid: [f64; 3], solid: (f64, [f64; 3], f64), domain: FluidDomain, procedure: Procedure) -> Model {
    let [rho, nu, alpha] = fluid;
    let mut m = Model::default();
    m.add_material("Fluid", hot_fluid(rho, nu, alpha)).unwrap();
    add_solid_material(&mut m, "Solid", solid.0, solid.1, solid.2);
    m.add_fluid_domain("Flow", domain).unwrap();
    m.add_step(Step::new("Run", procedure), None).unwrap();
    m
}

/// A column of `n` cells across y ∈ [0, 1], periodic in x and z, the faces at y = 0 and 1 held
/// at `t[0]` and `t[1]`, obstacle `Slab` conducting with `power`.
fn slab_column(n: usize, t: [f64; 2], power: f64) -> FluidDomain {
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Slab".into()];
    let mut th = FluidThermal::default();
    th.faces.y_min = ThermalCondition::Temperature { value: kelvin(t[0]) };
    th.faces.y_max = ThermalCondition::Temperature { value: kelvin(t[1]) };
    th.solids = vec![ConductingSolid { body: "Slab".into(), material: "Solid".into(), power }];
    d.thermal = Some(th);
    d
}

fn slab_obstacle(slab: [f64; 2]) -> fieldcraft_lbm::Obstacle {
    obstacle("Slab", cuboid([-1.0, slab[0], -1.0], [2.0, slab[1], 2.0]), &CUBOID_FACES)
}

/// The steady temperature through layers `(thickness, conductivity)` from y = 0 at `t0`, with
/// heat flux `q` in the −y direction (from the y = 1 side to the y = 0 side) and no sources.
fn layered(y: f64, t0: f64, q: f64, layers: &[(f64, f64)]) -> f64 {
    let (mut base, mut t) = (0.0, t0);
    for &(thick, k) in layers {
        if y <= base + thick {
            return t + q * (y - base) / k;
        }
        t += q * thick / k;
        base += thick;
    }
    t
}

/// Fluid, a conducting slab and fluid between held plates, no flow: the steady temperature is
/// linear in each layer with the flux continuous across both interfaces, `q = ΔT / Σ (L / k)`.
/// The interface temperature from the series resistance of the half cells makes the
/// discrete solution exact, so it matches to round-off, on two grids, with the slab in the
/// middle, against the cold plate, and with a conductivity that differs along the other axes.
#[test]
fn a_composite_wall_conducts_in_series() {
    let (alpha, kf) = (0.05, 0.05);
    // (cells, slab, conductivity, the fluid's lattice diffusivity: 1/8 is τ_T = 1, 0.005 is
    // τ_T = 0.52, where the interface links couple most strongly per unit of conductance)
    for (n, slab, k, lattice) in [
        (16, [0.25, 0.5], [0.2; 3], 0.125),
        (32, [0.25, 0.5], [0.2; 3], 0.125),
        (16, [-1.0, 0.25], [0.2; 3], 0.125),
        (8, [-1.0, 0.25], [0.2; 3], 0.005),
        (16, [0.25, 0.5], [2.0, 0.01, 0.5], 0.125),
    ] {
        let dx = 1.0 / n as f64;
        let dt = lattice * dx * dx / alpha;
        let m = cht_model([1.0, 0.1, alpha], (1.0, k, 1.0), slab_column(n, [0.0, 1.0], 0.0), step(8.0 / alpha, dt, None));
        let obstacles = [slab_obstacle(slab)];
        let (setup, mut flow) = prepare(&m, &obstacles);
        run_all(&setup, &mut flow);
        let lo = slab[0].max(0.0);
        let layers = [(lo, kf), (slab[1] - lo, k[1]), (1.0 - slab[1], kf)];
        let q = 1.0 / layers.iter().map(|(l, k)| l / k).sum::<f64>();
        let mut err: f64 = 0.0;
        for j in 0..n {
            let y = (j as f64 + 0.5) * dx;
            err = err.max((flow.temperature(j).unwrap() - layered(y, 0.0, q, &layers)).abs());
        }
        let solid = flow.thermal().unwrap().solid().unwrap();
        // every plate passes the flux q over its area, the slab's too when it touches one
        let hs = setup.thermal.as_ref().unwrap().heat_scale;
        let (cold, hot) = (flow.heat()[BoxFace::YMin.index()] * hs, flow.heat()[BoxFace::YMax.index()] * hs);
        println!(
            "composite wall, {n} cells, slab {slab:?}, k {k:?}, lattice α {lattice}: max error {err:.2e}, {} sub-steps; plates {cold:.6e} {hot:.6e} (q A = {:.6e})",
            solid.substeps,
            q * dx * dx
        );
        assert!(err < 1e-9, "{n} {slab:?} {k:?}: {err}");
        assert!((cold - q * dx * dx).abs() < 1e-9 * q * dx * dx && (cold + hot).abs() < 1e-9 * cold, "{cold} {hot}");
        assert_eq!(solid.cells.len(), (n as f64 * (slab[1] - lo)).round() as usize);
    }
}

/// A slab generating heat between two plates held at zero: every plate takes half the power,
/// the fluid layers are linear and the slab parabolic. The cell-centred conduction holds the
/// parabola exactly; the half cell at each interface assumes a linear profile, an error of
/// `s Δx² / (8 k_s)` that falls with the square of the spacing.
#[test]
fn a_heat_generating_slab_balances_and_converges() {
    let (alpha, kf, ks) = (0.05, 0.05, 0.5);
    let slab = [0.375, 0.625];
    let power_per_area = 2.0;
    let mut errors = Vec::new();
    for n in [16usize, 32] {
        let dx = 1.0 / n as f64;
        let dt = dx * dx / (8.0 * alpha);
        // the slab spans the cell's x and z extent: its power per area of the plates
        let power = power_per_area * dx * dx;
        let m = cht_model([1.0, 0.1, alpha], (1.0, [ks; 3], 1.0), slab_column(n, [0.0, 0.0], power), step(8.0 / alpha, dt, None));
        let (setup, mut flow) = prepare(&m, &[slab_obstacle(slab)]);
        run_all(&setup, &mut flow);
        let thick = slab[1] - slab[0];
        let s = power_per_area / thick;
        let q = 0.5 * power_per_area;
        let t_i = q * slab[0] / kf;
        let exact = |y: f64| {
            if y < slab[0] {
                q * y / kf
            } else if y > slab[1] {
                q * (1.0 - y) / kf
            } else {
                t_i + s / (2.0 * ks) * ((0.5 * thick).powi(2) - (y - 0.5).powi(2))
            }
        };
        let peak = exact(0.5);
        let mut err: f64 = 0.0;
        for j in 0..n {
            let y = (j as f64 + 0.5) * dx;
            err = err.max((flow.temperature(j).unwrap() - exact(y)).abs() / peak);
        }
        // the plates (box faces y_min, y_max) receive the whole power at steady state
        let hs = setup.thermal.as_ref().unwrap().heat_scale;
        let heat = flow.heat();
        let plates = (heat[BoxFace::YMin.index()] + heat[BoxFace::YMax.index()]) * hs;
        println!("heat-generating slab, {n} cells: max error {err:.2e} of the peak, plates take {plates:.6} of {power:.6}");
        assert!((plates - power).abs() < 1e-6 * power, "{n}: {plates} vs {power}");
        errors.push(err);
    }
    let order = orders(&errors, 2.0);
    println!("orders {order:?}");
    assert!(errors[0] < 0.02, "{errors:?}");
    assert!(order[0] > 1.8, "{order:?}");
}

/// A channel along x, periodic, its floor a solid wall y ∈ [0, ¼] generating heat (adiabatic
/// underneath), a Poiseuille flow driven by a body force over y ∈ [¼, 1] and the roof held at
/// zero. The flow runs along the isotherms, so the steady temperature is the conduction one:
/// linear across the fluid, all the wall's heat crossing it, and parabolic in the wall with its
/// peak on the adiabatic face. The interface sees a moving fluid beside it and must not care.
#[test]
fn a_heated_wall_under_a_channel_flow_conducts_across_it() {
    let (alpha, nu, kf, ks) = (0.05, 0.05, 0.05, 0.25);
    let n = 16usize;
    let dx = 1.0 / n as f64;
    let dt = dx * dx / (8.0 * alpha);
    let flux = 2.0;
    let wall = 0.25;
    let mut d = slab_column(n, [0.0, 0.0], flux * dx * dx);
    let th = d.thermal.as_mut().unwrap();
    th.faces.y_min = ThermalCondition::Adiabatic {};
    let (gap, u_peak) = (1.0 - wall, 0.4);
    d.acceleration = [8.0 * nu * u_peak / (gap * gap), 0.0, 0.0];
    let m = cht_model([1.0, nu, alpha], (1.0, [ks; 3], 1.0), d, step(8.0 / alpha, dt, None));
    let (setup, mut flow) = prepare(&m, &[slab_obstacle([-1.0, wall])]);
    run_all(&setup, &mut flow);
    let s = flux / wall;
    let t_i = flux * gap / kf;
    let exact = |y: f64| if y > wall { flux * (1.0 - y) / kf } else { t_i + s / (2.0 * ks) * (wall * wall - y * y) };
    let peak = exact(0.0);
    let (mut err, mut speed): (f64, f64) = (0.0, 0.0);
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        err = err.max((flow.temperature(j).unwrap() - exact(y)).abs() / peak);
        if let Some((_, u)) = flow.moments(j) {
            speed = speed.max(u[0] * setup.velocity_scale());
        }
    }
    println!("heated wall under a channel flow: peak speed {speed:.4} (Poiseuille {u_peak}), max error {err:.2e} of the peak temperature");
    assert!((speed / u_peak - 1.0).abs() < 0.02, "{speed}");
    assert!(err < 2e-3, "{err}");
}

/// A heated cube in a still 16³ box held at zero, conducting 8 times the fluid (plastic in air)
/// at τ_T = 0.52, where each fluid link takes `b = 2 w k_s / (k_s + k_f)` ≈ ¼ of the solid's
/// temperature per step: an explicit link went unstable at the steady run's scaled capacity and
/// at a heat capacity of 0.03 of the fluid's. With the links implicit both run, and the cube's
/// heat balances every step.
#[test]
fn a_low_capacity_solid_with_strong_links_stays_stable() {
    use fieldcraft_lbm::Limits;
    let (n, alpha) = (16usize, 1.0);
    let dx = 1.0 / n as f64;
    let dt = 0.005 * dx * dx / alpha;
    let steps = 800;
    for (capacity, steady) in [(1.0, true), (0.03, false)] {
        let mut d = domain([0.0; 3], [1.0; 3], dx);
        d.obstacles = vec!["Cube".into()];
        let mut th = FluidThermal::default();
        for face in BoxFace::ALL {
            *th.faces.get_mut(face) = ThermalCondition::Temperature { value: kelvin(0.0) };
        }
        th.solids = vec![ConductingSolid { body: "Cube".into(), material: "Solid".into(), power: 1.0 }];
        d.thermal = Some(th);
        let rule = steady.then(|| fieldcraft_model::Steady { monitors: vec![], tolerance: 1e-9, window: None, min_time: None, divergence: None });
        let procedure = Procedure::Fluid {
            domain: "Flow".into(),
            end_time: Some(steps as f64 * dt),
            steps: None,
            time_step: Some(dt),
            frames: 1,
            magic: None,
            collision: None,
            steady: rule,
        };
        let m = cht_model([1.0, alpha, alpha], (capacity, [8.0 * alpha; 3], 1.0), d, procedure);
        let cube = [obstacle("Cube", cuboid([0.25; 3], [0.75; 3]), &CUBOID_FACES)];
        let setup = fieldcraft_lbm::Setup::new(&m, None, &cube, &Limits::default()).unwrap();
        let ts = setup.thermal.as_ref().unwrap();
        assert!((ts.tau - 0.52).abs() < 1e-12, "{}", ts.tau);
        let mut flow = Flow::new(&setup, &cube).unwrap();
        let content = |f: &Flow| f.thermal().unwrap().solid().unwrap().heat_content();
        let power = ts.solids[0].power;
        let mut worst: f64 = 0.0;
        for _ in 0..steps {
            let before = content(&flow);
            flow.step().unwrap();
            let faces: f64 = flow.groups().iter().zip(flow.heat()).filter(|(g, _)| g.body.is_some()).map(|(_, h)| h).sum();
            worst = worst.max((content(&flow) - before - power - faces).abs() / power);
        }
        let peak = flow.thermal().unwrap().solid().unwrap().temperature.iter().fold(0.0, |m: f64, t| m.max(*t));
        println!(
            "low-capacity cube, capacity {capacity} (scaled to {:.3}), steady {steady}: {steps} steps, peak {peak:.4e}, balance {worst:.1e} of the power",
            ts.solids[0].capacity
        );
        assert!(peak.is_finite() && peak > 0.0, "{peak}");
        assert!(worst < 1e-9, "{worst}");
    }
}

/// A powered solid thinner than the spacing covers no cell centre: it is refused (its power would
/// vanish while the report still showed it), while an unpowered one passes with the generic
/// warning.
#[test]
fn a_powered_solid_that_covers_no_cell_is_refused() {
    let thin = [0.44, 0.46];
    let n = 16;
    let dx = 1.0 / n as f64;
    let dt = dx * dx / (8.0 * 0.05);
    for (power, refused) in [(1.0, true), (0.0, false)] {
        let m = cht_model([1.0, 0.1, 0.05], (1.0, [0.2; 3], 1.0), slab_column(n, [0.0, 1.0], power), step(1.0, dt, None));
        let obstacles = [slab_obstacle(thin)];
        let setup = fieldcraft_lbm::Setup::new(&m, None, &obstacles, &fieldcraft_lbm::Limits::default()).unwrap();
        match Flow::new(&setup, &obstacles) {
            Err(e) => {
                assert!(refused, "{e}");
                assert!(e.to_string().contains("covers no cell centre"), "{e}");
            }
            Ok(f) => assert!(!refused && f.empty_obstacles() == ["Slab".to_string()]),
        }
    }
}

/// A heated block in a channel flow: whatever the block's faces take from the fluid plus its
/// power is what its cells gain, step by step, to round-off; at steady state its faces give
/// the fluid its whole power.
#[test]
fn a_heated_block_in_a_flow_conserves_heat() {
    let (n, alpha, nu) = (8usize, 0.02, 0.02);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [3.0, 1.0, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Block".into()];
    let mut th = FluidThermal { initial_temperature: kelvin(1.0), ..FluidThermal::default() };
    th.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    th.solids = vec![ConductingSolid { body: "Block".into(), material: "Solid".into(), power: 0.01 }];
    d.thermal = Some(th);
    let dt = 0.25 * dx * dx / nu;
    let m = cht_model([1.0, nu, alpha], (2.0, [0.5; 3], 1.0), d, step(400.0, dt, None));
    let block = obstacle("Block", cuboid([0.75, 0.375, -1.0], [1.25, 0.625, 1.0]), &CUBOID_FACES);
    let (setup, mut flow) = prepare(&m, &[block]);
    let ts = setup.thermal.as_ref().unwrap();
    let power = ts.solids[0].power;
    let content = |f: &Flow| f.thermal().unwrap().solid().unwrap().heat_content();
    let block_heat = |f: &Flow| f.groups().iter().zip(f.heat()).filter(|(g, _)| g.body.as_deref() == Some("Block")).map(|(_, h)| h).sum::<f64>();
    let mut worst: f64 = 0.0;
    for k in 0..setup.steps {
        let before = content(&flow);
        flow.step().unwrap();
        let gained = content(&flow) - before;
        let expected = power + block_heat(&flow);
        worst = worst.max((gained - expected).abs() / power);
        if k == 10 {
            assert!(block_heat(&flow).abs() > 0.0, "the block exchanges heat with the fluid");
        }
    }
    let out = -block_heat(&flow);
    println!("heated block: {} steps, per-step balance within {worst:.1e} of the power; faces give {out:.6e} of {power:.6e}", setup.steps);
    assert!(worst < 1e-9, "{worst}");
    assert!((out - power).abs() < 1e-3 * power, "{out} vs {power}");
}

/// The heated block of [`a_heated_block_in_a_flow_conserves_heat`] with solid heat capacity
/// `capacity` (times the fluid's), run for `end_time` or to steady state.
fn heated_block(capacity: f64, end_time: f64, steady: Option<fieldcraft_model::Steady>) -> (Model, Vec<fieldcraft_lbm::Obstacle>) {
    let (n, alpha, nu) = (8usize, 0.02, 0.02);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [3.0, 1.0, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Block".into()];
    let mut th = FluidThermal { initial_temperature: kelvin(1.0), ..FluidThermal::default() };
    th.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    th.solids = vec![ConductingSolid { body: "Block".into(), material: "Solid".into(), power: 0.01 }];
    d.thermal = Some(th);
    let dt = 0.25 * dx * dx / nu;
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(dt),
        frames: 2,
        magic: None,
        collision: None,
        steady,
    };
    let m = cht_model([1.0, nu, alpha], (capacity, [0.5; 3], 1.0), d, procedure);
    (m, vec![obstacle("Block", cuboid([0.75, 0.375, -1.0], [1.25, 0.625, 1.0]), &CUBOID_FACES)])
}

/// A block a thousand times slower to heat than the fluid around it (aluminium in air is two
/// thousand): run to steady state, its heat capacity is scaled down, it stops in a small
/// fraction of its own time constant, and its temperatures are those a fast block reaches in a
/// long fixed run (the steady temperatures do not depend on the heat capacity).
#[test]
fn a_steady_run_scales_a_slow_solid_and_lands_on_the_same_temperatures() {
    use fieldcraft_lbm::{Limits, StopReason, solve_fluid};
    let rule = fieldcraft_model::Steady { monitors: vec![], tolerance: 1e-4, window: Some(25.0), min_time: None, divergence: None };
    let (slow, obstacles) = heated_block(2000.0, 20000.0, Some(rule));
    let setup = fieldcraft_lbm::Setup::new(&slow, None, &obstacles, &Limits::default()).unwrap();
    let solid = &setup.thermal.as_ref().unwrap().solids[0];
    assert!(solid.capacity_scale < 0.05, "{}", solid.capacity_scale);
    assert!(setup.warnings.iter().any(|w| w.contains("scaled")), "{:?}", setup.warnings);
    let run = solve_fluid(&slow, None, &obstacles, &Limits::default()).unwrap();
    let report = run.steady.as_ref().unwrap();
    assert_eq!(report.reason, StopReason::Steady, "{:?}", report.message);
    // the slow block's own time constant (its capacity over its conductance to the fluid) is
    // about 500 time units, so it would take thousands to settle unscaled
    assert!(run.stats.end_time < 800.0, "{}", run.stats.end_time);
    let (fast, obstacles) = heated_block(2.0, 800.0, None);
    let reference = solve_fluid(&fast, None, &obstacles, &Limits::default()).unwrap();
    let (a, b) = (run.frames.last().unwrap(), reference.frames.last().unwrap());
    let rise = b.temperature.iter().map(|t| t - 1.0).fold(0.0, f32::max) as f64;
    let diff = a.temperature.iter().zip(&b.temperature).map(|(x, y)| (x - y).abs()).fold(0.0, f32::max) as f64;
    println!(
        "steady block: scaled by {:.3e}, stopped at t = {:.1} ({} steps); peak rise {rise:.4}, differs from the long run by {:.2e} of it",
        solid.capacity_scale,
        run.stats.end_time,
        run.stats.steps,
        diff / rise
    );
    assert!(diff < 1e-3 * rise, "{diff} of {rise}");
}
