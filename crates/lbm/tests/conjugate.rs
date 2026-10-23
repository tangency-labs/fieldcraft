//! Conjugate heat transfer: heat conducting through solid obstacles, coupled to the temperature
//! lattice of the fluid, against closed-form solutions and the heat balance.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{BoxFace, FaceCondition, FluidDomain, FluidThermal, Model, Outflow, Procedure, Profile, ThermalCondition};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

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
    th.solids = vec![conducting("Slab", power)];
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

/// A column across y like [`slab_column`] with obstacles `Die` (material `Solid`) and `Lid`
/// (material `Other`) conducting, joined by `contact`, the plates held at 0 and 1.
fn contact_column(n: usize, power: f64, contact: fieldcraft_model::ThermalContact) -> FluidDomain {
    let mut d = slab_column(n, [0.0, 1.0], power);
    d.obstacles = vec!["Die".into(), "Lid".into()];
    let th = d.thermal.as_mut().unwrap();
    th.solids =
        vec![conducting("Die", power), fieldcraft_model::ConductingSolid { body: "Lid".into(), material: "Other".into(), power: 0.0, limit: None }];
    th.contacts = vec![contact];
    d
}

/// The [`contact_column`]'s model: fluid diffusivity `alpha`, `Die` of conductivity `k[0]`,
/// `Lid` of `k[1]`, both of the fluid's heat capacity.
fn contact_model(alpha: f64, k: [f64; 2], d: FluidDomain, procedure: Procedure) -> Model {
    let mut m = Model::default();
    m.add_material("Fluid", hot_fluid(1.0, 0.1, alpha)).unwrap();
    add_solid_material(&mut m, "Solid", 1.0, [k[0]; 3], 1.0);
    add_solid_material(&mut m, "Other", 1.0, [k[1]; 3], 1.0);
    m.add_fluid_domain("Flow", d).unwrap();
    m.add_step(fieldcraft_model::Step::new("Run", procedure), None).unwrap();
    m
}

fn paste(resistance: Option<f64>, layer: Option<(f64, f64)>) -> fieldcraft_model::ThermalContact {
    fieldcraft_model::ThermalContact {
        bodies: ["Die".into(), "Lid".into()],
        resistance,
        thickness: layer.map(|l| l.0),
        conductivity: layer.map(|l| l.1),
    }
}

/// A two-layer wall (a die and its lid) with a thermal interface layer between them, between
/// plates held at 0 and 1, with fluid layers outside or filling the box: the steady temperature
/// is linear in each layer and jumps by `q R` across the contact, `q = ΔT / (Σ L / k + R)`
/// (Incropera et al. 2011, §3.1.4). The face between the two bodies carries the series of the
/// two half cells and the contact, so the discrete solution is exact: it matches to round-off,
/// on two grids, with the resistance given as such or as a layer's thickness and conductivity
/// (16 cells, R = 1.25: q = 0.076923, the drop 0.096154 = q R to 1e-14).
#[test]
fn a_contact_resistance_drops_the_temperature_by_q_r() {
    let (alpha, kf) = (0.05, 0.05);
    let k = [0.2, 0.5];
    // (cells, Die, Lid, contact)
    for (n, die, lid, contact) in [
        (16, [0.25, 0.5], [0.5, 0.75], paste(Some(1.25), None)),
        (16, [-1.0, 0.5], [0.5, 2.0], paste(None, Some((0.05, 0.1)))),
        (32, [-1.0, 0.375], [0.375, 2.0], paste(Some(0.1), None)),
    ] {
        let r = contact.resistance();
        let dx = 1.0 / n as f64;
        let dt = 0.125 * dx * dx / alpha;
        let m = contact_model(alpha, k, contact_column(n, 0.0, contact), step(12.0 / alpha, dt, None));
        let obstacles = [
            obstacle("Die", cuboid([-1.0, die[0], -1.0], [2.0, die[1], 2.0]), &CUBOID_FACES),
            obstacle("Lid", cuboid([-1.0, lid[0], -1.0], [2.0, lid[1], 2.0]), &CUBOID_FACES),
        ];
        let (setup, mut flow) = prepare(&m, &obstacles);
        run_all(&setup, &mut flow);
        let (lo, hi) = (die[0].max(0.0), lid[1].min(1.0));
        let layers = [(lo, kf), (die[1] - lo, k[0]), (hi - lid[0], k[1]), (1.0 - hi, kf)];
        let q = 1.0 / (layers.iter().map(|(l, k)| l / k).sum::<f64>() + r);
        // linear in every layer, the jump q R at the contact (y = die[1])
        let exact = |y: f64| layered(y, 0.0, q, &layers) + if y > die[1] { q * r } else { 0.0 };
        let mut err: f64 = 0.0;
        for j in 0..n {
            let y = (j as f64 + 0.5) * dx;
            err = err.max((flow.temperature(j).unwrap() - exact(y)).abs());
        }
        // the drop across the contact: each side's profile extrapolated to the interface
        let b = (die[1] / dx).round() as usize;
        let a = b - 1;
        let side_a = flow.temperature(a).unwrap() + q * 0.5 * dx / k[0];
        let side_b = flow.temperature(b).unwrap() - q * 0.5 * dx / k[1];
        let drop = side_b - side_a;
        let hs = setup.thermal.as_ref().unwrap().heat_scale;
        let (cold, hot) = (flow.heat()[BoxFace::YMin.index()] * hs, flow.heat()[BoxFace::YMax.index()] * hs);
        println!(
            "two-layer wall with a contact, {n} cells, Die {die:?}, Lid {lid:?}, R {r}: q {q:.6}, drop {drop:.9} vs q R {:.9}, max error {err:.2e}; plates {cold:.6e} {hot:.6e} (q A = {:.6e})",
            q * r,
            q * dx * dx
        );
        assert!(err < 1e-9, "{n} {die:?}: {err}");
        assert!((drop - q * r).abs() < 1e-9, "{drop} vs {}", q * r);
        assert!((cold - q * dx * dx).abs() < 1e-9 * q * dx * dx && (cold + hot).abs() < 1e-9 * cold, "{cold} {hot}");
    }
}

/// The die powered, the contact between it and the lid, run from rest: every step the solids
/// gain what the die generates less what the plates and the fluid take, to round-off; a contact
/// between solids that do not touch is refused, naming it.
#[test]
fn heat_crosses_a_contact_conserved_and_a_contact_must_touch() {
    let (n, alpha) = (16usize, 0.05);
    let dx = 1.0 / n as f64;
    let dt = 0.125 * dx * dx / alpha;
    let power = 2.0 * dx * dx;
    let m = contact_model(alpha, [0.2, 0.5], contact_column(n, power, paste(Some(0.5), None)), step(400.0 * dt, dt, None));
    let obstacles = |gap: f64| {
        [
            obstacle("Die", cuboid([-1.0, 0.25, -1.0], [2.0, 0.5, 2.0]), &CUBOID_FACES),
            obstacle("Lid", cuboid([-1.0, 0.5 + gap, -1.0], [2.0, 0.75, 2.0]), &CUBOID_FACES),
        ]
    };
    let (setup, mut flow) = prepare(&m, &obstacles(0.0));
    let content = |f: &Flow| f.thermal().unwrap().solid().unwrap().heat_content();
    let into_solids = |f: &Flow| f.groups().iter().zip(f.heat()).filter(|(g, _)| g.body.is_some()).map(|(_, h)| h).sum::<f64>();
    let lattice_power = setup.thermal.as_ref().unwrap().solids[0].power;
    let mut worst: f64 = 0.0;
    for _ in 0..setup.steps {
        let before = content(&flow);
        flow.step().unwrap();
        worst = worst.max((content(&flow) - before - lattice_power - into_solids(&flow)).abs() / lattice_power);
    }
    let t = &flow.thermal().unwrap().solid().unwrap().temperature;
    println!(
        "powered die across a contact: {} steps, balance within {worst:.1e} of the power, solids {:.4}..{:.4}",
        setup.steps,
        t.iter().fold(f64::INFINITY, |a, b| a.min(*b)),
        t.iter().fold(0.0, |a: f64, b| a.max(*b))
    );
    assert!(worst < 1e-9, "{worst}");
    let setup = fieldcraft_lbm::Setup::new(&m, None, &obstacles(0.125), &fieldcraft_lbm::Limits::default()).unwrap();
    let e = Flow::new(&setup, &obstacles(0.125)).err().map(|e| e.to_string()).unwrap_or_default();
    assert!(e.contains("thermal contact `Die–Lid`") && e.contains("share no cell face"), "{e}");
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
        th.solids = vec![conducting("Cube", 1.0)];
        d.thermal = Some(th);
        let rule = steady.then(|| fieldcraft_model::Steady {
            monitors: vec![],
            tolerance: 1e-9,
            window: None,
            min_time: None,
            divergence: None,
            averaging: None,
        });
        let procedure = Procedure::Fluid {
            domain: "Flow".into(),
            end_time: Some(steps as f64 * dt),
            steps: None,
            time_step: Some(dt),
            frames: 1,
            magic: None,
            collision: None,
            turbulence: None,
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
    th.solids = vec![conducting("Block", 0.01)];
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
    th.solids = vec![conducting("Block", 0.01)];
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
        turbulence: None,
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
    let rule = fieldcraft_model::Steady { monitors: vec![], tolerance: 1e-4, window: Some(25.0), min_time: None, divergence: None, averaging: None };
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

/// The de Vahl Davis cavity (hot wall x = 0 at 1, cold wall x = H at 0, adiabatic top and
/// bottom, Pr 0.71, `n` cells across, one cell deep and periodic in z) with a centred square
/// body of side H / 2 conducting `kappa` times the fluid, its heat capacity the fluid's.
/// Returns the Nusselt number of the hot wall and the heats into the hot and cold walls.
fn conjugate_cavity(ra: f64, kappa: f64, n: usize, u_max: f64, diffusive_times: f64) -> (f64, f64, f64) {
    use fieldcraft_lbm::{Limits, run};
    let (h, pr, alpha) = (1.0, 0.71, 1.0);
    let dx = h / n as f64;
    let nu = pr * alpha;
    let g_beta = ra * nu * alpha / (h * h * h);
    let mut d = domain([0.0; 3], [h, h, dx], dx);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Body".into()];
    let mut t = FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -g_beta, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    t.solids = vec![conducting("Body", 0.0)];
    d.thermal = Some(t);
    let dt = (0.05 * dx / (u_max * alpha / h)).min(dx * dx / (4.0 * alpha));
    let m = cht_model([1.0, nu, alpha], (1.0, [kappa * alpha; 3], 1.0), d, step(diffusive_times * h * h / alpha, dt, None));
    let body = [obstacle("Body", cuboid([0.25, 0.25, -1.0], [0.75, 0.75, 1.0]), &CUBOID_FACES)];
    let setup = fieldcraft_lbm::Setup::new(&m, None, &body, &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &body).unwrap();
    let r = run(&setup, &mut flow, &Limits::default()).unwrap();
    let steady = |name: &str| {
        let w = r.heat.iter().find(|w| w.name == name && w.body.is_none()).unwrap();
        let m = w.power.len();
        let tail = &w.power[m - m / 10..];
        tail.iter().sum::<f64>() / tail.len() as f64
    };
    let (hot, cold) = (steady("x_min"), steady("x_max"));
    let nusselt = -hot * h / (alpha * h * dx);
    println!("conjugate cavity Ra {ra:.0e}, κ {kappa}, {n} cells: Nu = {nusselt:.4}, hot {hot:.6e}, cold {cold:.6e}, {} steps", setup.steps);
    (nusselt, hot, cold)
}

/// The conduction limit of the House et al. cavity (no flow), every run: the hot wall's Nusselt
/// number against the paper's 0.7069, 1 and 1.4138 for conductivities 0.2, 1 and 5 times the
/// fluid's, and the heat balances. (On 32 cells the low ratio is within 0.06 % and the high one
/// within 0.5 %: the heat crowds into the body's corners.)
#[test]
fn a_conducting_body_in_a_cavity_without_flow_matches_house_et_al() {
    for (kappa, published) in [(0.2, 0.7069), (1.0, 1.0), (5.0, 1.4138)] {
        let (nu, hot, cold) = conjugate_cavity(0.0, kappa, 16, 1.0, 1.0);
        assert!((nu / published - 1.0).abs() < 0.015, "κ {kappa}: Nu {nu} vs {published}");
        assert!((hot + cold).abs() < 1e-6 * cold, "κ {kappa}: {hot} vs {cold}");
    }
}

/// House, Beckermann and Smith, Numer. Heat Transfer A 18 (1990) 213–225, Table 3: natural
/// convection at Ra = 10⁵, Pr = 0.71, round a centred body of half the cavity's side conducting
/// 0.2, 1 and 5 times the fluid: the published hot-wall Nusselt numbers are 4.6257, 4.5061 and
/// 4.3221. The grid converges at second order from above: on 48 cells 4.6511, 4.5299 and 4.3433
/// (+0.5 %), on 64 cells 4.6356, 4.5151 and 4.3297 (+0.2 %); extrapolated from those two,
/// 4.6157, 4.4961 and 4.3122, −0.22 % from the published values. The heat into the hot wall is the heat out of the cold one
/// to 1e-7. The test runs 48 cells within 0.7 %.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn conjugate_cavity_matches_house_et_al() {
    for (kappa, published) in [(0.2, 4.6257), (1.0, 4.5061), (5.0, 4.3221)] {
        let (nu, hot, cold) = conjugate_cavity(1e5, kappa, N_HEAVY, 68.6, 1.0);
        assert!((nu / published - 1.0).abs() < 0.007, "κ {kappa}: Nu {nu} vs {published}");
        assert!((hot + cold).abs() < 1e-4 * cold, "κ {kappa}: {hot} vs {cold}");
    }
}

/// Cells across the heavy cavity (48: 2.5 minutes in release on one core).
const N_HEAVY: usize = 48;

/// Rayleigh's square array of cylinders, conduction only: the effective conductivity of a square
/// array (period 1) of cylinders of area fraction `f` conducting `kappa` times the matrix, from
/// the lattice-sum series of Y. A. Godin, J. Math. Phys. 53:063703 (2012) (arXiv:1201.1419; the
/// f⁸ pole term from arXiv:1301.0110), whose square-array lattice
/// sums `S₂ = 3.15121`, `S₄ = 4.25577` give the coefficients of W. T. Perrins, D. R. McKenzie,
/// R. C. McPhedran, Proc. R. Soc. Lond. A 369:207–225 (1979):
/// `1 + 2fρ / (1 − fρ − 0.305827 ρ² f⁴ / (1 − 1.402958 ρ² f⁸) − 0.013362 ρ² f⁸)`,
/// `ρ = (κ − 1)/(κ + 1)`.
fn square_array(f: f64, kappa: f64) -> f64 {
    let r = (kappa - 1.0) / (kappa + 1.0);
    1.0 + 2.0 * f * r / (1.0 - f * r - 0.305827 * r * r * f.powi(4) / (1.0 - 1.402958 * r * r * f.powi(8)) - 0.013362 * r * r * f.powi(8))
}

/// One period of the array, a cylinder of radius 0.3 at its centre (`n` cells across, one deep):
/// the faces between cylinders across the gradient are isotherms and the faces along it carry no
/// heat (mirror symmetry), so x = 0 held at 1, x = 1 at 0 and adiabatic y faces are exact.
/// Returns the effective conductivity over the matrix's, and the series at the voxelised area
/// fraction (the staircase changes the area; that change is factored out, so what remains is
/// the interface's own error).
fn cylinder_array(n: usize, kappa: f64) -> (f64, f64) {
    let (alpha, r) = (1.0, 0.3);
    let dx = 1.0 / n as f64;
    let mut d = domain([0.0; 3], [1.0, 1.0, dx], dx);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.obstacles = vec!["Rod".into()];
    let mut th = FluidThermal::default();
    th.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    th.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    th.solids = vec![conducting("Rod", 0.0)];
    d.thermal = Some(th);
    // τ_T = 1; the rods' heat capacity at most the fluid's (a larger one only slows the settling:
    // the steady state does not depend on it), and no less than κ times it, so they diffuse at
    // least as fast as the fluid
    let dt = dx * dx / (8.0 * alpha);
    let m = cht_model([1.0, 0.1, alpha], (kappa.min(1.0), [kappa * alpha; 3], 1.0), d, step(2.0 / alpha, dt, None));
    let rod = [obstacle("Rod", cylinder([0.5, 0.5], r, -1.0, 1.0, 256), &["side", "bottom", "top"])];
    let (setup, mut flow) = prepare(&m, &rod);
    run_all(&setup, &mut flow);
    let hs = setup.thermal.as_ref().unwrap().heat_scale;
    let (hot, cold) = (-flow.heat()[BoxFace::XMin.index()] * hs, flow.heat()[BoxFace::XMax.index()] * hs);
    let effective = cold / (alpha * dx);
    let f = flow.thermal().unwrap().solid().unwrap().cells.len() as f64 / (n * n) as f64;
    let series = square_array(f, kappa);
    println!(
        "cylinder array, {n} cells, κ {kappa}: k_eff {effective:.6} vs {series:.6} at f {f:.5} ({:+.2} %), heat in {hot:.6e} out {cold:.6e}",
        100.0 * (effective / series - 1.0)
    );
    assert!((hot - cold).abs() < 1e-6 * cold, "{hot} vs {cold}");
    (effective, series)
}

/// The array on 16 cells, every run: conducting rods (κ 10) within 1 % of the series, insulating
/// ones (κ 0.1) within 6 %: for them the staircase is a first-order error (the heavy test follows
/// it down). Measured with the thermal lattice's Λ_T = 1/12: +0.55 % and −5.66 %. (With a single
/// relaxation time these depended on τ_T: at this test's τ_T = 1, Λ_T = 1/4, +0.18 % and −3.85 %;
/// at τ_T = 0.6, Λ_T = 0.01, +0.88 % and −12.4 %.) The series itself against its value at the round area fraction, and Keller's
/// reciprocity (a square array's k(κ) k(1/κ) = 1).
#[test]
fn a_square_array_of_cylinders_conducts_as_rayleigh_s_series() {
    let f = std::f64::consts::PI * 0.09;
    assert!((square_array(f, 10.0) - 1.602942).abs() < 1e-6 && (square_array(f, 0.1) - 0.623853).abs() < 1e-6);
    assert!((square_array(f, 10.0) * square_array(f, 0.1) - 1.0).abs() < 1e-12);
    let (k10, s10) = cylinder_array(16, 10.0);
    assert!((k10 / s10 - 1.0).abs() < 0.01, "{k10} vs {s10}");
    let (k01, s01) = cylinder_array(16, 0.1);
    assert!((k01 / s01 - 1.0).abs() < 0.06, "{k01} vs {s01}");
}

/// Insulating rods (κ 0.1) on 32 and 64 cells: −2.99 % and −1.21 % of the series (Λ_T = 1/12),
/// the staircase interface converging at first order or better; conducting rods (κ 10) on 32
/// cells within 0.6 % (+0.54 %). With a single relaxation time at this test's τ_T = 1
/// (Λ_T = 1/4) they were −2.2 %, −0.85 % and +0.31 %, at τ_T = 0.6 (Λ_T = 0.01) −4.7 % and +0.77 %
/// on 32 cells.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn insulating_cylinders_converge_at_first_order() {
    let (k10, s10) = cylinder_array(32, 10.0);
    assert!((k10 / s10 - 1.0).abs() < 0.006, "{k10} vs {s10}");
    let errors: Vec<f64> = [32, 64]
        .iter()
        .map(|n| {
            let (k, s) = cylinder_array(*n, 0.1);
            (k / s - 1.0).abs()
        })
        .collect();
    let order = orders(&errors, 2.0);
    println!("insulating cylinders: errors {errors:?}, orders {order:?}");
    assert!(errors[1] < 0.013, "{errors:?}");
    assert!(order.iter().all(|o| *o > 0.9), "{order:?}");
}

/// The heat-generating slab of [`a_heat_generating_slab_balances_and_converges`] (16 cells,
/// plates at 0 °C), started at 10 °C everywhere, against temperature limits either side of its
/// steady peak (7.625 °C at its centre). Run to steady state its peak is the steady value, the
/// cell centres 1/32 off the centre 0.0078 lower and the scheme within 2 % of it; a short
/// transient run reports the highest over its frames, at or above the 10 °C it starts at (its
/// power heats it before the plates draw the heat out).
#[test]
fn a_heated_slab_passes_and_fails_limits_either_side_of_its_peak() {
    use fieldcraft_lbm::{Limits, StopReason, solve_fluid};
    let (n, alpha, ks) = (16usize, 0.05, 0.5);
    let dx = 1.0 / n as f64;
    let dt = dx * dx / (8.0 * alpha);
    let exact = 7.625 - 8.0 / (2.0 * ks) * (0.5 * dx).powi(2);
    let solve = |limit: f64, steady: bool| {
        let mut d = slab_column(n, [0.0, 0.0], 2.0 * dx * dx);
        let th = d.thermal.as_mut().unwrap();
        th.initial_temperature = kelvin(10.0);
        th.solids[0].limit = Some(kelvin(limit));
        let mut procedure = step(if steady { 8.0 / alpha } else { 1.0 }, dt, None);
        if let Procedure::Fluid { steady: rule, .. } = &mut procedure
            && steady
        {
            *rule = Some(fieldcraft_model::Steady {
                monitors: vec![],
                tolerance: 1e-5,
                window: Some(2.0),
                min_time: None,
                divergence: None,
                averaging: None,
            });
        }
        let m = cht_model([1.0, 0.1, alpha], (1.0, [ks; 3], 1.0), d, procedure);
        solve_fluid(&m, None, &[slab_obstacle([0.375, 0.625])], &Limits::default()).unwrap()
    };
    for (limit, pass) in [(exact * 1.05, true), (exact * 0.95, false)] {
        let r = solve(limit, true);
        assert_eq!(r.steady.as_ref().map(|s| s.reason), Some(StopReason::Steady));
        let c = &r.components()[0];
        let peak = c.peak.unwrap();
        assert!(peak.steady && peak.time == r.frames.last().unwrap().time, "{peak:?}");
        let rise = peak.temperature - kelvin(0.0);
        assert!((rise - exact).abs() < 0.02 * exact, "peak {rise} vs {exact}");
        assert!((c.margin.unwrap() - (kelvin(limit) - peak.temperature)).abs() < 1e-9, "{c:?}");
        assert_eq!((c.ok, r.components_ok()), (Some(pass), Some(pass)), "limit {limit}: {c:?}");
    }
    // transient: the hottest of the frames, above a limit the steady slab keeps under
    let r = solve(exact * 1.05, false);
    let c = &r.components()[0];
    let peak = c.peak.unwrap();
    let hottest = |f: &fieldcraft_lbm::Frame| r.solids[0].cells.iter().map(|k| f.temperature[*k as usize]).fold(f32::MIN, f32::max);
    assert!(!peak.steady && r.frames.iter().all(|f| f64::from(hottest(f)) <= peak.temperature - kelvin(0.0) + 1e-9), "{peak:?}");
    assert!(r.frames.iter().any(|f| f.time == peak.time) && peak.temperature >= kelvin(10.0), "{peak:?}");
    assert_eq!((c.ok, r.components_ok()), (Some(false), Some(false)));
    // without a limit: a peak, no verdict
    let mut r = r;
    r.solids[0].limit = None;
    let c = &r.components()[0];
    assert!(c.peak.is_some() && c.margin.is_none() && c.ok.is_none() && r.components_ok().is_none(), "{c:?}");
}
