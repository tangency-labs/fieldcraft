//! Trusses, point masses, springs and dashpots in explicit dynamics against closed-form
//! solutions. Run with `--nocapture` to see the errors.
//!
//! - A mass on a spring started with a velocity: x = (v₀/ω) sin ωt.
//! - The same with a dashpot: the amplitude decays as e^{−ξωt}, the period is 2π/ω_d with
//!   ω_d = ω√(1 − ξ²); a mass on a dashpot alone slows as v₀ e^{−ct/m}.
//! - Two masses on a two-point spring along its line: ω² = k(1/m₁ + 1/m₂), momentum conserved.
//! - A truss bar driven at one end: a wave front at c t with stress −ρ c v behind it.

use std::f64::consts::PI;

use fieldcraft_explicit::{ExplicitResult, Limits, solve_explicit};
use fieldcraft_fem::truss::axial;
use fieldcraft_model::*;

fn model(step: ExplicitDynamic) -> Model {
    let mut m = Model::default();
    m.add_step(Step::new("Dynamic", Procedure::ExplicitDynamic(ExplicitDynamic::new(step.end_time))), None).unwrap();
    // The step's settings name sets; they are put in once the sets exist (see `solve`).
    PENDING.with(|p| *p.borrow_mut() = Some(step));
    m
}

thread_local! {
    static PENDING: std::cell::RefCell<Option<ExplicitDynamic>> = const { std::cell::RefCell::new(None) };
}

/// Put the step's settings in, keeping what is active in it, and solve.
fn solve(m: &mut Model) -> ExplicitResult {
    if let Some(settings) = PENDING.with(|p| p.borrow_mut().take()) {
        let mut step = m.step("Dynamic").unwrap().clone();
        step.procedure = Procedure::ExplicitDynamic(settings);
        m.replace_step(step).unwrap();
    }
    solve_explicit(m, None, &Limits::default()).unwrap()
}

fn set(m: &mut Model, name: &str, nodes: &[u32]) -> Target {
    m.add_set(name, Set { nodes: nodes.to_vec(), ..Set::default() }).unwrap();
    Target::Set(name.into())
}

fn hold(m: &mut Model, name: &str, target: Target, values: [Option<f64>; 6]) {
    m.add_boundary_condition(name, BoundaryCondition::Displacement { target, values }).unwrap();
    m.set_boundary_condition_active("Dynamic", name, true).unwrap();
}

fn step(end: f64, scale: f64, moving: &str, velocity: [f64; 3]) -> ExplicitDynamic {
    let mut s = ExplicitDynamic::new(end);
    s.dt_scale = scale;
    s.history_nodes = Some(moving.into());
    s.history_interval = Some(end / 20_000.0);
    s.initial_velocities = vec![InitialVelocity { set: moving.into(), velocity }];
    s.bulk_viscosity = BulkViscosity { linear: 0.0, quadratic: 0.0 };
    s
}

/// The history of one direction of the first history node.
fn trace(r: &ExplicitResult, k: usize) -> Vec<(f64, f64)> {
    r.history.time.iter().zip(&r.history.displacement).map(|(t, d)| (*t, d[0][k])).collect()
}

/// Local maxima of a sampled signal (parabolic interpolation through the three samples).
fn peaks(x: &[(f64, f64)]) -> Vec<(f64, f64)> {
    x.windows(3)
        .filter(|w| w[1].1 > w[0].1 && w[1].1 >= w[2].1)
        .map(|w| {
            let (h, a, b, c) = (w[1].0 - w[0].0, w[0].1, w[1].1, w[2].1);
            let d = 0.5 * (a - c) / (a - 2.0 * b + c);
            (w[1].0 + d * h, b - 0.25 * (a - c) * d)
        })
        .collect()
}

/// A mass on a spring to the ground (and optionally a dashpot), started at velocity v₀.
fn oscillator(k: f64, m: f64, damping: Option<f64>, periods: f64, scale: f64) -> ExplicitResult {
    let omega = (k / m).sqrt();
    let v0 = 0.3;
    let mut model = model(step(periods * 2.0 * PI / omega, scale, "n", [v0, 0.0, 0.0]));
    model.set_node(1, [0.0, 0.0, 0.0]).unwrap();
    let n = set(&mut model, "n", &[1]);
    model.add_spring("k", Spring::Spring { first: n.clone(), second: None, direction: Direction::Dof(1), stiffness: k }).unwrap();
    if let Some(c) = damping {
        model.add_spring("c", Spring::Dashpot { first: n.clone(), second: None, direction: Direction::Dof(1), coefficient: c }).unwrap();
    }
    model.add_mass("m", PointMass { target: n.clone(), mass: Some(Mass::Isotropic(m)), inertia: None }).unwrap();
    hold(&mut model, "yz", n, [None, Some(0.0), Some(0.0), None, None, None]);
    solve(&mut model)
}

#[test]
fn spring_mass_oscillation() {
    let (k, m) = (200.0f64, 0.5);
    let r = oscillator(k, m, None, 3.0, 0.02);
    let omega = (k / m).sqrt();
    let amplitude = 0.3 / omega;
    let worst = trace(&r, 0).iter().map(|(t, x)| (x - amplitude * (omega * t).sin()).abs()).fold(0.0, f64::max);
    let p = peaks(&trace(&r, 0));
    let period = (p[2].0 - p[0].0) / 2.0;
    eprintln!(
        "spring–mass: worst error {:.2e} of the amplitude, period {period:.6} vs {:.6}, Δt {:.3e}, energy error {:.1e}",
        worst / amplitude,
        2.0 * PI / omega,
        r.initial_increment,
        r.energy.relative_error(r.initial_kinetic)
    );
    assert!(worst < 2e-3 * amplitude, "second-order phase error");
    assert!((period - 2.0 * PI / omega).abs() < 1e-4 * period);
    assert!(r.energy.relative_error(r.initial_kinetic) < 1e-3);
    assert!(r.warnings.iter().any(|w| w.contains("limits the time increment")), "{:?}", r.warnings);
    // The increment is the scale times the spring's 2/ω.
    assert!((r.initial_increment - 0.02 * 2.0 / omega).abs() < 1e-12 * r.initial_increment);
}

#[test]
fn damped_decay_rate() {
    let (k, m, xi) = (200.0f64, 0.5, 0.05);
    let omega = (k / m).sqrt();
    let c = 2.0 * xi * m * omega;
    let r = oscillator(k, m, Some(c), 6.0, 0.01);
    let p = peaks(&trace(&r, 0));
    let wd = omega * (1.0 - xi * xi).sqrt();
    // Decay rate from the first and last peaks, period from their spacing.
    let (first, last) = (p[0], p[p.len() - 1]);
    let rate = (first.1 / last.1).ln() / (last.0 - first.0);
    let period = (last.0 - first.0) / (p.len() - 1) as f64;
    let rate_err = (rate - xi * omega).abs() / (xi * omega);
    let period_err = (period - 2.0 * PI / wd).abs() / (2.0 * PI / wd);
    eprintln!("damped oscillator: decay rate {rate:.6} vs ξω {:.6} ({rate_err:.2e}), period error {period_err:.2e}", xi * omega);
    assert!(rate_err < 0.01, "{rate} vs {}", xi * omega);
    assert!(period_err < 1e-3);
    // The dashpot dissipates what the oscillator loses: the balance holds.
    assert!(r.energy.viscous > 0.9 * r.initial_kinetic && r.energy.relative_error(r.initial_kinetic) < 0.01, "{:?}", r.energy);

    // A dashpot alone: v = v₀ e^{−ct/m}.
    let mut model = model(step(3.0 * m / c, 0.001, "n", [0.3, 0.0, 0.0]));
    model.set_node(1, [0.0; 3]).unwrap();
    let n = set(&mut model, "n", &[1]);
    model.add_spring("c", Spring::Dashpot { first: n.clone(), second: None, direction: Direction::Vector([2.0, 0.0, 0.0]), coefficient: c }).unwrap();
    model.add_mass("m", PointMass { target: n.clone(), mass: Some(Mass::Isotropic(m)), inertia: None }).unwrap();
    hold(&mut model, "yz", n, [None, Some(0.0), Some(0.0), None, None, None]);
    let r = solve(&mut model);
    let (t, v) = (r.history.time.last().unwrap(), r.history.velocity.last().unwrap()[0][0]);
    let exact = 0.3 * (-c * t / m).exp();
    eprintln!("dashpot alone: v {v:.6e} vs {exact:.6e} at t = {t:.4}");
    assert!((v - exact).abs() < 0.01 * exact);
    // Travel: x(∞) = v₀ m / c.
    let x = r.history.displacement.last().unwrap()[0][0];
    let travel = 0.3 * m / c * (1.0 - (-c * t / m).exp());
    assert!((x - travel).abs() < 0.01 * travel, "{x} vs {travel}");
}

#[test]
fn two_masses_on_a_line_spring() {
    let (k, m1, m2) = (50.0f64, 1.0, 3.0);
    let omega = (k * (1.0 / m1 + 1.0 / m2)).sqrt();
    let dir = [0.6, 0.0, 0.8];
    let mut model = model(step(4.0 * 2.0 * PI / omega, 0.02, "a", [0.0; 3]));
    model.set_node(1, [0.0; 3]).unwrap();
    model.set_node(2, dir.map(|c| 2.0 * c)).unwrap();
    let a = set(&mut model, "a", &[1]);
    let b = set(&mut model, "b", &[2]);
    model.add_spring("k", Spring::Spring { first: a.clone(), second: Some(b.clone()), direction: Direction::Line, stiffness: k }).unwrap();
    model.add_mass("m1", PointMass { target: a, mass: Some(Mass::Isotropic(m1)), inertia: None }).unwrap();
    model.add_mass("m2", PointMass { target: b, mass: Some(Mass::Isotropic(m2)), inertia: None }).unwrap();
    // Node 1 starts towards node 2, node 2 at rest: the centre of mass drifts, the spring rings.
    PENDING.with(|p| p.borrow_mut().as_mut().unwrap().initial_velocities = vec![InitialVelocity { set: "a".into(), velocity: dir.map(|c| 0.5 * c) }]);
    let r = solve(&mut model);
    // Momentum is conserved exactly (internal forces are equal and opposite).
    let p0 = r.history.momentum[0];
    let drift = r.history.momentum.iter().map(|p| (0..3).map(|k| (p[k] - p0[k]).abs()).fold(0.0, f64::max)).fold(0.0, f64::max);
    // The relative stretch along the line oscillates at ω.
    let x1 = trace(&r, 0);
    let along: Vec<(f64, f64)> = r
        .history
        .time
        .iter()
        .zip(&r.history.displacement)
        .map(|(t, d)| (*t, d[0][0] * dir[0] + d[0][2] * dir[2] - 0.5 * m1 / (m1 + m2) * t))
        .collect();
    let p = peaks(&along);
    let period = (p[3].0 - p[0].0) / 3.0;
    eprintln!("two masses: period {period:.6} vs {:.6}, momentum drift {drift:.1e}", 2.0 * PI / omega);
    assert!((period - 2.0 * PI / omega).abs() < 1e-4 * period);
    assert!(drift < 1e-14);
    assert!(x1.iter().all(|(_, x)| x.is_finite()));
}

/// A rod2 bar along x, held across, its end pushed at velocity v: a wave front at c t and
/// σ = −ρ c v behind it, carried at the bar's own critical increment L/c.
#[test]
fn velocity_wave_in_a_truss_bar() {
    let (young, rho, area, length, cells) = (1.0e4f64, 1.0, 0.01, 10.0, 200usize);
    let c = (young / rho).sqrt();
    let t = 0.6 * length / c;
    let v = 0.01;
    for kind in [ElementKind::Rod2, ElementKind::Rod3] {
        let mut s = ExplicitDynamic::new(t);
        s.output_interval = Some(t / 2.0);
        // A rod3 is small-displacement only; the rod2 runs corotational (the default).
        s.nonlinear_geometry = kind == ElementKind::Rod2;
        let mut m = model(s);
        let h = length / cells as f64;
        for i in 0..=cells {
            m.set_node(1 + i as u32, [i as f64 * h, 0.0, 0.0]).unwrap();
        }
        for i in 0..cells {
            let mut nodes = vec![1 + i as u32, 2 + i as u32];
            if kind == ElementKind::Rod3 {
                m.set_node(1000 + i as u32, [(i as f64 + 0.5) * h, 0.0, 0.0]).unwrap();
                nodes.push(1000 + i as u32);
            }
            m.set_element(1 + i as u32, kind, nodes).unwrap();
        }
        m.add_set("bar", Set { elements: (1..=cells as u32).collect(), ..Set::default() }).unwrap();
        m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young, poisson: 0.0 }), density: Some(rho), ..Material::default() })
            .unwrap();
        m.add_section(
            "bar",
            Section {
                orientation: None,
                set: "bar".into(),
                material: "Mat".into(),
                kind: SectionKind::Rod { area },
                formulation: Formulation::default(),
            },
        )
        .unwrap();
        let all: Vec<u32> = m.nodes().keys().copied().collect();
        let line = set(&mut m, "line", &all);
        hold(&mut m, "line", line, [None, Some(0.0), Some(0.0), None, None, None]);
        let start = set(&mut m, "start", &[1]);
        hold(&mut m, "push", start, [Some(v * t), None, None, None, None, None]);
        let r = solve(&mut m);
        let f = r.frames.last().unwrap();
        // σ behind the front (x < 0.8 c t) and the front where it falls to half.
        let stress: Vec<(f64, f64)> = (0..cells).map(|i| ((i as f64 + 0.5) * h, axial(&f.stress[i], [1.0, 0.0, 0.0], false))).collect();
        let behind: Vec<f64> = stress.iter().filter(|(x, _)| *x > 0.1 * c * t && *x < 0.8 * c * t).map(|s| s.1).collect();
        let sigma = behind.iter().sum::<f64>() / behind.len() as f64;
        let front =
            stress.windows(2).find(|w| w[0].1 / (-rho * c * v) >= 0.5 && w[1].1 / (-rho * c * v) < 0.5).map(|w| (w[0].0 + w[1].0) / 2.0).unwrap();
        let expected_dt = if kind == ElementKind::Rod2 { 0.9 * h / c } else { r.initial_increment };
        eprintln!(
            "{kind:?} bar wave: σ {sigma:.6} vs −ρcv {:.6}, front {front:.4} vs c t {:.4}, Δt {:.3e}, energy error {:.1e}",
            -rho * c * v,
            c * t,
            r.initial_increment,
            r.energy.relative_error(r.initial_kinetic)
        );
        assert!((sigma + rho * c * v).abs() < 0.01 * rho * c * v);
        assert!((front - c * t).abs() < 3.0 * h);
        assert!((r.initial_increment - expected_dt).abs() < 1e-12 * expected_dt);
        assert!(r.energy.relative_error(r.initial_kinetic) < 0.01, "{:?}", r.energy);
    }
}

#[test]
fn what_explicit_cannot_take_is_named() {
    // A node with a spring and no mass.
    let mut m = model(ExplicitDynamic::new(1.0));
    m.set_node(1, [0.0; 3]).unwrap();
    let n = set(&mut m, "n", &[1]);
    m.add_spring("k", Spring::Spring { first: n.clone(), second: None, direction: Direction::Dof(1), stiffness: 1.0 }).unwrap();
    let e = solve_explicit(&m, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("node 1 has no mass"), "{e}");
    // Rotary inertia, and a rotational spring.
    m.add_reference_point("RP", ReferencePoint { xyz: [0.0; 3], label: 0 }).unwrap();
    m.add_mass("j", PointMass { target: Target::ReferencePoint("RP".into()), mass: None, inertia: Some([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]) }).unwrap();
    let e = solve_explicit(&m, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("rotary inertia"), "{e}");
}

/// A mass whirling on a stiff rod2 tethered to a fixed point (nonlinear geometry): the rod turns
/// with it and keeps its length, carrying the centripetal force m v²/L. A small-displacement rod
/// would let the mass fly off along the tangent.
#[test]
fn whirling_mass_on_a_corotational_rod() {
    let (l, mass, v, axial) = (1.0f64, 1.0, 1.0, 1.0e4);
    let t = std::f64::consts::FRAC_PI_2 * l / v;
    let mut s = step(t, 0.9, "tip", [0.0, v, 0.0]);
    s.bulk_viscosity = BulkViscosity { linear: 0.0, quadratic: 0.0 };
    let mut m = model(s);
    m.set_node(1, [0.0; 3]).unwrap();
    m.set_node(2, [l, 0.0, 0.0]).unwrap();
    m.set_element(1, ElementKind::Rod2, vec![1, 2]).unwrap();
    m.add_set("rod", Set { elements: vec![1], ..Set::default() }).unwrap();
    // A light rod (E A = 1e4, ρ A = 1e-4) with the mass at its end.
    m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: axial, poisson: 0.0 }), density: Some(1e-4), ..Material::default() })
        .unwrap();
    m.add_section(
        "rod",
        Section {
            orientation: None,
            set: "rod".into(),
            material: "Mat".into(),
            kind: SectionKind::Rod { area: 1.0 },
            formulation: Formulation::default(),
        },
    )
    .unwrap();
    let pivot = set(&mut m, "pivot", &[1]);
    m.add_boundary_condition("pivot", BoundaryCondition::Fixed { target: pivot }).unwrap();
    m.set_boundary_condition_active("Dynamic", "pivot", true).unwrap();
    let tip = set(&mut m, "tip", &[2]);
    m.add_mass("m", PointMass { target: tip.clone(), mass: Some(Mass::Isotropic(mass)), inertia: None }).unwrap();
    hold(&mut m, "plane", tip, [None, None, Some(0.0), None, None, None]);
    let r = solve(&mut m);
    let radius = r.history.displacement.iter().map(|d| ((l + d[0][0]).powi(2) + d[0][1].powi(2)).sqrt()).fold(0.0f64, |w, x| w.max((x - l).abs()));
    let end = r.history.displacement.last().unwrap()[0];
    eprintln!(
        "whirling mass: largest change of radius {radius:.3e} (centripetal stretch m v²/EA = {:.1e}), end at ({:.4}, {:.4})",
        mass * v * v / axial,
        l + end[0],
        end[1]
    );
    assert!(radius < 5.0 * mass * v * v / axial, "the rod keeps its length to its elastic stretch");
    // A quarter turn: the mass ends near (0, L).
    assert!((l + end[0]).abs() < 1e-2 && (end[1] - l).abs() < 1e-2, "{end:?}");
}

/// Mass scaling to twice a uniform bar's stable increment quadruples its density: the added mass
/// is three times the bar's, the stable increment doubles with the halved wave speed
/// c = √(E/ρ), and the kinetic energy of the added mass is three quarters of the total.
#[test]
fn mass_scaling_slows_the_waves_as_the_added_density_says() {
    use fieldcraft_explicit::{Limits, solve_explicit};
    use fieldcraft_model::*;
    let build = |scaled: Option<f64>| {
        let mut m = Model::default();
        let n = 40;
        for i in 0..=n {
            for (j, (y, z)) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)].iter().enumerate() {
                m.set_node((4 * i + j + 1) as u32, [i as f64, *y, *z]).unwrap();
            }
        }
        for i in 0..n {
            let a = |j: usize| (4 * i + j + 1) as u32;
            let b = |j: usize| (4 * (i + 1) + j + 1) as u32;
            m.set_element(i as u32 + 1, ElementKind::Hex8, vec![a(0), b(0), b(1), a(1), a(3), b(3), b(2), a(2)]).unwrap();
        }
        m.add_set("all", Set { elements: (1..=n as u32).collect(), ..Set::default() }).unwrap();
        m.add_material("Mat", Material { elastic: Some(Elastic::Isotropic { young: 1.0, poisson: 0.0 }), density: Some(1.0), ..Material::default() })
            .unwrap();
        let formulation = Formulation { integration: Integration::Reduced, ..Formulation::default() };
        m.add_section("bar", Section { orientation: None, set: "all".into(), material: "Mat".into(), kind: SectionKind::Solid {}, formulation })
            .unwrap();
        let mut step = ExplicitDynamic::new(10.0);
        step.initial_velocities = vec![InitialVelocity { set: "all".into(), velocity: [0.01, 0.0, 0.0] }];
        step.mass_scaling = scaled.map(|t| MassScaling { target_increment: t, set: None });
        m.add_step(Step::new("Run", Procedure::ExplicitDynamic(step)), None).unwrap();
        solve_explicit(&m, None, &Limits::default()).unwrap()
    };
    let plain = build(None);
    let scaled = build(Some(2.0 * plain.initial_increment));
    let s = scaled.mass_scaling.as_ref().unwrap();
    eprintln!(
        "mass scaling: {s:?}; increments {} vs {}; Δt {:.4e} vs {:.4e}",
        scaled.increments, plain.increments, scaled.initial_increment, plain.initial_increment
    );
    assert_eq!(s.elements, 40);
    assert!((s.added_mass / plain.total_mass - 3.0).abs() < 1e-9 && (s.fraction - 0.75).abs() < 1e-9, "{s:?}");
    assert_eq!(s.sections.len(), 1);
    assert!((scaled.initial_increment / plain.initial_increment - 2.0).abs() < 1e-9);
    assert!((scaled.total_mass / plain.total_mass - 4.0).abs() < 1e-9);
    assert!((scaled.energy.added_kinetic / scaled.energy.kinetic - 0.75).abs() < 1e-9, "{:?}", scaled.energy);
    assert!(scaled.warnings.iter().any(|w| w.contains("mass scaling")));
    assert!(
        scaled.controlling.first().is_some_and(|c| (0.9 * c.increment / scaled.initial_increment - 1.0).abs() < 1e-6),
        "{:?}",
        scaled.controlling
    );
}
