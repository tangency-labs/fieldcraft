//! Heat transfer with temperature-dependent properties (`fieldcraft_model::curves`), against the
//! Kirchhoff transform: across a slab held at T₀ and T₁, θ(T) = ∫k dT is linear in x, so with
//! k = a + bT the temperature is the root of a T + b T²/2 = θ(x) (H. S. Carslaw, J. C. Jaeger,
//! *Conduction of Heat in Solids*, 2nd ed., Oxford 1959, §1.7). The conductivity is the cited
//! EN 1993-1-2 carbon-steel curve of the material library.

mod common;

use common::*;
use fieldcraft_model::library::{self, Property};
use fieldcraft_model::units::CELSIUS_ZERO;
use fieldcraft_model::*;
use fieldcraft_solver::{HeatResult, Limits, solve_heat, solve_nonlinear, solve_static};

fn hold(g: &mut Grid, name: &str, value: f64, f: impl Fn([f64; 3]) -> bool) {
    g.node_set(name, f);
    g.bc(name, BoundaryCondition::Temperature { target: Target::Set(name.into()), value });
}

fn temperature_at(r: &HeatResult, label: u32) -> f64 {
    r.temperature()[r.nodes.binary_search(&label).unwrap()]
}

/// A slab of `kind` from x = 0 (held at `cold`) to x = `length` (held at `hot`), of `material`.
fn slab(kind: ElementKind, length: f64, material: Material, cold: f64, hot: f64, procedure: Procedure) -> Grid {
    let mut g = grid([8, 1, 1], kind, |p| [length * p[0], 0.01 * p[1], 0.01 * p[2]]);
    g.define(iso(2e11, 0.3), material.density, full(), procedure);
    g.model.replace_material("Mat", material).unwrap();
    hold(&mut g, "cold", cold, |p| near(p[0], 0.0));
    hold(&mut g, "hot", hot, |p| near(p[0], 1.0));
    g
}

fn steel() -> Material {
    library::entry("EN10025-S355").unwrap().material(library::DEFAULT_TEMPERATURE, library::DEFAULT_TEMPERATURE).unwrap()
}

#[test]
fn a_conductivity_falling_with_the_temperature_bows_the_slab_profile_as_kirchhoff_says() {
    let (cold, hot, length) = (CELSIUS_ZERO + 20.0, CELSIUS_ZERO + 700.0, 0.2);
    // EN 1993-1-2: λ = 54 − 3.33e-2 θ (θ in °C), so k = a + b T in kelvin.
    let b = -3.33e-2;
    let a = 54.0 - b * CELSIUS_ZERO;
    let theta = |t: f64| a * t + b * t * t / 2.0;
    let exact = |x: f64| {
        let target = theta(cold) + (theta(hot) - theta(cold)) * x / length;
        // b T²/2 + a T − θ = 0, the root between the end temperatures.
        (-a + (a * a + 2.0 * b * target).sqrt()) / b
    };
    let curve = &library::entry("EN10025-S355").unwrap().curves[&Property::Conductivity];
    let linear = curve.points.0.iter().filter(|p| p[0] < CELSIUS_ZERO + 800.0);
    assert!(
        linear.clone().count() > 5 && linear.clone().all(|&[t, k]| (k - (a + b * t)).abs() < 1e-12 * k),
        "the cited curve is linear below 800 °C"
    );
    for kind in [ElementKind::Hex8, ElementKind::Hex20] {
        let g = slab(kind, length, steel(), cold, hot, Procedure::HeatTransfer { transient: None });
        let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
        assert!(r.warnings.is_empty(), "{kind:?}: {:?}", r.warnings);
        // Linear elements with the element-mean conductivity are exact at the nodes; quadratic ones
        // take one conductivity across a curved profile, off by up to 1.7e-3 (relative, in
        // kelvin) at their midside nodes here.
        let tolerance = if kind == ElementKind::Hex8 { 1e-9 } else { 2e-3 };
        let (mut bow, mut worst): (f64, f64) = (0.0, 0.0);
        for (l, p) in &g.param {
            let (t, want) = (temperature_at(&r, *l), exact(length * p[0]));
            worst = worst.max((t - want).abs() / want);
            assert!((t - want).abs() < tolerance * want, "{kind:?} x = {}: {t} vs {want}", length * p[0]);
            bow = bow.max(cold + (hot - cold) * p[0] - t);
        }
        println!("{kind:?}: largest error {worst:.2e} of the temperature");
        // The dependence acts: the constant-conductivity profile is straight, this one sags below
        // it (the cold side conducts better, so its gradient is gentler).
        println!("{kind:?}: the profile sags {bow:.2} K below the straight line");
        assert!(bow > 25.0, "{kind:?}: {bow}");
        // The heat through the slab is (θ(T₁) − θ(T₀)) / L per area.
        let flow = (theta(hot) - theta(cold)) / length * 0.01 * 0.01;
        let into: f64 = g.nodes_where(|p| near(p[0], 1.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()]).sum();
        assert!((into - flow).abs() < 10.0 * tolerance * flow, "{kind:?}: {into} vs {flow}");
    }
}

#[test]
fn a_steel_bar_radiating_from_its_end_settles_both_nonlinearities_in_one_loop() {
    // Held at T₀ on x = 0 and radiating from x = L: the heat q through the bar is what leaves
    // the end, q = ε σ (T_L⁴ − T_amb⁴), and the Kirchhoff transform carries it across,
    // θ(T₀) − θ(T_L) = q L. Both ends of the balance follow the temperature, so the bar's end
    // temperature is the root of one equation (bisection here), and the nodes between follow
    // θ(x) = θ(T₀) − q x, exact for linear elements as above.
    let (hot, length, eps, ambient) = (CELSIUS_ZERO + 700.0, 0.2, 0.8, CELSIUS_ZERO + 20.0);
    let b = -3.33e-2;
    let a = 54.0 - b * CELSIUS_ZERO;
    let theta = |t: f64| a * t + b * t * t / 2.0;
    let temperature = |target: f64| (-a + (a * a + 2.0 * b * target).sqrt()) / b;
    let out = |t: f64| eps * STEFAN_BOLTZMANN * (t.powi(4) - ambient.powi(4));
    let (mut lo, mut hi) = (ambient, hot);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if out(mid) - (theta(hot) - theta(mid)) / length > 0.0 { hi = mid } else { lo = mid }
    }
    let end = 0.5 * (lo + hi);
    let q = out(end);
    let mut g = grid([8, 1, 1], ElementKind::Hex8, |p| [length * p[0], 0.01 * p[1], 0.01 * p[2]]);
    let material = steel();
    g.define(iso(2e11, 0.3), material.density, full(), Procedure::HeatTransfer { transient: None });
    g.model.replace_material("Mat", material).unwrap();
    hold(&mut g, "hot", hot, |p| near(p[0], 0.0));
    g.face_set("end", |p| near(p[0], 1.0));
    g.load("glow", Load::Radiation { set: "end".into(), emissivity: eps, ambient });
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    for (l, p) in &g.param {
        let (t, want) = (temperature_at(&r, *l), temperature(theta(hot) - q * length * p[0]));
        assert!((t - want).abs() < 1e-9 * want, "x = {}: {t} vs {want}", length * p[0]);
    }
    // Both act: the end sits well below the held end and well above the surroundings, and the
    // heat the radiation takes is what the held end supplies.
    assert!(end < hot - 100.0 && end > ambient + 300.0, "{end}");
    let glow = r.loads.iter().find(|f| f.name == "glow").unwrap();
    let area = 0.01 * 0.01;
    assert!(
        (glow.heat + q * area).abs() < 1e-9 * q * area && (r.reaction_total - q * area).abs() < 1e-9 * q * area,
        "{} {}",
        glow.heat,
        r.reaction_total
    );
}

#[test]
fn constant_tables_reproduce_the_constant_material_exactly() {
    let base = Material {
        density: Some(7850.0),
        thermal: Some(Thermal { conductivity: Some(Directional::Isotropic(45.0)), specific_heat: Some(460.0), ..Thermal::default() }),
        ..Material::default()
    };
    let flat = |v: f64| Some(TemperatureTable(vec![[200.0, v], [900.0, v]]));
    let tabled = Material {
        temperature_dependent: Some(TemperatureDependent {
            conductivity: flat(45.0),
            specific_heat: flat(460.0),
            density: flat(7850.0),
            viscosity: None,
        }),
        ..base.clone()
    };
    let transient = Transient { end_time: 50.0, increment: 5.0, theta: 1.0, initial_temperature: 300.0 };
    for procedure in [Procedure::HeatTransfer { transient: None }, Procedure::HeatTransfer { transient: Some(transient) }] {
        let solve = |m: &Material| {
            solve_heat(&slab(ElementKind::Hex8, 0.1, m.clone(), 300.0, 400.0, procedure.clone()).model, None, &Limits::default()).unwrap()
        };
        let (constant, variable) = (solve(&base), solve(&tabled));
        assert_eq!(constant.frames, variable.frames);
        assert_eq!(constant.reaction, variable.reaction);
        assert!(variable.warnings.is_empty(), "{:?}", variable.warnings);
    }
}

#[test]
fn a_heat_capacity_rising_with_the_temperature_slows_the_heating() {
    // A transient slab heated from one face: with ρc rising with the temperature (400 at 300 K,
    // 800 at 600 K, the line extended past both so the early undershoot of the consistent
    // capacity stays inside) it lags the slab whose ρc stays at its cold value, and leads the
    // one at its hot value.
    let with = |c: Option<TemperatureTable>, c0: f64| Material {
        density: Some(7850.0),
        thermal: Some(Thermal { conductivity: Some(Directional::Isotropic(45.0)), specific_heat: Some(c0), ..Thermal::default() }),
        temperature_dependent: c.map(|c| TemperatureDependent { specific_heat: Some(c), ..TemperatureDependent::default() }),
        ..Material::default()
    };
    let transient = Transient { end_time: 200.0, increment: 2.0, theta: 1.0, initial_temperature: 300.0 };
    let middle = |m: Material| {
        let g = slab(ElementKind::Hex8, 0.1, m, 300.0, 600.0, Procedure::HeatTransfer { transient: Some(transient) });
        let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
        temperature_at(&r, g.nodes_where(|p| near(p[0], 0.75) && near(p[1], 0.0) && near(p[2], 0.0))[0])
    };
    let rising = middle(with(Some(TemperatureTable(vec![[250.0, 400.0 - 200.0 / 3.0], [650.0, 800.0 + 200.0 / 3.0]])), 400.0));
    let (cold, hot) = (middle(with(None, 400.0)), middle(with(None, 800.0)));
    assert!(hot < rising && rising < cold, "{hot} < {rising} < {cold}");
}

#[test]
fn earlier_transient_frames_report_the_flux_of_their_own_conductivity() {
    // Two increments of the steel slab heated from 20 °C to 700 °C at x = L: in frame 1 (not the
    // last, so recovered by `flux_at`) the node at x = L belongs to one element, whose profile
    // is linear in x, so its flux is −k(T̄) (T_L − T_{L−h}) / h exactly, k at the element's mean.
    let (cold, hot, length) = (CELSIUS_ZERO + 20.0, CELSIUS_ZERO + 700.0, 0.2);
    let transient = Transient { end_time: 200.0, increment: 100.0, theta: 1.0, initial_temperature: cold };
    let g = slab(ElementKind::Hex8, length, steel(), cold, hot, Procedure::HeatTransfer { transient: Some(transient) });
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    assert_eq!(r.frames.len(), 3);
    let flux = r.flux_at(&g.model, 1).unwrap();
    let at = |x: f64| {
        let label = g.nodes_where(|p| near(p[0], x) && near(p[1], 0.0) && near(p[2], 0.0))[0];
        let i = r.nodes.binary_search(&label).unwrap();
        (i, r.frames[1][i])
    };
    let ((end, t_end), (_, t_before)) = (at(1.0), at(7.0 / 8.0));
    let k = TemperatureTable(library::entry("EN10025-S355").unwrap().curves[&Property::Conductivity].points.0.clone());
    let conductivity = k.value(0.5 * (t_end + t_before));
    let want = -conductivity * (t_end - t_before) / (length / 8.0);
    assert!((flux[end][0] - want).abs() < 1e-9 * want.abs(), "{} vs {want}", flux[end][0]);
    // The material's constant (at 20 °C) would have overstated it by a third or more.
    assert!(k.value(cold) > 1.3 * conductivity, "{conductivity}");
}

#[test]
fn temperatures_outside_a_table_are_held_at_its_end_and_reported() {
    let g = slab(ElementKind::Hex8, 0.2, steel(), CELSIUS_ZERO + 20.0, CELSIUS_ZERO + 1500.0, Procedure::HeatTransfer { transient: None });
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    // Reported as data in kelvin, for the caller to word in its display system.
    let [range] = r.table_ranges.as_slice() else { panic!("{:?}", r.table_ranges) };
    assert_eq!(
        (range.material.as_str(), range.tables.as_slice()),
        ("Mat", &[("conductivity".to_string(), [CELSIUS_ZERO + 20.0, CELSIUS_ZERO + 1200.0])][..])
    );
    assert!(range.reached[1] > CELSIUS_ZERO + 1200.0 && range.reached[1] < CELSIUS_ZERO + 1500.0, "{range:?}");
    assert!(range.to_string().contains("outside its conductivity (293.15 K to 1473.15 K) table range"), "{range}");
    // Above 1200 °C the conductivity is held at its 27.3 W/(m K).
    assert!(r.temperature().iter().all(|t| t.is_finite()));
}

#[test]
fn density_and_heat_capacity_tables_keep_an_orthotropic_conductivity() {
    // Heat along Y through a plate of conductivity [45, 30, 20]: with only density and specific
    // heat tables the conductivity keeps its axes (k_y, not k_x), steady and transient.
    let orthotropic = |tables: Option<TemperatureDependent>| Material {
        density: Some(7850.0),
        thermal: Some(Thermal { conductivity: Some(Directional::Orthotropic([45.0, 30.0, 20.0])), specific_heat: Some(460.0), ..Thermal::default() }),
        temperature_dependent: tables,
        ..Material::default()
    };
    let flat = |v: f64| Some(TemperatureTable(vec![[200.0, v], [900.0, v]]));
    let tables = TemperatureDependent { specific_heat: flat(460.0), density: flat(7850.0), ..TemperatureDependent::default() };
    let transient = Transient { end_time: 50.0, increment: 5.0, theta: 1.0, initial_temperature: 300.0 };
    for procedure in [Procedure::HeatTransfer { transient: None }, Procedure::HeatTransfer { transient: Some(transient) }] {
        let solve = |m: Material| {
            let mut g = grid([1, 4, 1], ElementKind::Hex8, |p| [0.01 * p[0], 0.1 * p[1], 0.01 * p[2]]);
            g.define(iso(2e11, 0.3), m.density, full(), procedure.clone());
            g.model.replace_material("Mat", m).unwrap();
            hold(&mut g, "cold", 300.0, |p| near(p[1], 0.0));
            hold(&mut g, "hot", 400.0, |p| near(p[1], 1.0));
            solve_heat(&g.model, None, &Limits::default()).unwrap()
        };
        let (constant, variable) = (solve(orthotropic(None)), solve(orthotropic(Some(tables.clone()))));
        assert_eq!(constant.frames, variable.frames);
        assert_eq!(constant.reaction, variable.reaction);
    }
}

#[test]
fn a_transient_that_starts_below_a_table_is_reported_even_if_it_ends_inside() {
    // Initially at 250 K, below the table's 280 K, then held at 300 K and 320 K: the early
    // increments used the held end value, so the run says so although it ends inside.
    let c = TemperatureTable(vec![[280.0, 460.0], [400.0, 500.0]]);
    let m = Material {
        density: Some(7850.0),
        thermal: Some(Thermal { conductivity: Some(Directional::Isotropic(45.0)), specific_heat: Some(460.0), ..Thermal::default() }),
        temperature_dependent: Some(TemperatureDependent { specific_heat: Some(c), ..TemperatureDependent::default() }),
        ..Material::default()
    };
    let transient = Transient { end_time: 2000.0, increment: 20.0, theta: 1.0, initial_temperature: 250.0 };
    let g = slab(ElementKind::Hex8, 0.05, m, 300.0, 320.0, Procedure::HeatTransfer { transient: Some(transient) });
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    assert!(r.temperature().iter().all(|t| *t > 280.0), "the step ends inside the table");
    assert!(r.table_ranges.iter().any(|t| t.tables.iter().any(|(what, _)| what == "specific heat") && t.reached[0] == 250.0), "{:?}", r.table_ranges);
}

#[test]
fn a_coupled_step_carries_the_heat_solution_s_warnings() {
    // The slab of `temperatures_outside_a_table_are_held_at_its_end_and_reported`, as a coupled
    // temperature–displacement step: its thermal strain takes temperatures outside the table, and
    // the static result says so as the heat result does.
    // A nonlinear static step with the same thermal conditions takes the same heat solution.
    for (procedure, nonlinear) in
        [(Procedure::CoupledTemperatureDisplacement {}, false), (Procedure::NonlinearStatic(NonlinearStatic::default()), true)]
    {
        let mut g = slab(ElementKind::Hex8, 0.2, steel(), CELSIUS_ZERO + 20.0, CELSIUS_ZERO + 1500.0, procedure);
        g.symmetry("x0", Symmetry::X, |p| near(p[0], 0.0));
        g.symmetry("y0", Symmetry::Y, |p| near(p[1], 0.0));
        g.symmetry("z0", Symmetry::Z, |p| near(p[2], 0.0));
        let r =
            if nonlinear { solve_nonlinear(&g.model, None, &Limits::default()) } else { solve_static(&g.model, None, &Limits::default()) }.unwrap();
        let [range] = r.table_ranges.as_slice() else { panic!("{nonlinear}: {:?}", r.table_ranges) };
        assert_eq!(range.tables[0].0, "conductivity", "{nonlinear}");
    }
}
