//! Analytic verification of heat transfer and thermal stress (docs/fem.md, docs/elements.md):
//! steady conduction through a slab (prescribed temperatures, a heat source, convection),
//! radial conduction through a thick cylinder, transient conduction into a semi-infinite solid
//! (erfc), free and constrained thermal expansion, a bimetal strip (Timoshenko 1925) and the
//! coupled temperature–displacement step.

mod common;

use common::*;
use fieldcraft_model::*;
use fieldcraft_solver::{Limits, solve_heat, solve_static};

const K: f64 = 45.0;

/// Replace the grid's material with `elastic` plus thermal properties.
fn thermal(g: &mut Grid, elastic: Elastic, density: Option<f64>, t: Thermal) {
    g.model.replace_material("Mat", Material { elastic: Some(elastic), density, thermal: Some(t), ..Material::default() }).unwrap();
}

fn conductor(k: f64) -> Thermal {
    Thermal { conductivity: Some(Directional::Isotropic(k)), ..Thermal::default() }
}

fn hold(g: &mut Grid, name: &str, value: f64, f: impl Fn([f64; 3]) -> bool) {
    g.node_set(name, f);
    g.bc(name, BoundaryCondition::Temperature { target: Target::Set(name.into()), value });
}

fn temperature_at(r: &fieldcraft_solver::HeatResult, label: u32) -> f64 {
    r.temperature()[r.nodes.binary_search(&label).unwrap()]
}

/// Abramowitz & Stegun 7.1.26 is too coarse here; erfc by its continued fraction for large x and
/// the Taylor series for small x, both to ~1e-15.
fn erfc(x: f64) -> f64 {
    if x < 2.0 {
        // erf x = 2/√π Σ (−1)ⁿ x^(2n+1) / (n! (2n + 1)).
        let mut sum = 0.0;
        let mut term = x;
        for n in 0..200 {
            let t = term / (2 * n + 1) as f64;
            sum += t;
            term *= -x * x / (n + 1) as f64;
            if t.abs() < 1e-17 {
                break;
            }
        }
        1.0 - 2.0 / std::f64::consts::PI.sqrt() * sum
    } else {
        // Lentz evaluation of erfc x = e^(−x²)/√π · 1/(x + ½/(x + 1/(x + 3/2/(x + …)))).
        let mut f = x;
        for n in (1..200).rev() {
            f = x + n as f64 / 2.0 / f;
        }
        (-x * x).exp() / std::f64::consts::PI.sqrt() / f
    }
}

#[test]
fn erfc_matches_reference_values() {
    for (x, want) in [(0.0, 1.0), (0.5, 0.479_500_122_186_953_5), (1.0, 0.157_299_207_050_285_1), (2.5, 4.069_520_174_449_59e-4)] {
        assert!((erfc(x) - want).abs() < 1e-13 * want.max(1e-3), "{x}: {} vs {want}", erfc(x));
    }
}

#[test]
fn slab_steady_conduction_is_linear_with_the_fourier_flux() {
    // A 4 × 1 × 1 slab, 20 at x = 0 and 100 at x = 4: T linear, q_x = −k ΔT / L everywhere.
    for (kind, f, name) in variants() {
        let mut g = grid([4, 1, 1], kind, |p| [4.0 * p[0], p[1], p[2]]);
        g.define(iso(2e5, 0.3), None, f, Procedure::HeatTransfer { transient: None });
        thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
        hold(&mut g, "cold", 20.0, |p| near(p[0], 0.0));
        hold(&mut g, "hot", 100.0, |p| near(p[0], 1.0));
        let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
        for (l, p) in &g.param {
            let t = temperature_at(&r, *l);
            assert!((t - (20.0 + 80.0 * p[0])).abs() < 1e-9, "{name} node {l}: {t}");
        }
        for e in &r.elements {
            for p in &e.points {
                assert!((p.flux[0] + K * 20.0).abs() < 1e-8 && p.flux[1].abs() < 1e-8 && p.flux[2].abs() < 1e-8, "{name}: {:?}", p.flux);
            }
        }
        // The heat entering at the hot face leaves at the cold one: q A = k ΔT / L · 1.
        let hot: f64 = g.nodes_where(|p| near(p[0], 1.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()]).sum();
        assert!((hot - K * 20.0).abs() < 1e-8, "{name}: {hot}");
        assert!(r.reaction_total.abs() < 1e-8 && r.warnings.iter().all(|w| !w.contains("differ")), "{name}");
    }
}

#[test]
fn slab_with_a_heat_source_and_convection() {
    // Uniform source Q in a slab of length L held at 0 at x = 0, convecting to T∞ at x = L:
    // T(x) = −Q x²/(2k) + a x with k a − Q L = h (T∞ − T(L)) from the energy balance at x = L,
    // a = (Q L + h T∞ + h Q L²/(2k)) / (k + h L). Quadratic in x: exact at the nodes of the
    // quadratic elements.
    let (q, h, sink, l) = (300.0, 25.0, 60.0, 2.0);
    let a = (q * l + h * sink + h * q * l * l / (2.0 * K)) / (K + h * l);
    let exact = |x: f64| -q * x * x / (2.0 * K) + a * x;
    for (kind, f, name) in variants().into_iter().filter(|v| matches!(v.0, ElementKind::Hex20 | ElementKind::Tet10)) {
        let mut g = grid([4, 1, 1], kind, |p| [l * p[0], 0.5 * p[1], 0.5 * p[2]]);
        g.define(iso(2e5, 0.3), None, f, Procedure::HeatTransfer { transient: None });
        thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
        hold(&mut g, "base", 0.0, |p| near(p[0], 0.0));
        g.face_set("end", |p| near(p[0], 1.0));
        g.load("film", Load::Convection { set: "end".into(), film: h, sink });
        g.load("source", Load::HeatSource { set: "all".into(), power: q });
        let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
        for (lb, p) in &g.param {
            let t = temperature_at(&r, *lb);
            assert!((t - exact(l * p[0])).abs() < 1e-9 * exact(l).abs(), "{name} x = {}: {t} vs {}", l * p[0], exact(l * p[0]));
        }
        // Steady balance: what the source and the film supply leaves through the held face,
        // (Q L + h (T∞ − T(L))) A.
        let supplied = (q * l + h * (sink - exact(l))) * 0.25;
        assert!((r.reaction_total + supplied).abs() < 1e-9 * supplied, "{name}: {} vs {supplied}", -r.reaction_total);
        assert!((r.applied_total - supplied).abs() < 1e-9 * supplied, "{name}");
    }
    // A flux on a face instead of convection: T = −Q x²/(2k) + (Q L + q″) x / k.
    let qs = 700.0;
    let mut g = grid([4, 1, 1], ElementKind::Hex20, |p| [l * p[0], 0.5 * p[1], 0.5 * p[2]]);
    g.define(iso(2e5, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
    hold(&mut g, "base", 0.0, |p| near(p[0], 0.0));
    g.face_set("end", |p| near(p[0], 1.0));
    g.load("flux", Load::SurfaceFlux { set: "end".into(), flux: qs });
    g.load("source", Load::HeatSource { set: "all".into(), power: q });
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    for (lb, p) in &g.param {
        let x = l * p[0];
        let want = -q * x * x / (2.0 * K) + (q * l + qs) * x / K;
        assert!((temperature_at(&r, *lb) - want).abs() < 1e-9 * want.abs().max(1.0));
    }
}

#[test]
fn thick_cylinder_radial_conduction() {
    // A quarter of a hollow cylinder r = 1..2 held at 100 inside and 0 outside:
    // T = 100 (1 − ln r / ln 2). The cut faces are insulated by symmetry.
    for (kind, n, tol) in [(ElementKind::Hex20, [6, 6, 1], 2e-4), (ElementKind::Hex8, [16, 16, 1], 2e-3)] {
        let mut g = grid(n, kind, |p| {
            let (r, th) = (1.0 + p[0], std::f64::consts::FRAC_PI_2 * p[1]);
            [r * th.cos(), r * th.sin(), 0.3 * p[2]]
        });
        g.define(iso(2e5, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
        thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
        hold(&mut g, "inner", 100.0, |p| near(p[0], 0.0));
        hold(&mut g, "outer", 0.0, |p| near(p[0], 1.0));
        let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
        let mut worst: f64 = 0.0;
        for (l, p) in &g.param {
            let rad = 1.0 + p[0];
            worst = worst.max((temperature_at(&r, *l) - 100.0 * (1.0 - rad.ln() / 2f64.ln())).abs());
        }
        assert!(worst < tol * 100.0, "{kind:?}: {worst}");
        // The heat flow through the quarter: Q = (π/2) k L ΔT / ln 2.
        let q_exact = std::f64::consts::FRAC_PI_2 * K * 0.3 * 100.0 / 2f64.ln();
        let q_in: f64 = g.nodes_where(|p| near(p[0], 0.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()]).sum::<f64>().abs();
        assert!((q_in - q_exact).abs() < 10.0 * tol * q_exact, "{kind:?}: {q_in} vs {q_exact}");
    }
}

#[test]
fn transient_conduction_into_a_semi_infinite_solid() {
    // A long bar at 0 whose end x = 0 is held at 1 from t = 0: T = erfc(x / (2√(κt))), κ = k/(ρc),
    // until the front reaches the far end (here 2√(κt) ≤ L/4).
    let (rho, c, length): (f64, f64, f64) = (7.8, 0.46, 10.0);
    let kappa = K / (rho * c);
    let end = (length / 8.0).powi(2) / kappa;
    for (theta, increments) in [(0.5, 200usize), (1.0, 800)] {
        let mut g = grid([100, 1, 1], ElementKind::Hex8, |p| [length * p[0], 0.1 * p[1], 0.1 * p[2]]);
        let tr = Transient { end_time: end, increment: end / increments as f64, theta, initial_temperature: 0.0 };
        g.define(iso(2e5, 0.3), Some(rho), full(), Procedure::HeatTransfer { transient: Some(tr) });
        thermal(&mut g, iso(2e5, 0.3), Some(rho), Thermal { specific_heat: Some(c), ..conductor(K) });
        hold(&mut g, "end", 1.0, |p| near(p[0], 0.0));
        let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
        assert_eq!(*r.times.last().unwrap(), end);
        let mut worst: f64 = 0.0;
        for (l, p) in &g.param {
            let x = length * p[0];
            worst = worst.max((temperature_at(&r, *l) - erfc(x / (2.0 * (kappa * end).sqrt()))).abs());
        }
        assert!(worst < 5e-3, "θ = {theta}: {worst}");
        // Frames start at the initial temperature.
        assert!(r.frames[0].iter().all(|t| *t == 0.0) && r.times[0] == 0.0);
    }
}

/// A 4 × 1 × 1 bar of `kind` heated by `dt` through a temperature load, held as `hold` says.
fn heated_bar(kind: ElementKind, f: Formulation, alpha: Directional, dt: f64, poisson: f64) -> Grid {
    let mut g = grid([4, 1, 1], kind, |p| [4.0 * p[0], p[1], p[2]]);
    g.define(iso(2e5, poisson), None, f, Procedure::Static {});
    thermal(&mut g, iso(2e5, poisson), None, Thermal { expansion: Some(alpha), reference_temperature: 20.0, ..Thermal::default() });
    g.node_set("everything", |_| true);
    g.load("heat", Load::Temperature { target: Target::Set("everything".into()), value: 20.0 + dt });
    g
}

#[test]
fn free_thermal_expansion_is_stress_free() {
    // Held only against rigid motion (3-2-1 by symmetry planes through a corner): ΔL = α L ΔT
    // along each axis, no stress.
    let alpha = [1.2e-5, 2.0e-5, 0.8e-5];
    for (kind, f, name) in variants() {
        let mut g = heated_bar(kind, f, Directional::Orthotropic(alpha), 80.0, 0.3);
        g.symmetry("x0", Symmetry::X, |p| near(p[0], 0.0));
        g.symmetry("y0", Symmetry::Y, |p| near(p[1], 0.0));
        g.symmetry("z0", Symmetry::Z, |p| near(p[2], 0.0));
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        for (i, l) in r.nodes.iter().enumerate() {
            let x = g.model.node(*l).unwrap().xyz;
            for k in 0..3 {
                assert!((r.displacement[i][k] - alpha[k] * 80.0 * x[k]).abs() < 1e-12, "{name} node {l}");
            }
            assert!(r.stress[i].iter().all(|s| s.abs() < 1e-6), "{name}: {:?}", r.stress[i]);
        }
        assert!(r.reaction.iter().flatten().all(|v| v.abs() < 1e-6), "{name}");
    }
}

#[test]
fn constrained_bar_carries_minus_e_alpha_dt() {
    // Both ends held axially, the sides free: σ_xx = −E α ΔT, the others zero, for any ν.
    for (kind, f, name) in variants() {
        let mut g = heated_bar(kind, f, Directional::Isotropic(1.2e-5), 50.0, 0.3);
        g.symmetry("x0", Symmetry::X, |p| near(p[0], 0.0));
        g.symmetry("x1", Symmetry::X, |p| near(p[0], 1.0));
        g.symmetry("y0", Symmetry::Y, |p| near(p[1], 0.0));
        g.symmetry("z0", Symmetry::Z, |p| near(p[2], 0.0));
        let r = solve_static(&g.model, None, &Limits::default()).unwrap();
        let want = -2e5 * 1.2e-5 * 50.0;
        for s in &r.stress {
            assert!((s[0] - want).abs() < 1e-8 * want.abs(), "{name}: {s:?}");
            assert!(s[1..].iter().all(|v| v.abs() < 1e-8 * want.abs()), "{name}: {s:?}");
        }
        // The ends push back with σ A.
        let rx: f64 = g.nodes_where(|p| near(p[0], 1.0)).iter().map(|l| r.reaction[r.nodes.binary_search(l).unwrap()][0]).sum();
        assert!((rx - want).abs() < 1e-8 * want.abs(), "{name}: {rx}");
    }
}

#[test]
fn bimetal_strip_bends_as_timoshenko_says() {
    // Two equal layers (thickness t each, same E, ν = 0) with α₁, α₂, heated by ΔT, clamped at
    // x = 0: curvature κ = 3 Δα ΔT / (2 h) with h = 2t (Timoshenko 1925, m = n = 1), tip rise
    // κ L² / 2.
    let (len, t, dt, a1, a2) = (20.0, 0.5, 100.0, 1.0e-5, 2.0e-5);
    let mut g = grid([20, 1, 4], ElementKind::Hex20, |p| [len * p[0], 1.0 * p[1], 2.0 * t * p[2]]);
    let all: Vec<u32> = g.model.elements().keys().copied().collect();
    let (bottom, top): (Vec<u32>, Vec<u32>) =
        all.into_iter().partition(|e| g.model.element(*e).unwrap().nodes.iter().all(|n| g.param[n][2] <= 0.5 + 1e-12));
    g.model.add_set("bottom", Set { elements: bottom, ..Set::default() }).unwrap();
    g.model.add_set("top", Set { elements: top, ..Set::default() }).unwrap();
    for (name, alpha) in [("Low", a1), ("High", a2)] {
        g.model
            .add_material(
                name,
                Material {
                    elastic: Some(iso(2e5, 0.0)),
                    density: None,
                    thermal: Some(Thermal { expansion: Some(Directional::Isotropic(alpha)), reference_temperature: 0.0, ..Thermal::default() }),
                    ..Material::default()
                },
            )
            .unwrap();
    }
    g.model.add_section("low", Section { set: "bottom".into(), material: "Low".into(), kind: SectionKind::Solid {}, formulation: full() }).unwrap();
    g.model.add_section("high", Section { set: "top".into(), material: "High".into(), kind: SectionKind::Solid {}, formulation: full() }).unwrap();
    g.model.add_step(Step::new("Step", Procedure::Static {}), None).unwrap();
    g.fixed("root", |p| near(p[0], 0.0));
    g.node_set("everything", |_| true);
    g.load("heat", Load::Temperature { target: Target::Set("everything".into()), value: dt });
    let r = solve_static(&g.model, None, &Limits::default()).unwrap();
    let kappa = 1.5 * (a2 - a1) * dt / (2.0 * t);
    let want = -kappa * len * len / 2.0;
    // The tip's mid-height nodes (the high-α layer on top lengthens, so the strip curls down).
    let tip: Vec<f64> = g
        .param
        .iter()
        .filter(|(_, p)| near(p[0], 1.0) && near(p[2], 0.5))
        .map(|(l, _)| r.displacement[r.nodes.binary_search(l).unwrap()][2])
        .collect();
    assert!(!tip.is_empty());
    for w in tip {
        assert!((w - want).abs() < 0.02 * want.abs(), "{w} vs {want}");
    }
}

#[test]
fn coupled_step_matches_heat_then_static() {
    // A slab held at 20 and 100 at its ends, free to expand (3-2-1): the coupled step gives the
    // heat step's temperatures and the static step's displacements for that field.
    let (alpha, tref) = (1.2e-5, 20.0);
    let build = |procedure: Procedure| {
        let mut g = grid([4, 1, 1], ElementKind::Hex20, |p| [4.0 * p[0], p[1], p[2]]);
        g.define(iso(2e5, 0.3), None, full(), procedure);
        let t = Thermal { expansion: Some(Directional::Isotropic(alpha)), reference_temperature: tref, ..conductor(K) };
        thermal(&mut g, iso(2e5, 0.3), None, t);
        g
    };
    let mut c = build(Procedure::CoupledTemperatureDisplacement {});
    hold(&mut c, "cold", 20.0, |p| near(p[0], 0.0));
    hold(&mut c, "hot", 100.0, |p| near(p[0], 1.0));
    c.symmetry("x0", Symmetry::X, |p| near(p[0], 0.0));
    c.symmetry("y0", Symmetry::Y, |p| near(p[1], 0.0));
    c.symmetry("z0", Symmetry::Z, |p| near(p[2], 0.0));
    let r = solve_static(&c.model, None, &Limits::default()).unwrap();
    let temps = r.temperature.clone().unwrap();
    for (i, l) in r.nodes.iter().enumerate() {
        let x = c.model.node(*l).unwrap().xyz;
        assert!((temps[i] - (20.0 + 20.0 * x[0])).abs() < 1e-9);
    }
    // The heat step alone gives the same temperatures.
    let mut h = build(Procedure::HeatTransfer { transient: None });
    hold(&mut h, "cold", 20.0, |p| near(p[0], 0.0));
    hold(&mut h, "hot", 100.0, |p| near(p[0], 1.0));
    let hr = solve_heat(&h.model, None, &Limits::default()).unwrap();
    assert_eq!(hr.temperature(), temps.as_slice());
    // A static step with the same field as temperature loads (one per plane x = const) gives the
    // same displacements and stresses.
    let mut s = build(Procedure::Static {});
    s.symmetry("x0", Symmetry::X, |p| near(p[0], 0.0));
    s.symmetry("y0", Symmetry::Y, |p| near(p[1], 0.0));
    s.symmetry("z0", Symmetry::Z, |p| near(p[2], 0.0));
    for k in 0..=8 {
        let name = format!("plane{k}");
        s.node_set(&name, |p| near(p[0], k as f64 / 8.0));
        s.load(&name, Load::Temperature { target: Target::Set(name.clone()), value: 20.0 + 20.0 * 4.0 * k as f64 / 8.0 });
    }
    let sr = solve_static(&s.model, None, &Limits::default()).unwrap();
    let scale = sr.displacement.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(scale > 0.0);
    for (a, b) in r.displacement.iter().flatten().zip(sr.displacement.iter().flatten()) {
        assert!((a - b).abs() < 1e-12 * scale, "{a} vs {b}");
    }
    for (a, b) in r.stress.iter().flatten().zip(sr.stress.iter().flatten()) {
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
    }
    assert!(sr.temperature.is_none());
}

#[test]
fn heat_steps_refuse_what_they_cannot_solve() {
    // Nothing holds the temperature level.
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(2e5, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
    g.face_set("end", |p| near(p[0], 1.0));
    g.load("flux", Load::SurfaceFlux { set: "end".into(), flux: 10.0 });
    let e = solve_heat(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("undetermined"), "{e}");
    // No conductivity: the model check says so.
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(2e5, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    hold(&mut g, "end", 1.0, |p| near(p[0], 0.0));
    let e = solve_heat(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("conductivity"), "{e}");
    // Transient without a specific heat.
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(2e5, 0.3), Some(1.0), full(), Procedure::HeatTransfer { transient: Some(Transient::new(1.0, 0.1)) });
    thermal(&mut g, iso(2e5, 0.3), Some(1.0), conductor(K));
    hold(&mut g, "end", 1.0, |p| near(p[0], 0.0));
    let e = solve_heat(&g.model, None, &Limits::default()).unwrap_err().to_string();
    assert!(e.contains("specific heat"), "{e}");
    // Conflicting prescribed temperatures.
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(2e5, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
    hold(&mut g, "a", 1.0, |p| near(p[0], 0.0));
    hold(&mut g, "b", 2.0, |p| near(p[0], 0.0));
    assert!(solve_heat(&g.model, None, &Limits::default()).is_err());
    // A hard cap on the unknowns.
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| p);
    g.define(iso(2e5, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    thermal(&mut g, iso(2e5, 0.3), None, conductor(K));
    hold(&mut g, "a", 1.0, |p| near(p[0], 0.0));
    let tight = Limits { max_dofs: 3, ..Limits::default() };
    assert!(solve_heat(&g.model, None, &tight).unwrap_err().to_string().contains("limit"));
    // A static step is not a heat step.
    let g = heated_bar(ElementKind::Hex8, full(), Directional::Isotropic(1e-5), 10.0, 0.3);
    assert!(solve_heat(&g.model, None, &Limits::default()).is_err());
}

/// One rod3 or two rod2 along X over [0, `len`], heated by `dt`, fixed at x = 0, the other
/// nodes held in Y and Z, and the end held along X too when `held`.
fn heated_rod(kind: ElementKind, len: f64, alpha: [f64; 3], dt: f64, held: bool) -> fieldcraft_solver::StaticResult {
    let mut m = Model::default();
    for (label, x) in [(1, 0.0), (2, len), (3, 0.5 * len)] {
        m.set_node(label, [x, 0.0, 0.0]).unwrap();
    }
    // One rod3, or two rod2 in a row.
    let elements = if kind == ElementKind::Rod3 {
        m.set_element(1, kind, vec![1, 2, 3]).unwrap();
        vec![1]
    } else {
        m.set_element(1, kind, vec![1, 3]).unwrap();
        m.set_element(2, kind, vec![3, 2]).unwrap();
        vec![1, 2]
    };
    m.add_set("rod", Set { elements, nodes: vec![1, 2, 3], ..Set::default() }).unwrap();
    m.add_set("a", Set { nodes: vec![1], ..Set::default() }).unwrap();
    m.add_set("b", Set { nodes: vec![2], ..Set::default() }).unwrap();
    m.add_set("c", Set { nodes: vec![3], ..Set::default() }).unwrap();
    let thermal = Thermal { expansion: Some(Directional::Orthotropic(alpha)), reference_temperature: 0.0, ..Thermal::default() };
    m.add_material("Mat", Material { elastic: Some(iso(2e5, 0.3)), thermal: Some(thermal), ..Material::default() }).unwrap();
    m.add_section("rod", Section { set: "rod".into(), material: "Mat".into(), kind: SectionKind::Rod { area: 4.0 }, formulation: full() }).unwrap();
    let mut step = Step::new("Step", Procedure::Static {});
    let lateral = [None, Some(0.0), Some(0.0)];
    for (name, v) in [("a", [Some(0.0); 3]), ("b", if held { [Some(0.0); 3] } else { lateral }), ("c", lateral)] {
        let values = [v[0], v[1], v[2], None, None, None];
        m.add_boundary_condition(name, BoundaryCondition::Displacement { target: Target::Set(name.into()), values }).unwrap();
        step.boundary_conditions.insert(name.into());
    }
    m.add_load("warm", Load::Temperature { target: Target::Set("rod".into()), value: dt }).unwrap();
    step.loads.insert("warm".into());
    m.add_step(step, None).unwrap();
    solve_static(&m, None, &Limits::default()).unwrap()
}

#[test]
fn trusses_take_thermal_strain_along_their_axis() {
    // Along X the member's thermal strain is αₓ ΔT (the projection Σ αₖ tₖ² is checked on the
    // kernel): free it lengthens by αₓ ΔT L stress-free; held at both ends it carries −E αₓ ΔT.
    let alpha = [1.0e-5, 2.0e-5, 3.0e-5];
    let (e, area, dt, len) = (2e5, 4.0, 50.0, 3.0);
    let sigma = -e * alpha[0] * dt;
    for kind in [ElementKind::Rod2, ElementKind::Rod3] {
        let r = heated_rod(kind, len, alpha, dt, true);
        for p in r.elements.iter().flat_map(|el| &el.points) {
            assert!((p.stress[0] - sigma).abs() < 1e-9 * sigma.abs(), "{kind:?}: {:?}", p.stress);
        }
        let reaction = r.reaction[r.nodes.binary_search(&2).unwrap()][0];
        assert!((reaction - sigma * area).abs() < 1e-9 * (sigma * area).abs(), "{kind:?}: {reaction}");
        let r = heated_rod(kind, len, alpha, dt, false);
        let u = r.displacement[r.nodes.binary_search(&2).unwrap()][0];
        assert!((u - alpha[0] * dt * len).abs() < 1e-12, "{kind:?}: {u}");
        assert!(r.elements.iter().flat_map(|el| &el.points).all(|p| p.stress[0].abs() < 1e-9));
    }
}

#[test]
fn heat_factor_cap_reports_value_bytes_before_allocation() {
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| [2.0 * p[0], p[1], p[2]]);
    g.define(iso(1000.0, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    thermal(&mut g, iso(1000.0, 0.3), None, conductor(K));
    hold(&mut g, "cold", 0.0, |p| near(p[0], 0.0));
    let before = g.model.clone();
    let tight = Limits { max_factor_entries: 1, ..Limits::default() };
    let e = solve_heat(&g.model, None, &tight).unwrap_err().to_string();
    assert!(e.contains("bytes of factor values alone") && e.contains("limit is 1 entries"), "{e}");
    assert_eq!(g.model, before);
    assert!(solve_heat(&g.model, None, &Limits::default()).is_ok());
}
