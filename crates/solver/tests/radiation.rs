//! Radiation to an ambient against closed forms (docs/manual/structural/heat-transfer.md): a
//! plate heated on one face that loses its heat from the other by radiation alone (to warm
//! surroundings and to 0 K) and by convection and radiation together (the steady energy balance
//! of the cooled face), and a lumped body cooling by radiation alone against the integrated
//! T(t), second order in time under Crank–Nicolson; and the bounds of the Newton iterations.

mod common;

use common::*;
use fieldcraft_model::*;
use fieldcraft_solver::heat::MAX_RADIATION_ITERATIONS;
use fieldcraft_solver::{Error, HeatResult, Limits, solve_heat};

const K: f64 = 15.0;

fn conductor(g: &mut Grid, density: Option<f64>, specific_heat: Option<f64>) {
    let thermal = Thermal { conductivity: Some(Directional::Isotropic(K)), specific_heat, ..Thermal::default() };
    g.model.replace_material("Mat", Material { elastic: Some(iso(2e11, 0.3)), density, thermal: Some(thermal), ..Material::default() }).unwrap();
}

fn temperature_at(r: &HeatResult, label: u32) -> f64 {
    r.temperature()[r.nodes.binary_search(&label).unwrap()]
}

/// The root of the increasing `f` between `lo` and `hi`, by bisection to round-off.
fn root(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 { hi = mid } else { lo = mid }
    }
    0.5 * (lo + hi)
}

/// A plate of thickness `l` (1 × 1 across) heated by `q` on x = 0, cooled on x = l by radiation
/// (ε, T_amb) and, when `film` is given, convection (h, T∞); the other faces insulated.
fn plate(kind: ElementKind, f: Formulation, l: f64, q: f64, emissivity: f64, ambient: f64, film: Option<(f64, f64)>) -> (Grid, HeatResult) {
    let mut g = grid([4, 1, 1], kind, |p| [l * p[0], p[1], p[2]]);
    g.define(iso(2e11, 0.3), None, f, Procedure::HeatTransfer { transient: None });
    conductor(&mut g, None, None);
    g.face_set("heated", |p| near(p[0], 0.0));
    g.face_set("cooled", |p| near(p[0], 1.0));
    g.load("heater", Load::SurfaceFlux { set: "heated".into(), flux: q });
    g.load("glow", Load::Radiation { set: "cooled".into(), emissivity, ambient });
    if let Some((h, sink)) = film {
        g.load("air", Load::Convection { set: "cooled".into(), film: h, sink });
    }
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    (g, r)
}

#[test]
fn plate_cooled_by_radiation_alone_balances_its_heat() {
    // Steady: all of q leaves the cooled face, q = ε σ (T_L⁴ − T_amb⁴), and conduction carries
    // it across, T = T_L + q (L − x) / k: linear, so exact at the nodes of every element.
    let (l, q, eps) = (0.02, 2000.0, 0.8);
    for (kind, f, name) in variants() {
        for ambient in [293.15, 0.0] {
            let (g, r) = plate(kind, f, l, q, eps, ambient, None);
            let tl = (q / (eps * STEFAN_BOLTZMANN) + ambient.powi(4)).powf(0.25);
            for (label, p) in &g.param {
                let want = tl + q * l * (1.0 - p[0]) / K;
                let t = temperature_at(&r, *label);
                assert!((t - want).abs() < 1e-9 * want, "{name}, T_amb {ambient}, x {}: {t} vs {want}", l * p[0]);
            }
            // The heat balance: the flux brings q A in, the radiation takes it out.
            let glow = r.loads.iter().find(|f| f.name == "glow").unwrap();
            assert_eq!(glow.kind, "radiation");
            assert!((glow.heat + q).abs() < 1e-9 * q, "{name}: radiated {}", glow.heat);
            assert!(r.applied_total.abs() < 1e-9 * q && r.reaction_total == 0.0, "{name}: {}", r.applied_total);
        }
    }
}

#[test]
fn plate_with_convection_and_radiation_matches_the_hand_calculation() {
    // The cooled face's balance q = h (T_L − T∞) + ε σ (T_L⁴ − T_amb⁴), solved by hand here by
    // bisection: in still air radiation carries a share comparable to convection.
    let (l, q, eps, ambient, h, sink): (f64, f64, f64, f64, f64, f64) = (0.01, 1500.0, 0.9, 300.0, 10.0, 295.0);
    let balance = |t: f64| h * (t - sink) + eps * STEFAN_BOLTZMANN * (t.powi(4) - ambient.powi(4)) - q;
    let tl = root(balance, ambient, 2000.0);
    let radiated = eps * STEFAN_BOLTZMANN * (tl.powi(4) - ambient.powi(4));
    assert!(radiated > 0.3 * q && radiated < 0.7 * q, "radiation's share {radiated} of {q}");
    for (kind, f, name) in variants() {
        let (g, r) = plate(kind, f, l, q, eps, ambient, Some((h, sink)));
        for (label, p) in &g.param {
            let want = tl + q * l * (1.0 - p[0]) / K;
            let t = temperature_at(&r, *label);
            assert!((t - want).abs() < 1e-9 * want, "{name}: {t} vs {want}");
        }
        let heat = |n: &str| r.loads.iter().find(|f| f.name == n).unwrap().heat;
        assert!((heat("glow") + radiated).abs() < 1e-8 * q, "{name}: {} vs {radiated}", heat("glow"));
        assert!((heat("air") + q - radiated).abs() < 1e-8 * q, "{name}");
        assert!(r.applied_total.abs() < 1e-9 * q, "{name}");
    }
}

/// ∫ dT / (T⁴ − a⁴) = ln((T − a)/(T + a)) / (4a³) − atan(T/a) / (2a³), for T > a > 0.
fn antiderivative(t: f64, a: f64) -> f64 {
    ((t - a) / (t + a)).ln() / (4.0 * a.powi(3)) - (t / a).atan() / (2.0 * a.powi(3))
}

/// A cube of side `s` whose six faces radiate, from `t0`: one hex8 keeps a uniform temperature
/// (each node's capacity and radiating area are the same eighth), so it is the lumped body
/// ρ c V dT/dt = −ε σ A (T⁴ − T_amb⁴) whatever its conductivity.
fn cooling_cube(s: f64, eps: f64, ambient: f64, t0: f64, end: f64, increments: usize, theta: f64) -> HeatResult {
    cooling_cube_with(s, eps, ambient, t0, end, increments, theta, &Limits::default())
}

/// [`cooling_cube`] under `limits`.
#[allow(clippy::too_many_arguments)]
fn cooling_cube_with(s: f64, eps: f64, ambient: f64, t0: f64, end: f64, increments: usize, theta: f64, limits: &Limits) -> HeatResult {
    let mut g = grid([1, 1, 1], ElementKind::Hex8, |p| p.map(|c| s * c));
    let tr = Transient { end_time: end, increment: end / increments as f64, theta, initial_temperature: t0 };
    g.define(iso(2e11, 0.3), Some(7800.0), full(), Procedure::HeatTransfer { transient: Some(tr) });
    conductor(&mut g, Some(7800.0), Some(460.0));
    let skin: Vec<ElementFace> = (0..3).flat_map(|a| [0.0, 1.0].map(|v| g.faces_where(|p| near(p[a], v)))).flatten().collect();
    g.model.add_set("skin", Set { faces: skin, ..Set::default() }).unwrap();
    g.load("glow", Load::Radiation { set: "skin".into(), emissivity: eps, ambient });
    solve_heat(&g.model, None, limits).unwrap()
}

#[test]
fn lumped_body_cools_by_radiation_as_the_closed_form() {
    // Time to cool from T₀ to T: t = (F(T₀) − F(T)) / b, b = ε σ A / (ρ c V) = 6 ε σ / (ρ c s),
    // F the antiderivative; the temperature at time t is the root of that.
    let (s, eps, ambient, t0) = (0.01, 0.7, 300.0, 1000.0);
    let b = 6.0 * eps * STEFAN_BOLTZMANN / (7800.0 * 460.0 * s);
    let exact = |t: f64| root(|temp| t - (antiderivative(t0, ambient) - antiderivative(temp, ambient)) / b, ambient + 1e-9, t0);
    // ...and the 0 K limit: T = T₀ (1 + 3 b T₀³ t)^(−1/3). Long enough to lose most of the heat.
    let end = 600.0;
    let t_end = exact(end);
    assert!(t_end < 0.75 * t0, "{t_end}");
    let mut errors = Vec::new();
    for increments in [80, 160, 320] {
        let r = cooling_cube(s, eps, ambient, t0, end, increments, 0.5);
        let t = r.temperature();
        let spread = t.iter().fold(0.0_f64, |m, v| m.max((v - t[0]).abs()));
        assert!(spread < 1e-9 * t0, "not uniform by {spread}");
        errors.push((t[0] - t_end).abs());
        // The whole history follows the closed form.
        for (time, frame) in r.times.iter().zip(&r.frames) {
            assert!((frame[0] - exact(*time)).abs() < 1e-3 * t0, "{increments} increments, at {time}: {} vs {}", frame[0], exact(*time));
        }
    }
    println!("lumped radiation, Crank–Nicolson, 80/160/320 increments: T(600 s) = {t_end:.6} K, errors {errors:?}");
    // Crank–Nicolson with the radiation weighted by θ too: the error falls fourfold per halving.
    assert!(errors[2] < 1e-4 * t0, "{errors:?}");
    for w in errors.windows(2) {
        let order = (w[0] / w[1]).log2();
        assert!((1.8..2.2).contains(&order), "order {order} from {errors:?}");
    }
    // Backward Euler converges too, first order.
    let be: Vec<f64> = [160, 320].map(|n| (cooling_cube(s, eps, ambient, t0, end, n, 1.0).temperature()[0] - t_end).abs()).to_vec();
    println!("backward Euler, 160/320 increments: errors {be:?}");
    assert!((0.8..1.2).contains(&(be[0] / be[1]).log2()), "{be:?}");
    // To 0 K: the closed form T₀ (1 + 3 b T₀³ t)^(−1/3).
    let r = cooling_cube(s, eps, 0.0, t0, end, 320, 0.5);
    let want = t0 * (1.0 + 3.0 * b * t0.powi(3) * end).powf(-1.0 / 3.0);
    println!("to 0 K: {:.6} K vs {want:.6} K", r.temperature()[0]);
    assert!((r.temperature()[0] - want).abs() < 1e-4 * want, "{} vs {want}", r.temperature()[0]);
}

/// The plate of [`plate`] with only the heater (none at `flux` 0) and the radiation, as a model.
fn radiating_plate(flux: f64, ambient: f64) -> Grid {
    let mut g = grid([2, 1, 1], ElementKind::Hex8, |p| [0.02 * p[0], p[1], p[2]]);
    g.define(iso(2e11, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    conductor(&mut g, None, None);
    g.face_set("heated", |p| near(p[0], 0.0));
    g.face_set("cooled", |p| near(p[0], 1.0));
    if flux != 0.0 {
        g.load("heater", Load::SurfaceFlux { set: "heated".into(), flux });
    }
    g.load("glow", Load::Radiation { set: "cooled".into(), emissivity: 0.5, ambient });
    g
}

#[test]
fn radiation_newton_iterations_are_bounded() {
    // The iteration cap stops the solve with an error, never a hang. (With convection beside the
    // radiation the lumped start is not the answer, which it is for radiation alone here.)
    let mut g = radiating_plate(2000.0, 300.0);
    g.load("air", Load::Convection { set: "cooled".into(), film: 10.0, sink: 295.0 });
    let limits = Limits { max_newton_iterations: 1, ..Limits::default() };
    let e = solve_heat(&g.model, None, &limits).unwrap_err();
    assert!(matches!(e, Error::NotConverged(_)), "{e}");
    // Enough iterations: a handful from the lumped start.
    let limits = Limits { max_newton_iterations: 12, ..Limits::default() };
    assert!(solve_heat(&g.model, None, &limits).is_ok());
}

#[test]
fn a_radiation_solve_stops_after_its_most_newton_iterations() {
    // One increment of 10⁸ s from 1 K under surroundings at 10⁴ K: the tangent at 1 K is nearly
    // nothing, so the first iteration overshoots by eight orders of magnitude and Newton's
    // method on T⁴ comes back by a quarter per iteration, more than MAX_RADIATION_ITERATIONS
    // of them. The solve stops there, naming the increment.
    let mut g = grid([1, 1, 1], ElementKind::Hex8, |p| p.map(|c| 0.01 * c));
    let tr = Transient { end_time: 1e8, increment: 1e8, theta: 1.0, initial_temperature: 1.0 };
    g.define(iso(2e11, 0.3), Some(7800.0), full(), Procedure::HeatTransfer { transient: Some(tr) });
    conductor(&mut g, Some(7800.0), Some(460.0));
    g.face_set("skin", |p| near(p[0], 1.0));
    g.load("glow", Load::Radiation { set: "skin".into(), emissivity: 0.7, ambient: 1e4 });
    let e = solve_heat(&g.model, None, &Limits::default()).unwrap_err();
    let text = e.to_string();
    assert!(matches!(e, Error::NotConverged(_)), "{e}");
    assert!(text.contains("increment 1") && text.contains(&format!("{MAX_RADIATION_ITERATIONS} Newton iterations")), "{text}");
}

#[test]
fn every_increment_has_its_own_newton_iterations() {
    // A cap of 5 Newton iterations is per increment, not over the step: 200 increments of the
    // cooling cube take some 400 iterations in all and finish.
    let (s, eps, ambient, t0, end) = (0.01, 0.7, 300.0, 1000.0, 600.0);
    let r = cooling_cube_with(s, eps, ambient, t0, end, 200, 0.5, &Limits { max_newton_iterations: 5, ..Limits::default() });
    let b = 6.0 * eps * STEFAN_BOLTZMANN / (7800.0 * 460.0 * s);
    let exact = root(|temp| end - (antiderivative(t0, ambient) - antiderivative(temp, ambient)) / b, ambient + 1e-9, t0);
    assert!((r.temperature()[0] - exact).abs() < 1e-4 * t0, "{} vs {exact}", r.temperature()[0]);
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_long_radiating_transient_finishes() {
    // An hour at Δt = 0.05 s: 72 000 increments, past what one step-wide count of
    // Limits::MAX_NEWTON_ITERATIONS would allow at two or more iterations each.
    let (s, eps, ambient, t0, end, increments) = (0.01, 0.7, 300.0, 1000.0, 3600.0, 72_000);
    assert!(2 * increments > Limits::MAX_NEWTON_ITERATIONS);
    let r = cooling_cube(s, eps, ambient, t0, end, increments, 0.5);
    let b = 6.0 * eps * STEFAN_BOLTZMANN / (7800.0 * 460.0 * s);
    let exact = root(|temp| end - (antiderivative(t0, ambient) - antiderivative(temp, ambient)) / b, ambient + 1e-9, t0);
    println!("an hour in {increments} increments: T = {:.9} K vs {exact:.9} K", r.temperature()[0]);
    assert!((r.temperature()[0] - exact).abs() < 1e-6 * t0, "{} vs {exact}", r.temperature()[0]);
}

#[test]
fn a_steady_solve_with_nothing_above_0_k_is_refused() {
    // Surroundings at 0 K, no heat supplied and no temperature held: every temperature is 0 K,
    // where radiation has no tangent. Refused with the reason, not as a singular matrix.
    let g = radiating_plate(0.0, 0.0);
    let e = solve_heat(&g.model, None, &Limits::default()).unwrap_err();
    let text = e.to_string();
    assert!(matches!(e, Error::Invalid(_)) && text.contains("0 K") && text.contains("raise the ambient"), "{text}");
    // Heat supplied, or warmer surroundings, solve.
    assert!(solve_heat(&radiating_plate(10.0, 0.0).model, None, &Limits::default()).is_ok());
    assert!(solve_heat(&radiating_plate(0.0, 300.0).model, None, &Limits::default()).is_ok());
}

#[test]
fn a_badly_conditioned_plate_still_converges() {
    // The radiating plate with its first element a nanometre thick beside three of 6.7 mm: the
    // conduction matrix spans eighteen orders of magnitude, the temperatures carry round-off of
    // some 10⁻⁷ of their size and the heat balance next to the thin element cannot be computed
    // more closely than its terms' round-off. Each equation is held to its own round-off, so the
    // solve converges, on the exact linear field to what the mesh can resolve.
    let (l, q, eps, ambient) = (0.02, 2000.0, 0.8, 293.15);
    let x = move |p: f64| if p <= 0.25 { p / 0.25 * 1e-9 } else { 1e-9 + (p - 0.25) / 0.75 * (l - 1e-9) };
    let mut g = grid([4, 1, 1], ElementKind::Hex8, move |p| [x(p[0]), p[1], p[2]]);
    g.define(iso(2e11, 0.3), None, full(), Procedure::HeatTransfer { transient: None });
    conductor(&mut g, None, None);
    g.face_set("heated", |p| near(p[0], 0.0));
    g.face_set("cooled", |p| near(p[0], 1.0));
    g.load("heater", Load::SurfaceFlux { set: "heated".into(), flux: q });
    g.load("glow", Load::Radiation { set: "cooled".into(), emissivity: eps, ambient });
    let r = solve_heat(&g.model, None, &Limits::default()).unwrap();
    let tl = (q / (eps * STEFAN_BOLTZMANN) + ambient.powi(4)).powf(0.25);
    for (label, p) in &g.param {
        let want = tl + q * (l - x(p[0])) / K;
        let t = temperature_at(&r, *label);
        assert!((t - want).abs() < 1e-6 * want, "x {}: {t} vs {want}", x(p[0]));
    }
}
