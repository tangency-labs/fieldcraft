//! Shear-dependent viscosity (`fieldcraft_lbm::rheology`) against closed forms: plane Poiseuille
//! flow of power-law fluids, shear-thinning and shear-thickening, with a resolution check; the
//! Newtonian limits of the power law (n = 1) and of Carreau–Yasuda (η₀ = η∞) against the
//! Newtonian solver; Carreau–Yasuda's two plateaus; and Cross-WLF's temperature, against the
//! Cross fluid of the same zero-shear viscosity.
//!
//! Each test prints its errors (`cargo test -p fieldcraft-lbm --test rheology -- --nocapture`).

mod common;

use common::*;
use fieldcraft_lbm::{Flow, Setup};
use fieldcraft_model::{Collision, Directional, FaceCondition, Fluid, FluidThermal, Material, Model, Thermal};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// Plates at y = 0 and y = H = 1; one cell across x and z, periodic.
const H: f64 = 1.0;

/// The peak speed the cases are built for (model units).
const U_MAX: f64 = 0.1;

/// A power-law fluid of index `n` whose plane Poiseuille flow peaks at [`U_MAX`] under the
/// acceleration it returns, with the viscosity `wall_nu` at the walls: `(fluid, G)`, with G = ρ a
/// (ρ = 1).
///
/// For η = K γ̇^(n−1) between plates at y = ±h the closed form is
/// `u(y) = n/(n+1) (G/K)^(1/n) (h^((n+1)/n) − |y|^((n+1)/n))` (Bird, Armstrong & Hassager 1987,
/// §4.2), with the wall shear rate `γ̇_w = (G h / K)^(1/n)`.
fn power_law(n: f64, wall_nu: f64) -> (Fluid, f64) {
    let h = H / 2.0;
    // (G/K)^(1/n) from the peak, then K from the wall viscosity K γ̇_w^(n−1).
    let ratio_root = U_MAX * (n + 1.0) / n / h.powf((n + 1.0) / n);
    let wall_rate = ratio_root * h.powf(1.0 / n);
    let k = wall_nu / wall_rate.powf(n - 1.0);
    let g = k * ratio_root.powf(n);
    (Fluid::PowerLaw { consistency: k, index: n, reference_shear_rate: None }, g)
}

/// The closed-form power-law profile at `y` (from the lower plate).
fn exact(fluid: &Fluid, g: f64, y: f64) -> f64 {
    let Fluid::PowerLaw { consistency: k, index: n, .. } = *fluid else { return f64::NAN };
    let h = H / 2.0;
    let e = (n + 1.0) / n;
    n / (n + 1.0) * (g / k).powf(1.0 / n) * (h.powf(e) - (y - h).abs().powf(e))
}

/// A channel of `cells` across with `fluid` (ρ = 1) driven by acceleration `g`, the time step
/// `Δx²` (the lattice viscosity is the model's), run for `time`.
fn channel(fluid: Fluid, g: f64, cells: usize, time: f64, collision: Collision) -> (Setup, Flow) {
    channel_at(fluid, None, g, cells, time, collision)
}

/// [`channel`], with heat in the fluid at a uniform `kelvin` when given (adiabatic walls, no
/// buoyancy: the temperature stays where it starts).
fn channel_at(fluid: Fluid, kelvin: Option<f64>, g: f64, cells: usize, time: f64, collision: Collision) -> (Setup, Flow) {
    let dx = H / cells as f64;
    let m = channel_model(fluid, kelvin, g, cells, collide(time, dx * dx, None, collision)).unwrap();
    let (setup, mut flow) = prepare(&m, &[]);
    run_all(&setup, &mut flow);
    (setup, flow)
}

/// The model of a channel of `cells` across with `fluid` (ρ = 1), heat at `kelvin` when given,
/// driven by acceleration `g`, solved by `procedure`.
fn channel_model(fluid: Fluid, kelvin: Option<f64>, g: f64, cells: usize, procedure: fieldcraft_model::Procedure) -> fieldcraft_model::Result<Model> {
    let dx = H / cells as f64;
    let mut d = domain([0.0; 3], [dx, H, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.acceleration = [g, 0.0, 0.0];
    d.thermal = kelvin.map(|t| FluidThermal { initial_temperature: t, ..FluidThermal::default() });
    let heat = kelvin.map(|_| Thermal { conductivity: Some(Directional::Isotropic(0.1)), specific_heat: Some(1.0), ..Thermal::default() });
    let mut m = Model::default();
    m.add_material("Fluid", Material { density: Some(1.0), fluid: Some(fluid), thermal: heat, ..Material::default() })?;
    m.add_fluid_domain("Flow", d)?;
    m.add_step(fieldcraft_model::Step::new("Run", procedure), None)?;
    Ok(m)
}

/// The largest difference between two profiles.
fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)
}

/// The streamwise velocity of every cell across the channel (model units).
fn profile(setup: &Setup, flow: &Flow) -> Vec<f64> {
    (0..setup.grid.cells[1]).map(|j| flow.moments(j).unwrap().1[0] * setup.velocity_scale()).collect()
}

/// The largest error of a power-law channel's profile over its peak, and the clamped share of
/// the last step's cells.
fn power_law_error(n: f64, cells: usize, collision: Collision) -> (f64, f64) {
    // The wall viscosity 1/6 is τ = 1 there; the time is several viscous times of the wall's
    // viscosity (the larger viscosity elsewhere settles sooner).
    let (fluid, g) = power_law(n, 1.0 / 6.0);
    let (setup, flow) = channel(fluid, g, cells, 8.0 * H * H * 6.0, collision);
    let dx = H / cells as f64;
    let u = profile(&setup, &flow);
    let err = u.iter().enumerate().map(|(j, v)| (v - exact(&fluid, g, (j as f64 + 0.5) * dx)).abs()).fold(0.0, f64::max) / U_MAX;
    let clamps = flow.viscosity_clamps().unwrap();
    let last = (clamps.last[0] + clamps.last[1]) as f64 / flow.fluid_cells() as f64;
    (err, last)
}

/// Power-law Poiseuille flow, shear-thinning (n = 0.5) and shear-thickening (n = 1.5), at 32
/// cells: 1.450e-3 and 1.278e-3 of the closed form's peak (the CPU is bit for bit the same on every
/// platform), held within about 10 % of that so a regression shows; closer at 64 (heavy).
#[test]
fn power_law_poiseuille_matches_the_closed_form() {
    for (n, bound) in [(0.5, 1.6e-3), (1.5, 1.4e-3)] {
        let (e, clamped) = power_law_error(n, 32, Collision::Trt);
        println!("power-law Poiseuille n = {n}, 32 cells: max error {e:.3e} of the peak, {:.1}% of cells clamped", 100.0 * clamped);
        assert!(e < bound, "n = {n}: {e}");
    }
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn power_law_poiseuille_converges_with_resolution() {
    for collision in COLLISIONS {
        for n in [0.5, 1.5] {
            let errors: Vec<f64> = [16, 32, 64].iter().map(|&c| power_law_error(n, c, collision).0).collect();
            println!("power-law Poiseuille {collision:?} n = {n}, 16 → 64 cells: errors {errors:?}");
            // 64 cells: 4.0e-4 and 4.1e-4 (TRT), 8.5e-4 and 1.4e-4 (regularized).
            assert!(errors[2] < errors[1] && errors[1] < errors[0], "n = {n}: {errors:?}");
            assert!(errors[2] < 1e-3, "n = {n}: {errors:?}");
        }
    }
}

/// The regularized and multiple-relaxation-time collisions take each cell's own rate the same way:
/// n = 0.5 at 32 cells, 1.455e-3 and 2.528e-3 of the peak, held within about 10 %.
#[test]
fn power_law_poiseuille_with_the_other_collisions() {
    for (collision, bound) in [(Collision::Regularized, 1.6e-3), (Collision::Mrt, 2.8e-3)] {
        let (e, _) = power_law_error(0.5, 32, collision);
        println!("power-law Poiseuille n = 0.5, {collision:?}, 32 cells: max error {e:.3e}");
        assert!(e < bound, "{collision:?}: {e}");
    }
}

/// n = 1, and Carreau–Yasuda with η₀ = η∞, are Newtonian: the same flow as the Newtonian solver
/// to round-off with each collision, with nothing clamped.
#[test]
fn newtonian_limits_reproduce_the_newtonian_solver() {
    let nu = 1.0 / 6.0;
    let g = 8.0 * nu * U_MAX / (H * H);
    let time = 2.0 * H * H / nu;
    for collision in [Collision::Trt, Collision::Regularized, Collision::Mrt] {
        let (setup, flow) = channel(Fluid::Newtonian { viscosity: nu }, g, 16, time, collision);
        let reference = profile(&setup, &flow);
        assert!(flow.viscosity_clamps().is_none() && setup.rheology.is_none());
        for fluid in [
            Fluid::PowerLaw { consistency: nu, index: 1.0, reference_shear_rate: None },
            Fluid::CarreauYasuda { zero_shear_viscosity: nu, infinite_shear_viscosity: nu, time_constant: 10.0, a: 2.0, index: 0.4 },
        ] {
            let (s, f) = channel(fluid, g, 16, time, collision);
            assert!((s.tau_plus - setup.tau_plus).abs() < 1e-12, "{} vs {}", s.tau_plus, setup.tau_plus);
            let d = max_diff(&profile(&s, &f), &reference) / U_MAX;
            println!("{collision:?} {}: largest difference from the Newtonian solver {d:.2e} of the peak", fluid.name());
            assert!(d < 1e-12, "{fluid:?}: {d}");
            assert_eq!(f.viscosity_clamps().unwrap().total, [0, 0]);
        }
    }
}

/// Carreau–Yasuda between its plateaus: at shear rates far below 1/λ it is the Newtonian fluid of
/// η₀, far above (with n close to 0) that of η∞.
#[test]
fn carreau_yasuda_plateaus_are_newtonian() {
    let (low, high) = (1.0 / 6.0, 1.0 / 12.0);
    let time = 2.0 * H * H / high;
    let newtonian = |nu: f64| {
        let g = 8.0 * nu * U_MAX / (H * H);
        let (s, f) = channel(Fluid::Newtonian { viscosity: nu }, g, 16, time, Collision::Trt);
        (g, profile(&s, &f))
    };
    // The wall shear rate is 4 U / H = 0.4; λ γ̇ ≤ 4e-7 keeps the zero-shear plateau.
    let (g, reference) = newtonian(low);
    let slow = Fluid::CarreauYasuda { zero_shear_viscosity: low, infinite_shear_viscosity: high, time_constant: 1e-6, a: 2.0, index: 0.3 };
    let (s, f) = channel(slow, g, 16, time, Collision::Trt);
    let d = max_diff(&profile(&s, &f), &reference) / U_MAX;
    println!("Carreau–Yasuda at λγ̇ ≤ 4e-7: {d:.2e} of the peak from Newtonian η₀");
    assert!(d < 1e-6, "{d}");
    // λ γ̇ ≥ 1e6 except within a millionth of the middle: η − η∞ ≤ (η₀ − η∞) (λγ̇)^(n−1) ≈ 1e-6.
    let (g, reference) = newtonian(high);
    let fast = Fluid::CarreauYasuda { zero_shear_viscosity: low, infinite_shear_viscosity: high, time_constant: 1e7, a: 2.0, index: 0.05 };
    let (s, f) = channel(fast, g, 16, time, Collision::Trt);
    let d = max_diff(&profile(&s, &f), &reference) / U_MAX;
    println!("Carreau–Yasuda at λγ̇ ≫ 1: {d:.2e} of the peak from Newtonian η∞");
    assert!(d < 1e-3, "{d}");
}

/// Cross-WLF takes each cell's temperature: in a melt held at a uniform temperature it is the
/// Cross fluid of that temperature's zero-shear viscosity (λ = η₀ / τ*), to round-off; a hotter
/// melt is thinner and flows faster. Without heat in the fluid it is refused by name.
#[test]
fn cross_wlf_follows_the_temperature() {
    let melt = |d1: f64| Fluid::CrossWlf { d1, d2: 373.0, a1: 28.0, a2: 51.6, critical_stress: 0.1, index: 0.5 };
    let kelvin = 400.0;
    // D1 for η₀ = 1/4 at 400 K: λ γ̇ ≈ 1 at the walls (γ̇_w ≈ 4 U / H), τ ≤ 1.25 everywhere.
    let shift = melt(1.0).viscosity(0.0, kelvin);
    let fluid = melt(0.25 / shift);
    let zero: f64 = fluid.viscosity(0.0, kelvin);
    assert!((zero - 0.25).abs() < 1e-12, "{zero}");
    let cross = Fluid::Cross { zero_shear_viscosity: zero, time_constant: zero / 0.1, index: 0.5 };
    let g = 8.0 * 0.15 * U_MAX / (H * H);
    let time = 8.0 * H * H * 6.0;
    let (s, f) = channel_at(fluid, Some(kelvin), g, 16, time, Collision::Trt);
    let (rs, rf) = channel_at(cross, Some(kelvin), g, 16, time, Collision::Trt);
    let (u, reference) = (profile(&s, &f), profile(&rs, &rf));
    let peak = reference.iter().copied().fold(0.0, f64::max);
    let d = max_diff(&u, &reference) / peak;
    println!("Cross-WLF at a uniform {kelvin} K: {d:.2e} of the peak from the Cross fluid of its η₀ (peak {peak:.4})");
    assert!(peak > 0.5 * U_MAX && d < 1e-10, "{peak} {d}");
    assert_eq!(f.viscosity_clamps().unwrap().total, [0, 0]);
    // 2 K hotter: η₀ falls by WLF (to 0.63 of it), the melt flows faster.
    let (hs, hf) = channel_at(fluid, Some(kelvin + 2.0), g, 16, time, Collision::Trt);
    let hot = profile(&hs, &hf).iter().copied().fold(0.0, f64::max);
    println!("Cross-WLF at {} K: peak {hot:.4}", kelvin + 2.0);
    assert!(hot > 1.2 * peak, "{hot} vs {peak}");
    // No heat in the fluid domain: refused by name.
    let e = channel_model(fluid, None, g, 16, step(1.0, 1.0 / 256.0, None)).unwrap_err().to_string();
    assert!(e.contains("Cross-WLF") && e.contains("thermal"), "{e}");
}

/// A turbulence model and porous zones with a shear-dependent fluid are refused by name; so is a
/// fluid with no velocity scale to take its viscosity at, and a solve reports clamping.
#[test]
fn unsupported_combinations_are_refused_and_clamping_is_reported() {
    let (fluid, g) = power_law(0.5, 1.0 / 6.0);
    let base = |f: &dyn Fn(&mut Model)| {
        let mut m = channel_model(fluid, None, g, 8, step(0.5, (H / 8.0).powi(2), None)).unwrap();
        f(&mut m);
        m
    };
    let turbulent = base(&|m| {
        let mut s = m.step("Run").unwrap().clone();
        if let fieldcraft_model::Procedure::Fluid { turbulence, .. } = &mut s.procedure {
            *turbulence = Some(fieldcraft_model::Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None });
        }
        m.replace_step(s).unwrap();
    });
    let e = Setup::new(&turbulent, None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("power-law") && e.contains("turbulence"), "{e}");
    let vreman = base(&|m| {
        let mut s = m.step("Run").unwrap().clone();
        if let fieldcraft_model::Procedure::Fluid { turbulence, .. } = &mut s.procedure {
            *turbulence = Some(fieldcraft_model::Turbulence::Vreman { constant: 0.18, prandtl: 0.9, wall_model: None });
        }
        m.replace_step(s).unwrap();
    });
    let e = Setup::new(&vreman, None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("power-law") && e.contains("turbulence"), "{e}");
    let still = base(&|m| {
        let mut d = m.fluid_domain("Flow").unwrap().clone();
        d.acceleration = [0.0; 3];
        m.replace_fluid_domain("Flow", d).unwrap();
    });
    let e = Setup::new(&still, None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("velocity scale"), "{e}");
    // The middle of a shear-thinning channel is at rest: its viscosity is capped, and said so.
    let r = fieldcraft_lbm::solve_fluid(&base(&|_| {}), None, &[], &fieldcraft_lbm::Limits::default()).unwrap();
    assert!(r.warnings.iter().any(|w| w.contains("capped")), "{:?}", r.warnings);
}

/// τ⁺ at the characteristic shear rate below τ_min is not refused when clamping handles it: a
/// shear-thinning fluid whose wall viscosity is under τ_min's runs, the floored cells counted and
/// the characteristic τ reported as such and held in the cells' range. A fluid floored at every
/// shear rate (Carreau–Yasuda with η₀ under τ_min's viscosity) is refused.
#[test]
fn a_characteristic_tau_below_the_floor_is_refused_only_when_clamping_cannot_help() {
    let cells = 8;
    let dx = H / cells as f64;
    // The lattice viscosity is the model's (Δt = Δx²): τ = 0.5 + 3ν, 0.502 at the wall.
    let (thinning, g) = power_law(0.5, 6e-4);
    let m = channel_model(thinning, None, g, cells, step(20.0 * dx * dx, dx * dx, None)).unwrap();
    let setup = Setup::new(&m, None, &[], &fieldcraft_lbm::Limits::default()).unwrap();
    assert!(setup.warnings.iter().any(|w| w.contains("characteristic shear rate is below")), "{:?}", setup.warnings);
    assert_eq!(setup.tau_plus, fieldcraft_lbm::setup::MIN_TAU);
    let r = fieldcraft_lbm::solve_fluid(&m, None, &[], &fieldcraft_lbm::Limits::default()).unwrap();
    let t = r.stats.shear_dependent_tau.unwrap();
    assert_eq!(t.range, [fieldcraft_lbm::setup::MIN_TAU, fieldcraft_lbm::rheology::MAX_TAU]);
    assert!(t.characteristic < fieldcraft_lbm::setup::MIN_TAU && Some(t) == setup.shear_dependent_tau, "{t:?}");
    assert!(r.warnings.iter().any(|w| w.contains("floored")), "{:?}", r.warnings);
    let thin = Fluid::CarreauYasuda { zero_shear_viscosity: 6e-4, infinite_shear_viscosity: 0.0, time_constant: 1.0, a: 2.0, index: 0.5 };
    let m = channel_model(thin, None, g, cells, step(20.0 * dx * dx, dx * dx, None)).unwrap();
    let e = Setup::new(&m, None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("Carreau–Yasuda") && e.contains("every shear rate"), "{e}");
}

/// A Newtonian fluid with a viscosity table takes each cell's temperature: held at a uniform
/// temperature it is the Newtonian fluid of the table's viscosity there (the n = 1 power law,
/// same collision path), to round-off; a hotter one, thinner in the table, flows faster. Without
/// heat in the fluid the constant holds, with a warning.
#[test]
fn a_viscosity_table_follows_the_temperature() {
    let table = fieldcraft_model::TemperatureTable(vec![[380.0, 0.3], [400.0, 0.2], [420.0, 0.15]]);
    let tabled = |kelvin: Option<f64>, procedure: fieldcraft_model::Procedure| {
        let mut m = channel_model(Fluid::Newtonian { viscosity: 0.2 }, kelvin, 8.0 * 0.15 * U_MAX / (H * H), 16, procedure).unwrap();
        let mut fluid = m.material("Fluid").unwrap().clone();
        fluid.temperature_dependent = Some(fieldcraft_model::TemperatureDependent { viscosity: Some(table.clone()), ..Default::default() });
        m.replace_material("Fluid", fluid).unwrap();
        m
    };
    let g = 8.0 * 0.15 * U_MAX / (H * H);
    let time = 8.0 * H * H * 6.0;
    let dx = H / 16.0;
    let run = |m: &Model| {
        let (setup, mut flow) = prepare(m, &[]);
        run_all(&setup, &mut flow);
        (profile(&setup, &flow), setup)
    };
    for kelvin in [390.0, 410.0] {
        let mu = table.value(kelvin);
        let (u, setup) = run(&tabled(Some(kelvin), collide(time, dx * dx, None, Collision::Trt)));
        assert!(setup.rheology.as_ref().is_some_and(|r| r.table.is_some()));
        let reference = Fluid::PowerLaw { consistency: mu, index: 1.0, reference_shear_rate: None };
        let (rs, rf) = channel_at(reference, Some(kelvin), g, 16, time, Collision::Trt);
        let expected = profile(&rs, &rf);
        let peak = expected.iter().copied().fold(0.0, f64::max);
        let d = max_diff(&u, &expected) / peak;
        println!("viscosity table at {kelvin} K (μ = {mu}): {d:.2e} of the peak from the constant fluid (peak {peak:.4})");
        assert!(peak > 0.3 * U_MAX && d < 1e-10, "{peak} {d}");
    }
    let peak = |kelvin| run(&tabled(Some(kelvin), collide(time, dx * dx, None, Collision::Trt))).0.iter().copied().fold(0.0, f64::max);
    assert!(peak(410.0) > 1.2 * peak(390.0));
    // No heat in the fluid: the constant viscosity, said so.
    let m = tabled(None, step(1.0, 1.0 / 256.0, None));
    let setup = Setup::new(&m, None, &[], &fieldcraft_lbm::Limits::default()).unwrap();
    assert!(setup.rheology.is_none());
    assert!(setup.warnings.iter().any(|w| w.contains("viscosity is held at its constant")), "{:?}", setup.warnings);
    // A table thinner somewhere than the grid resolves is refused as a constant that thin would
    // be, although its constant (0.2) is resolved.
    let mut thin = tabled(Some(400.0), collide(time, dx * dx, None, Collision::Trt));
    let mut fluid = thin.material("Fluid").unwrap().clone();
    fluid.temperature_dependent = Some(fieldcraft_model::TemperatureDependent {
        viscosity: Some(fieldcraft_model::TemperatureTable(vec![[380.0, 0.3], [420.0, 1e-6]])),
        ..Default::default()
    });
    thin.replace_material("Fluid", fluid).unwrap();
    let e = Setup::new(&thin, None, &[], &fieldcraft_lbm::Limits::default()).unwrap_err().to_string();
    assert!(e.contains("thinnest viscosity"), "{e}");
}

/// With a viscosity table, the solver's own time step and the velocity scale of an
/// acceleration-driven flow come from the table's thinnest viscosity, not the constant: a table
/// fifty times thinner than its constant somewhere still runs at the step the solver chose,
/// every cell's τ⁺ above the floor, and the Mach check sees the faster flow.
#[test]
fn a_viscosity_table_sizes_the_automatic_step_from_its_thinnest_viscosity() {
    let (g, cells) = (8.0 * 0.015 * U_MAX / (H * H), 16);
    let auto = fieldcraft_model::Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(1.0),
        steps: None,
        time_step: None,
        frames: 2,
        magic: None,
        sigma: None,
        collision: Some(Collision::Trt),
        turbulence: None,
        steady: None,
    };
    let mut m = channel_model(Fluid::Newtonian { viscosity: 0.2 }, Some(400.0), g, cells, auto).unwrap();
    let mut fluid = m.material("Fluid").unwrap().clone();
    fluid.temperature_dependent = Some(fieldcraft_model::TemperatureDependent {
        viscosity: Some(fieldcraft_model::TemperatureTable(vec![[380.0, 0.3], [420.0, 0.004]])),
        ..Default::default()
    });
    m.replace_material("Fluid", fluid).unwrap();
    let setup = Setup::new(&m, None, &[], &fieldcraft_lbm::Limits::default()).unwrap();
    let r = setup.rheology.as_ref().unwrap();
    let thinnest = 0.5 + 3.0 * r.table.as_ref().unwrap().range()[0] * r.viscosity_scale;
    assert!(thinnest >= fieldcraft_lbm::setup::MIN_TAU, "{thinnest}");
    // The box is one cell thick along x, the smallest extent the velocity scale takes.
    let dx = H / cells as f64;
    let at = |nu: f64| (g * dx * dx / (8.0 * nu)).min((2.0 * g * H).sqrt());
    assert_eq!(setup.velocity, at(0.004));
    assert!(setup.velocity > at(0.2));
    let constant_step = (fieldcraft_lbm::setup::AUTO_TAU_MAX - 0.5) * dx * dx / (3.0 * 0.2);
    assert!(setup.time_step > 10.0 * constant_step, "{} vs {constant_step}", setup.time_step);
}
