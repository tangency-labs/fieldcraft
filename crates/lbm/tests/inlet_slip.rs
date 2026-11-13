//! #820: boundary-local TRT inlet repair and unchanged collision/side controls.
mod common;
use common::inlet::*;
use fieldcraft_lbm::Backend;
use fieldcraft_model::Collision;

#[test]
fn uniform_inlet_zero_seed_is_an_equilibrium() {
    // MRT also amplifies roundoff here; it is not a stable control for this inlet.
    for collision in [Collision::Trt, Collision::Regularized, Collision::Hrr] {
        let mut flow = box_flow(collision, 0.505, "slip", 0.0, [8; 3], 0, 1.0);
        flow.advance(100).unwrap();
        let peak = departures(&mut flow, 0, 1.0);
        assert!(peak.into_iter().all(|v| v <= 1e-12), "{collision:?}: {peak:?}");
    }
}

/// A two-dimensional oblique seed failed at step2574 before the repair. This small case
/// reaches beyond that failure in the ordinary gate; the full cube sweeps remain heavy.
#[test]
fn inlet_oblique_seed_decays() {
    let mut flow = box_flow(Collision::Trt, 0.505, "slip", 0.0, [16, 16, 1], 0, 1.0);
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| (1.0, [0.05 + 1e-6 * (std::f64::consts::PI * 11.0 * (grid.coords(cell)[1] as f64 + 0.5) / 16.0).cos(), 0.0, 0.0]));
    flow.advance(4000).unwrap();
    let peak = departures(&mut flow, 0, 1.0);
    println!("#820 oblique final {peak:?}");
    assert!(peak.into_iter().all(|v| v < 1e-6), "{peak:?}");
}

/// Original plus/minus checkerboard, broad-band seeds at two viscosities, all normal axes
/// and both inlet signs. Check decreasing late windows as well as the final departure.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn trt_inlet_perturbations_decay_for_thirty_thousand_steps() {
    for (tau, seed, random, axis, sign) in [
        (0.505, 1e-6, false, 0, 1.0),
        (0.505, -1e-6, false, 0, 1.0),
        (0.505, 1e-6, true, 0, 1.0),
        (0.518, 1e-6, true, 0, 1.0),
        (0.505, 1e-6, true, 0, -1.0),
        (0.505, 1e-6, true, 1, 1.0),
        (0.505, 1e-6, true, 1, -1.0),
        (0.505, 1e-6, true, 2, 1.0),
        (0.505, 1e-6, true, 2, -1.0),
    ] {
        let mut flow = box_flow(Collision::Trt, tau, "slip", seed, [16; 3], axis, sign);
        if random {
            broadband(&mut flow, axis, sign, seed);
        }
        let mut windows = [0.0f64; 6];
        for step in 0..600 {
            flow.advance(50).unwrap();
            windows[step / 100] = windows[step / 100].max(departures(&mut flow, axis, sign)[0]);
        }
        let peak = departures(&mut flow, axis, sign);
        println!("#820 tau={tau} seed={seed} random={random} axis={axis} sign={sign}: final={peak:?}, density windows={windows:?}");
        assert!(peak.into_iter().all(|v| v <= 1e-6), "{peak:?}");
        assert!(windows[2..].windows(2).all(|w| w[1] <= w[0] + 1e-12), "{windows:?}");
    }
}

/// The original checkerboard controls. A repair also needs non-checkerboard transverse seeds:
/// periodic sides alone do not rule out a growing inlet mode (#820 investigation).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn inlet_controls_decay() {
    for (collision, tau, sides, seed, random) in [
        (Collision::Trt, 0.505, "periodic", 1e-6, false),
        (Collision::Trt, 0.505, "closed", 1e-6, false),
        (Collision::Trt, 0.505, "slip", 0.0, false),
        (Collision::Trt, 0.512, "slip", 1e-6, false),
        (Collision::Regularized, 0.505, "slip", 1e-6, false),
        (Collision::Hrr, 0.505, "slip", 1e-6, false),
        (Collision::Trt, 0.505, "periodic", 1e-6, true),
        (Collision::Trt, 0.505, "closed", 1e-6, true),
    ] {
        let mut flow = box_flow(collision, tau, sides, seed, [16; 3], 0, 1.0);
        if random {
            broadband(&mut flow, 0, 1.0, seed);
        }
        flow.advance(30_000).unwrap();
        let peak = departures(&mut flow, 0, 1.0);
        println!("#820 control: {collision:?} tau={tau} sides={sides} seed={seed} random={random}: final density/velocity {peak:?}");
        assert!(peak.into_iter().all(|v| v <= 1e-6), "{peak:?}");
    }
}

/// Below the laminar floor, reachable only with a turbulence model, the TRT velocity inlet is not
/// stable at #820's Mach number: the setup says so, and only where a TRT velocity inlet exists.
#[test]
fn a_trt_velocity_inlet_below_the_floor_is_warned_about() {
    use common::*;
    use fieldcraft_lbm::{Limits, Setup};
    use fieldcraft_model::{FaceCondition, Outflow, Procedure, Profile, Turbulence};
    let setup = |tau: f64, collision: Collision, inlet: bool| {
        let mut d = domain([0.0; 3], [8.0; 3], 1.0);
        if inlet {
            d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
            d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
        } else {
            d.faces.x_min = FaceCondition::Periodic {};
            d.faces.x_max = FaceCondition::Periodic {};
        }
        let mut p = collide(10.0, 1.0, None, collision);
        if let Procedure::Fluid { turbulence, .. } = &mut p {
            *turbulence = Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None });
        }
        Setup::new(&model(1.0, (tau - 0.5) / 3.0, d, p), None, &[], &Limits::default()).unwrap().warnings
    };
    let warned = |w: Vec<String>| w.iter().any(|w| w.contains("velocity inlet is not stable"));
    assert!(warned(setup(0.5005, Collision::Trt, true)));
    assert!(!warned(setup(0.505, Collision::Trt, true)));
    assert!(!warned(setup(0.5005, Collision::Trt, false)));
    assert!(!warned(setup(0.5005, Collision::Hrr, true)));
}

/// Mass balance of the corrected inlet in a walled duct at steady state: the total mass stays put
/// and the mass entering through the inlet plane crosses every interior plane. The density
/// correction is high-passed along the face, so it adds no net mass
/// (`the_density_correction_adds_no_mass_over_a_face` in `inlet.rs`): the inflow is Zou–He's.
#[test]
fn a_duct_with_the_corrected_inlet_balances_its_mass_at_steady_state() {
    use common::*;
    use fieldcraft_model::{FaceCondition, Outflow, Profile};
    // About 1.3 viscous times H²/ν to settle, then a window to measure over.
    let (n, u, tau, settle, window) = ([24, 8, 8], 0.05, 0.52, 12_000, 2_000);
    let mut d = domain([0.0; 3], n.map(|n| n as f64), 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
    d.initial_velocity = [u, 0.0, 0.0];
    let (_, mut flow) = prepare(&model(1.0, (tau - 0.5) / 3.0, d, collide(1e5, 1.0, None, Collision::Trt)), &[]);
    let grid = *flow.grid();
    let measure = |flow: &mut fieldcraft_lbm::Flow| {
        let f = flow.fields().unwrap();
        let flux = |x: usize| -> f64 {
            (0..n[1]).flat_map(|y| (0..n[2]).map(move |z| grid.index([x, y, z]))).map(|c| f.density[c] * f.velocity[c][0]).sum()
        };
        let inflow = flux(0);
        let imbalance = (1..n[0] - 1).map(|x| ((flux(x) - inflow) / inflow).abs()).fold(0.0, f64::max);
        (f.density.iter().sum::<f64>(), inflow, imbalance)
    };
    flow.advance(settle).unwrap();
    let (m0, ..) = measure(&mut flow);
    flow.advance(window).unwrap();
    let (m1, inflow, imbalance) = measure(&mut flow);
    // Measured 7e-20 per step and 5e-14 (round-off).
    let drift = (m1 - m0).abs() / m0 / window as f64;
    println!("#820 duct: mass drift {drift:.2e} per step, inflow {inflow:.6}, largest plane imbalance {imbalance:.2e}");
    assert!(drift < 1e-15, "{drift}");
    assert!(imbalance < 1e-12, "{imbalance}");
}
