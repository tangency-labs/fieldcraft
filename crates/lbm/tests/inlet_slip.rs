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
