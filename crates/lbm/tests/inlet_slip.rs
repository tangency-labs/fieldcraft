//! #820: characterize the oscillatory inlet mode before repairing its closure.
//! The historical TRT diagnostic below is deliberately not an accuracy acceptance test.
//! A repair must replace its growth expectation with the 30,000-step decay sweep.
mod common;

use common::*;
use fieldcraft_lbm::{Backend, Flow};
use fieldcraft_model::{Collision, FaceCondition, Outflow, Profile};

fn box_flow(collision: Collision, tau: f64, sides: &str, seed: f64, size: usize) -> Flow {
    let mut d = domain([0.0; 3], [size as f64; 3], 1.0);
    let side = if sides == "slip" { FaceCondition::Slip {} } else { FaceCondition::Periodic {} };
    d.faces.y_min = side;
    d.faces.y_max = side;
    d.faces.z_min = side;
    d.faces.z_max = side;
    if sides == "closed" {
        d.faces.x_min = FaceCondition::Periodic {};
        d.faces.x_max = FaceCondition::Periodic {};
    } else {
        d.faces.x_min = FaceCondition::Velocity { velocity: [0.05, 0.0, 0.0], profile: Profile::Uniform };
        d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective };
    }
    d.initial_velocity = [0.05, 0.0, 0.0];
    let (_, mut flow) = prepare(&model(1.0, (tau - 0.5) / 3.0, d, collide(30_000.0, 1.0, None, collision)), &[]);
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| (1.0, [0.05 + if grid.coords(cell)[1] % 2 == 0 { seed } else { -seed }, 0.0, 0.0]));
    flow
}

fn departures(flow: &Flow) -> [f64; 2] {
    let mut peak = [0.0f64; 2];
    for cell in 0..flow.grid().len() {
        let (rho, u) = flow.moments(cell).unwrap();
        peak[0] = peak[0].max((rho - 1.0).abs());
        peak[1] = peak[1].max((u[0] - 0.05).abs()).max(u[1].abs()).max(u[2].abs());
    }
    peak
}

#[test]
fn uniform_inlet_zero_seed_is_an_equilibrium() {
    // MRT also amplifies roundoff here; it is not a stable control for this inlet.
    for collision in [Collision::Trt, Collision::Regularized, Collision::Hrr] {
        let mut flow = box_flow(collision, 0.505, "slip", 0.0, 8);
        flow.advance(100).unwrap();
        let peak = departures(&flow);
        assert!(peak.into_iter().all(|v| v <= 1e-12), "{collision:?}: {peak:?}");
    }
}

/// Historical evidence only. Once the repair is present this must assert decay, rather than
/// preserve the old guard failure. The sign-reversed disturbance is the same linear mode.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn characterize_the_unrepaired_trt_inlet_mode() {
    for seed in [1e-6, -1e-6] {
        let mut flow = box_flow(Collision::Trt, 0.505, "slip", seed, 16);
        let mut density_at_1500 = 0.0;
        let mut growth = 0.0;
        let mut failure = None;
        for step in 1..=4000 {
            if let Err(error) = flow.step() {
                assert!(matches!(error, fieldcraft_lbm::Error::Unstable(_)), "{error}");
                failure = Some(step);
                break;
            }
            if step == 1500 {
                density_at_1500 = departures(&flow)[0];
            } else if step == 2000 {
                growth = (departures(&flow)[0] / density_at_1500).ln() / 500.0;
            }
        }
        println!("historical #820: seed={seed}, guard={failure:?}, density growth={growth:.6}/step");
        assert!(failure.is_some_and(|step| (2500..3500).contains(&step)), "{failure:?}");
        assert!((0.004..0.006).contains(&growth), "{growth}");
    }
}

/// The original checkerboard controls. A repair also needs non-checkerboard transverse seeds:
/// periodic sides alone do not rule out a growing inlet mode (#820 investigation).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn inlet_checkerboard_controls_decay() {
    for (collision, tau, sides, seed) in [
        (Collision::Trt, 0.505, "periodic", 1e-6),
        (Collision::Trt, 0.505, "closed", 1e-6),
        (Collision::Trt, 0.505, "slip", 0.0),
        (Collision::Trt, 0.512, "slip", 1e-6),
        (Collision::Regularized, 0.505, "slip", 1e-6),
        (Collision::Hrr, 0.505, "slip", 1e-6),
    ] {
        let mut flow = box_flow(collision, tau, sides, seed, 16);
        flow.advance(30_000).unwrap();
        let peak = departures(&flow);
        println!("#820 control: {collision:?} tau={tau} sides={sides} seed={seed}: final density/velocity {peak:?}");
        assert!(peak.into_iter().all(|v| v <= 1e-6), "{peak:?}");
    }
}
