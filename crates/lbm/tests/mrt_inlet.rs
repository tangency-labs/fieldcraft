//! A velocity inlet under the MRT collision stays stable. The Zou-He closure leaves the
//! non-hydrodynamic part of the inlet cell's non-equilibrium unconstrained; the regularised
//! collisions project it every step, MRT relaxes it raw, so without the projection of the inlet
//! cell (`zou_he`, `project`) a plate-with-slot channel at low viscosity went unstable in
//! about 80 steps, in the first columns behind the inlet. This channel (uniform inlet at 0.05 in
//! lattice units with a cross-flow of 0.02, pressure outlet, no-slip walls on y) runs 3000 steps at
//! tau - 1/2 = 1e-2, 6e-3 and 5e-3 and must stay finite with a bounded speed (on develop's inlet
//! the first goes unstable at step 855). The inlet's accuracy is checked against plane Poiseuille
//! flow, and its stability against the seeded box of `inlet_slip.rs`.

mod common;

use common::inlet::{box_flow, broadband, departures};
use common::*;
use fieldcraft_lbm::Backend;
use fieldcraft_model::{Collision, FaceCondition, Outflow, Profile};

fn flow(tau_minus_half: f64, collision: Collision, steps: u64) -> fieldcraft_lbm::Flow {
    let nu = tau_minus_half / 3.0;
    let u = 0.05;
    let (nx, ny, nz) = (40usize, 12usize, 4usize);
    let mut d = domain([0.0; 3], [nx as f64, ny as f64, nz as f64], 1.0);
    for f in [&mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.02, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    d.initial_velocity = [u, 0.0, 0.0];
    let p = collide(steps as f64, 1.0, None, collision);
    let (_, mut flow) = prepare(&model(1.0, nu, d, p), &[]);
    for _ in 0..steps {
        flow.step().unwrap();
    }
    flow
}

fn peak_speed(flow: &fieldcraft_lbm::Flow) -> f64 {
    let mut peak = 0.0f64;
    for c in 0..flow.grid().len() {
        if flow.is_fluid(c)
            && let Some((_, v)) = flow.moments(c)
        {
            peak = peak.max((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt());
        }
    }
    peak
}

#[test]
fn a_velocity_inlet_under_mrt_stays_stable_at_low_viscosity() {
    for tau in [1e-2, 6e-3, 5e-3] {
        let peak = peak_speed(&flow(tau, Collision::Mrt, 3000));
        assert!(peak.is_finite() && peak < 0.2, "tau - 1/2 {tau}: peak speed {peak}");
    }
}

/// A seeded 16^3 slip box with a uniform inlet (the setup of `inlet_slip.rs`, where the TRT inlet
/// holds 30 000 steps), MRT, broadband seed 1e-6 on top of U = 0.05, at the viscosities where
/// develop's inlet table records Zou-He failing (tau - 1/2 = 5e-3 to 1e-2). Survival means finite
/// fields and a density departure below 0.1, the criterion that table stops on.
#[test]
fn a_seeded_box_under_mrt_survives_at_low_viscosity() {
    for tau in [0.505, 0.507, 0.51] {
        let mut flow = box_flow(Collision::Mrt, tau, "slip", 1e-6, [16; 3], 0, 1.0);
        broadband(&mut flow, 0, 1.0, 1e-6);
        let mut peak = [0.0f64; 2];
        for _ in 0..60 {
            flow.advance(50).unwrap();
            peak = departures(&mut flow, 0, 1.0);
            assert!(peak.into_iter().all(|v| v.is_finite() && v < 0.1), "tau {tau}, step {}: {peak:?}", flow.steps_done());
        }
        println!("MRT seeded box, tau {tau}: departures {peak:?} after {} steps", flow.steps_done());
    }
}

/// Plane Poiseuille flow in a channel of `h` cells (no-slip y walls, periodic z), a parabolic
/// inlet of peak U = 0.05 and a pressure outlet, tau = 0.6 under MRT, run to steady state
/// (10 h^2 / nu steps). The parabola u(y) = 4 U y (h - y) / h^2 is the exact solution; returns the
/// error relative to U in the first column behind the inlet and in the middle of the channel.
fn poiseuille_error(h: usize) -> [f64; 2] {
    let (u, nu) = (0.05, 0.1 / 3.0);
    let nx = 3 * h;
    let mut d = domain([0.0; 3], [nx as f64, h as f64, 1.0], 1.0);
    for f in [&mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d.faces.x_min = FaceCondition::Velocity { velocity: [u, 0.0, 0.0], profile: Profile::Parabolic };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    let steps = (10.0 * (h * h) as f64 / nu) as u64;
    let (_, mut flow) = prepare(&model(1.0, nu, d, collide(steps as f64, 1.0, None, Collision::Mrt)), &[]);
    flow.advance(steps).unwrap();
    let grid = *flow.grid();
    let column = |x: usize| {
        (0..h)
            .map(|j| {
                let y = j as f64 + 0.5;
                let exact = 4.0 * u * y * (h as f64 - y) / (h * h) as f64;
                let got = flow.moments(grid.index([x, j, 0])).map_or(f64::NAN, |(_, v)| v[0]);
                (got - exact).abs() / u
            })
            .fold(0.0, f64::max)
    };
    [column(1), column(nx / 2)]
}

/// The inlet reproduces the exact parabola next to it, and the error falls at second order with
/// the resolution (halfway bounce-back and Zou-He are second order): the observed order between
/// h = 16 and 32 is at least 1.8. The mid-channel column is printed, not asserted: it carries the
/// pressure outlet's error as well, which is not the inlet's.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_inlet_reproduces_plane_poiseuille_flow_at_second_order() {
    let (coarse, fine) = (poiseuille_error(16), poiseuille_error(32));
    println!("MRT Poiseuille mid channel: error {:.3e} (h = 16), {:.3e} (h = 32)", coarse[1], fine[1]);
    let order = (coarse[0] / fine[0]).log2();
    println!("MRT Poiseuille next to the inlet: error {:.3e} (h = 16), {:.3e} (h = 32), order {order:.2}", coarse[0], fine[0]);
    assert!(order > 1.8, "order {order}");
}
