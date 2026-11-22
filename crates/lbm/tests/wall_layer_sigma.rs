//! The hybrid collision's opt-in wall treatment (docs/manual/fluid/collision.md, "Wall layer"): σ = 1 in the layer of the resolved walls
//! and one trace relaxation rate in every cell. This file holds the two gates it is judged by, on the collision as it is:
//!
//! - a slot between two plates (common/slot_bed.rs), driven through at a held mean speed, for the cell-to-cell mode along x that
//!   pure σ = 1 in the layer carries and that the plain hybrid collision at σ 0.99 already carries on a slot of 6 cells;
//! - a periodic square duct against the exact laminar series (Boussinesq), the flow rate at a given body force.

mod common;

use common::slot_bed::{LX, PLATE, THICK, build};
use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure};

/// What a slot run reads: the largest second difference along x, downstream of the plate, of the wall row's u_x and of the slot's flux,
/// over the slot's mean speed. `setup` is applied to the flow before it runs.
fn slot_alternation(slot: usize, sigma: f64, setup: impl Fn(&mut Flow)) -> [f64; 2] {
    let tau = 2.4e-3 * slot as f64 / 8.0;
    let u_lat = 0.047;
    let bed = build(slot, tau, u_lat, sigma).expect("the slot bed's setup");
    let mut flow = Flow::new(&bed.setup, &bed.obstacles).expect("the slot's flow");
    setup(&mut flow);
    let grid = *flow.grid();
    let u_at = |flow: &Flow, x: usize, y: usize| flow.moments(grid.index([x, y, 0])).map_or(f64::NAN, |(_, u)| u[0]);
    let rows = slot..2 * slot;
    let mut a = bed.g_lat;
    for _ in 0..60 {
        Backend::advance(&mut flow, 100).expect("100 steps");
        let mean = rows.clone().map(|y| u_at(&flow, PLATE, y)).sum::<f64>() / slot as f64;
        assert!(mean.is_finite(), "the slot flow left the finite numbers");
        // hold the slot's mean speed: a bounded correction of the force each 100 steps
        a *= (1.0 + 0.5 * (u_lat / mean.max(1e-3) - 1.0)).clamp(0.8, 1.25);
        Backend::set_acceleration(&mut flow, [a, 0.0, 0.0]).expect("the body force");
    }
    let mean = rows.clone().map(|y| u_at(&flow, PLATE, y)).sum::<f64>() / slot as f64;
    let wall_row = |x: usize| u_at(&flow, x, slot);
    let flux = |x: usize| rows.clone().map(|y| u_at(&flow, x, y)).sum::<f64>();
    let alt = |f: &dyn Fn(usize) -> f64, norm: f64| {
        (PLATE + THICK + 1..LX - 1).map(|x| (f(x + 1) - 2.0 * f(x) + f(x - 1)).abs()).fold(0.0, f64::max) / norm
    };
    [alt(&wall_row, mean), alt(&flux, mean * slot as f64)]
}

/// The plain hybrid collision at σ 0.99 is free of the cell-to-cell mode on slots of 6 and 8 cells: the baseline the wall treatment is gated
/// against (the alternation stays below `2 × baseline + 1e-4`).
#[test]
fn the_plain_hybrid_collision_has_no_mode_on_slots_of_6_and_8_cells() {
    for slot in [6, 8] {
        let [wall, flux] = slot_alternation(slot, 0.99, |_| {});
        println!("slot {slot}, plain hybrid, sigma 0.99: wall-row alternation {wall:.3e}, flux alternation {flux:.3e}");
        assert!(wall < 1e-2 && flux < 1e-2, "slot {slot}: the baseline carries a mode: {wall:e}, {flux:e}");
    }
}

/// The exact developed laminar flow in a square duct of side `w` (Boussinesq): U = c G w² / ν with the coefficient c below.
fn duct_coefficient() -> f64 {
    let pi = std::f64::consts::PI;
    let sum: f64 = (0..200).map(|m| 2 * m + 1).map(|n| ((n as f64) * pi / 2.0).tanh() / (n as f64).powi(5)).sum();
    (1.0 - 192.0 / pi.powi(5) * sum) / 12.0
}

/// The mean speed of a periodic square duct of `w` cells at τ − ½ 6e-3 (hybrid, σ 0.99), over the exact series', run to a steady state.
fn duct_ratio(w: usize, setup: impl Fn(&mut Flow)) -> f64 {
    let tau = 6e-3;
    let nu = tau / 3.0;
    let dx = 1.0 / w as f64;
    let dt = nu * dx * dx;
    let mut d = domain([0.0; 3], [4.0 * dx, 1.0, 1.0], dx);
    d.faces.x_min = FaceCondition::Periodic {};
    d.faces.x_max = FaceCondition::Periodic {};
    // the lattice force that drives U = 0.04 in the exact duct
    let u_target = 0.04;
    let g_lat = u_target * nu / (duct_coefficient() * (w * w) as f64);
    d.acceleration = [g_lat * dx / (dt * dt), 0.0, 0.0];
    let mut p = collide(1e6 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { sigma, .. } = &mut p {
        *sigma = Some(HrrSigma::Value(0.99));
    }
    let mut flow = prepare(&model(1.0, 1.0, d, p), &[]).1;
    setup(&mut flow);
    let grid = *flow.grid();
    let mean = |flow: &Flow| {
        let (mut s, mut n) = (0.0, 0);
        for y in 0..w {
            for z in 0..w {
                if let Some((_, u)) = flow.moments(grid.index([0, y, z])) {
                    s += u[0];
                    n += 1;
                }
            }
        }
        s / n as f64
    };
    let (mut last, mut steps) = (0.0, 0u64);
    while steps < 600_000 {
        Backend::advance(&mut flow, 5000).expect("5000 steps");
        steps += 5000;
        let now = mean(&flow);
        assert!(now.is_finite(), "the duct flow left the finite numbers");
        if steps > 20_000 && (now - last).abs() < 1e-7 * now.abs() {
            break;
        }
        last = now;
    }
    mean(&flow) / u_target
}

/// The plain hybrid collision's duct flow rate against the series, at 12 cells.
#[test]
fn the_plain_hybrid_duct_is_within_five_percent_of_the_series() {
    let r = duct_ratio(12, |_| {});
    println!("periodic duct, 12 cells, plain hybrid, sigma 0.99: flow rate / exact = {r:.5}");
    assert!((r - 1.0).abs() < 0.05, "the duct's flow rate is {r} of the exact series'");
}
