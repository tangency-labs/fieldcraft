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

/// The hash of the 200-step channel run on the hybrid collision before the layer existed.
const DEVELOP_HASH: u64 = 0xe68c_caf8_0bc4_26d7;

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

// ---- the wall layer (`Flow::set_wall_layer_sigma_one`) ----

const CHANNEL_CELLS: usize = 8;
const CHANNEL_STEPS: usize = 200;

/// A 16-row channel as a one-cell column on the hybrid collision (σ 0.99), between two box-face walls or between two obstacle slabs (the
/// same fluid).
fn channel(slabs: bool) -> Flow {
    let dx = 1.0 / CHANNEL_CELLS as f64;
    let (lo, hi) = if slabs { (-dx, 2.0 + dx) } else { (0.0, 2.0) };
    let mut d = domain([0.0, lo, 0.0], [dx, hi, dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    let obstacles = if slabs {
        d.faces.y_min = FaceCondition::Periodic {};
        d.faces.y_max = FaceCondition::Periodic {};
        d.obstacles = vec!["Floor".into(), "Roof".into()];
        vec![
            obstacle("Floor", cuboid([-1.0, -dx, -1.0], [1.0 + dx, 0.0, 1.0 + dx]), &CUBOID_FACES),
            obstacle("Roof", cuboid([-1.0, 2.0, -1.0], [1.0 + dx, 2.0 + dx, 1.0 + dx]), &CUBOID_FACES),
        ]
    } else {
        d.faces.y_min = FaceCondition::Wall {};
        d.faces.y_max = FaceCondition::Wall {};
        Vec::new()
    };
    d.acceleration = [1.0, 0.0, 0.0];
    d.initial_velocity = [12.0, 0.0, 0.0];
    // lattice viscosity 0.0125 and a bulk speed near 0.1 (τ − ½ about 0.04)
    let dt = 0.1 * dx / 12.0;
    let mut p = collide(CHANNEL_STEPS as f64 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { sigma, .. } = &mut p {
        *sigma = Some(HrrSigma::Value(0.95));
    }
    prepare(&model(1.0, 1.0 / 20.0, d, p), &obstacles).1
}

/// ρ and u of every cell, as bits, folded into one hash.
fn hash(flow: &Flow) -> u64 {
    (0..flow.grid().len())
        .flat_map(|c| flow.moments(c).map_or([f64::NAN; 4], |(r, u)| [r, u[0], u[1], u[2]]))
        .fold(0xcbf2_9ce4_8422_2325u64, |h, x| (h ^ x.to_bits()).wrapping_mul(0x0000_0100_0000_01b3))
}

fn run_channel(slabs: bool, layer: Option<bool>) -> Flow {
    let mut flow = channel(slabs);
    if let Some(on) = layer {
        flow.set_wall_layer_sigma_one(on);
    }
    for _ in 0..CHANNEL_STEPS {
        flow.step().expect("a step");
    }
    flow
}

/// The layer is off by default, and off is the collision as it was before the layer existed (the hash the same run gave on develop).
#[test]
fn the_wall_layer_is_off_by_default_and_off_is_the_collision_before_it() {
    let flow = channel(true);
    assert!(!flow.wall_layer_sigma_one(), "the layer is off by default");
    assert!(flow.wall_sigma_layer_mask().iter().any(|&m| m), "the slabs' walls have a layer");
    let (default, off, on) = (hash(&run_channel(true, None)), hash(&run_channel(true, Some(false))), hash(&run_channel(true, Some(true))));
    println!("hashes: default {default:#018x} off {off:#018x} on {on:#018x}");
    assert_eq!(default, off, "the default is the layer off");
    assert_ne!(on, off, "the layer acts on a half-way wall");
    assert_eq!(off, DEVELOP_HASH, "off is the collision on develop before the layer");
}

/// The layer is what the walls are, not how they are described: the channel between box faces and between slabs of obstacle cells differs
/// by round-off with the layer off and with it on alike (one rule for both kinds of wall).
#[test]
fn a_box_face_wall_and_a_slab_wall_agree_with_the_layer_on_and_off() {
    let gap = |on: bool| {
        let (a, b) = (run_channel(false, Some(on)), run_channel(true, Some(on)));
        let (ga, gb) = (*a.grid(), *b.grid());
        let (mut worst, mut peak) = (0.0f64, 0.0f64);
        for y in 0..2 * CHANNEL_CELLS {
            let (p, q) = (a.moments(ga.index([0, y, 0])), b.moments(gb.index([0, y + 1, 0])));
            let (Some((ra, ua)), Some((rb, ub))) = (p, q) else { panic!("row {y} is not fluid") };
            worst = worst.max((ra - rb).abs());
            for k in 0..3 {
                worst = worst.max((ua[k] - ub[k]).abs());
            }
            peak = peak.max(ua[0]);
        }
        (worst, peak)
    };
    let (off, peak) = gap(false);
    let (on, _) = gap(true);
    println!("box face against slab, lattice units: layer off {off:.3e}, on {on:.3e} (peak speed {peak:.3e})");
    assert!(peak > 1e-3, "the channel did not accelerate: {peak:e}");
    assert!(off <= 1e-10 * peak.max(1.0), "the same wall differs by {off:e} with the layer off");
    assert!(on <= 1e-10 * peak.max(1.0), "the same wall differs by {on:e} with the layer on");
}

/// A symmetry plane and a slip face are not no-slip walls: a channel with a wall at the bottom and one of them at the top has the layer
/// in its bottom row only.
#[test]
fn a_symmetry_plane_and_a_slip_face_are_not_in_the_wall_layer() {
    for top in [FaceCondition::Symmetry {}, FaceCondition::Slip {}] {
        let dx = 1.0 / CHANNEL_CELLS as f64;
        let mut d = domain([0.0; 3], [dx, 1.0, dx], dx);
        for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
            *f = FaceCondition::Periodic {};
        }
        d.faces.y_min = FaceCondition::Wall {};
        d.faces.y_max = top;
        let mut p = collide(1e-3, 0.1 * dx / 12.0, None, Collision::Hrr);
        if let Procedure::Fluid { sigma, .. } = &mut p {
            *sigma = Some(HrrSigma::Value(0.99));
        }
        let flow = prepare(&model(1.0, 1.0 / 20.0, d, p), &[]).1;
        let grid = *flow.grid();
        let mask = flow.wall_sigma_layer_mask();
        let rows: Vec<bool> = (0..grid.cells[1]).map(|y| mask[grid.index([0, y, 0])]).collect();
        assert!(rows[0], "{top:?}: the wall's row is in the layer");
        assert!(rows[1..].iter().all(|&m| !m), "{top:?}: only the wall's row is in the layer: {rows:?}");
    }
}

/// σ = 1 in the layer alone, without the trace relaxation of the next commit, carries the cell-to-cell mode on a slot of 6 cells.
#[test]
fn the_wall_layer_alone_carries_a_mode_on_a_slot_of_6_cells() {
    let [wall, flux] = slot_alternation(6, 0.99, |f| f.set_wall_layer_sigma_one(true));
    println!("slot 6, layer only: wall-row alternation {wall:.3e}, flux alternation {flux:.3e}");
    assert!(wall > 1e-2 || flux > 1e-2, "no mode: {wall:e}, {flux:e}");
}

// ---- the trace relaxation (`Flow::set_trace_bulk_one`) with the layer ----

fn layer_and_trace(f: &mut Flow) {
    f.set_wall_layer_sigma_one(true);
    f.set_trace_bulk_one(true);
}

/// With σ = 1 in the layer and the trace relaxed at ω_b = 1 the slots of 6 and 8 cells carry no more alternation than the plain hybrid
/// collision (the layer alone has 50 times as much).
#[test]
fn the_layer_with_the_trace_relaxation_carries_no_mode_on_slots_of_6_and_8_cells() {
    for slot in [6, 8] {
        let plain = slot_alternation(slot, 0.99, |_| {});
        let both = slot_alternation(slot, 0.99, layer_and_trace);
        println!("slot {slot}: plain {plain:?}, layer and trace {both:?}");
        for k in 0..2 {
            assert!(both[k] < 2.0 * plain[k] + 1e-4, "slot {slot}: alternation {k} {:.3e} against {:.3e} for the plain collision", both[k], plain[k]);
        }
    }
}

/// The periodic duct's flow rate against the series with the layer and the trace relaxation: within 1 % at 8, 12 and 16 cells where the
/// plain collision is 2 to 6 % off, and closer than the plain collision's at every width.
#[test]
fn the_duct_flow_rate_with_the_layer_and_the_trace_relaxation_is_within_one_percent_of_the_series() {
    for w in [8usize, 12, 16] {
        let on = duct_ratio(w, layer_and_trace);
        let plain = duct_ratio(w, |_| {});
        println!("duct {w} cells: flow rate / exact {on:.5} with the layer and the trace relaxation, {plain:.5} plain");
        assert!((on - 1.0).abs() < (plain - 1.0).abs() / 2.0, "{w} cells: {on} is not twice as close to 1 as the plain collision's {plain}");
        assert!((on - 1.0).abs() < 0.01, "{w} cells: {on}");
    }
}
