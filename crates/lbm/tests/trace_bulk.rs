//! The hybrid collision's opt-in trace relaxation (the step's `trace_bulk_one`, docs/manual/fluid/collision.md, "Trace relaxation"): the
//! non-equilibrium density mode relaxed at one rate, ω_b = 1, in every hybrid cell. Two gates it is judged by, on the collision as it is:
//!
//! - a slot between two plates (common/slot_bed.rs), driven through at a held mean speed, for a cell-to-cell mode along x: the plain hybrid
//!   collision at σ 0.99 has none on slots of 6 and 8 cells, and the trace relaxation must not add one;
//! - a periodic square duct against the exact laminar series (Boussinesq), the flow rate at a given body force, for the 8, 12 and 16 cell
//!   ducts. The flow rate is read as the mean of the cell-centre speeds, which is the midpoint rule on the exact profile and so reads high by
//!   a known amount ([`midpoint_bias`]: 1.85, 0.82, 0.46 % at 8, 12, 16 cells); the corrected ratio divides it out.

mod common;

use common::slot_bed::{LX, PLATE, THICK, build};
use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_model::{Collision, FaceCondition, HrrSigma, Procedure};

/// What a slot run reads: the largest second difference along x, downstream of the plate, of the wall row's u_x and of the slot's flux,
/// over the slot's mean speed. `on` is the step's `trace_bulk_one`.
///
/// The controller below only holds the slot's mean speed, so the alternation measure does not depend on its constants: it moves the
/// uniform body force, which has no cell-to-cell component along x. A correction every 100 steps is about two acoustic crossings of the 32-cell
/// box (32 / c_s = 55 steps), so the speed it reads is settled; the factor 0.5 and the 60 corrections only set how
/// closely the mean speed is held; the clamp (0.8, 1.25) is reciprocal, so one correction changes the force by 25 % at most.
fn slot_alternation(slot: usize, sigma: f64, on: bool) -> [f64; 2] {
    let tau = 2.4e-3 * slot as f64 / 8.0;
    let u_lat = 0.047;
    let bed = build(slot, tau, u_lat, sigma, on).expect("the slot bed's setup");
    let mut flow = Flow::new(&bed.setup, &bed.obstacles).expect("the slot's flow");
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

/// The bar of a cell-to-cell mode: one per cent of the mean speed in the second difference along x. It is a regression bar, not a theory
/// bar: no discretisation-error estimate of the slot case gives it, and it is not fitted to the relaxation's output. The plain collision
/// measures 4e-3 at 6 cells (docs/manual/fluid/collision.md), so the bar sits above that baseline and says the baseline is itself clean. The
/// gate of the relaxation is the comparison with the plain collision ([`ADDED_SHARE`]).
///
/// The controller in `slot_alternation` (gain 0.5, the correction every 100 steps, the clamp (0.8, 1.25)) only holds the mean speed; the
/// measure is taken relative to that mean, so its constants do not enter the claim.
const MODE_BAR: f64 = 1e-2;

/// How much the relaxation may add to the plain collision's alternation: 10 % of it, with 1e-4 for a plain value near round-off (the
/// documented run has 4.1e-3 against 4.0e-3, 2.5 %).
const ADDED_SHARE: f64 = 0.1;
const ADDED_FLOOR: f64 = 1e-4;

/// The plain hybrid collision at σ 0.99 is free of the cell-to-cell mode on slots of 6 and 8 cells: the baseline the trace relaxation is gated
/// against.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_plain_hybrid_collision_has_no_mode_on_slots_of_6_and_8_cells() {
    for slot in [6, 8] {
        let [wall, flux] = slot_alternation(slot, 0.99, false);
        println!("slot {slot}, plain hybrid, sigma 0.99: wall-row alternation {wall:.3e}, flux alternation {flux:.3e}");
        assert!(wall < MODE_BAR && flux < MODE_BAR, "slot {slot}: the baseline carries a mode: {wall:e}, {flux:e}");
    }
}

/// With the trace relaxation on, the slots of 6 and 8 cells carry no more alternation than the plain hybrid collision, up to [`ADDED_SHARE`].
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_trace_relaxation_adds_no_mode_on_slots_of_6_and_8_cells() {
    for slot in [6, 8] {
        let plain = slot_alternation(slot, 0.99, false);
        let on = slot_alternation(slot, 0.99, true);
        println!("slot {slot}: plain {plain:?}, trace relaxation {on:?}");
        for k in 0..2 {
            assert!(
                on[k] <= (1.0 + ADDED_SHARE) * plain[k] + ADDED_FLOOR,
                "slot {slot}: alternation {k} {:.3e} against {:.3e} for the plain collision",
                on[k],
                plain[k]
            );
        }
    }
}

/// The mean speed the duct's body force is chosen to give the exact series.
const DUCT_SPEED: f64 = 0.04;

/// The exact developed laminar flow in a square duct of side `w` (Boussinesq): U = c G w² / ν with the coefficient below.
fn duct_coefficient() -> f64 {
    let pi = std::f64::consts::PI;
    let sum: f64 = (0..200).map(|m| 2 * m + 1).map(|n| ((n as f64) * pi / 2.0).tanh() / (n as f64).powi(5)).sum();
    (1.0 - 192.0 / pi.powi(5) * sum) / 12.0
}

/// The exact profile of the laminar square duct of side 1 at unit G/ν, Boussinesq's series: with odd n,
/// u = (4/π³) Σ n⁻³ (1 − cosh(nπ(z − ½)) / cosh(nπ/2)) sin(nπ y).
fn duct_profile(y: f64, z: f64) -> f64 {
    let pi = std::f64::consts::PI;
    (0..60)
        .map(|m| 2 * m + 1)
        .map(|n| {
            let k = n as f64 * pi;
            (1.0 - (k * (z - 0.5)).cosh() / (k / 2.0).cosh()) * (k * y).sin() / (n as f64).powi(3)
        })
        .sum::<f64>()
        * 4.0
        / pi.powi(3)
}

/// How far the mean of the cell-centre speeds of a `w` × `w` duct reads above the exact mean (the midpoint rule on the exact profile).
fn midpoint_bias(w: usize) -> f64 {
    let n = w as f64;
    let mean = (0..w).flat_map(|j| (0..w).map(move |k| ((j as f64 + 0.5) / n, (k as f64 + 0.5) / n))).map(|(y, z)| duct_profile(y, z)).sum::<f64>()
        / (n * n);
    mean / duct_coefficient() - 1.0
}

/// A periodic square duct of `w` cells at τ − ½ 6e-3 (hybrid, σ 0.99), driven by the lattice force that gives the exact duct a mean speed of
/// 0.04.
fn duct_flow(w: usize, on: bool) -> Flow {
    let tau = 6e-3;
    let nu = tau / 3.0;
    let dx = 1.0 / w as f64;
    let dt = nu * dx * dx;
    let mut d = domain([0.0; 3], [4.0 * dx, 1.0, 1.0], dx);
    d.faces.x_min = FaceCondition::Periodic {};
    d.faces.x_max = FaceCondition::Periodic {};
    let g_lat = DUCT_SPEED * nu / (duct_coefficient() * (w * w) as f64);
    d.acceleration = [g_lat * dx / (dt * dt), 0.0, 0.0];
    let mut p = collide(1e6 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { sigma, trace_bulk_one, .. } = &mut p {
        *sigma = Some(HrrSigma::Value(0.99));
        *trace_bulk_one = on;
    }
    prepare(&model(1.0, 1.0, d, p), &[]).1
}

/// The mean speed of the duct's cell centres over the exact series', run to a steady state: not yet corrected for the midpoint rule.
fn duct_ratio(w: usize, on: bool) -> f64 {
    let mut flow = duct_flow(w, on);
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
    mean(&flow) / DUCT_SPEED
}

/// The midpoint rule's bias of the cell-centre mean on the exact profile: 1.85, 0.82 and 0.46 % at 8, 12 and 16 cells, close to its leading
/// term 1/(24 c N²) with the series' own coefficient c = 0.03514.
#[test]
fn the_cell_centre_mean_of_the_exact_duct_reads_high_by_the_midpoint_rule() {
    let c = duct_coefficient();
    for (w, expected) in [(8usize, 1.848e-2), (12, 0.822e-2), (16, 0.463e-2)] {
        let b = midpoint_bias(w);
        let leading = 1.0 / (24.0 * c * (w * w) as f64);
        println!("duct {w} cells: midpoint bias {:.3} % (leading term {:.3} %)", 100.0 * b, 100.0 * leading);
        assert!((b - expected).abs() < 2e-5, "{w} cells: {b}");
    }
}

/// The plain hybrid collision's duct flow rate against the series, at 12 cells.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_plain_hybrid_duct_is_within_five_percent_of_the_series() {
    let r = duct_ratio(12, false) / (1.0 + midpoint_bias(12));
    println!("periodic duct, 12 cells, plain hybrid, sigma 0.99: flow rate / exact = {r:.5} (midpoint bias removed)");
    assert!((r - 1.0).abs() < 0.05, "the duct's flow rate is {r} of the exact series'");
}

/// The periodic duct's flow rate with the trace relaxation, over the series with the midpoint bias removed, at 8, 12 and 16 cells, beside the
/// plain collision's. A laminar duct has no non-equilibrium density mode to remove, so the relaxation leaves the flow rate where the plain
/// collision has it (within 0.1 %): the deficit of the hybrid collision there (3.7 % low at 12 cells) is the blend's, not the trace's.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_trace_relaxation_leaves_the_laminar_duct_flow_rate_where_the_plain_collision_has_it() {
    for w in [8usize, 12, 16] {
        let bias = 1.0 + midpoint_bias(w);
        let (on, plain) = (duct_ratio(w, true) / bias, duct_ratio(w, false) / bias);
        println!("duct {w} cells: flow rate / exact (midpoint bias removed) {on:.5} with the trace relaxation, {plain:.5} plain");
        assert!((on - plain).abs() < 1e-3, "{w} cells: {on} against {plain} for the plain collision");
    }
}

/// Lattice shear viscosity of the wave box.
const WAVE_NU: f64 = 0.01;

/// A periodic box of `n` cells along x (4 across) at lattice viscosity [`WAVE_NU`], hybrid collision at `sigma`, holding a plane density wave of
/// amplitude 1e-4 and wavelength `n` cells, at rest.
fn wave_box(n: usize, sigma: f64, on: bool) -> Flow {
    let dx = 1.0 / n as f64;
    let dt = WAVE_NU * dx * dx;
    let mut d = domain([0.0; 3], [n as f64 * dx, 4.0 * dx, 4.0 * dx], dx);
    d.faces.x_min = FaceCondition::Periodic {};
    d.faces.x_max = FaceCondition::Periodic {};
    d.faces.y_min = FaceCondition::Periodic {};
    d.faces.y_max = FaceCondition::Periodic {};
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    let mut p = collide(1e6 * dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { sigma: s, trace_bulk_one, .. } = &mut p {
        *s = Some(HrrSigma::Value(sigma));
        *trace_bulk_one = on;
    }
    let mut flow = prepare(&model(1.0, 1.0, d, p), &[]).1;
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let grid = *flow.grid();
    flow.set_equilibrium(|c| (1.0 + WAVE_AMPLITUDE * (k * (grid.centre(c)[0] * n as f64)).cos(), [0.0; 3]));
    flow
}

const WAVE_AMPLITUDE: f64 = 1e-4;

/// The decay rate of the wave's envelope, per step, by a least-squares line through ln of sqrt(a² + (b/c_s)²) (a the density wave's cosine
/// component, b the velocity's sine component: the standing wave's two quadratures) over the steps the amplitude is between 2 % and
/// 0.3 % of its start, sampled every 5 steps.
fn wave_decay(n: usize, sigma: f64, on: bool) -> f64 {
    let mut flow = wave_box(n, sigma, on);
    let grid = *flow.grid();
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let cs = (1.0f64 / 3.0).sqrt();
    let envelope = |flow: &Flow| {
        let (mut a, mut b) = (0.0, 0.0);
        for i in 0..n {
            if let Some((rho, u)) = flow.moments(grid.index([i, 0, 0])) {
                let x = k * (i as f64 + 0.5);
                a += (rho - 1.0) * x.cos();
                b += u[0] * x.sin();
            }
        }
        let (a, b) = (2.0 * a / n as f64, 2.0 * b / n as f64);
        (a * a + (b / cs) * (b / cs)).sqrt()
    };
    let start = envelope(&flow);
    assert!((start - WAVE_AMPLITUDE).abs() < 1e-3 * WAVE_AMPLITUDE, "the wave starts at {start}");
    let mut samples: Vec<(f64, f64)> = Vec::new();
    let mut step = 0u64;
    while step < 200_000 {
        Backend::advance(&mut flow, 5).expect("5 steps");
        step += 5;
        let e = envelope(&flow);
        assert!(e.is_finite(), "the wave left the finite numbers");
        if e < 0.02 * start && e > 0.003 * start {
            samples.push((step as f64, e.ln()));
        }
        if e <= 0.003 * start {
            break;
        }
    }
    assert!(samples.len() > 20, "{} samples of the decay", samples.len());
    let m = samples.len() as f64;
    let (st, sy) = samples.iter().fold((0.0, 0.0), |(a, b), (t, y)| (a + t, b + y));
    let (stt, sty) = samples.iter().fold((0.0, 0.0), |(a, b), (t, y)| (a + t * t, b + t * y));
    -(m * sty - st * sy) / (m * stt - st * st)
}

/// The error of the wave's measured decay rate against Γ = (k²/2)(4ν/3 + ζ), the bulk viscosity ζ added to the shear's share, and the bar it
/// is judged by. The hybrid collision's trace relaxes at the shear rate, so its bulk viscosity is 2ν/3 and Γ = ν k²; with the trace relaxation
/// on (ω_b = 1) ζ = (2/9)(1/ω_b − ½) = 1/9 lattice units. The bar was set before the run: a relative error of Γ at most (2π/N)², the lattice's own
/// second-order error in k.
fn wave_error(n: usize, sigma: f64, on: bool) -> f64 {
    let k = 2.0 * std::f64::consts::PI / n as f64;
    let zeta = if on { 1.0 / 9.0 } else { 2.0 * WAVE_NU / 3.0 };
    let theory = 0.5 * k * k * (4.0 * WAVE_NU / 3.0 + zeta);
    let got = wave_decay(n, sigma, on);
    let err = (got / theory - 1.0).abs();
    println!(
        "wave N {n}, sigma {sigma}, trace relaxation {}: decay {got:.4e} per step, theory {theory:.4e}, error {:.3} % (bar {:.2} %)",
        if on { "on" } else { "off" },
        100.0 * err,
        100.0 * k * k
    );
    assert!(err <= k * k, "N {n}, sigma {sigma}, on {on}: error {err} above the bar {}", k * k);
    err
}

/// The trace relaxation sets the bulk viscosity: on 32 cells at σ 0.99 the plane density wave decays at the theory's rate with the relaxation off
/// (ζ = 2ν/3) and on (ζ = 1/9, six times the rate).
#[test]
fn the_trace_relaxation_sets_the_bulk_viscosity_the_wave_decays_by() {
    let off = wave_error(32, 0.99, false);
    let on = wave_error(32, 0.99, true);
    println!("wave, 32 cells, sigma 0.99: errors {off:.3e} off, {on:.3e} on");
}

/// The same at σ 0.99 and 1, on 32 and 64 cells, with the observed order of the error in the wavelength (about 2 where the error is above round-off;
/// the plain collision at σ 1 is exact to 1e-6, where an order has no meaning).
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn the_wave_decay_matches_the_bulk_viscosity_on_32_and_64_cells_at_sigma_099_and_1() {
    for sigma in [0.99, 1.0] {
        for on in [false, true] {
            let (e32, e64) = (wave_error(32, sigma, on), wave_error(64, sigma, on));
            println!("sigma {sigma}, on {on}: observed order {:.2} (errors {e32:.3e}, {e64:.3e})", (e32 / e64).log2());
            if on {
                assert!((e32 / e64).log2() >= 1.5, "sigma {sigma}: order {}", (e32 / e64).log2());
            }
        }
    }
}

/// FNV-1a over the bits of every fluid cell's density and velocity.
fn fingerprint(flow: &Flow) -> u64 {
    let grid = *flow.grid();
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for c in 0..grid.len() {
        if let Some((rho, u)) = flow.moments(c) {
            for v in [rho, u[0], u[1], u[2]] {
                for byte in v.to_bits().to_le_bytes() {
                    h = (h ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
                }
            }
        }
    }
    h
}

/// Off, the collision is bit for bit what it was before the relaxation existed: the 8-cell duct's fields after 500 steps hash to the value
/// the collision gave then (taken from the tree without the relaxation), and with the relaxation on they do not.
#[test]
fn the_trace_relaxation_off_leaves_the_collision_bit_for_bit_as_it_was() {
    let run = |on: bool| {
        let mut flow = duct_flow(8, on);
        Backend::advance(&mut flow, 500).expect("500 steps");
        fingerprint(&flow)
    };
    let off = run(false);
    println!("duct, 8 cells, 500 steps: fingerprint {off:#018x} off, {:#018x} on", run(true));
    assert_eq!(off, PINNED_OFF, "the collision with the relaxation off changed: {off:#018x}");
    assert_ne!(off, run(true), "the relaxation did not act");
}

/// Taken on the tree without the relaxation, by the same test body.
const PINNED_OFF: u64 = 0xc412_164b_fcc4_e9c5;

/// A periodic box of 8 × 4 × 4 cells at τ − ½ 6e-3 (hybrid, σ 0.99) in uniform motion under a uniform body force, and the trace
/// Σ f_q (|c_q|² − 1) of its populations after one step. Uniform fields do not change in streaming and the finite-difference strain is zero,
/// so the populations after a step are the collision's output.
fn trace_after_one_step(on: bool) -> (f64, f64) {
    let tau = 6e-3;
    let nu = tau / 3.0;
    let dx = 1.0 / 8.0;
    let dt = nu * dx * dx;
    let mut d = domain([0.0; 3], [1.0, 4.0 * dx, 4.0 * dx], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    // The setup bounds the speed a body force reaches over the box, sqrt(2 a L) with L = 1 = 8 cells: lattice Mach sqrt(2 a_lat 8)/c_s
    // must stay below 0.3, so a_lat = 1e-3 (Mach 0.22).
    let a_lat = 1e-3;
    d.acceleration = [a_lat * dx / (dt * dt), 0.0, 0.0];
    let mut p = collide(dt, dt, None, Collision::Hrr);
    if let Procedure::Fluid { sigma, trace_bulk_one, .. } = &mut p {
        *sigma = Some(HrrSigma::Value(0.99));
        *trace_bulk_one = on;
    }
    let mut flow = prepare(&model(1.0, 1.0, d, p), &[]).1;
    let u0 = [0.05, 0.02, -0.01];
    flow.set_equilibrium(|_| (1.0, u0));
    let cell = 0;
    let (rho, u) = flow.moments(cell).expect("a fluid cell");
    let force = flow.acceleration().map(|a| rho * a);
    let expected = rho * u.iter().map(|v| v * v).sum::<f64>() + (0..3).map(|k| u[k] * force[k]).sum::<f64>();
    Backend::advance(&mut flow, 1).expect("one step");
    let trace: f64 = (0..fieldcraft_lbm::lattice::Q)
        .map(|q| {
            let c = fieldcraft_lbm::lattice::C[q];
            flow.population(q, cell) * f64::from(c[0] * c[0] + c[1] * c[1] + c[2] * c[2] - 1)
        })
        .sum();
    (trace, expected)
}

/// Theory: with the relaxation on, a uniformly accelerated periodic box leaves the collision with the trace ρ u² + u·F, the equilibrium's ρ u²
/// and the force's share at ω_b = 1, `(1 − ω_b/2) · 2 u·F = u·F`, with F = ρ a. The test checks this identity, which does not depend on the
/// code's own formula. It does not separate the relaxation from the plain collision on this bed: the plain hybrid collision's trace is already
/// within 1e-6 of the same sum here (printed), so the gate is the identity, not the difference.
#[test]
fn the_accelerated_periodic_box_leaves_the_collision_with_the_trace_rho_u_squared_plus_u_dot_f() {
    let (on, expected) = trace_after_one_step(true);
    let (off, _) = trace_after_one_step(false);
    println!("trace after one step: {on:.12e} with the relaxation, {off:.12e} without, ρu² + u·F = {expected:.12e}");
    assert!((on - expected).abs() < 1e-14, "{on} against {expected}");
}
