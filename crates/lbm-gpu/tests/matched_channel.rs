//! The matched-wall channel against Lee & Moser on the GPU (#1108): the CPU study of
//! `crates/lbm/tests/turbulence.rs` run on the device, with the same geometry, protocol, stop rule and
//! pass bars. The flow is read back every `sampling_interval` steps instead of every step, so that
//! the device is not stalled by a read-back of every step.
//!
//! Sub-sampling: the block statistic is the mean of each block's samples. Taking every m-th step of a
//! correlated series inside a block leaves the block mean an unbiased estimate of the same quantity, so
//! the stop rule keeps its meaning; the CPU study samples every step.

use std::f64::consts::PI;
use std::time::Instant;

#[path = "../../lbm/tests/common/mod.rs"]
mod common;
#[path = "verification/device.rs"]
mod device;

use common::lee_moser::{LM_1000_CHANNEL, LM_5200_CHANNEL, LeeMoserChannel, at_re_tau};
use common::*;
use device::gpu;
use fieldcraft_lbm::backend::Backend;
use fieldcraft_lbm::bounce::MATCHING_ROWS;
use fieldcraft_lbm::statistics::Blocks;
use fieldcraft_lbm::{Flow, Limits, Setup};
use fieldcraft_lbm_gpu::Streaming;
use fieldcraft_model::{Collision, FaceCondition, Procedure, Turbulence, WallModel};

/// The environment variable that names the grid (cells per δ) of the GPU study.
const GRID: &str = "FIELDCRAFT_MATCHED_CHANNEL_GPU";
/// The environment variable that names the friction Reynolds number: 1000, 2000 or 5200 (1000 when unset).
const RE_TAU: &str = "FIELDCRAFT_MATCHED_RE_TAU";
/// The read-back count per block: a block (one δ/u_τ) is read back at least this many times.
const SAMPLES_PER_BLOCK: usize = 32;
/// The periodic face condition of the channel's x and z faces.
const PERIODIC: FaceCondition = FaceCondition::Periodic {};
/// The bulk U_b⁺ of the CPU bounded regression (`a_coarse_matched_wall_channel_at_re_tau_1000_stays_stable_over_a_short_run`,
/// `crates/lbm/tests/turbulence.rs`, 8 cells per δ, the bounded protocol), as that test printed it on the CPU.
const CPU_BOUNDED_BULK: f64 = 19.792;

/// The read-back interval in steps for a block of `block` steps: every `block / SAMPLES_PER_BLOCK`
/// steps (at least one), so each block has at least `SAMPLES_PER_BLOCK` samples.
fn sampling_interval(block: usize) -> usize {
    (block / SAMPLES_PER_BLOCK).max(1)
}

#[test]
fn the_bulk_is_read_back_at_least_32_times_per_block() {
    // A block of 3200 steps is read back every 100 steps: 32 samples per δ/u_τ.
    assert_eq!(sampling_interval(3200), 100);
    // A block shorter than 32 steps is read back every step.
    assert_eq!(sampling_interval(10), 1);
    for block in [32usize, 33, 3198, 8320] {
        let m = sampling_interval(block);
        assert!(block / m >= SAMPLES_PER_BLOCK, "block {block} read back every {m}");
    }
}

#[test]
fn the_setup_refuses_re_tau_5200_at_16_cells_per_delta_and_accepts_re_tau_1000() {
    // The relaxation time τ − 1/2 = 0.3 n / (Re_τ U_b⁺) must reach 0.0001 (τ ≥ 0.5001): 16 cells per δ is too coarse at Re_τ 5200.
    let refused = matched_setup(16, &LM_5200_CHANNEL, &MATCHED_STUDY);
    assert!(refused.as_ref().is_err_and(|e| e.contains("relaxation time")), "{:?}", refused.err());
    assert!(matched_setup(16, &LM_1000_CHANNEL, &MATCHED_STUDY).is_ok());
}

/// The run's times in δ/u_τ: start-up discarded, the least time averaged, and the cap (a run that
/// reaches it without an accepted mean is inconclusive, not a pass). The CPU study's.
struct MatchedProtocol {
    discard: f64,
    average_min: f64,
    cap: f64,
}

/// The study's protocol (the CPU study's, `crates/lbm/tests/turbulence.rs`).
const MATCHED_STUDY: MatchedProtocol = MatchedProtocol { discard: 20.0, average_min: 50.0, cap: 400.0 };
/// The bounded regression's protocol (the CPU study's, in the heavy job).
const MATCHED_BOUNDED: MatchedProtocol = MatchedProtocol { discard: 2.0, average_min: 2.0, cap: 4.0 };
/// The block length of the stop rule, in δ/u_τ (the CPU study's).
const MATCHED_BLOCK: f64 = 1.0;
/// The accuracy bar against Lee & Moser (the CPU study's, fixed before the runs): C_f and every row's U⁺.
const MATCHED_BAR: f64 = 0.05;
/// The lattice sound speed, 1/√3 (c_s² = 1/3).
const LATTICE_SOUND_SPEED: f64 = 0.577_350_269_189_625_8;

/// The time step in lattice units for a grid of `n` cells per δ: dt = 0.1 Δx / U_b⁺ (the CPU study's).
fn time_step(n: usize, lm: &LeeMoserChannel) -> f64 {
    0.1 * (1.0 / n as f64) / lm.bulk
}

/// The matched-wall channel's set-up on `n` cells per δ at `lm`'s Re_τ: the box, the walls, the driving,
/// Smagorinsky's model with the matched wall, and the start-up field (the CPU study's). The set-up refuses
/// a grid whose relaxation time is below the stability bound.
fn matched_setup(n: usize, lm: &LeeMoserChannel, protocol: &MatchedProtocol) -> Result<(Setup, Flow), String> {
    let dx = 1.0 / n as f64;
    let (lx, lz) = ((2.0 * PI * n as f64).round() * dx, (PI * n as f64).round() * dx);
    let mut d = domain([0.0; 3], [lx, 2.0 + dx, lz], dx);
    for f in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    let dt = time_step(n, lm);
    let mut p = collide(protocol.cap, dt, None, Collision::Regularized);
    if let Procedure::Fluid { turbulence, .. } = &mut p {
        *turbulence =
            Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: Some(WallModel::SpaldingBounceMatched) });
    }
    // The protocol's run is longer than the default update cap on the finer grids (as in the CPU study).
    let box_cells = (2.0 * PI * n as f64).round() as u64 * (2 * n as u64 + 1) * (PI * n as f64).round() as u64;
    let steps_cap = (protocol.cap / dt).round() as u64;
    let limits = Limits { max_cell_updates: box_cells * steps_cap, ..Limits::default() };
    let setup = Setup::new(&model(1.0, 1.0 / lm.re_tau as f64, d, p), None, &[], &limits).map_err(|e| e.to_string())?;
    let mut flow = Flow::new(&setup, &[]).map_err(|e| e.to_string())?;
    let grid = *flow.grid();
    let to_lattice = 1.0 / setup.velocity_scale();
    let (kx, kz) = (2.0 * PI / lx, 2.0 * PI / lz);
    flow.set_equilibrium(|cell| {
        let [x, y, z] = grid.centre(cell);
        let eta = y - 1.0 - 0.5 * dx;
        let wall = (1.0 - eta * eta).max(0.0);
        // libm's functions, so the start-up field's bits do not depend on the platform (as the CPU study).
        let u = 1.15 * lm.bulk * (1.0 - libm::pow(eta.abs().min(1.0), 8.0)) + 2.0 * libm::sin(2.0 * kz * z) * libm::cos(kx * x) * wall;
        let v = 1.5 * wall * libm::cos(2.0 * kz * z) * libm::sin(3.0 * kx * x);
        let w = 1.5 * wall * eta * libm::sin(kx * x) * libm::cos(kz * z);
        (1.0, [u * to_lattice, v * to_lattice, w * to_lattice])
    });
    Ok((setup, flow))
}

/// One grid of the GPU study: the time-averaged profile, its bulk speed, the stop and the cost.
struct MatchedRun {
    /// U⁺ per fluid row j = 1 … 2n − 1 (u_τ = 1 in SI), and the bulk U_b⁺ over the full 2δ.
    rows: Vec<f64>,
    bulk: f64,
    stop: String,
    steps: u64,
    averaged: f64,
    /// The largest local Mach number over the read-backs (not over every step).
    max_mach: f64,
    readbacks: u64,
    interval: usize,
    /// The span of one block in δ/u_τ: `interval` × the block's samples, not exactly one δ/u_τ.
    block_span: f64,
    wall: f64,
}

/// The matched-wall channel at `lm`'s Re_τ on `n` cells per δ, on the GPU (`Streaming::InPlace`), with the
/// set-up of [`matched_setup`] and the stop rule of the CPU study (`matched_channel_run` in
/// `crates/lbm/tests/turbulence.rs`). The flow is read back every [`sampling_interval`] steps; a block (one
/// δ/u_τ) gives `block / interval` samples, and the stop is checked at each block's end. `Ok(None)` when no
/// GPU adapter is present.
fn matched_channel_gpu(n: usize, lm: &LeeMoserChannel, protocol: &MatchedProtocol) -> Result<Option<MatchedRun>, String> {
    let started = Instant::now();
    let (setup, flow) = matched_setup(n, lm, protocol)?;
    let dt = time_step(n, lm);
    let grid = *flow.grid();
    let ny = grid.cells[1];
    // The fluid mask and each cell's row, from the CPU flow: the GPU flow does not expose them.
    let fluid: Vec<bool> = (0..grid.len()).map(|c| flow.is_fluid(c)).collect();
    let row: Vec<usize> = (0..grid.len()).map(|c| grid.coords(c)[1]).collect();
    let Some(mut device) = gpu(&flow, Streaming::InPlace) else { return Ok(None) };
    let discard = (protocol.discard / dt).round() as u64;
    let cap = (protocol.cap / dt).round() as u64;
    let block = (MATCHED_BLOCK / dt).round().max(1.0) as usize;
    let interval = sampling_interval(block);
    let block_samples = block / interval;
    let block_span = (block_samples * interval) as f64 * dt;
    // The burn-in counts the samples taken within the discarded start-up.
    let mut blocks = Blocks::new(block_samples, discard / interval as u64).map_err(|e| e.to_string())?;
    let mut row_sum = vec![0.0; ny];
    // The per-read-back row sums, reused across read-backs.
    let mut sum = vec![0.0; ny];
    let mut count = vec![0usize; ny];
    let mut samples = 0u64;
    let mut averaged = 0u64;
    let mut max_u = 0.0f64;
    let mut stop = String::from("time cap (inconclusive)");
    let mut steps = 0u64;
    let mut readbacks = 0u64;
    // The cap in whole read-back intervals.
    let cap_steps = cap / interval as u64 * interval as u64;
    while steps < cap_steps {
        device.advance(interval as u64).map_err(|e| e.to_string())?;
        steps += interval as u64;
        readbacks += 1;
        let fields = device.fields().map_err(|e| e.to_string())?;
        sum.fill(0.0);
        count.fill(0);
        for cell in (0..grid.len()).filter(|c| fluid[*c]) {
            let (Some(u), Some(j)) = (fields.velocity.get(cell), row.get(cell)) else { continue };
            // The mean is in SI (u_τ = 1); the lattice speed stays for the Mach number.
            sum[*j] += u[0] * setup.velocity_scale();
            count[*j] += 1;
            max_u = max_u.max((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt());
        }
        let mut bulk_sample = 0.0;
        for j in 1..ny - 1 {
            let mean = sum[j] / count[j].max(1) as f64;
            bulk_sample += mean;
            if steps > discard {
                row_sum[j] += mean;
            }
        }
        blocks.push(bulk_sample / (ny - 1) as f64).map_err(|e| e.to_string())?;
        samples += 1;
        if steps > discard {
            averaged += 1;
        }
        if samples.is_multiple_of(block_samples as u64)
            && averaged as f64 * interval as f64 * dt >= protocol.average_min
            && blocks.acceptance().accepted
        {
            stop = String::from("accepted by the block statistics");
            break;
        }
    }
    let avg = averaged.max(1) as f64;
    let rows: Vec<f64> = (1..ny - 1).map(|j| row_sum[j] / avg).collect();
    let bulk = rows.iter().sum::<f64>() / (ny - 1) as f64;
    Ok(Some(MatchedRun {
        rows,
        bulk,
        stop,
        steps,
        averaged: averaged as f64 * interval as f64 * dt,
        max_mach: max_u / LATTICE_SOUND_SPEED,
        readbacks,
        interval,
        block_span,
        wall: started.elapsed().as_secs_f64(),
    }))
}

/// The Re_τ the study runs: `FIELDCRAFT_MATCHED_RE_TAU` when set, which must be 1000, 2000 or 5200 (any
/// other value is refused, never run as 1000); 1000 when unset.
fn re_tau_from_env() -> &'static LeeMoserChannel {
    match std::env::var(RE_TAU) {
        Err(_) => &LM_1000_CHANNEL,
        Ok(v) => match v.parse::<u32>().ok().and_then(at_re_tau) {
            Some(lm) => lm,
            None => panic!("{RE_TAU}={v:?}: the Lee & Moser profiles are at Re_τ 1000, 2000 and 5200"),
        },
    }
}

/// The matched-wall channel against Lee & Moser on the GPU, on one grid (cells per δ from
/// `FIELDCRAFT_MATCHED_CHANNEL_GPU`) at the friction Reynolds number `FIELDCRAFT_MATCHED_RE_TAU`. Pass on a grid: C_f
/// within 5 % of Lee & Moser, (U_b⁺_LM / U_b⁺)² − 1, and U⁺ within 5 % on every row from the sample row (y = MATCHING_ROWS Δx
/// from the nearer wall) to the centreline. A time cap is inconclusive, never a pass. The bar is not adjusted to pass. Hours
/// long on the GPU (docs/verification.md, "Long studies"); f32 on the device.
#[test]
#[ignore = "hours: set FIELDCRAFT_MATCHED_CHANNEL_GPU=<cells per δ> on a GPU (docs/verification.md, \"Long studies\")"]
fn matched_wall_channel_on_the_gpu_against_lee_and_moser() {
    let Some(n) = std::env::var(GRID).ok().and_then(|v| v.parse::<usize>().ok()) else {
        println!("the GPU matched-wall channel runs only with {GRID}=<cells per δ> (hours long): skipped");
        return;
    };
    let lm = re_tau_from_env();
    let Some(run) = matched_channel_gpu(n, lm, &MATCHED_STUDY).expect("the matched-wall channel runs on the GPU") else {
        println!("no GPU adapter: skipped");
        return;
    };
    let cf_error = (lm.bulk / run.bulk).powi(2) - 1.0;
    println!(
        "matched wall on the GPU, Re_τ {}, {n} cells per δ: U_b⁺ {:.3} (Lee & Moser {:.3}), C_f error {:+.2} %",
        lm.re_tau,
        run.bulk,
        lm.bulk,
        100.0 * cf_error
    );
    println!(
        "stop: {} after {} steps, {:.1} δ/u_τ averaged; largest local Mach {:.3} over {} read-backs, every {} steps (a block spans {:.3} δ/u_τ); wall {:.1} s",
        run.stop, run.steps, run.averaged, run.max_mach, run.readbacks, run.interval, run.block_span, run.wall
    );
    // Kawai & Larsson's second condition (§III.B): h_wm about 0.1δ in outer units.
    let h_wm = MATCHING_ROWS as f64 / n as f64;
    println!("h_wm/δ = {MATCHING_ROWS}/{n} = {h_wm:.3}{}", if h_wm > 0.1 { " (above 0.1: outside Kawai & Larsson's condition)" } else { "" });
    let mut worst = 0.0f64;
    for (k, u) in run.rows.iter().enumerate() {
        let j = k + 1;
        // The distance to the nearer wall, in cells: the walls are the rows 0 and 2n.
        let cells = j.min(2 * n - j);
        let y = cells as f64 / n as f64;
        let reference = lm.at(y);
        let err = u / reference - 1.0;
        if cells >= MATCHING_ROWS {
            worst = worst.max(err.abs());
        }
        println!("  y/δ {y:.4}  U⁺ {u:.3}  Lee & Moser {reference:.3}  {:+.2} %", 100.0 * err);
    }
    assert!(run.stop.starts_with("accepted"), "inconclusive: {}", run.stop);
    assert!(cf_error.abs() < MATCHED_BAR, "C_f error {:+.2} %", 100.0 * cf_error);
    assert!(worst < MATCHED_BAR, "the largest row error from the sample row is {:.2} %", 100.0 * worst);
}

/// The matched-wall channel's bounded regression on the GPU, in the heavy job: 8 cells per δ at Re_τ 1000, 2 δ/u_τ
/// of start-up discarded and 2 averaged (the cap), the CPU bounded regression's code and settings, prints its U_b⁺ beside
/// the CPU's (`CPU_BOUNDED_BULK`). It asserts a stable run (largest local Mach below 0.3), finite positive row means and
/// U_b⁺ within 30 % of Lee & Moser's; it is not the accuracy bar.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn a_coarse_matched_wall_channel_on_the_gpu_stays_stable_over_a_short_run() {
    let lm = &LM_1000_CHANNEL;
    let Some(run) = matched_channel_gpu(8, lm, &MATCHED_BOUNDED).expect("the bounded matched-wall channel runs on the GPU") else {
        println!("no GPU adapter: skipped");
        return;
    };
    println!(
        "matched wall on the GPU, Re_τ 1000, 8 cells per δ, bounded: U_b⁺ {:.3} (the CPU bounded regression {CPU_BOUNDED_BULK:.3}; Lee & Moser {:.3}), {} steps, {:.1} δ/u_τ averaged, largest Mach {:.3}",
        run.bulk, lm.bulk, run.steps, run.averaged, run.max_mach
    );
    assert!(run.rows.iter().all(|u| u.is_finite() && *u > 0.0), "{:?}", run.rows);
    assert!(run.max_mach < 0.3, "largest local Mach {}", run.max_mach);
    assert!((run.bulk / lm.bulk - 1.0).abs() < 0.3, "U_b⁺ {} against {}", run.bulk, lm.bulk);
}
