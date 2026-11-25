//! `inlet_ramp` on a bluff body: the far-field pressure rms of a 2-D cylinder started at rest, with no ramp, and with
//! the inlet ramped over one, two and two and a half periods T₀ of the box. These reproduce the manual's table
//! (`docs/manual/fluid/boundary-conditions.md`, "Raising the inlet"). Heavy: each run is 38 400 steps, several minutes
//! on one CPU core.
//!
//! Source: the expert probe of #785, CPU, f64, its runs of 2026-11-24 on develop 655d92fe with the probe's ramp: the
//! box 24 D × 12 D in one cell of span, 16 cells per D, the Hrr collision with Smagorinsky (C = 0.1) at Re 3900, halfway
//! walls, the inlet at lattice speed 0.05, the run to 120 D/U. The far-field probe is 3 D upstream of the centre and 3 D
//! above it, its density sampled every 90 steps; the rms is over the second half of the samples. The probe logged
//! 20.190, 7.001 and 0.049; this code prints 20.190, 7.020 (one period) and 0.049, the 0.3 % gap not isolated.
//!
//! The probe set its initial field to rest after giving the domain an initial velocity of 1; here the domain's initial
//! velocity is zero, so the start is the same rest field and the inlet is the only ramped quantity.
//!
//! The printed rms values are frozen: the run is bit-for-bit deterministic, so they must print as this code printed
//! them. Beside each ramped row the test prints the references of `common::ramp_tube`: the linear formula (the limit
//! Ma → 0), the ideal tube (the same ramp by characteristics, which keeps the Mach-number effects of the period and of
//! the round trip, sampled as the probe is), and the ideal tube times this run's own decay of the fundamental (printed,
//! not asserted, since that decay is read from the run it would check). The ideal tube's own error at this grid is the
//! Richardson bound on 384, 768 and 1536 intervals, printed with its observed order.
//!
//! Nothing is asserted against the ideal tube here. The rms carries the shedding's pressure variance at Re 3900, which
//! the one-dimensional tube does not have, and this test does not bound it; so the rows are regression values, and the
//! comparison to the ideal tube is printed. The sharp check of the solver against the ideal tube is
//! `inlet_ramp_tube.rs`, with no body in the tube.

mod common;

use common::ramp_tube::{RampTube, linear_rms, rest_period, richardson, rms};
use common::*;
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile, Turbulence, WallTreatment};

/// The box and the flow, in D = 1 units: 24 D long, 12 D high, 6 D of the cylinder's centre from the inlet.
const LENGTH_D: f64 = 24.0;
const HEIGHT_D: f64 = 12.0;
const UP_D: f64 = 6.0;
const CELLS_PER_D: f64 = 16.0;
const RE: f64 = 3900.0;
const LATTICE_SPEED: f64 = 0.05;
/// The end of the run in D/U (the model's time, D = U = 1).
const END_D_PER_U: f64 = 120.0;
/// The Smagorinsky constant of the LES study's cylinder case.
const SMAGORINSKY: f64 = 0.1;
/// Steps between far-field samples (the probe's).
const SAMPLE_EVERY: u64 = 90;
/// The ramp unit of the manual's rows, in steps: 4 · 384 · √3 rounded, the rest period of a tube as long as the box.
const PERIOD_STEPS: u64 = 2660;
/// The far-field probe's cell from the inlet: 3 D at 16 cells per D.
const PROBE_CELL: f64 = 48.0;

/// The box's cells along the flow.
fn cells() -> f64 {
    LENGTH_D * CELLS_PER_D
}

/// The far-field pressure coefficient of the probe after every step, C_p = (ρ − 1) / 3 / (½ U²) with c_s² = 1/3.
/// `ramp_steps` is the inlet's ramp in steps (none: no ramp).
fn far_field(ramp_steps: Option<u64>) -> Vec<f64> {
    let dx = 1.0 / CELLS_PER_D;
    let span = dx;
    let (xc, yc) = (UP_D, HEIGHT_D / 2.0);
    let dt = LATTICE_SPEED * dx;
    let mut d = domain([0.0; 3], [LENGTH_D, HEIGHT_D, span], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective, sigma: None };
    d.faces.y_min = FaceCondition::Slip {};
    d.faces.y_max = FaceCondition::Slip {};
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    d.obstacles = vec!["Body".into()];
    let body = obstacle("Body", cylinder([xc, yc], 0.5, -dx, span + dx, 1024), &["side", "bottom", "top"]);
    let procedure = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(END_D_PER_U),
        steps: None,
        time_step: Some(dt),
        frames: 1,
        magic: None,
        sigma: None,
        collision: Some(Collision::Hrr),
        turbulence: Some(Turbulence::Smagorinsky { constant: SMAGORINSKY, prandtl: 0.9, van_driest: None, wall_model: None }),
        steady: None,
        walls: WallTreatment::Halfway,
        inlet_ramp: ramp_steps.map(|s| s as f64 * dt),
    };
    let m = model(1.0, 1.0 / RE, d, procedure);
    let obstacles = [body];
    let (setup, mut flow) = prepare(&m, &obstacles);
    assert_eq!(flow.inlet_ramp_steps(), ramp_steps.unwrap_or(0), "the ramp is not the requested number of steps");
    let grid = *flow.grid();
    let probe = grid.cell_at([(xc - 3.0).max(0.5), yc + 3.0, span / 2.0]).expect("the far-field probe is in the box");
    assert_eq!(grid.cell_at([(PROBE_CELL + 0.5) * dx, yc + 3.0, span / 2.0]), Some(probe), "the probe is not cell {PROBE_CELL} from the inlet");
    let mut cp = Vec::with_capacity(setup.steps as usize);
    for _ in 0..setup.steps {
        flow.step().unwrap();
        let (rho, _) = flow.moments(probe).expect("the probe is a fluid cell");
        cp.push((rho - 1.0) / 3.0 / (0.5 * LATTICE_SPEED * LATTICE_SPEED));
    }
    cp
}

/// The probe's rms: a sample every `SAMPLE_EVERY` steps, over the second half of the samples.
fn sampled_rms(cp: &[f64]) -> f64 {
    let samples: Vec<f64> = (1..=cp.len() / SAMPLE_EVERY as usize).map(|k| cp[k * SAMPLE_EVERY as usize - 1]).collect();
    rms(&samples[samples.len() / 2..])
}

/// The pressure coefficient per unit of density: c_s² / (½ U²) with c_s² = 1/3.
fn cp_per_density() -> f64 {
    1.0 / 3.0 / (0.5 * LATTICE_SPEED * LATTICE_SPEED)
}

/// The ideal tube for this box, on `intervals` cells. Its acoustic length is between two bounds: the faces set their
/// cells' moments at the cell centres, 383 cells apart (`placement` 0), and a face acting at the box's face, half a cell
/// outside its centre, would add half a cell at each end, 384 (`placement` 1). The reference is the middle, 383.5; the
/// probe moves with the inlet's half cell.
fn ideal_tube(placement: f64, ramp_steps: u64, intervals: usize) -> RampTube {
    RampTube {
        length: cells() - 1.0 + placement,
        speed: LATTICE_SPEED,
        ramp_steps: ramp_steps as f64,
        probe: PROBE_CELL + 0.5 * placement,
        intervals,
    }
}

/// The ideal tube's rms, sampled as the probe is, on `intervals` cells.
fn ideal_rms(placement: f64, ramp_steps: u64, steps: u64, intervals: usize) -> f64 {
    let rho = ideal_tube(placement, ramp_steps, intervals).density(steps);
    sampled_rms(&rho) * cp_per_density()
}

/// The ideal tube's error at this box's grid, relative: the Richardson bound on its 384-interval rms from 768 and 1536
/// intervals, with the observed order.
fn reference_bound(ramp_steps: u64, steps: u64) -> (f64, f64) {
    let c = |intervals| ideal_rms(0.5, ramp_steps, steps, intervals);
    let (order, bound) = richardson(c(cells() as usize), c(2 * cells() as usize), c(4 * cells() as usize));
    (order, bound / c(cells() as usize))
}

/// The linear formula's rms for this box and probe (Ma → 0), the manual's fundamental.
fn linear(ramp_steps: u64) -> f64 {
    let tube = ideal_tube(0.5, ramp_steps, cells() as usize);
    linear_rms(tube.length, LATTICE_SPEED, ramp_steps as f64, tube.probe, 1)
}

/// The period of the fundamental with the inlet at its speed, T₀ / (1 − Ma²): the wave goes down at c_s + U and
/// comes back at c_s − U.
fn flow_period() -> f64 {
    rest_period(cells() - 0.5) / (1.0 - 3.0 * LATTICE_SPEED * LATTICE_SPEED)
}

/// The decay rate of the fundamental, per step: the least-squares slope of the logarithm of its amplitude over
/// consecutive windows of two periods from the ramp's end (the projection on the sine and cosine of `flow_period`).
fn decay_rate(cp: &[f64], ramp_steps: u64) -> f64 {
    let omega = 2.0 * std::f64::consts::PI / flow_period();
    let window = (2.0 * flow_period()).round() as usize;
    let mut logs = Vec::new();
    let mut start = ramp_steps as usize;
    while start + window <= cp.len() {
        let w = &cp[start..start + window];
        let mean = w.iter().sum::<f64>() / window as f64;
        let (mut s, mut c) = (0.0, 0.0);
        for (i, v) in w.iter().enumerate() {
            s += (v - mean) * (omega * i as f64).sin();
            c += (v - mean) * (omega * i as f64).cos();
        }
        logs.push((2.0 * (s * s + c * c).sqrt() / window as f64).ln());
        start += window;
    }
    let n = logs.len() as f64;
    let (mean_k, mean_log) = ((n - 1.0) / 2.0, logs.iter().sum::<f64>() / n);
    let slope = logs.iter().enumerate().map(|(k, l)| (k as f64 - mean_k) * (l - mean_log)).sum::<f64>()
        / logs.iter().enumerate().map(|(k, _)| (k as f64 - mean_k).powi(2)).sum::<f64>();
    -slope / window as f64
}

/// The references of a ramped row, printed: the linear formula, the ideal tube, and the ideal tube times this run's
/// decay from the ramp's end to the middle of the rms window. Returns the ideal tube's rms.
fn print_references(label: &str, cp: &[f64], ramp_steps: u64) -> f64 {
    let steps = cp.len() as u64;
    let ideal = ideal_rms(0.5, ramp_steps, steps, cells() as usize);
    let rate = decay_rate(cp, ramp_steps);
    // The rms window is the second half of the samples: its middle is three quarters of the run.
    let elapsed = 0.75 * steps as f64 - ramp_steps as f64;
    let (order, bound) = reference_bound(ramp_steps, steps);
    println!("{label}: linear formula rms {:.3}", linear(ramp_steps));
    println!("{label}: ideal tube rms {ideal:.3}; reference error {:.2e} relative at observed order {order:.2} (Richardson)", bound);
    println!(
        "{label}: decay of the fundamental {rate:.2e} per step; ideal tube times the decay {:.3} (printed, not asserted)",
        ideal * (-rate * elapsed).exp()
    );
    ideal
}

/// The rms must print as this code printed it, to its three decimals (a regression value, frozen).
fn assert_printed(label: &str, rms: f64, printed: &str) {
    println!("{label}: far-field C_p rms {rms:.3}");
    assert_eq!(format!("{rms:.3}"), printed, "{label}");
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn far_field_rms_without_a_ramp_is_a_regression_value() {
    // No reference: without a ramp the ideal tube's wave steepens into shocks, which the characteristics do not capture.
    assert_printed("no ramp", sampled_rms(&far_field(None)), "20.190");
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn far_field_rms_with_two_and_a_half_period_ramp_is_a_regression_value() {
    let ramp = PERIOD_STEPS * 5 / 2;
    let cp = far_field(Some(ramp));
    assert_printed("ramp 2.5 T0", sampled_rms(&cp), "0.049");
    // The ideal tube's value at the two bounds of the acoustic length: it moves by about a third between them, which
    // is the spread the 2.5 T0 comparison is printed against.
    let steps = cp.len() as u64;
    let (short, long) = (ideal_rms(0.0, ramp, steps, cells() as usize), ideal_rms(1.0, ramp, steps, cells() as usize));
    println!("ramp 2.5 T0: ideal tube at the bounds {short:.3} and {long:.3}, a spread of {:.0} %", 100.0 * (long - short).abs() / short);
    print_references("ramp 2.5 T0", &cp, ramp);
}

/// The one- and two-period rms, printed against their references. The ratio is printed next to the linear formula's
/// (5) and the ideal tube's.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn one_and_two_period_ramps_print_the_references() {
    let cp = far_field(Some(PERIOD_STEPS));
    let one = sampled_rms(&cp);
    assert_printed("ramp 1 T0", one, "7.020");
    let ideal_one = print_references("ramp 1 T0", &cp, PERIOD_STEPS);
    let cp_two = far_field(Some(2 * PERIOD_STEPS));
    let two = sampled_rms(&cp_two);
    assert_printed("ramp 2 T0", two, "1.350");
    let ideal_two = print_references("ramp 2 T0", &cp_two, 2 * PERIOD_STEPS);
    println!(
        "ratio rms(1 T0) / rms(2 T0) = {:.3}; linear formula {:.3}; ideal tube {:.3}",
        one / two,
        linear(PERIOD_STEPS) / linear(2 * PERIOD_STEPS),
        ideal_one / ideal_two
    );
}
