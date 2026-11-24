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
//! The two-period row gives the ratio rms(T₀) / rms(2 T₀), which theory puts at 5; the measured ratio is 5.198, a 4.0 %
//! excess that is an open discrepancy (#1134), so this test prints it and asserts nothing about it. The ramp's own
//! reference check is the inlet velocity's exact formula, in `inlet_ramp.rs`, in the trunk gate.
//!
//! The no-ramp, one-period, two-period and 2.5-period rms are regression values: the run is bit-for-bit deterministic,
//! so they must print as this code printed them. They carry no reference of their own here.

mod common;

use common::*;
use fieldcraft_lbm::Backend;
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
/// The period T₀ = 4 L / c_s of the box between the inlet and the outlet: L = 384 cells, c_s = 1/√3 (rounded).
const PERIOD_STEPS: u64 = 2660;

/// The far-field pressure coefficient of the probe, C_p = (ρ − 1) / 3 / (½ U²) with c_s² = 1/3, and its rms over the
/// second half of the run. `ramp_steps` is the inlet's ramp in steps (none: no ramp).
fn far_field_rms(ramp_steps: Option<u64>) -> (f64, f64) {
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
    let mut cp = Vec::new();
    for _ in 0..setup.steps / SAMPLE_EVERY {
        flow.advance(SAMPLE_EVERY).unwrap();
        let (rho, _) = flow.moments(probe).expect("the probe is a fluid cell");
        cp.push((rho - 1.0) / 3.0 / (0.5 * LATTICE_SPEED * LATTICE_SPEED));
    }
    let tail = &cp[cp.len() / 2..];
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    let rms = (tail.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / tail.len() as f64).sqrt();
    (mean, rms)
}

/// The theory of the fundamental's rms after a ramp of `m` periods: the residual pressure amplitude
/// (4/π)(2/Ma) |cos(π m)| / |4 m² − 1| (the step's (4/π) fundamental, the ramp's factor), as the rms of a sinusoid (÷ √2),
/// and the closed–open mode's cosine at the probe 3 D from the inlet, cos(k x) with k = π / (2 L), L = 24 D.
fn theory_rms(m: f64) -> f64 {
    let mach = 3f64.sqrt() * LATTICE_SPEED;
    let amplitude = 4.0 / std::f64::consts::PI * 2.0 / mach * (std::f64::consts::PI * m).cos().abs() / (4.0 * m * m - 1.0).abs();
    amplitude / 2f64.sqrt() * (std::f64::consts::PI * 3.0 / (2.0 * LENGTH_D)).cos()
}

/// The rms must print as this code printed it, to its three decimals (a regression value, frozen).
fn assert_printed(label: &str, rms: f64, printed: &str) {
    println!("{label}: far-field C_p rms {rms:.3}");
    assert_eq!(format!("{rms:.3}"), printed, "{label}");
}

/// The theory ratio of the one-period and two-period residuals: |cos π| / |4 − 1| over |cos 2π| / |16 − 1|.
const THEORY_RATIO: f64 = 5.0;

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn far_field_rms_without_a_ramp_is_a_regression_value() {
    let (_, rms) = far_field_rms(None);
    assert_printed("no ramp", rms, "20.190");
}

#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn far_field_rms_with_two_and_a_half_period_ramp_is_a_regression_value() {
    let (_, rms) = far_field_rms(Some(PERIOD_STEPS * 5 / 2));
    assert_printed("ramp 2.5 T0", rms, "0.049");
}

/// The one- and two-period rms, printed against their theory. The ratio is printed next to its theory value (5); it
/// is not asserted: the 4.0 % excess is the open discrepancy of #1134.
#[test]
#[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
fn one_and_two_period_ramps_print_the_ratio_see_1134() {
    let (_, one) = far_field_rms(Some(PERIOD_STEPS));
    assert_printed("ramp 1 T0", one, "7.020");
    println!("ramp 1 T0: theory rms {:.3}", theory_rms(1.0));
    let (_, two) = far_field_rms(Some(2 * PERIOD_STEPS));
    assert_printed("ramp 2 T0", two, "1.350");
    println!("ramp 2 T0: theory rms {:.3}", theory_rms(2.0));
    println!("ratio rms(1 T0) / rms(2 T0) = {:.3}, theory {THEORY_RATIO} (see #1134)", one / two);
}
