//! `inlet_ramp`: a velocity face rises from the initial field over a time instead of stepping to its velocity at the first step.
//! Theory (plane acoustics, Joukowsky): a velocity step U at the inlet of a tube sends p' = rho c_s U, so the density
//! rises by U / c_s = sqrt(3) U in lattice units. With the outlet held at rho = 1 the tube's fundamental period is
//! T0 = 4 L / c_s; a ramp of T = m T0 leaves the fundamental with |cos(pi m)| / |4 m^2 - 1| of the step's amplitude.
mod common;

use common::*;
use fieldcraft_lbm::{Error, Flow, Limits, Setup, Target, inlet::ramp};
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile};

/// Cells in the tube, and its lattice speed (Mach number 0.017).
const CELLS: usize = 64;
const SPEED: f64 = 0.01;
/// The tube's fundamental period T0 = 4 L / c_s in steps, with c_s = 1 / sqrt(3) (the rounded value).
const T0: u64 = 443;

/// The procedure of `collide`, with its `inlet_ramp` set to `ramp` seconds.
fn procedure(steps: u64, collision: Collision, ramp: Option<f64>) -> Procedure {
    let mut p = collide(steps as f64, 1.0, None, collision);
    if let Procedure::Fluid { inlet_ramp, .. } = &mut p {
        *inlet_ramp = ramp;
    }
    p
}

/// The tube: a velocity face at `x_min`, a uniform pressure face at `x_max`, periodic in the other directions.
fn tube(cells: usize, speed: f64) -> fieldcraft_model::FluidDomain {
    let mut d = domain([0.0; 3], [cells as f64, 1.0, 1.0], 1.0);
    d.faces.x_min = FaceCondition::Velocity { velocity: [speed, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None };
    for f in [&mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *f = FaceCondition::Periodic {};
    }
    d
}

/// The largest |rho - 1| at the cell next to the inlet over `window` (steps), in a tube started at rest and run
/// for `steps` steps, with the inlet ramped over `ramp` seconds.
fn pressure_wave(steps: u64, window: std::ops::Range<u64>, collision: Collision, ramp: Option<f64>) -> f64 {
    let m = model(1.0, 0.1 / 3.0, tube(CELLS, SPEED), procedure(steps, collision, ramp)); // tau+ 0.6
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &[]).unwrap();
    let mut worst = 0.0f64;
    for k in 1..=steps {
        flow.step().unwrap();
        if window.contains(&k) {
            worst = worst.max((flow.moments(1).unwrap().0 - 1.0).abs());
        }
    }
    worst
}

/// The modal bound B(m) of the step's wave in units of J = sqrt(3) U: (4/pi) sum_n 1 / ((2n+1) |4 m^2 (2n+1)^2 - 1|),
/// over the first three modes of the closed-open tube (the step's modal amplitudes are 4 / (pi (2n+1))).
fn modal_bound(m: f64) -> f64 {
    let sum: f64 = (0..3)
        .map(|n| {
            let k = (2 * n + 1) as f64;
            1.0 / (k * (4.0 * m * m * k * k - 1.0).abs())
        })
        .sum();
    4.0 / std::f64::consts::PI * sum
}

#[test]
fn the_ramp_is_the_raised_cosine() {
    assert_eq!(ramp(0, 100), 0.0);
    assert!((ramp(50, 100) - 0.5).abs() < 1e-15);
    assert!((ramp(25, 100) - 0.5 * (1.0 - (std::f64::consts::PI / 4.0).cos())).abs() < 1e-15);
    assert_eq!(ramp(100, 100), 1.0);
    assert_eq!(ramp(7, 0), 1.0);
    let mut last = 0.0;
    for k in 0..=100 {
        let s = ramp(k, 100);
        assert!(s >= last, "the ramp falls at step {k}");
        last = s;
    }
}

#[test]
fn a_start_from_rest_sends_joukowsky_s_wave_and_a_ramp_removes_it() {
    let t0 = T0;
    assert_eq!((4.0 * CELLS as f64 * 3f64.sqrt()).round() as u64, t0);
    let joukowsky = 3f64.sqrt() * SPEED;
    for collision in [Collision::Trt, Collision::Regularized] {
        // No ramp: linear acoustics gives the step's wave at its first quarter period, (1 - U / c_s) J at least.
        let first = pressure_wave(t0, 1..t0 / 4, collision, None);
        let bar = (1.0 - SPEED * 3f64.sqrt()) * joukowsky; // U / c_s = sqrt(3) U with c_s = 1 / sqrt(3)
        assert!(first >= bar, "{collision:?}: the step's wave is {first} < {bar}");
        // A ramp of one and two periods: the modal bound.
        for (m, ramp_steps) in [(1.0, t0), (2.0, 2 * t0)] {
            let seconds = ramp_steps as f64;
            let after = pressure_wave(ramp_steps + 3 * t0, ramp_steps.max(1)..ramp_steps + 3 * t0, collision, Some(seconds));
            let limit = modal_bound(m) * joukowsky;
            assert!(after <= limit, "{collision:?}, m {m}: |rho - 1| {after} > bound {limit}");
        }
    }
}

#[test]
fn a_ramp_needs_rest_and_a_velocity_face() {
    fn refusal(m: &fieldcraft_model::Model) -> String {
        match Setup::new(m, None, &[], &Limits::default()) {
            Err(Error::Invalid(text)) => text,
            Err(e) => panic!("unexpected error: {e}"),
            Ok(_) => panic!("the ramp was accepted"),
        }
    }
    // No velocity face.
    let mut d = tube(CELLS, SPEED);
    d.faces.x_min = FaceCondition::Periodic {};
    d.faces.x_max = FaceCondition::Periodic {};
    let m = model(1.0, 0.1 / 3.0, d, procedure(T0, Collision::Trt, Some(100.0)));
    assert!(refusal(&m).contains("velocity face"));
    // The ramp is longer than the run.
    let m = model(1.0, 0.1 / 3.0, tube(CELLS, SPEED), procedure(100, Collision::Trt, Some(1000.0)));
    assert!(refusal(&m).contains("longer than the step"));
}

#[test]
fn a_ramp_starts_from_the_initial_field() {
    // The domain's initial field is 0.005 and the inlet's target 0.01: u(t) = u0 + r(t) (u - u0) at every step.
    let mut d = tube(CELLS, SPEED);
    d.initial_velocity = [0.005, 0.0, 0.0];
    let m = model(1.0, 0.1 / 3.0, d, procedure(100, Collision::Trt, Some(100.0)));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &[]).unwrap();
    let full: Vec<Target> = flow.inlet_targets();
    let from = setup.initial_velocity;
    assert_ne!(from, [0.0; 3], "the test needs a non-zero initial field");
    for k in 1..=30u64 {
        flow.step().unwrap();
        let share = ramp(k, 100);
        for (oc, target) in flow.open_cells().iter().zip(&full) {
            if let (Target::Velocity(u), Target::Velocity(v)) = (oc.target, target) {
                let expected: [f64; 3] = std::array::from_fn(|c| from[c] + share * (v[c] - from[c]));
                assert_eq!(u, expected, "step {k}");
            }
        }
    }
}

#[test]
fn the_imposed_inlet_velocity_is_the_ramp_formula() {
    // The velocity the inlet cell imposes, read from its moments after k steps, is u0 + r(k) (SPEED − u0), with the
    // raised cosine written out here: r(k) = ½ (1 − cos(π k / T)) for k < T, and 1 from T on.
    const RAMP: u64 = 100;
    const INITIAL: f64 = 0.005;
    let mut d = tube(CELLS, SPEED);
    d.initial_velocity = [INITIAL, 0.0, 0.0];
    let m = model(1.0, 0.1 / 3.0, d, procedure(200, Collision::Trt, Some(RAMP as f64)));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &[]).unwrap();
    let inlet = flow.grid().index([0, 0, 0]);
    // The moments sum the 19 populations of the inlet cell, each of magnitude at most about 1 at density 1, so the
    // rounding of the imposed velocity is at most 19 ε (f64::EPSILON); Zou–He's reconstruction is exact in exact
    // arithmetic. The tolerance is that bound, not a fitted one.
    let tolerance = 19.0 * f64::EPSILON;
    for k in 1..=200u64 {
        flow.step().unwrap();
        if [1, 25, 50, 100, 101, 200].contains(&k) {
            let share = if k >= RAMP { 1.0 } else { 0.5 * (1.0 - (std::f64::consts::PI * k as f64 / RAMP as f64).cos()) };
            let expected = INITIAL + share * (SPEED - INITIAL);
            let (_, u) = flow.moments(inlet).unwrap();
            assert!((u[0] - expected).abs() <= tolerance, "after {k} steps: imposed {} against {expected}", u[0]);
        }
    }
}

#[test]
fn the_ramp_counts_seconds_of_the_time_step() {
    // A time step of 0.5 s: a 25 s ramp is 50 steps, and after 25 steps the share is ramp(25, 50).
    let mut p = collide(50.0, 0.5, None, Collision::Trt);
    if let Procedure::Fluid { inlet_ramp, .. } = &mut p {
        *inlet_ramp = Some(25.0);
    }
    let m = model(1.0, 0.1 / 3.0, tube(CELLS, SPEED), p);
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &[]).unwrap();
    assert_eq!(flow.inlet_ramp_steps(), 50);
    let full = flow.inlet_targets();
    for _ in 0..25 {
        flow.step().unwrap();
    }
    let share = ramp(25, 50);
    for (oc, target) in flow.open_cells().iter().zip(&full) {
        if let (Target::Velocity(u), Target::Velocity(v)) = (oc.target, target) {
            for c in 0..3 {
                assert_eq!(u[c], share * v[c], "component {c}");
            }
        }
    }
}

#[test]
fn without_a_ramp_nothing_changes() {
    let m = model(1.0, 0.1 / 3.0, tube(CELLS, SPEED), procedure(10, Collision::Trt, None));
    let setup = Setup::new(&m, None, &[], &Limits::default()).unwrap();
    let mut flow = Flow::new(&setup, &[]).unwrap();
    assert_eq!(flow.inlet_ramp_steps(), 0);
    let before: Vec<Target> = flow.open_cells().iter().map(|c| c.target).collect();
    for _ in 0..10 {
        flow.step().unwrap();
    }
    let after: Vec<Target> = flow.open_cells().iter().map(|c| c.target).collect();
    assert_eq!(format!("{before:?}"), format!("{after:?}"));
}
