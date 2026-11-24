//! A recorded wall force is the mean of two consecutive steps (Ladd 1994): the history of a square cylinder
//! near tau+ = 1/2 does not alternate with the parity of the sampled step.
mod common;

use common::*;
use fieldcraft_lbm::{Limits, solve_fluid};
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile, Turbulence};

/// The parity check's tolerance: |odd − even| as a share of the mean. Odd and even entries of the pair-mean history
/// are one step apart; for a smooth signal their means differ by the drift over the window, which is under 1/20 of
/// the mean. A single-step history gives about 0.9 of the mean here (the test's failure output).
const PARITY_TOLERANCE: f64 = 0.05;

/// The drag (SI force, x component summed over the body's faces) at every step of a square cylinder (8 cells across,
/// Re 8000, lattice speed 0.05: tau+ 0.50015) in a slip-walled box, `steps` steps.
fn drag_history(collision: Collision, steps: u64) -> Vec<f64> {
    let dx = 1.0 / 8.0;
    let mut d = domain([0.0; 3], [12.0, 8.0, dx], dx);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Convective, sigma: None };
    d.faces.y_min = FaceCondition::Slip {};
    d.faces.y_max = FaceCondition::Slip {};
    d.faces.z_min = FaceCondition::Periodic {};
    d.faces.z_max = FaceCondition::Periodic {};
    d.obstacles = vec!["Body".into()];
    d.initial_velocity = [1.0, 0.0, 0.0];
    let body = obstacle("Body", cuboid([3.5, 3.5, -dx], [4.5, 4.5, 2.0 * dx]), &CUBOID_FACES);
    let p = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: None,
        steps: Some(steps),
        time_step: Some(0.05 * dx),
        frames: 1,
        magic: None,
        sigma: None,
        collision: Some(collision),
        turbulence: Some(Turbulence::Smagorinsky { constant: 0.1, prandtl: 0.9, van_driest: None, wall_model: None }),
        steady: None,
        walls: Default::default(),
        inlet_ramp: None,
        trace_bulk_one: false,
    };
    let m = model(1.0, 1.0 / 8000.0, d, p);
    let r = solve_fluid(&m, None, &[body], &Limits::default()).unwrap();
    let n = r.history_times.len();
    (0..n).map(|i| r.forces.iter().filter(|h| h.body.as_deref() == Some("Body")).map(|h| h.force[i][0]).sum()).collect()
}

/// A wall force is recorded as the mean of two consecutive steps (Ladd 1994): the entries of a step-by-step history
/// do not depend on whether the step is odd or even. On single-step readings the two halves differ by about the
/// mean itself.
#[test]
fn the_force_history_does_not_depend_on_the_parity_of_the_sampled_step() {
    let steps = 3000;
    let f = drag_history(Collision::Hrr, steps);
    assert_eq!(f.len(), steps as usize, "one history entry per step (the parity statistic reads consecutive steps)");
    let last = &f[f.len() - 1000..];
    let mean = last.iter().sum::<f64>() / 1000.0;
    let odd = last.iter().skip(1).step_by(2).sum::<f64>() / 500.0;
    let even = last.iter().step_by(2).sum::<f64>() / 500.0;
    assert!((odd - even).abs() <= PARITY_TOLERANCE * mean.abs(), "odd {odd} even {even} mean {mean}");
}
