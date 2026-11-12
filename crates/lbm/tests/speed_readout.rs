//! The run reports the true largest speed over the fluid cells, and warns above Mach 0.1, whatever the
//! setup's estimate from the characteristic velocity said.
mod common;

use common::*;
use fieldcraft_lbm::{Limits, solve_fluid};
use fieldcraft_model::{Collision, FaceCondition, Outflow, Procedure, Profile};

/// A channel of 32 x 8 x 4 cells with a uniform inlet and two slabs that leave a throat a quarter of
/// the width: the inlet sets the setup's estimate (Mach 0.087), the throat runs several times faster.
fn throat(gap: f64, slabs: bool) -> (fieldcraft_model::Model, [fieldcraft_lbm::Obstacle; 2]) {
    let mut d = domain([0.0; 3], [1.6, 0.4, 0.2], 0.05);
    d.faces.x_min = FaceCondition::Velocity { velocity: [1.0, 0.0, 0.0], profile: Profile::Uniform };
    d.faces.x_max = FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform };
    d.initial_velocity = [1.0, 0.0, 0.0];
    if slabs {
        d.obstacles = vec!["Lo".into(), "Hi".into()];
    }
    let mid = 0.2;
    let obstacles = [
        obstacle("Lo", cuboid([0.6, -1.0, -1.0], [1.0, mid - gap / 2.0, 1.0]), &["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"]),
        obstacle("Hi", cuboid([0.6, mid + gap / 2.0, -1.0], [1.0, 2.0, 1.0]), &["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"]),
    ];
    let run = Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(1.0),
        steps: None,
        time_step: None,
        frames: 2,
        magic: None,
        sigma: None,
        collision: Some(Collision::Trt),
        turbulence: None,
        steady: None,
        walls: Default::default(),
    };
    (model(1.0, 0.01, d, run), obstacles)
}

#[test]
fn a_throat_faster_than_the_setup_estimate_is_reported_and_warned_about() {
    let (m, obstacles) = throat(0.1, true);
    let r = solve_fluid(&m, None, &obstacles, &Limits::default()).unwrap_or_else(|e| panic!("{e}"));
    let peak = r.stats.validity.as_ref().and_then(|v| v.peak_mach()).map_or(0.0, |s| s.mach_max);
    eprintln!("estimate {:.3}, true largest Mach {peak:.3}; warnings {:?}", r.stats.mach, r.warnings);
    assert!(r.stats.mach <= 0.1, "the setup estimate is the inlet's: {}", r.stats.mach);
    assert!(peak > 0.1, "the throat runs faster than the inlet: {peak}");
    assert!(r.warnings.iter().any(|w| w.contains("largest local Mach number")), "{:?}", r.warnings);
}

#[test]
fn a_uniform_stream_at_the_estimate_gets_no_such_warning() {
    let (mut m, _) = throat(0.1, false);
    let mut d = domain([0.0; 3], [0.8, 0.2, 0.2], 0.05);
    for face in [&mut d.faces.x_min, &mut d.faces.x_max, &mut d.faces.y_min, &mut d.faces.y_max, &mut d.faces.z_min, &mut d.faces.z_max] {
        *face = FaceCondition::Periodic {};
    }
    d.initial_velocity = [1.0, 0.0, 0.0];
    m.replace_fluid_domain("Flow", d).unwrap_or_else(|e| panic!("{e}"));
    let r = solve_fluid(&m, None, &[], &Limits::default()).unwrap_or_else(|e| panic!("{e}"));
    let peak = r.stats.validity.as_ref().and_then(|v| v.peak_mach()).map_or(0.0, |s| s.mach_max);
    assert!((peak - r.stats.mach).abs() < 1e-6, "true {peak} against estimate {}", r.stats.mach);
    assert!(!r.warnings.iter().any(|w| w.contains("largest local Mach number")), "{:?}", r.warnings);
}
