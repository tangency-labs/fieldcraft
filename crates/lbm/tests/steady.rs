//! Runs to steady state: plane Poiseuille flow and the de Vahl Davis cavity stop on their own and
//! match a fixed long run; an unstable flow and a divergence limit stop as diverged; faces a
//! monitor cannot find are refused.

mod common;

use common::*;
use fieldcraft_lbm::{Error, FluidResult, Limits, StopReason, solve_fluid};
use fieldcraft_model::{BoxFace, FaceCondition, FlowMonitor, FluidBoundary, FluidThermal, Model, Procedure, Steady, ThermalCondition};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

fn rule(tolerance: f64) -> Steady {
    Steady { monitors: vec![], tolerance, window: None, min_time: None, divergence: None, averaging: None }
}

fn fluid(end_time: f64, time_step: f64, steady: Option<Steady>) -> Procedure {
    Procedure::Fluid {
        domain: "Flow".into(),
        end_time: Some(end_time),
        steps: None,
        time_step: Some(time_step),
        frames: 4,
        magic: None,
        collision: None,
        turbulence: None,
        steady,
    }
}

fn solve(m: &Model) -> Result<FluidResult, Error> {
    solve_fluid(m, None, &[], &Limits::default())
}

/// Plane Poiseuille flow between walls at y = 0 and 1 driven by a body force, periodic in x and
/// z, 16 cells across, τ = 1: the peak velocity 0.1 settles on the viscous time h² / (π² ν).
fn poiseuille(steady: Option<Steady>) -> Model {
    poiseuille_for(4.0 * 6.0, steady)
}

fn poiseuille_for(end_time: f64, steady: Option<Steady>) -> Model {
    let (h, nu, umax) = (1.0, 1.0 / 6.0, 0.1);
    let dx = h / 16.0;
    let mut d = domain([0.0; 3], [dx, h, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d.acceleration = [8.0 * nu * umax / (h * h), 0.0, 0.0];
    model(1.0, nu, d, fluid(end_time, dx * dx, steady))
}

/// The largest difference between the last frames' velocities, over the largest speed.
fn velocity_difference(a: &FluidResult, b: &FluidResult) -> f64 {
    let (fa, fb) = (a.frames.last().unwrap(), b.frames.last().unwrap());
    let peak = fb.velocity.iter().map(|v| v[0].abs()).fold(0.0, f32::max);
    fa.velocity.iter().zip(&fb.velocity).map(|(u, v)| (0..3).map(|k| (u[k] - v[k]).abs()).fold(0.0, f32::max)).fold(0.0, f32::max) as f64
        / peak as f64
}

#[test]
fn poiseuille_flow_stops_at_steady_state_with_the_long_run_s_profile() {
    let fixed = solve(&poiseuille(None)).unwrap();
    assert!(fixed.steady.is_none());
    let walls = || {
        let mut r = rule(1e-5);
        r.monitors = vec![FlowMonitor::Force { faces: vec![FluidBoundary::Box(BoxFace::YMin), FluidBoundary::Box(BoxFace::YMax)] }];
        r
    };
    let run = solve(&poiseuille(Some(walls()))).unwrap();
    let report = run.steady.as_ref().unwrap();
    println!(
        "Poiseuille: steady after {} of {} steps ({:.0}% of the time), {} samples; velocity difference {:.2e}",
        run.stats.steps,
        fixed.stats.steps,
        100.0 * run.stats.steps as f64 / fixed.stats.steps as f64,
        report.samples.len(),
        velocity_difference(&run, &fixed)
    );
    assert_eq!(report.reason, StopReason::Steady, "{report:?}");
    assert_eq!(report.quantities, ["velocity", "force y_min, y_max"]);
    assert!(run.stats.steps < fixed.stats.steps / 2, "it saves at least half the time: {} of {}", run.stats.steps, fixed.stats.steps);
    assert_eq!(run.frames.last().unwrap().step, run.stats.steps, "a frame at the stop");
    assert!(run.frames.windows(2).all(|w| w[0].step < w[1].step));
    assert!(velocity_difference(&run, &fixed) < 1e-3, "{}", velocity_difference(&run, &fixed));
    // The walls take the driving force ρ a V between them: the monitored total is that.
    let last = report.samples.last().unwrap();
    let drive = 8.0 / 6.0 * 0.1 * (1.0 / 16.0) * (1.0 / 16.0);
    assert!((last.values[1][0] / drive - 1.0).abs() < 1e-3, "{last:?}");
    assert!(last.changes.iter().all(|c| c.unwrap() < 1e-5));
    // The same decisions every run.
    assert_eq!(solve(&poiseuille(Some(walls()))).unwrap(), run);
}

/// The de Vahl Davis cavity at Ra 10³ on 16 cells, as in the thermal tests: hot wall at 1 on
/// x_min, cold wall at 0 on x_max, Pr 0.71.
fn cavity(steady: Option<Steady>) -> Model {
    let (h, pr, ra, u_max) = (1.0, 0.71, 1e3, 3.65);
    let dx = h / 16.0;
    let (alpha, nu) = (1.0, pr);
    let mut d = domain([0.0; 3], [h, h, dx], dx);
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    let mut t =
        FluidThermal { initial_temperature: kelvin(0.5), gravity: [0.0, -ra * nu * alpha, 0.0], expansion: Some(1.0), ..FluidThermal::default() };
    t.faces.x_min = ThermalCondition::Temperature { value: kelvin(1.0) };
    t.faces.x_max = ThermalCondition::Temperature { value: kelvin(0.0) };
    d.thermal = Some(t);
    let dt = (0.05 * dx / u_max).min(dx * dx / 4.0);
    hot_model(1.0, nu, alpha, d, fluid(4.0, dt, steady))
}

fn hot_wall(r: &FluidResult) -> f64 {
    *r.heat.iter().find(|w| w.name == "x_min").unwrap().power.last().unwrap()
}

#[test]
fn the_de_vahl_davis_cavity_stops_at_steady_state_with_the_long_run_s_heat() {
    let fixed = solve(&cavity(None)).unwrap();
    let mut r = rule(1e-4);
    r.monitors = vec![
        FlowMonitor::Heat { faces: vec![FluidBoundary::Box(BoxFace::XMin)] },
        FlowMonitor::Heat { faces: vec![FluidBoundary::Box(BoxFace::XMax)] },
        FlowMonitor::MeanTemperature {},
    ];
    r.window = Some(0.05);
    let run = solve(&cavity(Some(r))).unwrap();
    let report = run.steady.as_ref().unwrap();
    let (q, q0) = (hot_wall(&run), hot_wall(&fixed));
    println!(
        "cavity: steady after {} of {} steps ({:.0}% of the time); hot wall {q:.6e} vs {q0:.6e}; velocity difference {:.2e}",
        run.stats.steps,
        fixed.stats.steps,
        100.0 * run.stats.steps as f64 / fixed.stats.steps as f64,
        velocity_difference(&run, &fixed)
    );
    assert_eq!(report.reason, StopReason::Steady, "{report:?}");
    assert_eq!(report.quantities, ["velocity", "temperature", "heat x_min", "heat x_max", "mean_temperature"]);
    assert!(run.stats.steps < fixed.stats.steps / 2, "{} of {}", run.stats.steps, fixed.stats.steps);
    assert!((q / q0 - 1.0).abs() < 1e-3, "{q} vs {q0}");
    assert!(velocity_difference(&run, &fixed) < 1e-2, "{}", velocity_difference(&run, &fixed));
    // The mean temperature of the symmetric cavity stays at ½, up to the lattice's own breaking
    // of the symmetry T → 1 − T through its weakly compressible flow, which falls with the grid:
    // 1.20e-3 low on these 16 cells, 3.8e-4 on 32 (the thermal lattice's Λ_T = 1/12; with a single
    // relaxation time, at this τ_T 1.38, 4.5e-4 and 1.3e-4).
    let mean = report.samples.last().unwrap().values[4][0];
    assert!((mean - kelvin(0.5)).abs() < 1.5e-3, "{mean}");
}

/// A periodic box under a body force with no walls accelerates without bound: the stability
/// check stops it, and with a stop rule the run ends as diverged with the frames saved so far.
#[test]
fn an_unstable_flow_stops_as_diverged() {
    let dx = 0.1;
    let mut d = domain([0.0; 3], [dx; 3], dx);
    for f in BoxFace::ALL {
        *d.faces.get_mut(f) = PERIODIC;
    }
    d.acceleration = [1.0, 0.0, 0.0];
    // a Δt² / Δx = 0.009: below the acceleration limit; the lattice speed passes 0.5 after 56 steps.
    let dt = (0.009 * dx).sqrt();
    let mut r = rule(1e-6);
    r.window = Some(10.0 * dt);
    let unstable = |steady| {
        model(
            1.0,
            0.01,
            d.clone(),
            Procedure::Fluid {
                domain: "Flow".into(),
                end_time: None,
                steps: Some(200),
                time_step: Some(dt),
                frames: 20,
                magic: None,
                collision: None,
                turbulence: None,
                steady,
            },
        )
    };
    assert!(matches!(solve(&unstable(None)), Err(Error::Unstable(_))), "without a rule it is an error");
    let run = solve(&unstable(Some(r.clone()))).unwrap();
    let report = run.steady.as_ref().unwrap();
    assert_eq!(report.reason, StopReason::Diverged);
    assert!(report.message.as_deref().unwrap().contains("speed"), "{report:?}");
    assert!(run.stats.steps < 60 && run.stats.steps >= 50, "{}", run.stats.steps);
    assert!(run.frames.iter().all(|f| f.step <= run.stats.steps && f.velocity.iter().all(|v| v[0].is_finite())));
    assert!(run.warnings.iter().any(|w| w.contains("diverged")), "{:?}", run.warnings);

    // A divergence limit stops it sooner: the speed keeps growing by more than the limit.
    r.divergence = Some(0.05);
    r.min_time = Some(20.0 * dt);
    let run = solve(&unstable(Some(r))).unwrap();
    let report = run.steady.as_ref().unwrap();
    assert_eq!(report.reason, StopReason::Diverged);
    assert_eq!(run.stats.steps, 20, "the first sample at or after min_time: {report:?}");
    assert!(report.message.as_deref().unwrap().contains("divergence limit"), "{report:?}");
}

#[test]
fn a_flow_that_does_not_settle_runs_its_length_and_says_so() {
    // Two time units: about three viscous times h² / (π² ν), still settling.
    let mut r = rule(1e-3);
    r.window = Some(0.125);
    let run = solve(&poiseuille_for(2.0, Some(r))).unwrap();
    let report = run.steady.as_ref().unwrap();
    assert_eq!(report.reason, StopReason::MaxTime, "{report:?}");
    assert_eq!(run.stats.steps, solve(&poiseuille_for(2.0, None)).unwrap().stats.steps);
    assert!(run.warnings.iter().any(|w| w.contains("did not reach steady state")), "{:?}", run.warnings);
    // 2 time units in windows of ⅛, and the initial sample.
    assert_eq!(report.samples.len(), 17);
    assert!(report.samples[0].changes.iter().all(Option::is_none));
}

#[test]
fn monitors_the_flow_cannot_measure_are_refused() {
    let with = |monitors: Vec<FlowMonitor>| {
        let mut r = rule(1e-3);
        r.monitors = monitors;
        solve(&poiseuille(Some(r)))
    };
    let face = |body: &str, face: &str| FluidBoundary::Surface(fieldcraft_model::GeometryRef::face(body, face).unwrap());
    assert!(matches!(with(vec![FlowMonitor::Force { faces: vec![face("Ball", "skin")] }]), Err(Error::NotFound(_))));
    assert!(matches!(with(vec![FlowMonitor::Heat { faces: vec![FluidBoundary::Box(BoxFace::YMin)] }]), Err(Error::Invalid(_))));
    assert!(matches!(with(vec![FlowMonitor::MeanTemperature {}]), Err(Error::Invalid(_))));
    // A mass flow through a periodic face is the flow along the channel: fine.
    let run = with(vec![FlowMonitor::MassFlow { face: BoxFace::XMax }, FlowMonitor::MassFlow { face: BoxFace::XMin }]).unwrap();
    let last = run.steady.unwrap().samples.pop().unwrap();
    // ρ ū H Δx with the mean of the parabola ⅔ of its peak 0.1 (cell-averaged).
    let expected = 2.0 / 3.0 * 0.1 * 1.0 / 16.0;
    assert!((last.values[1][0] / expected - 1.0).abs() < 1e-2, "{last:?}");
    assert!((last.values[1][0] + last.values[2][0]).abs() < 1e-12, "{last:?}");
}
