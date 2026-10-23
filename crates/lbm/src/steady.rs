//! Running a fluid step to steady state: the monitored quantities and the stop rule
//! ([`fieldcraft_model::Steady`]).
//!
//! The run samples every monitored quantity [`SUB_SAMPLES`] times per `window` time steps, at
//! the same chunk boundaries on every backend (a device backend reads its fields back a few times
//! per window, never per step), so CPU and GPU take their decisions at the same steps. Each
//! quantity `q` (a whole field for the residuals, a vector or a number otherwise) has, over each
//! window, the relative spread `‖max q − min q‖ / max ‖q‖`: the componentwise extremes over the
//! window's samples (both ends included), the Euclidean norm, and the largest norm over every
//! sample so far; it is 0 while `q` has always been 0. The spread, unlike the difference of the
//! window's ends, sees an oscillation whose period divides the window. Temperatures count from
//! the initial temperature, so an offset (kelvin, say) does not hide a change.
//!
//! The run stops as steady at the end of the second consecutive full window, at or after
//! `min_time`, in which every relative spread is below the tolerance; as diverged when one
//! exceeds the divergence limit after `min_time`, is not a finite number, or the backend reports
//! the flow unstable. A last window cut short by the end of the step takes no decision.
//!
//! With averaging ([`fieldcraft_model::Averaging`]), for a turbulent flow that never stops
//! fluctuating, the decision is taken on each quantity's running mean of its reported value (the
//! root mean square speed for the velocity field) from the averaging start on: its relative
//! spread over a window counts only for a window that begins at or after the start, so the
//! start-up never enters a decision. The report adds each quantity's mean, its fluctuation (the
//! root mean square about the mean) and its extremes over the samples averaged.

use fieldcraft_model::{BoxFace, FlowMonitor, FluidBoundary, GeometryRef, Steady};

use crate::backend::{Backend, Fields};
use crate::setup::Setup;
use crate::{Error, Result};

/// Why a run to steady state stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum StopReason {
    /// Every monitored quantity stopped changing.
    Steady,
    /// The step's full length passed first.
    MaxTime,
    /// A monitored quantity blew up, or the flow went unstable.
    Diverged,
}

impl StopReason {
    pub fn name(self) -> &'static str {
        match self {
            StopReason::Steady => "steady",
            StopReason::MaxTime => "max_time",
            StopReason::Diverged => "diverged",
        }
    }
}

/// Samples per window.
pub const SUB_SAMPLES: u64 = 8;

/// Consecutive steady windows before the run stops as steady.
pub const STEADY_WINDOWS: u32 = 2;

/// The stop rule in time steps.
#[derive(Debug, Clone, PartialEq)]
pub struct SteadySetup {
    pub rule: Steady,
    /// Time steps per window: a whole number of samples.
    pub window: u64,
    /// Time steps between samples.
    pub sample: u64,
    /// The first step a decision may be taken at.
    pub min_steps: u64,
    /// With averaging: the first step whose sample enters the running means.
    pub average_from: Option<u64>,
}

impl SteadySetup {
    /// The rule for a step of `steps` time steps of `dt`. The window is at least one step and
    /// one [`crate::solve::MAX_HISTORY`]-th of the step, so the history stays bounded.
    pub fn new(rule: &Steady, steps: u64, dt: f64) -> Result<SteadySetup> {
        rule.validate()?;
        let floor = steps.div_ceil(crate::solve::MAX_HISTORY).max(1);
        let window = match rule.window {
            Some(w) => {
                let n = (w / dt).round();
                if !n.is_finite() {
                    return Err(Error::Invalid(format!("the steady-state window {w} is not a number of time steps of {dt}")));
                }
                // Saturating: a window beyond the step's length samples only at its end.
                (n as u64).clamp(1, steps.max(1))
            }
            None => steps.div_ceil(20).max(1),
        }
        .max(floor);
        let sample = (window / SUB_SAMPLES).max(1);
        let window = sample * (window / sample);
        let min_steps = rule.min_time.map_or(0, |t| ((t / dt).ceil().max(0.0) as u64).min(steps));
        let average_from = rule.averaging.map(|a| ((a.start / dt).ceil().max(0.0) as u64).min(steps));
        Ok(SteadySetup { rule: rule.clone(), window, sample, min_steps, average_from })
    }
}

/// One sample of every monitored quantity.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MonitorSample {
    pub step: u64,
    pub time: f64,
    /// Per quantity: its value in model units (the root mean square over the fluid cells for a
    /// residual's field, the vector for a force, the number otherwise).
    pub values: Vec<Vec<f64>>,
    /// Per quantity: its relative spread over the window that ends here (`None` at the start).
    pub changes: Vec<Option<f64>>,
}

/// How a run to steady state went.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SteadyReport {
    pub reason: StopReason,
    pub tolerance: f64,
    pub window: f64,
    /// The monitored quantities: `velocity`, `temperature` (with heat), then the rule's
    /// monitors by [`FlowMonitor::label`].
    pub quantities: Vec<String>,
    pub samples: Vec<MonitorSample>,
    /// Why it diverged.
    pub message: Option<String>,
    /// With averaging: each quantity's statistics over the samples averaged (in the order of
    /// `quantities`; empty when no sample was averaged).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub averages: Vec<Average>,
}

/// A quantity's statistics over the samples a statistically steady run averaged.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Average {
    /// The flow time the averaging began at.
    pub from: f64,
    pub samples: u64,
    /// Componentwise, in model units: the mean, the root mean square about it, and the extremes.
    pub mean: Vec<f64>,
    pub rms: Vec<f64>,
    pub min: Vec<f64>,
    pub max: Vec<f64>,
}

#[derive(Debug, Clone)]
enum Kind {
    Velocity,
    Temperature,
    MassFlow(BoxFace),
    Force(Vec<usize>),
    Heat(Vec<usize>),
    MeanTemperature,
}

#[derive(Debug, Clone)]
struct Quantity {
    kind: Kind,
    /// The componentwise extremes over the current window's samples so far.
    low: Vec<f64>,
    high: Vec<f64>,
    peak: f64,
    /// With averaging: the running sums of the reported value and its square, its extremes, and
    /// the samples averaged.
    sum: Vec<f64>,
    squares: Vec<f64>,
    least: Vec<f64>,
    most: Vec<f64>,
    count: u64,
}

/// The monitor of a run: its quantities and their history.
#[derive(Debug, Clone)]
pub struct Monitor {
    pub setup: SteadySetup,
    /// Consecutive steady windows so far.
    steady_windows: u32,
    quantities: Vec<Quantity>,
    report: SteadyReport,
}

/// The force groups a monitor's faces name; every face must name at least one.
fn groups(flow: &dyn Backend, faces: &[FluidBoundary]) -> Result<Vec<usize>> {
    let mut out = Vec::new();
    for face in faces {
        let before = out.len();
        for (i, g) in flow.groups().iter().enumerate() {
            let hit = match (face, &g.body) {
                (FluidBoundary::Box(b), None) => g.name == b.name(),
                (FluidBoundary::Surface(r), Some(body)) => GeometryRef::face(body.as_str(), g.name.as_str()).is_ok_and(|tag| r.selects(&tag)),
                _ => false,
            };
            if hit && !out.contains(&i) {
                out.push(i);
            }
        }
        if out.len() == before {
            let name = match face {
                FluidBoundary::Box(b) => format!("box face {}", b.name()),
                FluidBoundary::Surface(r) => r.to_string(),
            };
            return Err(Error::NotFound(format!("the monitored {name} is not a wall of the domain")));
        }
    }
    Ok(out)
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

impl Monitor {
    /// The monitor of `setup`'s rule on `flow`, its faces resolved; `None` without a rule.
    pub fn new(setup: &Setup, flow: &dyn Backend) -> Result<Option<Monitor>> {
        let Some(steady) = &setup.steady else { return Ok(None) };
        let thermal = setup.thermal.is_some();
        let mut labels = vec!["velocity".to_string()];
        let mut kinds = vec![Kind::Velocity];
        if thermal {
            labels.push("temperature".into());
            kinds.push(Kind::Temperature);
        }
        for m in &steady.rule.monitors {
            let kind = match m {
                FlowMonitor::MassFlow { face } => Kind::MassFlow(*face),
                FlowMonitor::Force { faces } => Kind::Force(groups(flow, faces)?),
                FlowMonitor::Heat { faces } => Kind::Heat(groups(flow, faces)?),
                FlowMonitor::MeanTemperature {} => Kind::MeanTemperature,
            };
            if matches!(kind, Kind::Heat(_) | Kind::MeanTemperature) && !thermal {
                return Err(Error::Invalid(format!("the monitored {} needs heat in the fluid domain", m.label())));
            }
            labels.push(m.label());
            kinds.push(kind);
        }
        Ok(Some(Monitor {
            quantities: kinds
                .into_iter()
                .map(|kind| Quantity {
                    kind,
                    low: Vec::new(),
                    high: Vec::new(),
                    peak: 0.0,
                    sum: Vec::new(),
                    squares: Vec::new(),
                    least: Vec::new(),
                    most: Vec::new(),
                    count: 0,
                })
                .collect(),
            report: SteadyReport {
                reason: StopReason::MaxTime,
                tolerance: steady.rule.tolerance,
                window: steady.window as f64 * setup.time_step,
                quantities: labels,
                samples: Vec::new(),
                message: None,
                averages: Vec::new(),
            },
            setup: steady.clone(),
            steady_windows: 0,
        }))
    }

    /// Whether step `k` is a sample point.
    pub fn samples_at(&self, k: u64) -> bool {
        k.is_multiple_of(self.setup.sample)
    }

    /// Sample every quantity at step `k` (a sample point) and, at the end of a window, decide:
    /// `Some` reason to stop now.
    pub fn sample(&mut self, setup: &Setup, flow: &mut dyn Backend, k: u64) -> Result<Option<StopReason>> {
        let Fields { density, velocity, temperature } = flow.fields()?;
        let fluid: Vec<usize> = (0..setup.grid.len()).filter(|c| flow.is_fluid(*c)).collect();
        let n = fluid.len().max(1) as f64;
        let vs = setup.velocity_scale();
        let t0 = setup.thermal.as_ref().map_or(0.0, |t| t.initial_temperature);
        let zero = setup.thermal.as_ref().map_or(0.0, |t| t.zero);
        let forces = if self.quantities.iter().any(|q| matches!(q.kind, Kind::Force(_))) { flow.forces()? } else { Vec::new() };
        let heat = if self.quantities.iter().any(|q| matches!(q.kind, Kind::Heat(_))) { flow.heat()? } else { Vec::new() };
        let end = k.is_multiple_of(self.setup.window);
        let mut values = Vec::with_capacity(self.quantities.len());
        let mut changes = Vec::with_capacity(self.quantities.len());
        for q in &mut self.quantities {
            // The monitored vector (temperatures from the initial one) and the reported value.
            let (state, value): (Vec<f64>, Vec<f64>) = match &q.kind {
                Kind::Velocity => {
                    let s: Vec<f64> = fluid.iter().flat_map(|c| velocity.get(*c).copied().unwrap_or([0.0; 3]).map(|u| u * vs)).collect();
                    let rms = norm(&s) / n.sqrt();
                    (s, vec![rms])
                }
                Kind::Temperature => {
                    // the fluid's cells and the conducting solids' (the slowest to settle)
                    let cells = fluid.iter().copied().chain(flow.conducting().iter().flatten().map(|c| *c as usize));
                    let s: Vec<f64> = cells.map(|c| temperature.get(c).copied().unwrap_or(t0) - t0).collect();
                    let rms = norm(&s) / (s.len().max(1) as f64).sqrt();
                    (s, vec![rms])
                }
                Kind::MeanTemperature => {
                    let mean = fluid.iter().map(|c| temperature.get(*c).copied().unwrap_or(t0) - t0).sum::<f64>() / n;
                    // The lattice is above `zero`; the report in kelvin.
                    (vec![mean], vec![fieldcraft_model::units::kelvin_from(mean + t0, zero)])
                }
                Kind::MassFlow(face) => {
                    let axis = face.axis();
                    let (layer, sign) = match face {
                        BoxFace::XMin | BoxFace::YMin | BoxFace::ZMin => (0, -1.0),
                        _ => (setup.grid.cells.get(axis).copied().unwrap_or(1).saturating_sub(1), 1.0),
                    };
                    let mut flux = 0.0;
                    for &c in &fluid {
                        if setup.grid.coords(c).get(axis) == Some(&layer) {
                            let rho = density.get(c).copied().unwrap_or(1.0);
                            let u = velocity.get(c).and_then(|u| u.get(axis).copied()).unwrap_or(0.0);
                            flux += sign * rho * u;
                        }
                    }
                    let m = flux * vs * setup.density * setup.grid.spacing * setup.grid.spacing;
                    (vec![m], vec![m])
                }
                Kind::Force(ids) => {
                    let fs = setup.force_scale();
                    let mut f = [0.0; 3];
                    for i in ids {
                        if let Some(g) = forces.get(*i) {
                            for a in 0..3 {
                                f[a] += g[a] * fs;
                            }
                        }
                    }
                    (f.to_vec(), f.to_vec())
                }
                Kind::Heat(ids) => {
                    let hs = setup.thermal.as_ref().map_or(0.0, |t| t.heat_scale);
                    let p: f64 = ids.iter().filter_map(|i| heat.get(*i)).map(|h| h * hs).sum();
                    (vec![p], vec![p])
                }
            };
            // With averaging, what must settle is the running mean of the reported value, over
            // windows that begin at or after the averaging start.
            let (state, counts) = match self.setup.average_from {
                Some(from) if k >= from => {
                    if q.count == 0 {
                        q.sum = vec![0.0; value.len()];
                        q.squares = vec![0.0; value.len()];
                        q.least.clone_from(&value);
                        q.most.clone_from(&value);
                    }
                    for (i, x) in value.iter().enumerate() {
                        if let (Some(s), Some(s2), Some(lo), Some(hi)) =
                            (q.sum.get_mut(i), q.squares.get_mut(i), q.least.get_mut(i), q.most.get_mut(i))
                        {
                            *s += x;
                            *s2 += x * x;
                            *lo = lo.min(*x);
                            *hi = hi.max(*x);
                        }
                    }
                    q.count += 1;
                    let n = q.count as f64;
                    let mean: Vec<f64> = q.sum.iter().map(|s| s / n).collect();
                    let started = k.saturating_sub(self.setup.window) >= from;
                    // An instantaneous value that is not a number still ends the run.
                    let mean = if value.iter().any(|x| !x.is_finite()) { vec![f64::NAN; mean.len()] } else { mean };
                    (mean, started)
                }
                // Before the start only a value that is not a number counts (it ends the run).
                Some(_) => (state, false),
                None => (state, true),
            };
            if counts || self.setup.average_from.is_none() {
                q.peak = q.peak.max(norm(&state));
            }
            if q.low.len() != state.len() {
                q.low.clone_from(&state);
                q.high.clone_from(&state);
            }
            for ((lo, hi), x) in q.low.iter_mut().zip(q.high.iter_mut()).zip(&state) {
                *lo = lo.min(*x);
                *hi = hi.max(*x);
                // A NaN is never steady.
                if x.is_nan() {
                    *hi = f64::INFINITY;
                }
            }
            let blown = q.high.contains(&f64::INFINITY);
            let change = (end && k > 0 && (counts || blown)).then(|| {
                let spread = q.low.iter().zip(&q.high).map(|(a, b)| (b - a) * (b - a)).sum::<f64>().sqrt();
                if q.peak > 0.0 && q.peak.is_finite() {
                    spread / q.peak
                } else if spread == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                }
            });
            if end {
                // The next window starts from this sample.
                q.low.clone_from(&state);
                q.high = state;
            }
            values.push(value);
            changes.push(change);
        }
        if !end {
            return Ok(None);
        }
        let rule = &self.setup.rule;
        let late = k >= self.setup.min_steps;
        let mut decision = None;
        if let Some((i, c)) = changes.iter().enumerate().find_map(|(i, c)| c.filter(|c| !c.is_finite()).map(|c| (i, c))) {
            decision = Some(StopReason::Diverged);
            self.report.message = Some(format!("the {} changed by {c} over the window", self.label(i)));
        } else if let (true, Some(limit)) = (late, rule.divergence)
            && let Some((i, c)) = changes.iter().enumerate().find_map(|(i, c)| c.filter(|c| *c > limit).map(|c| (i, c)))
        {
            decision = Some(StopReason::Diverged);
            self.report.message = Some(format!("the {} changed by {c:.3e} over the window, above the divergence limit {limit:e}", self.label(i)));
        } else if changes.iter().all(|c| c.is_some_and(|c| c < rule.tolerance)) {
            self.steady_windows += 1;
            if late && self.steady_windows >= STEADY_WINDOWS {
                decision = Some(StopReason::Steady);
            }
        } else {
            self.steady_windows = 0;
        }
        self.report.samples.push(MonitorSample { step: k, time: k as f64 * setup.time_step, values, changes });
        if let Some(r) = decision {
            self.report.reason = r;
        }
        Ok(decision)
    }

    fn label(&self, i: usize) -> &str {
        self.report.quantities.get(i).map_or("quantity", String::as_str)
    }

    /// The backend stopped with an error: an unstable flow ends the run as diverged.
    pub fn unstable(&mut self, message: String) {
        self.report.reason = StopReason::Diverged;
        self.report.message = Some(message);
    }

    pub fn report(mut self) -> SteadyReport {
        if let Some(from) = self.setup.average_from {
            let start = self.setup.rule.averaging.map_or(from as f64, |a| a.start);
            self.report.averages = self
                .quantities
                .iter()
                .filter(|q| q.count > 0)
                .map(|q| {
                    let n = q.count as f64;
                    let mean: Vec<f64> = q.sum.iter().map(|s| s / n).collect();
                    let rms = q.squares.iter().zip(&mean).map(|(s2, m)| (s2 / n - m * m).max(0.0).sqrt()).collect();
                    Average { from: start, samples: q.count, mean, rms, min: q.least.clone(), max: q.most.clone() }
                })
                .collect();
        }
        self.report
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flow::Group;
    use crate::{Limits, run};
    use fieldcraft_model::{BoxFaces, Fluid, FluidDomain, Material, Model, Procedure, Step};

    /// A stand-in backend whose uniform velocity is `u(step)`.
    struct Scripted {
        k: u64,
        cells: usize,
        u: fn(u64) -> f64,
    }

    impl Backend for Scripted {
        fn name(&self) -> String {
            "cpu".into()
        }
        fn advance(&mut self, steps: u64) -> Result<()> {
            self.k += steps;
            Ok(())
        }
        fn steps_done(&self) -> u64 {
            self.k
        }
        fn fields(&mut self) -> Result<Fields> {
            Ok(Fields { density: vec![1.0; self.cells], velocity: vec![[(self.u)(self.k), 0.0, 0.0]; self.cells], temperature: Vec::new() })
        }
        fn mass(&mut self) -> Result<f64> {
            Ok(self.cells as f64)
        }
        fn forces(&mut self) -> Result<Vec<[f64; 3]>> {
            Ok(Vec::new())
        }
        fn heat(&mut self) -> Result<Vec<f64>> {
            Ok(Vec::new())
        }
        fn groups(&self) -> &[Group] {
            &[]
        }
        fn fluid_cells(&self) -> usize {
            self.cells
        }
        fn is_fluid(&self, _: usize) -> bool {
            true
        }
        fn empty_obstacles(&self) -> &[String] {
            &[]
        }
        fn conducting(&self) -> &[Vec<u32>] {
            &[]
        }
    }

    /// 2000 steps, a window of 50 steps.
    fn setup() -> Setup {
        let mut m = Model::default();
        m.add_material("Fluid", Material { density: Some(1.0), fluid: Some(Fluid::Newtonian { viscosity: 0.1 }), ..Material::default() }).unwrap();
        let domain = FluidDomain {
            material: "Fluid".into(),
            min: [0.0; 3],
            max: [1.0; 3],
            spacing: 0.5,
            faces: BoxFaces::default(),
            obstacles: vec![],
            surfaces: vec![],
            initial_velocity: [0.0; 3],
            acceleration: [0.0; 3],
            tolerance: None,
            thermal: None,
            porous: vec![],
            fans: vec![],
        };
        m.add_fluid_domain("Flow", domain).unwrap();
        let rule = Steady { monitors: vec![], tolerance: 1e-3, window: Some(50.0), min_time: None, divergence: None, averaging: None };
        let p = Procedure::Fluid {
            domain: "Flow".into(),
            end_time: None,
            steps: Some(2000),
            time_step: Some(1.0),
            frames: 1,
            magic: None,
            collision: None,
            turbulence: None,
            steady: Some(rule),
        };
        m.add_step(Step::new("Run", p), None).unwrap();
        Setup::new(&m, None, &[], &Limits::default()).unwrap()
    }

    fn stop(u: fn(u64) -> f64) -> (StopReason, u64) {
        stop_with(setup(), u)
    }

    fn stop_with(s: Setup, u: fn(u64) -> f64) -> (StopReason, u64) {
        let mut b = Scripted { k: 0, cells: s.grid.len(), u };
        let r = run(&s, &mut b, &Limits::default()).unwrap();
        (r.steady.unwrap().reason, r.stats.steps)
    }

    #[test]
    fn an_oscillation_whose_period_divides_the_window_is_not_steady() {
        use std::f64::consts::PI;
        assert_eq!(setup().steady.unwrap().sample, 6);
        // 50 steps make 8 samples of 6 steps: the window is 48 steps.
        assert_eq!(setup().steady.unwrap().window, 48);
        assert_eq!(stop(|k| 0.05 + 0.02 * (2.0 * PI * k as f64 / 48.0).sin()), (StopReason::MaxTime, 2000));
        assert_eq!(stop(|k| 0.05 + 0.02 * (2.0 * PI * k as f64 / 24.0).sin()).0, StopReason::MaxTime);
        assert_eq!(stop(|k| 0.05 + 0.02 * (2.0 * PI * k as f64 / 50.0).sin()).0, StopReason::MaxTime);
        // A constant flow is steady after two windows; a settling one once it has settled.
        assert_eq!(stop(|_| 0.05), (StopReason::Steady, 96));
        let (reason, steps) = stop(|k| 0.05 * (1.0 - (-(k as f64) / 100.0).exp()));
        assert_eq!(reason, StopReason::Steady);
        assert!((600..800).contains(&steps), "{steps}");
    }

    #[test]
    fn a_window_cut_short_by_the_end_takes_no_decision() {
        // Windows of 1496 steps in 2000: one full window, then 504 steps that do not count.
        let mut s = setup();
        let rule = s.steady.as_ref().unwrap().rule.clone();
        s.steady = Some(SteadySetup::new(&Steady { window: Some(1500.0), ..rule }, 2000, 1.0).unwrap());
        assert_eq!(stop_with(s, |_| 0.05), (StopReason::MaxTime, 2000));
    }
    fn averaged(start: f64, tolerance: f64) -> Setup {
        let mut s = setup();
        let rule = s.steady.as_ref().unwrap().rule.clone();
        let rule = Steady { averaging: Some(fieldcraft_model::Averaging { start }), tolerance, ..rule };
        s.steady = Some(SteadySetup::new(&rule, 2000, 1.0).unwrap());
        s
    }

    #[test]
    fn a_fluctuating_flow_stops_once_its_running_mean_settles() {
        use std::f64::consts::PI;
        // A start-up far from the mean, then a fluctuation that never dies out.
        let u: fn(u64) -> f64 = |k| if k < 200 { 0.5 } else { 0.05 + 0.02 * (2.0 * PI * k as f64 / 48.0).sin() };
        let s = averaged(200.0, 1e-2);
        assert_eq!(s.steady.as_ref().unwrap().average_from, Some(200));
        let mut b = Scripted { k: 0, cells: s.grid.len(), u };
        let r = run(&s, &mut b, &Limits::default()).unwrap();
        let report = r.steady.unwrap();
        assert_eq!(report.reason, StopReason::Steady, "{:?}", report.message);
        assert!((300..2000).contains(&r.stats.steps), "{}", r.stats.steps);
        // The statistics leave the start-up out.
        let a = report.averages.first().expect("the velocity's average");
        assert_eq!(a.from, 200.0);
        // In model units: a lattice speed of 1 is a cell per time step.
        let v = s.grid.spacing / s.time_step;
        assert!((a.mean[0] - 0.05 * v).abs() < 0.05 * v * 1e-2, "{:?}", a.mean);
        assert!((a.rms[0] - 0.02 * v / 2f64.sqrt()).abs() < 1e-3 * v, "{:?}", a.rms);
        assert!(a.min[0] >= 0.03 * v - 1e-9 && a.max[0] <= 0.07 * v + 1e-9, "{a:?}");
        // Without averaging the fluctuation alone never stops.
        assert_eq!(stop(|k| 0.05 + 0.02 * (2.0 * PI * k as f64 / 48.0).sin()).0, StopReason::MaxTime);
    }

    #[test]
    fn a_flow_that_blows_up_before_the_averaging_start_still_diverges() {
        let s = averaged(1000.0, 1e-2);
        let (reason, steps) = stop_with(s, |k| if k < 300 { 0.05 } else { f64::NAN });
        assert_eq!(reason, StopReason::Diverged);
        assert!(steps < 1000, "{steps}");
    }
}
