//! Analysis steps and the step × load/boundary-condition table.
//!
//! Steps run in order. Each step lists the boundary conditions and loads active in it, by name;
//! nothing carries over from one step to the next unless it is listed again, so the table a user
//! sees is the whole truth (no propagation, no hidden modified or inactive states). The state
//! before the first step is the unloaded, unconstrained model; a boundary condition active in
//! every step can be written at model level by formats that have one.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::{Error, Friction, Model, Result, check_name, context, finite, positive, require};

/// Most steps in a model.
pub const MAX_STEPS: usize = 1000;

/// Most modes a frequency step extracts.
pub const MAX_MODES: u32 = 100_000;

/// Most time increments an explicit dynamic step may ask for.
pub const MAX_INCREMENTS: u64 = 1_000_000_000;

/// Most field frames an explicit dynamic step writes (end time over output interval).
pub const MAX_FRAMES: f64 = 10_000.0;

/// Most history samples an explicit dynamic step records (end time over history interval).
pub const MAX_HISTORY_SAMPLES: f64 = 1_000_000.0; // unit-factors: a count of samples

/// Most initial velocities or rigid walls an explicit dynamic step lists.
pub const MAX_STEP_ITEMS: usize = 10_000;

/// Most time steps a fluid step asks for (the solver caps a run lower by default).
pub const MAX_FLUID_STEPS: u64 = 1_000_000_000;

/// Most output frames a fluid step writes.
pub const MAX_FLUID_FRAMES: u32 = 1000;

fn default_frames() -> u32 {
    10
}

/// What a step solves.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum Procedure {
    /// Linear static equilibrium under the active loads and boundary conditions.
    Static {},
    /// Natural frequencies and mode shapes: the lowest `modes` modes, optionally only those
    /// between `min_hz` and `max_hz` (cycles per time unit of the unit system).
    Frequency {
        modes: u32,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        min_hz: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max_hz: Option<f64>,
    },
    /// Linear buckling: the load factors λ at which the step's loads (and prescribed motions),
    /// applied as a linear preload, make K + λ K_σ singular, and their mode shapes; the lowest
    /// `modes` positive ones, or those nearest above `shift` when one is given.
    Buckle {
        modes: u32,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        shift: Option<f64>,
    },
    /// Transient dynamics by explicit time integration up to `end_time`.
    ExplicitDynamic(ExplicitDynamic),
    /// Time-dependent incompressible flow in a fluid domain, from rest (or the domain's initial
    /// velocity) to `end_time`, or for `steps` time steps of `time_step` (exactly one of the two
    /// lengths). The solver picks the time step when `time_step` is left out. `frames` velocity
    /// and pressure fields are written, evenly spaced, after the initial one; forces on every
    /// wall are recorded throughout. `collision` picks the lattice Boltzmann collision operator
    /// ([`Collision`], two relaxation times by default); `magic` is the two-relaxation-time
    /// parameter Λ (default 3/16, which puts bounce-back walls exactly halfway between cells).
    /// With `steady`, the run stops on its own once the flow stops changing ([`Steady`]); its
    /// length is then the most it may run. `turbulence` adds a subgrid eddy viscosity
    /// ([`Turbulence`]; none by default).
    Fluid {
        domain: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        end_time: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        steps: Option<u64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        time_step: Option<f64>,
        #[serde(default = "default_frames")]
        frames: u32,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        magic: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        collision: Option<Collision>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        turbulence: Option<Turbulence>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        steady: Option<Steady>,
    },
    /// Heat conduction: the steady temperature field, or its history over time when `transient`
    /// is given. Temperature boundary conditions and heat loads (flux, convection, sources)
    /// act; mechanical ones are ignored.
    HeatTransfer {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        transient: Option<Transient>,
    },
    /// Steady heat conduction and the static equilibrium of the thermal strains it causes,
    /// together: every boundary condition and load of the step acts.
    CoupledTemperatureDisplacement {},
    /// Static equilibrium followed in load increments with Newton–Raphson iterations: large
    /// displacements, rotations and strains (`nonlinear_geometry`), hyperelastic and plastic
    /// materials and follower pressure. Loads and prescribed displacements are the step's load
    /// factor times their values: a ramp from 0 to 1, or the step's load path.
    NonlinearStatic(NonlinearStatic),
}

/// The collision operator of a lattice Boltzmann fluid step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Collision {
    /// Two relaxation times: the symmetric and antisymmetric parts of each pair of opposite
    /// populations relax separately (the default).
    #[default]
    Trt,
    /// Recursive regularised: before relaxing, the non-equilibrium part is rebuilt from the
    /// momentum flux alone (and the third-order moments it implies), which removes the
    /// non-hydrodynamic modes that make the two-relaxation-time scheme unstable as τ approaches ½.
    Regularized,
}

impl Collision {
    /// `trt` or `regularized`, as written in the model.
    pub fn name(self) -> &'static str {
        match self {
            Collision::Trt => "trt",
            Collision::Regularized => "regularized",
        }
    }
}

/// A fluid step's turbulence model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Turbulence {
    /// Smagorinsky's subgrid eddy viscosity ν_t = (C Δx)² |S|, with `constant` C (0.1 to 0.2 is
    /// usual; 0.17 is Lilly's for isotropic turbulence) and |S| = √(2 S:S) the local strain
    /// rate. With heat, the turbulent diffusivity α_t = ν_t / Pr_t adds to the fluid's, with
    /// `prandtl` the turbulent Prandtl number Pr_t (0.85 to 0.9 is usual; default 0.9).
    Smagorinsky {
        constant: f64,
        #[serde(default = "default_turbulent_prandtl")]
        prandtl: f64,
    },
}

fn default_turbulent_prandtl() -> f64 {
    0.9
}

/// Most quantities a steady-state rule monitors besides the residual.
pub const MAX_FLOW_MONITORS: usize = 100;

/// Most faces one monitored quantity sums over.
pub const MAX_MONITOR_FACES: usize = 1000;

/// The relative change per window below which a flow is steady, when none is given.
const STEADY_TOLERANCE: f64 = 1e-3; // unit-factors: a relative change

fn default_steady_tolerance() -> f64 {
    STEADY_TOLERANCE
}

/// The stop rule of a fluid step run to steady state.
///
/// The run samples its monitored quantities 8 times per `window` of flow time: the velocity
/// field (and the temperature field when the flow carries heat) and every quantity of
/// `monitors`. `tolerance` is the relative change per window: the spread (largest minus smallest
/// sample) of each quantity over the window, relative to the largest magnitude it has had. The
/// run stops as steady after two consecutive full windows, ending at or after `min_time`, with
/// every relative change below `tolerance`; when `divergence` is given, a relative change above
/// it after `min_time` stops the run as diverged. Otherwise it runs the step's full length.
///
/// A turbulent flow never stops fluctuating, so its fields never settle within a tolerance. With
/// `averaging`, each quantity's running mean from `averaging.start` on is what must settle
/// instead: its spread over a window, relative to the largest mean it has had, below `tolerance`
/// for two consecutive windows that both begin at or after the start. The report then gives each
/// quantity's mean, its fluctuation (the root mean square about the mean) and its extremes over
/// the averaging period.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Steady {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub monitors: Vec<FlowMonitor>,
    #[serde(default = "default_steady_tolerance")]
    pub tolerance: f64,
    /// Flow time the change is measured over (default: a twentieth of the step's length).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_time: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub divergence: Option<f64>,
    /// Stop on running means instead of the instantaneous values (a turbulent flow).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub averaging: Option<Averaging>,
}

/// How a statistically steady run averages ([`Steady::averaging`]).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Averaging {
    /// Flow time discarded as start-up before the running means begin.
    pub start: f64,
}

impl Steady {
    pub fn validate(&self) -> Result<()> {
        if !(self.tolerance > 0.0 && self.tolerance < 1.0) {
            return Err(Error::Invalid(format!("the steady-state tolerance is a relative change above 0 and below 1, got {}", self.tolerance)));
        }
        if let Some(w) = self.window {
            positive("steady-state window", w)?;
        }
        if let Some(t) = self.min_time
            && !(t >= 0.0 && t.is_finite())
        {
            return Err(Error::Invalid(format!("the steady-state minimum time is not negative, got {t}")));
        }
        if let Some(d) = self.divergence
            && !(d > self.tolerance && d.is_finite())
        {
            return Err(Error::Invalid(format!("the divergence limit is a relative change above the tolerance, got {d}")));
        }
        if let Some(a) = self.averaging
            && !(a.start >= 0.0 && a.start.is_finite())
        {
            return Err(Error::Invalid(format!("the averaging start is a flow time, not negative, got {}", a.start)));
        }
        if self.monitors.len() > MAX_FLOW_MONITORS {
            return Err(Error::Invalid(format!("a steady-state rule monitors at most {MAX_FLOW_MONITORS} quantities, got {}", self.monitors.len())));
        }
        for m in &self.monitors {
            if let FlowMonitor::Force { faces } | FlowMonitor::Heat { faces } = m
                && !(1..=MAX_MONITOR_FACES).contains(&faces.len())
            {
                return Err(Error::Invalid(format!("a monitored force or heat sums 1 to {MAX_MONITOR_FACES} faces, got {}", faces.len())));
            }
        }
        Ok(())
    }
}

/// A quantity a steady-state rule monitors besides the residual.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum FlowMonitor {
    /// The mass flow out through a box face.
    MassFlow { face: crate::BoxFace },
    /// The total force of the fluid on walls (box faces or obstacle faces).
    Force { faces: Vec<crate::FluidBoundary> },
    /// The total heat the fluid gives walls.
    Heat { faces: Vec<crate::FluidBoundary> },
    /// The mean temperature of the fluid.
    MeanTemperature {},
}

impl FlowMonitor {
    /// A short label for reports: `mass_flow x_max`, `force x_min, Body1/F3`, `mean_temperature`.
    pub fn label(&self) -> String {
        let faces = |faces: &[crate::FluidBoundary]| {
            faces
                .iter()
                .map(|f| match f {
                    crate::FluidBoundary::Box(b) => b.name().to_string(),
                    crate::FluidBoundary::Surface(r) => format!("{}/{}", r.body, r.entity.name().unwrap_or_default()),
                })
                .collect::<Vec<_>>()
                .join(", ")
        };
        match self {
            FlowMonitor::MassFlow { face } => format!("mass_flow {}", face.name()),
            FlowMonitor::Force { faces: f } => format!("force {}", faces(f)),
            FlowMonitor::Heat { faces: f } => format!("heat {}", faces(f)),
            FlowMonitor::MeanTemperature {} => "mean_temperature".into(),
        }
    }
}

/// Most load increments a nonlinear static step may take.
pub const MAX_NONLINEAR_INCREMENTS: u32 = 100_000;

/// Most Newton iterations a nonlinear static step may take in one increment.
pub const MAX_NONLINEAR_ITERATIONS: u32 = 1000;

/// The settings of a nonlinear static step.
///
/// The step's load factor λ runs from 0 to 1: every load is λ times its value and every
/// prescribed displacement λ times its value (a linear ramp). Increments are fractions of the
/// step: the first is `initial_increment`; an increment that converges quickly lets the next
/// grow (up to `max_increment`), one that diverges, fails to converge in `max_iterations` or
/// inverts an element is retried with a quarter of its size, and the step fails, naming the
/// load factor and the reason, when an increment would fall below `min_increment` or the step
/// has taken `max_increments` increments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct NonlinearStatic {
    /// Large displacements, rotations and strains (total Lagrangian elements, follower
    /// pressure); `false` keeps the kinematics of small displacements (only the materials are
    /// nonlinear).
    pub nonlinear_geometry: bool,
    /// The first increment, as a fraction of the step.
    pub initial_increment: f64,
    /// The smallest increment a cutback may reach.
    pub min_increment: f64,
    /// The largest increment.
    pub max_increment: f64,
    /// Cap on the increments (cutbacks included).
    pub max_increments: u32,
    /// Cap on the Newton iterations of one increment.
    pub max_iterations: u32,
    /// Convergence tolerances.
    pub tolerances: Tolerances,
    /// A line search along each Newton correction (off by default: full Newton steps converge
    /// best for elastic large-rotation problems; a search helps where the material softens).
    pub line_search: bool,
    /// The load factor against the step fraction: piecewise linear through these
    /// `[fraction, factor]` points, from `[0, 0]` to fraction 1 (`None`: the linear ramp, the
    /// factor equal to the fraction). Loads and prescribed displacements are the factor times
    /// their values, so a path that comes down unloads and one that turns negative reverses
    /// them (cyclic loading). Increments end at every point of the path.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub load_path: Option<Vec<[f64; 2]>>,
    /// A geometric imperfection: the mode shape of a buckling step, scaled, added to the node
    /// coordinates before this step is solved (the step and its results are on that geometry).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub imperfection: Option<Imperfection>,
    /// Arc-length (Riks–Crisfield) control: the load factor becomes an unknown, so the step can
    /// follow its equilibrium path through limit points (snap-through, snap-back). The first
    /// increment's load factor change is `initial_increment`; the increment sizes become arc
    /// lengths in proportion. `None`: load control.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub arc_length: Option<ArcLength>,
}

/// The constraint of an arc-length step.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArcConstraint {
    /// ‖Δu‖ = Δl (Crisfield's cylindrical arc).
    #[default]
    Cylindrical,
    /// ‖Δu‖² + Δλ² ‖u_f‖² = Δl², u_f the first increment's tangent displacement of the unit
    /// load (the spherical arc, scaled to displacement units).
    Spherical,
}

/// When an arc-length step ends: at `max_load_factor` (default 1), when the displacement of
/// `node` along `direction` (1, 2, 3) reaches `displacement` in magnitude, or at the step's
/// increment cap, whichever comes first.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ArcLength {
    pub constraint: ArcConstraint,
    pub max_load_factor: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub displacement_limit: Option<DisplacementLimit>,
}

impl Default for ArcLength {
    fn default() -> ArcLength {
        ArcLength { constraint: ArcConstraint::Cylindrical, max_load_factor: 1.0, displacement_limit: None }
    }
}

/// A displacement that ends an arc-length step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DisplacementLimit {
    pub node: u32,
    pub direction: u8,
    pub displacement: f64,
}

/// Mode `mode` (from 1) of buckling step `step`, scaled so its largest component is `scale`
/// (a length), as an initial imperfection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Imperfection {
    pub step: String,
    pub mode: u32,
    pub scale: f64,
}

/// Most points of a nonlinear step's load path.
pub const MAX_LOAD_PATH_POINTS: usize = 1000;

/// Load fractions this close to each other (relatively, or absolutely near 0) are the same.
const SAME_FRACTION: f64 = 1e-12; // unit-factors: a load fraction

/// When an increment's Newton iterations have converged: all three hold.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Tolerances {
    /// The largest residual force over the largest force (applied, internal or reaction) the
    /// step has seen, and the largest residual moment over the moment scale (the largest moment,
    /// at least that force times the model's size).
    pub force: f64,
    /// The largest translation of the correction over the largest translation of the increment,
    /// and its largest rotation over the largest rotation of the increment (at least that
    /// translation over the model's size).
    pub displacement: f64,
    /// The work of the last correction against the residual over that of the first.
    pub energy: f64,
}

impl Default for Tolerances {
    fn default() -> Tolerances {
        Tolerances { force: 1e-6, displacement: 1e-5, energy: 1e-10 } // unit-factors: relative tolerances
    }
}

impl Default for NonlinearStatic {
    fn default() -> NonlinearStatic {
        NonlinearStatic {
            nonlinear_geometry: true,
            initial_increment: 0.1,
            min_increment: 1e-5,
            max_increment: 1.0,
            max_increments: 100,
            max_iterations: 16,
            tolerances: Tolerances::default(),
            line_search: false,
            load_path: None,
            imperfection: None,
            arc_length: None,
        }
    }
}

impl NonlinearStatic {
    pub fn validate(&self) -> Result<()> {
        let fraction = |what: &str, v: f64| {
            if v > 0.0 && v <= 1.0 { Ok(v) } else { Err(Error::Invalid(format!("the {what} is a fraction of the step in (0, 1], got {v}"))) }
        };
        if let Some(arc) = &self.arc_length {
            if !(arc.max_load_factor.is_finite() && arc.max_load_factor > 0.0) {
                return Err(Error::Invalid(format!("an arc-length step's maximum load factor is a positive number, got {}", arc.max_load_factor)));
            }
            if let Some(l) = &arc.displacement_limit {
                if !(1..=3).contains(&l.direction) {
                    return Err(Error::Invalid(format!("a displacement limit's direction is 1, 2 or 3, got {}", l.direction)));
                }
                if !(l.displacement.is_finite() && l.displacement != 0.0) {
                    return Err(Error::Invalid(format!("a displacement limit is a nonzero number, got {}", l.displacement)));
                }
            }
            if self.load_path.is_some() {
                return Err(Error::Invalid("an arc-length step finds its own load factor; it takes no load path".into()));
            }
        }
        if let Some(imp) = &self.imperfection {
            if !(1..=MAX_MODES).contains(&imp.mode) {
                return Err(Error::Invalid(format!("an imperfection takes mode 1 to {MAX_MODES}, got {}", imp.mode)));
            }
            finite("imperfection scale", imp.scale)?;
        }
        let initial = fraction("initial increment", self.initial_increment)?;
        let min = fraction("minimum increment", self.min_increment)?;
        let max = fraction("maximum increment", self.max_increment)?;
        if !(min <= initial && initial <= max) {
            return Err(Error::Invalid(format!("the increments must satisfy minimum ≤ initial ≤ maximum, got {min}, {initial} and {max}")));
        }
        if !(1..=MAX_NONLINEAR_INCREMENTS).contains(&self.max_increments) {
            return Err(Error::Invalid(format!("a nonlinear step takes 1 to {MAX_NONLINEAR_INCREMENTS} increments, got {}", self.max_increments)));
        }
        if !(1..=MAX_NONLINEAR_ITERATIONS).contains(&self.max_iterations) {
            return Err(Error::Invalid(format!(
                "a nonlinear step takes 1 to {MAX_NONLINEAR_ITERATIONS} iterations per increment, got {}",
                self.max_iterations
            )));
        }
        let t = self.tolerances;
        for (what, v) in [("force", t.force), ("displacement", t.displacement), ("energy", t.energy)] {
            if !(v > 0.0 && v < 1.0) {
                return Err(Error::Invalid(format!("the {what} tolerance is in (0, 1), got {v}")));
            }
        }
        if let Some(path) = &self.load_path {
            if !(2..=MAX_LOAD_PATH_POINTS).contains(&path.len()) {
                return Err(Error::Invalid(format!("a load path has 2 to {MAX_LOAD_PATH_POINTS} points, got {}", path.len())));
            }
            if path.first() != Some(&[0.0, 0.0]) || path.last().map(|p| p[0]) != Some(1.0) {
                return Err(Error::Invalid("a load path starts at [0, 0] (the unloaded start of the step) and ends at fraction 1".into()));
            }
            if path.iter().flatten().any(|v| !v.is_finite()) {
                return Err(Error::Invalid("the points of a load path must be finite".into()));
            }
            if path.windows(2).any(|w| w.get(1).map(|b| b[0]) <= w.first().map(|a| a[0])) {
                return Err(Error::Invalid("the fractions of a load path must increase".into()));
            }
        }
        Ok(())
    }

    /// The load factor at step fraction `t` (in [0, 1]).
    pub fn load_factor(&self, t: f64) -> f64 {
        let Some(path) = &self.load_path else { return t };
        for w in path.windows(2) {
            if let [[t0, f0], [t1, f1]] = *w
                && t <= t1
            {
                return if t1 > t0 { f0 + (f1 - f0) * ((t - t0) / (t1 - t0)).clamp(0.0, 1.0) } else { f1 };
            }
        }
        path.last().map_or(t, |p| p[1])
    }

    /// The first point of the load path after fraction `t` (1 without a path), where an
    /// increment from `t` must end.
    pub fn next_turn(&self, t: f64) -> f64 {
        self.load_path.as_ref().and_then(|path| path.iter().map(|p| p[0]).find(|f| *f > t * (1.0 + SAME_FRACTION) + SAME_FRACTION)).unwrap_or(1.0)
    }
}

/// Most time increments of a transient heat transfer step.
pub const MAX_HEAT_INCREMENTS: f64 = 100_000.0;

/// An increment count this close above a whole number is that number: what is left of dividing
/// the end time by the increment.
const INCREMENT_ROUND_OFF: f64 = 1e-9; // unit-factors: a fraction of an increment

/// Time integration of a transient heat transfer step: the generalized trapezoidal (θ) method
/// with fixed increments, from a uniform initial temperature.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Transient {
    pub end_time: f64,
    pub increment: f64,
    /// 1 is backward Euler (the default), ½ Crank–Nicolson; θ ≥ ½ is unconditionally stable.
    #[serde(default = "default_theta")]
    pub theta: f64,
    /// The uniform temperature at the start (kelvin; 0 °C when left out).
    #[serde(default = "default_initial_temperature")]
    pub initial_temperature: f64,
}

fn default_theta() -> f64 {
    1.0
}

fn default_initial_temperature() -> f64 {
    crate::units::CELSIUS_ZERO
}

impl Transient {
    /// Backward Euler from 0.
    pub fn new(end_time: f64, increment: f64) -> Transient {
        Transient { end_time, increment, theta: default_theta(), initial_temperature: default_initial_temperature() }
    }

    /// The number of increments (the last one may be shorter).
    pub fn increments(&self) -> usize {
        (self.end_time / self.increment - INCREMENT_ROUND_OFF).ceil().max(1.0) as usize
    }

    pub fn validate(&self) -> Result<()> {
        let end = positive("end time", self.end_time)?;
        let dt = positive("time increment", self.increment)?;
        if end / dt > MAX_HEAT_INCREMENTS {
            return Err(Error::Invalid(format!("a transient heat step takes at most {MAX_HEAT_INCREMENTS} increments, got {}", end / dt)));
        }
        if !(0.5..=1.0).contains(&self.theta) {
            return Err(Error::Invalid(format!("θ is in [0.5, 1] (unconditionally stable), got {}", self.theta)));
        }
        finite("initial temperature", self.initial_temperature)?;
        Ok(())
    }
}

/// The settings of an explicit dynamic step.
///
/// Loads act in full from the start of the step; a displacement boundary condition with a value
/// is reached at constant velocity over the step (value / `end_time`), a zero value holds the
/// node. The step starts from the undeformed model at rest, apart from `initial_velocities`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExplicitDynamic {
    /// Duration of the step (time units of the unit system).
    pub end_time: f64,
    /// Fraction of the stable time increment taken (0 < scale ≤ 1).
    #[serde(default = "default_dt_scale")]
    pub dt_scale: f64,
    /// Cap on the time increments; the run stops early (and says so) when it is reached.
    #[serde(default = "default_max_increments")]
    pub max_increments: u64,
    /// Time between field frames (displacement, velocity, stress); by default a twentieth of the
    /// step.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub output_interval: Option<f64>,
    /// Time between history samples (energies, momentum, constraint forces and the histories of
    /// `history_nodes`); by default a thousandth of the step.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub history_interval: Option<f64>,
    /// A set whose nodes' displacement and velocity are recorded as histories.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub history_nodes: Option<String>,
    /// Velocities of the nodes of sets at the start of the step (a node in several sets takes
    /// the last one listed).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub initial_velocities: Vec<InitialVelocity>,
    /// Fixed rigid planes the nodes of sets cannot pass through, with Coulomb friction.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rigid_walls: Vec<RigidWall>,
    /// Hourglass control of one-point hexahedra.
    #[serde(default)]
    pub hourglass: Hourglass,
    /// Artificial bulk viscosity.
    #[serde(default)]
    pub bulk_viscosity: BulkViscosity,
    /// Large deformation (updated Lagrangian: gradients of the current configuration, stress
    /// updated in rate form with an objective stress rate); `false` keeps the gradients of the
    /// undeformed mesh (small displacements and rotations).
    #[serde(default = "yes")]
    pub nonlinear_geometry: bool,
    /// The objective stress rate of large deformation.
    #[serde(default)]
    pub stress_rate: StressRate,
    /// Opt-in mass scaling: the elements whose stable increment is below a target get the
    /// density that brings them to it (the added mass is reported).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mass_scaling: Option<MassScaling>,
}

/// Mass scaling to a target time increment: every element (of `set`, or of the whole model) whose
/// stable increment, times the step's scale, is below `target_increment` has its density raised
/// by (target / increment)², which brings it to the target. It changes the dynamics: the added
/// mass is reported, per section, with its share of the kinetic energy.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MassScaling {
    pub target_increment: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub set: Option<String>,
}

fn yes() -> bool {
    true
}

/// The objective stress rate a large-deformation stress update uses.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StressRate {
    /// Jaumann (co-rotational with the spin), integrated with the Hughes–Winget rotation.
    #[default]
    Jaumann,
    /// Green–Naghdi (co-rotational with the rotation of the polar decomposition), integrated as
    /// Flanagan and Taylor do.
    GreenNaghdi,
}

fn default_dt_scale() -> f64 {
    0.9
}

fn default_max_increments() -> u64 {
    10_000_000
}

impl ExplicitDynamic {
    /// A step of `end_time` with every other setting at its default.
    pub fn new(end_time: f64) -> ExplicitDynamic {
        ExplicitDynamic {
            end_time,
            dt_scale: default_dt_scale(),
            max_increments: default_max_increments(),
            output_interval: None,
            history_interval: None,
            history_nodes: None,
            initial_velocities: Vec::new(),
            rigid_walls: Vec::new(),
            hourglass: Hourglass::default(),
            bulk_viscosity: BulkViscosity::default(),
            nonlinear_geometry: true,
            stress_rate: StressRate::default(),
            mass_scaling: None,
        }
    }

    /// The time between field frames.
    pub fn output_interval(&self) -> f64 {
        self.output_interval.unwrap_or(self.end_time / 20.0)
    }

    /// The time between history samples.
    pub fn history_interval(&self) -> f64 {
        self.history_interval.unwrap_or(self.end_time / DEFAULT_HISTORY_SAMPLES)
    }

    /// The sets the step refers to.
    pub fn sets(&self) -> impl Iterator<Item = &str> {
        self.history_nodes
            .iter()
            .map(String::as_str)
            .chain(self.initial_velocities.iter().map(|v| v.set.as_str()))
            .chain(self.rigid_walls.iter().map(|w| w.set.as_str()))
            .chain(self.mass_scaling.iter().filter_map(|m| m.set.as_deref()))
    }

    pub fn validate(&self) -> Result<()> {
        let end = positive("end time", self.end_time)?;
        if !(self.dt_scale > 0.0 && self.dt_scale <= 1.0) {
            return Err(Error::Invalid(format!("the time increment scale is in (0, 1], got {}", self.dt_scale)));
        }
        if !(1..=MAX_INCREMENTS).contains(&self.max_increments) {
            return Err(Error::Invalid(format!("an explicit step takes 1 to {MAX_INCREMENTS} increments, got {}", self.max_increments)));
        }
        for (what, interval, most) in
            [("output interval", self.output_interval, MAX_FRAMES), ("history interval", self.history_interval, MAX_HISTORY_SAMPLES)]
        {
            if let Some(i) = interval {
                positive(what, i)?;
                if end / i > most {
                    return Err(Error::Invalid(format!("the {what} {i} gives more than {most} samples over the step")));
                }
            }
        }
        if self.initial_velocities.len() > MAX_STEP_ITEMS || self.rigid_walls.len() > MAX_STEP_ITEMS {
            return Err(Error::Limit(format!("a step lists at most {MAX_STEP_ITEMS} initial velocities and rigid walls")));
        }
        for v in &self.initial_velocities {
            for c in v.velocity {
                finite("initial velocity", c)?;
            }
        }
        for w in &self.rigid_walls {
            for c in w.point.iter().chain(&w.normal) {
                finite("rigid wall", *c)?;
            }
            if !(w.normal.iter().map(|c| c * c).sum::<f64>() > 0.0) {
                return Err(Error::Invalid("a rigid wall's normal must not be zero".into()));
            }
            w.friction.validate().map_err(|e| Error::Invalid(format!("rigid wall on set `{}`: {e}", w.set)))?;
        }
        if let Some(m) = &self.mass_scaling {
            positive("mass scaling target increment", m.target_increment)?;
            if m.target_increment > end {
                return Err(Error::Invalid(format!("the mass scaling target increment {} is longer than the step", m.target_increment)));
            }
        }
        self.hourglass.validate()?;
        self.bulk_viscosity.validate()
    }
}

/// A velocity given to the nodes of a set at the start of an explicit step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InitialVelocity {
    pub set: String,
    pub velocity: [f64; 3],
}

/// A fixed rigid plane through `point`. Its `normal` points to the side the nodes of `set` stay
/// on; nodes sliding on it feel `friction` (frictionless by default).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RigidWall {
    pub set: String,
    pub point: [f64; 3],
    pub normal: [f64; 3],
    #[serde(default, skip_serializing_if = "Friction::is_frictionless")]
    pub friction: Friction,
}

/// How one-point hexahedra resist their hourglass (zero-energy) modes.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HourglassKind {
    /// Forces proportional to the hourglass velocities.
    #[default]
    Viscous,
    /// Forces proportional to the hourglass displacements.
    Stiffness,
}

/// Hourglass control: the kind and a dimensionless coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hourglass {
    #[serde(rename = "type", default)]
    pub kind: HourglassKind,
    #[serde(default = "default_hourglass_coefficient")]
    pub coefficient: f64,
}

/// Largest hourglass coefficient.
pub const MAX_HOURGLASS_COEFFICIENT: f64 = 0.25;

/// History samples of an explicit dynamic step when no interval is given.
const DEFAULT_HISTORY_SAMPLES: f64 = 1000.0; // unit-factors: a count of samples

fn default_hourglass_coefficient() -> f64 {
    0.1
}

impl Default for Hourglass {
    fn default() -> Hourglass {
        Hourglass { kind: HourglassKind::Viscous, coefficient: default_hourglass_coefficient() }
    }
}

impl Hourglass {
    pub fn validate(&self) -> Result<()> {
        if !(0.0..=MAX_HOURGLASS_COEFFICIENT).contains(&self.coefficient) {
            return Err(Error::Invalid(format!("the hourglass coefficient is in 0..={MAX_HOURGLASS_COEFFICIENT}, got {}", self.coefficient)));
        }
        Ok(())
    }
}

/// Artificial bulk viscosity in compression: a pressure `ρ L (quadratic L ε̇ᵥ² − linear c ε̇ᵥ)`
/// for a volumetric strain rate ε̇ᵥ < 0, element size L and wave speed c.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BulkViscosity {
    #[serde(default = "default_linear_viscosity")]
    pub linear: f64,
    #[serde(default = "default_quadratic_viscosity")]
    pub quadratic: f64,
}

/// Largest bulk viscosity coefficient.
pub const MAX_BULK_VISCOSITY: f64 = 10.0;

fn default_linear_viscosity() -> f64 {
    0.06
}

fn default_quadratic_viscosity() -> f64 {
    1.5
}

impl Default for BulkViscosity {
    fn default() -> BulkViscosity {
        BulkViscosity { linear: default_linear_viscosity(), quadratic: default_quadratic_viscosity() }
    }
}

impl BulkViscosity {
    pub fn validate(&self) -> Result<()> {
        for (what, v) in [("linear", self.linear), ("quadratic", self.quadratic)] {
            if !(0.0..=MAX_BULK_VISCOSITY).contains(&v) {
                return Err(Error::Invalid(format!("the {what} bulk viscosity coefficient is in 0..={MAX_BULK_VISCOSITY}, got {v}")));
            }
        }
        Ok(())
    }
}

impl Procedure {
    /// The procedure's type, as the model file names it.
    pub fn name(&self) -> &'static str {
        match self {
            Procedure::Static {} => "static",
            Procedure::Frequency { .. } => "frequency",
            Procedure::Buckle { .. } => "buckle",
            Procedure::ExplicitDynamic(_) => "explicit_dynamic",
            Procedure::HeatTransfer { .. } => "heat_transfer",
            Procedure::CoupledTemperatureDisplacement {} => "coupled_temperature_displacement",
            Procedure::Fluid { .. } => "fluid",
            Procedure::NonlinearStatic(_) => "nonlinear_static",
        }
    }

    /// The fluid domain a fluid step solves in.
    pub fn fluid_domain(&self) -> Option<&str> {
        match self {
            Procedure::Fluid { domain, .. } => Some(domain),
            _ => None,
        }
    }

    /// Whether the procedure solves for the motion of the solids (not pure heat transfer or
    /// fluid flow).
    pub fn is_structural(&self) -> bool {
        !matches!(self, Procedure::HeatTransfer { .. } | Procedure::Fluid { .. })
    }

    /// Whether the procedure solves for temperatures.
    pub fn is_thermal(&self) -> bool {
        matches!(self, Procedure::HeatTransfer { .. } | Procedure::CoupledTemperatureDisplacement {})
    }

    /// Whether the procedure is transient dynamics (it needs mass; a free body is legitimate).
    pub fn is_dynamic(&self) -> bool {
        matches!(self, Procedure::ExplicitDynamic(_))
    }

    /// The sets the procedure's own settings refer to.
    pub fn sets(&self) -> impl Iterator<Item = &str> {
        let explicit = match self {
            Procedure::ExplicitDynamic(e) => Some(e),
            _ => None,
        };
        explicit.into_iter().flat_map(ExplicitDynamic::sets)
    }

    pub fn validate(&self) -> Result<()> {
        if let Procedure::ExplicitDynamic(e) = self {
            return e.validate();
        }
        if let Procedure::NonlinearStatic(n) = self {
            return n.validate();
        }
        if let Procedure::Fluid { end_time, steps, time_step, frames, magic, ref steady, .. } = *self {
            if let Some(st) = steady {
                st.validate()?;
            }
            match (end_time, steps) {
                (Some(t), None) => {
                    crate::positive("end time", t)?;
                }
                (None, Some(n)) if (1..=MAX_FLUID_STEPS).contains(&n) => {}
                (None, Some(n)) => return Err(Error::Invalid(format!("a fluid step runs 1 to {MAX_FLUID_STEPS} time steps, got {n}"))),
                _ => return Err(Error::Invalid("a fluid step needs exactly one of `end_time` and `steps`".into())),
            }
            if let Some(dt) = time_step {
                crate::positive("time step", dt)?;
            }
            if !(1..=MAX_FLUID_FRAMES).contains(&frames) {
                return Err(Error::Invalid(format!("a fluid step writes 1 to {MAX_FLUID_FRAMES} frames, got {frames}")));
            }
            if let Some(m) = magic
                && !(m > 0.0 && m <= 1.0)
            {
                return Err(Error::Invalid(format!("the two-relaxation-time parameter Λ must be above 0 and at most 1, got {m}")));
            }
        }
        if let Procedure::HeatTransfer { transient: Some(t) } = self {
            return t.validate();
        }
        if let Procedure::Buckle { modes, shift } = *self {
            if !(1..=MAX_MODES).contains(&modes) {
                return Err(Error::Invalid(format!("a buckling step extracts 1 to {MAX_MODES} modes, got {modes}")));
            }
            if let Some(s) = shift {
                finite("buckling shift", s)?;
            }
        }
        if let Procedure::Frequency { modes, min_hz, max_hz } = *self {
            if !(1..=MAX_MODES).contains(&modes) {
                return Err(Error::Invalid(format!("a frequency step extracts 1 to {MAX_MODES} modes, got {modes}")));
            }
            for f in [min_hz, max_hz].into_iter().flatten() {
                if finite("frequency", f)? < 0.0 {
                    return Err(Error::Invalid(format!("frequencies are not negative, got {f}")));
                }
            }
            if let (Some(lo), Some(hi)) = (min_hz, max_hz)
                && lo >= hi
            {
                return Err(Error::Invalid(format!("the frequency range {lo}..{hi} is empty")));
            }
        }
        Ok(())
    }

    /// The output a step of this procedure requests unless told otherwise.
    pub fn default_output(&self) -> BTreeSet<OutputVariable> {
        match self {
            Procedure::Static {} | Procedure::NonlinearStatic(_) => OutputVariable::ALL.into_iter().collect(),
            Procedure::Frequency { .. } | Procedure::Buckle { .. } => [OutputVariable::Displacement].into_iter().collect(),
            Procedure::ExplicitDynamic(_) => [OutputVariable::Displacement, OutputVariable::Stress].into_iter().collect(),
            // Velocity, pressure and wall forces, always.
            Procedure::Fluid { .. } => BTreeSet::new(),
            Procedure::HeatTransfer { .. } => [OutputVariable::Temperature, OutputVariable::HeatFlux].into_iter().collect(),
            Procedure::CoupledTemperatureDisplacement {} => {
                OutputVariable::ALL.into_iter().chain([OutputVariable::Temperature, OutputVariable::HeatFlux]).collect()
            }
        }
    }
}

/// A result a step writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OutputVariable {
    /// Nodal displacements (mode shapes in a frequency step).
    Displacement,
    /// Reaction forces at constrained nodes.
    ReactionForce,
    Stress,
    Strain,
    /// Nodal temperatures (heat transfer).
    Temperature,
    /// Heat flux at the integration points (heat transfer).
    HeatFlux,
}

impl OutputVariable {
    /// The structural variables.
    pub const ALL: [OutputVariable; 4] =
        [OutputVariable::Displacement, OutputVariable::ReactionForce, OutputVariable::Stress, OutputVariable::Strain];
}

/// An analysis step.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Step {
    pub name: String,
    pub procedure: Procedure,
    pub output: BTreeSet<OutputVariable>,
    /// Boundary conditions active in this step.
    #[serde(default)]
    pub boundary_conditions: BTreeSet<String>,
    /// Loads active in this step.
    #[serde(default)]
    pub loads: BTreeSet<String>,
}

impl Step {
    /// A step with the procedure's default output and nothing active.
    pub fn new(name: impl Into<String>, procedure: Procedure) -> Step {
        Step { name: name.into(), output: procedure.default_output(), procedure, boundary_conditions: BTreeSet::new(), loads: BTreeSet::new() }
    }

    /// Change the procedure. Moving between a structural and a fluid procedure starts from the
    /// new procedure's default output and, for a fluid step, nothing active (a fluid step takes
    /// its conditions from its domain); between structural procedures the output and table stay.
    pub fn set_procedure(&mut self, procedure: Procedure) {
        if procedure.fluid_domain().is_some() != self.procedure.fluid_domain().is_some() {
            self.output = procedure.default_output();
            if procedure.fluid_domain().is_some() {
                self.boundary_conditions.clear();
                self.loads.clear();
            }
        }
        self.procedure = procedure;
    }
}

/// What solving one step reads from a [`Model`] ([`Model::step_inputs`]), borrowed: two models
/// whose step inputs serialize the same give the step the same result.
#[derive(Debug, Serialize)]
pub struct StepInputs<'a> {
    nodes: &'a BTreeMap<u32, crate::Node>,
    elements: &'a BTreeMap<u32, crate::Element>,
    classification: &'a crate::Classification,
    sets: BTreeMap<&'a str, &'a crate::Set>,
    materials: &'a BTreeMap<String, crate::Material>,
    sections: &'a BTreeMap<String, crate::Section>,
    reference_points: &'a BTreeMap<String, crate::ReferencePoint>,
    couplings: &'a BTreeMap<String, crate::Coupling>,
    contacts: &'a BTreeMap<String, crate::Contact>,
    ties: &'a BTreeMap<String, crate::Tie>,
    masses: &'a BTreeMap<String, crate::PointMass>,
    springs: &'a BTreeMap<String, crate::Spring>,
    boundary_conditions: BTreeMap<&'a str, &'a crate::BoundaryCondition>,
    loads: BTreeMap<&'a str, &'a crate::Load>,
    fluid_domains: &'a BTreeMap<String, crate::FluidDomain>,
    step: &'a Step,
}

impl Model {
    /// What solving step `name` reads: the mesh, the materials, sections, interactions and every
    /// definition that acts in all steps, the step itself and the boundary conditions, loads and
    /// sets it uses. Other steps and what only they use are left out (steps are solved each on
    /// its own), and so are the meshing inputs (datums, partitions, seeds: the mesh they made is
    /// in). `None` when there is no such step.
    pub fn step_inputs(&self, name: &str) -> Option<StepInputs<'_>> {
        let step = self.step(name)?;
        Some(StepInputs {
            nodes: &self.nodes,
            elements: &self.elements,
            classification: &self.classification,
            sets: self.sets.iter().filter(|(n, _)| self.step_uses_set(step, n)).map(|(n, s)| (n.as_str(), s)).collect(),
            materials: &self.materials,
            sections: &self.sections,
            reference_points: &self.reference_points,
            couplings: &self.couplings,
            contacts: &self.contacts,
            ties: &self.ties,
            masses: &self.masses,
            springs: &self.springs,
            boundary_conditions: self
                .boundary_conditions
                .iter()
                .filter(|(n, _)| step.boundary_conditions.contains(*n))
                .map(|(n, b)| (n.as_str(), b))
                .collect(),
            loads: self.loads.iter().filter(|(n, _)| step.loads.contains(*n)).map(|(n, l)| (n.as_str(), l)).collect(),
            fluid_domains: &self.fluid_domains,
            step,
        })
    }

    /// The steps in analysis order.
    pub fn steps(&self) -> &[Step] {
        &self.steps
    }

    pub fn step(&self, name: &str) -> Option<&Step> {
        self.steps.iter().find(|s| s.name == name)
    }

    /// The step named `name`, or the only step whose procedure is `wanted` when no name is
    /// given; it must be such a step (`what` names the kind in messages, e.g. "static").
    pub fn select_step(&self, name: Option<&str>, wanted: fn(&Procedure) -> bool, what: &str) -> Result<&Step> {
        let step = match name {
            Some(n) => self.step(n).ok_or_else(|| Error::Invalid(format!("step `{n}` does not exist")))?,
            None => {
                let matching: Vec<&Step> = self.steps.iter().filter(|s| wanted(&s.procedure)).collect();
                match matching.as_slice() {
                    [one] => *one,
                    [] => return Err(Error::Invalid(format!("the model has no {what} step"))),
                    many => {
                        let names: Vec<&str> = many.iter().map(|s| s.name.as_str()).collect();
                        return Err(Error::Invalid(format!("name the step: the model has {what} steps {}", names.join(", "))));
                    }
                }
            }
        };
        if !wanted(&step.procedure) {
            return Err(Error::Invalid(format!("step `{}` is a {} step, not a {what} one", step.name, step.procedure.name())));
        }
        Ok(step)
    }

    fn step_index(&self, name: &str) -> Result<usize> {
        self.steps.iter().position(|s| s.name == name).ok_or_else(|| Error::NotFound(format!("step `{name}`")))
    }

    fn step_mut(&mut self, name: &str) -> Result<&mut Step> {
        self.steps.iter_mut().find(|s| s.name == name).ok_or_else(|| Error::NotFound(format!("step `{name}`")))
    }

    fn check_step(&self, step: &Step) -> Result<()> {
        step.procedure.validate().map_err(|e| context("step", &step.name, e))?;
        for set in step.procedure.sets() {
            require(&self.sets, "set", set).map_err(|e| context("step", &step.name, e))?;
        }
        if let Some(domain) = step.procedure.fluid_domain() {
            require(&self.fluid_domains, "fluid domain", domain).map_err(|e| context("step", &step.name, e))?;
            if !(step.boundary_conditions.is_empty() && step.loads.is_empty() && step.output.is_empty()) {
                return Err(context(
                    "step",
                    &step.name,
                    Error::Invalid("a fluid step takes its conditions from its fluid domain and always writes velocity, pressure and wall forces: list no boundary conditions, loads or output".into()),
                ));
            }
        }
        for b in &step.boundary_conditions {
            require(&self.boundary_conditions, "boundary condition", b).map_err(|e| context("step", &step.name, e))?;
        }
        for l in &step.loads {
            require(&self.loads, "load", l).map_err(|e| context("step", &step.name, e))?;
        }
        Ok(())
    }

    /// Add a step at the end, or just before the step named `before`.
    pub fn add_step(&mut self, step: Step, before: Option<&str>) -> Result<()> {
        check_name("step", &step.name)?;
        if self.step(&step.name).is_some() {
            return Err(Error::Invalid(format!("step `{}` already exists", step.name)));
        }
        if self.steps.len() >= MAX_STEPS {
            return Err(Error::Limit(format!("a model holds at most {MAX_STEPS} steps")));
        }
        self.check_step(&step)?;
        let at = match before {
            Some(b) => self.step_index(b)?,
            None => self.steps.len(),
        };
        self.steps.insert(at, step);
        Ok(())
    }

    /// Replace a step's definition (procedure, output and table entries), keeping its place.
    pub fn replace_step(&mut self, step: Step) -> Result<()> {
        self.check_step(&step)?;
        let name = step.name.clone();
        *self.step_mut(&name)? = step;
        Ok(())
    }

    pub fn remove_step(&mut self, name: &str) -> Result<()> {
        let i = self.step_index(name)?;
        self.steps.remove(i);
        Ok(())
    }

    /// Make a boundary condition active or inactive in a step.
    pub fn set_boundary_condition_active(&mut self, step: &str, name: &str, active: bool) -> Result<()> {
        require(&self.boundary_conditions, "boundary condition", name)?;
        switch(&mut self.step_mut(step)?.boundary_conditions, name, active);
        Ok(())
    }

    /// Make a load active or inactive in a step.
    pub fn set_load_active(&mut self, step: &str, name: &str, active: bool) -> Result<()> {
        require(&self.loads, "load", name)?;
        switch(&mut self.step_mut(step)?.loads, name, active);
        Ok(())
    }
}

fn switch(active: &mut BTreeSet<String>, name: &str, on: bool) {
    if on {
        active.insert(name.to_string());
    } else {
        active.remove(name);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BoundaryCondition, Set, Target};

    #[test]
    fn steps_keep_order_and_an_explicit_table() {
        let mut m = Model::default();
        m.add_set("base", Set { nodes: vec![1], ..Set::default() }).unwrap();
        m.add_boundary_condition("fix", BoundaryCondition::Fixed { target: Target::Set("base".into()) }).unwrap();
        m.add_step(Step::new("Load", Procedure::Static {}), None).unwrap();
        m.add_step(Step::new("Modes", Procedure::Frequency { modes: 10, min_hz: None, max_hz: None }), Some("Load")).unwrap();
        assert_eq!(m.steps().iter().map(|s| s.name.as_str()).collect::<Vec<_>>(), ["Modes", "Load"]);
        assert!(m.add_step(Step::new("Load", Procedure::Static {}), None).is_err());
        assert!(m.add_step(Step::new("X", Procedure::Static {}), Some("Nope")).is_err());
        assert!(m.add_step(Step::new("Y", Procedure::Frequency { modes: 0, min_hz: None, max_hz: None }), None).is_err());
        assert!(m.add_step(Step::new("Y", Procedure::Frequency { modes: 5, min_hz: Some(10.0), max_hz: Some(1.0) }), None).is_err());
        assert!(serde_json::from_value::<Procedure>(serde_json::json!({"type": "static", "modes": 3})).is_err(), "typos are errors");

        m.set_boundary_condition_active("Load", "fix", true).unwrap();
        assert!(m.step("Load").unwrap().boundary_conditions.contains("fix"));
        assert!(m.step("Modes").unwrap().boundary_conditions.is_empty(), "no propagation");
        assert!(m.set_boundary_condition_active("Load", "nope", true).is_err());
        assert!(m.set_load_active("Load", "fix", true).is_err(), "loads and boundary conditions are separate");

        // Nothing still referred to can be deleted: the set by the boundary condition, the
        // boundary condition by the step that lists it.
        assert!(matches!(m.remove_set("base"), Err(Error::InUse(_))));
        assert!(m.remove_boundary_condition("fix").unwrap_err().to_string().contains("step `Load`"));
        m.set_boundary_condition_active("Load", "fix", false).unwrap();
        m.remove_boundary_condition("fix").unwrap();
        m.remove_set("base").unwrap();
        m.remove_step("Modes").unwrap();
        assert_eq!(m.steps().len(), 1);
    }

    #[test]
    fn explicit_steps_validate_and_hold_their_sets() {
        let mut m = Model::default();
        m.add_set("front", Set { nodes: vec![1], ..Set::default() }).unwrap();
        let json = serde_json::json!({
            "type": "explicit_dynamic", "end_time": 1e-3, "history_nodes": "front",
            "initial_velocities": [{"set": "front", "velocity": [-10, 0, 0]}],
            "rigid_walls": [{"set": "front", "point": [0, 0, 0], "normal": [1, 0, 0], "friction": {"static": 0.2}}],
            "hourglass": {"type": "stiffness"}
        });
        let p: Procedure = serde_json::from_value(json).unwrap();
        let Procedure::ExplicitDynamic(e) = &p else { panic!("{p:?}") };
        assert_eq!((e.dt_scale, e.hourglass.kind, e.hourglass.coefficient), (0.9, HourglassKind::Stiffness, 0.1));
        assert_eq!(e.bulk_viscosity, BulkViscosity { linear: 0.06, quadratic: 1.5 });
        assert!((e.output_interval() - 5e-5).abs() < 1e-18 && (e.history_interval() - 1e-6).abs() < 1e-18);
        let back: Procedure = serde_json::from_value(serde_json::to_value(&p).unwrap()).unwrap();
        assert_eq!(back, p);
        m.add_step(Step::new("Impact", p.clone()), None).unwrap();
        assert!(matches!(m.remove_set("front"), Err(Error::InUse(_))), "the step refers to the set");
        assert!(m.add_step(Step::new("Impact", p), None).is_err(), "names stay unique");
        let missing = serde_json::json!({"type": "explicit_dynamic", "end_time": 1.0, "history_nodes": "nope"});
        assert!(m.add_step(Step::new("Bad", serde_json::from_value(missing).unwrap()), None).is_err());

        let base = ExplicitDynamic::new(1.0);
        let bad = [
            ExplicitDynamic { end_time: 0.0, ..base.clone() },
            ExplicitDynamic { end_time: f64::NAN, ..base.clone() },
            ExplicitDynamic { dt_scale: 1.5, ..base.clone() },
            ExplicitDynamic { dt_scale: 0.0, ..base.clone() },
            ExplicitDynamic { max_increments: 0, ..base.clone() },
            ExplicitDynamic { output_interval: Some(1e-9), ..base.clone() },
            ExplicitDynamic { history_interval: Some(-1.0), ..base.clone() },
            ExplicitDynamic {
                rigid_walls: vec![RigidWall { set: "front".into(), point: [0.0; 3], normal: [0.0; 3], friction: Friction::NONE }],
                ..base.clone()
            },
            ExplicitDynamic {
                rigid_walls: vec![RigidWall { set: "front".into(), point: [0.0; 3], normal: [1.0, 0.0, 0.0], friction: Friction::coulomb(f64::NAN) }],
                ..base.clone()
            },
            ExplicitDynamic {
                initial_velocities: vec![InitialVelocity { set: "front".into(), velocity: [f64::INFINITY, 0.0, 0.0] }],
                ..base.clone()
            },
            ExplicitDynamic { hourglass: Hourglass { kind: HourglassKind::Viscous, coefficient: 2.0 }, ..base.clone() },
            ExplicitDynamic { bulk_viscosity: BulkViscosity { linear: -1.0, quadratic: 1.0 }, ..base.clone() },
            ExplicitDynamic { mass_scaling: Some(MassScaling { target_increment: 0.0, set: None }), ..base.clone() },
            ExplicitDynamic { mass_scaling: Some(MassScaling { target_increment: f64::NAN, set: None }), ..base.clone() },
            ExplicitDynamic { mass_scaling: Some(MassScaling { target_increment: 2.0, set: None }), ..base.clone() },
        ];
        for b in bad {
            assert!(Procedure::ExplicitDynamic(b.clone()).validate().is_err(), "{b:?}");
        }
        assert!(serde_json::from_value::<Procedure>(serde_json::json!({"type": "explicit_dynamic", "end_time": 1, "dt": 2})).is_err());
    }

    #[test]
    fn nonlinear_static_steps_validate_their_controls() {
        let p: Procedure = serde_json::from_value(serde_json::json!({"type": "nonlinear_static"})).unwrap();
        assert_eq!(p, Procedure::NonlinearStatic(NonlinearStatic::default()));
        assert_eq!(p.name(), "nonlinear_static");
        assert!(p.validate().is_ok() && p.is_structural() && !p.is_dynamic());
        assert_eq!(p.default_output(), Procedure::Static {}.default_output());
        let p: Procedure = serde_json::from_value(serde_json::json!({
            "type": "nonlinear_static", "nonlinear_geometry": false, "initial_increment": 0.5, "min_increment": 0.01,
            "max_increments": 20, "max_iterations": 8, "tolerances": {"force": 1e-4}, "line_search": true
        }))
        .unwrap();
        let Procedure::NonlinearStatic(n) = &p else { panic!("{p:?}") };
        assert_eq!((n.nonlinear_geometry, n.max_increment, n.tolerances.force, n.tolerances.energy), (false, 1.0, 1e-4, 1e-10));
        let back: Procedure = serde_json::from_value(serde_json::to_value(&p).unwrap()).unwrap();
        assert_eq!(back, p);
        let base = NonlinearStatic::default();
        let bad = [
            NonlinearStatic { initial_increment: 0.0, ..base.clone() },
            NonlinearStatic { initial_increment: f64::NAN, ..base.clone() },
            NonlinearStatic { initial_increment: 2.0, max_increment: 2.0, ..base.clone() },
            NonlinearStatic { min_increment: 0.5, initial_increment: 0.1, ..base.clone() },
            NonlinearStatic { max_increment: 0.05, ..base.clone() },
            NonlinearStatic { max_increments: 0, ..base.clone() },
            NonlinearStatic { max_increments: MAX_NONLINEAR_INCREMENTS + 1, ..base.clone() },
            NonlinearStatic { max_iterations: 0, ..base.clone() },
            NonlinearStatic { max_iterations: u32::MAX, ..base.clone() },
            NonlinearStatic { tolerances: Tolerances { force: 0.0, ..Tolerances::default() }, ..base.clone() },
            NonlinearStatic { tolerances: Tolerances { energy: f64::INFINITY, ..Tolerances::default() }, ..base.clone() },
            NonlinearStatic { tolerances: Tolerances { displacement: 1.0, ..Tolerances::default() }, ..base.clone() },
        ];
        for b in bad {
            assert!(Procedure::NonlinearStatic(b.clone()).validate().is_err(), "{b:?}");
        }
        assert!(serde_json::from_value::<Procedure>(serde_json::json!({"type": "nonlinear_static", "nlgeom": true})).is_err());
        // A load path: up, down through zero, back to zero.
        let cyclic = NonlinearStatic { load_path: Some(vec![[0.0, 0.0], [0.25, 1.0], [0.75, -1.0], [1.0, 0.0]]), ..base.clone() };
        assert!(Procedure::NonlinearStatic(cyclic.clone()).validate().is_ok());
        for (t, f) in [(0.0, 0.0), (0.125, 0.5), (0.25, 1.0), (0.5, 0.0), (0.75, -1.0), (0.875, -0.5), (1.0, 0.0)] {
            assert!((cyclic.load_factor(t) - f).abs() < 1e-15, "{t}: {}", cyclic.load_factor(t));
        }
        assert_eq!((cyclic.next_turn(0.0), cyclic.next_turn(0.25), cyclic.next_turn(0.3), cyclic.next_turn(0.75)), (0.25, 0.75, 0.75, 1.0));
        assert_eq!((base.load_factor(0.3), base.next_turn(0.3)), (0.3, 1.0));
        for path in [
            vec![[0.0, 0.0]],
            vec![[0.0, 0.5], [1.0, 1.0]],
            vec![[0.0, 0.0], [0.9, 1.0]],
            vec![[0.0, 0.0], [0.5, 1.0], [0.5, 2.0], [1.0, 0.0]],
            vec![[0.0, 0.0], [0.5, f64::NAN], [1.0, 0.0]],
            vec![[0.0, 0.0]; MAX_LOAD_PATH_POINTS + 1],
        ] {
            assert!(Procedure::NonlinearStatic(NonlinearStatic { load_path: Some(path.clone()), ..base.clone() }).validate().is_err(), "{path:?}");
        }
    }

    #[test]
    fn heat_transfer_procedures() {
        let steady: Procedure = serde_json::from_value(serde_json::json!({"type": "heat_transfer"})).unwrap();
        assert_eq!(steady, Procedure::HeatTransfer { transient: None });
        assert!(steady.is_thermal() && steady.validate().is_ok());
        assert_eq!(steady.default_output(), [OutputVariable::Temperature, OutputVariable::HeatFlux].into_iter().collect());
        let t: Procedure =
            serde_json::from_value(serde_json::json!({"type": "heat_transfer", "transient": {"end_time": 10, "increment": 0.3}})).unwrap();
        let Procedure::HeatTransfer { transient: Some(tr) } = t else { panic!() };
        // From 0 °C by default.
        assert_eq!((tr.theta, tr.initial_temperature, tr.increments()), (1.0, 273.15, 34));
        assert_eq!(Transient::new(1.0, 0.25).increments(), 4);
        for bad in [
            Transient { theta: 0.4, ..Transient::new(1.0, 0.1) },
            Transient { theta: f64::NAN, ..Transient::new(1.0, 0.1) },
            Transient::new(1.0, 0.0),
            Transient::new(-1.0, 0.1),
            Transient::new(1e300, 1e-300),
            Transient { initial_temperature: f64::INFINITY, ..Transient::new(1.0, 0.1) },
        ] {
            assert!(Procedure::HeatTransfer { transient: Some(bad) }.validate().is_err(), "{bad:?}");
        }
        let coupled: Procedure = serde_json::from_value(serde_json::json!({"type": "coupled_temperature_displacement"})).unwrap();
        assert!(coupled.is_thermal() && coupled.default_output().contains(&OutputVariable::Stress));
        assert_eq!(coupled.name(), "coupled_temperature_displacement");
    }
}
