//! The plane-pulse check of the characteristic outlet's reflection ([`outlet_reflection`]).
//!
//! A duct of [`CELLS`] × 1 × 1 cells (periodic sides, TRT τ [`TAU`]) has a velocity inlet at
//! rest at x_min and a uniform pressure outlet at x_max. A right-going Gaussian pulse
//! `ρ′ = A exp(−(x − x₀)² / 2w²)`, `u′ = c_s ρ′ / ρ₀`, starts at x₀ = [`START`], and a probe at
//! [`PROBE`] records the density each step: first the incident pulse, then what the outlet sends
//! back. The duct runs twice, with the uniform outlet and with the characteristic outlet of σ
//! ([`crate::Flow::characteristic_outlet`]; T. J. Poinsot, S. K. Lele, J. Comput. Phys.
//! 101:104–129, 1992, with σ ≈ 0.25 after D. H. Rudy, J. C. Strikwerda, J. Comput. Phys.
//! 36:55–70, 1980).
//!
//! **The theory.** The relaxed outlet reflects `R(ω) = −1 / (1 − 2iω/K)` of a normal wave, so
//! `|R| = K / √(K² + 4ω²)` (L. Selle, F. Nicoud, T. Poinsot, *Actual impedance of nonreflecting
//! boundary conditions: implications for computation of resonators*, AIAA J. 42(5):958–964,
//! 2004). In time, that is a one-pole low-pass of rate K/2 applied to the pulse a fully
//! reflecting outlet returns: the uniform outlet's reflection, filtered by
//! `y[t] = a y[t−1] + (1 − a) u[t]`, `a = e^(−K/2)`, predicts the characteristic outlet's. The
//! check reports the largest departure of the measured reflection from that prediction and the
//! decay rate of its tail over K/2, both up to the step the inlet's own reflection would return.
//!
//! **The windowed spectral ratio** ([`OutletReflection::windowed`]) is the ratio of the discrete
//! Fourier amplitudes of the reflected and the incident pulse, each over a rectangular window of
//! ±5w/c_s steps. The window cuts off the relaxation's slow tail (time constant 2/K), so it is
//! not the outlet's |R| at those frequencies: it is this duct's number, kept so that a change to
//! the outlet shows.
//!
//! **The layer check** ([`layer_reflection`]) runs the same pulse in a duct of [`LAYER_CELLS`] cells
//! with and without an absorbing layer ([`crate::absorbing`]) at its pressure end. The layer's
//! stated reflection R is that of the round trip through the layer alone; the duct's own pressure
//! face returns the rest, so the layered run is compared with the plain one and the ratio of the
//! two returned pulses is the layer's reflection. Its difference from R is the viscous part the
//! setup reports and the grid's own part, which theory does not give and the manual reports.

use fieldcraft_model::{AbsorbingLayer, BoxFace, BoxFaces, Collision, FaceCondition, Fluid, FluidDomain, Material, Model, Outflow, Procedure, Profile, Step};

use crate::absorbing::LayerReport;
use crate::lattice::CS2;
use crate::{Error, Flow, Limits, Result, Setup};

/// The duct's length in cells.
pub const CELLS: usize = 2048;
/// The pulse's starting centre (cells from x_min).
pub const START: f64 = 1000.0; // unit-factors: a position in cells
/// The probe's cell.
pub const PROBE: usize = 1300;
/// The pulse's width w (cells).
pub const WIDTH: f64 = 40.0;
/// The pulse's density amplitude over the reference density.
pub const AMPLITUDE: f64 = 1e-4;
/// The relaxation time.
pub const TAU: f64 = 0.51;
/// The steps each run takes: past the end of the record.
pub const STEPS: usize = 7400;
/// The most frequencies one check reports.
pub const MAX_FREQUENCIES: usize = 16;

/// A windowed spectral ratio at one angular frequency (radians per step).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Windowed {
    pub omega: f64,
    /// The measured ratio.
    pub ratio: f64,
}

/// What [`outlet_reflection`] measures.
#[derive(Debug, Clone, PartialEq)]
pub struct OutletReflection {
    pub sigma: f64,
    /// The relaxation rate K = σ c_s / L (per step).
    pub rate: f64,
    /// The steps of the record the check uses (before the inlet's reflection returns).
    pub record_steps: usize,
    /// The uniform outlet's reflected over incident peak: the control (a fully reflecting end).
    pub uniform_peak_ratio: f64,
    /// The characteristic outlet's incident pulse peak at the probe (density over ρ₀).
    pub incident_peak: f64,
    /// The largest departure of its reflection from the one-pole prediction, over the incident
    /// peak.
    pub largest_departure: f64,
    /// Its reflection's tail level after the pulse has passed, over the incident peak.
    pub tail_level: f64,
    /// The decay rate fitted to that tail, over K/2.
    pub fitted_rate: f64,
    /// The characteristic outlet's windowed spectral ratio (includes the truncated relaxation
    /// tail).
    pub windowed: Vec<Windowed>,
    /// The same window applied to the one-pole prediction.
    pub predicted_windowed: Vec<Windowed>,
    /// The uniform outlet's windowed spectral ratio.
    pub uniform_windowed: Vec<Windowed>,
}

/// Run the duct with the uniform outlet and with the characteristic outlet of `sigma`, and
/// compare the second's reflection with the one-pole response of the first ([module](self)).
/// `omegas` are the angular frequencies (radians per step) of the windowed spectral ratios.
pub fn outlet_reflection(sigma: f64, omegas: &[f64]) -> Result<OutletReflection> {
    if !(sigma > 0.0 && sigma <= 1.0) {
        return Err(Error::Invalid(format!("σ is {sigma}: give one above 0 and at most 1")));
    }
    if omegas.len() > MAX_FREQUENCIES {
        return Err(Error::Limit(format!("at most {MAX_FREQUENCIES} frequencies")));
    }
    if let Some(w) = omegas.iter().find(|w| !(w.is_finite() && **w > 0.0 && **w <= std::f64::consts::PI)) {
        return Err(Error::Invalid(format!("the frequency {w} is not above 0 and at most π radians per step")));
    }
    let cs = CS2.sqrt();
    let n = CELLS as f64;
    let rate = sigma * cs / n;
    let uniform = duct_series(&outlet_duct(None))?.0;
    let characteristic = duct_series(&outlet_duct(Some(sigma)))?.0;

    // The incident pulse passes the probe at t_inc, its reflection at t_ref; the windows span
    // ±5w/c_s steps; the record ends before the inlet's reflection of the reflection returns.
    let probe = PROBE as f64;
    let (t_inc, t_ref) = ((probe - START) / cs, ((n - START) + (n - probe)) / cs);
    let half = (5.0 * WIDTH / cs) as usize;
    let split = (0.5 * (t_inc + t_ref)) as usize;
    let end = ((t_ref + 2.0 * probe / cs) as usize).saturating_sub(half);
    let at = |v: &[f64], t: usize| v.get(t).copied().unwrap_or(0.0);
    if end >= characteristic.len() || split >= end {
        return Err(Error::Invalid("the record is shorter than the check needs".into()));
    }

    // The one-pole prediction from the uniform outlet's reflection.
    let a = libm::exp(-rate / 2.0);
    let mut predicted = vec![0.0; uniform.len()];
    let mut y = 0.0;
    for (t, p) in predicted.iter_mut().enumerate().skip(split) {
        y = a * y + (1.0 - a) * at(&uniform, t);
        *p = y;
    }
    let incident_peak = characteristic.iter().take(split).fold(0.0_f64, |m, v| m.max(v.abs()));
    if !(incident_peak > 0.0 && incident_peak.is_finite()) {
        return Err(Error::Unstable("the incident pulse did not reach the probe".into()));
    }
    let worst = (split..end).map(|t| (at(&characteristic, t) - at(&predicted, t)).abs()).fold(0.0_f64, f64::max);

    // The tail's decay rate: least squares of ln |r| on t, after the pulse has passed.
    let first = (t_ref as usize).saturating_add(half);
    let (mut m, mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for t in first..end {
        let r = at(&characteristic, t).abs();
        if !(r > 0.0 && r.is_finite()) {
            return Err(Error::Unstable(format!("the reflection's tail crosses zero at step {t}")));
        }
        let (x, l) = (t as f64, libm::log(r));
        m += 1.0;
        sx += x;
        sy += l;
        sxx += x * x;
        sxy += x * l;
    }
    let fitted = -(m * sxy - sx * sy) / (m * sxx - sx * sx);

    let window = |centre: f64| {
        let c0 = centre as usize;
        c0.saturating_sub(half)..(c0 + half).min(characteristic.len())
    };
    let (incident, reflected) = (window(t_inc), window(t_ref));
    let amplitude = |v: &[f64], span: std::ops::Range<usize>, omega: f64| {
        let (mut re, mut im) = (0.0, 0.0);
        for t in span {
            let (x, wt) = (at(v, t), omega * t as f64);
            re += x * libm::cos(wt);
            im -= x * libm::sin(wt);
        }
        (re * re + im * im).sqrt()
    };
    let ratios = |v: &[f64], of: &[f64]| -> Vec<Windowed> {
        omegas
            .iter()
            .map(|&omega| Windowed { omega, ratio: amplitude(v, reflected.clone(), omega) / amplitude(of, incident.clone(), omega) })
            .collect()
    };
    let peak = |v: &[f64], span: std::ops::Range<usize>| span.map(|t| at(v, t).abs()).fold(0.0_f64, f64::max);
    Ok(OutletReflection {
        sigma,
        rate,
        record_steps: end,
        uniform_peak_ratio: peak(&uniform, reflected.clone()) / peak(&uniform, incident.clone()),
        incident_peak,
        largest_departure: worst / incident_peak,
        tail_level: at(&characteristic, first).abs() / incident_peak,
        fitted_rate: fitted / (rate / 2.0),
        windowed: ratios(&characteristic, &characteristic),
        predicted_windowed: ratios(&predicted, &characteristic),
        uniform_windowed: ratios(&uniform, &uniform),
    })
}

/// The duct of the layer check: its length in cells, the pulse's starting centre and the probe, in cells.
pub const LAYER_CELLS: usize = 1024;
pub const LAYER_START: f64 = 500.0;
pub const LAYER_PROBE: usize = 700;
/// The largest stream the layer check takes, in lattice speed: the round trip's exponent is sound well below
/// the speed of sound.
const MAX_STREAM: f64 = 0.2;

/// What [`layer_reflection`] measures.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerReflection {
    /// The layer's setup report (its stated reflection, the predicted returns and the peak rate).
    pub report: LayerReport,
    /// The duct without the layer: the pressure face's own returned over incident peak.
    pub plain: f64,
    /// The layer's measured reflection: the layered run's returned over incident peak, over `plain`.
    pub measured: f64,
}

/// Run the duct without and with an absorbing layer of `cells` cells and stated `reflection` at its
/// pressure end, in a stream of lattice speed `stream` along +x, and compare the returned pulse. The
/// pulse runs right from [`LAYER_START`] and is read at [`LAYER_PROBE`]; the record ends after its
/// return and the window of its peaks (±5 w / c_s steps, as the outlet check's).
pub fn layer_reflection(collision: Collision, cells: usize, reflection: f64, stream: f64) -> Result<LayerReflection> {
    if !(stream.is_finite() && stream.abs() < MAX_STREAM) {
        return Err(Error::Invalid(format!("the stream is {stream}: give one finite and below {MAX_STREAM} in lattice speed")));
    }
    if cells == 0 || cells >= LAYER_PROBE / 2 {
        return Err(Error::Invalid(format!("the layer has {cells} cells: give 1 to {} (the probe's half of the duct)", LAYER_PROBE / 2 - 1)));
    }
    let cs = CS2.sqrt();
    let n = LAYER_CELLS as f64;
    let probe = LAYER_PROBE as f64;
    let (forward, backward) = (cs + stream, cs - stream);
    let t_inc = (probe - LAYER_START) / forward;
    let t_ref = t_inc + (n - probe) * (1.0 / forward + 1.0 / backward);
    let half = (5.0 * WIDTH / cs) as usize;
    let steps = t_ref as usize + half + 1;
    let duct = |layer: Option<(usize, f64)>| Duct {
        cells: LAYER_CELLS,
        start: LAYER_START,
        probe: LAYER_PROBE,
        steps,
        collision,
        stream,
        outlet_sigma: None,
        layer,
    };
    let (plain_series, _) = duct_series(&duct(None))?;
    let (layered_series, setup) = duct_series(&duct(Some((cells, reflection))))?;
    let report = setup.layers.first().map(|l| l.report).ok_or_else(|| Error::Invalid("the layer check's duct has no layer".into()))?;
    let incident_plain = peak_about(&plain_series, t_inc, half);
    let incident_layered = peak_about(&layered_series, t_inc, half);
    for (incident, what) in [(incident_plain, "plain"), (incident_layered, "layered")] {
        if !(incident > 0.0 && incident.is_finite()) {
            return Err(Error::Unstable(format!("the {what} incident pulse did not reach the probe")));
        }
    }
    let plain = peak_about(&plain_series, t_ref, half) / incident_plain;
    let layered = peak_about(&layered_series, t_ref, half) / incident_layered;
    Ok(LayerReflection { report, plain, measured: layered / plain })
}

/// The largest |ρ − 1| of a series within `half` steps of `centre`.
fn peak_about(series: &[f64], centre: f64, half: usize) -> f64 {
    let c0 = centre as usize;
    let last = series.len().saturating_sub(1);
    (c0.saturating_sub(half)..=c0.saturating_add(half).min(last)).filter_map(|t| series.get(t)).fold(0.0_f64, |m, v| m.max(v.abs()))
}

/// One run of the duct: its cells, the pulse's starting centre and the probe (cells), the steps, the collision,
/// the stream along +x (lattice speed), the characteristic outlet's σ at x_max (none: uniform), and an absorbing
/// layer at x_max as (cells, reflection) (none: no layer).
struct Duct {
    cells: usize,
    start: f64,
    probe: usize,
    steps: usize,
    collision: Collision,
    stream: f64,
    outlet_sigma: Option<f64>,
    layer: Option<(usize, f64)>,
}

/// The duct of the outlet check: TRT, at rest, with the uniform outlet or the characteristic one of σ.
fn outlet_duct(sigma: Option<f64>) -> Duct {
    Duct {
        cells: CELLS,
        start: START,
        probe: PROBE,
        steps: STEPS,
        collision: Collision::Trt,
        stream: 0.0,
        outlet_sigma: sigma,
        layer: None,
    }
}

/// The probe's density over ρ₀, less one, after each step of the duct, and the setup it ran with.
fn duct_series(duct: &Duct) -> Result<(Vec<f64>, Setup)> {
    let cs = CS2.sqrt();
    let stream = [duct.stream, 0.0, 0.0];
    let faces = BoxFaces {
        x_min: FaceCondition::Velocity { velocity: stream, profile: Profile::Uniform },
        x_max: FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None },
        y_min: FaceCondition::Periodic {},
        y_max: FaceCondition::Periodic {},
        z_min: FaceCondition::Periodic {},
        z_max: FaceCondition::Periodic {},
    };
    let layers = duct
        .layer
        .map(|(cells, reflection)| AbsorbingLayer { face: BoxFace::XMax, width: cells as f64, reflection, reference_velocity: Some(stream) })
        .into_iter()
        .collect();
    let domain = FluidDomain {
        material: "Fluid".into(),
        min: [0.0; 3],
        max: [duct.cells as f64, 1.0, 1.0],
        spacing: 1.0,
        faces,
        obstacles: vec![],
        surfaces: vec![],
        initial_velocity: stream,
        acceleration: [0.0; 3],
        tolerance: None,
        thermal: None,
        porous: vec![],
        fans: vec![],
        baffles: vec![],
        free_surface: None,
        layers,
    };
    let mut model = Model::default();
    let fluid = Material { density: Some(1.0), fluid: Some(Fluid::Newtonian { viscosity: (TAU - 0.5) / 3.0 }), ..Material::default() };
    model.add_material("Fluid", fluid)?;
    model.add_fluid_domain("Duct", domain)?;
    let procedure = Procedure::Fluid {
        domain: "Duct".into(),
        end_time: Some(duct.steps as f64),
        steps: None,
        time_step: Some(1.0),
        frames: 1,
        magic: None,
        sigma: None,
        collision: Some(duct.collision),
        turbulence: None,
        steady: None,
        walls: Default::default(),
        inlet_ramp: None,
    };
    model.add_step(Step::new("Run", procedure), None)?;
    let setup = Setup::new(&model, None, &[], &Limits::default())?;
    let mut flow = Flow::new(&setup, &[])?;
    if let Some(s) = duct.outlet_sigma {
        flow.characteristic_outlet(BoxFace::XMax, s)?;
    }
    let grid = *flow.grid();
    let (start, stream_speed) = (duct.start, duct.stream);
    flow.set_equilibrium(|cell| {
        let x = grid.centre(cell).first().copied().unwrap_or(0.0);
        let r = AMPLITUDE * libm::exp(-(x - start).powi(2) / (2.0 * WIDTH * WIDTH));
        (1.0 + r, [stream_speed + cs * r, 0.0, 0.0])
    });
    let probe = grid.index([duct.probe, 0, 0]);
    let mut series = Vec::with_capacity(duct.steps);
    for _ in 0..duct.steps {
        fieldcraft_model::progress::check()?;
        flow.step()?;
        let (rho, _) = flow.moments(probe).ok_or_else(|| Error::Invalid("the probe is not a fluid cell".into()))?;
        series.push(rho - 1.0);
    }
    Ok((series, setup))
}
