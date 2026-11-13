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

use fieldcraft_model::{BoxFace, BoxFaces, Collision, FaceCondition, Fluid, FluidDomain, Material, Model, Outflow, Procedure, Profile, Step};

use crate::lattice::CS2;
use crate::{Error, Flow, Limits, Result, Setup};

/// The duct's length in cells.
pub const CELLS: usize = 2048;
/// The pulse's starting centre (cells from x_min).
pub const START: f64 = 1000.0;
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
    let uniform = probe_series(None)?;
    let characteristic = probe_series(Some(sigma))?;

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
    let a = (-rate / 2.0).exp();
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
        let (x, l) = (t as f64, r.ln());
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
            re += x * wt.cos();
            im -= x * wt.sin();
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

/// The probe's density over ρ₀, less one, after each of [`STEPS`] steps of the duct with the
/// uniform outlet or with the characteristic outlet of σ.
fn probe_series(sigma: Option<f64>) -> Result<Vec<f64>> {
    let cs = CS2.sqrt();
    let faces = BoxFaces {
        x_min: FaceCondition::Velocity { velocity: [0.0; 3], profile: Profile::Uniform },
        x_max: FaceCondition::Pressure { pressure: 0.0, outflow: Outflow::Uniform, sigma: None },
        y_min: FaceCondition::Periodic {},
        y_max: FaceCondition::Periodic {},
        z_min: FaceCondition::Periodic {},
        z_max: FaceCondition::Periodic {},
    };
    let domain = FluidDomain {
        material: "Fluid".into(),
        min: [0.0; 3],
        max: [CELLS as f64, 1.0, 1.0],
        spacing: 1.0,
        faces,
        obstacles: vec![],
        surfaces: vec![],
        initial_velocity: [0.0; 3],
        acceleration: [0.0; 3],
        tolerance: None,
        thermal: None,
        porous: vec![],
        fans: vec![],
        baffles: vec![],
        free_surface: None,
    };
    let mut model = Model::default();
    let fluid = Material { density: Some(1.0), fluid: Some(Fluid::Newtonian { viscosity: (TAU - 0.5) / 3.0 }), ..Material::default() };
    model.add_material("Fluid", fluid)?;
    model.add_fluid_domain("Duct", domain)?;
    let procedure = Procedure::Fluid {
        domain: "Duct".into(),
        end_time: Some(STEPS as f64),
        steps: None,
        time_step: Some(1.0),
        frames: 1,
        magic: None,
        sigma: None,
        collision: Some(Collision::Trt),
        turbulence: None,
        steady: None,
    };
    model.add_step(Step::new("Run", procedure), None)?;
    let setup = Setup::new(&model, None, &[], &Limits::default())?;
    let mut flow = Flow::new(&setup, &[])?;
    if let Some(s) = sigma {
        flow.characteristic_outlet(BoxFace::XMax, s)?;
    }
    let grid = *flow.grid();
    flow.set_equilibrium(|cell| {
        let x = grid.centre(cell).first().copied().unwrap_or(0.0);
        let r = AMPLITUDE * (-(x - START).powi(2) / (2.0 * WIDTH * WIDTH)).exp();
        (1.0 + r, [cs * r, 0.0, 0.0])
    });
    let probe = grid.index([PROBE, 0, 0]);
    let mut series = Vec::with_capacity(STEPS);
    for _ in 0..STEPS {
        flow.step()?;
        let (rho, _) = flow.moments(probe).ok_or_else(|| Error::Invalid("the probe is not a fluid cell".into()))?;
        series.push(rho - 1.0);
    }
    Ok(series)
}
