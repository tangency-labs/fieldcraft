//! Sampled numerical validity, separate from completion, statistical convergence and validation.
//!
//! Density departure and local Mach follow the weakly compressible scaling in Krüger et al.,
//! *The Lattice Boltzmann Method*, Springer (2017). The 2% / 0.2 envelope is the initial
//! low-Mach study policy of #785/#793, not a universal stability limit. Force peaks are
//! descriptive indicators: a real transient may have a large peak.

use crate::{Error, Fields, Result};

/// Maximum regular intervals in a runtime field window. Initial and early-stop/final samples
/// are included. No fields are retained: at most 258 small records.
pub const MAX_INTERVALS: u64 = 256;
pub const MAX_SAMPLES: usize = MAX_INTERVALS as usize + 2;
pub const DENSITY_ENVELOPE: f64 = 0.02;
pub const MACH_ENVELOPE: f64 = 0.2;

// Fixed operation order, with scaling to avoid squaring large input-derived values.
fn magnitude(v: [f64; 3]) -> f64 {
    if !v.iter().all(|x| x.is_finite()) {
        return f64::NAN;
    }
    let scale = v.iter().map(|x| x.abs()).fold(0.0, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    let u = v.map(|x| x / scale);
    (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt() * scale
}

/// Extrema over every included fluid cell, including inlet and outlet layers. Densities are
/// ratios to the reference supplied to `measure` (1 in production lattice units).
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct Sample {
    pub step: u64,
    /// The backend returned unusable fields. Extrema are zero placeholders in this record;
    /// the existing backend guard/steady monitor remains responsible for divergence handling.
    pub invalid: bool,
    pub density_min: f64,
    pub density_max: f64,
    pub density_departure: f64,
    pub mach_min: f64,
    pub mach_max: f64,
}

/// Measure without downcasting to output-frame precision. A spatial reference allows studies
/// to account for an intended pressure gradient. `include` must exclude solid and gas cells.
pub fn measure(step: u64, fields: &Fields, include: impl Fn(usize) -> bool, reference: impl Fn(usize) -> f64) -> Result<Sample> {
    if fields.density.len() != fields.velocity.len() {
        return Err(Error::Invalid("validity field lengths differ".into()));
    }
    let mut sample = Sample { step, density_min: f64::MAX, mach_min: f64::MAX, ..Sample::default() };
    let mut count = 0usize;
    for (cell, (&rho, u)) in fields.density.iter().zip(&fields.velocity).enumerate() {
        if !include(cell) {
            continue;
        }
        let reference = reference(cell);
        let ratio = rho / reference;
        let mach = magnitude(*u) * 3.0_f64.sqrt();
        if !reference.is_finite() || reference <= 0.0 || !rho.is_finite() || rho <= 0.0 || !ratio.is_finite() || ratio <= 0.0 || !mach.is_finite() {
            return Err(Error::Invalid("validity needs positive finite density/reference and finite velocity".into()));
        }
        sample.density_min = sample.density_min.min(ratio);
        sample.density_max = sample.density_max.max(ratio);
        sample.density_departure = sample.density_departure.max((ratio - 1.0).abs());
        sample.mach_min = sample.mach_min.min(mach);
        sample.mach_max = sample.mach_max.max(mach);
        count += 1;
    }
    if count == 0 {
        return Err(Error::Invalid("validity window has no included fluid cells".into()));
    }
    Ok(sample)
}

/// A peak relative to the norm of the signed mean force, not the mean magnitude. A zero mean
/// has no defined ratio. Opposite forces are never cancelled across wall groups.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct ForceIndicator {
    pub name: String,
    pub body: Option<String>,
    pub samples: usize,
    pub peak_over_mean: Option<f64>,
    pub largest_jump_over_mean: Option<f64>,
}

impl ForceIndicator {
    pub fn from_history(history: &crate::ForceHistory) -> Self {
        let n = history.force.len();
        let mut mean = [0.0; 3];
        let mut peak = 0.0_f64;
        let mut jump = 0.0_f64;
        let mut previous = None;
        for &f in &history.force {
            for (m, v) in mean.iter_mut().zip(f) {
                *m += v / n as f64;
            }
            peak = peak.max(magnitude(f));
            if let Some(p) = previous {
                let p: [f64; 3] = p;
                let difference = magnitude([f[0] - p[0], f[1] - p[1], f[2] - p[2]]);
                jump = if difference.is_finite() { jump.max(difference) } else { f64::INFINITY };
            }
            previous = Some(f);
        }
        let scale = magnitude(mean);
        let ratio = |value| {
            let ratio: f64 = value / scale;
            (scale > 0.0 && scale.is_finite() && ratio.is_finite() && history.force.iter().flatten().all(|x| x.is_finite())).then_some(ratio)
        };
        Self { name: history.name.clone(), body: history.body.clone(), samples: n, peak_over_mean: ratio(peak), largest_jump_over_mean: ratio(jump) }
    }
}

/// Numerical diagnostics only. `within_study_envelope` is about sampled fields, never a
/// validation certificate; excursions between samples remain unknown.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Report {
    pub field_every_steps: u64,
    pub force_every_steps: u64,
    pub max_field_samples: usize,
    pub max_field_sample_bytes: usize,
    pub density_reference: &'static str,
    pub invalid_field_samples: usize,
    pub extrema: Option<Sample>,
    pub within_study_envelope: Option<bool>,
    pub samples: Vec<Sample>,
    pub forces: Vec<ForceIndicator>,
}

impl Report {
    pub fn new(samples: Vec<Sample>, field_every_steps: u64, force_every_steps: u64, forces: &[crate::ForceHistory]) -> Result<Self> {
        if samples.len() > MAX_SAMPLES || field_every_steps == 0 || force_every_steps == 0 {
            return Err(Error::Invalid("validity report exceeds its sample cap or has zero cadence".into()));
        }
        if samples.windows(2).any(|s| s.first().zip(s.get(1)).is_some_and(|(a, b)| a.step >= b.step)) || samples.iter().any(|s| !s.valid()) {
            return Err(Error::Invalid("validity report has invalid samples".into()));
        }
        let mut usable = samples.iter().filter(|s| !s.invalid);
        let extrema = usable.next().cloned().map(|mut s| {
            for next in usable {
                s.step = next.step;
                s.density_min = s.density_min.min(next.density_min);
                s.density_max = s.density_max.max(next.density_max);
                s.density_departure = s.density_departure.max(next.density_departure);
                s.mach_min = s.mach_min.min(next.mach_min);
                s.mach_max = s.mach_max.max(next.mach_max);
            }
            s
        });
        let invalid_field_samples = samples.iter().filter(|s| s.invalid).count();
        let within_study_envelope =
            extrema.as_ref().filter(|_| invalid_field_samples == 0).map(|s| s.density_departure <= DENSITY_ENVELOPE && s.mach_max <= MACH_ENVELOPE);
        Ok(Self {
            field_every_steps,
            force_every_steps,
            max_field_samples: MAX_SAMPLES,
            max_field_sample_bytes: MAX_SAMPLES * std::mem::size_of::<Sample>(),
            density_reference: "material density (lattice density 1); intended pressure gradients are not subtracted",
            invalid_field_samples,
            extrema,
            within_study_envelope,
            samples,
            forces: forces.iter().map(ForceIndicator::from_history).collect(),
        })
    }
}

impl Sample {
    pub(crate) fn valid(&self) -> bool {
        if self.invalid {
            return [self.density_min, self.density_max, self.density_departure, self.mach_min, self.mach_max].iter().all(|x| *x == 0.0);
        }
        [self.density_min, self.density_max, self.density_departure, self.mach_min, self.mach_max].iter().all(|x| x.is_finite())
            && self.density_min > 0.0
            && self.density_max >= self.density_min
            && self.density_departure >= (self.density_min - 1.0).abs().max((self.density_max - 1.0).abs())
            && self.mach_min >= 0.0
            && self.mach_max >= self.mach_min
    }
}

pub(crate) fn capture(samples: &mut Vec<Sample>, flow: &mut dyn crate::Backend, step: u64) -> Result<()> {
    if samples.last().is_some_and(|s| s.step == step) {
        return Ok(());
    }
    if samples.len() >= MAX_SAMPLES {
        return Err(Error::Limit("validity sample cap reached".into()));
    }
    let fields = flow.fields()?;
    // Fully gaseous free-surface domains have no meaningful fluid-density envelope.
    let include = |cell| flow.is_fluid(cell) && fields.fill.get(cell).is_none_or(|f| *f > 0.0);
    if !(0..fields.density.len()).any(include) {
        return Ok(());
    }
    let sample = match measure(step, &fields, include, |_| 1.0) {
        Ok(sample) => sample,
        Err(Error::Invalid(_)) => Sample { step, invalid: true, ..Sample::default() },
        Err(e) => return Err(e),
    };
    samples.try_reserve_exact(1).map_err(|_| Error::Limit("not enough memory for validity samples".into()))?;
    samples.push(sample);
    Ok(())
}
