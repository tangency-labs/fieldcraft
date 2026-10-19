//! Mesh-wide quality: statistics per metric, and the validity gate meshers run on their output.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{Element, ElementKind, Model, check_coordinates};

use super::{ElementQuality, JacobianCheck, Metric, Validity, element_quality, element_validity, is_measured};
use crate::{Error, Result};

/// The most histogram bins a report may ask for.
pub const MAX_BINS: usize = 1000;

/// Statistics of one metric over the measured elements.
#[derive(Debug, Clone, PartialEq)]
pub struct MetricStats {
    pub metric: Metric,
    /// Elements that have this metric.
    pub count: usize,
    pub min: f64,
    pub max: f64,
    /// Mean of the finite values (`None` when there are none).
    pub mean: Option<f64>,
    /// The element with the worst value (the first one on ties).
    pub worst: u32,
    /// Counts over [`Metric::range`] in equal bins; values outside the range fall in the end bins.
    pub histogram: Vec<usize>,
}

/// Quality of a set of elements.
#[derive(Debug, Clone, PartialEq)]
pub struct QualityReport {
    /// Elements measured (solids and shells).
    pub measured: usize,
    /// Elements of kinds with no quality metrics (rods, beams), by kind.
    pub skipped: Vec<(ElementKind, usize)>,
    /// Solids proven valid.
    pub valid: usize,
    pub invalid: Vec<u32>,
    pub undetermined: Vec<u32>,
    /// Total signed volume of the solids and total area of the shells.
    pub volume: f64,
    pub area: f64,
    /// One entry per metric present, in [`Metric`] declaration order.
    pub metrics: Vec<MetricStats>,
}

fn coordinates(model: &Model, label: u32, element: &Element) -> Result<Vec<[f64; 3]>> {
    element
        .nodes
        .iter()
        .map(|n| {
            model.node(*n).map(|node| node.xyz).ok_or_else(|| Error::Model(fieldcraft_model::Error::NotFound(format!("element {label}: node {n}"))))
        })
        .collect()
}

/// The elements to look at: all of them, or the given labels (sorted, without repeats).
fn selection<'a>(model: &'a Model, labels: Option<&[u32]>) -> Result<Vec<(u32, &'a Element)>> {
    match labels {
        None => Ok(model.elements().iter().map(|(l, e)| (*l, e)).collect()),
        Some(labels) => labels
            .iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .map(|l| model.element(*l).map(|e| (*l, e)).ok_or_else(|| Error::Model(fieldcraft_model::Error::NotFound(format!("element {l}")))))
            .collect(),
    }
}

/// Quality statistics of the model's elements, or of the given ones, with `bins` histogram bins
/// (1..=[`MAX_BINS`]) per metric. Deterministic: elements are visited in label order.
pub fn quality_report(model: &Model, labels: Option<&[u32]>, bins: usize) -> Result<QualityReport> {
    if !(1..=MAX_BINS).contains(&bins) {
        return Err(Error::Invalid(format!("bins must be in 1..={MAX_BINS}, got {bins}")));
    }
    let mut report = QualityReport {
        measured: 0,
        skipped: Vec::new(),
        valid: 0,
        invalid: Vec::new(),
        undetermined: Vec::new(),
        volume: 0.0,
        area: 0.0,
        metrics: Vec::new(),
    };
    let mut skipped: BTreeMap<ElementKind, usize> = BTreeMap::new();
    let mut stats: BTreeMap<Metric, (MetricStats, f64, usize)> = BTreeMap::new();
    for (label, element) in selection(model, labels)? {
        if !is_measured(element.kind) {
            *skipped.entry(element.kind).or_insert(0) += 1;
            continue;
        }
        let q: ElementQuality = element_quality(element.kind, &coordinates(model, label, element)?)?;
        report.measured += 1;
        match q.validity() {
            Some(Validity::Valid) => report.valid += 1,
            Some(Validity::Invalid) => report.invalid.push(label),
            Some(Validity::Undetermined) => report.undetermined.push(label),
            None => {}
        }
        if element.kind.dimension() == 3 {
            report.volume += q.size;
        } else {
            report.area += q.size;
        }
        for (metric, v) in q.metrics {
            let (s, sum, finite) = stats.entry(metric).or_insert_with(|| {
                (
                    MetricStats { metric, count: 0, min: f64::INFINITY, max: f64::NEG_INFINITY, mean: None, worst: label, histogram: vec![0; bins] },
                    0.0,
                    0,
                )
            });
            if v.is_nan() {
                continue;
            }
            let worse = if metric.higher_is_better() { v < s.min } else { v > s.max };
            if worse || s.count == 0 {
                s.worst = label;
            }
            s.count += 1;
            s.min = s.min.min(v);
            s.max = s.max.max(v);
            if v.is_finite() {
                *sum += v;
                *finite += 1;
            }
            let (lo, hi) = metric.range();
            // `as` saturates: −∞ lands in bin 0, +∞ in the last bin.
            let bin = (((v - lo) / (hi - lo)) * bins as f64).floor().max(0.0) as usize;
            if let Some(b) = s.histogram.get_mut(bin.min(bins - 1)) {
                *b += 1;
            }
        }
    }
    report.skipped = skipped.into_iter().collect();
    report.metrics = stats
        .into_values()
        .map(|(mut s, sum, finite)| {
            s.mean = (finite > 0).then(|| sum / finite as f64);
            s
        })
        .collect();
    Ok(report)
}

/// An element that fails a quality limit, or is not proven valid.
#[derive(Debug, Clone, PartialEq)]
pub struct Marked {
    pub label: u32,
    pub kind: ElementKind,
    /// The mean of its nodes.
    pub centroid: [f64; 3],
    /// The certified validity, for solids.
    pub validity: Option<Validity>,
    /// Each limit it fails: the metric, its value and the limit.
    pub failed: Vec<(Metric, f64, f64)>,
    /// Its scaled Jacobian, when it has one (the worst-first order).
    pub scaled_jacobian: Option<f64>,
}

impl Marked {
    /// Not proven valid (inverted or undetermined).
    pub fn is_invalid(&self) -> bool {
        matches!(self.validity, Some(Validity::Invalid | Validity::Undetermined))
    }
}

/// The elements (all, or the given labels) that fail a limit: a value below the limit of a
/// metric where higher is better ([`Metric::higher_is_better`]), above it otherwise, or no value
/// (NaN). Solids not proven valid are always marked. A metric an element does not have does not
/// apply to it. Returns how many elements were measured and the marked ones, worst first: not
/// valid, then by scaled Jacobian (lowest first; none counts as 1), then by label.
pub fn mark_elements(model: &Model, labels: Option<&[u32]>, limits: &[(Metric, f64)]) -> Result<(usize, Vec<Marked>)> {
    if let Some((m, v)) = limits.iter().find(|(_, v)| !v.is_finite()) {
        return Err(Error::Invalid(format!("the limit of {} must be a finite number, got {v}", m.name())));
    }
    let mut measured = 0;
    let mut marked = Vec::new();
    for (label, element) in selection(model, labels)? {
        if !is_measured(element.kind) {
            continue;
        }
        let xyz = coordinates(model, label, element)?;
        let q = element_quality(element.kind, &xyz)?;
        measured += 1;
        let failed: Vec<(Metric, f64, f64)> = limits
            .iter()
            .filter_map(|&(metric, limit)| {
                let v = q.get(metric)?;
                let fails = v.is_nan() || if metric.higher_is_better() { v < limit } else { v > limit };
                fails.then_some((metric, v, limit))
            })
            .collect();
        let validity = q.validity();
        let invalid = matches!(validity, Some(Validity::Invalid | Validity::Undetermined));
        if failed.is_empty() && !invalid {
            continue;
        }
        let n = xyz.len().max(1) as f64;
        let centroid = [0, 1, 2].map(|k| xyz.iter().map(|p| p[k]).sum::<f64>() / n);
        marked.push(Marked { label, kind: element.kind, centroid, validity, failed, scaled_jacobian: q.get(Metric::ScaledJacobian) });
    }
    let key = |m: &Marked| (!m.is_invalid(), m.scaled_jacobian.filter(|v| !v.is_nan()).unwrap_or(1.0));
    marked.sort_by(|a, b| {
        let (ka, kb) = (key(a), key(b));
        ka.0.cmp(&kb.0).then(ka.1.total_cmp(&kb.1)).then(a.label.cmp(&b.label))
    });
    Ok((measured, marked))
}

/// What the validity gate checks besides element validity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GateOptions {
    /// The volume the mesh must fill: the boundary's enclosed volume ([`enclosed_volume`]).
    pub expected_volume: Option<f64>,
    /// Relative tolerance on the volume. Curved (quadratic) meshes of a faceted boundary need a
    /// looser one than linear meshes.
    pub volume_tolerance: f64,
}

impl Default for GateOptions {
    fn default() -> Self {
        Self { expected_volume: None, volume_tolerance: 1e-9 }
    }
}

/// The result of the validity gate.
#[derive(Debug, Clone, PartialEq)]
pub struct GateReport {
    /// Solid elements checked (other kinds are ignored).
    pub checked: usize,
    pub invalid: Vec<u32>,
    pub undetermined: Vec<u32>,
    /// Total signed volume of the solids.
    pub volume: f64,
    pub options: GateOptions,
}

impl GateReport {
    /// Whether the volume matches the expected one (true when none is expected).
    pub fn volume_matches(&self) -> bool {
        self.options.expected_volume.is_none_or(|v| (self.volume - v).abs() <= self.options.volume_tolerance * v.abs())
    }

    /// Every solid is proven valid and the volume matches.
    pub fn passed(&self) -> bool {
        self.invalid.is_empty() && self.undetermined.is_empty() && self.volume_matches()
    }

    /// The report, or [`Error::Failed`] describing why the gate failed.
    pub fn into_result(self) -> Result<Self> {
        if self.passed() {
            return Ok(self);
        }
        let mut why: Vec<crate::message::Message> = Vec::new();
        let list = |v: &[u32]| v.iter().take(10).map(u32::to_string).collect::<Vec<_>>().join(", ");
        if !self.invalid.is_empty() {
            why.push(crate::message![format!("{} inverted or degenerate elements ({})", self.invalid.len(), list(&self.invalid))]);
        }
        if !self.undetermined.is_empty() {
            why.push(crate::message![format!(
                "{} elements whose validity could not be proven ({})",
                self.undetermined.len(),
                list(&self.undetermined)
            )]);
        }
        if let (false, Some(v)) = (self.volume_matches(), self.options.expected_volume) {
            why.push(crate::message![
                "mesh volume ",
                crate::message::volume(self.volume),
                " differs from the boundary's ",
                crate::message::volume(v)
            ]);
        }
        let mut message = crate::message!["mesh validity gate: "];
        for (i, w) in why.into_iter().enumerate() {
            if i > 0 {
                message.0.push("; ".into());
            }
            message.0.extend(w.0);
        }
        Err(Error::failed(message))
    }
}

/// Accumulates the gate's report one element at a time.
struct Gate(GateReport);

impl Gate {
    fn new(options: &GateOptions) -> Result<Self> {
        if !(options.volume_tolerance >= 0.0 && options.volume_tolerance.is_finite()) || options.expected_volume.is_some_and(|v| !v.is_finite()) {
            return Err(Error::Invalid("gate options must be finite, with a non-negative tolerance".into()));
        }
        Ok(Gate(GateReport { checked: 0, invalid: Vec::new(), undetermined: Vec::new(), volume: 0.0, options: *options }))
    }

    fn add(&mut self, label: u32, kind: ElementKind, xyz: &[[f64; 3]]) -> Result<()> {
        if kind.dimension() != 3 {
            return Ok(());
        }
        self.add_check(label, &element_validity(kind, xyz)?);
        Ok(())
    }

    fn add_check(&mut self, label: u32, check: &JacobianCheck) {
        let r = &mut self.0;
        r.checked += 1;
        r.volume += check.volume;
        match check.validity {
            Validity::Valid => {}
            Validity::Invalid => r.invalid.push(label),
            Validity::Undetermined => r.undetermined.push(label),
        }
    }
}

/// [`validity_gate_elements`] over certified checks already made, for a mesher that measured its
/// elements with [`super::element_quality`] (which certifies them too) and should not certify
/// them twice: `(label, check)` per solid element.
pub fn validity_gate_checks<'a, I>(checks: I, options: &GateOptions) -> Result<GateReport>
where
    I: IntoIterator<Item = (u32, &'a JacobianCheck)>,
{
    let mut gate = Gate::new(options)?;
    for (label, check) in checks {
        gate.add_check(label, check);
    }
    Ok(gate.0)
}

/// The validity gate every mesher runs on its output before it reaches the model: each solid
/// element's Jacobian determinant is proven positive everywhere, and the elements' total volume
/// equals the volume the boundary encloses. Elements come as `(label, kind, node coordinates in
/// node order)`, straight from the mesher's arrays; other than solids are ignored. Call
/// [`GateReport::into_result`] to turn a failure into an error.
pub fn validity_gate_elements<I, E>(elements: I, options: &GateOptions) -> Result<GateReport>
where
    I: IntoIterator<Item = (u32, ElementKind, E)>,
    E: AsRef<[[f64; 3]]>,
{
    let mut gate = Gate::new(options)?;
    for (label, kind, xyz) in elements {
        gate.add(label, kind, xyz.as_ref())?;
    }
    Ok(gate.0)
}

/// [`validity_gate_elements`] over a model's elements.
pub fn validity_gate(model: &Model, options: &GateOptions) -> Result<GateReport> {
    let mut gate = Gate::new(options)?;
    for (label, element) in model.elements() {
        gate.add(*label, element.kind, &coordinates(model, *label, element)?)?;
    }
    Ok(gate.0)
}

/// The signed volume enclosed by a closed triangle surface (positive when the triangles' right-
/// hand normals point outwards), by the divergence theorem.
pub fn enclosed_volume(points: &[[f64; 3]], triangles: &[[usize; 3]]) -> Result<f64> {
    check_coordinates(points.as_flattened())?;
    let o = points.first().copied().unwrap_or([0.0; 3]);
    let mut six_v = 0.0;
    for t in triangles {
        let p = |i: usize| points.get(i).map(|p| super::shape::sub(*p, o)).ok_or_else(|| Error::Invalid(format!("triangle {t:?}: no point {i}")));
        six_v += super::shape::det(p(t[0])?, p(t[1])?, p(t[2])?);
    }
    Ok(six_v / 6.0)
}
