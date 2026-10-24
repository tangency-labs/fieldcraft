//! Read-only solid-field path sampling in SI. The spatial index retains only element labels
//! and conservative bounds, not copies of mesh coordinates or results. One cached index is
//! shared by session snapshots and reused while the exact result and its definition match.

use std::sync::{Arc, Weak};

use fieldcraft_fem::locate::{Bounds, Work, bounds, locate};
use fieldcraft_model::{ElementKind, Model};

use crate::results::{Field, Pick};
use crate::{EngineError, Result, StepResult};

pub const MAX_POINTS: usize = 1000;
pub const MAX_SAMPLES: usize = 10_000;
const MAX_ELEMENTS: usize = 500_000;
const MAX_RESULT_NODES: usize = 2_000_000;
const MAX_INDEX_BYTES: usize = 96 * 1024 * 1024;
const MAX_CANDIDATES: usize = 4096;
const MAX_VISITS: usize = 2_000_000;
const MAX_WORK: usize = 4_000_000;
const MAX_BUILD_WORK: usize = 14_000_000;
const LEAF: usize = 8;

fn error(why: impl Into<String>) -> EngineError {
    EngineError::invalid("points", why)
}

#[derive(Debug, Clone)]
struct Entry {
    label: u32,
    bounds: Bounds,
}

#[derive(Debug, Clone)]
struct Branch {
    bounds: Bounds,
    start: usize,
    end: usize,
    children: Option<(usize, usize)>,
}

#[derive(Debug)]
pub struct Index {
    entries: Vec<Entry>,
    branches: Vec<Branch>,
}

/// One cache entry keeps at most 96 MiB of checked index allocations. A result's existing
/// nodal vectors are borrowed; queries separately cap result nodes, samples and total work.
#[derive(Debug, Clone, Default)]
pub struct Cache(Option<Cached>);

#[derive(Debug, Clone)]
struct Cached {
    result: Weak<StepResult>,
    definition: [u8; 32],
    index: Arc<Index>,
}

impl Cache {
    pub fn index(&mut self, model: &Model, result: &Arc<StepResult>, definition: [u8; 32]) -> Result<Arc<Index>> {
        if let Some(cached) = &self.0
            && cached.definition == definition
            && cached.result.upgrade().is_some_and(|r| Arc::ptr_eq(&r, result))
        {
            return Ok(cached.index.clone());
        }
        let index = Arc::new(Index::build(model, result.nodes())?);
        self.0 = Some(Cached { result: Arc::downgrade(result), definition, index: index.clone() });
        Ok(index)
    }
}

/// The last successful query, in SI and reference geometry. It is session view state,
/// not a model edit or journal entry. The exact result identity guards stale overlays.
#[derive(Debug, Clone)]
pub struct Plot {
    pub step: String,
    pub field: Field,
    pub pick: Pick,
    /// Original query mode, None for non-frequency results even after the step is replaced.
    pub mode: Option<usize>,
    pub points: Vec<[f64; 3]>,
    pub distance: Vec<f64>,
    pub xyz: Vec<[f64; 3]>,
    pub values: Vec<Option<f64>>,
    pub elements: Vec<Option<u32>>,
    /// At most 10000 triples; the exact accepted location, for deformation of the overlay.
    pub natural: Vec<Option<[f64; 3]>>,
    pub result: Weak<StepResult>,
}

impl Plot {
    pub fn current(&self, session: &crate::Session) -> bool {
        session.result(&self.step).is_some_and(|r| self.result.upgrade().is_some_and(|a| Arc::ptr_eq(&a, r)))
    }

    /// The stored reference station plus the displayed frame/mode's interpolated displacement.
    /// `factor` is the existing viewport deformation scale * deformation_unit * animation swing.
    /// No guessed extrapolation is drawn across gaps; unrelated/stale views have no overlay.
    pub fn displayed_point(&self, session: &crate::Session, index: usize, display_step: &str, pick: Pick, factor: f64) -> Option<[f64; 3]> {
        if display_step != self.step || !self.current(session) || !factor.is_finite() {
            return None;
        }
        let point = *self.xyz.get(index)?;
        if factor == 0.0 {
            return Some(point);
        }
        let result = session.result(&self.step)?;
        let displacement = crate::results::displacements(result, pick).ok()?;
        let label = self.elements.get(index).copied().flatten()?;
        let natural = self.natural.get(index).copied().flatten()?;
        let element = session.model.element(label)?;
        let shape = fieldcraft_fem::Basis::of(element.kind).eval(natural);
        let mut offset = [0.0; 3];
        for (n, node) in shape.values().iter().zip(&element.nodes) {
            if *n == 0.0 {
                continue;
            }
            let k = result.nodes().binary_search(node).ok()?;
            let value = displacement.get(k)?;
            if value.iter().any(|v| !v.is_finite()) {
                return None;
            }
            for a in 0..3 {
                offset[a] += n * value[a];
            }
        }
        let shown = [0, 1, 2].map(|k| point[k] + factor * offset[k]);
        shown.iter().all(|v| v.is_finite()).then_some(shown)
    }
}

/// Map a genuine displayed face-triangle hit back to its element's reference surface.
/// Barycentric interpolation is in *natural* coordinates, then the full FE map is evaluated:
/// quadratic surface picks stay on the FE face instead of its piecewise-linear drawing chord.
pub fn reference_surface_point(model: &Model, label: u32, triangle_nodes: [u32; 3], barycentric: [f64; 3]) -> Option<[f64; 3]> {
    let tolerance = 256.0 * f64::EPSILON;
    if barycentric.iter().any(|v| !v.is_finite() || *v < -tolerance || *v > 1.0 + tolerance)
        || (barycentric.iter().sum::<f64>() - 1.0).abs() > tolerance
    {
        return None;
    }
    let element = model.element(label)?;
    let basis = fieldcraft_fem::Basis::of(element.kind);
    if basis.dim() != 3 || element.nodes.len() != basis.count() {
        return None;
    }
    let mut natural = [0.0; 3];
    for (node, weight) in triangle_nodes.iter().zip(barycentric) {
        let i = element.nodes.iter().position(|n| n == node)?;
        let xi = basis.nodes().get(i)?;
        for k in 0..3 {
            natural[k] += weight * xi[k];
        }
    }
    let shape = basis.eval(natural);
    let origin = model.node(*element.nodes.first()?)?.xyz;
    let mut offset = [0.0; 3];
    for (n, node) in shape.values().iter().zip(&element.nodes) {
        let p = model.node(*node)?.xyz;
        for k in 0..3 {
            offset[k] += n * (p[k] - origin[k]);
        }
    }
    let point = [0, 1, 2].map(|k| origin[k] + offset[k]);
    point.iter().all(|v| v.is_finite()).then_some(point)
}

fn geometry(model: &Model, label: u32) -> Result<(ElementKind, [[f64; 3]; 27], usize)> {
    let element = model.elements().get(&label).ok_or_else(|| error("path index element is missing; rebuild the result"))?;
    let mut xyz = [[0.0; 3]; 27];
    if element.nodes.len() > xyz.len() {
        return Err(error("path element node count exceeds the supported basis"));
    }
    for (out, node) in xyz.iter_mut().zip(&element.nodes) {
        *out = model.nodes().get(node).ok_or_else(|| error("path element refers to a missing mesh node"))?.xyz;
    }
    Ok((element.kind, xyz, element.nodes.len()))
}

impl Index {
    fn build(model: &Model, nodes: &[u32]) -> Result<Self> {
        fieldcraft_model::progress::check()?;
        if model.elements().len() > MAX_ELEMENTS || nodes.len() > MAX_RESULT_NODES {
            return Err(error(format!(
                "path index supports at most {MAX_ELEMENTS} mesh elements and {MAX_RESULT_NODES} result nodes; use a smaller model"
            )));
        }
        if nodes.windows(2).any(|w| matches!(w,[a,b] if a>=b)) {
            return Err(error("path result node labels must be strictly increasing"));
        }
        let count = model.elements().len();
        // Balanced splits of a range larger than 8 leave at least 4 elements per leaf.
        let max_branches = count.div_ceil(LEAF / 2).saturating_mul(2).max(1);
        let bytes = count
            .checked_mul(std::mem::size_of::<Entry>())
            .and_then(|b| max_branches.checked_mul(std::mem::size_of::<Branch>()).and_then(|n| n.checked_add(b)));
        if bytes.is_none_or(|b| b > MAX_INDEX_BYTES) {
            return Err(error("path index memory limit exceeded; use a smaller model"));
        }
        let mut entries = Vec::new();
        entries.try_reserve_exact(count).map_err(|_| error("not enough memory for the bounded path index"))?;
        fieldcraft_model::progress::stage("Path index bounds", Some(count as f64), "elements");
        let mut remaining_build_work = MAX_BUILD_WORK;
        for (i, (&label, element)) in model.elements().iter().enumerate() {
            if i.is_multiple_of(128) {
                fieldcraft_model::progress::check()?;
                fieldcraft_model::progress::count(i as f64);
            }
            if !element.nodes.iter().all(|n| nodes.binary_search(n).is_ok()) {
                continue;
            }
            let (kind, xyz, len) = geometry(model, label)?;
            remaining_build_work = remaining_build_work
                .checked_sub(if kind == ElementKind::Hex20 { 28 } else { 1 })
                .ok_or_else(|| error("path index construction work limit reached; use a smaller model"))?;
            let enclosure =
                bounds(kind, xyz.get(..len).ok_or_else(|| error("invalid path node count"))?).map_err(|e| error(format!("element {label}: {e}")))?;
            entries.push(Entry { label, bounds: enclosure });
        }
        if entries.is_empty() {
            return Err(error("the solved step has no supported solid elements with complete nodal data"));
        }
        let mut index = Self { entries, branches: Vec::new() };
        index.branches.try_reserve_exact(max_branches).map_err(|_| error("not enough memory for the bounded path index"))?;
        fieldcraft_model::progress::stage("Path index hierarchy", Some(index.entries.len() as f64), "elements");
        index.split(0, index.entries.len(), 0)?;
        fieldcraft_model::progress::check()?;
        Ok(index)
    }

    fn split(&mut self, start: usize, end: usize, depth: usize) -> Result<usize> {
        fieldcraft_model::progress::check()?;
        if depth > 20 {
            return Err(error("path index depth limit reached"));
        }
        let range = self.entries.get(start..end).filter(|e| !e.is_empty()).ok_or_else(|| error("invalid path index range"))?;
        let first = range.first().ok_or_else(|| error("empty path index range"))?;
        let enclosure = range.iter().skip(1).fold(first.bounds, |a, b| a.union(b.bounds));
        let id = self.branches.len();
        self.branches.push(Branch { bounds: enclosure, start, end, children: None });
        if end - start > LEAF {
            let axis = (0..3).max_by(|a, b| (enclosure.max[*a] - enclosure.min[*a]).total_cmp(&(enclosure.max[*b] - enclosure.min[*b]))).unwrap_or(0);
            let mid = start + (end - start) / 2;
            let range = self.entries.get_mut(start..end).ok_or_else(|| error("invalid path index split"))?;
            let center = |entry: &Entry| entry.bounds.min[axis] + (entry.bounds.max[axis] - entry.bounds.min[axis]) * 0.5;
            range.select_nth_unstable_by(mid - start, |a, b| center(a).total_cmp(&center(b)).then_with(|| a.label.cmp(&b.label)));
            let left = self.split(start, mid, depth + 1)?;
            let right = self.split(mid, end, depth + 1)?;
            self.branches.get_mut(id).ok_or_else(|| error("invalid path index branch"))?.children = Some((left, right));
        } else {
            fieldcraft_model::progress::count(end as f64);
        }
        Ok(id)
    }

    fn candidates(&self, point: [f64; 3], visits: &mut usize) -> Result<Vec<u32>> {
        let mut stack = Vec::with_capacity(22);
        stack.push(0);
        let mut candidates = Vec::new();
        while let Some(id) = stack.pop() {
            if (*visits).is_multiple_of(256) {
                fieldcraft_model::progress::check()?;
            }
            *visits = visits
                .checked_add(1)
                .filter(|v| *v <= MAX_VISITS)
                .ok_or_else(|| error("path index traversal limit reached; reduce samples or split the path"))?;
            let branch = self.branches.get(id).ok_or_else(|| error("invalid path index branch"))?;
            if !branch.bounds.contains_near(point) {
                continue;
            }
            if let Some((left, right)) = branch.children {
                stack.push(right);
                stack.push(left);
            } else {
                for entry in self.entries.get(branch.start..branch.end).ok_or_else(|| error("invalid path index leaf"))? {
                    if entry.bounds.contains_near(point) {
                        if candidates.len() >= MAX_CANDIDATES {
                            return Err(error("too many overlapping path elements; refine the model or move the path"));
                        }
                        candidates.push(entry.label);
                    }
                }
            }
        }
        // First containing candidate is the lowest label, independent of BVH traversal.
        candidates.sort_unstable();
        Ok(candidates)
    }

    fn station(
        &self,
        model: &Model,
        nodes: &[u32],
        values: &[f64],
        point: [f64; 3],
        work: &mut Work,
        visits: &mut usize,
    ) -> Result<(Option<u32>, Option<f64>, Option<[f64; 3]>)> {
        for label in self.candidates(point, visits)? {
            fieldcraft_model::progress::check()?;
            let (kind, xyz, len) = geometry(model, label)?;
            if let Some(location) = locate(kind, xyz.get(..len).ok_or_else(|| error("invalid path node count"))?, point, work)
                .map_err(|e| error(format!("element {label}: {e}")))?
            {
                let element = model.elements().get(&label).ok_or_else(|| error("path element missing"))?;
                let mut nodal = [0.0; 27];
                for (out, n) in nodal.iter_mut().zip(&element.nodes) {
                    *out = nodes.binary_search(n).ok().and_then(|k| values.get(k)).copied().unwrap_or(f64::NAN);
                }
                return Ok((
                    Some(label),
                    location.interpolate(nodal.get(..len).ok_or_else(|| error("invalid path field count"))?),
                    Some(location.natural),
                ));
            }
        }
        Ok((None, None, None))
    }

    pub fn sample(
        &self,
        model: &Model,
        result: &Arc<StepResult>,
        field: Field,
        pick: Pick,
        points: Vec<[f64; 3]>,
        samples: usize,
        values: &[f64],
    ) -> Result<Plot> {
        if values.len() != result.nodes().len() {
            return Err(error("the selected field does not have one value per result node"));
        }
        let (distance, xyz) = stations(&points, samples)?;
        let mut plot = Plot {
            step: result.step().into(),
            field,
            pick,
            mode: matches!(result.kind(), crate::results::Kind::Frequency(_)).then_some(pick.mode),
            points,
            distance,
            values: Vec::with_capacity(samples),
            elements: Vec::with_capacity(samples),
            natural: Vec::with_capacity(samples),
            xyz,
            result: Arc::downgrade(result),
        };
        let mut work = Work::new(MAX_WORK);
        let mut visits = 0;
        fieldcraft_model::progress::stage("Sampling path", Some(samples as f64), "stations");
        for (i, &point) in plot.xyz.iter().enumerate() {
            fieldcraft_model::progress::check()?;
            if i.is_multiple_of(16) {
                fieldcraft_model::progress::count(i as f64);
            }
            let (element, value, natural) = self.station(model, result.nodes(), values, point, &mut work, &mut visits)?;
            plot.elements.push(element);
            plot.values.push(value);
            plot.natural.push(natural);
        }
        fieldcraft_model::progress::check()?;
        fieldcraft_model::progress::count(samples as f64);
        Ok(plot)
    }
}

/// Equally spaced reference-coordinate arc-length stations. Zero-length segments are skipped;
/// an entirely coincident path is refused. Endpoint coordinates are preserved exactly.
pub fn stations(points: &[[f64; 3]], samples: usize) -> Result<(Vec<f64>, Vec<[f64; 3]>)> {
    if !(2..=MAX_POINTS).contains(&points.len()) || !(2..=MAX_SAMPLES).contains(&samples) || points.iter().flatten().any(|v| !v.is_finite()) {
        return Err(error(format!("use 2 to {MAX_POINTS} finite polyline points and 2 to {MAX_SAMPLES} samples")));
    }
    let mut length = Vec::with_capacity(points.len());
    length.push(0.0);
    let mut total = 0.0;
    for pair in points.windows(2) {
        if let [a, b] = pair {
            total += (b[0] - a[0]).hypot(b[1] - a[1]).hypot(b[2] - a[2]);
            length.push(total);
        }
    }
    if !(total > 0.0 && total.is_finite()) {
        return Err(error("path must have a finite, nonzero arc length"));
    }
    let mut distance = Vec::with_capacity(samples);
    let mut xyz = Vec::with_capacity(samples);
    let mut segment = 1;
    for i in 0..samples {
        let d = if i + 1 == samples { total } else { total * (i as f64 / (samples - 1) as f64) };
        while segment + 1 < length.len() && length.get(segment).is_some_and(|v| *v <= d) {
            segment += 1;
        }
        let a = points.get(segment - 1).ok_or_else(|| error("invalid path segment"))?;
        let b = points.get(segment).ok_or_else(|| error("invalid path segment"))?;
        let start = length.get(segment - 1).copied().ok_or_else(|| error("invalid path length"))?;
        let end = length.get(segment).copied().ok_or_else(|| error("invalid path length"))?;
        let t = if end > start { (d - start) / (end - start) } else { 0.0 };
        let point = if i == 0 {
            *points.first().ok_or_else(|| error("empty path"))?
        } else if i + 1 == samples {
            *points.last().ok_or_else(|| error("empty path"))?
        } else {
            [0, 1, 2].map(|k| a[k] * (1.0 - t) + b[k] * t)
        };
        if point.iter().any(|v| !v.is_finite()) {
            return Err(error("path station coordinates overflow"));
        }
        distance.push(d);
        xyz.push(point);
    }
    Ok((distance, xyz))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polyline_stations_follow_arc_length_and_preserve_endpoints() {
        let points = [[0.0; 3], [0.0; 3], [3.0, 0.0, 0.0], [3.0, 4.0, 0.0], [3.0, 4.0, 0.0]];
        let (distance, xyz) = stations(&points, 8).unwrap();
        assert_eq!(distance, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
        assert_eq!(xyz[0], points[0]);
        assert_eq!(xyz[3], points[2]);
        assert_eq!(xyz[7], points[4]);
        assert_eq!(xyz[5], [3.0, 2.0, 0.0]);
        assert!(stations(&[[0.0; 3]; 2], 5).is_err());
        assert!(stations(&points, MAX_SAMPLES + 1).is_err());
        assert!(stations(&vec![[0.0; 3]; MAX_POINTS + 1], 5).is_err());
        assert!(stations(&[[0.0; 3], [f64::NAN; 3]], 5).is_err());
        assert!(stations(&[[-f64::MAX; 3], [f64::MAX; 3]], 5).is_err());
    }

    #[test]
    fn deterministic_overlap_undefined_values_and_bounded_candidate_work() {
        let mut model = Model::default();
        for (i, p) in fieldcraft_fem::Basis::Hex8.nodes().iter().enumerate() {
            model.set_node(i as u32 + 1, *p).unwrap();
        }
        let nodes: Vec<_> = (1..=8).collect();
        model.set_element(17, ElementKind::Hex8, nodes.clone()).unwrap();
        model.set_element(3, ElementKind::Hex8, nodes.clone()).unwrap();
        let index = Index::build(&model, &nodes).unwrap();
        let found = index.station(&model, &nodes, &[2.0; 8], [0.0; 3], &mut Work::new(100), &mut 0).unwrap();
        assert_eq!(found, (Some(3), Some(2.0), Some([0.0; 3])));
        let unknown = index.station(&model, &nodes, &[f64::NAN; 8], [0.0; 3], &mut Work::new(100), &mut 0).unwrap();
        assert_eq!(unknown, (Some(3), None, Some([0.0; 3])));
        assert!(index.station(&model, &nodes, &[2.0; 8], [0.0; 3], &mut Work::new(0), &mut 0).is_err());
        let mut visits = MAX_VISITS;
        assert!(index.candidates([0.0; 3], &mut visits).is_err());
        for label in 1..=MAX_CANDIDATES as u32 + 1 {
            model.set_element(label, ElementKind::Hex8, nodes.clone()).unwrap();
        }
        let index = Index::build(&model, &nodes).unwrap();
        assert!(index.candidates([0.0; 3], &mut 0).unwrap_err().to_string().contains("overlapping"));
        // Unsupported dimensional families are a named refusal, not silent gaps.
        model.set_element(1, ElementKind::Tet4, vec![1, 2, 3, 5]).unwrap();
        model.set_element(2, ElementKind::Tri3, vec![1, 2, 3]).unwrap();
        assert!(Index::build(&model, &nodes).unwrap_err().to_string().contains("solid elements"));
        assert!(Index::build(&model, &[2, 1]).unwrap_err().to_string().contains("increasing"));
    }

    #[test]
    fn curved_display_face_pick_maps_natural_coordinates_not_the_drawing_chord() {
        let kind = ElementKind::Hex20;
        let basis = fieldcraft_fem::Basis::of(kind);
        let mut model = Model::default();
        for (i, [r, s, t]) in basis.nodes().iter().enumerate() {
            model.set_node(i as u32 + 1, [*r, s + r - (1.0 - r * r), *t]).unwrap();
        }
        model.set_element(1, kind, (1..=20).collect()).unwrap();
        let bary = [1.0 / 3.0; 3];
        let point = reference_surface_point(&model, 1, [1, 9, 12], bary).unwrap();
        let xi = [-2.0 / 3.0, -2.0 / 3.0, -1.0];
        let expected = [xi[0], xi[1] + xi[0] - (1.0 - xi[0] * xi[0]), xi[2]];
        assert!((0..3).all(|k| (point[k] - expected[k]).abs() < 1e-12));
        let chord = [1, 9, 12].map(|n| model.node(n).unwrap().xyz[1]).iter().sum::<f64>() / 3.0;
        assert!((point[1] - chord).abs() > 0.1);
        assert!(reference_surface_point(&model, 1, [1, 9, 999], bary).is_none());
        assert!(reference_surface_point(&model, 1, [1, 9, 12], [f64::NAN; 3]).is_none());
        assert!(reference_surface_point(&model, 1, [1, 9, 12], [2.0, -1.0, 0.0]).is_none());
    }
}
