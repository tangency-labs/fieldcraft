//! Extruding a quadrilateral section up to a target surface.
//!
//! [`extrude_to_surface`] stacks layers of hexahedra on a quad mesh of a section (a face of a
//! body, a partition face or a 2D mesh) and ends them on one or more B-rep faces, the *target*,
//! instead of on a copy of the section: a boss meeting a curved wall, an air column up to a
//! domed lid, a fillet region capped by a curved face. The region need not be sweepable
//! between two parallel faces.
//!
//! **Paths.** Every source node has its own path: a straight line along the given direction, or
//! along the section's normal at the node. The path's end is where the line meets the target,
//! found by Newton's method on the tangent plane of the target's closest point (each step
//! intersects the line with the plane through the closest point, normal to the gap; it is exact
//! after one step on a plane and converges quadratically on smooth faces). The layers divide each
//! path by the same parameters, uniform or graded geometrically by `bias` (the last layer's
//! thickness over the first's).
//!
//! **Side faces.** When the section's sides are given, the boundary nodes of every layer are
//! projected onto them (onto the curve between two sides at a corner), and the boundary's last
//! layer is found by alternating projections onto the sides and the target (von Neumann's
//! alternating projections converge to the intersection of surfaces meeting transversally), so
//! the layers follow side faces that are not straight along the path. Interior nodes follow
//! their straight paths plus the boundary's displacement from its straight paths interpolated by
//! inverse-distance weights (the boundary-residual correction of Knupp, *Next-generation sweep
//! tool: a method for generating all-hex meshes on two-and-one-half dimensional geometries*,
//! Proc. 7th International Meshing Roundtable, 1998), and their last layer is placed back on the
//! target along the path direction.
//!
//! The source nodes are kept as given (the first nodes of the result, in order), so a section
//! taken from a neighbouring mesh conforms with it. Faces are projected onto through
//! [`crate::project::Geometry`], the trimmed B-rep projection of projection blocking: a point
//! beyond a face's boundary goes to the face, never to its surface's extension, so a path that
//! misses a target face does not converge on it. The result passes the certified validity
//! gate (every element's Jacobian determinant proven positive) and the tiling check (the
//! elements' volume equals the volume their boundary encloses).

use std::collections::BTreeMap;

use fieldcraft_model::ElementKind;

use super::face_flux;
use super::geom::{add, dist, scale, unit};
use super::sweep::shepard_weights;
use crate::project::Geometry;
use crate::quality::{GateOptions, GateReport, Metric, element_quality, validity_gate_checks};
use crate::vec3::{P, cross, dot, norm, sub};
use crate::{Error, Result};

/// The most layers.
pub const MAX_LAYERS: usize = 10_000;
/// Newton steps per path end, and alternating-projection rounds per node.
const MAX_NEWTON: usize = 60;
const MAX_ALTERNATING: usize = 100;
/// The most interior-node × boundary-node × layer products the residual correction may spend.
const MAX_CORRECTION_WORK: u64 = 2_000_000_000;
/// The largest layer bias (last over first thickness) and its inverse.
pub const MAX_BIAS: f64 = 1000.0;

/// The path of each source node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExtrudePath {
    /// Every node along this direction (need not be unit length).
    Direction(P),
    /// Every node along the section's normal there, turned toward the target.
    Normal,
}

/// How many layers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExtrudeLayers {
    Count(usize),
    /// Layers of about this thickness (along the mean path length).
    Size(f64),
}

/// What [`extrude_to_surface`] extrudes and where to.
#[derive(Debug, Clone)]
pub struct ExtrudeSpec<'a> {
    /// The section's nodes and quadrilaterals (indices into `nodes`).
    pub nodes: &'a [P],
    pub quads: &'a [[u32; 4]],
    pub path: ExtrudePath,
    /// The target faces and the side faces (face ids of the geometry).
    pub targets: &'a [u32],
    pub sides: &'a [u32],
    pub layers: ExtrudeLayers,
    /// The last layer's thickness over the first's (1: uniform).
    pub bias: f64,
    /// How far a node may be from a face and still be on it (the geometry's accuracy).
    pub tolerance: f64,
    /// The most hexahedra ([`Error::Limit`] beyond).
    pub max_elements: usize,
}

/// What a boundary face of the extruded mesh lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtrudeRole {
    /// The section.
    Source,
    /// A target face (its id).
    Target(u32),
    /// The side face it lies on, when sides were given.
    Side(Option<u32>),
}

/// A boundary face of the extruded mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtrudeBoundary {
    /// Element index (0-based).
    pub element: u32,
    /// Face number (1-based, as in [`ElementKind::faces`]).
    pub face: u8,
    pub role: ExtrudeRole,
}

/// The extruded hex8 mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Extruded {
    /// The source nodes first (unchanged, in order), then each layer's nodes.
    pub nodes: Vec<P>,
    pub hexes: Vec<[u32; 8]>,
    pub boundary: Vec<ExtrudeBoundary>,
    pub layers: usize,
    pub min_scaled_jacobian: f64,
    pub mean_scaled_jacobian: f64,
    /// The largest distance of a last-layer node from the target.
    pub max_target_gap: f64,
    pub gate: GateReport,
}

/// Extrude the section `spec` up to its target surface. See the module docs.
pub fn extrude_to_surface(spec: &ExtrudeSpec, geometry: &dyn Geometry) -> Result<Extruded> {
    let n_src = spec.nodes.len();
    check(spec)?;
    let tol = spec.tolerance;
    let ex = Extrusion { g: geometry, targets: spec.targets, tol };
    let boundary = section_boundary(spec.quads)?;

    // Path directions.
    let dirs: Vec<P> = match spec.path {
        ExtrudePath::Direction(d) => {
            let d = unit(d).ok_or_else(|| Error::Invalid("the direction has no length".into()))?;
            vec![d; n_src]
        }
        ExtrudePath::Normal => ex.normals(spec)?,
    };

    // Each side node's side faces.
    let mut side_faces: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    if !spec.sides.is_empty() {
        for &v in &boundary.nodes {
            let p = at(spec.nodes, v)?;
            let on: Vec<u32> = spec.sides.iter().copied().filter(|&f| geometry.project(f, p).is_some_and(|q| dist(p, q) <= 10.0 * tol)).collect();
            if on.is_empty() {
                return Err(Error::invalid(crate::message![
                    "the section's boundary node at ",
                    crate::message::point(p),
                    " is on none of the side faces"
                ]));
            }
            side_faces.insert(v, on);
        }
    }

    // Straight ends, and the boundary's ends on the sides.
    let mut straight = Vec::with_capacity(n_src);
    for (i, (&s, &d)) in spec.nodes.iter().zip(&dirs).enumerate() {
        let (e, _) = ex.hit(s, d, true).ok_or_else(|| {
            Error::failed(crate::message![
                format!("the path of section node {i} from "),
                crate::message::point(s),
                " along ",
                format!("({:.4}, {:.4}, {:.4})", d[0], d[1], d[2]),
                " does not reach the target surface"
            ])
        })?;
        straight.push(e);
    }
    let mut ends: BTreeMap<u32, P> = BTreeMap::new();
    for (&v, faces) in &side_faces {
        let e = at(&straight, v)?;
        ends.insert(
            v,
            ex.corner(faces, e)
                .ok_or_else(|| Error::failed(crate::message!["the side faces and the target surface do not meet near ", crate::message::point(e)]))?,
        );
    }

    // Layer parameters.
    let mean_length = spec.nodes.iter().zip(&straight).map(|(s, e)| dist(*s, *e)).sum::<f64>() / n_src as f64;
    let n = match spec.layers {
        ExtrudeLayers::Count(n) => n,
        ExtrudeLayers::Size(h) => ((mean_length / h).round() as usize).max(1),
    };
    if n == 0 || n > MAX_LAYERS {
        return Err(Error::Invalid(format!("{n} layers; an extrusion has 1 to {MAX_LAYERS}")));
    }
    let count = spec.quads.len().saturating_mul(n);
    if count > spec.max_elements {
        return Err(Error::Limit(format!("the extrusion would have {count} hexahedra; the limit is {}", spec.max_elements)));
    }
    let ts = parameters(n, spec.bias);

    // Interior correction weights, in the plane across the mean direction.
    let correct = !side_faces.is_empty();
    let ring: Vec<u32> = side_faces.keys().copied().collect();
    let interior: Vec<u32> = (0..n_src as u32).filter(|v| !side_faces.contains_key(v)).collect();
    let work = (interior.len() as u64).saturating_mul(ring.len() as u64).saturating_mul(n as u64);
    if correct && work > MAX_CORRECTION_WORK {
        return Err(Error::Limit(format!("placing the interior nodes would take {work} steps; the limit is {MAX_CORRECTION_WORK}")));
    }
    let weights: Vec<Vec<(usize, f64)>> = if correct {
        let mean_dir = unit(dirs.iter().fold([0.0; 3], |a, d| add(a, *d))).unwrap_or([0.0, 0.0, 1.0]);
        let (u, w) = across(mean_dir);
        let flat = |p: P| [dot(p, u), dot(p, w)];
        let ring2: Vec<[f64; 2]> = ring.iter().map(|&v| at(spec.nodes, v).map(flat)).collect::<Result<_>>()?;
        interior.iter().map(|&v| at(spec.nodes, v).map(|p| shepard_weights(flat(p), &ring2))).collect::<Result<_>>()?
    } else {
        Vec::new()
    };

    let mut nodes: Vec<P> = Vec::with_capacity(n_src * (n + 1));
    nodes.extend_from_slice(spec.nodes);
    let lerp = |s: P, e: P, t: f64| add(s, scale(sub(e, s), t));
    for (k, &t) in ts.iter().enumerate().skip(1) {
        let last = k == n;
        let mut layer: Vec<P> = spec.nodes.iter().zip(&straight).map(|(s, e)| lerp(*s, *e, t)).collect();
        // Boundary: onto the sides; their displacement from the straight paths.
        let mut residual = Vec::with_capacity(ring.len());
        for &v in &ring {
            let i = v as usize;
            let straight_at = layer.get(i).copied().unwrap_or([0.0; 3]);
            let x = if last {
                ends.get(&v).copied().unwrap_or(straight_at)
            } else {
                let faces = side_faces.get(&v).map_or(&[][..], Vec::as_slice);
                let s = at(spec.nodes, v)?;
                let target = lerp(s, ends.get(&v).copied().unwrap_or(straight_at), t);
                ex.onto(faces, target).ok_or_else(|| {
                    Error::failed(crate::message!["cannot place a layer node on the side faces near ", crate::message::point(target)])
                })?
            };
            residual.push(sub(x, straight_at));
            if let Some(slot) = layer.get_mut(i) {
                *slot = x;
            }
        }
        // Interior: the straight path plus the interpolated displacement; the last layer back on
        // the target along the path.
        for (j, &v) in interior.iter().enumerate() {
            let i = v as usize;
            let Some(&y0) = layer.get(i) else { continue };
            let correction = weights
                .get(j)
                .map_or([0.0; 3], |w| w.iter().fold([0.0; 3], |acc, &(r, wr)| add(acc, scale(residual.get(r).copied().unwrap_or([0.0; 3]), wr))));
            let mut y = add(y0, correction);
            if last && correct {
                let d = dirs.get(i).copied().unwrap_or([0.0, 0.0, 1.0]);
                y = ex.hit(y, d, false).map(|(p, _)| p).ok_or_else(|| {
                    Error::failed(crate::message!["the path near ", crate::message::point(y), " does not reach the target surface"])
                })?;
            }
            if let Some(slot) = layer.get_mut(i) {
                *slot = y;
            }
        }
        nodes.extend(layer);
    }
    fieldcraft_model::check_coordinates(nodes.as_flattened())?;

    // Hexahedra, each quad turned so its normal runs along the path.
    let mut hexes = Vec::with_capacity(count);
    let mut quads = Vec::with_capacity(spec.quads.len());
    for q in spec.quads {
        let [a, b, c, d] = q.map(|v| at(spec.nodes, v).unwrap_or([0.0; 3]));
        let normal = cross(sub(c, a), sub(d, b));
        let along = q.iter().fold([0.0; 3], |acc, &v| {
            let s = at(spec.nodes, v).unwrap_or([0.0; 3]);
            add(acc, sub(nodes.get(n_src + v as usize).copied().unwrap_or(s), s))
        });
        quads.push(if dot(normal, along) < 0.0 { [q[0], q[3], q[2], q[1]] } else { *q });
    }
    let stride = n_src as u32;
    for q in &quads {
        for k in 0..n as u32 {
            let lo = q.map(|v| v + k * stride);
            let hi = q.map(|v| v + (k + 1) * stride);
            hexes.push([lo[0], lo[1], lo[2], lo[3], hi[0], hi[1], hi[2], hi[3]]);
        }
    }

    // Boundary faces: the section, the last layer and the sides.
    let centre = |e: usize, face: usize| -> P {
        let h = hexes.get(e).copied().unwrap_or([0; 8]);
        let local = ElementKind::Hex8.faces().get(face).copied().unwrap_or(&[]);
        let pts: Vec<P> = local.iter().filter_map(|&i| h.get(i)).filter_map(|&v| nodes.get(v as usize).copied()).collect();
        scale(pts.iter().fold([0.0; 3], |a, p| add(a, *p)), 1.0 / pts.len().max(1) as f64)
    };
    let nearest = |faces: &[u32], p: P| -> Option<u32> {
        faces
            .iter()
            .filter_map(|&f| geometry.project(f, p).map(|q| (f, dist(p, q))))
            .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))
            .map(|(f, _)| f)
    };
    let mut bfaces = Vec::new();
    for (qi, q) in quads.iter().enumerate() {
        let first = qi * n;
        let last = first + n - 1;
        bfaces.push(ExtrudeBoundary { element: first as u32, face: 1, role: ExtrudeRole::Source });
        let f = nearest(spec.targets, centre(last, 1)).ok_or_else(|| Error::Failed("a last-layer face is on no target face".into()))?;
        bfaces.push(ExtrudeBoundary { element: last as u32, face: 2, role: ExtrudeRole::Target(f) });
        for m in 0..4 {
            let (a, b) = (q[m], q[(m + 1) % 4]);
            if !boundary.edges.contains(&(a.min(b), a.max(b))) {
                continue;
            }
            for e in first..=last {
                let side = if spec.sides.is_empty() { None } else { nearest(spec.sides, centre(e, 2 + m)) };
                bfaces.push(ExtrudeBoundary { element: e as u32, face: 3 + m as u8, role: ExtrudeRole::Side(side) });
            }
        }
    }

    // Quality and the certified gate.
    let kind = ElementKind::Hex8;
    let coords = |e: usize| -> Result<Vec<P>> { hexes.get(e).map_or(Ok(Vec::new()), |h| h.iter().map(|&v| at(&nodes, v)).collect()) };
    let mut enclosed = 0.0;
    for b in &bfaces {
        let x = coords(b.element as usize)?;
        let local = kind.faces().get(b.face as usize - 1).copied().unwrap_or(&[]);
        let face: Vec<P> = local.iter().filter_map(|&i| x.get(i).copied()).collect();
        enclosed += face_flux(&face);
    }
    let (mut min_sj, mut sum_sj, mut worst) = (f64::INFINITY, 0.0, 0usize);
    let mut checks = Vec::with_capacity(hexes.len());
    for e in 0..hexes.len() {
        let q = element_quality(kind, &coords(e)?)?;
        let sj = q.get(Metric::ScaledJacobian).unwrap_or(-1.0);
        if sj < min_sj {
            (min_sj, worst) = (sj, e);
        }
        sum_sj += sj;
        checks.push((e as u32 + 1, q.jacobian.ok_or_else(|| Error::Failed("an element was not certified".into()))?));
    }
    let gate = validity_gate_checks(checks.iter().map(|(l, c)| (*l, c)), &GateOptions { expected_volume: Some(enclosed), volume_tolerance: 1e-9 })?
        .into_result()
        .map_err(|e| {
            Error::failed(e.message().then(crate::message![
                format!(" (worst element {} near ", worst + 1),
                crate::message::point(centre(worst, 0)),
                ": a path crosses another or a side; try other layers, a smaller bias or the side faces)"
            ]))
        })?;
    let max_target_gap = nodes
        .iter()
        .skip(n * n_src)
        .map(|&p| spec.targets.iter().filter_map(|&f| geometry.project(f, p)).map(|q| dist(p, q)).fold(f64::INFINITY, f64::min))
        .fold(0.0, f64::max);
    Ok(Extruded {
        nodes,
        hexes,
        boundary: bfaces,
        layers: n,
        min_scaled_jacobian: min_sj,
        mean_scaled_jacobian: sum_sj / count.max(1) as f64,
        max_target_gap,
        gate,
    })
}

fn check(spec: &ExtrudeSpec) -> Result<()> {
    if spec.nodes.is_empty() || spec.quads.is_empty() {
        return Err(Error::Invalid("the section has no quadrilaterals".into()));
    }
    fieldcraft_model::check_coordinates(spec.nodes.as_flattened())?;
    if let Some(v) = spec.quads.iter().flatten().find(|&&v| v as usize >= spec.nodes.len()) {
        return Err(Error::Invalid(format!("section node {v} is out of range")));
    }
    if spec.quads.iter().any(|q| (0..4).any(|i| (i + 1..4).any(|j| q[i] == q[j]))) {
        return Err(Error::Invalid("a section quadrilateral repeats a node".into()));
    }
    if spec.targets.is_empty() {
        return Err(Error::Invalid("no target surface".into()));
    }
    if !(spec.tolerance > 0.0 && spec.tolerance.is_finite()) {
        return Err(Error::Invalid("the tolerance must be a positive number".into()));
    }
    if !(spec.bias >= 1.0 / MAX_BIAS && spec.bias <= MAX_BIAS) {
        return Err(Error::Invalid(format!("the bias must be from {} to {MAX_BIAS}", 1.0 / MAX_BIAS)));
    }
    match spec.layers {
        ExtrudeLayers::Count(n) if n == 0 || n > MAX_LAYERS => Err(Error::Invalid(format!("an extrusion has 1 to {MAX_LAYERS} layers"))),
        ExtrudeLayers::Size(h) if !(h > 0.0 && h.is_finite()) => Err(Error::Invalid("the layer size must be a positive number".into())),
        _ => Ok(()),
    }
}

/// The section's boundary: edges on one quadrilateral, and their nodes.
struct SectionBoundary {
    edges: Vec<(u32, u32)>,
    nodes: Vec<u32>,
}

fn section_boundary(quads: &[[u32; 4]]) -> Result<SectionBoundary> {
    let mut uses: BTreeMap<(u32, u32), u32> = BTreeMap::new();
    for q in quads {
        for m in 0..4 {
            let (a, b) = (q[m], q[(m + 1) % 4]);
            *uses.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    if let Some((&(a, b), _)) = uses.iter().find(|(_, n)| **n > 2) {
        return Err(Error::Invalid(format!("the section is not a surface: the edge between nodes {a} and {b} is on more than two quadrilaterals")));
    }
    let edges: Vec<(u32, u32)> = uses.into_iter().filter(|(_, n)| *n == 1).map(|(e, _)| e).collect();
    let mut nodes: Vec<u32> = edges.iter().flat_map(|&(a, b)| [a, b]).collect();
    nodes.sort_unstable();
    nodes.dedup();
    Ok(SectionBoundary { edges, nodes })
}

/// The layer parameters `0 = t_0 < … < t_n = 1`, thicknesses in geometric progression with the
/// last over the first equal to `bias`.
pub(crate) fn parameters(n: usize, bias: f64) -> Vec<f64> {
    let r = if n > 1 { libm::pow(bias, 1.0 / (n - 1) as f64) } else { 1.0 };
    let mut out = Vec::with_capacity(n + 1);
    let (mut sum, mut h) = (0.0, 1.0);
    out.push(0.0);
    for _ in 0..n {
        sum += h;
        out.push(sum);
        h *= r;
    }
    let total = sum;
    for t in &mut out {
        *t /= total;
    }
    if let Some(t) = out.last_mut() {
        *t = 1.0;
    }
    out
}

/// The geometry queries of one extrusion.
struct Extrusion<'a> {
    g: &'a dyn Geometry,
    targets: &'a [u32],
    tol: f64,
}

impl Extrusion<'_> {
    /// Where the line `o + t d` meets the target: the hit with the smallest `t` ahead of `o`
    /// (`ahead`), or nearest to `o` either way. Returns the point and its target face.
    fn hit(&self, o: P, d: P, ahead: bool) -> Option<(P, u32)> {
        let mut best: Option<(P, u32, f64)> = None;
        for &f in self.targets {
            let Some((p, t)) = self.hit_face(f, o, d) else { continue };
            if ahead && t <= self.tol {
                continue;
            }
            let key = if ahead { t } else { t.abs() };
            if best.is_none_or(|b| key < b.2) {
                best = Some((p, f, key));
            }
        }
        best.map(|(p, f, _)| (p, f))
    }

    /// Newton's method for the line `o + t d` on face `f`: each step meets the tangent plane at
    /// the closest point. `None` when it does not converge onto the face (the line misses it,
    /// or grazes it).
    fn hit_face(&self, f: u32, o: P, d: P) -> Option<(P, f64)> {
        let q0 = self.g.project(f, o)?;
        let mut t = dot(sub(q0, o), d);
        let eps = 1e-10 * (1.0 + norm(o).max(dist(o, q0)));
        for _ in 0..MAX_NEWTON {
            let x = add(o, scale(d, t));
            let q = self.g.project(f, x)?;
            let gap = sub(x, q);
            let g2 = dot(gap, gap);
            if g2.sqrt() <= eps {
                return Some((q, t));
            }
            let den = dot(d, gap);
            if den.abs() <= 1e-6 * g2.sqrt() {
                return None;
            }
            t -= g2 / den;
            if !t.is_finite() {
                return None;
            }
        }
        // A free-form face's facets can make the last steps cycle within its accuracy.
        let x = add(o, scale(d, t));
        let q = self.g.project(f, x)?;
        (dist(x, q) <= self.tol).then_some((q, t))
    }

    /// `p` on all of `faces` (one face: its closest point; several: alternating projections).
    fn onto(&self, faces: &[u32], p: P) -> Option<P> {
        let mut x = p;
        if let [f] = faces {
            return self.g.project(*f, x);
        }
        for _ in 0..MAX_ALTERNATING {
            let before = x;
            for &f in faces {
                x = self.g.project(f, x)?;
            }
            if dist(before, x) <= 1e-12 * (1.0 + norm(x)) {
                break;
            }
        }
        faces.iter().all(|&f| self.g.project(f, x).is_some_and(|q| dist(q, x) <= self.tol)).then_some(x)
    }

    /// The point near `p` on the side `faces` and the target (alternating projections, ending on
    /// the target exactly).
    fn corner(&self, faces: &[u32], p: P) -> Option<P> {
        let mut x = p;
        let mut face = self.nearest_target(x)?;
        for _ in 0..MAX_ALTERNATING {
            let before = x;
            for &f in faces {
                x = self.g.project(f, x)?;
            }
            face = self.nearest_target(x)?;
            x = self.g.project(face, x)?;
            if dist(before, x) <= 1e-12 * (1.0 + norm(x)) {
                break;
            }
        }
        let on_sides = faces.iter().all(|&f| self.g.project(f, x).is_some_and(|q| dist(q, x) <= self.tol));
        let on_target = self.g.project(face, x).is_some_and(|q| dist(q, x) <= self.tol);
        (on_sides && on_target).then_some(x)
    }

    fn nearest_target(&self, p: P) -> Option<u32> {
        self.targets
            .iter()
            .filter_map(|&f| self.g.project(f, p).map(|q| (f, dist(p, q))))
            .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))
            .map(|(f, _)| f)
    }

    /// Unit normals of the section at its nodes, turned toward the target.
    fn normals(&self, spec: &ExtrudeSpec) -> Result<Vec<P>> {
        let mut sum = vec![[0.0; 3]; spec.nodes.len()];
        for q in spec.quads {
            let [a, b, c, d] = q.map(|v| at(spec.nodes, v).unwrap_or([0.0; 3]));
            let n = cross(sub(c, a), sub(d, b));
            for &v in q {
                if let Some(s) = sum.get_mut(v as usize) {
                    *s = add(*s, n);
                }
            }
        }
        let mut out = Vec::with_capacity(sum.len());
        for (i, s) in sum.iter().enumerate() {
            out.push(unit(*s).ok_or_else(|| Error::Invalid(format!("the section has no normal at node {i} (it folds there)")))?);
        }
        // Toward the target: the side of the section its nearest point is on.
        let centre = scale(spec.nodes.iter().fold([0.0; 3], |a, p| add(a, *p)), 1.0 / spec.nodes.len() as f64);
        let mean = out.iter().fold([0.0; 3], |a, n| add(a, *n));
        let f = self.nearest_target(centre).ok_or_else(|| Error::Failed("no target face to extrude toward".into()))?;
        let q = self.g.project(f, centre).ok_or_else(|| Error::Failed("no target face to extrude toward".into()))?;
        if dot(sub(q, centre), mean) < 0.0 {
            for n in &mut out {
                *n = scale(*n, -1.0);
            }
        }
        Ok(out)
    }
}

/// Two unit vectors across `d`.
fn across(d: P) -> (P, P) {
    let t = if d[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let u = unit(cross(d, t)).unwrap_or([1.0, 0.0, 0.0]);
    (u, cross(d, u))
}

fn at(nodes: &[P], v: u32) -> Result<P> {
    nodes.get(v as usize).copied().ok_or_else(|| Error::Failed(format!("node {v} out of range")))
}

#[cfg(test)]
mod tests;
