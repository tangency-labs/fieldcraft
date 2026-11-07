//! Strong relaxation: an optimisation of the node positions of a whole region toward elements of
//! scaled Jacobian 1 (`mesh.node.relax` with `"mode": "strong"`).
//!
//! The unknowns are the positions of the region's corner nodes: interior nodes in space, nodes on
//! a geometry face in the face's parameters ([`Constraints::slide_face`]); nodes on geometry edges
//! and vertices are held, so feature lines stay as they are. The objective is
//!
//! ```text
//! F = Σₑ φ(qₑ) + ½ ln²(Vₑ/V₀ₑ)   for a valid element (q its scaled Jacobian, V its size),
//! φ(q) = (1/q − 1)²               phase 1 (smooth): every element toward 1, the worst most,
//! φ(q) = q⁻¹⁶ − 1                 phase 2 (max–min): a smooth stand-in for the smallest q,
//! B (2 − q)                        for an invalid element (B large),
//! ```
//!
//! (on a hexahedron, linear or quadratic, φ is the mean of φ over its nine Jacobian samples, smooth
//! where φ of the smallest has a kink),
//! φ continuing as its tangent line below q = 0.3, so it rises steeply but stays finite as an
//! element degenerates, while the invalid branch falls as an inverted element untangles and lies
//! above every valid element's term (a barrier: no step that lowers `F` inverts an element unless
//! it untangles worse ones; and passing through q → 0⁺ never costs more than staying inverted). The size term is there because the scaled Jacobian does not see an
//! element flattening along one direction: without it the optimiser squeezed elements against
//! the boundary to square up their neighbours (the Pi 5 base's largest edge ratio went from 5 to
//! 34). For the same reason an element longer than the region's largest edge ratio before costs
//! as much as an invalid one.
//!
//! `F` is minimised by block coordinate descent (nonlinear Gauss–Seidel): each sweep visits the
//! movable nodes in label order and lowers the part of `F` that depends on the node, first by
//! trying the centroid of its edge neighbours (Laplacian), then by a gradient step (central
//! differences, backtracking line search, the step length adapted per node). A move is kept when
//! it lowers the local objective, so individual elements may get worse on the way. L. Freitag,
//! P. Plassmann, *Local optimization-based simplicial mesh untangling and improvement*, IJNME 49,
//! 2000; P. Knupp, *Hexahedral and tetrahedral mesh untangling*, Eng. Comput. 17, 2001 (the
//! barrier objective on the Jacobian); J. Nocedal, S. Wright, *Numerical Optimization*, ch. 3
//! (backtracking).
//!
//! Phase 1 has up to half the sweeps, phase 2 (from phase 1's best state) the rest; a phase stops
//! when a sweep lowers its objective by less than a relative 1e-5, and after the first sweep only
//! nodes whose neighbourhood moved are visited. The sweeps keep the best state they pass through (highest smallest quality, then highest
//! mean). The result is that state if it has no invalid element and its smallest quality is not
//! below the starting one; otherwise nothing moves. Deterministic: no randomness, fixed order.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::Model;

use super::{Constraints, Elem, Freedom, P, ROUND_OFF, RegionQuality, State, add, dist, dot, measure, region_quality, scale, sub};
use crate::{Error, Result};

/// The smallest fraction of a step tried.
const MIN_STEP: f64 = 1e-6; // unit-factors: a fraction of a step

/// The relative change of the local objective over an edge length below which a difference
/// gradient's direction is not trusted (its round-off is about 1e-11 of that).
const TRUSTED_GRADIENT: f64 = 1e-3; // unit-factors: a relative change of the objective

/// Weight of an invalid element in the objective: above every valid element's term (the shape
/// term is at most [`shape`] at 0, a few 1e9, plus the size term), so no step that untangles an
/// element ever looks worse than staying inverted.
const INVALID: f64 = 1e11;

/// Below this quality the shape term continues as its tangent line, so it stays finite (and
/// below [`INVALID`]) as an element degenerates instead of growing without bound.
const KNEE: f64 = 0.3;

/// Weight of the size term.
const SIZE_WEIGHT: f64 = 0.5;

/// Exponent of the max–min phase.
const MIN_MAX_POWER: i32 = 16;

/// What a phase minimises.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum Phase {
    /// `(1/q − 1)²`: every element toward 1, the worst weighted most.
    #[default]
    Smooth,
    /// `q⁻¹⁶`: a smooth stand-in for the smallest quality, which it lifts.
    MinMax,
}

/// A value of the objective, `INVALID · barrier + valid`, kept in its two parts: summed as one
/// number, the barrier terms of invalid and stretched elements (1e11 each) rounded the valid
/// elements' terms away, so a move that left the barrier as it was compared round-off, and the
/// difference gradient of a node next to an invalid element pointed where round-off did (#612).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Objective {
    /// The invalid and stretched elements' terms, in units of [`INVALID`].
    barrier: f64,
    /// The valid elements' terms.
    valid: f64,
}

impl Objective {
    /// `self − other`, the parts subtracted first.
    fn minus(self, other: Objective) -> f64 {
        INVALID * (self.barrier - other.barrier) + (self.valid - other.valid)
    }

    fn plus(self, other: Objective) -> Objective {
        Objective { barrier: self.barrier + other.barrier, valid: self.valid + other.valid }
    }

    /// The size of the value, for its round-off.
    fn magnitude(self) -> f64 {
        INVALID * self.barrier.abs() + self.valid.abs()
    }

    fn is_zero(self) -> bool {
        self.barrier == 0.0 && self.valid == 0.0
    }
}

/// The objective's term for one element of quality `q` and shape term `shaped` whose size went
/// from `v0` to `v`.
fn phi(q: f64, shaped: f64, valid: bool, v: f64, v0: f64) -> Objective {
    if valid && q > 0.0 {
        let g = if v > 0.0 && v0 > 0.0 { libm::log(v / v0) } else { 0.0 };
        Objective { barrier: 0.0, valid: shaped + SIZE_WEIGHT * g * g }
    } else {
        Objective { barrier: 2.0 - q.min(1.0), valid: 0.0 }
    }
}

/// The shape term of a valid element of quality `q > 0`: `φ(q)` from [`KNEE`] up, its tangent
/// line below. Unbounded, `q⁻¹⁶` reached 1e96 as an element neared flat, so far above an
/// invalid element's term that the descent kept elements inverted rather than pass through
/// q → 0⁺ (#438: the distorted plate stayed tangled).
fn shape(phase: Phase, q: f64) -> f64 {
    // φ and its derivative at `x`.
    let at = |x: f64| match phase {
        Phase::Smooth => ((1.0 / x - 1.0).powi(2), -2.0 * (1.0 / x - 1.0) / (x * x)),
        Phase::MinMax => (x.powi(-MIN_MAX_POWER) - 1.0, -f64::from(MIN_MAX_POWER) * x.powi(-MIN_MAX_POWER - 1)),
    };
    if q >= KNEE {
        at(q).0
    } else {
        let (v, d) = at(KNEE);
        v + d * (q - KNEE)
    }
}

/// A linear hexahedron's scaled Jacobian at its corners and centre (the samples
/// [`crate::quality::element_quality`] takes), straight from the trilinear map: the sweeps'
/// inner loop, where the certified check would dominate. A state is certified before it is kept.
/// The Jacobian at the centre and the scaled Jacobian of every sample.
fn hex8_scaled_jacobian(x: &[P]) -> (f64, [f64; 9]) {
    let mut samples = [0.0; 9];
    let corners = &crate::quality::jacobian::HEX_CORNERS;
    let at = |k: usize| x.get(k).copied().unwrap_or([0.0; 3]);
    let mut centre = 0.0;
    for (sample, u) in samples.iter_mut().zip(corners.iter().map(|c| c.map(|v| v as f64)).chain(std::iter::once([0.5; 3]))) {
        let mut col = [[0.0; 3]; 3];
        for (k, bits) in corners.iter().enumerate() {
            for (d, c) in col.iter_mut().enumerate() {
                let w: f64 = (0..3)
                    .map(|e| {
                        let b = bits[e] as f64;
                        if e == d { 2.0 * b - 1.0 } else { b * u[e] + (1.0 - b) * (1.0 - u[e]) }
                    })
                    .product();
                *c = add(*c, scale(at(k), w));
            }
        }
        (centre, *sample) = jacobian_sample(col);
    }
    (centre, samples)
}

/// A quadratic (20-node serendipity) hexahedron's scaled Jacobian at the same samples as
/// [`hex8_scaled_jacobian`], from the shape functions' derivatives: the Jacobian at the centre and
/// the scaled Jacobian of every sample.
pub(super) fn hex20_scaled_jacobian(x: &[P]) -> (f64, [f64; 9]) {
    use crate::quality::jacobian::{HEX_EDGES, hex_corner};
    let at = |k: usize| x.get(k).copied().unwrap_or([0.0; 3]);
    let corners: [P; 8] = std::array::from_fn(hex_corner);
    let mut samples = [0.0; 9];
    let mut centre = 0.0;
    for (sample, u) in samples.iter_mut().zip(corners.iter().copied().chain(std::iter::once([0.0; 3]))) {
        // ∂x/∂uᵈ: the shape functions' derivatives at `u` times the nodes.
        let mut col = [[0.0; 3]; 3];
        for (k, c) in corners.iter().enumerate() {
            let f = [0, 1, 2].map(|e| 1.0 + u[e] * c[e]);
            for (d, column) in col.iter_mut().enumerate() {
                let (o1, o2) = ((d + 1) % 3, (d + 2) % 3);
                let w = 0.125 * c[d] * f[o1] * f[o2] * (2.0 * u[d] * c[d] + u[o1] * c[o1] + u[o2] * c[o2] - 1.0);
                *column = add(*column, scale(at(k), w));
            }
        }
        for (i, &(a, b)) in HEX_EDGES.iter().enumerate() {
            // A mid-side node's factor along each axis and its derivative: 1 − u² along its edge,
            // 1 + u·r across it.
            let (f, df): ([f64; 3], [f64; 3]) = {
                let pair = [0, 1, 2].map(|e| {
                    let r = 0.5 * (corners[a][e] + corners[b][e]);
                    if r == 0.0 { (1.0 - u[e] * u[e], -2.0 * u[e]) } else { (1.0 + u[e] * r, r) }
                });
                (pair.map(|p| p.0), pair.map(|p| p.1))
            };
            for (d, column) in col.iter_mut().enumerate() {
                *column = add(*column, scale(at(8 + i), 0.25 * df[d] * f[(d + 1) % 3] * f[(d + 2) % 3]));
            }
        }
        (centre, *sample) = jacobian_sample(col);
    }
    (centre, samples)
}

/// A Jacobian sample from its columns: the determinant and the scaled Jacobian (−1 when a
/// column is degenerate).
fn jacobian_sample(col: [P; 3]) -> (f64, f64) {
    let l = dot(col[0], col[0]).sqrt() * dot(col[1], col[1]).sqrt() * dot(col[2], col[2]).sqrt();
    let det = dot(col[0], super::cross(col[1], col[2]));
    (det, if l > 0.0 && l.is_finite() { det / l } else { -1.0 })
}

/// A hexahedron's Jacobian samples ([`hex8_scaled_jacobian`], [`hex20_scaled_jacobian`]); `None`
/// for another kind.
fn hex_scaled_jacobian(kind: fieldcraft_model::ElementKind, x: &[P]) -> Option<(f64, [f64; 9])> {
    match kind {
        fieldcraft_model::ElementKind::Hex8 => Some(hex8_scaled_jacobian(x)),
        fieldcraft_model::ElementKind::Hex20 => Some(hex20_scaled_jacobian(x)),
        _ => None,
    }
}

/// An element's size for the size term: a hexahedron's Jacobian at its centre (its volume, near
/// enough, up to a constant factor), else the size [`crate::quality::element_quality`] reports.
pub(super) fn size(kind: fieldcraft_model::ElementKind, x: &[P]) -> f64 {
    match hex_scaled_jacobian(kind, x) {
        Some((centre, _)) => centre,
        None => crate::quality::element_quality(kind, x).map_or(0.0, |q| q.size),
    }
}

/// An element's objective term, quality and validity; `certified` takes the certified check
/// for every kind.
fn term_of(phase: Phase, e: &Elem, s: &State<'_>, certified: bool) -> (Objective, f64, bool) {
    let x = s.coords(e);
    let fast = if certified { None } else { hex_scaled_jacobian(e.kind, &x) };
    let (q, valid, v, shaped) = if let Some((v, samples)) = fast {
        // The shape term of every sample, averaged: the smallest one's had a kink where the
        // worst sample changes, over which the difference gradient turned within the difference
        // step and a round-off change of the position turned the step with it (#612).
        let q = samples.iter().copied().fold(f64::INFINITY, f64::min);
        (q, q > 0.0, v, samples.iter().map(|&c| shape(phase, c)).sum::<f64>() / samples.len() as f64)
    } else {
        let (q, valid) = measure(e, &x);
        (q, valid, size(e.kind, &x), shape(phase, q))
    };
    // An element stretched past the cap is as bad as an invalid one.
    let ratio = super::edge_ratio(e.kind, &x);
    let stretch = if ratio > s.ratio_cap * (1.0 + ROUND_OFF) { ratio / s.ratio_cap } else { 0.0 };
    (phi(q, shaped, valid, v, e.size0).plus(Objective { barrier: stretch, valid: 0.0 }), q, valid)
}

/// The change of the objective from the terms `from` to the terms `to` of the same elements,
/// element by element: an element that did not change adds exactly nothing, where a difference
/// of sums would add the round-off of its term, which (an invalid element's 1e11, or a bad
/// element's q⁻¹⁶) can be far larger than the change of the elements that did (#612).
fn change(to: &[Objective], from: &[Objective]) -> f64 {
    to.iter().zip(from).map(|(a, b)| a.minus(*b)).sum()
}

impl State<'_> {
    /// The objective's terms of the elements at `n`.
    fn local_objective(&self, phase: Phase, n: u32) -> Vec<Objective> {
        let around = self.around.get(&n).map(Vec::as_slice).unwrap_or(&[]);
        around.iter().filter_map(|&i| self.elems.get(i)).map(|e| term_of(phase, e, self, false).0).collect()
    }

    /// `n` moved by `d` under `freedom`, kept when the local objective drops below `f0`; the new
    /// local objective when kept.
    fn try_lower(&mut self, phase: Phase, n: u32, freedom: Freedom, d: P, f0: &[Objective]) -> Option<Vec<Objective>> {
        let p = self.at(n);
        let q = self.moved(freedom, p, d)?;
        // A move of round-off is no move: whether it lowers the objective is round-off too.
        if dist(p, q) <= ROUND_OFF * self.local_length(n) {
            return None;
        }
        let undo = self.place(n, q);
        let f = self.local_objective(phase, n);
        if change(&f, f0) < 0.0 {
            Some(f)
        } else {
            self.restore(undo);
            None
        }
    }

    /// One optimisation of node `n`; `step` is its adaptive step (relative to the local edge
    /// length). Whether it moved.
    fn optimize_node(&mut self, phase: Phase, n: u32, step: &mut f64) -> bool {
        let freedom = self.freedom(n);
        if freedom == Freedom::Fixed {
            return false;
        }
        let mut f0 = self.local_objective(phase, n);
        if f0.iter().all(|f| f.is_zero()) {
            return false;
        }
        let mut moved = false;
        // Laplacian: the centroid of the edge neighbours, kept when it lowers the objective.
        let nb: Vec<u32> = self.neighbours.get(&n).map(|s| s.iter().copied().collect()).unwrap_or_default();
        if !nb.is_empty() {
            let c = scale(nb.iter().fold([0.0; 3], |s, &m| add(s, self.at(m))), 1.0 / nb.len() as f64);
            let d = sub(c, self.at(n));
            for k in [1.0, 0.5, 0.25] {
                if let Some(f) = self.try_lower(phase, n, freedom, scale(d, k), &f0) {
                    f0 = f;
                    moved = true;
                    break;
                }
            }
        }
        // Gradient step (central differences in space; on a face the slide maps the difference
        // onto the surface, so the gradient is the surface gradient).
        let h = self.local_length(n);
        if !(h > 0.0) {
            return moved;
        }
        let here = self.at(n);
        let eps = 1e-5 * h;
        let mut g = [0.0; 3];
        for k in 0..3 {
            let mut f = [Vec::new(), Vec::new()];
            for (j, sign) in [1.0, -1.0].into_iter().enumerate() {
                let mut d = [0.0; 3];
                d[k] = sign * eps;
                f[j] = match self.moved(freedom, here, d) {
                    Some(q) => {
                        let undo = self.place(n, q);
                        let v = self.local_objective(phase, n);
                        self.restore(undo);
                        v
                    }
                    None => f0.clone(),
                };
            }
            g[k] = change(&f[0], &f[1]) / (2.0 * eps);
        }
        let gn = dot(g, g).sqrt();
        if !(gn > 0.0 && gn.is_finite()) {
            return moved;
        }
        // A gradient that changes the objective by less than a relative `TRUSTED_GRADIENT` over an
        // edge length is mostly round-off in direction: the step shrinks with it.
        let magnitude: f64 = f0.iter().map(|o| o.magnitude()).sum();
        let trust = (gn * h / (TRUSTED_GRADIENT * magnitude)).min(1.0);
        let dir = scale(g, -trust / gn);
        let mut t = *step;
        for attempt in 0..8 {
            if let Some(_f) = self.try_lower(phase, n, freedom, scale(dir, t * h), &f0) {
                // Grow after a first-try success, keep the shortened step otherwise.
                *step = if attempt == 0 { (t * 2.0).min(0.5) } else { t };
                return true;
            }
            t *= 0.5;
        }
        *step = (t * 2.0).max(MIN_STEP);
        moved
    }
}

/// When [`morph`] keeps a cluster's morph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Guards {
    /// Kept whenever no element it touches inverts (more invalid elements than before).
    #[default]
    Valid,
    /// Kept (at full or a halved scale) only when, besides, the worst element it touches is no
    /// worse, the smooth objective there drops and no element gets longer than the longest was.
    Strict,
}

/// Default reach of [`morph`], in mean edge lengths of the bad elements.
pub const DEFAULT_REACH_LENGTHS: f64 = 5.0;

/// What [`morph`] does.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MorphOptions {
    /// Elements below this quality are the bad ones.
    pub below: f64,
    /// The radius round the bad nodes whose nodes follow them; `None`:
    /// [`DEFAULT_REACH_LENGTHS`] times the bad elements' mean edge length.
    pub reach: Option<f64>,
    /// Most rounds (bad nodes moved, the region morphed after them).
    pub rounds: usize,
    /// When a cluster's morph is kept.
    pub guards: Guards,
}

impl Default for MorphOptions {
    fn default() -> Self {
        MorphOptions { below: 0.7, reach: None, rounds: 20, guards: Guards::Valid }
    }
}

/// One round of [`morph`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Round {
    /// Elements below the threshold, and their movable corner nodes (the ones moved first).
    pub bad_elements: usize,
    pub bad_nodes: usize,
    /// Nodes within the reach of the clusters (the ones that follow), summed.
    pub region_nodes: usize,
    pub reach: f64,
    /// Clusters of bad nodes (within the reach of each other), and those whose morph was kept.
    pub clusters: usize,
    pub kept: usize,
    /// Over all the elements, after the round.
    pub min_quality: f64,
    pub mean_quality: f64,
}

/// What [`morph`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Morphed {
    /// The new position of every node that moved.
    pub positions: BTreeMap<u32, P>,
    /// Corner nodes that moved.
    pub moved: usize,
    pub max_displacement: f64,
    pub before: RegionQuality,
    pub after: RegionQuality,
    pub mean_before: f64,
    pub mean_after: f64,
    pub rounds: Vec<Round>,
}

/// A certified look at some elements ([`State::certified`]).
struct Checked {
    invalid: usize,
    worst: f64,
    objective: Objective,
    ratio: f64,
}

impl Checked {
    /// Whether `self` (after a move) is better than `was`: fewer invalid elements, or as many
    /// with the worst valid one no worse and the objective lower.
    fn improves(&self, was: &Checked) -> bool {
        self.invalid < was.invalid
            || (self.invalid == was.invalid && self.worst >= was.worst - ROUND_OFF && self.objective.minus(was.objective) < 0.0)
    }
}

/// Sweeps of the bad nodes' own optimisation (their neighbours held).
const TARGET_SWEEPS: usize = 5;
/// Most Gauss–Seidel sweeps of the harmonic solve, and its tolerance (of the largest prescribed
/// displacement).
const HARMONIC_SWEEPS: usize = 2000;
const HARMONIC_TOL: f64 = 1e-6; // unit-factors: relative to the largest prescribed displacement
/// Most rounds of [`morph`].
pub const MAX_ROUNDS: usize = 200;

/// Nodes within `r` of `centres`, from a uniform grid of cell `r` over `points` (a k-d tree's
/// fixed-radius query, deterministic: the result is a sorted set).
fn within(points: &BTreeMap<u32, P>, centres: &[P], r: f64) -> BTreeSet<u32> {
    let cell = |p: P| p.map(|x| (x / r).floor() as i64);
    let mut grid: BTreeMap<[i64; 3], Vec<(u32, P)>> = BTreeMap::new();
    for (&n, &p) in points {
        grid.entry(cell(p)).or_default().push((n, p));
    }
    let mut out = BTreeSet::new();
    for &c in centres {
        let [i, j, k] = cell(c);
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    for &(n, p) in grid.get(&[i + di, j + dj, k + dk]).map(Vec::as_slice).unwrap_or(&[]) {
                        if dist(p, c) <= r {
                            out.insert(n);
                        }
                    }
                }
            }
        }
    }
    out
}

/// `bad` split into clusters: nodes within `reach` of each other (transitively) are one.
fn clusters(s: &State<'_>, bad: &BTreeSet<u32>, reach: f64) -> Vec<BTreeSet<u32>> {
    let points: BTreeMap<u32, P> = bad.iter().map(|&n| (n, s.at(n))).collect();
    let mut left = bad.clone();
    let mut out = Vec::new();
    while let Some(first) = left.pop_first() {
        let mut cluster = BTreeSet::from([first]);
        let mut front = vec![first];
        while let Some(n) = front.pop() {
            for m in within(&points, &[s.at(n)], reach) {
                if left.remove(&m) {
                    cluster.insert(m);
                    front.push(m);
                }
            }
        }
        out.push(cluster);
    }
    out
}

impl State<'_> {
    /// Morph round one cluster of bad nodes: their own optimisation (the rest held) gives their
    /// displacements, the harmonic extension carries them to the nodes within `reach`, applied at
    /// full or a halved scale when that improves the elements it touches. The region's size and
    /// whether it was kept. With [`Guards::Valid`] the full morph is kept unless an element
    /// it touches inverts; with [`Guards::Strict`] only a scale that improves the elements is.
    fn morph_cluster(&mut self, bad: &BTreeSet<u32>, reach: f64, guards: Guards, movable: &dyn Fn(&State<'_>, u32) -> bool) -> Result<(usize, bool)> {
        let start: BTreeMap<u32, P> = bad.iter().map(|&n| (n, self.at(n))).collect();
        let mut undo = Vec::new();
        let mut steps: BTreeMap<u32, f64> = bad.iter().map(|&n| (n, 0.1)).collect();
        for _ in 0..TARGET_SWEEPS {
            for &n in bad {
                let before = self.at(n);
                let mut step = steps.get(&n).copied().unwrap_or(0.1);
                if self.optimize_node(Phase::MinMax, n, &mut step) {
                    undo.push((n, before));
                }
                steps.insert(n, step);
            }
        }
        let targets: BTreeMap<u32, P> = start.iter().map(|(&n, p0)| (n, sub(self.at(n), *p0))).collect();
        for (n, p) in undo.into_iter().rev() {
            self.place(n, p);
        }
        let points: BTreeMap<u32, P> = self.around.keys().map(|&n| (n, self.at(n))).collect();
        let centres: Vec<P> = start.values().copied().collect();
        let region = within(&points, &centres, reach);
        let free: Vec<u32> = region.iter().copied().filter(|&n| !bad.contains(&n) && movable(self, n)).collect();
        let u = self.harmonic(&free, &targets)?;
        let touched: BTreeSet<usize> = u.keys().flat_map(|n| self.around.get(n).into_iter().flatten().copied()).collect();
        let was = self.certified(&touched);
        for t in [1.0, 0.5, 0.25, 0.125] {
            let mut undo = Vec::new();
            for (&n, d) in &u {
                if dot(*d, *d) == 0.0 {
                    continue;
                }
                let p0 = self.at(n);
                if let Some(p) = self.moved(self.freedom(n), p0, scale(*d, t)) {
                    undo.extend(self.place(n, p));
                }
            }
            let now = self.certified(&touched);
            let keep = match guards {
                Guards::Valid => now.invalid <= was.invalid,
                Guards::Strict => now.improves(&was) && now.ratio <= self.ratio_cap.max(was.ratio) * (1.0 + ROUND_OFF),
            };
            if keep {
                return Ok((region.len(), true));
            }
            self.restore(undo);
        }
        Ok((region.len(), false))
    }

    /// Every element's fast quality, its smallest and mean.
    fn qualities(&self) -> (Vec<f64>, f64, f64) {
        let q: Vec<f64> = self.elems.iter().map(|e| term_of(Phase::Smooth, e, self, false).1).collect();
        if q.is_empty() {
            return (q, 0.0, 0.0);
        }
        let min = q.iter().copied().fold(f64::INFINITY, f64::min);
        let mean = q.iter().sum::<f64>() / q.len() as f64;
        (q, min, mean)
    }

    /// Over elements `elems`, certified: how many are invalid, the worst quality of the valid
    /// ones, the smooth objective's sum and the largest edge ratio.
    fn certified(&self, elems: &BTreeSet<usize>) -> Checked {
        let mut c = Checked { invalid: 0, worst: f64::INFINITY, objective: Objective::default(), ratio: 0.0 };
        for &i in elems {
            let Some(e) = self.elems.get(i) else { continue };
            let x = self.coords(e);
            let (q, valid) = measure(e, &x);
            if valid {
                c.worst = c.worst.min(q);
            } else {
                c.invalid += 1;
            }
            c.objective = c.objective.plus(phi(q, shape(Phase::Smooth, q), valid, size(e.kind, &x), e.size0));
            c.ratio = c.ratio.max(super::edge_ratio(e.kind, &x));
        }
        c
    }

    /// The harmonic extension over `free` of the displacements `fixed` (zero at every other
    /// node): Gauss–Seidel on the graph Laplacian of the element edges.
    pub(super) fn harmonic(&self, free: &[u32], fixed: &BTreeMap<u32, P>) -> Result<BTreeMap<u32, P>> {
        let mut u: BTreeMap<u32, P> = fixed.clone();
        for &n in free {
            u.insert(n, [0.0; 3]);
        }
        let scale_of = fixed.values().map(|d| dot(*d, *d).sqrt()).fold(0.0, f64::max);
        for sweep in 0..HARMONIC_SWEEPS {
            if sweep % 64 == 0 {
                fieldcraft_model::progress::check()?;
            }
            let mut change = 0.0f64;
            for &n in free {
                let Some(nb) = self.neighbours.get(&n) else { continue };
                if nb.is_empty() {
                    continue;
                }
                let sum = nb.iter().fold([0.0; 3], |s, m| add(s, u.get(m).copied().unwrap_or([0.0; 3])));
                let v = scale(sum, 1.0 / nb.len() as f64);
                let old = u.insert(n, v).unwrap_or([0.0; 3]);
                change = change.max(dist(old, v));
            }
            if change <= HARMONIC_TOL * scale_of {
                break;
            }
        }
        Ok(u)
    }
}

/// Morph the mesh round its bad elements: each round finds the elements of quality below
/// `options.below` and, cluster by cluster (bad nodes within the reach of each other), moves
/// their movable corner nodes to better places (a few sweeps of their own
/// optimisation, the rest of the mesh held), and carries that displacement to the nodes within
/// the reach by harmonic interpolation (the graph Laplacian solved with the bad nodes'
/// displacements and zero beyond the reach and at held nodes). Face nodes are put back on their
/// face, edge and vertex nodes do not move, nodes outside the reach stay exactly as they are. A
/// cluster's morph is kept as long as no element it touches inverts ([`Guards::Valid`]), or, with
/// [`Guards::Strict`], at full or a halved scale only if besides no element gets longer than the
/// longest one was, the worst element it touches is no worse and the smooth objective there drops.
/// The result is the round with the highest smallest quality (never below the start). Only nodes
/// of `nodes` move. The model is not changed: the result has the new positions.
pub fn morph(model: &Model, nodes: &BTreeSet<u32>, options: MorphOptions, geometry: &dyn Constraints) -> Result<Morphed> {
    if !(options.below > 0.0 && options.below <= 1.0) {
        return Err(Error::Invalid(format!("below must be in (0, 1], got {}", options.below)));
    }
    if options.reach.is_some_and(|r| !(r > 0.0 && r.is_finite())) {
        return Err(Error::Invalid("the reach must be positive".into()));
    }
    if options.rounds == 0 || options.rounds > MAX_ROUNDS {
        return Err(Error::Invalid(format!("rounds must be 1 to {MAX_ROUNDS}, got {}", options.rounds)));
    }
    if nodes.len() > super::MAX_NODES {
        return Err(Error::Limit(format!("more than {} nodes to relax", super::MAX_NODES)));
    }
    let mut state = State::new(model, nodes, geometry, true)?;
    let original = state.pos.clone();
    let before = region_quality(&state);
    state.ratio_cap = before.max_edge_ratio.max(1.0);
    let (_, _, mean_before) = state.qualities();
    let movable = |s: &State<'_>, n: u32| nodes.contains(&n) && s.around.contains_key(&n) && s.freedom(n) != Freedom::Fixed;
    let mut rounds = Vec::new();
    // The best state of the run: the highest smallest quality, then the fewest elements below
    // `below`, then the highest mean; the start until a round beats it.
    let score = |s: &State<'_>| {
        let (q, min, mean) = s.qualities();
        (min, q.iter().filter(|v| **v < options.below).count(), mean)
    };
    // Qualities within round-off of each other are the same (#612).
    let better = |a: (f64, usize, f64), b: (f64, usize, f64)| {
        a.0 > b.0 + ROUND_OFF || (a.0 >= b.0 - ROUND_OFF && (a.1 < b.1 || (a.1 == b.1 && a.2 > b.2 + ROUND_OFF)))
    };
    let mut best = (score(&state), original.clone());
    for round in 0..options.rounds {
        fieldcraft_model::progress::stage(&format!("morph round {}", round + 1), Some(options.rounds as f64), "rounds");
        fieldcraft_model::progress::count(round as f64);
        fieldcraft_model::progress::check()?;
        let (q, _, _) = state.qualities();
        let bad_elems: Vec<usize> = (0..q.len()).filter(|&i| q.get(i).is_some_and(|v| *v < options.below)).collect();
        let mut bad: BTreeSet<u32> = BTreeSet::new();
        let mut edge_sum = (0.0, 0usize);
        for &i in &bad_elems {
            let Some(e) = state.elems.get(i) else { continue };
            bad.extend(e.nodes.iter().take(e.kind.corner_count()).copied().filter(|&n| movable(&state, n)));
            let x = state.coords(e);
            for edge in e.kind.edges() {
                if let (Some(a), Some(b)) = (edge.first().and_then(|&k| x.get(k)), edge.get(1).and_then(|&k| x.get(k))) {
                    edge_sum = (edge_sum.0 + dist(*a, *b), edge_sum.1 + 1);
                }
            }
        }
        if bad.is_empty() {
            break;
        }
        let reach = options.reach.unwrap_or(if edge_sum.1 > 0 { DEFAULT_REACH_LENGTHS * edge_sum.0 / edge_sum.1 as f64 } else { 0.0 });
        if !(reach > 0.0) {
            break;
        }
        // Bad nodes within the reach of each other are one cluster, morphed and checked together.
        let clusters = clusters(&state, &bad, reach);
        let mut round_out = Round { bad_elements: bad_elems.len(), bad_nodes: bad.len(), reach, clusters: clusters.len(), ..Round::default() };
        for cluster in &clusters {
            fieldcraft_model::progress::check()?;
            let (region, kept) = state.morph_cluster(cluster, reach, options.guards, &|s: &State<'_>, n| movable(s, n))?;
            round_out.region_nodes += region;
            round_out.kept += usize::from(kept);
        }
        let now = score(&state);
        (round_out.min_quality, round_out.mean_quality) = (now.0, now.2);
        if better(now, best.0) {
            best = (now, state.pos.clone());
        }
        let any = round_out.kept > 0;
        rounds.push(round_out);
        if !any {
            break;
        }
    }
    state.pos = best.1;
    let mut after = region_quality(&state);
    let (_, _, mut mean_after) = state.qualities();
    // Never worse than the start.
    if after.invalid > before.invalid || after.min_quality < before.min_quality - ROUND_OFF {
        state.pos = original.clone();
        after = before.clone();
        mean_after = mean_before;
    }
    let mut out = Morphed { before, after, mean_before, mean_after, rounds, ..Default::default() };
    for (n, p) in &state.pos {
        let Some(o) = original.get(n) else { continue };
        let d = dist(*o, *p);
        if d > 0.0 {
            out.max_displacement = out.max_displacement.max(d);
            out.moved += usize::from(nodes.contains(n));
            out.positions.insert(*n, *p);
        }
    }
    Ok(out)
}
