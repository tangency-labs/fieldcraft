//! Node relaxation: moving a few nodes of a finished mesh to fix its worst elements, the nodes on
//! the geometry staying on it (`mesh.node.relax`).
//!
//! Each node to move carries a weight in `(0, 1]` (the strength, times a falloff that blends the
//! edited zone into the mesh around it). A sweep visits the nodes in label order and tries, for
//! each, two moves:
//!
//! 1. **Smart Laplacian.** The weighted step toward the centroid of the corner nodes it shares an
//!    element edge with, kept only when the worst element around the node does not get worse and
//!    no valid element becomes invalid; otherwise the step is halved (three times at most).
//!    L. Freitag, *On combining Laplacian and optimization-based mesh smoothing techniques*, AMD
//!    220, 1997 (the quality-guarded Laplacian); S. Canann, J. Tristano, M. Staten, *An approach
//!    to combined Laplacian and optimization-based smoothing for triangular, quadrilateral, and
//!    quad-dominant meshes*, IMR 1998 (smart Laplacian: a move is accepted only if the local
//!    quality does not drop).
//! 2. **Optimization.** A step up the gradient of the worst element's quality around the node
//!    (central differences), with a backtracking line search, kept only when the worst quality
//!    rises: the max–min formulation of L. Freitag, C. Ollivier-Gooch, *Tetrahedral mesh
//!    improvement using swapping and smoothing*, IJNME 40, 1997, and P. Knupp, *Achieving finite
//!    element mesh quality via optimization of the Jacobian matrix norm and associated
//!    quantities*, IJNME 48, 2000, restricted to the active (worst) element.
//!
//! An element's quality is its scaled Jacobian (a triangle's: its smallest angle over 60°), from
//! [`crate::quality::element_quality`]; a solid is valid when its certified Jacobian check is,
//! a shell when its normal keeps its side. So no move inverts an element, and none makes the
//! worst element around a node worse, and none makes an element's edge ratio exceed the
//! region's largest before the relaxation.
//!
//! **Constraints** ([`Freedom`]): an interior node moves freely; a node on a geometry face moves
//! in the face's parameters and is evaluated back on it ([`Constraints::slide_face`]); a node on a
//! geometry edge slides along it ([`Constraints::slide_curve`]); a fixed node (a geometry vertex,
//! or a skin node whose geometry is unknown) does not move. Mid-side nodes of quadratic elements
//! follow their edge's corners, keeping their offset from the edge's midpoint (the curvature),
//! and are put back on their face or edge.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::{ElementKind, Model, RelaxationCell};

use crate::quality::{self, Metric, Validity};
use crate::{Error, Result};

/// Qualities this close are equal.
const SAME_QUALITY: f64 = 1e-9; // unit-factors: a difference of qualities
/// The difference step in an element's parameters.
const STEP: f64 = 1e-3; // unit-factors: a step in an element's parameters

type P = [f64; 3];

/// Most sweeps over the nodes.
pub const MAX_ITERATIONS: usize = 100;
/// Most nodes one relaxation moves.
pub const MAX_NODES: usize = 200_000;
/// Elements listed as the worst of the region.
pub const WORST_LISTED: usize = 5;

/// How a node may move.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Freedom {
    /// Anywhere (an interior node).
    Free,
    /// Not at all.
    Fixed,
    /// On the face with this id (the caller's).
    Face(u32),
    /// Along the curve with this id (the caller's).
    Curve(u32),
}

/// Where nodes may go: the caller's faces and curves.
pub trait Constraints {
    /// How `node` may move; `None` when the caller does not know the node's geometry: it then
    /// moves freely inside the mesh and is held on its skin.
    fn freedom(&self, node: u32) -> Option<Freedom>;
    /// `p` (on face `face`) moved by `d` in the face's parameters; `None` refuses the move.
    fn slide_face(&self, face: u32, p: P, d: P) -> Option<P>;
    /// `p` (on curve `curve`) moved by `d` along the curve; `None` refuses the move.
    fn slide_curve(&self, curve: u32, p: P, d: P) -> Option<P>;
}

/// No geometry: interior nodes move, every node on the skin is held.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoGeometry;

impl Constraints for NoGeometry {
    fn freedom(&self, _: u32) -> Option<Freedom> {
        None
    }
    fn slide_face(&self, _: u32, _: P, _: P) -> Option<P> {
        None
    }
    fn slide_curve(&self, _: u32, _: P, _: P) -> Option<P> {
        None
    }
}

/// How the weight falls from the centre of a sphere of influence to its radius.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Falloff {
    /// `(1 − t²)²` of the relative distance `t`: 1 at the centre, 0 with a zero slope at the
    /// radius, so the edited zone blends into the mesh around it.
    #[default]
    Smooth,
    /// `1 − t`.
    Linear,
    /// The full strength everywhere inside.
    None,
}

impl Falloff {
    pub fn from_name(name: &str) -> Option<Falloff> {
        match name {
            "smooth" => Some(Falloff::Smooth),
            "linear" => Some(Falloff::Linear),
            "none" => Some(Falloff::None),
            _ => None,
        }
    }

    /// The weight at relative distance `t` in `[0, 1]`.
    pub fn weight(self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        match self {
            Falloff::Smooth => (1.0 - t * t).powi(2),
            Falloff::Linear => 1.0 - t,
            Falloff::None => 1.0,
        }
    }
}

/// The quality of the elements around the moved nodes.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RegionQuality {
    pub elements: usize,
    /// Smallest element quality (scaled Jacobian; a triangle's smallest angle over 60°).
    pub min_quality: f64,
    /// Largest ratio of an element's longest to its shortest edge.
    pub max_edge_ratio: f64,
    /// Elements whose Jacobian (or a shell's normal) is not valid.
    pub invalid: usize,
    /// The worst elements and their quality, worst first.
    pub worst: Vec<(u32, f64)>,
}

/// What [`relax`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Relaxed {
    /// The new position of every node that moved (corners and the mid-side nodes following them).
    pub positions: BTreeMap<u32, P>,
    /// Corner nodes that moved.
    pub moved: usize,
    /// Largest distance a node moved.
    pub max_displacement: f64,
    pub before: RegionQuality,
    pub after: RegionQuality,
    /// Sweeps made (fewer than asked when a sweep moved nothing).
    pub iterations: usize,
}

/// The corner nodes of measured elements (solids and shells) inside the sphere, weighted by
/// `strength` × the falloff; only nodes of `only` when given.
pub fn sphere_weights(
    model: &Model,
    center: P,
    radius: f64,
    strength: f64,
    falloff: Falloff,
    only: Option<&BTreeSet<u32>>,
) -> Result<BTreeMap<u32, f64>> {
    if !(radius > 0.0 && radius.is_finite()) || !center.iter().all(|c| c.is_finite()) {
        return Err(Error::Invalid("the sphere needs a finite centre and a positive radius".into()));
    }
    check_strength(strength)?;
    let mut out = BTreeMap::new();
    for n in corner_nodes(model) {
        if only.is_some_and(|o| !o.contains(&n)) {
            continue;
        }
        let Some(node) = model.node(n) else { continue };
        let t = dist(node.xyz, center) / radius;
        if t < 1.0 {
            let w = strength * falloff.weight(t);
            if w > 0.0 {
                out.insert(n, w);
            }
        }
        if out.len() > MAX_NODES {
            return Err(Error::Limit(format!("more than {MAX_NODES} nodes inside the sphere; give a smaller radius")));
        }
    }
    Ok(out)
}

/// The corner nodes among `nodes` (nodes of measured elements that are a corner of one), each
/// weighted `strength`.
pub fn uniform_weights(model: &Model, nodes: &BTreeSet<u32>, strength: f64) -> Result<BTreeMap<u32, f64>> {
    check_strength(strength)?;
    let corners = corner_nodes(model);
    let out: BTreeMap<u32, f64> = nodes.iter().filter(|n| corners.contains(n)).map(|&n| (n, strength)).collect();
    if out.len() > MAX_NODES {
        return Err(Error::Limit(format!("more than {MAX_NODES} nodes to relax")));
    }
    Ok(out)
}

fn check_strength(strength: f64) -> Result<()> {
    if strength > 0.0 && strength <= 1.0 { Ok(()) } else { Err(Error::Invalid(format!("strength must be in (0, 1], got {strength}"))) }
}

/// Corner nodes of the solid and shell elements.
fn corner_nodes(model: &Model) -> BTreeSet<u32> {
    let mut out = BTreeSet::new();
    for e in model.elements().values().filter(|e| quality::is_measured(e.kind)) {
        out.extend(e.nodes.iter().take(e.kind.corner_count()).copied());
    }
    out
}

/// Nodes on the skin of the mesh: on a solid face no other solid shares, or on a shell.
fn skin_nodes(model: &Model, near: &BTreeSet<u32>) -> BTreeSet<u32> {
    let mut faces: BTreeMap<Vec<u32>, (usize, Vec<u32>)> = BTreeMap::new();
    let mut skin = BTreeSet::new();
    for e in model.elements().values() {
        if !e.nodes.iter().any(|n| near.contains(n)) {
            continue;
        }
        match e.kind.dimension() {
            3 => {
                for f in e.kind.faces() {
                    let nodes: Vec<u32> = f.iter().filter_map(|&i| e.nodes.get(i).copied()).collect();
                    let mut key: Vec<u32> = nodes.iter().take(corner_len(f.len())).copied().collect();
                    key.sort_unstable();
                    faces.entry(key).or_insert((0, nodes)).0 += 1;
                }
            }
            2 => skin.extend(e.nodes.iter().copied()),
            _ => {}
        }
    }
    for (count, nodes) in faces.into_values() {
        if count == 1 {
            skin.extend(nodes);
        }
    }
    skin
}

/// Corners of a face with `n` nodes (3, 4 linear; 6, 8 quadratic).
fn corner_len(n: usize) -> usize {
    match n {
        6 => 3,
        8 => 4,
        n => n,
    }
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: P, b: P) -> P {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: P, s: f64) -> P {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: P, b: P) -> P {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dist(a: P, b: P) -> f64 {
    dot(sub(a, b), sub(a, b)).sqrt()
}

/// One element of the region.
struct Elem {
    label: u32,
    kind: ElementKind,
    nodes: Vec<u32>,
    /// A shell's normal as it was: its side, which a move must keep.
    normal: Option<P>,
    /// Its size as it was ([`strong::size`]): the strong mode keeps sizes near it (only the strong
    /// mode measures it, a certified check per element; 0 otherwise).
    size0: f64,
}

/// The shell normal from its corners (`None` for solids).
fn shell_normal(kind: ElementKind, x: &[P]) -> Option<P> {
    if kind.dimension() != 2 {
        return None;
    }
    let c = |i: usize| x.get(i).copied().unwrap_or([0.0; 3]);
    Some(match kind.corner_count() {
        3 => cross(sub(c(1), c(0)), sub(c(2), c(0))),
        _ => cross(sub(c(2), c(0)), sub(c(3), c(1))),
    })
}

/// An element's quality and whether it is valid.
fn measure(e: &Elem, x: &[P]) -> (f64, bool) {
    let Ok(q) = quality::element_quality(e.kind, x) else { return (-1.0, false) };
    let score = q.get(Metric::ScaledJacobian).or_else(|| q.get(Metric::MinAngle).map(|a| a / 60.0)).unwrap_or(-1.0);
    let valid = match (q.validity(), e.normal) {
        (Some(v), _) => v == Validity::Valid,
        (None, Some(n0)) => shell_normal(e.kind, x).is_some_and(|n| dot(n, n0) > 0.0) && q.size > 0.0,
        (None, None) => q.size > 0.0,
    };
    (score, valid)
}

/// Longest over shortest corner edge.
fn edge_ratio(kind: ElementKind, x: &[P]) -> f64 {
    let (mut lo, mut hi) = (f64::INFINITY, 0.0f64);
    for e in kind.edges() {
        if let (Some(&a), Some(&b)) = (e.first().and_then(|&i| x.get(i)), e.get(1).and_then(|&i| x.get(i))) {
            let l = dist(a, b);
            lo = lo.min(l);
            hi = hi.max(l);
        }
    }
    if lo > 0.0 { hi / lo } else { f64::INFINITY }
}

/// The working state: positions, the region's elements and who depends on whom.
struct State<'a> {
    pos: BTreeMap<u32, P>,
    elems: Vec<Elem>,
    /// Corner node → its elements (indices into `elems`).
    around: BTreeMap<u32, Vec<usize>>,
    /// Corner node → corner nodes it shares an element edge with.
    neighbours: BTreeMap<u32, BTreeSet<u32>>,
    /// Mid-side node → its edge's corners and its offset from their midpoint.
    mids: BTreeMap<u32, (u32, u32, P)>,
    /// Corner node → the mid-side nodes on its edges.
    mids_of: BTreeMap<u32, Vec<u32>>,
    skin: BTreeSet<u32>,
    geometry: &'a dyn Constraints,
    /// Hold nodes on geometry edges (the strong mode keeps feature lines as they are).
    pin_curves: bool,
    /// The largest edge ratio the strong mode may make.
    ratio_cap: f64,
}

impl State<'_> {
    fn at(&self, n: u32) -> P {
        self.pos.get(&n).copied().unwrap_or([0.0; 3])
    }

    fn coords(&self, e: &Elem) -> Vec<P> {
        e.nodes.iter().map(|&n| self.at(n)).collect()
    }

    /// How `n` may move: a skin node the caller left free is held.
    fn freedom(&self, n: u32) -> Freedom {
        match self.geometry.freedom(n) {
            Some(Freedom::Curve(_)) if self.pin_curves => Freedom::Fixed,
            Some(f) => f,
            None if self.skin.contains(&n) => Freedom::Fixed,
            None => Freedom::Free,
        }
    }

    /// `p` moved by `d` under `freedom`.
    fn moved(&self, freedom: Freedom, p: P, d: P) -> Option<P> {
        let q = match freedom {
            Freedom::Fixed => return None,
            Freedom::Free => add(p, d),
            Freedom::Face(f) => self.geometry.slide_face(f, p, d)?,
            Freedom::Curve(c) => self.geometry.slide_curve(c, p, d)?,
        };
        q.iter().all(|x| x.is_finite() && x.abs() <= fieldcraft_model::MAX_COORD).then_some(q)
    }

    /// The worst quality around `n` and the validity of each element there.
    fn local(&self, n: u32) -> (f64, Vec<bool>) {
        let mut worst = f64::INFINITY;
        let mut valid = Vec::new();
        for &i in self.around.get(&n).map(Vec::as_slice).unwrap_or(&[]) {
            let Some(e) = self.elems.get(i) else { continue };
            let (q, v) = measure(e, &self.coords(e));
            worst = worst.min(q);
            valid.push(v);
        }
        (worst, valid)
    }

    /// The largest edge ratio of the elements around `n`.
    fn local_ratio(&self, n: u32) -> f64 {
        let around = self.around.get(&n).map(Vec::as_slice).unwrap_or(&[]);
        around.iter().filter_map(|&i| self.elems.get(i)).map(|e| edge_ratio(e.kind, &self.coords(e))).fold(0.0, f64::max)
    }

    /// The worst element around `n`: its index and quality.
    fn worst_around(&self, n: u32) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;
        for &i in self.around.get(&n)? {
            let Some(e) = self.elems.get(i) else { continue };
            let (q, _) = measure(e, &self.coords(e));
            if best.is_none_or(|(_, b)| q < b) {
                best = Some((i, q));
            }
        }
        best
    }

    /// Put `n` at `p` with its mid-side nodes following; returns what to restore.
    fn place(&mut self, n: u32, p: P) -> Vec<(u32, P)> {
        let mut undo = vec![(n, self.at(n))];
        self.pos.insert(n, p);
        for m in self.mids_of.get(&n).cloned().unwrap_or_default() {
            let Some(&(a, b, bulge)) = self.mids.get(&m) else { continue };
            let old = self.at(m);
            let guess = add(scale(add(self.at(a), self.at(b)), 0.5), bulge);
            let q = match self.freedom(m) {
                Freedom::Free => Some(guess),
                Freedom::Fixed => None,
                f => self.moved(f, old, sub(guess, old)),
            };
            if let Some(q) = q {
                undo.push((m, old));
                self.pos.insert(m, q);
            }
        }
        undo
    }

    fn restore(&mut self, undo: Vec<(u32, P)>) {
        for (n, p) in undo.into_iter().rev() {
            self.pos.insert(n, p);
        }
    }

    /// Try moving `n` by `d` (under its freedom): kept when no valid element around it becomes
    /// invalid and its worst quality does not drop (rises by more than `gain`).
    fn attempt(&mut self, n: u32, freedom: Freedom, d: P, before: &(f64, Vec<bool>), gain: f64) -> bool {
        let p = self.at(n);
        let Some(q) = self.moved(freedom, p, d) else { return false };
        if dist(p, q) == 0.0 {
            return false;
        }
        let undo = self.place(n, q);
        let (worst, valid) = self.local(n);
        let ok = worst >= before.0 + gain && before.1.iter().zip(&valid).all(|(was, is)| !was || *is) && self.local_ratio(n) <= self.ratio_cap;
        if !ok {
            self.restore(undo);
        }
        ok
    }

    /// The mean length of the element edges at `n`.
    fn local_length(&self, n: u32) -> f64 {
        let p = self.at(n);
        let nb = self.neighbours.get(&n);
        let (sum, count) = nb.into_iter().flatten().fold((0.0, 0usize), |(s, c), &m| (s + dist(p, self.at(m)), c + 1));
        if count > 0 { sum / count as f64 } else { 0.0 }
    }

    /// One relaxation of `n` with weight `w`; whether it moved.
    fn relax_node(&mut self, n: u32, w: f64) -> bool {
        let freedom = self.freedom(n);
        if freedom == Freedom::Fixed {
            return false;
        }
        let p = self.at(n);
        let mut moved = false;
        // 1. Smart Laplacian.
        let nb: Vec<u32> = self.neighbours.get(&n).map(|s| s.iter().copied().collect()).unwrap_or_default();
        if !nb.is_empty() {
            let centroid = scale(nb.iter().fold([0.0; 3], |s, &m| add(s, self.at(m))), 1.0 / nb.len() as f64);
            let before = self.local(n);
            let mut d = scale(sub(centroid, p), w);
            for _ in 0..4 {
                if self.attempt(n, freedom, d, &before, 0.0) {
                    moved = true;
                    break;
                }
                d = scale(d, 0.5);
            }
        }
        // 2. Optimization of the worst element around the node.
        let h = self.local_length(n);
        if !(h > 0.0) {
            return moved;
        }
        let Some((worst, q0)) = self.worst_around(n) else { return moved };
        let Some(e) = self.elems.get(worst) else { return moved };
        let slot = e.nodes.iter().position(|&m| m == n);
        let Some(slot) = slot else { return moved };
        let mut x = self.coords(e);
        let here = self.at(n);
        let eps = 1e-4 * h;
        let mut grad = [0.0; 3];
        for k in 0..3 {
            let mut f = [0.0; 2];
            for (j, sign) in [1.0, -1.0].into_iter().enumerate() {
                if let Some(c) = x.get_mut(slot) {
                    c[k] = here[k] + sign * eps;
                }
                f[j] = measure(e, &x).0;
            }
            if let Some(c) = x.get_mut(slot) {
                c[k] = here[k];
            }
            grad[k] = (f[0] - f[1]) / (2.0 * eps);
        }
        let g = dot(grad, grad).sqrt();
        if !(g > 0.0 && g.is_finite()) || q0 >= 1.0 - SAME_QUALITY {
            return moved;
        }
        let before = self.local(n);
        let mut step = 0.25 * w * h;
        for _ in 0..5 {
            if self.attempt(n, freedom, scale(grad, step / g), &before, SAME_QUALITY) {
                return true;
            }
            step *= 0.5;
        }
        moved
    }
}

impl<'a> State<'a> {
    /// The elements with a corner among `movers`, their nodes' positions and adjacency; with
    /// `pin_curves`, nodes on geometry edges are held.
    fn new(model: &Model, movers: &BTreeSet<u32>, geometry: &'a dyn Constraints, pin_curves: bool) -> Result<State<'a>> {
        let mut state = State {
            pos: BTreeMap::new(),
            elems: Vec::new(),
            around: BTreeMap::new(),
            neighbours: BTreeMap::new(),
            mids: BTreeMap::new(),
            mids_of: BTreeMap::new(),
            skin: skin_nodes(model, movers),
            geometry,
            pin_curves,
            ratio_cap: f64::INFINITY,
        };
        for (&label, e) in model.elements() {
            if !quality::is_measured(e.kind) || e.nodes.len() != e.kind.node_count() {
                continue;
            }
            let corners = e.kind.corner_count();
            if !e.nodes.iter().take(corners).any(|n| movers.contains(n)) {
                continue;
            }
            let mut x = Vec::with_capacity(e.nodes.len());
            for &n in &e.nodes {
                let p = model.node(n).map(|n| n.xyz).ok_or_else(|| Error::Invalid(format!("element {label} refers to missing node {n}")))?;
                state.pos.insert(n, p);
                x.push(p);
            }
            let index = state.elems.len();
            for &n in e.nodes.iter().take(corners) {
                if movers.contains(&n) {
                    state.around.entry(n).or_default().push(index);
                }
            }
            for edge in e.kind.edges() {
                let (Some(&a), Some(&b)) = (edge.first().and_then(|&i| e.nodes.get(i)), edge.get(1).and_then(|&i| e.nodes.get(i))) else { continue };
                for (u, v) in [(a, b), (b, a)] {
                    if movers.contains(&u) {
                        state.neighbours.entry(u).or_default().insert(v);
                    }
                }
                if let Some(&m) = edge.get(2).and_then(|&i| e.nodes.get(i)) {
                    let (pa, pb, pm) = (x.get(edge[0]).copied(), x.get(edge[1]).copied(), x.get(edge[2]).copied());
                    if let (Some(pa), Some(pb), Some(pm)) = (pa, pb, pm) {
                        state.mids.entry(m).or_insert((a, b, sub(pm, scale(add(pa, pb), 0.5))));
                        for c in [a, b] {
                            let list = state.mids_of.entry(c).or_default();
                            if !list.contains(&m) {
                                list.push(m);
                            }
                        }
                    }
                }
            }
            state.elems.push(Elem {
                label,
                kind: e.kind,
                normal: shell_normal(e.kind, &x),
                size0: if pin_curves { strong::size(e.kind, &x) } else { 0.0 },
                nodes: e.nodes.clone(),
            });
        }
        Ok(state)
    }
}

/// Quality over the region.
fn region_quality(s: &State<'_>) -> RegionQuality {
    let mut r = RegionQuality { elements: s.elems.len(), min_quality: f64::INFINITY, max_edge_ratio: 0.0, invalid: 0, worst: Vec::new() };
    let mut all = Vec::with_capacity(s.elems.len());
    for e in &s.elems {
        let x = s.coords(e);
        let (q, valid) = measure(e, &x);
        r.min_quality = r.min_quality.min(q);
        r.max_edge_ratio = r.max_edge_ratio.max(edge_ratio(e.kind, &x));
        r.invalid += usize::from(!valid);
        all.push((e.label, q));
    }
    all.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    all.truncate(WORST_LISTED);
    r.worst = all;
    if r.elements == 0 {
        r.min_quality = 0.0;
    }
    r
}

/// Relax the nodes of `weights` (node → weight in `(0, 1]`) for `iterations` sweeps, nodes
/// moving as `geometry` allows. The model is not changed: the result has the new positions.
pub fn relax(model: &Model, weights: &BTreeMap<u32, f64>, iterations: usize, geometry: &dyn Constraints) -> Result<Relaxed> {
    if iterations == 0 || iterations > MAX_ITERATIONS {
        return Err(Error::Invalid(format!("iterations must be 1 to {MAX_ITERATIONS}, got {iterations}")));
    }
    if weights.len() > MAX_NODES {
        return Err(Error::Limit(format!("more than {MAX_NODES} nodes to relax")));
    }
    if let Some((n, w)) = weights.iter().find(|(_, w)| !(**w > 0.0 && **w <= 1.0)) {
        return Err(Error::Invalid(format!("node {n} has weight {w}; weights are in (0, 1]")));
    }
    let movers: BTreeSet<u32> = weights.keys().copied().collect();
    let mut state = State::new(model, &movers, geometry, false)?;
    let original = state.pos.clone();
    let before = region_quality(&state);
    // A move may not stretch an element past the region's longest one (the scaled Jacobian alone
    // lets an element flatten: the crank's worst edge ratio went from 10.5 to 2759).
    state.ratio_cap = before.max_edge_ratio.max(1.0);
    let mut done = 0;
    for _ in 0..iterations {
        fieldcraft_model::progress::stage("relax sweep", Some(iterations as f64), "sweeps");
        fieldcraft_model::progress::count(done as f64);
        fieldcraft_model::progress::check()?;
        done += 1;
        let mut any = false;
        for (&n, &w) in weights {
            if state.around.contains_key(&n) {
                any |= state.relax_node(n, w);
            }
        }
        if !any {
            break;
        }
    }
    let after = region_quality(&state);
    let mut out = Relaxed { before, after, iterations: done, ..Default::default() };
    for (n, p) in &state.pos {
        let Some(o) = original.get(n) else { continue };
        let d = dist(*o, *p);
        if d > 0.0 {
            out.max_displacement = out.max_displacement.max(d);
            out.moved += usize::from(movers.contains(n));
            out.positions.insert(*n, *p);
        }
    }
    Ok(out)
}

// ---- Relaxation fields ----

/// The linear kind of a quadratic element's corners.
fn linear(kind: ElementKind) -> Option<ElementKind> {
    Some(match kind {
        ElementKind::Tri3 | ElementKind::Tri6 => ElementKind::Tri3,
        ElementKind::Quad4 | ElementKind::Quad8 => ElementKind::Quad4,
        ElementKind::Tet4 | ElementKind::Tet10 => ElementKind::Tet4,
        ElementKind::Wedge6 => ElementKind::Wedge6,
        ElementKind::Hex8 | ElementKind::Hex20 => ElementKind::Hex8,
        _ => return None,
    })
}

/// The relaxation as a field ([`fieldcraft_model::Relaxation`]'s cells): every element with a
/// moved node, its corners as they were in `model` and as `relaxed` left them.
pub fn field_cells(model: &Model, relaxed: &Relaxed) -> Vec<RelaxationCell> {
    let mut cells = Vec::new();
    for e in model.elements().values() {
        let Some(kind) = linear(e.kind) else { continue };
        if !e.nodes.iter().any(|n| relaxed.positions.contains_key(n)) {
            continue;
        }
        let corners: Vec<u32> = e.nodes.iter().take(kind.node_count()).copied().collect();
        let at: Vec<P> = corners.iter().filter_map(|n| model.node(*n).map(|x| x.xyz)).collect();
        if at.len() != corners.len() {
            continue;
        }
        let moved = corners.iter().zip(&at).map(|(n, a)| relaxed.positions.get(n).copied().unwrap_or(*a)).collect();
        cells.push(RelaxationCell { kind, at, moved });
    }
    cells
}

/// Linear shape functions of a cell at parametric `xi`, and the domain's dimension.
fn shape(kind: ElementKind, xi: P) -> (Vec<f64>, usize) {
    let [u, v, w] = xi;
    match kind {
        ElementKind::Tri3 => (vec![1.0 - u - v, u, v], 2),
        ElementKind::Quad4 => (vec![(1.0 - u) * (1.0 - v), u * (1.0 - v), u * v, (1.0 - u) * v], 2),
        ElementKind::Tet4 => (vec![1.0 - u - v - w, u, v, w], 3),
        ElementKind::Wedge6 => {
            let t = [1.0 - u - v, u, v];
            (t.iter().map(|a| a * (1.0 - w)).chain(t.iter().map(|a| a * w)).collect(), 3)
        }
        _ => {
            let b = [(1.0 - u) * (1.0 - v), u * (1.0 - v), u * v, (1.0 - u) * v];
            (b.iter().map(|a| a * (1.0 - w)).chain(b.iter().map(|a| a * w)).collect(), 3)
        }
    }
}

/// How far `xi` lies outside the cell's parametric domain (0 inside).
fn outside(kind: ElementKind, xi: P) -> f64 {
    let [u, v, w] = xi;
    let box01 = |x: f64| (-x).max(x - 1.0).max(0.0);
    match kind {
        ElementKind::Tri3 => (-u).max(-v).max(u + v - 1.0).max(0.0),
        ElementKind::Tet4 => (-u).max(-v).max(-w).max(u + v + w - 1.0).max(0.0),
        ElementKind::Quad4 => box01(u).max(box01(v)),
        ElementKind::Wedge6 => (-u).max(-v).max(u + v - 1.0).max(0.0).max(box01(w)),
        _ => box01(u).max(box01(v)).max(box01(w)),
    }
}

/// The nearest point of the parametric domain to `xi`.
fn clamp_domain(kind: ElementKind, xi: P) -> P {
    let c = |x: f64| x.clamp(0.0, 1.0);
    let simplex = |u: f64, v: f64, w: f64| {
        let (u, v, w) = (u.max(0.0), v.max(0.0), w.max(0.0));
        let s = u + v + w;
        if s > 1.0 { (u / s, v / s, w / s) } else { (u, v, w) }
    };
    match kind {
        ElementKind::Tri3 => {
            let (u, v, _) = simplex(xi[0], xi[1], 0.0);
            [u, v, 0.0]
        }
        ElementKind::Tet4 => {
            let (u, v, w) = simplex(xi[0], xi[1], xi[2]);
            [u, v, w]
        }
        ElementKind::Wedge6 => {
            let (u, v, _) = simplex(xi[0], xi[1], 0.0);
            [u, v, c(xi[2])]
        }
        ElementKind::Quad4 => [c(xi[0]), c(xi[1]), 0.0],
        _ => [c(xi[0]), c(xi[1]), c(xi[2])],
    }
}

fn combine(n: &[f64], x: &[P]) -> P {
    n.iter().zip(x).fold([0.0; 3], |s, (w, p)| add(s, scale(*p, *w)))
}

/// The parametric coordinates of `x` in the cell (Gauss–Newton on the linear map; for a shell,
/// of the nearest point of its surface), and the distance left.
fn locate(c: &RelaxationCell, x: P) -> (P, f64) {
    let (_, dim) = shape(c.kind, [0.0; 3]);
    let mut xi = match c.kind {
        ElementKind::Tri3 | ElementKind::Tet4 => [0.25; 3],
        ElementKind::Wedge6 => [1.0 / 3.0, 1.0 / 3.0, 0.5],
        _ => [0.5; 3],
    };
    if dim == 2 {
        xi[2] = 0.0;
    }
    let map = |xi: P| combine(&shape(c.kind, xi).0, &c.at);
    for _ in 0..20 {
        let r = sub(x, map(xi));
        // The Jacobian's columns (the map is at most trilinear: a unit-step difference of the
        // two sides of `xi` is exact along each axis).
        let h = STEP;
        let mut j = [[0.0; 3]; 3];
        for (k, col) in j.iter_mut().enumerate().take(dim) {
            let (mut a, mut b) = (xi, xi);
            a[k] += h;
            b[k] -= h;
            *col = scale(sub(map(a), map(b)), 0.5 / h);
        }
        // Normal equations JᵀJ δ = Jᵀ r.
        let mut m = [[0.0; 3]; 3];
        let mut rhs = [0.0; 3];
        for a in 0..dim {
            rhs[a] = dot(j[a], r);
            for b in 0..dim {
                m[a][b] = dot(j[a], j[b]);
            }
        }
        if dim == 2 {
            m[2][2] = 1.0;
        }
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if !(det.abs() > 1e-300) || !det.is_finite() {
            break;
        }
        let solve = |col: usize| {
            let mut t = m;
            for row in 0..3 {
                t[row][col] = rhs[row];
            }
            (t[0][0] * (t[1][1] * t[2][2] - t[1][2] * t[2][1]) - t[0][1] * (t[1][0] * t[2][2] - t[1][2] * t[2][0])
                + t[0][2] * (t[1][0] * t[2][1] - t[1][1] * t[2][0]))
                / det
        };
        let delta = [solve(0), solve(1), if dim == 3 { solve(2) } else { 0.0 }];
        xi = add(xi, delta);
        if !xi.iter().all(|v| v.is_finite()) || dot(delta, delta) < 1e-28 {
            break;
        }
    }
    (xi, dist(x, map(clamp_domain(c.kind, xi))))
}

/// The cell's size: its longest corner-to-corner span.
fn cell_size(c: &RelaxationCell) -> f64 {
    let mut s = 0.0f64;
    for a in &c.at {
        for b in &c.at {
            s = s.max(dist(*a, *b));
        }
    }
    s
}

/// Parametric slack accepted around a cell (a node of a finer mesh on a curved face lies off the
/// coarse cell's flat side by its chord error).
const SLACK: f64 = 0.25;

/// The field's displacement at `x`: from the cell holding it (or the nearest within the slack),
/// with that cell's shape functions; `None` outside the field.
pub fn field_displacement(cells: &[RelaxationCell], x: P) -> Option<P> {
    let mut best: Option<(f64, &RelaxationCell, P)> = None;
    for c in cells {
        let size = cell_size(c);
        let lo = c.at.iter().fold([f64::INFINITY; 3], |m, p| [m[0].min(p[0]), m[1].min(p[1]), m[2].min(p[2])]);
        let hi = c.at.iter().fold([f64::NEG_INFINITY; 3], |m, p| [m[0].max(p[0]), m[1].max(p[1]), m[2].max(p[2])]);
        let pad = SLACK * size;
        if (0..3).any(|k| x[k] < lo[k] - pad || x[k] > hi[k] + pad) {
            continue;
        }
        let (xi, gap) = locate(c, x);
        let badness = outside(c.kind, xi).max(gap / size.max(f64::MIN_POSITIVE));
        if badness <= SLACK && best.is_none_or(|(b, _, _)| badness < b) {
            best = Some((badness, c, xi));
        }
    }
    let (_, c, xi) = best?;
    let n = shape(c.kind, clamp_domain(c.kind, xi)).0;
    let d: Vec<P> = c.moved.iter().zip(&c.at).map(|(m, a)| sub(*m, *a)).collect();
    Some(combine(&n, &d))
}

/// Apply a relaxation field to the mesh: each node of `nodes` (default every node of a measured
/// element) inside the field moves by the interpolated displacement as `geometry` allows (on its
/// face or edge); a move that would make a valid element invalid is undone; then one guarded
/// smoothing sweep (strength `polish`, 0 for none) evens out what the interpolation left. The
/// model is not changed.
pub fn apply_field(
    model: &Model,
    cells: &[RelaxationCell],
    nodes: Option<&BTreeSet<u32>>,
    geometry: &dyn Constraints,
    polish: f64,
) -> Result<Relaxed> {
    let candidates: BTreeSet<u32> = match nodes {
        Some(n) => n.clone(),
        None => model.elements().values().filter(|e| quality::is_measured(e.kind)).flat_map(|e| e.nodes.iter().copied()).collect(),
    };
    if candidates.len() > fieldcraft_model::MAX_NODES {
        return Err(Error::Limit("too many nodes".into()));
    }
    let skin = skin_nodes(model, &candidates);
    let mut work = model.clone();
    let mut moved: BTreeMap<u32, P> = BTreeMap::new();
    for &n in &candidates {
        let Some(x) = model.node(n).map(|p| p.xyz) else { continue };
        let Some(d) = field_displacement(cells, x) else { continue };
        if !(dot(d, d) > 0.0) {
            continue;
        }
        let freedom = match geometry.freedom(n) {
            Some(f) => f,
            None if skin.contains(&n) => Freedom::Fixed,
            None => Freedom::Free,
        };
        let q = match freedom {
            Freedom::Fixed => None,
            Freedom::Free => Some(add(x, d)),
            Freedom::Face(f) => geometry.slide_face(f, x, d),
            Freedom::Curve(c) => geometry.slide_curve(c, x, d),
        };
        if let Some(q) = q.filter(|q| q.iter().all(|v| v.is_finite() && v.abs() <= fieldcraft_model::MAX_COORD)) {
            work.set_node(n, q)?;
            moved.insert(n, q);
        }
    }
    // Undo the moves around any element they would invalidate.
    let touched: Vec<(u32, ElementKind, Vec<u32>)> = model
        .elements()
        .iter()
        .filter(|(_, e)| quality::is_measured(e.kind) && e.nodes.iter().any(|n| moved.contains_key(n)))
        .map(|(l, e)| (*l, e.kind, e.nodes.clone()))
        .collect();
    let at = |m: &Model, nodes: &[u32]| -> Vec<P> { nodes.iter().map(|n| m.node(*n).map_or([0.0; 3], |p| p.xyz)).collect() };
    for _ in 0..10 {
        let mut reverted = false;
        for (label, kind, nodes) in &touched {
            let x0 = at(model, nodes);
            let e = Elem { label: *label, kind: *kind, nodes: nodes.clone(), normal: shell_normal(*kind, &x0), size0: 0.0 };
            if measure(&e, &x0).1 && !measure(&e, &at(&work, nodes)).1 {
                for n in nodes {
                    if moved.remove(n).is_some()
                        && let Some(p) = model.node(*n)
                    {
                        work.set_node(*n, p.xyz)?;
                        reverted = true;
                    }
                }
            }
        }
        if !reverted {
            break;
        }
    }
    let mut out = Relaxed::default();
    let corners = corner_nodes(model);
    if polish > 0.0 {
        let weights = uniform_weights(&work, &moved.keys().copied().filter(|n| corners.contains(n)).collect(), polish.min(1.0))?;
        if !weights.is_empty() {
            let r = relax(&work, &weights, 1, geometry)?;
            for (n, p) in &r.positions {
                work.set_node(*n, *p)?;
                moved.insert(*n, *p);
            }
            out.after = r.after;
        }
    }
    for (n, p) in moved {
        let Some(o) = model.node(n) else { continue };
        let d = dist(o.xyz, p);
        if d > 0.0 {
            out.max_displacement = out.max_displacement.max(d);
            out.moved += usize::from(corners.contains(&n));
            out.positions.insert(n, p);
        }
    }
    out.iterations = 1;
    Ok(out)
}

mod moving;
mod strong;
pub use moving::{
    DEFAULT_MIN_QUALITY, ElementChange, GUARD_STEPS, Guard, Influence, MIN_RELATIVE_RADIUS, MoveRun, Moved, Propagate, Unconstrained, centroid,
    move_nodes, move_nodes_with, move_weights, nearest_corner,
};
pub use strong::{DEFAULT_REACH_LENGTHS, Guards, MAX_ROUNDS, MorphOptions, Morphed, Round, morph};

#[cfg(test)]
mod tests;
