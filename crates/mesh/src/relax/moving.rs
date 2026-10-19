//! Moving nodes by hand (`mesh.node.move`): grabbed nodes are displaced by a vector, and the
//! nodes round them follow, so an analyst can drag a node of a finished mesh and see what the
//! move does to the elements around it.
//!
//! **Who follows.** Every node gets a share `w` in `[0, 1]` of the displacement ([`move_weights`]):
//!
//! - *Falloff*: the grabbed nodes move by all of it, a node within the radius of them by
//!   `strength × falloff(t)` (`t` its distance to the nearest grabbed node over the radius,
//!   [`Falloff`]); with a sphere and no grabbed nodes, every node inside it by
//!   `strength × falloff(t)` of its distance to the centre (a soft selection).
//! - *Harmonic*: the share is the harmonic function on the graph of element edges that is 1 at
//!   the grabbed nodes and 0 beyond the radius and at the nodes that may not move (Gauss–Seidel
//!   on the graph Laplacian, as the strong relaxation's morph: [`super::strong`]). The
//!   displacement then spreads through the mesh as a membrane would carry it, following the
//!   mesh's connectivity rather than straight-line distance (M. Botsch, O. Sorkine, *On linear
//!   variational surface deformation methods*, IEEE TVCG 14, 2008, §3: harmonic deformation
//!   fields with fixed and handle regions).
//!
//! **Where a node may go** is the relaxation's [`Constraints`]: a node on a geometry face slides
//! in the face's parameters, a node on a geometry edge along it, a node on a vertex stays; with
//! the geometry off ([`Unconstrained`]) every node moves in space. Mid-side nodes follow their
//! edge's corners, as in a relaxation.
//!
//! **Guard.** With [`Guard::NoInversion`] the displacement is scaled back (by bisection, at
//! most [`GUARD_STEPS`] halvings) to the largest share of it with which no valid element becomes
//! invalid or drops below the guard's smallest scaled Jacobian (default
//! [`DEFAULT_MIN_QUALITY`]; or its own, when that was lower); the scale is reported.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_model::Model;

use super::{Constraints, Falloff, Freedom, P, RegionQuality, State, add, corner_nodes, dist, dot, measure, scale, sub};
use crate::{Error, Result};

/// A node whose falloff weight is below this does not move.
const MIN_WEIGHT: f64 = 1e-12; // unit-factors: a weight from 0 to 1

/// How the displacement reaches the nodes round the grabbed ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Propagate {
    /// `strength × falloff` of the distance.
    #[default]
    Falloff,
    /// The harmonic function of the element-edge graph.
    Harmonic,
}

impl Propagate {
    pub fn from_name(name: &str) -> Option<Propagate> {
        match name {
            "falloff" => Some(Propagate::Falloff),
            "harmonic" => Some(Propagate::Harmonic),
            _ => None,
        }
    }
}

/// What keeps a move from breaking elements.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Guard {
    /// The move as asked, whatever it does to the elements.
    None,
    /// Scaled back until no valid element becomes invalid or falls below `min_quality` (scaled
    /// Jacobian; an element already below keeps at least its own).
    NoInversion { min_quality: f64 },
}

impl Default for Guard {
    fn default() -> Self {
        Guard::NoInversion { min_quality: DEFAULT_MIN_QUALITY }
    }
}

/// Halvings of the guard's bisection: the scale is found to within 2⁻²⁰.
pub const GUARD_STEPS: usize = 20;

/// The default smallest scaled Jacobian the guard lets a move make: with 0 the bisection stops at
/// the very edge of inversion, with an element of quality 0.
pub const DEFAULT_MIN_QUALITY: f64 = 0.05;

/// The smallest radius of influence, relative to the mesh's size (its bounding box diagonal):
/// below it no node but the grabbed ones is within reach, and the search grid would overflow.
pub const MIN_RELATIVE_RADIUS: f64 = 1e-6; // unit-factors: relative to the mesh's size

/// Every node moves in space, on the skin too: a move with the geometry off.
#[derive(Debug, Clone, Copy, Default)]
pub struct Unconstrained;

impl Constraints for Unconstrained {
    fn freedom(&self, _: u32) -> Option<Freedom> {
        Some(Freedom::Free)
    }
    fn slide_face(&self, _: u32, _: P, _: P) -> Option<P> {
        None
    }
    fn slide_curve(&self, _: u32, _: P, _: P) -> Option<P> {
        None
    }
}

/// Which nodes follow the grabbed ones, and how.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Influence {
    /// The sphere's centre: shares by the distance to it rather than to the grabbed nodes.
    pub centre: Option<P>,
    /// Nodes this close to the grabbed nodes (or the centre) follow; `None`: only the grabbed
    /// nodes move.
    pub radius: Option<f64>,
    pub propagate: Propagate,
    /// Falloff only: the largest share of a node that is not grabbed, in `(0, 1]`.
    pub strength: f64,
    pub falloff: Falloff,
}

impl Default for Influence {
    fn default() -> Self {
        Influence { centre: None, radius: None, propagate: Propagate::Falloff, strength: 1.0, falloff: Falloff::Smooth }
    }
}

/// For each candidate, its distance to the nearest of `from` when that is at most `r` (a uniform
/// grid of cell `r` over `from`; deterministic).
pub(super) fn nearest_within(from: &[P], candidates: &BTreeMap<u32, P>, r: f64) -> BTreeMap<u32, f64> {
    // Saturating casts and additions: no coordinate over radius overflows the index.
    let cell = |p: P| p.map(|x| (x / r).floor() as i64);
    let mut grid: BTreeMap<[i64; 3], Vec<P>> = BTreeMap::new();
    for &p in from {
        grid.entry(cell(p)).or_default().push(p);
    }
    let mut out = BTreeMap::new();
    for (&n, &p) in candidates {
        let [i, j, k] = cell(p);
        let mut best = f64::INFINITY;
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    for &q in grid.get(&[i.saturating_add(di), j.saturating_add(dj), k.saturating_add(dk)]).map(Vec::as_slice).unwrap_or(&[]) {
                        best = best.min(dist(p, q));
                    }
                }
            }
        }
        if best <= r {
            out.insert(n, best);
        }
    }
    out
}

/// The share of the displacement each node takes (node → share in `(0, 1]`): the grabbed
/// corner nodes `handles` all of it, the others as `influence` says. Nodes that may not move
/// under `geometry` (a vertex, a skin node of unknown geometry) are anchors of the harmonic
/// function and get no share. Only corner nodes of measured elements count (mid-side nodes
/// follow them).
pub fn move_weights(model: &Model, handles: &BTreeSet<u32>, influence: &Influence, geometry: &dyn Constraints) -> Result<BTreeMap<u32, f64>> {
    let Influence { centre, radius, propagate, strength, falloff } = *influence;
    if !(strength > 0.0 && strength <= 1.0) {
        return Err(Error::Invalid(format!("strength must be in (0, 1], got {strength}")));
    }
    if radius.is_some_and(|r| !(r > 0.0 && r.is_finite())) {
        return Err(Error::Invalid("the radius must be positive".into()));
    }
    if centre.is_some_and(|c| !c.iter().all(|x| x.is_finite())) {
        return Err(Error::Invalid("the centre must be finite".into()));
    }
    let corners = corner_nodes(model);
    let handles: BTreeSet<u32> = handles.iter().copied().filter(|n| corners.contains(n)).collect();
    if handles.len() > super::MAX_NODES {
        return Err(Error::Limit(format!("more than {} nodes to move", super::MAX_NODES)));
    }
    let positions: BTreeMap<u32, P> = corners.iter().filter_map(|&n| model.node(n).map(|p| (n, p.xyz))).collect();
    if let Some(r) = radius {
        let (lo, hi) = positions.values().fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |(lo, hi), p| {
            ([lo[0].min(p[0]), lo[1].min(p[1]), lo[2].min(p[2])], [hi[0].max(p[0]), hi[1].max(p[1]), hi[2].max(p[2])])
        });
        let size = if positions.is_empty() { 0.0 } else { dist(lo, hi) };
        let min = MIN_RELATIVE_RADIUS * size;
        if r < min {
            return Err(Error::invalid(crate::message![
                "the radius must be at least ",
                crate::message::length(min),
                " (a millionth of the mesh's size), got ",
                crate::message::length(r)
            ]));
        }
    }
    // The nodes within the radius, with their distance to the nearest handle (or the centre).
    let near = match (radius, centre) {
        (None, _) => BTreeMap::new(),
        (Some(r), Some(c)) => nearest_within(&[c], &positions, r),
        (Some(r), None) => {
            let from: Vec<P> = handles.iter().filter_map(|n| positions.get(n).copied()).collect();
            nearest_within(&from, &positions, r)
        }
    };
    if near.len() > super::MAX_NODES {
        return Err(Error::Limit(format!("more than {} nodes within the radius; give a smaller one", super::MAX_NODES)));
    }
    let mut out: BTreeMap<u32, f64> = handles.iter().map(|&n| (n, 1.0)).collect();
    match propagate {
        Propagate::Falloff => {
            let r = radius.unwrap_or(1.0);
            for (&n, &d) in &near {
                let w = strength * falloff.weight(d / r);
                if w > 0.0 && !handles.contains(&n) {
                    out.insert(n, w);
                }
            }
        }
        Propagate::Harmonic => {
            if radius.is_none() {
                return Err(Error::Invalid("a harmonic move needs a radius".into()));
            }
            if handles.is_empty() {
                return Err(Error::Invalid("a harmonic move needs nodes to grab".into()));
            }
            let movers: BTreeSet<u32> = near.keys().chain(&handles).copied().collect();
            let state = State::new(model, &movers, geometry, false)?;
            let free: Vec<u32> = near.keys().copied().filter(|n| !handles.contains(n) && state.freedom(*n) != Freedom::Fixed).collect();
            let fixed: BTreeMap<u32, P> = handles.iter().map(|&n| (n, [1.0, 0.0, 0.0])).collect();
            let u = state.harmonic(&free, &fixed)?;
            for n in free {
                let w = u.get(&n).map_or(0.0, |v| v[0]).clamp(0.0, 1.0);
                if w > MIN_WEIGHT {
                    out.insert(n, w);
                }
            }
        }
    }
    Ok(out)
}

/// One element's quality before and after a move.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElementChange {
    pub label: u32,
    /// Scaled Jacobian (a triangle's smallest angle over 60°).
    pub before: f64,
    pub after: f64,
    pub valid_before: bool,
    pub valid_after: bool,
}

impl ElementChange {
    /// Its quality dropped, or it became invalid.
    pub fn worse(&self) -> bool {
        self.after < self.before - SAME || (self.valid_before && !self.valid_after)
    }
    /// Its quality rose.
    pub fn better(&self) -> bool {
        self.after > self.before + SAME && !self.worse()
    }
}

/// What [`move_nodes`] did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Moved {
    /// The new position of every node that moved (corners and the mid-side nodes following them).
    pub positions: BTreeMap<u32, P>,
    /// Corner nodes that moved.
    pub moved: usize,
    pub max_displacement: f64,
    /// The quality of the elements with a corner that may move, before and after.
    pub before: RegionQuality,
    pub after: RegionQuality,
    pub mean_before: f64,
    pub mean_after: f64,
    /// Elements whose quality dropped, and rose.
    pub worse: usize,
    pub better: usize,
    /// The share of the asked displacement applied (below 1 when the guard scaled it back).
    pub scale: f64,
    /// Every element of the region, by label.
    pub elements: Vec<ElementChange>,
}

/// A change in quality smaller than this is no change.
const SAME: f64 = 1e-9; // unit-factors: a difference of qualities

/// Move each node of `weights` by its share of `d`, as `geometry` lets it, the displacement
/// scaled back under `guard`. The model is not changed: the result has the new positions.
pub fn move_nodes(model: &Model, weights: &BTreeMap<u32, f64>, d: P, guard: Guard, geometry: &dyn Constraints) -> Result<Moved> {
    move_nodes_with(model, weights, d, guard, geometry, &MoveRun::default())
}

/// How [`move_nodes_with`] runs: what a live preview reuses and how finely it guards.
#[derive(Debug, Clone, Default)]
pub struct MoveRun {
    /// The region's elements' quality and validity before any move, in label order (from an
    /// earlier run's [`Moved::elements`] with the same weights); measured again when absent or
    /// of another length.
    pub baseline: Option<Vec<(f64, bool)>>,
    /// Halvings of the guard's bisection; 0 is [`GUARD_STEPS`].
    pub guard_steps: usize,
}

/// Each element's quality and validity at the state's positions, measured in parallel (fixed
/// blocks: the same result on any thread count).
fn measure_all(elems: &[super::Elem], pos: &BTreeMap<u32, P>) -> Vec<(f64, bool)> {
    fieldcraft_parallel::map(elems.len(), |i| {
        let Some(e) = elems.get(i) else { return (-1.0, false) };
        let x: Vec<P> = e.nodes.iter().map(|n| pos.get(n).copied().unwrap_or([0.0; 3])).collect();
        measure(e, &x)
    })
}

/// The region's quality from its elements' measures.
fn summary(state: &State<'_>, values: &[(f64, bool)]) -> RegionQuality {
    let mut r = RegionQuality { elements: state.elems.len(), min_quality: f64::INFINITY, max_edge_ratio: 0.0, invalid: 0, worst: Vec::new() };
    let mut all = Vec::with_capacity(values.len());
    for (e, &(q, valid)) in state.elems.iter().zip(values) {
        r.min_quality = r.min_quality.min(q);
        r.max_edge_ratio = r.max_edge_ratio.max(super::edge_ratio(e.kind, &state.coords(e)));
        r.invalid += usize::from(!valid);
        all.push((e.label, q));
    }
    all.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    all.truncate(super::WORST_LISTED);
    r.worst = all;
    if r.elements == 0 {
        r.min_quality = 0.0;
    }
    r
}

/// [`move_nodes`], reusing what `run` carries.
pub fn move_nodes_with(model: &Model, weights: &BTreeMap<u32, f64>, d: P, guard: Guard, geometry: &dyn Constraints, run: &MoveRun) -> Result<Moved> {
    if !d.iter().all(|x| x.is_finite()) || dot(d, d).sqrt() > fieldcraft_model::MAX_COORD {
        return Err(Error::Invalid("the displacement must be finite".into()));
    }
    if weights.len() > super::MAX_NODES {
        return Err(Error::Limit(format!("more than {} nodes to move", super::MAX_NODES)));
    }
    let floor = match guard {
        Guard::None => 0.0,
        Guard::NoInversion { min_quality } if (0.0..=1.0).contains(&min_quality) => min_quality,
        Guard::NoInversion { min_quality } => return Err(Error::Invalid(format!("the smallest quality must be in [0, 1], got {min_quality}"))),
    };
    if let Some((n, w)) = weights.iter().find(|(_, w)| !(**w > 0.0 && **w <= 1.0)) {
        return Err(Error::Invalid(format!("node {n} has share {w}; shares are in (0, 1]")));
    }
    let steps = if run.guard_steps == 0 { GUARD_STEPS } else { run.guard_steps.min(GUARD_STEPS) };
    let movers: BTreeSet<u32> = weights.keys().copied().collect();
    let mut state = State::new(model, &movers, geometry, false)?;
    let original = state.pos.clone();
    let start = match &run.baseline {
        Some(b) if b.len() == state.elems.len() => b.clone(),
        _ => measure_all(&state.elems, &state.pos),
    };
    let before = summary(&state, &start);
    let freedom: BTreeMap<u32, Freedom> = movers.iter().map(|&n| (n, state.freedom(n))).collect();
    // The state with the displacement applied at scale `t`.
    let apply = |state: &mut State<'_>, t: f64| {
        state.pos.clone_from(&original);
        for (&n, &w) in weights {
            let (Some(&f), Some(&p0)) = (freedom.get(&n), original.get(&n)) else { continue };
            if !state.around.contains_key(&n) {
                continue;
            }
            if let Some(q) = state.moved(f, p0, scale(d, w * t)) {
                state.place(n, q);
            }
        }
    };
    // Whether every element valid before is still valid, and no worse than the guard's floor
    // (or its own quality, when that was lower).
    let sound = |now: &[(f64, bool)]| now.iter().zip(&start).all(|(&(q, valid), &(q0, was))| !was || (valid && q >= floor.min(q0)));
    apply(&mut state, 1.0);
    let mut now = measure_all(&state.elems, &state.pos);
    let mut t = 1.0;
    if guard != Guard::None && !sound(&now) {
        let (mut lo, mut hi) = (0.0, 1.0);
        let mut at_lo = None;
        for _ in 0..steps {
            fieldcraft_model::progress::check()?;
            let mid = 0.5 * (lo + hi);
            apply(&mut state, mid);
            let m = measure_all(&state.elems, &state.pos);
            if sound(&m) {
                lo = mid;
                at_lo = Some(m);
            } else {
                hi = mid;
            }
        }
        t = lo;
        apply(&mut state, t);
        now = at_lo.unwrap_or_else(|| start.clone());
    }
    let after = summary(&state, &now);
    let mut out = Moved { before, after, scale: t, ..Moved::default() };
    let (mut sum_before, mut sum_after) = (0.0, 0.0);
    for (e, (&(q0, v0), &(q1, v1))) in state.elems.iter().zip(start.iter().zip(&now)) {
        sum_before += q0;
        sum_after += q1;
        let change = ElementChange { label: e.label, before: q0, after: q1, valid_before: v0, valid_after: v1 };
        out.worse += usize::from(change.worse());
        out.better += usize::from(change.better());
        out.elements.push(change);
    }
    if !state.elems.is_empty() {
        out.mean_before = sum_before / state.elems.len() as f64;
        out.mean_after = sum_after / state.elems.len() as f64;
    }
    for (n, p) in &state.pos {
        let Some(o) = original.get(n) else { continue };
        let moved = dist(*o, *p);
        if moved > 0.0 {
            out.max_displacement = out.max_displacement.max(moved);
            out.moved += usize::from(movers.contains(n));
            out.positions.insert(*n, *p);
        }
    }
    Ok(out)
}

/// The centroid of `nodes`' positions (`None` when none has one): where a move `to` a point
/// starts from.
pub fn centroid(model: &Model, nodes: &BTreeSet<u32>) -> Option<P> {
    let pts: Vec<P> = nodes.iter().filter_map(|n| model.node(*n).map(|p| p.xyz)).collect();
    if pts.is_empty() {
        return None;
    }
    Some(scale(pts.iter().fold([0.0; 3], |s, p| add(s, *p)), 1.0 / pts.len() as f64))
}

/// The corner node of a measured element nearest `p` (the node a sphere's harmonic move grabs).
pub fn nearest_corner(model: &Model, p: P) -> Option<u32> {
    let mut best: Option<(f64, u32)> = None;
    for n in corner_nodes(model) {
        let Some(x) = model.node(n) else { continue };
        let d = dot(sub(x.xyz, p), sub(x.xyz, p));
        if best.is_none_or(|(b, _)| d < b) {
            best = Some((d, n));
        }
    }
    best.map(|(_, n)| n)
}
