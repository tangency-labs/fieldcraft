//! Moving nodes by hand: `mesh.node.move` displaces grabbed nodes (listed, on listed B-rep faces,
//! or a sphere's) by a vector or to a point, the nodes round them following by a falloff or by
//! harmonic interpolation ([`fieldcraft_mesh::relax::move_weights`]), nodes on the geometry
//! staying on it as Relax Nodes moves them ([`super::relax`]), the move scaled back until no
//! element inverts; it reports what the move did to the elements' quality.
//!
//! The interactive tool (drag a node in the view) works out the same move live: [`prepare`] does
//! the part that does not depend on the displacement once (which nodes follow, and the design's
//! faces they slide on), and [`Prepared::preview`] applies a displacement to it.

use std::collections::{BTreeMap, BTreeSet};

use fieldcraft_mesh::relax::{self, Constraints, Falloff, Guard, Influence, Moved, Propagate, Unconstrained};
use serde_json::{Value, json};

use super::relax::{DesignGeometry, design_geometry, face_refs, nodes_on_faces, point, quality_json, unmeshed_faces};
use crate::requirement::Requirement;
use crate::schema::ProblemKind;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const MOVE: CommandSpec = CommandSpec {
    id: "mesh.node.move",
    label: "Move Nodes",
    params: r#"{"nodes": [labels] (the nodes to grab) | "faces": [{"body": "Base", "face": "F3"}] (every node on these B-rep faces) | "center": [x, y, z] with "radius": length (a sphere: every node inside follows by the falloff of its distance to the centre), "by": [dx, dy, dz] | "to": [x, y, z] (where the grabbed nodes' centroid, or the centre, goes), "radius"?: length (with nodes or faces: the nodes this close to a grabbed one follow; default none, only the grabbed nodes move), "propagate"?: "falloff" (default) | "harmonic" (the displacement carried into the nodes within the radius by harmonic interpolation on the graph of element edges, 1 at the grabbed nodes and zero beyond the radius; needs a radius, and the grabbed nodes: a sphere grabs its nearest node), "falloff"?: "smooth" (default) | "linear" | "none", "strength"?: 1 (0..1]: the largest share of the displacement a node that is not grabbed takes (falloff only), "on_geometry"?: true (nodes on a B-rep face slide on it in its UV parameters, on a B-rep edge along it, on a vertex stay, and skin nodes of a mesh without geometry stay; false: every node moves in space), "guard"?: "no_inversion" (default: the move is scaled back until no element inverts or falls below `min_quality`) | "none", "min_quality"?: 0.05 (0..1, no_inversion only: the smallest scaled Jacobian the move may leave an element with; an element already below keeps at least its own)} — moves the nodes; mid-side nodes follow. Reports the nodes moved, the largest displacement, the displacement applied and its scale (below 1 when the guard scaled it back), the quality of the elements round the moved nodes before and after (smallest and mean scaled Jacobian, largest edge ratio, invalid elements, the worst) and how many got worse and better. One undo step; a move that moves nothing (a zero displacement, or one the guard scales back to nothing) changes nothing and records no undo step (`"moved": 0`, `"nothing_moved": true`)"#,
    form: Some(super::forms::NODE_MOVE),
    kind: CommandKind::Model,
    requires: &[Requirement::Mesh],
    run: move_nodes,
};

/// Where nodes may go: the design's faces and edges, or anywhere.
enum Geometry {
    Design(DesignGeometry),
    Free,
}

impl Geometry {
    fn constraints(&self) -> &dyn Constraints {
        match self {
            Geometry::Design(g) => g,
            Geometry::Free => &Unconstrained,
        }
    }
}

/// A move worked out up to its displacement: who follows and by how much, and where they may go.
pub struct Prepared {
    /// Node → its share of the displacement.
    weights: BTreeMap<u32, f64>,
    /// The grabbed nodes.
    handles: BTreeSet<u32>,
    /// Where a move `to` a point starts from: the grabbed nodes' centroid, or the centre.
    origin: Option<[f64; 3]>,
    /// The sphere of influence round the centre or the grabbed nodes, for the view.
    pub sphere: Option<([f64; 3], f64)>,
    guard: Guard,
    geometry: Geometry,
    /// The region's elements' quality before the move, measured by the first live preview.
    baseline: std::sync::OnceLock<Vec<(f64, bool)>>,
    /// The mesh the move works on, in millimetres (the meshers' system).
    model: std::sync::Arc<fieldcraft_model::Model>,
    /// For a tool's preview ([`prepare`]): the system its parameters are typed in. Its
    /// parameters are then brought to millimetres, and what it shows back to SI.
    typed_in: Option<fieldcraft_model::DisplaySystem>,
}

/// Halvings of the guard's bisection in a live preview (the scale to within 2⁻⁶; the command
/// itself finds it to 2⁻²⁰).
const PREVIEW_GUARD_STEPS: usize = 6;

impl Prepared {
    /// Nodes that take a share of the displacement.
    pub fn nodes(&self) -> usize {
        self.weights.len()
    }

    /// The grabbed nodes.
    pub fn handles(&self) -> &BTreeSet<u32> {
        &self.handles
    }

    /// The displacement `params` ask for (`by`, or `to` less the origin), in millimetres.
    pub fn displacement(&self, params: &Value) -> Result<[f64; 3]> {
        let converted;
        let params = match &self.typed_in {
            Some(input) => {
                converted = crate::units::convert_params(super::forms::NODE_MOVE, params, input, crate::units::geometry())?;
                &converted
            }
            None => params,
        };
        match (point(params, "by")?, point(params, "to")?) {
            (Some(_), Some(_)) => Err(EngineError::Params("give `by` or `to`, not both".into())),
            (Some(d), None) => Ok(d),
            (None, Some(to)) => {
                let o = self.origin.ok_or_else(|| EngineError::Params("no node to move `to` a point".into()))?;
                Ok([to[0] - o[0], to[1] - o[1], to[2] - o[2]])
            }
            (None, None) => Err(EngineError::Params("say how far to move: give `by` [dx, dy, dz] or `to` [x, y, z]".into())),
        }
    }

    /// The move by `d` ([`Prepared::displacement`]) on the mesh it was prepared on, without doing
    /// it; for a tool's preview, its positions and largest displacement in SI.
    pub fn preview(&self, _s: &Session, d: [f64; 3]) -> Result<Moved> {
        let m = relax::move_nodes(&self.model, &self.weights, d, self.guard, self.geometry.constraints())?;
        Ok(self.shown(m))
    }

    /// [`Prepared::preview`] for a drag frame: the region's quality before the move is measured
    /// once and reused, and the guard's scale is found more coarsely.
    pub fn preview_live(&self, _s: &Session, d: [f64; 3]) -> Result<Moved> {
        let run = relax::MoveRun { baseline: self.baseline.get().cloned(), guard_steps: PREVIEW_GUARD_STEPS };
        let m = relax::move_nodes_with(&self.model, &self.weights, d, self.guard, self.geometry.constraints(), &run)?;
        if self.baseline.get().is_none() {
            let _ = self.baseline.set(m.elements.iter().map(|c| (c.before, c.valid_before)).collect());
        }
        Ok(self.shown(m))
    }

    /// A move as the session (SI) has it, for a tool's preview.
    fn shown(&self, mut m: Moved) -> Moved {
        if self.typed_in.is_some() {
            for p in m.positions.values_mut() {
                *p = crate::units::geometry_point(*p);
            }
            m.max_displacement = crate::units::geometry_length(m.max_displacement);
        }
        m
    }
}

/// The corner node nearest `p` (SI, as the session stores the mesh): the node a drag in the
/// view grabs.
pub fn nearest_node(s: &Session, p: [f64; 3]) -> Option<u32> {
    relax::nearest_corner(&s.model, p)
}

/// A positive number parameter, when given.
fn positive(p: &Value, key: &str) -> Result<Option<f64>> {
    match p.get(key) {
        None => Ok(None),
        Some(_) => Some(super::optional_number(p, key, 0.0)?)
            .filter(|r| *r > 0.0 && r.is_finite())
            .ok_or_else(|| EngineError::param(key, ProblemKind::NotPositive))
            .map(Some),
    }
}

/// Everything of `mesh.node.move` but its displacement, for the tool's preview: `p` as typed (in
/// the display system), the session as stored (SI); what it shows comes back in SI.
pub fn prepare(s: &Session, p: &Value) -> Result<Prepared> {
    match prepare_capped(s, p, None)? {
        Capped::Ready(prep) => Ok(*prep),
        Capped::TooLarge(n) => Err(EngineError::Params(format!("{n} nodes to move"))),
    }
}

/// [`prepare`], or the count of the nodes that would follow when it is more than `cap`: the
/// costly part (the design's faces, the harmonic solve) is then not done.
pub enum Capped {
    Ready(Box<Prepared>),
    TooLarge(usize),
}

/// [`prepare`] that stops before the costly part when more than `cap` nodes would follow.
pub fn prepare_capped(s: &Session, p: &Value, cap: Option<usize>) -> Result<Capped> {
    let framed = crate::units::geometry_session(s)?;
    let params = crate::units::geometry_params(s, MOVE.id, p)?;
    Ok(match prepare_framed(&framed, &params, cap)? {
        Capped::Ready(mut prep) => {
            prep.typed_in = Some(s.input_units().clone());
            prep.sphere = prep.sphere.map(|(c, r)| (crate::units::geometry_point(c), crate::units::geometry_length(r)));
            Capped::Ready(prep)
        }
        c => c,
    })
}

/// [`prepare_capped`] on a session whose model is in millimetres, with parameters in them (the
/// command's own frame).
fn prepare_framed(s: &Session, p: &Value, cap: Option<usize>) -> Result<Capped> {
    let model = &s.model;
    let propagate = match super::opt::<String>(p, "propagate")? {
        None => Propagate::Falloff,
        Some(v) => Propagate::from_name(&v).ok_or_else(|| EngineError::param("propagate", ProblemKind::NotAChoice))?,
    };
    if propagate == Propagate::Harmonic && ["falloff", "strength"].iter().any(|k| p.get(*k).is_some()) {
        return Err(EngineError::Params("`falloff` and `strength` go with \"propagate\": \"falloff\"".into()));
    }
    let falloff = match super::opt::<String>(p, "falloff")? {
        None => Falloff::Smooth,
        Some(f) => Falloff::from_name(&f).ok_or_else(|| EngineError::param("falloff", ProblemKind::NotAChoice))?,
    };
    let strength = super::optional_number(p, "strength", 1.0)?;
    if !(strength > 0.0 && strength <= 1.0) {
        return Err(EngineError::param("strength", ProblemKind::NotInRange { min: 0.0, max: 1.0 }));
    }
    let min_quality = super::optional_number(p, "min_quality", relax::DEFAULT_MIN_QUALITY)?;
    if !(0.0..=1.0).contains(&min_quality) {
        return Err(EngineError::param("min_quality", ProblemKind::NotInRange { min: 0.0, max: 1.0 }));
    }
    let guard = match super::opt::<String>(p, "guard")?.as_deref() {
        None | Some("no_inversion") => Guard::NoInversion { min_quality },
        Some("none") if p.get("min_quality").is_some() => {
            return Err(EngineError::Params("`min_quality` goes with \"guard\": \"no_inversion\"".into()));
        }
        Some("none") => Guard::None,
        Some(_) => return Err(EngineError::param("guard", ProblemKind::NotAChoice)),
    };
    let on_geometry = super::opt::<bool>(p, "on_geometry")?.unwrap_or(true);
    let radius = positive(p, "radius")?;
    let centre = point(p, "center")?;
    let listed = ["nodes", "faces"].iter().filter(|k| p.get(**k).is_some()).count() + usize::from(centre.is_some());
    if listed != 1 {
        return Err(EngineError::Params("say which nodes to move: give one of `nodes`, `faces` or `center` (with `radius`)".into()));
    }
    let mut handles: BTreeSet<u32> = BTreeSet::new();
    if p.get("nodes").is_some() {
        for n in super::labels(p, "nodes", fieldcraft_model::MAX_NODES)? {
            if model.node(n).is_none() {
                return Err(EngineError::Params(format!("no node {n}")));
            }
            handles.insert(n);
        }
    }
    if p.get("faces").is_some() {
        handles.extend(nodes_on_faces(model, &face_refs(p)?)?);
    }
    let mut origin = relax::centroid(model, &handles);
    if let Some(c) = centre {
        if radius.is_none() {
            return Err(EngineError::Params("a sphere needs both `center` and `radius`".into()));
        }
        origin = Some(c);
        // A harmonic move needs something to carry: the sphere grabs the node nearest its centre.
        if propagate == Propagate::Harmonic {
            handles.extend(relax::nearest_corner(model, c));
        }
    }
    let influence = Influence { centre: if propagate == Propagate::Falloff { centre } else { None }, radius, propagate, strength, falloff };
    // The constraints of every node that might take a share: the grabbed ones and those within
    // the radius (all of the mesh's corner nodes when the region is not known yet).
    let candidates =
        relax::move_weights(model, &handles, &Influence { propagate: Propagate::Falloff, falloff: Falloff::None, ..influence }, &Unconstrained)?;
    if let Some(cap) = cap
        && candidates.len() > cap
    {
        return Ok(Capped::TooLarge(candidates.len()));
    }
    let geometry = if on_geometry {
        let mut nodes: BTreeSet<u32> = candidates.keys().copied().collect();
        // Their elements' nodes too: the mid-side nodes that follow.
        for e in model.elements().values() {
            if e.nodes.iter().any(|n| candidates.contains_key(n)) {
                nodes.extend(e.nodes.iter().copied());
            }
        }
        Geometry::Design(design_geometry(s, &nodes)?)
    } else {
        Geometry::Free
    };
    let weights = relax::move_weights(model, &handles, &influence, geometry.constraints())?;
    if weights.is_empty() {
        return Err(EngineError::Params("no mesh nodes to move: the nodes given are no element's corners, or none lies in the sphere".into()));
    }
    let sphere = radius.and_then(|r| centre.or(origin).map(|c| (c, r)));
    Ok(Capped::Ready(Box::new(Prepared {
        weights,
        handles,
        origin,
        sphere,
        guard,
        geometry,
        baseline: std::sync::OnceLock::new(),
        model: std::sync::Arc::new(model.clone()),
        typed_in: None,
    })))
}

/// The quality of the region as JSON, with its mean scaled Jacobian.
fn quality(q: &relax::RegionQuality, mean: f64) -> Value {
    let mut v = quality_json(q);
    if let Some(o) = v.as_object_mut() {
        o.insert("mean_scaled_jacobian".into(), json!(mean));
    }
    v
}

fn move_nodes(s: &mut Session, p: &Value) -> Result<Value> {
    let prep = match prepare_framed(s, p, None)? {
        Capped::Ready(prep) => *prep,
        Capped::TooLarge(n) => return Err(EngineError::Params(format!("{n} nodes to move"))),
    };
    let d = prep.displacement(p)?;
    let r = prep.preview(s, d)?;
    let asked = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if r.moved == 0 && asked > 0.0 && r.scale > 0.0 {
        return Err(EngineError::Params(
            "the grabbed nodes may not move: they are on B-rep vertices (or on the skin of a mesh without geometry); give \"on_geometry\": false to move them freely".into(),
        ));
    }
    if r.positions.is_empty() {
        // Nothing to do: the model is unchanged, so no undo step is recorded.
        let why = if asked == 0.0 {
            "the displacement is zero"
        } else {
            "the guard scaled the move back to nothing (any move would invert an element or go below min_quality)"
        };
        return Ok(json!({
            "moved": 0,
            "nodes_changed": 0,
            "nothing_moved": true,
            "scale": r.scale,
            "before": quality(&r.before, r.mean_before),
            "summary": format!("nothing moved: {why}"),
        }));
    }
    for (n, x) in &r.positions {
        s.model.set_node(*n, *x)?;
    }
    let summary = format!(
        "moved {} node{} (largest {}{}); scaled Jacobian min {:.3} to {:.3}, mean {:.3} to {:.3}; {} element{} worse, {} better",
        r.moved,
        if r.moved == 1 { "" } else { "s" },
        crate::units::shown_length(s, r.max_displacement),
        if r.scale < 1.0 { format!(", scaled to {:.3} so that no element inverts", r.scale) } else { String::new() },
        r.before.min_quality,
        r.after.min_quality,
        r.mean_before,
        r.mean_after,
        r.worse,
        if r.worse == 1 { "" } else { "s" },
        r.better,
    );
    Ok(json!({
        "moved": r.moved,
        "nodes_changed": r.positions.len(),
        "grabbed": prep.handles.len(),
        "max_displacement": r.max_displacement,
        "displacement": d.map(|x| x * r.scale),
        "scale": r.scale,
        "before": quality(&r.before, r.mean_before),
        "after": quality(&r.after, r.mean_after),
        "worse": r.worse,
        "better": r.better,
        "sphere": prep.sphere.map(|(c, r)| json!({"center": c, "radius": r})),
        "faces_without_mesh": unmeshed_faces(&s.model, p),
        "summary": summary,
    }))
}

#[cfg(test)]
mod tests;
