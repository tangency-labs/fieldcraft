//! Multi-level sweeping along a straight direction (extrusions) or around an axis (revolutions).
//!
//! **Recognition.** In a sweep frame every point has cross-section coordinates `q` and a sweep
//! coordinate `h` (height along the direction, or angle about the axis). A body is sweepable in
//! the frame when every face is a *cap* (a plane across the sweep: normal along the direction, or
//! a half-plane through the axis) or a *wall* (along the sweep: planes and cylinders parallel to
//! the direction; planes across the axis and coaxial cylinders, cones, spheres and tori), and
//! every edge runs either across the sweep (constant `h`) or along it (constant `q`). The caps
//! and the across edges sit at a finite set of *levels*.
//!
//! **Decomposition.** Projecting into the cross-section every edge that starts a wall upward
//! (the lower outline of each wall) and imprinting them on each other gives a planar arrangement
//! ([`super::arrangement`]) whose regions are the columns of the body: over each region the body
//! is a union of whole level intervals. This is many-to-many sweeping by imprinting all sources
//! and targets into one projection plane (Lai, Benzley, White, IJNME 49(1–2):261–275, 2000;
//! Knupp, *Next-generation sweep tool: a method for generating all-hex meshes on two-and-one-half
//! dimensional geometries*, Proc. 7th International Meshing Roundtable, 1998). Which intervals a
//! region fills is decided by the inside test at its centre ([`super::geom::Inside`]).
//!
//! **Blocks.** Regions no structured scheme fits are cut into blocks that fit
//! ([`super::decompose`]: holes apart, O-grids, cuts from concave corners and fillets), and the
//! arrangement is laid out again with the cuts, so they are curves like the others.
//!
//! **Meshing.** Interval assignment ([`super::interval`]) gives every curve its count so that
//! every region's scheme fits ([`super::section`]: four-sided regions mapped with opposite sides
//! equal, three-sided ones the triangle primitive, rings and disks O-grids, the rest even loops
//! for the free all-quad mesher); each level interval gets its own layer count. Each region's quads are
//! stacked through its intervals; shared curves and levels share nodes, so the hexahedra
//! conform. Nodes of curves on analytic walls are projected onto them; in a revolution, the
//! rotated nodes then lie exactly on every surface of revolution.

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_4, PI, TAU};

use super::arrangement::{self, Arrangement, Q};
use super::brep::{Brep, Face};
use super::decompose;
use super::geom::{Inside, Surface, add, dist, scale, unit};
use super::interval::Intervals;
use super::region::{self, Scheme};
use super::section;
use super::{Ctx, Raw};
use crate::message;
use crate::message::Message;
use crate::predicates::Sign;
use crate::seeds::EdgeSeed;
use crate::size::{SizeField, divide_polyline, edge_divisions, ideal_divisions};
use crate::vec3::{P, cross, dot, norm, sub};
use crate::{Error, Result};
use fieldcraft_model::SeedConstraint;

/// Angular tolerance (radians) of "parallel" and "perpendicular" for the recognition.
const ANGLE_TOL: f64 = 1e-6;
/// Distance (relative to the diagonal) within which two axes of revolution are one.
const RELATIVE_AXIS: f64 = 1e-6;
/// The most candidate frames tried per method.
pub(crate) const MAX_CANDIDATES: usize = 6;

/// A sweep frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Frame {
    /// Cross-section `(p·u, p·v)` relative to `o`, height `p·d`; `u × v = d`.
    Translate { o: P, u: P, v: P, d: P },
    /// Cross-section `(z, r)` (along the axis, from the axis), angle about `a` from `e1`
    /// counter-clockwise; `e2 = a × e1`.
    Rotate { o: P, a: P, e1: P, e2: P },
}

fn perpendicular(d: P) -> P {
    let k = if d[0].abs() <= d[1].abs() && d[0].abs() <= d[2].abs() {
        [1.0, 0.0, 0.0]
    } else if d[1].abs() <= d[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    unit(cross(d, k)).unwrap_or([1.0, 0.0, 0.0])
}

impl Frame {
    pub(crate) fn translate(o: P, d: P) -> Option<Frame> {
        let d = unit(d)?;
        let u = perpendicular(d);
        Some(Frame::Translate { o, u, v: cross(d, u), d })
    }

    pub(crate) fn rotate(o: P, a: P) -> Option<Frame> {
        let a = unit(a)?;
        let e1 = perpendicular(a);
        Some(Frame::Rotate { o, a, e1, e2: cross(a, e1) })
    }

    pub(crate) fn coords(&self, p: P) -> (Q, f64) {
        match *self {
            Frame::Translate { o, u, v, d } => {
                let w = sub(p, o);
                ([dot(w, u), dot(w, v)], dot(w, d))
            }
            Frame::Rotate { o, a, e1, e2 } => {
                let w = sub(p, o);
                let z = dot(w, a);
                let (x, y) = (dot(w, e1), dot(w, e2));
                let t = crate::vec3::atan2(y, x);
                ([z, libm::hypot(x, y)], if t < 0.0 { t + TAU } else { t })
            }
        }
    }

    pub(crate) fn point(&self, q: Q, h: f64) -> P {
        match *self {
            Frame::Translate { o, u, v, d } => add(o, add(add(scale(u, q[0]), scale(v, q[1])), scale(d, h))),
            Frame::Rotate { o, a, e1, e2 } => add(o, add(scale(a, q[0]), scale(add(scale(e1, libm::cos(h)), scale(e2, libm::sin(h))), q[1]))),
        }
    }

    fn periodic(&self) -> bool {
        matches!(self, Frame::Rotate { .. })
    }

    /// The signed sweep distance from `h0` to `h1` (angles wrapped to (−π, π]).
    fn dh(&self, h1: f64, h0: f64) -> f64 {
        let d = h1 - h0;
        match self {
            Frame::Translate { .. } => d,
            Frame::Rotate { .. } => {
                let w = d.rem_euclid(TAU);
                if w > PI { w - TAU } else { w }
            }
        }
    }

    /// The length a sweep distance `dh` spans at cross-section point `q`.
    fn span(&self, q: Q, dh: f64) -> f64 {
        match self {
            Frame::Translate { .. } => dh.abs(),
            Frame::Rotate { .. } => dh.abs() * q[1].abs(),
        }
    }

    pub(crate) fn describe(&self) -> Message {
        let f = |p: P| format!("({:.6}, {:.6}, {:.6})", p[0], p[1], p[2]);
        match *self {
            Frame::Translate { d, .. } => message![format!("extruded along {}", f(d))],
            Frame::Rotate { o, a, .. } => message!["revolved about the axis through ", crate::message::point(o), format!(" along {}", f(a))],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Role {
    Cap,
    Wall,
}

fn parallel(a: P, b: P) -> bool {
    norm(cross(a, b)) <= ANGLE_TOL
}

fn perpendicular_to(a: P, b: P) -> bool {
    dot(a, b).abs() <= ANGLE_TOL
}

/// Distance from `p` to the line through `o` along unit `a`.
fn line_distance(p: P, o: P, a: P) -> f64 {
    let w = sub(p, o);
    norm(sub(w, scale(a, dot(w, a))))
}

/// The face's role in the frame, or why it has none.
fn role(brep: &Brep, face: &Face, frame: &Frame, tol: f64) -> std::result::Result<Role, String> {
    let name = face.id;
    match *frame {
        Frame::Translate { d, .. } => {
            if let Some((_, n)) = face.plane {
                if parallel(n, d) {
                    return Ok(Role::Cap);
                }
                if perpendicular_to(n, d) {
                    return Ok(Role::Wall);
                }
                return Err(format!("planar face {name} is slanted to the direction"));
            }
            match face.surface {
                Some(Surface::Cylinder { axis, .. }) if parallel(axis, d) => Ok(Role::Wall),
                Some(Surface::Cylinder { .. }) => Err(format!("cylindrical face {name} is not parallel to the direction")),
                Some(_) => Err(format!("face {name} is curved across the direction")),
                None => {
                    // A free-form face is a wall when every triangle is parallel to the direction.
                    let ok = face.triangles.iter().all(|&t| {
                        let Some(tri) = brep.surface.triangles().get(t as usize) else { return false };
                        let [a, b, c] = tri.map(|v| brep.point(v));
                        unit(cross(sub(b, a), sub(c, a))).is_some_and(|n| perpendicular_to(n, d))
                    });
                    if ok { Ok(Role::Wall) } else { Err(format!("free-form face {name} is not parallel to the direction")) }
                }
            }
        }
        Frame::Rotate { o, a, .. } => {
            let coaxial = |origin: P, axis: P| parallel(axis, a) && line_distance(origin, o, a) <= tol;
            if let Some((p, n)) = face.plane {
                if perpendicular_to(n, a) && dot(sub(p, o), n).abs() <= tol {
                    return Ok(Role::Cap);
                }
                if parallel(n, a) {
                    return Ok(Role::Wall);
                }
                return Err(format!("planar face {name} neither contains nor crosses the axis"));
            }
            match face.surface {
                Some(Surface::Cylinder { origin, axis, .. })
                | Some(Surface::Cone { apex: origin, axis, .. })
                | Some(Surface::Torus { center: origin, axis, .. })
                    if coaxial(origin, axis) =>
                {
                    Ok(Role::Wall)
                }
                Some(Surface::Sphere { center, .. }) if line_distance(center, o, a) <= tol => Ok(Role::Wall),
                Some(_) => Err(format!("face {name} is not a surface of revolution about the axis")),
                None => Err(format!("free-form face {name}: a revolution needs analytic faces")),
            }
        }
    }
}

/// Candidate extrusion directions: the normals of planar faces and the axes of cylinders, most
/// common first, over the faces of every body given.
pub(crate) fn translations(breps: &[&Brep]) -> Vec<Frame> {
    let mut dirs: Vec<(P, f64)> = Vec::new();
    let mut push = |d: P, weight: f64| {
        let Some(mut d) = unit(d) else { return };
        // A canonical sign: the largest component positive.
        let k = (0..3).max_by(|&i, &j| d[i].abs().total_cmp(&d[j].abs()).then(j.cmp(&i))).unwrap_or(0);
        if d.get(k).copied().unwrap_or(0.0) < 0.0 {
            d = scale(d, -1.0);
        }
        match dirs.iter_mut().find(|(e, _)| parallel(*e, d)) {
            Some(slot) => slot.1 += weight,
            None => dirs.push((d, weight)),
        }
    };
    for f in breps.iter().flat_map(|b| &b.faces) {
        if let Some((_, n)) = f.plane {
            push(n, f.area);
        } else if let Some(Surface::Cylinder { axis, .. }) = f.surface {
            push(axis, 0.0);
        }
    }
    dirs.sort_by(|a, b| b.1.total_cmp(&a.1));
    let Some(o) = breps.first().map(|b| b.point(0)) else { return Vec::new() };
    dirs.into_iter().filter_map(|(d, _)| Frame::translate(o, d)).collect()
}

/// Candidate axes of revolution: the axes of the analytic surfaces of revolution of every body
/// given.
pub(crate) fn rotations(breps: &[&Brep]) -> Vec<Frame> {
    let mut axes: Vec<(P, P)> = Vec::new();
    let diagonal = breps.iter().map(|b| b.diagonal).fold(0.0, f64::max);
    for f in breps.iter().flat_map(|b| &b.faces) {
        let (o, a) = match f.surface {
            Some(Surface::Cylinder { origin, axis, .. })
            | Some(Surface::Cone { apex: origin, axis, .. })
            | Some(Surface::Torus { center: origin, axis, .. }) => (origin, axis),
            _ => continue,
        };
        if !axes.iter().any(|(p, b)| parallel(*b, a) && line_distance(o, *p, *b) <= RELATIVE_AXIS * diagonal) {
            axes.push((o, a));
        }
    }
    axes.into_iter().filter_map(|(o, a)| Frame::rotate(o, a)).collect()
}

/// The frame's classification of the bodies swept together: face roles (per body) and levels.
pub(crate) struct Setup {
    frame: Frame,
    roles: Vec<Vec<Role>>,
    /// Sorted levels (heights or angles in `[0, 2π)`).
    levels: Vec<f64>,
    /// Sweep coordinates closer than this are one level.
    h_tol: f64,
}

/// Classify faces and edges and find the levels of the bodies (one context each, sharing one
/// tolerance) swept together; an `Unsupported` error says why a body does not sweep in this
/// frame, naming it when there are several.
pub(crate) fn setup(ctxs: &[Ctx], frame: &Frame) -> Result<Setup> {
    let tol = ctxs.first().map_or(0.0, |c| c.tol);
    let r_max = ctxs.iter().flat_map(|c| c.brep.surface.points()).map(|&p| frame.coords(p).0[1].abs()).fold(0.0, f64::max);
    let h_tol = match frame {
        Frame::Translate { .. } => tol,
        Frame::Rotate { .. } => tol / r_max.max(tol),
    };
    let (mut roles, mut heights) = (Vec::with_capacity(ctxs.len()), Vec::new());
    for (c, ctx) in ctxs.iter().enumerate() {
        let named = |e: Error| match e {
            Error::Unsupported(m) if ctxs.len() > 1 => Error::Unsupported(format!("cell {c}: {m}")),
            e => e,
        };
        roles.push(classify(ctx, frame, tol, &mut heights).map_err(named)?);
    }
    heights.sort_by(f64::total_cmp);
    let mut levels: Vec<f64> = Vec::new();
    for h in heights {
        if levels.last().is_none_or(|&l| h - l > h_tol) {
            levels.push(h);
        }
    }
    if frame.periodic()
        && levels.len() > 1
        && let (Some(&first), Some(&last)) = (levels.first(), levels.last())
        && first + TAU - last <= h_tol
    {
        levels.pop();
    }
    let need = if frame.periodic() { 1 } else { 2 };
    if levels.len() < need {
        return Err(Error::Unsupported("no source and target faces across the sweep".into()));
    }
    Ok(Setup { frame: *frame, roles, levels, h_tol })
}

/// One body's face roles in the frame; the heights of its edges across the sweep and of its caps
/// are added to `heights`.
fn classify(ctx: &Ctx, frame: &Frame, tol: f64, heights: &mut Vec<f64>) -> Result<Vec<Role>> {
    let brep = ctx.brep;
    let mut roles = Vec::with_capacity(brep.faces.len());
    for f in &brep.faces {
        roles.push(role(brep, f, frame, tol).map_err(Error::Unsupported)?);
    }
    if frame.periodic() {
        // A revolved region touching the axis would collapse its elements there.
        if let Some(p) = brep.surface.points().iter().find(|&&p| frame.coords(p).0[1] <= tol) {
            return Err(Error::Measured(
                crate::Kind::Unsupported,
                crate::message![
                    "the body reaches the axis of revolution (at ",
                    crate::message::point(*p),
                    "); revolved hexahedra would collapse there"
                ],
            ));
        }
    }
    for (e, edge) in brep.edges.iter().enumerate() {
        let pts: Vec<(Q, f64)> = edge.points.iter().map(|&v| frame.coords(brep.point(v))).collect();
        let Some(&(q0, h0)) = pts.first() else { continue };
        let across = pts.iter().all(|&(q, h)| frame.span(q, frame.dh(h, h0)) <= tol);
        let along = pts.iter().all(|&(q, _)| dist([q[0], q[1], 0.0], [q0[0], q0[1], 0.0]) <= tol);
        if across {
            heights.push(h0);
        } else if !along {
            let [f0, f1] = edge.faces.map(|f| brep.face_tag(f));
            return Err(Error::Unsupported(format!("edge {e} between faces {f0} and {f1} runs neither along nor across the sweep")));
        }
    }
    for (f, r) in brep.faces.iter().zip(&roles) {
        if *r == Role::Cap
            && let Some(&v) = f.vertices.first()
        {
            heights.push(frame.coords(brep.point(v)).1);
        }
    }
    Ok(roles)
}

impl Setup {
    /// The number of levels the frame sweeps through (for ordering candidates).
    pub(crate) fn level_count(&self) -> usize {
        self.levels.len()
    }
}

/// A size field on the cross-section plane at sweep coordinate `h`.
pub(crate) struct PlaneField<'a> {
    pub field: &'a dyn SizeField,
    pub frame: Frame,
    pub h: f64,
}

impl SizeField for PlaneField<'_> {
    fn size(&self, p: [f64; 3]) -> f64 {
        self.field.size(self.frame.point([p[0], p[1]], self.h))
    }
}

/// The cross-section of a frame the bodies sweep in: the arrangement (cut into blocks),
/// the level intervals and which body fills each region in each of them. It does not depend on
/// the element size, so refinement retries reuse it.
pub(crate) struct Prepared {
    setup: Setup,
    arr: Arrangement,
    intervals: Vec<(f64, f64)>,
    filled: Vec<Vec<Option<u32>>>,
    used: Vec<bool>,
    /// The first arrangement tag of each body's faces (curve tags number the faces of all bodies
    /// one after the other).
    offsets: Vec<usize>,
}

impl Prepared {
    /// The body and face of an arrangement tag.
    fn face(&self, ctxs: &[Ctx], tag: usize) -> Option<(usize, usize)> {
        let c = self.offsets.iter().rposition(|&o| o <= tag)?;
        let f = tag - self.offsets.get(c)?;
        ctxs.get(c).filter(|ctx| f < ctx.brep.faces.len()).map(|_| (c, f))
    }
}

/// Lay out the cross-section of `setup`'s frame, for the bodies of `ctxs` (those `setup` was
/// made for) swept together: the walls of every body are imprinted on each other, so where two
/// bodies meet on a wall they get the same curves, and each region of each level interval is
/// filled by the body it lies in.
pub(crate) fn prepare(ctxs: &[Ctx], setup: Setup) -> Result<Prepared> {
    let ctx = ctxs.first().ok_or_else(|| Error::Failed("no body to sweep".into()))?;
    let tol = ctx.tol;
    let frame = &setup.frame;
    let levels = &setup.levels;
    // Lower outlines of the walls, and every vertex, of every body.
    let mut polylines: Vec<(Vec<Q>, usize)> = Vec::new();
    let mut vertices: Vec<Q> = Vec::new();
    let mut offsets = Vec::with_capacity(ctxs.len());
    let mut offset = 0;
    for (c, body) in ctxs.iter().enumerate() {
        let brep = body.brep;
        let roles = setup.roles.get(c).map(Vec::as_slice).unwrap_or(&[]);
        offsets.push(offset);
        for (e, edge) in brep.edges.iter().enumerate() {
            let pts: Vec<(Q, f64)> = edge.points.iter().map(|&v| frame.coords(brep.point(v))).collect();
            let Some(&(_, h0)) = pts.first() else { continue };
            if !pts.iter().all(|&(q, h)| frame.span(q, frame.dh(h, h0)) <= tol) {
                continue;
            }
            for &f in &edge.faces {
                if roles.get(f) != Some(&Role::Wall) {
                    continue;
                }
                // The wall's triangle on the edge's first segment: is its third vertex above?
                let Some(v) = brep.third_vertex(e, f) else { continue };
                let (qc, hc) = frame.coords(brep.point(v));
                if frame.dh(hc, h0) > 0.0 && frame.span(qc, frame.dh(hc, h0)) > tol {
                    polylines.push((pts.iter().map(|x| x.0).collect(), offset + f));
                }
            }
        }
        // Curve ends: the B-rep vertices and the ends of open edges. A closed edge without a
        // vertex on it (a hole's rim after its halves are merged, #512) starts wherever its chain
        // was picked up: that point is no end, and the cross-section's mesher places its nodes.
        vertices.extend(
            brep.edges
                .iter()
                .filter(|e| !e.closed())
                .flat_map(|e| [e.first(), e.last()])
                .chain(brep.is_vertex.iter().enumerate().filter(|(_, v)| **v).map(|(i, _)| i as u32))
                .map(|v| frame.coords(brep.point(v)).0),
        );
        offset += brep.faces.len();
    }
    let coincide = 2.0 * ctx.opts.chordal_tolerance.max(0.0) + 4.0 * tol;
    // Level intervals; a revolution's last one wraps round to the first level.
    let (Some(&first_level), Some(rest)) = (levels.first(), levels.get(1..)) else {
        return Err(Error::Unsupported("no levels".into()));
    };
    let ends: Vec<f64> = if frame.periodic() { rest.iter().copied().chain(std::iter::once(first_level + TAU)).collect() } else { rest.to_vec() };
    let intervals: Vec<(f64, f64)> = levels.iter().copied().zip(ends).collect();
    let mut inside = Occupants::new(ctxs, tol);
    let mut arr = arrangement::build(&vertices, &polylines, tol, coincide)?;
    let (mut filled, mut used) = occupancy(&arr, frame, &intervals, &mut inside)?;
    // Cut the filled regions into blocks the structured schemes fit, and lay the arrangement out
    // again with the cuts.
    let h_ref = levels.first().copied().unwrap_or(0.0);
    let plane = PlaneField { field: ctx.field, frame: *frame, h: h_ref };
    let settings = decompose::Settings { field: &plane, tol, opts: &ctx.opts.quad };
    let locked = seeded_curves(ctxs, frame, h_ref, &intervals, setup.h_tol, &arr);
    let cut_lines = decompose::cuts(&arr, &used, &locked, &settings)?;
    if !cut_lines.is_empty() {
        polylines.extend(cut_lines.into_iter().map(|c| (c, decompose::CUT)));
        arr = arrangement::build(&vertices, &polylines, tol, coincide)?;
        (filled, used) = occupancy(&arr, frame, &intervals, &mut inside)?;
    }
    Ok(Prepared { setup, arr, intervals, filled, used, offsets })
}

/// Mesh the prepared cross-section, with its curves `2^refine` times finer than the size asks;
/// with the body (index into `ctxs`) of each hexahedron.
pub(crate) fn mesh(ctxs: &[Ctx], prepared: &Prepared, refine: u32) -> Result<(Raw, Vec<u32>)> {
    let ctx = ctxs.first().ok_or_else(|| Error::Failed("no body to sweep".into()))?;
    let Prepared { setup, arr, intervals, filled, used, .. } = prepared;
    let (frame, levels) = (&setup.frame, &setup.levels);
    let h_ref = levels.first().copied().unwrap_or(0.0);
    // Interval assignment over the curves of filled regions, and the layer counts.
    let mut ia = Intervals::default();
    let plane_field = PlaneField { field: ctx.field, frame: *frame, h: h_ref };
    let curve3d = |c: usize| -> Vec<P> {
        arr.curves
            .get(c)
            .map(|cv| cv.points.iter().map(|&p| frame.point(arr.points.get(p).copied().unwrap_or([0.0; 2]), h_ref)).collect())
            .unwrap_or_default()
    };
    let fine = f64::from(1u32 << refine.min(8));
    for c in 0..arr.curves.len() {
        let goal = ideal_divisions(ctx.field, &curve3d(c))?;
        ia.var(goal * fine, 1);
    }
    let schemes = section::schemes(arr, used, &mut ia);
    let seeds = Seeds::of(ctxs, frame, h_ref, intervals, setup.h_tol, &(0..arr.curves.len()).map(curve3d).collect::<Vec<_>>())?;
    for (c, seed) in seeds.curves.iter().enumerate() {
        if let Some((seed, _)) = seed {
            ia.seed(c, seed);
        }
    }
    let counts = ia.solve()?;
    // A biased seed's spacing carries across the cross-section's mapped regions to the opposite
    // curves without seeds of their own (#257).
    let ties: Vec<[crate::seeds::Use; 2]> = arr
        .regions
        .iter()
        .zip(&schemes)
        .filter_map(|(region, &scheme)| region.loops.first().map(|lp| region::opposite_ties(lp, scheme)))
        .flatten()
        .collect();
    let carried = crate::seeds::carry(&seeds.curves, |s: &&EdgeSeed| s.bias.is_some(), &ties);
    let spacing: Vec<Option<(&EdgeSeed, bool)>> = seeds.curves.iter().zip(carried).map(|(own, carried)| own.or(carried)).collect();
    let r_max = arr.points.iter().map(|q| q[1].abs()).fold(0.0, f64::max);
    let mut layers = Vec::with_capacity(intervals.len());
    for &(lo, hi) in intervals.iter() {
        let path: Vec<P> = match frame {
            Frame::Translate { .. } => vec![frame.point([0.0, 0.0], lo), frame.point([0.0, 0.0], hi)],
            Frame::Rotate { .. } => (0..=16).map(|k| frame.point([0.0, r_max], lo + (hi - lo) * k as f64 / 16.0)).collect(),
        };
        let mut n = edge_divisions(ctx.field, &path)?;
        if frame.periodic() {
            n = n.max(((hi - lo) / FRAC_PI_4).ceil() as usize);
        }
        if let Some(Some(l)) = seeds.layers.get(layers.len()) {
            let seeded = l.fractions.len().saturating_sub(1).max(1);
            n = match l.seed.constraint {
                SeedConstraint::Exact | SeedConstraint::Free => seeded,
                SeedConstraint::Min => n.max(seeded),
            };
        }
        layers.push(n);
    }
    // The cross-section's nodes (curve interiors projected onto analytic walls) and quadrilaterals.
    let place = |c: usize, n: usize| -> Result<Vec<Q>> {
        let pts3 = match spacing.get(c).copied().flatten() {
            Some((seed, rev)) => crate::seeds::place(&curve3d(c), &seed.fractions(n, rev))?,
            None => divide_polyline(ctx.field, &curve3d(c), n)?,
        };
        let wall = arr
            .curves
            .get(c)
            .and_then(|curve| {
                curve.tags.iter().find_map(|&t| {
                    let (c, f) = prepared.face(ctxs, t)?;
                    ctxs.get(c)?.brep.faces.get(f)?.surface
                })
            })
            .filter(|s| !matches!(s, Surface::Plane { .. }));
        Ok(pts3.iter().map(|p| frame.coords(wall.and_then(|s| s.project(*p)).unwrap_or(*p)).0).collect())
    };
    let section::Section { nodes: nodes2, quads: region_quads, schemes } =
        section::mesh(arr, used, schemes, &counts, &place, &plane_field, &ctx.opts.quad)?;
    // Layer coordinates.
    let mut heights: Vec<f64> = Vec::new();
    let mut first_layer: Vec<usize> = Vec::with_capacity(intervals.len());
    for (k, &(lo, hi)) in intervals.iter().enumerate() {
        first_layer.push(heights.len());
        let n = layers.get(k).copied().unwrap_or(1);
        let fractions = match seeds.layers.get(k).and_then(Option::as_ref) {
            Some(l) if l.fractions.len() == n + 1 => l.fractions.clone(),
            _ => (0..=n).map(|j| j as f64 / n as f64).collect(),
        };
        for f in fractions.iter().take(n) {
            heights.push(lo + (hi - lo) * f);
        }
    }
    if !frame.periodic()
        && let Some(&(_, hi)) = intervals.last()
    {
        heights.push(hi);
    }
    let total_layers = heights.len();
    // Element budget.
    let mut count = 0usize;
    for (r, quads) in region_quads.iter().enumerate() {
        for (k, f) in filled.get(r).map(Vec::as_slice).unwrap_or(&[]).iter().enumerate() {
            if f.is_some() {
                count = count.saturating_add(quads.len().saturating_mul(layers.get(k).copied().unwrap_or(0)));
            }
        }
    }
    if count > ctx.opts.max_elements {
        return Err(Error::Limit(format!("the mesh would have {count} hexahedra; the limit is {}", ctx.opts.max_elements)));
    }
    // Hexahedra, with nodes made as first used.
    let mut nodes: Vec<P> = Vec::new();
    let mut node_at: HashMap<(u32, usize), u32> = HashMap::new();
    let mut node = |n2: u32, layer: usize, nodes: &mut Vec<P>| -> u32 {
        let layer = if frame.periodic() { layer % total_layers.max(1) } else { layer };
        *node_at.entry((n2, layer)).or_insert_with(|| {
            let q = nodes2.get(n2 as usize).copied().unwrap_or([0.0; 2]);
            nodes.push(frame.point(q, heights.get(layer).copied().unwrap_or(0.0)));
            nodes.len() as u32 - 1
        })
    };
    let mut hexes: Vec<[u32; 8]> = Vec::with_capacity(count);
    let mut cells: Vec<u32> = Vec::with_capacity(count);
    for (r, quads) in region_quads.iter().enumerate() {
        for (k, f) in filled.get(r).map(Vec::as_slice).unwrap_or(&[]).iter().enumerate() {
            let Some(cell) = *f else { continue };
            let (s, n) = (first_layer.get(k).copied().unwrap_or(0), layers.get(k).copied().unwrap_or(0));
            for q in quads {
                for j in s..s + n {
                    let b = q.map(|v| node(v, j, &mut nodes));
                    let t = q.map(|v| node(v, j + 1, &mut nodes));
                    hexes.push([b[0], b[1], b[2], b[3], t[0], t[1], t[2], t[3]]);
                    cells.push(cell);
                }
            }
        }
    }
    let filled_regions = used.iter().filter(|&&u| u).count();
    let structured = schemes.iter().zip(used.iter()).filter(|(s, u)| **u && **s != Scheme::Free).count();
    let bodies = if ctxs.len() > 1 { format!(" ({} cells together)", ctxs.len()) } else { String::new() };
    let description = frame.describe().then(message![format!(
        "{bodies}: {} level{}, {} region{} ({} structured, {} free), {} layer{}",
        levels.len(),
        if levels.len() == 1 { "" } else { "s" },
        filled_regions,
        if filled_regions == 1 { "" } else { "s" },
        structured,
        filled_regions - structured,
        total_layers.saturating_sub(usize::from(!frame.periodic())),
        if total_layers == 2 { "" } else { "s" },
    )]);
    Ok((Raw { nodes, hexes, description }, cells))
}

/// The edge seeds as a layered mesh takes them: on the cross-section's curves (seeds of edges
/// across the sweep, projected into the cross-section) and on the level intervals (seeds of
/// edges along the sweep, which set the interval's layers). Each with whether it runs against
/// the curve, or from the interval's top down.
struct Seeds<'s> {
    curves: Vec<Option<(&'s EdgeSeed, bool)>>,
    layers: Vec<Option<LayerSeed<'s>>>,
}

/// A seed's layers in one level interval: its nodes there, as fractions from the interval's
/// bottom (0) to its top (1).
#[derive(Clone)]
struct LayerSeed<'s> {
    seed: &'s EdgeSeed,
    fractions: Vec<f64>,
}

/// Split a seed's node positions (fractions from the bottom of the span it seeds, 0 to 1) over
/// the level intervals in it (`bounds`: the levels from the bottom of the span to its top, as
/// fractions): each inner level takes the nearest node still free, and the nodes between two
/// levels are stretched to fit between them. Every interval gets at least one division; `None`
/// when there are fewer nodes than that needs.
fn split_over(nodes: &[f64], bounds: &[f64]) -> Option<Vec<Vec<f64>>> {
    let (n, m) = (nodes.len().checked_sub(1)?, bounds.len().checked_sub(1)?);
    if m == 0 || n < m {
        return None;
    }
    // The node index at each level.
    let mut at = vec![0usize; m + 1];
    at[m] = n;
    for i in 1..m {
        let b = bounds.get(i).copied()?;
        let (lo, hi) = (at[i - 1] + 1, n - (m - i));
        at[i] = (lo..=hi).min_by(|&a, &c| (nodes[a] - b).abs().total_cmp(&(nodes[c] - b).abs()))?;
    }
    (0..m)
        .map(|i| {
            let (j0, j1) = (at[i], at[i + 1]);
            let (z0, z1) = (nodes.get(j0).copied()?, nodes.get(j1).copied()?);
            let span = z1 - z0;
            Some(
                (j0..=j1)
                    .map(|j| {
                        if j == j0 {
                            0.0
                        } else if j == j1 {
                            1.0
                        } else {
                            nodes.get(j).map_or(0.0, |z| (z - z0) / span)
                        }
                    })
                    .collect(),
            )
        })
        .collect()
}

impl<'s> Seeds<'s> {
    fn of(ctxs: &'s [Ctx], frame: &Frame, h_ref: f64, intervals: &[(f64, f64)], h_tol: f64, curves: &[Vec<P>]) -> Result<Seeds<'s>> {
        let mut out = Seeds { curves: vec![None; curves.len()], layers: vec![None; intervals.len()] };
        for ctx in ctxs {
            out.add(ctx, frame, h_ref, intervals, h_tol, curves)?;
        }
        Ok(out)
    }

    /// Add the seeds on the edges of one body.
    fn add(&mut self, ctx: &'s Ctx, frame: &Frame, h_ref: f64, intervals: &[(f64, f64)], h_tol: f64, curves: &[Vec<P>]) -> Result<()> {
        let brep = ctx.brep;
        let tol = ctx.tol;
        let out = self;
        let stol = super::seed_tolerance(ctx.opts, tol);
        let conflict = |a: &EdgeSeed, b: &EdgeSeed, why: &str| {
            Error::Invalid(format!(
                "the exact edge seeds {} and {} conflict: {why}; change a count or make a seed `min` or `free`",
                a.describe(),
                b.describe()
            ))
        };
        let unsupported = |seed: &EdgeSeed, why: &str| Error::Unsupported(format!("the exact edge seed {} {why}", seed.describe()));
        for e in 0..brep.edges.len() {
            let Some((seed, rev)) = ctx.seed_of(e) else { continue };
            let exact = seed.constraint == SeedConstraint::Exact;
            let pts: Vec<(Q, f64)> = brep.edge_points(e, true).into_iter().map(|p| frame.coords(p)).collect();
            let (Some(&(_, h0)), Some(&(_, h1))) = (pts.first(), pts.last()) else { continue };
            let across = pts.iter().all(|&(q, h)| frame.span(q, frame.dh(h, h0)) <= tol);
            if across {
                // The seed moved into the cross-section plane at the reference level.
                let projected = EdgeSeed { points: seed.points.iter().map(|&p| frame.point(frame.coords(p).0, h_ref)).collect(), ..seed.clone() };
                match curves.iter().position(|c| crate::seeds::matches(&projected, c, stol).is_some()) {
                    Some(c) => {
                        let r = crate::seeds::matches(&projected, curves.get(c).map(Vec::as_slice).unwrap_or(&[]), stol).unwrap_or(false);
                        let slot = out.curves.get_mut(c).ok_or_else(|| Error::Failed("curve index".into()))?;
                        if let Some((other, _)) = slot
                            && exact
                            && other.constraint == SeedConstraint::Exact
                            && other.count != seed.count
                        {
                            return Err(conflict(other, seed, "both lie on one curve of the swept cross-section"));
                        }
                        *slot = Some((seed, r));
                    }
                    None if exact => return Err(unsupported(seed, "is not one curve of the swept cross-section")),
                    None => {}
                }
                continue;
            }
            // Along the sweep: the interval from its lower to its upper end.
            let span = frame.dh(h1, h0);
            let (lo, hi, down) = if span >= 0.0 { (h0, h0 + span, false) } else { (h0 + span, h0, true) };
            let same = |a: f64, b: f64| frame.dh(a, b).abs() <= h_tol;
            let first = intervals.iter().position(|&(l, _)| same(l, lo));
            let last = intervals.iter().position(|&(_, h)| same(h, hi));
            let (Some(k0), Some(k1)) = (first, last) else {
                if exact {
                    return Err(unsupported(seed, "does not run from one level of the sweep to another"));
                }
                continue;
            };
            let spanned = intervals.get(k0..=k1).unwrap_or(&[]);
            let length = frame.dh(hi, lo).abs().max(f64::MIN_POSITIVE);
            let mut bounds: Vec<f64> = spanned.iter().map(|&(l, _)| frame.dh(l, lo) / length).collect();
            bounds.push(1.0);
            // The seed's nodes from the bottom: the seed runs down when the edge does and the
            // seed runs with it, or up against a downward edge.
            let nodes = seed.fractions(seed.count, rev != down);
            let Some(parts) = split_over(&nodes, &bounds) else {
                if exact {
                    return Err(unsupported(seed, &format!("has fewer divisions than the {} level intervals it spans", spanned.len())));
                }
                continue;
            };
            for (k, fractions) in (k0..=k1).zip(parts) {
                let slot = out.layers.get_mut(k).ok_or_else(|| Error::Failed("interval index".into()))?;
                if let Some(other) = slot
                    && exact
                    && other.seed.constraint == SeedConstraint::Exact
                    && other.fractions.len() != fractions.len()
                {
                    return Err(conflict(other.seed, seed, "both set the layers between the same two levels"));
                }
                *slot = Some(LayerSeed { seed, fractions });
            }
        }
        Ok(())
    }
}

/// The bodies swept together, for inside tests: each body's test and bounding box.
struct Occupants<'a> {
    bodies: Vec<(Inside<'a>, P, P)>,
}

impl<'a> Occupants<'a> {
    fn new(ctxs: &[Ctx<'a>], tol: f64) -> Self {
        let bodies = ctxs
            .iter()
            .map(|c| {
                let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
                for p in c.brep.surface.points() {
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k] - tol);
                        hi[k] = hi[k].max(p[k] + tol);
                    }
                }
                (Inside::new(c.brep.surface), lo, hi)
            })
            .collect();
        Occupants { bodies }
    }

    /// The first body `x` lies in.
    fn body_at(&mut self, x: P) -> Result<Option<u32>> {
        for (i, (inside, lo, hi)) in self.bodies.iter_mut().enumerate() {
            if (0..3).all(|k| x[k] >= lo[k] && x[k] <= hi[k]) && inside.contains(x)? {
                return Ok(u32::try_from(i).ok());
            }
        }
        Ok(None)
    }
}

/// Which body fills each region of the arrangement in each level interval (inside tests at the
/// region's centre, mid-interval), and which regions any body fills.
fn occupancy(arr: &Arrangement, frame: &Frame, intervals: &[(f64, f64)], inside: &mut Occupants) -> Result<(Vec<Vec<Option<u32>>>, Vec<bool>)> {
    if arr.regions.is_empty() {
        return Err(Error::Unsupported("the projected outline encloses no region".into()));
    }
    let mut filled: Vec<Vec<Option<u32>>> = Vec::with_capacity(arr.regions.len());
    for region in &arr.regions {
        let c = interior_point(arr, region)?;
        let mut row = Vec::with_capacity(intervals.len());
        for &(lo, hi) in intervals {
            row.push(inside.body_at(frame.point(c, 0.5 * (lo + hi)))?);
        }
        filled.push(row);
    }
    let used: Vec<bool> = filled.iter().map(|r| r.iter().any(Option::is_some)).collect();
    if !used.iter().any(|&u| u) {
        return Err(Error::Unsupported("no region of the cross-section lies inside the body".into()));
    }
    Ok((filled, used))
}

/// Which curves of `arr` carry an edge seed (a cut may not end inside them).
fn seeded_curves(ctxs: &[Ctx], frame: &Frame, h_ref: f64, intervals: &[(f64, f64)], h_tol: f64, arr: &Arrangement) -> Vec<bool> {
    let curves: Vec<Vec<P>> =
        arr.curves.iter().map(|cv| cv.points.iter().map(|&p| frame.point(arr.points.get(p).copied().unwrap_or([0.0; 2]), h_ref)).collect()).collect();
    match Seeds::of(ctxs, frame, h_ref, intervals, h_tol, &curves) {
        Ok(seeds) => seeds.curves.iter().map(Option::is_some).collect(),
        Err(_) => vec![false; arr.curves.len()],
    }
}

/// A point well inside a region: off the middle of one of its outer edges, inside the outer
/// loop, outside the holes and clear of every boundary segment.
pub(crate) fn interior_point(arr: &Arrangement, region: &arrangement::Region) -> Result<Q> {
    let at = |p: usize| arr.points.get(p).copied().unwrap_or([0.0; 2]);
    let polygon = |lp: &[(usize, bool)]| -> Vec<Q> {
        let mut out = Vec::new();
        for &(c, fwd) in lp {
            let mut v: Vec<usize> = arr.curves.get(c).map(|x| x.points.clone()).unwrap_or_default();
            if !fwd {
                v.reverse();
            }
            v.pop();
            out.extend(v.into_iter().map(at));
        }
        out
    };
    let loops: Vec<Vec<Q>> = region.loops.iter().map(|l| polygon(l)).collect();
    // Crossing-number winding with exact side tests (a probe point off the plane's domain
    // counts as outside).
    let winding = |p: Q, poly: &[Q]| -> i64 {
        let mut wn = 0;
        for k in 0..poly.len() {
            let (a, b) = (poly[k], poly[(k + 1) % poly.len()]);
            let side = arrangement::sign(a, b, p).unwrap_or(Sign::Zero);
            if a[1] <= p[1] {
                if b[1] > p[1] && side == Sign::Positive {
                    wn += 1;
                }
            } else if b[1] <= p[1] && side == Sign::Negative {
                wn -= 1;
            }
        }
        wn
    };
    let clearance = |p: Q| -> f64 {
        loops
            .iter()
            .flat_map(|poly| (0..poly.len()).map(move |k| (poly[k], poly[(k + 1) % poly.len()])))
            .map(|(a, b)| arrangement::seg_dist2(p, a, b).0.sqrt())
            .fold(f64::INFINITY, f64::min)
    };
    let outer = loops.first().ok_or_else(|| Error::Failed("region without a loop".into()))?;
    let mut edges: Vec<(f64, usize)> = (0..outer.len()).map(|k| (arrangement::d2(outer[k], outer[(k + 1) % outer.len()]), k)).collect();
    edges.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut best: Option<(f64, Q)> = None;
    for &(len2, k) in edges.iter().take(64) {
        let (a, b) = (outer[k], outer[(k + 1) % outer.len()]);
        let len = len2.sqrt();
        let (m, n) = ([(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5], [-(b[1] - a[1]) / len, (b[0] - a[0]) / len]);
        for f in [0.5, 0.25, 0.1, 0.03, 0.01, 1e-3] {
            let p = [m[0] + n[0] * f * len, m[1] + n[1] * f * len];
            if winding(p, outer) == 0 || loops.iter().skip(1).any(|h| winding(p, h) != 0) {
                continue;
            }
            let c = clearance(p);
            if best.is_none_or(|(bc, _)| c > bc) {
                best = Some((c, p));
            }
        }
    }
    best.filter(|(c, _)| *c > 0.0).map(|(_, p)| p).ok_or_else(|| Error::Failed("no point found inside a cross-section region".into()))
}

#[cfg(test)]
mod seed_tests {
    use super::split_over;

    #[test]
    fn a_seed_over_several_levels_is_split_at_them() {
        // Six even divisions over levels at 0, 0.4 and 1: the level takes the node at 1/3.
        let nodes: Vec<f64> = (0..=6).map(|j| j as f64 / 6.0).collect();
        let parts = split_over(&nodes, &[0.0, 0.4, 1.0]).unwrap();
        assert_eq!(parts.iter().map(|p| p.len() - 1).collect::<Vec<_>>(), vec![2, 4]);
        assert!(parts.iter().all(|p| p[0] == 0.0 && *p.last().unwrap() == 1.0 && p.windows(2).all(|w| w[0] < w[1])));
        // Every interval keeps at least one division; too few nodes is None.
        let parts = split_over(&[0.0, 0.5, 1.0], &[0.0, 0.1, 0.2, 1.0]);
        assert!(parts.is_none());
        let parts = split_over(&[0.0, 0.9, 0.95, 1.0], &[0.0, 0.1, 0.2, 1.0]).unwrap();
        assert_eq!(parts.iter().map(|p| p.len() - 1).collect::<Vec<_>>(), vec![1, 1, 1]);
    }
}
