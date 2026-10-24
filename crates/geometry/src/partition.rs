//! Partitions: a design's bodies cut into cells by planes, for meshing cell by cell.
//!
//! A partition cuts bodies of the open design by a plane with SolidCraft's kernel (the body's
//! intersections with a box on each side of the plane, as SolidCraft's split body does), in
//! memory: the design file is never changed. The partitions live in the FieldCraft model and are applied again, in order,
//! whenever the design is opened, so an edit of the design followed by a re-open replays them.
//!
//! Names are FieldCraft's own, so they do not depend on the order of the design's timeline:
//! - **Cells.** A body cut by a plane gives `<body>/a` (the piece on the plane's negative side)
//!   and `<body>/b`; several disjoint pieces on one side are `<body>/a1`, `<body>/a2`, … ordered
//!   by centroid.
//! - **Cut faces** pair across the cut: on a cell on the negative side (its outward normal along
//!   the plane's normal) the face is `partition:<name>:a`, on the positive side
//!   `partition:<name>:b`. Several disjoint faces on one cell (a cut through a hole) are `…#0`,
//!   `…#1`, ordered along the plane, which a set naming `partition:<name>:a` selects together.
//! - **Every other face** inherits the persistent name of the face it is a piece of (found on the
//!   same surface under a point of the piece). Where a cell has several pieces of one face they
//!   are `<name>#0`, `<name>#1`, …, which sets naming `<name>` select.

use solidcraft_geom::{Plane, Vec3};
use solidcraft_kernel::guard;

use crate::{Design, Error, Face, FaceSurface, FeatureError, Result, Shape, SurfaceMesh};

/// The largest coordinate a cut plane may have.
const MAX_COORD: f64 = 1e9; // unit-factors: a cap on hostile input
/// A cut plane's normal must be longer than this.
const MIN_NORMAL: f64 = 1e-12; // unit-factors: the length of a direction vector
/// Unit directions this close to parallel are one.
const SAME_DIRECTION: f64 = 1e-9; // unit-factors: a component of a unit vector
/// Unit directions whose cross product is shorter than this are parallel.
const PARALLEL_SINE: f64 = 1e-6; // unit-factors: the sine of an angle
/// Cone angles this close are one.
const SAME_ANGLE: f64 = 1e-9; // unit-factors: an angle in radians
/// A face lies on the cut within this fraction of the tessellation tolerance.
const ON_FACE: f64 = 1e-3; // unit-factors: relative to the tolerance

/// Most partitions a design takes.
pub const MAX_CUTS: usize = 256;
/// Most cells the partitions make in all.
pub const MAX_CELLS: usize = 4096;
/// A cut face lies on the plane within this fraction of its body's size.
const ON_PLANE: f64 = 1e-6; // unit-factors: relative to the body's size
/// Tessellation tolerance (fraction of a body's size) for matching and ordering faces.
const MATCH_TOLERANCE: f64 = 2e-3;

/// One partition: the bodies it cuts by a plane.
#[derive(Debug, Clone, PartialEq)]
pub struct Cut {
    /// Names the cut faces (`partition:<name>:a|b`).
    pub name: String,
    /// A point on the plane and its normal (any length but zero).
    pub point: [f64; 3],
    pub normal: [f64; 3],
    /// The bodies to cut, by name, as they are when this partition comes (the cells of earlier
    /// partitions included); a body earlier partitions cut stands for its cells (`B` for `B/a`,
    /// `B/b/a`, …); empty: every body the plane divides.
    pub bodies: Vec<String>,
}

/// What a partition did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutReport {
    pub name: String,
    /// The bodies it cut, each into two or more cells.
    pub split: Vec<String>,
    /// How many cut faces it named.
    pub faces: usize,
    /// What may make its cells hard to mesh: the plane lying on a face of a body it cut.
    pub warnings: Vec<String>,
}

fn vec3(p: [f64; 3]) -> Vec3 {
    Vec3 { x: p[0], y: p[1], z: p[2] }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn lexical(a: &[f64; 3], b: &[f64; 3]) -> std::cmp::Ordering {
    a.iter().zip(b).map(|(x, y)| x.total_cmp(y)).find(|o| o.is_ne()).unwrap_or(std::cmp::Ordering::Equal)
}

/// The plane of a cut, with its unit normal.
fn plane_of(cut: &Cut) -> Result<(Plane, [f64; 3])> {
    let finite = cut.point.iter().chain(&cut.normal).all(|v| v.is_finite() && v.abs() <= MAX_COORD);
    let len = dot(cut.normal, cut.normal).sqrt();
    if !finite || !(len > MIN_NORMAL) {
        return Err(Error::Invalid(format!("partition `{}`: the plane needs a finite point and a non-zero normal", cut.name)));
    }
    let n = cut.normal.map(|v| v / len);
    let plane = Plane::from_normal(vec3(cut.point), vec3(n))
        .ok_or_else(|| Error::Invalid(format!("partition `{}`: the plane needs a non-zero normal", cut.name)))?;
    Ok((plane, n))
}

impl Design {
    /// The design as opened, cut by `cuts` in order (the partitions replace any applied before).
    /// A cut that names a body the design does not have, that does not divide a body it names,
    /// or that the kernel fails to make on any body it applies to, is listed in [`Design::errors`] (feature `partition <name>`), as a feature of the
    /// design that failed to rebuild; the rest still apply.
    pub fn partitioned(&self, cuts: &[Cut]) -> Result<Design> {
        if cuts.len() > MAX_CUTS {
            return Err(Error::Invalid(format!("a design takes at most {MAX_CUTS} partitions")));
        }
        let planes: Vec<(Plane, [f64; 3])> = cuts.iter().map(plane_of).collect::<Result<_>>()?;
        let mut bodies: Vec<Shape> = self.source.as_ref().clone();
        let mut errors: Vec<FeatureError> = self.errors.iter().filter(|e| !crate::virtual_faces::derived(&e.feature)).cloned().collect();
        let mut reports = Vec::with_capacity(cuts.len());
        // Each cut is a kernel boolean per body it reaches, seconds each on a large part: the
        // job's progress counts them, and a cancelled job stops between them.
        fieldcraft_model::progress::stage("Partition cuts", Some(cuts.len() as f64), "cuts");
        for (k, (cut, (plane, n))) in cuts.iter().zip(&planes).enumerate() {
            fieldcraft_model::progress::count(k as f64);
            fieldcraft_model::progress::check()?;
            let feature = format!("partition {}", cut.name);
            // A listed name is the body of that name, or the cells earlier partitions cut it into.
            let names = |name: &str, b: &str| name == b || name.strip_prefix(b).is_some_and(|rest| rest.starts_with('/'));
            for b in &cut.bodies {
                if !bodies.iter().any(|x| names(&x.name, b)) {
                    errors.push(FeatureError { feature: feature.clone(), message: format!("no body `{b}` to cut"), bodies: vec![b.clone()] });
                }
            }
            let mut report = CutReport { name: cut.name.clone(), split: Vec::new(), faces: 0, warnings: Vec::new() };
            let mut next = Vec::with_capacity(bodies.len() + 1);
            // The listed names some body of which the plane divides.
            let mut divided = vec![false; cut.bodies.len()];
            for body in bodies {
                let listed: Vec<usize> = (0..cut.bodies.len()).filter(|&i| cut.bodies.get(i).is_some_and(|b| names(&body.name, b))).collect();
                if !(cut.bodies.is_empty() || !listed.is_empty()) {
                    next.push(body);
                    continue;
                }
                fieldcraft_model::progress::check()?;
                match split(&body, plane, *n, cut)? {
                    Split::Cells(cells, faces, warnings) => {
                        for &i in &listed {
                            if let Some(d) = divided.get_mut(i) {
                                *d = true;
                            }
                        }
                        report.split.push(body.name.clone());
                        report.faces += faces;
                        report.warnings.extend(warnings);
                        next.extend(cells);
                    }
                    Split::Misses => next.push(body),
                    // Listed or not: the plane may well divide the body, and a cut that silently
                    // leaves it whole would be mistaken for one that does not reach it (#283).
                    Split::Failed(e) => {
                        // Reported here; not again as a plane that divides nothing.
                        for &i in &listed {
                            if let Some(d) = divided.get_mut(i) {
                                *d = true;
                            }
                        }
                        errors.push(FeatureError {
                            feature: feature.clone(),
                            message: format!("the kernel could not cut `{}` ({e})", body.name),
                            bodies: vec![body.name.clone()],
                        });
                        next.push(body);
                    }
                }
                if next.len() > MAX_CELLS {
                    return Err(Error::Invalid(format!("the partitions make at most {MAX_CELLS} cells")));
                }
            }
            for (b, divided) in cut.bodies.iter().zip(divided) {
                // A body that is not there is reported above; one that is, and that the plane
                // misses (all of its cells), is an error: the cut would do nothing to it.
                if !divided && next.iter().any(|x| names(&x.name, b)) {
                    errors.push(FeatureError {
                        feature: feature.clone(),
                        message: format!("the plane does not divide `{b}`"),
                        bodies: vec![b.clone()],
                    });
                }
            }
            bodies = next;
            reports.push(report);
        }
        Ok(Design {
            name: self.name.clone(),
            bodies,
            errors,
            cuts: reports,
            source: std::sync::Arc::clone(&self.source),
            opened: std::sync::Arc::clone(&self.source),
            history: self.history.clone(),
            defeat: self.defeat.clone(),
        })
    }
}

/// A coordinate-axis half-space tool using the kernel's box primitive, as `solid.box` does.
/// Its planar face construction differs from extruding along the cutting plane's normal.
/// Only exact axis normals qualify: an almost-axis plane must keep its original orientation.
fn axis_box(plane: &Plane, n: [f64; 3], bbox: &solidcraft_geom::Aabb3, positive: bool) -> solidcraft_kernel::Result<Option<solidcraft_kernel::Body>> {
    let axis = match n {
        [x, 0.0, 0.0] if x.abs() == 1.0 => 0,
        [0.0, y, 0.0] if y.abs() == 1.0 => 1,
        [0.0, 0.0, z] if z.abs() == 1.0 => 2,
        _ => return Ok(None),
    };
    let p = crate::xyz(plane.origin);
    let c = crate::xyz(bbox.center());
    // Twice each bounding-box extent covers the body with half an extent of clearance on
    // every side. A diagonal-sized cube can be needlessly large on the thinner axes, making
    // the kernel's boolean require a tool nudge even for a well-positioned plane (#690).
    let extent = crate::xyz(bbox.size()).map(|v| v.max(LINEAR_REACH));
    let mut lo = std::array::from_fn(|i| c[i] - extent[i]);
    let mut hi = std::array::from_fn(|i| c[i] + extent[i]);
    if positive == (n[axis] > 0.0) {
        lo[axis] = p[axis];
    } else {
        hi[axis] = p[axis];
    }
    match solidcraft_kernel::box_solid(vec3(lo), vec3(hi)) {
        Ok(tool) => Ok(Some(tool)),
        // The primitive has stricter size limits than a general extrusion. Keep the existing
        // path for those tools, without enlarging either constructor's accepted domain.
        Err(solidcraft_kernel::KernelError::Invalid(_)) => Ok(None),
        Err(error) => Err(error),
    }
}

/// The parts of `body` on either side of `plane` (only the non-empty ones): its intersections
/// with two boxes, one on each side. Exact coordinate-axis planes use the body's individual
/// bounding-box dimensions; other planes keep square extrusions sized by the diagonal.
/// SolidCraft's `split_by_plane` uses boxes
/// several times the body's size, and its boolean then puts the cut up to a few micrometres off
/// the plane (1.8 µm on a 140 mm body), off the faces it should meet (#384); boxes of the body's own
/// size keep the cut on the plane.
///
/// A boolean the kernel could only make by moving one body against the other (its last-resort
/// tool nudge, micrometres and microradians) is refused: the cut would lie off the plane, next to
/// the face it was meant to run along, and leave a cell bounded by slivers no mesher takes (#427).
fn halves(body: &Shape, plane: &Plane, n: [f64; 3]) -> solidcraft_kernel::Result<Vec<solidcraft_kernel::Body>> {
    let bbox = solidcraft_kernel::measure(&body.body)?.bbox;
    let reach = bbox.diagonal().max(LINEAR_REACH);
    let c = plane.to_local(bbox.center());
    let square = solidcraft_geom::Region2 {
        outer: solidcraft_geom::Loop2::polygon(&[
            solidcraft_geom::Vec2::new(c.x - reach, c.y - reach),
            solidcraft_geom::Vec2::new(c.x + reach, c.y - reach),
            solidcraft_geom::Vec2::new(c.x + reach, c.y + reach),
            solidcraft_geom::Vec2::new(c.x - reach, c.y + reach),
        ]),
        holes: vec![],
    };
    let mut out = Vec::with_capacity(2);
    for (lo, hi) in [(0.0, reach), (-reach, 0.0)] {
        let side = match axis_box(plane, n, &bbox, lo == 0.0)? {
            Some(tool) => tool,
            None => {
                let Some(tool) = solidcraft_kernel::extrude(plane, std::slice::from_ref(&square), lo, hi)?.pop() else { continue };
                tool
            }
        };
        let (cut, report) =
            solidcraft_kernel::capture_boolean_diagnostics(|| solidcraft_kernel::boolean(&body.body, &side, solidcraft_kernel::BoolOp::Intersect));
        if let Some(moved) = report.warning() {
            return Err(solidcraft_kernel::KernelError::Failed(format!(
                "the cut is only made off the plane (SolidCraft's kernel: {moved}); the plane runs along a face or edge of the body: move it off (by a fraction of the element size)"
            )));
        }
        out.extend(cut?);
    }
    Ok(out)
}

/// Planar faces of `body` the cut's plane lies on (within [`COINCIDENT`] of the body's size),
/// as warnings: the kernel may leave micrometre slivers where a cut runs along a face (#384).
fn coincident(body: &Shape, parent: &Faces, cut: &Cut, n: [f64; 3], size: f64) -> Vec<String> {
    let reach = COINCIDENT * size.max(1.0);
    let mut out = Vec::new();
    for (i, surface) in parent.surfaces.iter().enumerate() {
        let Some(FaceSurface::Plane { normal, point }) = surface else { continue };
        let along = dot(*normal, n).abs();
        let distance = dot(sub(*point, cut.point), n).abs();
        if (along - 1.0).abs() < SAME_DIRECTION && distance <= reach {
            let name = parent.names.get(i).map_or_else(|| i.to_string(), |f| f.name.clone());
            out.push(format!(
                "partition `{}`: the plane lies on face `{name}` of `{}` ({:.1e} of its size from it); a cut along an existing face can leave micrometre slivers that no mesher takes: move the plane off the face (by a fraction of the element size)",
                cut.name,
                body.name,
                if size > 0.0 { distance / size } else { 0.0 }
            ));
        }
    }
    out
}

/// A plane within this fraction of a body's size of one of its faces lies on it.
const COINCIDENT: f64 = 1e-5;

/// The least reach of the boxes [`halves`] cuts with (mm), for bodies with no extent.
const LINEAR_REACH: f64 = 1.0;

/// What cutting one body gave.
enum Split {
    /// Its cells, named, how many cut faces they have, and warnings.
    Cells(Vec<Shape>, usize, Vec<String>),
    /// The plane does not divide it.
    Misses,
    /// The kernel failed to cut it: its message.
    Failed(String),
}

/// Cut one body.
fn split(body: &Shape, plane: &Plane, n: [f64; 3], cut: &Cut) -> Result<Split> {
    let sides = match guard("partition", || halves(body, plane, n)) {
        Ok(p) if p.len() >= 2 => p,
        Ok(_) => return Ok(Split::Misses),
        Err(e) => return Ok(Split::Failed(e.to_string())),
    };
    // Each side may fall apart into pieces that do not touch (#388): every piece is a cell.
    let mut parts = Vec::with_capacity(sides.len());
    for side in sides {
        match guard("partition", || side.lumps()) {
            Ok(pieces) if !pieces.is_empty() => parts.extend(pieces),
            Ok(_) => parts.push(side),
            Err(e) => return Ok(Split::Failed(e.to_string())),
        }
    }
    let size = body.size()?;
    let tol = (MATCH_TOLERANCE * size).max(1e-4);
    let parent = Faces::of(body, tol)?;
    let mut cells: Vec<(bool, [f64; 3], Shape)> = Vec::with_capacity(parts.len());
    for part in parts {
        let mut shape = Shape { name: String::new(), body: part, faces: Vec::new(), face_map: Vec::new(), cuts: Vec::new() };
        let m = shape.measure()?;
        let negative = dot(sub(m.centroid, cut.point), n) < 0.0;
        shape.faces = (0..m.faces).map(|_| Face { name: String::new(), exact: false, members: Vec::new() }).collect();
        cells.push((negative, m.centroid, shape));
    }
    cells.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| lexical(&a.1, &b.1)));
    let mut faces = 0;
    let count = |neg: bool| cells.iter().filter(|c| c.0 == neg).count();
    let (na, nb) = (count(true), count(false));
    let (mut ka, mut kb) = (0, 0);
    let mut out = Vec::with_capacity(cells.len());
    for (negative, _, mut shape) in cells {
        let (k, of) = if negative {
            ka += 1;
            (ka, na)
        } else {
            kb += 1;
            (kb, nb)
        };
        let side = if negative { "a" } else { "b" };
        shape.name = if of > 1 { format!("{}/{side}{k}", body.name) } else { format!("{}/{side}", body.name) };
        faces += name_faces(&mut shape, &parent, cut, n, size, tol)?;
        out.push(shape);
    }
    Ok(Split::Cells(out, faces, coincident(body, &parent, cut, n, size)))
}

/// A body's faces as matching needs them: each face's analytic surface, name and triangles.
struct Faces {
    surfaces: Vec<Option<FaceSurface>>,
    names: Vec<Face>,
    mesh: SurfaceMesh,
}

impl Faces {
    fn of(body: &Shape, tol: f64) -> Result<Faces> {
        Ok(Faces { surfaces: body.face_surfaces()?, names: body.faces.clone(), mesh: body.tessellate(tol)? })
    }

    /// The face nearest to `p` (by distance to its triangles), among those `keep` allows.
    fn nearest(&self, p: [f64; 3], keep: impl Fn(usize) -> bool) -> Option<usize> {
        let mut best: Option<(f64, usize)> = None;
        for (t, &f) in self.mesh.triangles.iter().zip(&self.mesh.tri_face) {
            let f = f as usize;
            if !keep(f) {
                continue;
            }
            let q = |k: usize| t.get(k).and_then(|&i| self.mesh.positions.get(i as usize)).copied().unwrap_or([0.0; 3]);
            let d = point_triangle(p, q(0), q(1), q(2));
            if best.is_none_or(|(b, _)| d < b) {
                best = Some((d, f));
            }
        }
        best.map(|(_, f)| f)
    }
}

impl Shape {
    /// For each point, the face (by index, the [`SurfaceMesh::tri_face`] tag) nearest to it,
    /// measured to a tessellation within `tolerance`; `None` for a body without faces.
    pub fn nearest_faces(&self, points: &[[f64; 3]], tolerance: f64) -> Result<Vec<Option<u32>>> {
        let faces = Faces { surfaces: Vec::new(), names: Vec::new(), mesh: self.tessellate(tolerance)? };
        Ok(points.iter().map(|p| faces.nearest(*p, |_| true).and_then(|f| u32::try_from(f).ok())).collect())
    }
}

/// Whether two analytic surfaces are the same surface (within `tol`); free-form faces match any
/// free-form face.
pub(crate) fn same_surface(a: Option<FaceSurface>, b: Option<FaceSurface>, tol: f64) -> bool {
    let close = |x: [f64; 3], y: [f64; 3]| dot(sub(x, y), sub(x, y)).sqrt() <= tol;
    let parallel = |x: [f64; 3], y: [f64; 3]| dot(cross(x, y), cross(x, y)).sqrt() <= PARALLEL_SINE;
    let on_line = |p: [f64; 3], o: [f64; 3], axis: [f64; 3]| {
        let c = cross(sub(p, o), axis);
        dot(c, c).sqrt() <= tol
    };
    match (a, b) {
        (None, None) => true,
        (Some(FaceSurface::Plane { normal: n1, point: p1 }), Some(FaceSurface::Plane { normal: n2, point: p2 })) => {
            dot(n1, n2) > 1.0 - SAME_DIRECTION && dot(sub(p2, p1), n1).abs() <= tol
        }
        (Some(FaceSurface::Cylinder { origin: o1, axis: a1, radius: r1 }), Some(FaceSurface::Cylinder { origin: o2, axis: a2, radius: r2 })) => {
            (r1 - r2).abs() <= tol && parallel(a1, a2) && on_line(o2, o1, a1)
        }
        (Some(FaceSurface::Cone { apex: p1, axis: a1, half_angle: h1 }), Some(FaceSurface::Cone { apex: p2, axis: a2, half_angle: h2 })) => {
            close(p1, p2) && dot(a1, a2) > 1.0 - SAME_DIRECTION && (h1 - h2).abs() <= SAME_ANGLE
        }
        (Some(FaceSurface::Sphere { center: c1, radius: r1 }), Some(FaceSurface::Sphere { center: c2, radius: r2 })) => {
            close(c1, c2) && (r1 - r2).abs() <= tol
        }
        (
            Some(FaceSurface::Torus { center: c1, axis: a1, major: m1, minor: n1 }),
            Some(FaceSurface::Torus { center: c2, axis: a2, major: m2, minor: n2 }),
        ) => close(c1, c2) && parallel(a1, a2) && (m1 - m2).abs() <= tol && (n1 - n2).abs() <= tol,
        _ => false,
    }
}

/// The distance from `p` to the triangle `a b c` (Ericson, Real-Time Collision Detection §5.1.5).
fn point_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    let at = |q: [f64; 3]| dot(sub(p, q), sub(p, q)).sqrt();
    let lerp = |o: [f64; 3], d: [f64; 3], t: f64| [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t];
    if d1 <= 0.0 && d2 <= 0.0 {
        return at(a);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return at(b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return at(lerp(a, ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return at(c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return at(lerp(a, ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return at(lerp(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))));
    }
    let denom = va + vb + vc;
    if denom.abs() < f64::MIN_POSITIVE {
        return at(a);
    }
    let (v, w) = (vb / denom, vc / denom);
    at([a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w])
}

/// The name without a trailing piece number (`F2:top#1` → `F2:top`).
pub(crate) fn base_name(name: &str) -> &str {
    match name.rsplit_once('#') {
        Some((b, k)) if !k.is_empty() && k.bytes().all(|c| c.is_ascii_digit()) => b,
        _ => name,
    }
}

/// Name a cell's faces: cut faces after the cut, the rest after the face of `parent` they are a
/// piece of. Returns how many cut faces it has.
fn name_faces(cell: &mut Shape, parent: &Faces, cut: &Cut, n: [f64; 3], size: f64, tol: f64) -> Result<usize> {
    let own = Faces::of(cell, tol)?;
    // A point on each face (the centroid of its largest triangle) and its area-weighted centroid.
    let count = cell.faces.len();
    let mut probe = vec![(0.0f64, [0.0; 3]); count];
    let mut centroid = vec![([0.0; 3], 0.0f64); count];
    for (t, &f) in own.mesh.triangles.iter().zip(&own.mesh.tri_face) {
        let q = |k: usize| t.get(k).and_then(|&i| own.mesh.positions.get(i as usize)).copied().unwrap_or([0.0; 3]);
        let (a, b, c) = (q(0), q(1), q(2));
        let cr = cross(sub(b, a), sub(c, a));
        let area = 0.5 * dot(cr, cr).sqrt();
        let mid = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0];
        if let Some(slot) = probe.get_mut(f as usize).filter(|s| area > s.0) {
            *slot = (area, mid);
        }
        if let Some((s, w)) = centroid.get_mut(f as usize) {
            for k in 0..3 {
                s[k] += area * mid[k];
            }
            *w += area;
        }
    }
    let centroid: Vec<[f64; 3]> = centroid.into_iter().map(|(s, w)| if w > 0.0 { s.map(|x| x / w) } else { s }).collect();
    let on_plane = ON_PLANE * size.max(1.0);
    let mut names: Vec<(String, bool)> = Vec::with_capacity(count); // (name, exact)
    let mut cut_faces = 0;
    for i in 0..count {
        let surface = own.surfaces.get(i).copied().flatten();
        if let Some(FaceSurface::Plane { normal, point }) = surface {
            let along = dot(normal, n);
            if (along.abs() - 1.0).abs() < SAME_DIRECTION && dot(sub(point, cut.point), n).abs() < on_plane.max(tol * ON_FACE) {
                // Outward along the normal: the cell is on the negative side.
                let side = if along > 0.0 { "a" } else { "b" };
                names.push((format!("partition:{}:{side}", cut.name), true));
                cut_faces += 1;
                continue;
            }
        }
        let p = probe.get(i).map_or([0.0; 3], |x| x.1);
        let from =
            parent.nearest(p, |f| same_surface(surface, parent.surfaces.get(f).copied().flatten(), tol)).or_else(|| parent.nearest(p, |_| true));
        let face = from.and_then(|f| parent.names.get(f));
        names.push(match face {
            Some(f) if !f.name.is_empty() => (f.name.clone(), f.exact),
            _ => (i.to_string(), false),
        });
    }
    // Several pieces of one face on the cell: numbered by position.
    let mut order: Vec<usize> = (0..count).collect();
    order.sort_by(|&a, &b| lexical(centroid.get(a).unwrap_or(&[0.0; 3]), centroid.get(b).unwrap_or(&[0.0; 3])));
    let bases: Vec<String> = names.iter().map(|(nm, _)| base_name(nm).to_string()).collect();
    let mut numbered = vec![None; count];
    for (i, base) in bases.iter().enumerate() {
        if numbered.get(i).is_some_and(Option::is_some) {
            continue;
        }
        let group: Vec<usize> = order.iter().copied().filter(|&j| bases.get(j) == Some(base)).collect();
        // Untouched pieces keep the numbers SolidCraft gave them; pieces the cut made (or a
        // name two pieces now share) are numbered again.
        let distinct = {
            let mut v: Vec<&str> = group.iter().filter_map(|&j| names.get(j).map(|x| x.0.as_str())).collect();
            v.sort_unstable();
            v.dedup();
            v.len() == group.len()
        };
        for (k, &j) in group.iter().enumerate() {
            if let Some(slot) = numbered.get_mut(j) {
                *slot = Some(if distinct { names.get(j).map(|x| x.0.clone()).unwrap_or_default() } else { format!("{base}#{k}") });
            }
        }
    }
    for ((face, name), (_, exact)) in cell.faces.iter_mut().zip(numbered).zip(&names) {
        face.name = name.unwrap_or_default();
        face.exact = *exact;
    }
    Ok(cut_faces)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::open_design;
    use crate::sample::{self, Plate};

    fn plate() -> Design {
        open_design(&sample::plate_design(&Plate::default()).unwrap()).unwrap()
    }

    fn box_design() -> Design {
        let mut shape = crate::block([0.0; 3], [40.0, 30.0, 10.0]).unwrap();
        shape.faces = (0..6).map(|i| Face { name: format!("box:{i}"), exact: true, members: Vec::new() }).collect();
        let bodies = vec![shape];
        Design {
            name: "Box".into(),
            source: std::sync::Arc::new(bodies.clone()),
            opened: std::sync::Arc::new(bodies.clone()),
            bodies,
            errors: Vec::new(),
            cuts: Vec::new(),
            history: None,
            defeat: None,
        }
    }

    fn cut(name: &str, point: [f64; 3], normal: [f64; 3]) -> Cut {
        Cut { name: name.into(), point, normal, bodies: Vec::new() }
    }

    fn names(b: &Shape) -> Vec<&str> {
        b.faces().iter().map(|f| f.name.as_str()).collect()
    }

    /// Coordinate-axis cuts retain the whole volume, paired faces on the requested plane,
    /// and the negative/positive cell names for either direction of the normal.
    #[test]
    fn coordinate_axis_tools_preserve_both_sides_and_face_names() {
        let d = box_design();
        let whole = d.bodies[0].measure().unwrap().volume;
        let original: Vec<&str> = names(&d.bodies[0]).into_iter().map(base_name).collect();
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let normal = std::array::from_fn(|i| if i == axis { sign } else { 0.0 });
                let at = [13.0, 11.0, 4.0];
                let p = d.partitioned(&[cut("Axis", at, normal)]).unwrap();
                assert!(p.errors.is_empty(), "axis {axis}, sign {sign}: {:?}", p.errors);
                assert_eq!(p.bodies.len(), 2);
                let volume: f64 = p.bodies.iter().map(|b| b.measure().unwrap().volume).sum();
                assert!((volume - whole).abs() < 1e-6 * whole, "axis {axis}: {volume} vs {whole}");
                for b in &p.bodies {
                    let negative = dot(sub(b.measure().unwrap().centroid, at), normal) < 0.0;
                    let side = if negative { "a" } else { "b" };
                    assert!(b.name.ends_with(&format!("/{side}")), "{}", b.name);
                    let mut cut_faces = 0;
                    for (face, surface) in b.faces().iter().zip(b.face_surfaces().unwrap()) {
                        if face.name.starts_with("partition:") {
                            assert!(face.name.starts_with(&format!("partition:Axis:{side}")), "{}", face.name);
                            let Some(FaceSurface::Plane { normal: outward, point }) = surface else { panic!("nonplanar cut face") };
                            assert!((point[axis] - at[axis]).abs() < 1e-9, "{point:?}");
                            assert!((dot(outward, normal) - if negative { 1.0 } else { -1.0 }).abs() < 1e-12);
                            cut_faces += 1;
                        } else {
                            assert!(original.contains(&base_name(&face.name)), "{}", face.name);
                        }
                        assert!(face.exact, "{}", face.name);
                    }
                    assert!(cut_faces > 0);
                    let mut all = names(b);
                    all.sort_unstable();
                    let before = all.len();
                    all.dedup();
                    assert_eq!(before, all.len());
                }
            }
        }
    }

    #[test]
    fn oblique_planes_keep_the_generic_tool_even_extremely_near_an_axis() {
        let d = box_design();
        let whole = d.bodies[0].measure().unwrap().volume;
        for normal in [[1.0, 1.0, 0.25], [1.0, 1e-12, 0.0]] {
            let cut = cut("Oblique", [13.0, 11.0, 4.0], normal);
            let (plane, n) = plane_of(&cut).unwrap();
            // In the second case n[0] rounds to exactly 1. Testing only that component would
            // incorrectly move the plane; the other nonzero component must exclude the box.
            assert!(n[1] != 0.0);
            assert!(axis_box(&plane, n, &solidcraft_kernel::measure(&d.bodies[0].body).unwrap().bbox, true).unwrap().is_none());
            let p = d.partitioned(&[cut]).unwrap();
            assert!(p.errors.is_empty(), "{normal:?}: {:?}", p.errors);
            assert_eq!(p.bodies.len(), 2);
            let volume: f64 = p.bodies.iter().map(|b| b.measure().unwrap().volume).sum();
            assert!((volume - whole).abs() < 1e-6 * whole, "{volume} vs {whole}");
        }
    }

    /// The issue's original SolidCraft primitive recipe, including its revolved relief groove.
    /// Building it and cutting all three planes belongs to the heavy CAD regression suite.
    #[cfg(feature = "scripts")]
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn a_revolved_groove_partitions_on_parallel_and_perpendicular_planes() {
        let recipe = r#"{"commands":[
            {"command":"solid.box","params":{"length":60,"width":30,"height":60,"corner":[-30,0,-30],"body_name":"Block"}},
            {"command":"solid.cylinder","params":{"diameter":30,"height":7.8,"base":[0,-0.01,0],"axis":[0,1,0],"operation":"cut","targets":["Block"]}},
            {"command":"solid.cylinder","params":{"diameter":14,"height":32,"base":[0,-1,0],"axis":[0,1,0],"operation":"cut","targets":["Block"]}},
            {"command":"sketch.create","params":{"plane":"XY","offset":"0"}},
            {"command":"sketch.line","params":{"points":[[16,7.5],[14.5,7.5],[14.5,9]],"infer":false}},
            {"command":"sketch.arc.three_point","params":{"start":[14.5,9],"end":[16,7.5],"through":[17.060660,10.060660]}},
            {"command":"sketch.finish","params":{}},
            {"command":"solid.revolve","params":{"axis":"y","angle":"360 deg","direction":"symmetric","operation":"cut","targets":["Block"]}}
        ]}"#;
        let original = open_design(&crate::script_design(recipe).unwrap()).unwrap();
        assert!(original.errors.is_empty(), "{:?}", original.errors);
        for (scale, shift) in [(1.0, [0.0; 3]), (2.0, [120.0, -70.0, 45.0])] {
            let mut d = original.clone();
            for body in &mut d.bodies {
                body.body = solidcraft_kernel::transform_matrix(
                    &body.body,
                    [[scale, 0.0, 0.0, 0.0], [0.0, scale, 0.0, 0.0], [0.0, 0.0, scale, 0.0], [shift[0], shift[1], shift[2], 1.0]],
                )
                .unwrap();
            }
            d.source = std::sync::Arc::new(d.bodies.clone());
            let original_bounds = original.bodies[0].measure().unwrap();
            let placed_bounds = d.bodies[0].measure().unwrap();
            for (axis, offset) in shift.iter().enumerate() {
                assert!((placed_bounds.bbox_min[axis] - (original_bounds.bbox_min[axis] * scale + offset)).abs() < 1e-6);
                assert!((placed_bounds.bbox_max[axis] - (original_bounds.bbox_max[axis] * scale + offset)).abs() < 1e-6);
            }
            // measure() uses adaptive tessellations, not exact integrals. The unchanged
            // perpendicular control already differs from the whole by 2.04 ppm. A 10 ppm
            // measurement bound tests coverage without demanding excessive tessellations;
            // the exact box control above retains its stricter 1 ppm bound.
            let volume = |body: &Shape| body.measure().unwrap().volume;
            let whole = volume(&d.bodies[0]);
            for (point, normal) in [([0.3, 0.0, 0.0], [1.0, 0.0, 0.0]), ([0.0, 0.0, -10.0], [0.0, 0.0, 1.0]), ([0.0, 20.0, 0.0], [0.0, 1.0, 0.0])] {
                let point = std::array::from_fn(|i| point[i] * scale + shift[i]);
                let p = d.partitioned(&[cut("Groove", point, normal)]).unwrap();
                assert!(p.errors.is_empty(), "scale {scale}, {point:?}: {:?}", p.errors);
                assert_eq!(p.bodies.len(), 2);
                let total: f64 = p.bodies.iter().map(volume).sum();
                assert!((total - whole).abs() < 1e-5 * whole, "scale {scale}, {point:?}: {total} vs {whole}");
                for b in &p.bodies {
                    assert!(b.body.validity().is_empty(), "{}: {:?}", b.name, b.body.validity());
                    let mut found = false;
                    for (face, surface) in b.faces().iter().zip(b.face_surfaces().unwrap()) {
                        if face.name.starts_with("partition:Groove:") {
                            let Some(FaceSurface::Plane { point: on, .. }) = surface else { panic!("nonplanar cut face") };
                            assert!(dot(sub(on, point), normal).abs() < 1e-9, "{point:?}: {on:?}");
                            assert!(face.exact);
                            found = true;
                        }
                    }
                    assert!(found, "{} has no cut face", b.name);
                }
            }
        }
    }

    #[test]
    fn a_plane_cuts_the_plate_into_cells_with_paired_cut_faces() {
        let d = plate();
        let whole = d.bodies[0].measure().unwrap().volume;
        // Through the hole's centre (x = 15): each half has two cut faces, either side of the hole.
        let p = d.partitioned(&[cut("Mid", [15.0, 0.0, 0.0], [1.0, 0.0, 0.0])]).unwrap();
        assert!(p.errors.is_empty(), "{:?}", p.errors);
        assert_eq!(p.bodies.len(), 2, "{:?}", p.bodies.iter().map(|b| &b.name).collect::<Vec<_>>());
        let volume: f64 = p.bodies.iter().map(|b| b.measure().unwrap().volume).sum();
        assert!((volume - whole).abs() < 1e-6 * whole, "{volume} vs {whole}");
        assert_eq!(p.cuts, vec![CutReport { name: "Mid".into(), split: vec![d.bodies[0].name.clone()], faces: 4, warnings: Vec::new() }]);
        let mut sides = Vec::new();
        for b in &p.bodies {
            let m = b.measure().unwrap();
            let side = if m.centroid[0] < 15.0 { "a" } else { "b" };
            let mut cut: Vec<&str> = names(b).into_iter().filter(|n| n.starts_with("partition:")).collect();
            cut.sort_unstable();
            assert_eq!(cut, vec![format!("partition:Mid:{side}#0"), format!("partition:Mid:{side}#1")], "{:?}", names(b));
            sides.push(side);
        }
        sides.sort_unstable();
        assert_eq!(sides, ["a", "b"]);
        // Faces the cut did not touch keep their names; the cut's faces pair by index (#0 with
        // #0) and lie on the same piece of the plane.
        let original: Vec<&str> = names(&d.bodies[0]).into_iter().map(base_name).collect();
        let body = &d.bodies[0].name;
        let mut cells: Vec<&str> = p.bodies.iter().map(|b| b.name.as_str()).collect();
        cells.sort_unstable();
        assert_eq!(cells, [format!("{body}/a"), format!("{body}/b")]);
        for b in &p.bodies {
            let n = names(b);
            assert!(n.iter().filter(|n| !n.starts_with("partition:")).all(|n| original.contains(&base_name(n))), "{n:?} {original:?}");
            assert!(b.faces().iter().all(|f| f.exact), "{:?}", b.faces());
            let mut unique = n.clone();
            unique.sort_unstable();
            unique.dedup();
            assert_eq!(unique.len(), n.len(), "{n:?}");
        }
        // The cut is replayed from the design as opened: partitioning again replaces it.
        let again = p.partitioned(&[]).unwrap();
        assert_eq!(again.bodies.len(), 1);
    }

    #[test]
    fn cuts_follow_each_other_and_name_only_the_bodies_they_divide() {
        let d = plate();
        let p = d.partitioned(&[cut("X", [30.0, 0.0, 0.0], [1.0, 0.0, 0.0]), cut("Z", [0.0, 0.0, 5.0], [0.0, 0.0, -1.0])]).unwrap();
        assert!(p.errors.is_empty(), "{:?}", p.errors);
        assert_eq!(p.bodies.len(), 4);
        assert_eq!(p.cuts.len(), 2);
        assert_eq!(p.cuts[1].split.len(), 2);
        for b in &p.bodies {
            let n = names(b);
            assert_eq!(n.iter().filter(|n| n.starts_with("partition:X:")).count(), 1, "{n:?}");
            assert_eq!(n.iter().filter(|n| n.starts_with("partition:Z:")).count(), 1, "{n:?}");
            // The normal points down: the cell above z = 5 is on its negative side.
            let above = b.measure().unwrap().centroid[2] > 5.0;
            assert!(n.contains(&if above { "partition:Z:a" } else { "partition:Z:b" }), "{n:?}");
        }
    }

    /// A body cut before stands for its cells: listing it cuts the cells the plane divides, and
    /// a plane that misses all of them is an error, as for a whole body.
    #[test]
    fn a_cut_body_listed_again_cuts_its_cells() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let on = |name: &str, point, normal| Cut { bodies: vec![body.clone()], ..cut(name, point, normal) };
        let p = d.partitioned(&[on("X", [30.0, 0.0, 0.0], [1.0, 0.0, 0.0]), on("Z", [0.0, 0.0, 5.0], [0.0, 0.0, 1.0])]).unwrap();
        assert!(p.errors.is_empty(), "{:?}", p.errors);
        assert_eq!(p.bodies.len(), 4);
        // Only one cell (x < 30) reaches x = 20.
        let p = d.partitioned(&[on("X", [30.0, 0.0, 0.0], [1.0, 0.0, 0.0]), on("Y", [20.0, 0.0, 0.0], [1.0, 0.0, 0.0])]).unwrap();
        assert!(p.errors.is_empty(), "{:?}", p.errors);
        assert_eq!(p.bodies.len(), 3);
        let p = d.partitioned(&[on("X", [30.0, 0.0, 0.0], [1.0, 0.0, 0.0]), on("Out", [100.0, 0.0, 0.0], [1.0, 0.0, 0.0])]).unwrap();
        assert_eq!(p.errors.iter().map(|e| e.feature.as_str()).collect::<Vec<_>>(), ["partition Out"]);
        assert_eq!(p.bodies.len(), 2);
    }

    #[test]
    fn a_plane_that_misses_or_a_missing_body_is_reported_not_fatal() {
        let d = plate();
        let body = d.bodies[0].name.clone();
        let miss = Cut { bodies: vec![body.clone()], ..cut("Out", [100.0, 0.0, 0.0], [1.0, 0.0, 0.0]) };
        let gone = Cut { bodies: vec!["Nope".into()], ..cut("Gone", [10.0, 0.0, 0.0], [1.0, 0.0, 0.0]) };
        let p = d.partitioned(&[miss, gone, cut("Free", [100.0, 0.0, 0.0], [1.0, 0.0, 0.0])]).unwrap();
        assert_eq!(p.bodies.len(), 1);
        let features: Vec<&str> = p.errors.iter().map(|e| e.feature.as_str()).collect();
        assert_eq!(features, ["partition Out", "partition Gone"], "{:?}", p.errors);
        assert_eq!(p.bodies[0].name, body);
        assert!(p.cuts.iter().all(|c| c.split.is_empty() && c.faces == 0), "{:?}", p.cuts);
    }

    #[test]
    fn the_nearest_face_of_a_point() {
        let d = plate();
        let b = &d.bodies[0];
        // Above the top (z = 10) and beside the side at y = 0.
        let f = b.nearest_faces(&[[30.0, 20.0, 10.5], [30.0, -0.5, 5.0]], 0.05).unwrap();
        let surfaces = b.face_surfaces().unwrap();
        let normal = |f: Option<u32>| match surfaces[f.unwrap() as usize] {
            Some(FaceSurface::Plane { normal, .. }) => normal,
            s => panic!("{s:?}"),
        };
        assert!((normal(f[0])[2] - 1.0).abs() < 1e-12);
        assert!((normal(f[1])[1] + 1.0).abs() < 1e-12);
    }

    #[test]
    fn hostile_planes_are_an_error() {
        let d = plate();
        for (point, normal) in [
            ([0.0; 3], [0.0; 3]),
            ([f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0; 3], [f64::INFINITY, 0.0, 0.0]),
            ([1e300, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ] {
            assert!(d.partitioned(&[cut("Bad", point, normal)]).is_err(), "{point:?} {normal:?}");
        }
        let many = vec![cut("M", [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]); MAX_CUTS + 1];
        assert!(d.partitioned(&many).is_err());
    }
}
