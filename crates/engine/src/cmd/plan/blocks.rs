//! The plan's block topology for the projection blocking mesher (#402 step 3, #399): the
//! translation of a mesh plan into one block topology over the whole body, so the rings round
//! its circles and the blend layer along its fillets, chamfers and shell walls come out of one
//! conforming mesh (no partitions, no ties, one division count per sheet of blocks).
//!
//! The body is read in a frame whose third axis is the plan's sweep direction (a principal axis):
//! - **Grid.** The planes of the body's planar faces square to each axis (its sides, a shell's
//!   inner walls and floor, a boss's top, a blind hole's bottom) and the frames round its
//!   circles make a grid of blocks over the body.
//! - **Rings.** The frame cell round each circle (concentric circles together) is O-gridded
//!   through every level, square inside square: the frame, then a square with its corners at
//!   45° on each circle, then a core. Corners on a circle sit on it before projection, so the
//!   ring's radial block edges meet the cylinder square on. A hole's frame is centred on it: a
//!   side of the body near the hole is the frame there and the far side is as far; a fillet
//!   round a corner of the outline that holds a hole (a screw hole in a rounded corner) is taken
//!   in by the hole's ring, its frame reaching the two sides and its outer corner projected onto
//!   the fillet's middle, so the corner needs no blend layer. A boss the body's walls run into
//!   (planes of the body crossing its circle outside it) is framed by the square round its
//!   whole circle, subdivided by those planes, every square inside it divided alike so the
//!   ring's blocks meet the walls node to node; on the walls' side the frame stops at their
//!   inner plane and its corner projects onto the inner fillet's middle (a boss just touching
//!   the fillet included).
//! - **Blend layer.** A side of the body a fillet or chamfer runs along, or that a shell's wall
//!   stands on, gets a layer of blocks between the grid (on the body shrunk by the blend's
//!   depth or the wall's thickness) and the side, one block per boundary face of the grid, the
//!   layer's vertices pushed out square to every side they lie on: diagonally at an edge and a
//!   corner. Projected, a corner vertex lands on the fillet (or chamfer) at its middle, so the
//!   two blocks meeting there share the blend's 180° as two right angles; on a shell's inner
//!   corner, as three 120°-ish ones. The layer's depth is the blend's measured depth (its
//!   stored radius or distance) or the wall's thickness; its divisions follow the element size.
//! - **Carving.** Every block whose middle lies outside the body (a cavity, the space round a
//!   boss, a hole's core) is removed, by the generalized winding number of the body's
//!   tessellation (solid angles, Van Oosterom and Strackee 1983).
//!
//! The projection mesher then projects the outer block faces onto the body's faces (each the
//! smooth group nearest it) and smooths the grid. Parallel block edges share one division count
//! across the whole topology, so neighbouring blocks always match: no shear bands, no ties.
//!
//! A body this layout does not fit (a face at an angle to the axes that is not a chamfer, a
//! curved face that is neither a circle's cylinder nor a blend, two circles in one frame) is
//! reported with the reason; the plan then keeps its frame cuts and swept cells.

use fieldcraft_geometry::{FaceSurface, Shape, SurfaceMesh};
use fieldcraft_mesh::project::Topology;
use fieldcraft_model::{PlanRing, RingKind};

type V3 = [f64; 3];

/// Most blocks a plan's layout makes before carving.
const MAX_LAYOUT_BLOCKS: usize = 20_000;
/// Most triangles the winding number classification sums over per point.
const MAX_CLASSIFY_TRIANGLES: usize = 400_000;

/// A body's block layout.
#[derive(Debug, Clone)]
pub(crate) struct Layout {
    pub(crate) topology: Topology,
    /// Circles ringed (concentric ones each).
    pub(crate) rings: usize,
    /// The sides with a blend layer: (side, depth, what).
    pub(crate) layers: Vec<(String, f64, &'static str)>,
    /// Grid planes per axis of the frame (across, across, along the sweep).
    pub(crate) planes: [usize; 3],
    /// The frame round each circle group: its axis point (the sweep's coordinate 0) and the
    /// smallest half-side of the frame.
    pub(crate) frames: Vec<(V3, f64)>,
    /// The ring-and-pillow layout ([`super::pillow`]): the pillow's least thickness along the walls.
    pub(crate) pillow: Option<f64>,
}

impl Layout {
    /// The half-side of the frame round the circle about the axis through `center` along the
    /// sweep `direction`.
    pub(crate) fn frame_of(&self, center: V3, direction: V3) -> Option<f64> {
        self.frames.iter().find_map(|(p, half)| {
            let off = sub(center, *p);
            let along = dot(off, direction);
            (dot(off, off) - along * along < 1e-12 * dot(center, center).max(1.0)).then_some(*half)
        })
    }

    pub(crate) fn summary(&self, show: super::Show) -> String {
        let length = |v: f64| show(fieldcraft_model::units::Quantity::Length, v);
        let layers: Vec<String> = self.layers.iter().map(|(s, d, w)| format!("{s} {} ({w})", length(*d))).collect();
        format!(
            "projection blocking: {} blocks over a {}×{}×{} grid, {} circle(s) ringed{}{}",
            self.topology.blocks.len(),
            self.planes[0].saturating_sub(1),
            self.planes[1].saturating_sub(1),
            self.planes[2].saturating_sub(1),
            self.rings,
            if layers.is_empty() { String::new() } else { format!(", a blend layer on {}", layers.join(", ")) },
            self.pillow.map_or_else(String::new, |t| format!(", a pillow layer from the carved grid to the walls ({} and more)", length(t)))
        )
    }
}

/// What the layout reads.
pub(crate) struct Input<'a> {
    pub(crate) shape: &'a Shape,
    pub(crate) surfaces: &'a [Option<FaceSurface>],
    pub(crate) mesh: &'a SurfaceMesh,
    pub(crate) direction: V3,
    pub(crate) rings: &'a [PlanRing],
    /// The faces fillets and chamfers made (indices on `shape`).
    pub(crate) blends: Vec<u32>,
    /// A shell's wall thickness.
    pub(crate) shell: Option<f64>,
    pub(crate) size: f64,
    /// How a message shows a length or an area.
    pub(crate) show: super::Show<'a>,
}

/// The local frame: local axis `a` is world axis `PERM[k][a]`, `k` the sweep's world axis (a
/// cyclic permutation, so blocks keep their handedness).
pub(super) const PERM: [[usize; 3]; 3] = [[1, 2, 0], [2, 0, 1], [0, 1, 2]];

const SIDE_NAMES: [[&str; 2]; 3] = [["-x", "+x"], ["-y", "+y"], ["-z", "+z"]];

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Twice the solid angle triangle `a b c` subtends at the origin, signed (Van Oosterom and
/// Strackee 1983: tan(Ω/2) = a·(b×c) / (|a||b||c| + (a·b)|c| + (a·c)|b| + (b·c)|a|)).
fn solid_angle(a: V3, b: V3, c: V3) -> f64 {
    let (la, lb, lc) = (norm(a), norm(b), norm(c));
    let num = dot(a, cross(b, c));
    let den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la;
    2.0 * num.atan2(den)
}

/// The generalized winding number of the closed surface `mesh` at `p`: 1 inside, 0 outside.
pub(super) fn winding(mesh: &SurfaceMesh, p: V3) -> f64 {
    let at = |i: u32| mesh.positions.get(i as usize).map(|q| sub(*q, p));
    let mut sum = 0.0;
    for t in mesh.triangles.iter().take(MAX_CLASSIFY_TRIANGLES) {
        if let (Some(a), Some(b), Some(c)) = (at(t[0]), at(t[1]), at(t[2])) {
            sum += solid_angle(a, b, c);
        }
    }
    sum / (4.0 * std::f64::consts::PI)
}

/// A circle's axis point and radius in the local frame, with its kind.
#[derive(Debug, Clone, Copy)]
struct Circle {
    c: [f64; 2],
    r: f64,
    kind: RingKind,
}

/// Concentric circles framed together: radii descending.
#[derive(Debug, Clone)]
struct Group {
    c: [f64; 2],
    circles: Vec<Circle>,
    /// The frame's extent per local axis (0, 1): [low, high].
    frame: [[f64; 2]; 2],
    /// The frame's half-side per local axis when a corner blend it takes in sets it: the
    /// distance from the centre to the body's side the blend runs along.
    forced: [Option<f64>; 2],
    /// A boss the body's walls run into: planes of the body cross its circle (outside it, where
    /// the walls join), so its frame is the square round the whole circle, subdivided by those
    /// planes, and the ring's blocks follow the subdivision.
    patch: bool,
    /// The extent per local axis of the concave blends a patch takes in (the inner fillet where
    /// the walls meet, in pieces where the boss touches it): the frame's side on the blends'
    /// side of the centre goes to their middle.
    blends: Option<[[f64; 2]; 2]>,
}

/// Plan the block topology of the body.
// Axes index several per-axis arrays at once (planes, frames, layers).
#[allow(clippy::needless_range_loop)]
pub(crate) fn layout(inp: &Input) -> Result<Layout, String> {
    let length = |v: f64| (inp.show)(fieldcraft_model::units::Quantity::Length, v);
    let d = inp.direction;
    let k = (0..3).find(|&k| (d[k].abs() - 1.0).abs() < 1e-9).ok_or("the sweep direction is not along an axis")?;
    let perm = PERM[k];
    let local = |p: V3| -> V3 { [p[perm[0]], p[perm[1]], p[perm[2]]] };
    let world = |l: V3| -> V3 {
        let mut w = [0.0; 3];
        for a in 0..3 {
            w[perm[a]] = l[a];
        }
        w
    };
    if inp.mesh.triangles.is_empty() || inp.mesh.triangles.len() > MAX_CLASSIFY_TRIANGLES {
        return Err(format!("the body's tessellation has {} triangles", inp.mesh.triangles.len()));
    }
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut face_points: Vec<Vec<V3>> = vec![Vec::new(); inp.surfaces.len()];
    for (t, &f) in inp.mesh.triangles.iter().zip(&inp.mesh.tri_face) {
        for &i in t {
            let Some(p) = inp.mesh.positions.get(i as usize).map(|p| local(*p)) else { continue };
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
            if let Some(list) = face_points.get_mut(f as usize) {
                list.push(p);
            }
        }
    }
    let diag = norm(sub(hi, lo));
    if !(diag.is_finite() && diag > 0.0) {
        return Err("the body has no extent".into());
    }
    let eps = 1e-6 * diag;
    let name = |f: usize| inp.shape.faces().get(f).map_or_else(|| format!("#{f}"), |x| x.name.clone());

    // Circles: the plan's rings along the sweep, in the local frame.
    let mut circles: Vec<Circle> = Vec::new();
    for r in inp.rings {
        if (dot(r.axis, d).abs() - 1.0).abs() > 1e-6 {
            return Err(format!("the ring of F{} is not square to the sweep", r.feature));
        }
        let c = local(r.center);
        circles.push(Circle { c: [c[0], c[1]], r: r.radius, kind: r.kind });
    }
    let on_circle = |origin: V3, radius: f64| {
        let o = local(origin);
        circles.iter().any(|c| (c.r - radius).abs() < 1e-6 * radius.max(1.0) && (c.c[0] - o[0]).hypot(c.c[1] - o[1]) < 1e-6 * radius.max(1.0))
    };

    // The faces: planes square to an axis (with their outward sign), the circles' cylinders, and
    // the blends; anything else is not this layout's.
    let mut planes_on: [Vec<(f64, f64)>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    // The planar faces square to an in-plane axis: (axis, offset, face).
    let mut side_faces: Vec<(usize, f64, usize)> = Vec::new();
    let mut blend_faces: Vec<usize> = Vec::new();
    for (f, s) in inp.surfaces.iter().enumerate() {
        let blend = inp.blends.contains(&(f as u32));
        match s {
            Some(FaceSurface::Plane { normal, point }) => {
                let n = local(*normal);
                match (0..3).find(|&a| (n[a].abs() - 1.0).abs() < 1e-6) {
                    Some(a) => {
                        planes_on[a].push((local(*point)[a], n[a].signum()));
                        if a < 2 {
                            side_faces.push((a, local(*point)[a], f));
                        }
                    }
                    None if blend => blend_faces.push(f),
                    None => return Err(format!("face `{}` is a plane at an angle to the axes", name(f))),
                }
            }
            // A ring's cylinder, even where a fillet's history names it (a boss over a filleted
            // corner).
            Some(FaceSurface::Cylinder { origin, axis, radius }) if (dot(*axis, d).abs() - 1.0).abs() < 1e-6 && on_circle(*origin, *radius) => {}
            Some(FaceSurface::Cylinder { .. }) if !blend => {
                return Err(format!("face `{}` is a cylinder that is not one of the plan's rings", name(f)));
            }
            _ if blend => blend_faces.push(f),
            _ => return Err(format!("face `{}` is curved and neither a ring's cylinder nor a fillet's", name(f))),
        }
    }

    // Circles, concentric ones grouped (radii descending).
    let mut groups: Vec<Group> = Vec::new();
    for c in &circles {
        match groups.iter_mut().find(|g| (g.c[0] - c.c[0]).hypot(g.c[1] - c.c[1]) < 1e-6 * c.r.max(1.0)) {
            Some(g) => {
                if !g.circles.iter().any(|x| (x.r - c.r).abs() < 1e-6 * c.r.max(1.0)) {
                    g.circles.push(*c);
                }
            }
            None => groups.push(Group { c: c.c, circles: vec![*c], frame: [[0.0; 2]; 2], forced: [None; 2], patch: false, blends: None }),
        }
    }
    for g in &mut groups {
        g.circles.sort_by(|a, b| b.r.total_cmp(&a.r));
        let Some(outer) = g.circles.first().copied() else { continue };
        let r0 = outer.r;
        let crossing: Vec<usize> = side_faces.iter().filter(|(a, v, _)| (v - g.c[*a]).abs() < r0 - eps).map(|(_, _, f)| *f).collect();
        if outer.kind != RingKind::Butterfly || crossing.is_empty() {
            continue;
        }
        // The walls join the boss outside its circle: no face of theirs inside it.
        for f in crossing {
            let pts = face_points.get(f).map_or(&[][..], Vec::as_slice);
            if pts.iter().any(|p| (p[0] - g.c[0]).hypot(p[1] - g.c[1]) < r0 * (1.0 - 1e-3)) {
                return Err(format!("face `{}` runs through the circle of radius {}", name(f), length(r0)));
            }
        }
        g.patch = true;
    }

    // The blend layer: per side, the depth of the blends along it, or the shell's wall. A blend
    // round a corner of the outline that a hole's ring can frame (a screw hole in a rounded
    // corner) is taken in by the ring instead.
    let mut ring = [[0.0f64; 2]; 3];
    let mut what = [[""; 2]; 3];
    for &f in &blend_faces {
        let pts = face_points.get(f).map_or(&[][..], Vec::as_slice);
        if pts.is_empty() {
            continue;
        }
        if let Some((g, half)) = corner_ring(&groups, pts, lo, hi, eps)
            && let Some(g) = groups.get_mut(g)
        {
            g.forced = half.map(Some);
            continue;
        }
        if let Some((g, extent)) = patch_blend(&groups, pts)
            && let Some(g) = groups.get_mut(g)
        {
            let union = g.blends.map_or(extent, |b| std::array::from_fn(|a| [b[a][0].min(extent[a][0]), b[a][1].max(extent[a][1])]));
            g.blends = Some(union);
            continue;
        }
        let mut touched = false;
        for a in 0..3 {
            for s in 0..2 {
                let depth = pts.iter().map(|p| if s == 0 { p[a] - lo[a] } else { hi[a] - p[a] }).fold(0.0f64, f64::max);
                if depth < 0.45 * (hi[a] - lo[a]) {
                    touched = true;
                    if depth > ring[a][s] {
                        ring[a][s] = depth;
                        what[a][s] = "blend";
                    }
                }
            }
        }
        if !touched {
            return Err(format!("the blend face `{}` does not run along a side of the body", name(f)));
        }
    }
    if let Some(t) = inp.shell.filter(|t| t.is_finite() && *t > 0.0) {
        // A wall: the inner face of the shell, `t` inside the side and facing away from it.
        for a in 0..3 {
            for s in 0..2 {
                let at = if s == 0 { lo[a] + t } else { hi[a] - t };
                let inward = if s == 0 { 1.0 } else { -1.0 };
                if planes_on[a].iter().any(|(v, sign)| (v - at).abs() < 1e-4 * t.max(1.0) && *sign == inward) {
                    ring[a][s] = t;
                    what[a][s] = "shell wall";
                }
            }
        }
    }
    let ilo: V3 = std::array::from_fn(|a| lo[a] + ring[a][0]);
    let ihi: V3 = std::array::from_fn(|a| hi[a] - ring[a][1]);
    if (0..3).any(|a| !(ihi[a] - ilo[a] > eps)) {
        return Err("the blend layers leave no room inside the body".into());
    }

    // The fixed planes per axis: faces square to it (a side with a layer moves in by its depth).
    let mut fixed: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
    for a in 0..3 {
        for &(v, _) in &planes_on[a] {
            let v = if ring[a][0] > 0.0 && (v - lo[a]).abs() < eps {
                ilo[a]
            } else if ring[a][1] > 0.0 && (v - hi[a]).abs() < eps {
                ihi[a]
            } else {
                v
            };
            if v < ilo[a] - eps || v > ihi[a] + eps {
                return Err(format!("a face square to {} lies inside the blend layer", SIDE_NAMES[perm[a]][0].trim_start_matches('-')));
            }
            fixed[a].push(v);
        }
        if ring[a][0] > 0.0 {
            fixed[a].push(ilo[a]);
        }
        if ring[a][1] > 0.0 {
            fixed[a].push(ihi[a]);
        }
        fixed[a].sort_by(f64::total_cmp);
        fixed[a].dedup_by(|x, y| (*x - *y).abs() < eps);
    }

    let snap = (0.3 * inp.size).min(0.25 * circles.iter().map(|c| c.r).fold(f64::INFINITY, f64::min));
    for g in &mut groups {
        let outer = g.circles.first().copied().ok_or("an empty circle group")?;
        let r0 = outer.r;
        for a in 0..2 {
            let c = g.c[a];
            if !g.patch && fixed[a].iter().any(|v| (v - c).abs() < r0 - eps) {
                return Err(format!("a face of the body runs through the circle of radius {}", length(r0)));
            }
        }
        let half = match outer.kind {
            // The square round the whole circle (its corners projected onto it at 45°), or
            // where a blend it takes in sets a side, at the walls' inner plane.
            RingKind::Butterfly if g.patch => std::array::from_fn(|a| {
                let mut h = [r0; 2];
                // The blends' near edge: the walls' inner plane, whose corner projects onto
                // the blend's middle.
                if let Some(b) = g.blends {
                    if 0.5 * (b[a][0] + b[a][1]) > g.c[a] {
                        h[1] = b[a][0] - g.c[a];
                    } else {
                        h[0] = g.c[a] - b[a][1];
                    }
                }
                h
            }),
            RingKind::Butterfly => [[r0 / 2f64.sqrt(); 2]; 2],
            RingKind::OGrid => o_grid_frame(g, &fixed, eps, inp.show)?,
        };
        for a in 0..2 {
            let c = g.c[a];
            for s in 0..2 {
                let h = half[a][s];
                if outer.kind == RingKind::Butterfly && !g.patch && fixed[a].iter().any(|v| (v - c).abs() < h - eps) {
                    return Err(format!("a face of the body runs through the square inside the circle of radius {}", length(r0)));
                }
                g.frame[a][s] = if s == 0 { c - h } else { c + h };
            }
        }
    }
    // Frames that overlap along an axis share their planes there (a plane may not cut another
    // frame); two frames that overlap along both are one cell, which the layout does not take.
    for _ in 0..=groups.len() * groups.len() {
        let mut changed = false;
        for a in 0..2 {
            for i in 0..groups.len() {
                for j in 0..groups.len() {
                    if i == j {
                        continue;
                    }
                    let (fi, fj) = (groups[i].frame[a], groups[j].frame[a]);
                    let overlap = fi[0] < fj[1] - eps && fj[0] < fi[1] - eps;
                    let same = (fi[0] - fj[0]).abs() < eps && (fi[1] - fj[1]).abs() < eps;
                    if overlap && !same {
                        let u = [fi[0].min(fj[0]), fi[1].max(fj[1])];
                        groups[i].frame[a] = u;
                        groups[j].frame[a] = u;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }
    for (i, g) in groups.iter().enumerate() {
        for h in groups.iter().skip(i + 1) {
            if (0..2).all(|a| g.frame[a][0] < h.frame[a][1] - eps && h.frame[a][0] < g.frame[a][1] - eps) {
                return Err(format!(
                    "the circles at ({:.4}, {:.4}) and ({:.4}, {:.4}) are too close for a frame each",
                    g.c[0], g.c[1], h.c[0], h.c[1]
                ));
            }
        }
    }
    // The planes: the fixed ones and the frames, a frame within `snap` of another plane moved
    // onto it (a fixed plane wins).
    let mut planes: [Vec<f64>; 3] = fixed.clone();
    for a in 0..2 {
        for g in &mut groups {
            for s in 0..2 {
                let v = g.frame[a][s];
                let near = fixed[a].iter().copied().filter(|p| (p - v).abs() < snap).min_by(|p, q| (p - v).abs().total_cmp(&(q - v).abs()));
                let v = match near {
                    Some(p) => p,
                    None => planes[a].iter().copied().find(|p| (p - v).abs() < snap).unwrap_or(v),
                };
                g.frame[a][s] = v;
                planes[a].push(v);
            }
        }
        planes[a].sort_by(f64::total_cmp);
        planes[a].dedup_by(|x, y| (*x - *y).abs() < eps);
    }
    for g in &groups {
        for a in 0..2 {
            // A patch's frame may stop short of its circle where the walls join (the ring's
            // faces there are projected out onto it), but holds the circles inside.
            let r0 = g.circles.first().map_or(0.0, |c| c.r);
            let hold = if g.patch { 1.2 * g.circles.get(1).map_or(0.0, |c| c.r) } else { 0.7 * r0 };
            if !(g.frame[a][0] < g.c[a] - hold && g.frame[a][1] > g.c[a] + hold) {
                return Err(format!("the frame round the circle of radius {} does not hold it", length(r0)));
            }
            if !g.patch && planes[a].iter().any(|p| *p > g.frame[a][0] + eps && *p < g.frame[a][1] - eps) {
                return Err(format!("a plane of the body runs through the frame round the circle of radius {}", length(r0)));
            }
        }
    }
    if planes.iter().any(|p| p.len() < 2) {
        return Err("the body has no planar sides or circles to lay a grid between".into());
    }
    let n = planes.each_ref().map(|p| p.len() - 1);
    if n.iter().product::<usize>() * (1 + 4 * circles.len()) > MAX_LAYOUT_BLOCKS {
        return Err(format!("a grid of {n:?} cells is too large for one layout"));
    }

    // The grid, in the local frame; each block's point that says whether it is material.
    let mut t = Topology::grid([&planes[0], &planes[1], &planes[2]]).map_err(|e| e.to_string())?;
    let mid = |p: &[f64], i: usize| 0.5 * (p[i] + p[i + 1]);
    let mut probe: Vec<V3> = t
        .blocks
        .iter()
        .enumerate()
        .map(|(b, _)| [mid(&planes[0], b % n[0]), mid(&planes[1], (b / n[0]) % n[1]), mid(&planes[2], b / (n[0] * n[1]))])
        .collect();
    let index = |p: &[f64], v: f64| p.iter().position(|x| (x - v).abs() < eps);
    // The grid cells the rings replace.
    let mut replaced: Vec<usize> = Vec::new();
    let mut rings = 0;
    for g in &groups {
        let (Some(i0), Some(i1), Some(j0), Some(j1)) =
            (index(&planes[0], g.frame[0][0]), index(&planes[0], g.frame[0][1]), index(&planes[1], g.frame[1][0]), index(&planes[1], g.frame[1][1]))
        else {
            return Err("a frame lost its planes".into());
        };
        let r0 = g.circles.first().map_or(0.0, |c| c.r);
        // The squares inside the frame: one per circle (corners at 45° on it), then the core.
        // Each: half-sides [-u, +u, -v, +v] and the radius it stands for per side (a patch's
        // frame stands for its circle, which the frame's corners are projected onto).
        let mut squares: Vec<([f64; 4], [f64; 4])> = Vec::new();
        let fr = [g.c[0] - g.frame[0][0], g.frame[0][1] - g.c[0], g.c[1] - g.frame[1][0], g.frame[1][1] - g.c[1]];
        squares.push((fr, if g.patch { [r0; 4] } else { fr }));
        for c in &g.circles {
            let h = c.r / 2f64.sqrt();
            if squares.last().is_some_and(|(s, _)| s.iter().all(|x| (x - h).abs() < eps)) || (g.patch && c.r == r0) {
                continue;
            }
            squares.push(([h; 4], [c.r; 4]));
        }
        rings += g.circles.len();
        let inner = g.circles.last().map_or(0.0, |c| c.r);
        squares.push(([inner / (2.0 * 2f64.sqrt()); 4], [0.0; 4]));
        if squares.windows(2).any(|w| (0..4).any(|i| !(w[1].0[i] < w[0].0[i] - eps))) {
            return Err(format!("the frame round the circle of radius {} is too tight for its ring", length(r0)));
        }
        // The frame's subdivision by the grid's planes, as fractions along each side: every
        // square inside is divided alike, so the ring's blocks meet node to node.
        let frac = |a: usize, lo: usize, hi: usize| -> Vec<f64> {
            let (p0, p1) = (planes[a][lo], planes[a][hi]);
            planes[a][lo..=hi].iter().map(|p| (p - p0) / (p1 - p0)).collect()
        };
        let (fu, fv) = (frac(0, i0, i1), frac(1, j0, j1));
        // Point (u, v) (fractions along the axes) of square `k`.
        let at = |k: usize, u: f64, v: f64| -> [f64; 2] {
            let h = squares.get(k).map_or([0.0; 4], |s| s.0);
            [g.c[0] - h[0] + u * (h[0] + h[1]), g.c[1] - h[2] + v * (h[2] + h[3])]
        };
        // The frame's boundary, counter-clockwise, as (u, v) fractions with the side each
        // segment that starts there lies on (0: -u, 1: +u, 2: -v, 3: +v).
        let mut around: Vec<([f64; 2], usize)> = Vec::new();
        around.extend(fu[..fu.len() - 1].iter().map(|&u| ([u, 0.0], 2)));
        around.extend(fv[..fv.len() - 1].iter().map(|&v| ([1.0, v], 1)));
        around.extend(fu[1..].iter().rev().map(|&u| ([u, 1.0], 3)));
        around.extend(fv[1..].iter().rev().map(|&v| ([0.0, v], 0)));
        for w in 0..n[2] {
            let (z0, z1) = (planes[2][w], planes[2][w + 1]);
            let zmid = 0.5 * (z0 + z1);
            for i in i0..i1 {
                for j in j0..j1 {
                    replaced.push((w * n[1] + j) * n[0] + i);
                }
            }
            let add = |q: [[f64; 2]; 4], pr: V3, t: &mut Topology, probe: &mut Vec<V3>| {
                let base = t.vertices.len() as u32;
                for z in [z0, z1] {
                    for p in q {
                        t.vertices.push([p[0], p[1], z]);
                    }
                }
                let mut b = [base, base + 1, base + 2, base + 3, base + 4, base + 5, base + 6, base + 7];
                let pos = b.map(|v| t.vertices.get(v as usize).copied().unwrap_or_default());
                if corner_det(&pos) < 0.0 {
                    b = [b[0], b[3], b[2], b[1], b[4], b[7], b[6], b[5]];
                }
                t.blocks.push(b);
                probe.push(pr);
            };
            for (k, pair) in squares.windows(2).enumerate() {
                let (outer_rho, inner_rho) = (pair[0].1, pair[1].1);
                for (s, &(p, side)) in around.iter().enumerate() {
                    let Some(&(q, _)) = around.get((s + 1) % around.len()) else { continue };
                    let quad = [at(k, p[0], p[1]), at(k, q[0], q[1]), at(k + 1, q[0], q[1]), at(k + 1, p[0], p[1])];
                    // The probe: between the two radii the squares stand for, toward the
                    // segment's middle.
                    let m = [0.5 * (quad[0][0] + quad[1][0]) - g.c[0], 0.5 * (quad[0][1] + quad[1][1]) - g.c[1]];
                    let len = m[0].hypot(m[1]).max(eps);
                    let rho = 0.5 * (outer_rho[side] + inner_rho[side]);
                    add(quad, [g.c[0] + m[0] / len * rho, g.c[1] + m[1] / len * rho, zmid], &mut t, &mut probe);
                }
            }
            // The core, divided as the frame.
            let last = squares.len() - 1;
            for u in fu.windows(2) {
                for v in fv.windows(2) {
                    let quad = [at(last, u[0], v[0]), at(last, u[1], v[0]), at(last, u[1], v[1]), at(last, u[0], v[1])];
                    let cx = 0.25 * quad.iter().map(|p| p[0]).sum::<f64>();
                    let cy = 0.25 * quad.iter().map(|p| p[1]).sum::<f64>();
                    add(quad, [cx, cy, zmid], &mut t, &mut probe);
                }
            }
        }
    }
    if probe.len() != t.blocks.len() {
        return Err("the layout lost track of its blocks".into());
    }
    replaced.sort_unstable();
    replaced.dedup();
    t.remove(&replaced).map_err(|e| e.to_string())?;
    let mut k = 0;
    probe.retain(|_| {
        k += 1;
        replaced.binary_search(&(k - 1)).is_err()
    });
    t.merge(eps);

    // The blend layer: a block on every boundary face of the grid on a side with a layer.
    let mut uses: std::collections::BTreeMap<[u32; 4], usize> = std::collections::BTreeMap::new();
    for block in &t.blocks {
        for f in hex_faces(block) {
            let mut key = f;
            key.sort_unstable();
            *uses.entry(key).or_insert(0) += 1;
        }
    }
    let mut moved: std::collections::BTreeMap<u32, u32> = std::collections::BTreeMap::new();
    let original = t.blocks.len();
    for b in 0..original {
        let Some(block) = t.blocks.get(b).copied() else { continue };
        for f in hex_faces(&block) {
            let mut key = f;
            key.sort_unstable();
            if uses.get(&key) != Some(&1) {
                continue;
            }
            let pts = f.map(|v| t.vertices.get(v as usize).copied().unwrap_or_default());
            let side = (0..3).flat_map(|a| [(a, 0), (a, 1)]).find(|&(a, s)| {
                let plane = if s == 0 { ilo[a] } else { ihi[a] };
                ring[a][s] > 0.0 && pts.iter().all(|p| (p[a] - plane).abs() < eps)
            });
            let Some((a, s)) = side else { continue };
            let mut outer = [0u32; 4];
            for (o, &v) in outer.iter_mut().zip(&f) {
                *o = match moved.get(&v) {
                    Some(&w) => w,
                    None => {
                        let p = t.vertices.get(v as usize).copied().unwrap_or_default();
                        let mut q = p;
                        for aa in 0..3 {
                            if ring[aa][0] > 0.0 && (p[aa] - ilo[aa]).abs() < eps {
                                q[aa] -= ring[aa][0];
                            }
                            if ring[aa][1] > 0.0 && (p[aa] - ihi[aa]).abs() < eps {
                                q[aa] += ring[aa][1];
                            }
                        }
                        let w = t.vertices.len() as u32;
                        t.vertices.push(q);
                        moved.insert(v, w);
                        w
                    }
                };
            }
            let mut nb = [f[0], f[1], f[2], f[3], outer[0], outer[1], outer[2], outer[3]];
            let at = |v: u32| t.vertices.get(v as usize).copied().unwrap_or_default();
            if corner_det(&nb.map(at)) < 0.0 {
                nb = [f[0], f[3], f[2], f[1], outer[0], outer[3], outer[2], outer[1]];
            }
            let centre = pts.iter().fold([0.0; 3], |acc, p| [acc[0] + p[0] / 4.0, acc[1] + p[1] / 4.0, acc[2] + p[2] / 4.0]);
            let mut pr = centre;
            pr[a] += if s == 0 { -0.5 * ring[a][s] } else { 0.5 * ring[a][s] };
            t.blocks.push(nb);
            probe.push(pr);
        }
    }
    if t.blocks.len() > MAX_LAYOUT_BLOCKS {
        return Err(format!("{} blocks: too many for one layout", t.blocks.len()));
    }

    // Carve: blocks whose probe lies outside the body.
    let gone: Vec<usize> = probe.iter().enumerate().filter(|(_, p)| winding(inp.mesh, world(**p)) < 0.5).map(|(b, _)| b).collect();
    if gone.len() == t.blocks.len() {
        return Err("every block of the layout lies outside the body".into());
    }
    if !gone.is_empty() {
        t.remove(&gone).map_err(|e| e.to_string())?;
    }
    for v in &mut t.vertices {
        *v = world(*v);
    }
    t.check().map_err(|e| format!("the layout's blocks do not fit: {e}"))?;
    let layers = (0..3)
        .flat_map(|a| [(a, 0), (a, 1)])
        .filter(|&(a, s)| ring[a][s] > 0.0)
        .map(|(a, s)| (SIDE_NAMES[perm[a]][s].to_string(), ring[a][s], what[a][s]))
        .collect();
    let frames = groups
        .iter()
        .map(|g| {
            let half = (0..2).flat_map(|a| [g.c[a] - g.frame[a][0], g.frame[a][1] - g.c[a]]).fold(f64::INFINITY, f64::min);
            (world([g.c[0], g.c[1], 0.0]), half)
        })
        .collect();
    Ok(Layout { topology: t, rings, layers, planes: planes.each_ref().map(Vec::len), frames, pillow: None })
}

/// The ring that takes in the blend with points `pts` (local frame): a blend round a corner of
/// the outline, touching the body's side along both in-plane axes, that lies inside the frame of
/// an O-grid ring reaching from the hole to those sides, the hole centred in it. Returns the
/// group and the frame's half-side per axis (the hole's distance to the side). A frame wider
/// than 6 radii (or narrower than 1.2) does not take a blend in: the hole is not the corner's.
fn corner_ring(groups: &[Group], pts: &[V3], lo: V3, hi: V3, eps: f64) -> Option<(usize, [f64; 2])> {
    let (mut bl, mut bh) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in pts {
        for a in 0..2 {
            bl[a] = bl[a].min(p[a]);
            bh[a] = bh[a].max(p[a]);
        }
    }
    let tol = 1e3 * eps;
    groups.iter().enumerate().find_map(|(i, g)| {
        let outer = g.circles.first()?;
        if outer.kind != RingKind::OGrid {
            return None;
        }
        let mut half = [0.0; 2];
        for a in 0..2 {
            let h = if (bl[a] - lo[a]).abs() < tol {
                g.c[a] - lo[a]
            } else if (bh[a] - hi[a]).abs() < tol {
                hi[a] - g.c[a]
            } else {
                return None;
            };
            if !(1.2 * outer.r..=6.0 * outer.r).contains(&h) || bl[a] < g.c[a] - h - tol || bh[a] > g.c[a] + h + tol {
                return None;
            }
            half[a] = h;
        }
        Some((i, half))
    })
}

/// The patch ring that takes in the concave blend with points `pts` (local frame): a blend
/// inside the square round a boss's whole circle and not inside the circle (the inner fillet
/// where the walls meet the boss, the boss touching it at most). Returns the group and, per
/// axis, the blend's extent.
fn patch_blend(groups: &[Group], pts: &[V3]) -> Option<(usize, [[f64; 2]; 2])> {
    let (mut bl, mut bh) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in pts {
        for a in 0..2 {
            bl[a] = bl[a].min(p[a]);
            bh[a] = bh[a].max(p[a]);
        }
    }
    groups.iter().enumerate().find_map(|(i, g)| {
        let r0 = g.circles.first()?.r;
        let tol = 1e-3 * r0;
        if !g.patch || pts.iter().any(|p| (p[0] - g.c[0]).hypot(p[1] - g.c[1]) < 0.9 * r0) {
            return None;
        }
        if (0..2).any(|a| bl[a] < g.c[a] - r0 - tol || bh[a] > g.c[a] + r0 + tol) {
            return None;
        }
        Some((i, [[bl[0], bh[0]], [bl[1], bh[1]]]))
    })
}

/// The half-sides of an O-grid ring's frame, per local axis and side: the hole centred in it.
/// A side of the body within 3 radii (or a corner blend the ring takes in) is the frame on that
/// side, and the frame's other side is as far from the centre unless the body's side is nearer
/// there too; an axis with no side near takes the half-side of the other axis's (a square
/// frame) when that is 1.5 to 3 radii, else 2 radii.
fn o_grid_frame(g: &Group, fixed: &[Vec<f64>; 3], eps: f64, show: super::Show) -> Result<[[f64; 2]; 2], String> {
    let length = |v: f64| show(fieldcraft_model::units::Quantity::Length, v);
    let r0 = g.circles.first().map_or(0.0, |c| c.r);
    let mut half: [Option<[f64; 2]>; 2] = [None; 2];
    for a in 0..2 {
        let c = g.c[a];
        let mut dist = [0.0; 2];
        for (s, d) in dist.iter_mut().enumerate() {
            *d = fixed[a].iter().map(|v| if s == 0 { c - v } else { v - c }).filter(|x| *x > 0.0).fold(f64::INFINITY, f64::min);
            if !d.is_finite() {
                return Err(format!("the hole of radius {} has no side of the body round it", length(r0)));
            }
            if *d < 1.05 * r0 {
                return Err(format!("the hole of radius {} is {} from a face: too close for a ring", length(r0), length(*d)));
            }
        }
        let near = dist[0].min(dist[1]);
        half[a] = match g.forced[a] {
            Some(h) => Some([h; 2]),
            // The side itself frames the ring when a strip narrower than the radius would be
            // left between them; the other side keeps the frame centred unless it is near too.
            None if near <= 3.0 * r0 + eps => Some(dist.map(|d| if d <= 3.0 * r0 + eps { d } else { near })),
            None => None,
        };
    }
    let square = |b: usize| half[b].map(|h| h[0].min(h[1])).filter(|h| (1.5 * r0 - eps..=3.0 * r0 + eps).contains(h)).unwrap_or(2.0 * r0);
    Ok([half[0].unwrap_or([square(1); 2]), half[1].unwrap_or([square(0); 2])])
}

/// The six faces of a hex8 block, as corner cycles.
fn hex_faces(b: &[u32; 8]) -> [[u32; 4]; 6] {
    [
        [b[0], b[1], b[2], b[3]],
        [b[4], b[5], b[6], b[7]],
        [b[0], b[1], b[5], b[4]],
        [b[1], b[2], b[6], b[5]],
        [b[2], b[3], b[7], b[6]],
        [b[3], b[0], b[4], b[7]],
    ]
}

/// The Jacobian determinant at corner 1 of a hex8 block (positive: right-handed).
pub(super) fn corner_det(c: &[V3; 8]) -> f64 {
    dot(sub(c[1], c[0]), cross(sub(c[3], c[0]), sub(c[4], c[0])))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn winding_of_a_cube_is_one_inside_and_zero_outside() {
        let shape = fieldcraft_geometry::block([0.0; 3], [2.0, 2.0, 2.0]).unwrap();
        let mesh = shape.tessellate(0.01).unwrap();
        assert!((winding(&mesh, [1.0, 1.0, 1.0]) - 1.0).abs() < 1e-6);
        assert!(winding(&mesh, [3.0, 1.0, 1.0]).abs() < 1e-6);
    }
}
