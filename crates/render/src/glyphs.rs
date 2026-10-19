//! Symbols for a model's analysis setup, drawn over the model by the desktop viewport and the
//! headless renderer alike, so both show the same thing. Every symbol is a wireframe: thin lines
//! and dots, no shaded faces, so it never hides the mesh or a contour plot under it.
//!
//! - **Boundary conditions:** a wire triangle per constrained translation, tip on the node,
//!   pointing along that axis from outside the part (where the set's faces tell which side that
//!   is); a constrained rotation adds a second, smaller triangle behind the first (the double
//!   mark), only where the node can rotate (shell and beam nodes, reference points). Fixed shows
//!   all three axes; a symmetry condition only the axis normal to its plane. Thermal conditions
//!   are dots.
//! - **Forces:** straight arrows with an open head, ending on their nodes. **Moments:** an axis
//!   line through the node and a curved arrow three quarters of the way around it, turning
//!   right-handed about the moment vector.
//! - **Pressures:** arrows normal to the loaded faces (onto them for a positive pressure, away for
//!   a suction); **tractions:** arrows onto the faces along the traction. 2D models get them on
//!   the loaded edges, in the plane.
//! - **Gravity, body forces and centrifugal loads:** one arrow beside the model (a centrifugal
//!   load: its axis and a curved arrow around it). **Initial velocities:** one arrow behind the
//!   moving set, along its motion.
//! - **Reference points:** a three-axis cross with the point's name; **couplings:** thin lines
//!   from the point to some of the coupled nodes; **contact and tie surfaces:** small squares on
//!   their faces; **rigid walls:** a hatched square.
//!
//! **Size on screen.** Symbols keep a constant size on screen whatever the zoom: [`Glyphs::draw`]
//! takes the symbol unit (a support triangle's length) in model length, which the caller works
//! out from the camera with [`symbol_unit`] (a size in pixels, [`DEFAULT_SIZE`] unless the user
//! chose another, times the model length a pixel spans). Line widths and dots are fixed in pixels
//! too. Loads are drawn in the same unit as supports, so the two stay in proportion.
//!
//! **Density.** Where a definition covers many nodes or faces, its symbols are thinned on a grid
//! a few symbol units across, so they are spread evenly at a density that follows the screen,
//! not the mesh; each definition draws at most [`MAX_PER_DEFINITION`] symbols ([`MAX_SYMBOLS`] in
//! all), and one with more places than that is listed in [`Drawn::notes`], so the viewer can say
//! so.
//!
//! **Depth.** Symbols are lines, depth-tested against the model with the small bias toward the
//! eye that keeps lines on the face they lie on: they show on their own surface and hide behind
//! other parts of the model. A concentrated load (at a reference point or a single node) and a
//! reference point's cross are [`Scene::on_top`], visible from any side.
//!
//! **Colour** is per kind ([`GlyphColors`], from the theme's tokens); the symbols of the
//! definition selected in the browser ([`GlyphView::highlight`]) are drawn in the highlight
//! colour, wider.
//!
//! [`Glyphs::from_model`] does the work that depends on the model only (resolving sets, the
//! places symbols may go); [`Glyphs::draw`] chooses among them for the symbol unit and places
//! them at the nodes' current positions (a results plot's deformed shape) in the given colours.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use fieldcraft_model::{BoundaryCondition, ContactSurfaces, ElementFace, Load, Model, Procedure, ResolvedSet, Target};

use crate::display::{centre, polygon_normal};
use crate::{Aabb3, Rgb, Scene, SceneFill, SceneLabel, SceneLine, ScenePoint, SceneTri, Vec3};

/// The symbol unit drawn when the caller gives none: a fraction of the model's bounding box
/// diagonal.
const UNIT_FRACTION: f64 = 0.025;
/// A symbol unit's default size on screen, in pixels (points on the desktop).
pub const DEFAULT_SIZE: f64 = 14.0;
/// The symbol sizes one may choose, in pixels.
pub const SIZE_RANGE: (f64, f64) = (6.0, 48.0);
/// Most symbols one definition draws.
pub const MAX_PER_DEFINITION: usize = 400;
/// Most symbols drawn in all.
pub const MAX_SYMBOLS: usize = 6000;
/// Most lines one coupling draws to its nodes.
const MAX_COUPLING_LINES: usize = 24;
/// How far symbols over a results plot fade toward the background ([`GlyphColors::over_results`]).
pub const MUTED: f64 = 0.55;
/// Most nodes or faces of one set looked at when choosing where symbols go (a larger set is
/// strided down to this first).
const MAX_SCAN: usize = 200_000;
/// Segments of a curved arrow.
const ARC_SEGMENTS: usize = 18;
/// Steps a symbol unit is rounded to per doubling ([`quantise_unit`]).
const UNIT_STEPS: f64 = 4.0;

/// Symbol line width, and the selected definition's, in output pixels.
const LINE_WIDTH: f64 = 1.5;
const HIGHLIGHT_WIDTH: f64 = 2.5;
/// Datum outline, datum axis and coupling line widths, in output pixels.
const OUTLINE_WIDTH: f64 = 1.5;
const SHAFT_WIDTH: f64 = 2.5;
const COUPLING_WIDTH: f64 = 1.0;
/// Reference point and thermal dot sizes, in output pixels.
const POINT_SIZE: f64 = 8.0;
const THERMAL_SIZE: f64 = 6.0;

/// The colours symbols are drawn in, by kind.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlyphColors {
    /// Boundary conditions and rigid walls.
    pub constraint: Rgb,
    /// Forces, moments, pressures, body loads and initial velocities.
    pub load: Rgb,
    /// Reference points and couplings.
    pub interaction: Rgb,
    /// Contact and tie surfaces.
    pub contact: Rgb,
    /// Temperatures, heat fluxes, convection and heat sources.
    pub thermal: Rgb,
    /// Reference point names.
    pub label: Rgb,
    /// Datum planes and axes.
    pub datum: Rgb,
    /// The symbols of the definition selected in the browser.
    pub highlight: Rgb,
}

impl GlyphColors {
    /// The colours moved a fraction `k` (0..=1) of the way to `toward`: symbols that stay out
    /// of the way (over a results plot). The highlight stays as it is.
    pub fn muted(self, toward: Rgb, k: f64) -> GlyphColors {
        let k = if k.is_finite() { k.clamp(0.0, 1.0) } else { 0.0 };
        let mix = |c: Rgb| {
            let m = |a: u8, b: u8| (f64::from(a) + (f64::from(b) - f64::from(a)) * k).round().clamp(0.0, 255.0) as u8;
            Rgb(m(c.0, toward.0), m(c.1, toward.1), m(c.2, toward.2))
        };
        GlyphColors {
            constraint: mix(self.constraint),
            load: mix(self.load),
            interaction: mix(self.interaction),
            contact: mix(self.contact),
            thermal: mix(self.thermal),
            label: mix(self.label),
            datum: mix(self.datum),
            highlight: self.highlight,
        }
    }

    /// The colours over a results plot drawn on `background`: [`MUTED`] of the way to it.
    pub fn over_results(self, background: Rgb) -> GlyphColors {
        self.muted(background, MUTED)
    }
}

/// Which groups of symbols are drawn (the browser's eye toggles).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct GlyphShow {
    /// Boundary conditions, loads and initial velocities.
    pub loads: bool,
    /// Reference points, couplings, contact, ties and rigid walls.
    pub interactions: bool,
}

impl GlyphShow {
    pub const ALL: GlyphShow = GlyphShow { loads: true, interactions: true };

    /// Whether a symbol of the Interactions group (`interaction`), or else of Loads & BCs, is shown.
    fn allows(self, interaction: bool) -> bool {
        if interaction { self.interactions } else { self.loads }
    }
}

/// The kind of definition a symbol draws (with its name, what the browser selects).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Definition {
    BoundaryCondition,
    Load,
    ReferencePoint,
    Coupling,
    Contact,
    Tie,
    /// A step's initial velocities and rigid walls.
    Step,
}

/// How [`Glyphs::draw`] draws: what is shown, in which colours, at which size, and which
/// definition is highlighted.
#[derive(Clone, Copy, Debug)]
pub struct GlyphView<'a> {
    pub show: GlyphShow,
    pub colors: &'a GlyphColors,
    /// One symbol unit (a support triangle's length) in model length, from [`symbol_unit`];
    /// when it is not a positive finite number, the model's own unit
    /// ([`Glyphs::model_unit`]).
    pub unit: f64,
    /// The definition selected in the browser.
    pub highlight: Option<(Definition, &'a str)>,
}

/// What [`Glyphs::draw`] drew.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Drawn {
    /// The box around the symbols, for framing and clipping.
    pub bounds: Aabb3,
    /// Definitions drawn with a sample of their symbols: `Clamp: 400 of 12000 nodes`.
    pub notes: Vec<String>,
}

/// The symbol unit for symbols `size` pixels across where one pixel spans `length_per_pixel`
/// of model length ([`crate::Camera::length_per_pixel`]); the size is held to [`SIZE_RANGE`]
/// ([`DEFAULT_SIZE`] when it is not a number). `None` when there is no such unit.
pub fn symbol_unit(length_per_pixel: f64, size: f64) -> Option<f64> {
    let size = if size.is_finite() { size.clamp(SIZE_RANGE.0, SIZE_RANGE.1) } else { DEFAULT_SIZE };
    let unit = length_per_pixel * size;
    (unit.is_finite() && unit > 0.0).then_some(unit)
}

/// The symbol unit for symbols `size_pt` points across in a view `height_px` pixels high at
/// `pixels_per_point` (1 for a headless image, the display's scale for the desktop), seen with
/// `cam`: the one helper the desktop, `ui.render` and headless renders share, so a render at
/// device resolution matches the window. A scale that is not a positive number counts as 1.
pub fn unit_for(cam: &crate::Camera, height_px: f64, size_pt: f64, pixels_per_point: f64) -> Option<f64> {
    let ppp = if pixels_per_point.is_finite() && pixels_per_point > 0.0 { pixels_per_point } else { 1.0 };
    symbol_unit(cam.length_per_pixel(height_px)? * ppp, size_pt)
}

/// A symbol unit rounded to a quarter of a doubling (within 9 %), and the step it is: a viewer
/// rebuilds its symbols only when zooming changes the step.
pub fn quantise_unit(unit: f64) -> Option<(i32, f64)> {
    if !(unit.is_finite() && unit > 0.0) {
        return None;
    }
    let step = (unit.log2() * UNIT_STEPS).round().clamp(-4000.0, 4000.0) as i32;
    Some((step, (f64::from(step) / UNIT_STEPS).exp2()))
}

/// Where a symbol sits: a mesh node, or a reference point (its own node label when a result
/// moves it, else its coordinates).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Anchor {
    Node(u32),
    Point { label: u32, xyz: Vec3 },
}

/// The corners of a loaded face (or the ends of a loaded edge), with the outward direction it
/// had on the undeformed model (the solver's convention: a face's right-hand normal, outward
/// for solids; an in-plane edge's normal to the right of its direction).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Face {
    corners: [u32; 4],
    count: usize,
    outward: Vec3,
}

/// A place a node symbol may go: its anchor, where it is on the undeformed model (for
/// thinning), the outward direction of the set's faces there (when it has faces) and, for a
/// constraint, the degrees of freedom held there.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Spot {
    at: Anchor,
    p: Vec3,
    normal: Option<Vec3>,
    dofs: u8,
}

/// What a group of symbols draws at each of its places.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Shape {
    /// The constrained translations (bits 0..3) and rotations (bits 3..6) of each spot.
    Constraint,
    Force(Vec3),
    Moment(Vec3),
    /// Heat generated in a set's elements: a dot at each node.
    Thermal,
    /// A prescribed temperature (a boundary condition or a temperature field): a stem off the
    /// part ending in a bulb, a thermometer.
    Temperature,
    /// A support along a unit direction that is no axis (a symmetry plane's normal): a
    /// constraint triangle along it, from outside the part.
    Support(Vec3),
    /// A prescribed motion along a unit direction: an arrow leaving the node from a short bar,
    /// telling it apart from a force (an arrow onto the node) and a support (triangles).
    Prescribed(Vec3),
    /// Convection to the surroundings: a wavy line off the face.
    Convection,
    /// A heat flux into the face: an arrow onto it.
    HeatFlux,
    /// A pressure: onto the face (positive) or away from it.
    Pressure {
        onto: bool,
    },
    /// A traction along a unit direction.
    Traction(Vec3),
    /// A contact or tie surface.
    Surface,
}

/// The symbols of one definition at many places, thinned when drawn.
#[derive(Clone, Debug, PartialEq)]
struct Group {
    owner: usize,
    shape: Shape,
    /// Node places (or face places, for face shapes), all of them up to [`MAX_SCAN`].
    nodes: Vec<Spot>,
    faces: Vec<(Face, Vec3)>,
    /// How many places the definition has.
    total: usize,
    /// The thinning grid's cell, in symbol units.
    cell: f64,
    /// A load at one place, or anything at a reference point: concentrated, drawn over the part
    /// (#78).
    on_top: bool,
}

/// A symbol drawn once.
#[derive(Clone, Debug, PartialEq)]
enum Mark {
    /// A fixed arrow over the model (gravity, body forces): the `slot`th beside `base`.
    Free {
        base: Vec3,
        slot: usize,
        dir: Vec3,
    },
    /// An initial velocity: an arrow along `dir` ending short of a set centred at `centre`
    /// reaching `reach` toward it.
    Velocity {
        centre: Vec3,
        reach: f64,
        dir: Vec3,
    },
    /// A centrifugal load: its axis through `centre` and a curved arrow around it.
    Spin {
        centre: Vec3,
        axis: Vec3,
        half_length: f64,
    },
    ReferencePoint {
        at: Anchor,
        name: String,
    },
    CouplingLine {
        from: Anchor,
        to: u32,
    },
    /// A rigid wall: a square through `point` normal to `normal`, backing a set of half
    /// diagonal `reach` (none: a small square), hatched on the side away from `normal`.
    Wall {
        point: Vec3,
        normal: Vec3,
        reach: Option<f64>,
    },
}

impl Mark {
    /// Whether the symbol belongs to the Interactions group (else Loads & BCs).
    fn interaction(&self) -> bool {
        matches!(self, Mark::ReferencePoint { .. } | Mark::CouplingLine { .. } | Mark::Wall { .. })
    }
}

/// The symbols of a model's analysis setup, ready to draw.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Glyphs {
    /// The definitions drawn: (kind, name).
    owners: Vec<(Definition, String)>,
    groups: Vec<Group>,
    marks: Vec<(usize, Mark)>,
    /// The model's own symbol unit ([`UNIT_FRACTION`] of its diagonal).
    unit: f64,
}

/// Every element of a long list, strided down to about [`MAX_SCAN`] of them.
fn scan<T>(items: &[T]) -> impl Iterator<Item = &T> {
    items.iter().step_by(items.len().div_ceil(MAX_SCAN).max(1))
}

/// The corner node labels of an element face, in the face's order.
fn face_corners(model: &Model, f: ElementFace) -> Option<Vec<u32>> {
    let kind = model.element(f.element)?.kind;
    let local = kind.faces().get(usize::from(f.face).checked_sub(1)?)?;
    let corners = local.iter().filter(|&&i| i < kind.corner_count()).count();
    let mut labels = model.face_nodes(f)?;
    labels.truncate(corners);
    Some(labels)
}

/// A node's coordinates, when finite.
fn xyz(model: &Model, label: u32) -> Option<Vec3> {
    model.node(label).map(|n| Vec3::from(n.xyz)).filter(|p| p.is_finite())
}

/// Two unit vectors perpendicular to `a` and to each other.
fn perpendiculars(a: Vec3) -> (Vec3, Vec3) {
    let helper = if a.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let u = a.cross(helper).normalized().unwrap_or(Vec3::Y);
    let v = a.cross(u).normalized().unwrap_or(Vec3::Z);
    (u, v)
}

/// Thin `items` at their positions to one per grid cell of side `cell` (the first in each),
/// then stride them down to `cap`.
fn thin<T: Copy>(items: impl IntoIterator<Item = (T, Vec3)>, cell: f64, cap: usize) -> Vec<T> {
    let cell = if cell.is_finite() && cell > 0.0 { cell } else { 1.0 };
    let key = |v: f64| (v / cell).floor().clamp(-1e15, 1e15) as i64;
    let mut seen = HashSet::with_capacity(cap.saturating_mul(4).min(MAX_SCAN));
    let mut kept: Vec<T> = Vec::new();
    for (item, p) in items {
        if p.is_finite() && seen.insert((key(p.x), key(p.y), key(p.z))) {
            kept.push(item);
        }
    }
    if kept.len() > cap {
        let stride = kept.len().div_ceil(cap.max(1));
        kept = kept.into_iter().step_by(stride).take(cap).collect();
    }
    kept
}

/// What [`Glyphs::from_model`] reads, with the sets resolved once each.
struct Builder<'m> {
    model: &'m Model,
    sets: BTreeMap<String, ResolvedSet>,
    /// The nodes that can rotate, worked out when a constraint needs them.
    rotating: Option<BTreeSet<u32>>,
    bounds: Aabb3,
    owners: Vec<(Definition, String)>,
    groups: Vec<Group>,
    marks: Vec<(usize, Mark)>,
    /// Free arrows placed so far (each beside the last).
    free: usize,
}

impl Builder<'_> {
    fn set(&mut self, name: &str) -> Option<&ResolvedSet> {
        if !self.sets.contains_key(name) {
            let r = self.model.resolve_set_partial(name)?;
            self.sets.insert(name.to_string(), r);
        }
        self.sets.get(name)
    }

    fn owner(&mut self, def: Definition, name: &str) -> usize {
        self.owners.push((def, name.to_string()));
        self.owners.len() - 1
    }

    fn group(&mut self, owner: usize, shape: Shape, nodes: Vec<Spot>, faces: Vec<(Face, Vec3)>, total: usize, cell: f64) {
        if !(nodes.is_empty() && faces.is_empty()) {
            let concentrated =
                matches!(shape, Shape::Force(_) | Shape::Moment(_)) || matches!(nodes.first(), Some(Spot { at: Anchor::Point { .. }, .. }));
            let on_top = total == 1 && faces.is_empty() && concentrated;
            self.groups.push(Group { owner, shape, nodes, faces, total, cell, on_top });
        }
    }

    /// A set's nodes (a sample of a huge one) with their coordinates, and how many it has.
    fn node_points(&mut self, set: &str) -> (Vec<(u32, Vec3)>, usize) {
        let model = self.model;
        let Some(r) = self.set(set) else { return (Vec::new(), 0) };
        (scan(&r.nodes).filter_map(|&n| Some((n, xyz(model, n)?))).collect(), r.nodes.len())
    }

    /// The places of a target (a sample of a huge set) and how many it has; a set's nodes carry
    /// the outward direction of its faces there, when it has faces.
    fn spots(&mut self, target: &Target) -> (Vec<Spot>, usize) {
        match target {
            Target::ReferencePoint(r) => {
                let spots: Vec<Spot> = self
                    .model
                    .reference_point(r)
                    .map(|p| {
                        let xyz = Vec3::from(p.xyz);
                        Spot { at: Anchor::Point { label: p.label, xyz }, p: xyz, normal: None, dofs: 0 }
                    })
                    .into_iter()
                    .collect();
                let n = spots.len();
                (spots, n)
            }
            Target::Set(s) => {
                let model = self.model;
                let (at, total) = self.node_points(s);
                let mut normals: BTreeMap<u32, Vec3> = BTreeMap::new();
                if let Some(r) = self.set(s) {
                    for &f in scan(&r.faces) {
                        let Some(labels) = face_corners(model, f) else { continue };
                        let Some(p) = labels.iter().map(|&n| xyz(model, n)).collect::<Option<Vec<_>>>() else { continue };
                        let Some(n) = polygon_normal(&p) else { continue };
                        for l in labels {
                            let sum = normals.entry(l).or_insert(Vec3::ZERO);
                            *sum = *sum + n;
                        }
                    }
                }
                let spots = at
                    .into_iter()
                    .map(|(n, p)| Spot { at: Anchor::Node(n), p, normal: normals.get(&n).and_then(|v| v.normalized()), dofs: 0 })
                    .collect();
                (spots, total)
            }
        }
    }

    /// A thermal condition's symbol (`shape`, at nodes) at its target's nodes.
    fn thermal(&mut self, owner: usize, target: &Target, shape: Shape) {
        let (spots, total) = self.spots(target);
        let cell = if shape == Shape::Temperature { 2.4 } else { 1.5 };
        self.group(owner, shape, spots, Vec::new(), total, cell);
    }

    /// The faces of a set (its element faces, else the edges of its 2D elements) with their
    /// centres, and how many it has.
    fn faces(&mut self, set: &str) -> (Vec<(Face, Vec3)>, usize) {
        let model = self.model;
        let Some(r) = self.set(set) else { return (Vec::new(), 0) };
        let face = |labels: &[u32], outward: Vec3| {
            let mut corners = [0u32; 4];
            for (slot, l) in corners.iter_mut().zip(labels) {
                *slot = *l;
            }
            Face { corners, count: labels.len().min(4), outward }
        };
        let corners_at = |labels: &[u32]| -> Option<Vec<Vec3>> { labels.iter().map(|&n| xyz(model, n)).collect() };
        let mut found: Vec<(Face, Vec3)> = Vec::new();
        for &f in scan(&r.faces) {
            let Some(labels) = face_corners(model, f) else { continue };
            let Some(p) = corners_at(&labels) else { continue };
            if let Some(n) = polygon_normal(&p) {
                found.push((face(&labels, n), centre(&p)));
            }
        }
        if !found.is_empty() {
            return (found, r.faces.len());
        }
        // A 2D model's loaded edges (in the XY plane).
        for &g in scan(&r.edges) {
            if model.element(g.element).is_none_or(|e| e.kind.dimension() != 2) {
                continue;
            }
            let Some(labels) = model.edge_nodes(g) else { continue };
            let labels = labels.get(..2).unwrap_or(&[]);
            let Some(p) = corners_at(labels) else { continue };
            let (Some(a), Some(b)) = (p.first(), p.get(1)) else { continue };
            let t = *b - *a;
            if let Some(n) = Vec3::new(t.y, -t.x, 0.0).normalized() {
                found.push((face(labels, n), centre(&p)));
            }
        }
        (found, r.edges.len())
    }

    /// The bounding box of a set's nodes.
    fn set_bounds(&mut self, set: &str) -> Aabb3 {
        let mut b = Aabb3::EMPTY;
        for (_, p) in self.node_points(set).0 {
            b.add(p);
        }
        b
    }

    /// A fixed arrow over the model (gravity, body forces).
    fn free_arrow(&mut self, owner: usize, dir: Vec3) {
        let Some(dir) = dir.normalized() else { return };
        let b = self.bounds;
        // Over the model's middle, in front of it, each one beside the last.
        let base = Vec3::new(b.center().x, b.min.y, b.max.z);
        self.marks.push((owner, Mark::Free { base, slot: self.free, dir }));
        self.free += 1;
    }

    fn boundary_condition(&mut self, name: &str, bc: &BoundaryCondition) {
        let owner = self.owner(Definition::BoundaryCondition, name);
        if bc.is_thermal() {
            return self.thermal(owner, bc.target(), Shape::Temperature);
        }
        let mut dofs = bc.dofs().bits();
        // A prescribed non-zero translation moves the nodes: an arrow along the motion; the
        // degrees of freedom it drives are no supports.
        if let BoundaryCondition::Displacement { values, .. } = bc {
            let moved = |i: usize| values.get(i).copied().flatten().filter(|v| v.is_finite() && *v != 0.0);
            let motion = Vec3::new(moved(0).unwrap_or(0.0), moved(1).unwrap_or(0.0), moved(2).unwrap_or(0.0));
            for (i, _) in values.iter().enumerate().filter(|(i, _)| moved(*i).is_some()) {
                dofs &= !(1u8 << i);
            }
            if let Some(dir) = motion.normalized() {
                let (spots, total) = self.spots(bc.target());
                self.group(owner, Shape::Prescribed(dir), spots, Vec::new(), total, 2.8);
            }
        }
        // A symmetry plane shows its normal only: the rotations it holds follow from it.
        if let BoundaryCondition::Symmetry { plane, .. } = bc {
            if let Some(n) = plane.inclined().and_then(|n| Vec3::from(n).normalized()) {
                let (spots, total) = self.spots(bc.target());
                return self.group(owner, Shape::Support(n), spots, Vec::new(), total, 2.4);
            }
            dofs &= 0b000_111;
        }
        if dofs == 0 {
            return;
        }
        let (spots, total) = self.spots(bc.target());
        if dofs & 0b111_000 != 0 && self.rotating.is_none() {
            self.rotating = Some(self.model.rotating_nodes());
        }
        let none = BTreeSet::new();
        let rotating = self.rotating.as_ref().unwrap_or(&none);
        let spots = spots
            .into_iter()
            .filter_map(|s| {
                let rotates = match s.at {
                    Anchor::Node(n) => rotating.contains(&n),
                    Anchor::Point { .. } => true,
                };
                let dofs = if rotates { dofs } else { dofs & 0b000_111 };
                (dofs != 0).then_some(Spot { dofs, ..s })
            })
            .collect();
        self.group(owner, Shape::Constraint, spots, Vec::new(), total, 2.4);
    }

    fn load(&mut self, name: &str, load: &Load) {
        let owner = self.owner(Definition::Load, name);
        match load {
            Load::Force { target, vector, .. } | Load::Moment { target, vector, .. } => {
                let Some(dir) = Vec3::from(*vector).normalized() else { return };
                let (shape, cell) = if matches!(load, Load::Moment { .. }) { (Shape::Moment(dir), 3.5) } else { (Shape::Force(dir), 2.8) };
                let (spots, total) = self.spots(target);
                self.group(owner, shape, spots, Vec::new(), total, cell);
            }
            Load::Pressure { set, magnitude } => {
                if !magnitude.is_finite() || *magnitude == 0.0 {
                    return;
                }
                let (faces, total) = self.faces(set);
                self.group(owner, Shape::Pressure { onto: *magnitude > 0.0 }, Vec::new(), faces, total, 1.8);
            }
            Load::Traction { set, vector } => {
                let Some(dir) = Vec3::from(*vector).normalized() else { return };
                let (faces, total) = self.faces(set);
                self.group(owner, Shape::Traction(dir), Vec::new(), faces, total, 1.8);
            }
            Load::Gravity { acceleration, .. } => self.free_arrow(owner, Vec3::from(*acceleration)),
            Load::BodyForce { vector, .. } => self.free_arrow(owner, Vec3::from(*vector)),
            Load::Centrifugal { point, axis, .. } => {
                let (Some(axis), point) = (Vec3::from(*axis).normalized(), Vec3::from(*point)) else { return };
                if !point.is_finite() {
                    return;
                }
                // The axis through the model's middle, as long as the model.
                let c = self.bounds.center();
                let centre = point + axis * axis.dot(c - point);
                self.marks.push((owner, Mark::Spin { centre, axis, half_length: 0.5 * self.bounds.diagonal() }));
            }
            Load::Temperature { target, .. } => self.thermal(owner, target, Shape::Temperature),
            Load::Convection { set, .. } | Load::SurfaceFlux { set, .. } => {
                let shape = if matches!(load, Load::Convection { .. }) { Shape::Convection } else { Shape::HeatFlux };
                let (faces, total) = self.faces(set);
                if faces.is_empty() {
                    self.thermal(owner, &Target::Set(set.clone()), Shape::Thermal);
                } else {
                    self.group(owner, shape, Vec::new(), faces, total, 2.2);
                }
            }
            Load::HeatSource { set, .. } => self.thermal(owner, &Target::Set(set.clone()), Shape::Thermal),
        }
    }

    /// Squares on the faces of contact and tie surfaces.
    fn surfaces(&mut self) {
        let model = self.model;
        let contacts = model.contacts().iter().map(|(name, c)| {
            let sets: Vec<String> = match &c.surfaces {
                ContactSurfaces::Pair { main, secondary } => vec![main.clone(), secondary.clone()],
                ContactSurfaces::General { surface } => surface.iter().cloned().collect(),
            };
            (Definition::Contact, name, sets)
        });
        let ties = model.ties().iter().map(|(name, t)| (Definition::Tie, name, vec![t.main.clone(), t.secondary.clone()]));
        for (def, name, sets) in contacts.chain(ties).collect::<Vec<_>>() {
            let owner = self.owner(def, name);
            for set in sets {
                let (faces, total) = self.faces(&set);
                self.group(owner, Shape::Surface, Vec::new(), faces, total, 2.5);
            }
        }
    }

    fn interactions(&mut self, unit: f64) {
        let model = self.model;
        for (name, p) in model.reference_points() {
            let xyz = Vec3::from(p.xyz);
            if xyz.is_finite() {
                let owner = self.owner(Definition::ReferencePoint, name);
                self.marks.push((owner, Mark::ReferencePoint { at: Anchor::Point { label: p.label, xyz }, name: name.clone() }));
            }
        }
        for (name, c) in model.couplings() {
            let Some(p) = model.reference_point(&c.reference_point) else { continue };
            let owner = self.owner(Definition::Coupling, name);
            let from = Anchor::Point { label: p.label, xyz: Vec3::from(p.xyz) };
            let (at, _) = self.node_points(&c.set);
            // Spread over the set: a grid of about four cells across its extent.
            let mut b = Aabb3::EMPTY;
            for (_, q) in &at {
                b.add(*q);
            }
            let cell = (b.diagonal() / 4.0).max(unit);
            for to in thin(at, cell, MAX_COUPLING_LINES) {
                self.marks.push((owner, Mark::CouplingLine { from, to }));
            }
        }
    }

    fn steps(&mut self) {
        let model = self.model;
        for step in model.steps() {
            let Procedure::ExplicitDynamic(e) = &step.procedure else { continue };
            let owner = self.owner(Definition::Step, &step.name);
            for v in &e.initial_velocities {
                let Some(dir) = Vec3::from(v.velocity).normalized() else { continue };
                let b = self.set_bounds(&v.set);
                if b.is_empty() {
                    continue;
                }
                // Behind the set, so it clears what the set moves toward (a rigid wall).
                let half = (b.max - b.min) * 0.5;
                let reach = (dir.x * half.x).abs() + (dir.y * half.y).abs() + (dir.z * half.z).abs();
                self.marks.push((owner, Mark::Velocity { centre: b.center(), reach, dir }));
            }
            for w in &e.rigid_walls {
                let (Some(normal), point) = (Vec3::from(w.normal).normalized(), Vec3::from(w.point)) else { continue };
                if !point.is_finite() {
                    continue;
                }
                // Behind the set's nodes and wide enough to back them.
                let b = self.set_bounds(&w.set);
                let (point, reach) = if b.is_empty() {
                    (point, None)
                } else {
                    let c = b.center();
                    (c - normal * normal.dot(c - point), Some(0.5 * b.diagonal()))
                };
                self.marks.push((owner, Mark::Wall { point, normal, reach }));
            }
        }
    }
}

impl Glyphs {
    /// The symbols of `model`'s boundary conditions, loads, contact and tie surfaces, reference
    /// points, couplings and explicit steps' initial velocities and rigid walls.
    pub fn from_model(model: &Model) -> Glyphs {
        let mut bounds = Aabb3::EMPTY;
        for n in model.nodes().values() {
            bounds.add(Vec3::from(n.xyz));
        }
        for p in model.reference_points().values() {
            bounds.add(Vec3::from(p.xyz));
        }
        let d = bounds.diagonal();
        let unit = if d.is_finite() && d > 0.0 { d * UNIT_FRACTION } else { 1.0 };
        let mut b =
            Builder { model, sets: BTreeMap::new(), rotating: None, bounds, owners: Vec::new(), groups: Vec::new(), marks: Vec::new(), free: 0 };
        for (name, bc) in model.boundary_conditions() {
            b.boundary_condition(name, bc);
        }
        for (name, load) in model.loads() {
            b.load(name, load);
        }
        b.steps();
        b.surfaces();
        b.interactions(unit);
        b.marks.truncate(MAX_SYMBOLS);
        Glyphs { owners: b.owners, groups: b.groups, marks: b.marks, unit }
    }

    /// [`Glyphs::from_model`] without the definitions `hidden` names (the view's hidden items).
    pub fn from_model_hiding(model: &Model, hidden: &dyn Fn(Definition, &str) -> bool) -> Glyphs {
        let mut g = Glyphs::from_model(model);
        let gone: Vec<bool> = g.owners.iter().map(|(d, n)| hidden(*d, n)).collect();
        let kept = |owner: usize| !gone.get(owner).copied().unwrap_or(false);
        g.groups.retain(|x| kept(x.owner));
        g.marks.retain(|(owner, _)| kept(*owner));
        g
    }

    /// The model's own symbol unit (a fortieth of its bounding box diagonal): the size symbols
    /// are drawn at when no camera gives one.
    pub fn model_unit(&self) -> f64 {
        self.unit
    }

    /// Draw the symbols into `scene` as `view` says, at the nodes' positions `pos(label)` (the
    /// undeformed model, or a results plot's deformed shape); a symbol at a node `pos` does not
    /// know is left out (a reference point falls back to its own coordinates).
    pub fn draw(&self, pos: impl Fn(u32) -> Option<Vec3>, view: &GlyphView, scene: &mut Scene) -> Drawn {
        let at = |a: Anchor| -> Option<Vec3> {
            match a {
                Anchor::Node(n) => pos(n),
                Anchor::Point { label, xyz } => Some(pos(label).unwrap_or(xyz)),
            }
            .filter(|p| p.is_finite())
        };
        let u = if view.unit.is_finite() && view.unit > 0.0 { view.unit } else { self.unit };
        let colors = view.colors;
        let mut d = Draw { scene, bounds: Aabb3::EMPTY, color: colors.load, width: LINE_WIDTH, on_top: false };
        let mut notes = Vec::new();
        let mut drawn = 0usize;
        // Whether an owner is the highlighted definition.
        let selected = |owner: usize| view.highlight.is_some_and(|(def, name)| self.owners.get(owner).is_some_and(|(d, n)| *d == def && n == name));
        for g in &self.groups {
            if !view.show.allows(g.shape == Shape::Surface) {
                continue;
            }
            let color = match g.shape {
                Shape::Constraint | Shape::Support(_) | Shape::Prescribed(_) => colors.constraint,
                Shape::Thermal | Shape::Temperature | Shape::Convection | Shape::HeatFlux => colors.thermal,
                Shape::Surface => colors.contact,
                _ => colors.load,
            };
            d.style(color, selected(g.owner), colors.highlight);
            let room = MAX_SYMBOLS.saturating_sub(drawn).min(MAX_PER_DEFINITION);
            let cell = g.cell * u;
            let kept = if g.faces.is_empty() {
                let spots = thin(g.nodes.iter().map(|s| (*s, s.p)), cell, room);
                d.on_top = g.on_top;
                for s in &spots {
                    let Some(p) = at(s.at) else { continue };
                    d.node_symbol(g.shape, p, s, u, g.total > 1);
                }
                d.on_top = false;
                spots.len()
            } else {
                let faces = thin(g.faces.iter().copied(), cell, room);
                for f in &faces {
                    d.face_symbol(g.shape, f, &pos, u);
                }
                faces.len()
            };
            drawn += kept;
            if kept < g.total && g.total > MAX_PER_DEFINITION {
                let name = self.owners.get(g.owner).map_or("?", |(_, n)| n.as_str());
                let what = if g.faces.is_empty() { "nodes" } else { "faces" };
                notes.push(format!("{name}: {kept} of {} {what}", g.total));
            }
        }
        for (owner, m) in &self.marks {
            if !view.show.allows(m.interaction()) || drawn >= MAX_SYMBOLS {
                continue;
            }
            drawn += 1;
            let color = match m {
                Mark::Wall { .. } => colors.constraint,
                Mark::ReferencePoint { .. } | Mark::CouplingLine { .. } => colors.interaction,
                _ => colors.load,
            };
            let picked = selected(*owner);
            d.style(color, picked, colors.highlight);
            match m {
                Mark::Free { base, slot, dir } => {
                    let middle = *base + Vec3::new(2.5 * *slot as f64, -4.0, 2.0) * u;
                    // The arrow's middle sits at the spot, whichever way it points.
                    d.arrow(middle - *dir * (2.5 * u), *dir, 5.0 * u, 1.2 * u);
                }
                Mark::Velocity { centre, reach, dir } => d.arrow(*centre - *dir * (reach + 6.0 * u), *dir, 5.0 * u, 1.2 * u),
                Mark::Spin { centre: c, axis, half_length } => {
                    let half = half_length.max(4.0 * u);
                    d.line(*c - *axis * half, *c + *axis * half);
                    d.curl(*c + *axis * (half * 0.8), *axis, 2.0 * u, u);
                }
                Mark::ReferencePoint { at: a, name } => {
                    let Some(p) = at(*a) else { continue };
                    d.on_top = true;
                    for axis in [Vec3::X, Vec3::Y, Vec3::Z] {
                        d.line(p - axis * (0.9 * u), p + axis * (0.9 * u));
                    }
                    // A dot: a zero-length line drawn round.
                    let width = d.width;
                    d.width = POINT_SIZE;
                    d.line(p, p);
                    (d.width, d.on_top) = (width, false);
                    let label = if picked { colors.highlight } else { colors.label };
                    d.scene.labels.push(SceneLabel { p, text: name.clone(), color: label });
                }
                Mark::CouplingLine { from, to } => {
                    if let (Some(a), Some(b)) = (at(*from), pos(*to)) {
                        if !picked {
                            d.width = COUPLING_WIDTH;
                        }
                        d.line(a, b);
                    }
                }
                Mark::Wall { point, normal, reach } => {
                    let half = reach.map_or(4.0 * u, |r| (r + u).clamp(2.0 * u, 40.0 * u));
                    d.wall(*point, *normal, half, u, picked);
                }
            }
        }
        Drawn { bounds: d.bounds, notes }
    }
}

/// The margin around a datum plane's rectangle and beyond a datum axis's ends, as a fraction
/// of the model's extent drawn.
const DATUM_MARGIN: f64 = 0.05;
/// A datum point's cross: each arm's length, as a fraction of the bounding box diagonal.
const DATUM_POINT: f64 = 0.02;
/// The smallest diagonal a datum is sized against, so an empty model still draws one.
const MIN_DIAGONAL: f64 = 1e-9; // unit-factors: a floor in metres (SI) against a zero-size box
/// A datum plane's fill: how opaque its translucent rectangle is.
pub const DATUM_FILL_ALPHA: f64 = 0.12;

/// Corner `k` (0..8, one bit per axis) of box `b`.
fn box_corner(b: &Aabb3, k: usize) -> Vec3 {
    let pick = |bit: usize, lo: f64, hi: f64| if k & bit == 0 { lo } else { hi };
    Vec3::new(pick(1, b.min.x, b.max.x), pick(2, b.min.y, b.max.y), pick(4, b.min.z, b.max.z))
}

/// The rectangle datum plane (`p`, unit normal `n`) is drawn as, in its own axes `(u, v)`
/// ([`perpendiculars`]): the model's extent within the plane (where the plane cuts `bounds`;
/// the box's shadow on the plane when it misses it) with a [`DATUM_MARGIN`] around it. Returns
/// the rectangle's centre on the plane and its half sizes along `u` and `v`.
fn plane_rect(bounds: &Aabb3, p: Vec3, n: Vec3, u: Vec3, v: Vec3) -> (Vec3, f64, f64) {
    if bounds.is_empty() || !bounds.diagonal().is_finite() {
        return (p, 1.0, 1.0);
    }
    let (lo, hi) = (bounds.min, bounds.max);
    let corner =
        |k: usize| Vec3::new(if k & 1 == 0 { lo.x } else { hi.x }, if k & 2 == 0 { lo.y } else { hi.y }, if k & 4 == 0 { lo.z } else { hi.z });
    let height = |q: Vec3| (q - p).dot(n);
    // The box's twelve edges join corners differing in one coordinate bit.
    let mut cut: Vec<Vec3> = Vec::with_capacity(12);
    for a in 0..8usize {
        for bit in [1usize, 2, 4].into_iter().filter(|bit| a & bit == 0) {
            let b = a | bit;
            let (qa, qb) = (corner(a), corner(b));
            let (ha, hb) = (height(qa), height(qb));
            if (ha <= 0.0 && hb >= 0.0) || (ha >= 0.0 && hb <= 0.0) {
                let t = if (ha - hb).abs() > 0.0 { ha / (ha - hb) } else { 0.0 };
                cut.push(qa + (qb - qa) * t);
            }
        }
    }
    if cut.is_empty() {
        cut = (0..8).map(|k| corner(k) - n * height(corner(k))).collect();
    }
    let (mut u0, mut u1, mut v0, mut v1) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for q in &cut {
        let (a, b) = ((*q - p).dot(u), (*q - p).dot(v));
        (u0, u1, v0, v1) = (u0.min(a), u1.max(a), v0.min(b), v1.max(b));
    }
    let size = (u1 - u0).max(v1 - v0);
    let margin = if size > 0.0 { DATUM_MARGIN * size } else { DATUM_MARGIN * bounds.diagonal().max(MIN_DIAGONAL) };
    let centre = p + u * ((u0 + u1) * 0.5) + v * ((v0 + v1) * 0.5);
    (centre, (u1 - u0) * 0.5 + margin, (v1 - v0) * 0.5 + margin)
}

/// The ends datum axis (`p`, unit direction `a`) is drawn between: the model's extent along it
/// (`bounds`' corners projected on it) with a [`DATUM_MARGIN`] beyond each end.
fn axis_ends(bounds: &Aabb3, p: Vec3, a: Vec3) -> (Vec3, Vec3) {
    if bounds.is_empty() || !bounds.diagonal().is_finite() {
        return (p - a, p + a);
    }
    let (mut t0, mut t1) = (f64::INFINITY, f64::NEG_INFINITY);
    for k in 0..8usize {
        let t = (box_corner(bounds, k) - p).dot(a);
        (t0, t1) = (t0.min(t), t1.max(t));
    }
    let span = t1 - t0;
    let margin = if span > 0.0 { DATUM_MARGIN * span } else { DATUM_MARGIN * bounds.diagonal().max(MIN_DIAGONAL) };
    (p + a * (t0 - margin), p + a * (t1 + margin))
}

/// Draw the model's datum planes (an outlined rectangle with a light translucent fill: the model's extent
/// within the plane, `bounds`, plus a margin), axes (a line along the model's extent) and points
/// (a small cross), each with its name; returns the box around them. Only the datums `shown`
/// takes are drawn, and `highlight` (the one selected in the browser) whether shown or not, in
/// the highlight colour.
pub fn datums(
    model: &Model,
    bounds: &Aabb3,
    shown: &dyn Fn(&str) -> bool,
    highlight: Option<&str>,
    colors: &GlyphColors,
    scene: &mut Scene,
) -> Aabb3 {
    let mut d = Draw { scene, bounds: Aabb3::EMPTY, color: colors.datum, width: OUTLINE_WIDTH, on_top: false };
    let diagonal = if bounds.is_empty() { 1.0 } else { bounds.diagonal().max(MIN_DIAGONAL) };
    // A suppressed datum is out of the model: drawn only while selected in the browser.
    let drawn = |n: &String| (shown(n) && !model.datum_suppressed(n)) || highlight == Some(n.as_str());
    for (name, datum) in model.datums().iter().filter(|(n, _)| drawn(n)) {
        let picked = highlight == Some(name.as_str());
        let color = if picked { colors.highlight } else { colors.datum };
        let outline = if picked { HIGHLIGHT_WIDTH } else { OUTLINE_WIDTH };
        // The one selected in the browser shows through the part.
        (d.color, d.on_top) = (color, picked);
        match *datum {
            fieldcraft_model::Datum::Plane { point, normal } => {
                let (p, n) = (Vec3::from(point), Vec3::from(normal));
                let (u, v) = perpendiculars(n);
                let (c, hu, hv) = plane_rect(bounds, p, n, u, v);
                let corner = |a: f64, b: f64| c + u * (a * hu) + v * (b * hv);
                let ring = [corner(-1.0, -1.0), corner(1.0, -1.0), corner(1.0, 1.0), corner(-1.0, 1.0)];
                d.width = outline;
                for (k, a) in ring.iter().enumerate() {
                    let b = ring.get((k + 1) % 4).copied().unwrap_or(*a);
                    d.line(*a, b);
                }
                let [a, b, c2, e] = ring;
                for p in [[a, b, c2], [a, c2, e]] {
                    d.scene.fills.push(SceneFill { tri: SceneTri { p, normal: n, color }, alpha: DATUM_FILL_ALPHA });
                }
                d.scene.labels.push(SceneLabel { p: corner(1.0, 1.0), text: name.clone(), color });
            }
            fieldcraft_model::Datum::Axis { point, direction } => {
                let (a, b) = axis_ends(bounds, Vec3::from(point), Vec3::from(direction));
                d.width = if picked { SHAFT_WIDTH + 1.0 } else { SHAFT_WIDTH };
                d.line(a, b);
                d.scene.labels.push(SceneLabel { p: b, text: name.clone(), color });
            }
            fieldcraft_model::Datum::Point { point } => {
                let p = Vec3::from(point);
                let arm = DATUM_POINT * diagonal;
                d.width = outline;
                for axis in [Vec3::X, Vec3::Y, Vec3::Z] {
                    d.line(p - axis * arm, p + axis * arm);
                }
                d.scene.labels.push(SceneLabel { p: p + (Vec3::X + Vec3::Y) * arm, text: name.clone(), color });
            }
        }
    }
    d.bounds
}

/// Symbol primitives into a scene, in the current colour and width, depth-tested or on top.
struct Draw<'s> {
    scene: &'s mut Scene,
    bounds: Aabb3,
    color: Rgb,
    width: f64,
    on_top: bool,
}

impl Draw<'_> {
    /// Draw in `color` at the line width, or in `highlight`, wider, when `selected`.
    fn style(&mut self, color: Rgb, selected: bool, highlight: Rgb) {
        (self.color, self.width) = if selected { (highlight, HIGHLIGHT_WIDTH) } else { (color, LINE_WIDTH) };
    }

    fn point(&mut self, p: Vec3, size: f64) {
        self.bounds.add(p);
        self.scene.points.push(ScenePoint { p, color: self.color, size });
    }

    fn line(&mut self, a: Vec3, b: Vec3) {
        if a.is_finite() && b.is_finite() {
            self.bounds.add(a);
            self.bounds.add(b);
            let l = SceneLine { a, b, color: self.color, width: self.width };
            if self.on_top { self.scene.on_top.push(l) } else { self.scene.lines.push(l) }
        }
    }

    /// Two crossed wire triangles with the tip at `tip`, pointing along unit `dir`, `len` long and
    /// `half` across at the base, in the planes of `dir` with `e1` and with `e2`: a triangle seen
    /// from any side.
    fn triangle(&mut self, tip: Vec3, dir: Vec3, len: f64, half: f64, e1: Vec3, e2: Vec3) {
        let base = tip - dir * len;
        for e in [e1, e2] {
            let (a, b) = (base + e * half, base - e * half);
            self.line(tip, a);
            self.line(tip, b);
            self.line(a, b);
        }
    }

    /// An open arrowhead: four lines back from `tip` against unit `dir`, `len` long.
    fn head(&mut self, tip: Vec3, dir: Vec3, len: f64) {
        let (e1, e2) = perpendiculars(dir);
        let back = tip - dir * len;
        for e in [e1, -e1, e2, -e2] {
            self.line(tip, back + e * (0.4 * len));
        }
    }

    /// An arrow from `tail` along unit `dir`, `len` long, its open head `head` long.
    fn arrow(&mut self, tail: Vec3, dir: Vec3, len: f64, head: f64) {
        let tip = tail + dir * len;
        self.line(tail, tip);
        self.head(tip, dir, head);
    }

    /// A curved arrow three quarters of the way around unit `axis` at `centre`, `radius` across,
    /// turning right-handed about the axis, its head `head` long and tangent at its end.
    fn curl(&mut self, centre: Vec3, axis: Vec3, radius: f64, head: f64) {
        let (u, v) = perpendiculars(axis);
        let at = |a: f64| centre + (u * a.cos() + v * a.sin()) * radius;
        let end = 1.5 * std::f64::consts::PI;
        for k in 0..ARC_SEGMENTS {
            let (a, b) = (end * k as f64 / ARC_SEGMENTS as f64, end * (k + 1) as f64 / ARC_SEGMENTS as f64);
            self.line(at(a), at(b));
        }
        // u × v = axis, so increasing angle turns right-handed about it.
        let tangent = (u * -end.sin() + v * end.cos()).normalized().unwrap_or(u);
        self.head(at(end), tangent, head);
    }

    /// The symbol of a node group's shape at `p` (`spread`: one of many, drawn smaller).
    fn node_symbol(&mut self, shape: Shape, p: Vec3, s: &Spot, u: f64, spread: bool) {
        match shape {
            Shape::Constraint => {
                let axes = [Vec3::X, Vec3::Y, Vec3::Z];
                for (i, axis) in axes.into_iter().enumerate() {
                    let (t, r) = (s.dofs & (1 << i) != 0, s.dofs & (1 << (i + 3)) != 0);
                    if !(t || r) {
                        continue;
                    }
                    // From outside the part: pointing against the faces' outward direction
                    // (along the axis where it lies in the faces).
                    let dir = match s.normal {
                        Some(n) if n.dot(axis) > 0.1 => -axis,
                        _ => axis,
                    };
                    let (e1, e2) = (axes[(i + 1) % 3], axes[(i + 2) % 3]);
                    if t {
                        self.triangle(p, dir, u, 0.45 * u, e1, e2);
                    }
                    if r {
                        // The double mark: a smaller triangle behind the translation's.
                        self.triangle(p - dir * u, dir, 0.8 * u, 0.36 * u, e1, e2);
                    }
                }
            }
            Shape::Support(n) => {
                let dir = match s.normal {
                    Some(f) if f.dot(n) > 0.1 => -n,
                    _ => n,
                };
                let (e1, e2) = perpendiculars(n);
                self.triangle(p, dir, u, 0.45 * u, e1, e2);
            }
            Shape::Force(dir) => {
                let (len, head) = if spread { (2.5 * u, 0.8 * u) } else { (4.0 * u, u) };
                self.arrow(p - dir * len, dir, len, head);
            }
            Shape::Moment(dir) => {
                let (radius, head) = if spread { (0.8 * u, 0.6 * u) } else { (1.3 * u, 0.8 * u) };
                self.line(p - dir * (1.6 * radius), p + dir * (1.6 * radius));
                self.curl(p, dir, radius, head);
            }
            Shape::Thermal => self.point(p, THERMAL_SIZE),
            Shape::Temperature => {
                // Off the part along its faces' outward direction (up, away from faces).
                let dir = s.normal.unwrap_or(Vec3::Z);
                let (stem, bulb) = if spread { (1.2 * u, 0.32 * u) } else { (1.8 * u, 0.45 * u) };
                self.line(p, p + dir * stem);
                let (e1, _) = perpendiculars(dir);
                let c = p + dir * (stem + bulb);
                for k in 0..ARC_SEGMENTS {
                    let at = |k: usize| {
                        let a = std::f64::consts::TAU * k as f64 / ARC_SEGMENTS as f64;
                        c + (dir * a.cos() + e1 * a.sin()) * bulb
                    };
                    self.line(at(k), at(k + 1));
                }
                self.point(p, THERMAL_SIZE * 0.7);
            }
            Shape::Prescribed(dir) => {
                let (len, head) = if spread { (2.2 * u, 0.7 * u) } else { (3.5 * u, 0.9 * u) };
                let (e1, _) = perpendiculars(dir);
                self.line(p - e1 * (0.45 * u), p + e1 * (0.45 * u));
                self.arrow(p, dir, len, head);
            }
            Shape::Pressure { .. } | Shape::Traction(_) | Shape::Surface | Shape::Convection | Shape::HeatFlux => {}
        }
    }

    /// The symbol of a face group's shape on `face` (at its corners' positions `pos`).
    fn face_symbol(&mut self, shape: Shape, face: &Face, pos: &impl Fn(u32) -> Option<Vec3>, u: f64) {
        let mut p = [Vec3::ZERO; 4];
        let corners = face.corners.iter().take(face.count).map(|&n| pos(n));
        if face.count == 0 || corners.zip(p.iter_mut()).any(|(q, slot)| q.map(|q| *slot = q).is_none()) {
            return;
        }
        let p = p.get(..face.count).unwrap_or(&[]);
        let c = centre(p);
        // The face's own normal now (it may be deformed), on the side it pointed.
        let outward = if face.count >= 3 { polygon_normal(p).map(|n| if n.dot(face.outward) < 0.0 { -n } else { n }) } else { None };
        let outward = outward.unwrap_or(face.outward);
        let (len, head) = (2.2 * u, 0.7 * u);
        match shape {
            Shape::Pressure { onto: true } => self.arrow(c + outward * len, -outward, len, head),
            Shape::Pressure { onto: false } => self.arrow(c, outward, len, head),
            Shape::Traction(dir) => self.arrow(c - dir * len, dir, len, head),
            Shape::HeatFlux => self.arrow(c + outward * len, -outward, len, head),
            Shape::Convection => {
                // A wave leaving the face.
                let (e1, _) = perpendiculars(outward);
                let steps = 6;
                let at = |k: usize| {
                    let along = len * k as f64 / steps as f64;
                    let side = if k == 0 { 0.0 } else { (k as f64 * std::f64::consts::FRAC_PI_2).sin() * 0.25 * u };
                    c + outward * along + e1 * side
                };
                for k in 0..steps {
                    self.line(at(k), at(k + 1));
                }
                self.head(at(steps), outward, head);
            }
            Shape::Surface => {
                // A small square on the face, lifted off it a little.
                let (e1, e2) = perpendiculars(outward);
                let at = |a: f64, b: f64| c + outward * (0.1 * u) + (e1 * a + e2 * b) * (0.35 * u);
                let ring = [at(-1.0, -1.0), at(1.0, -1.0), at(1.0, 1.0), at(-1.0, 1.0)];
                for (k, a) in ring.iter().enumerate() {
                    self.line(*a, ring.get((k + 1) % 4).copied().unwrap_or(*a));
                }
            }
            Shape::Constraint
            | Shape::Support(_)
            | Shape::Force(_)
            | Shape::Moment(_)
            | Shape::Thermal
            | Shape::Temperature
            | Shape::Prescribed(_) => {}
        }
    }

    /// A rigid wall: a square through `point` normal to `normal`, `half` across, with hatching
    /// on its back.
    fn wall(&mut self, point: Vec3, normal: Vec3, half: f64, unit: f64, selected: bool) {
        let (u, v) = perpendiculars(normal);
        let c = |a: f64, b: f64| point + u * (a * half) + v * (b * half);
        let corners = [c(-1.0, -1.0), c(1.0, -1.0), c(1.0, 1.0), c(-1.0, 1.0)];
        for (k, a) in corners.iter().enumerate() {
            self.line(*a, corners.get((k + 1) % 4).copied().unwrap_or(*a));
        }
        // Slanted ticks behind the square, on a grid about two symbol units apart.
        let ticks = ((half / unit).floor() as usize).clamp(2, 20);
        let slant = (u + v).normalized().unwrap_or(u);
        let width = self.width;
        self.width = if selected { width } else { COUPLING_WIDTH };
        for i in 0..=ticks {
            for j in 0..=ticks {
                let at = c(-1.0 + 2.0 * i as f64 / ticks as f64, -1.0 + 2.0 * j as f64 / ticks as f64);
                self.line(at, at - normal * (0.8 * unit) + slant * (0.4 * unit));
            }
        }
        self.width = width;
    }
}

#[cfg(test)]
mod tests {
    use fieldcraft_model::{Coupling, CouplingKind, Dofs, ElementKind, ReferencePoint, Set};

    use super::*;
    use crate::{Camera, StandardView};

    const COLORS: GlyphColors = GlyphColors {
        constraint: Rgb(0, 0, 255),
        load: Rgb(255, 0, 0),
        interaction: Rgb(0, 255, 0),
        contact: Rgb(200, 200, 0),
        thermal: Rgb(255, 128, 0),
        label: Rgb(0, 0, 0),
        datum: Rgb(0, 128, 128),
        highlight: Rgb(255, 0, 255),
    };

    /// A row of `n` unit hex8 cells along x.
    fn bar(n: u32) -> Model {
        let mut m = Model::default();
        let id = |i: u32, j: u32, k: u32| 1 + i * 4 + j * 2 + k;
        for i in 0..=n {
            for j in 0..2 {
                for k in 0..2 {
                    m.set_node(id(i, j, k), [i as f64, j as f64, k as f64]).unwrap();
                }
            }
        }
        for i in 0..n {
            let nodes = vec![id(i, 0, 0), id(i + 1, 0, 0), id(i + 1, 1, 0), id(i, 1, 0), id(i, 0, 1), id(i + 1, 0, 1), id(i + 1, 1, 1), id(i, 1, 1)];
            m.set_element(i + 1, ElementKind::Hex8, nodes).unwrap();
        }
        m.add_set("root", Set { nodes: vec![1, 2, 3, 4], ..Set::default() }).unwrap();
        m.add_set("tip", Set { nodes: (0..4).map(|k| 4 * n + 1 + k).collect(), ..Set::default() }).unwrap();
        m
    }

    fn view(show: GlyphShow, unit: f64) -> GlyphView<'static> {
        GlyphView { show, colors: &COLORS, unit, highlight: None }
    }

    fn drawn_at(g: &Glyphs, m: &Model, v: &GlyphView) -> (Scene, Drawn) {
        let mut s = Scene::default();
        let d = g.draw(|l| m.node(l).map(|n| Vec3::from(n.xyz)), v, &mut s);
        (s, d)
    }

    fn drawn(g: &Glyphs, m: &Model, show: GlyphShow) -> Scene {
        drawn_at(g, m, &view(show, g.model_unit())).0
    }

    /// The midpoints of the triangles' bases drawn with their tip at `tip` (the lines that do not
    /// touch it and whose ends are each joined to it).
    fn bases(s: &Scene, tip: Vec3) -> Vec<Vec3> {
        let joined = |p: Vec3| s.lines.iter().any(|l| (l.a == tip && l.b.dist(p) < 1e-9) || (l.b == tip && l.a.dist(p) < 1e-9));
        s.lines.iter().filter(|l| l.a != tip && l.b != tip && joined(l.a) && joined(l.b)).map(|l| (l.a + l.b) * 0.5).collect()
    }

    fn near(v: &[Vec3], p: Vec3) -> usize {
        v.iter().filter(|q| q.dist(p) < 1e-9).count()
    }

    #[test]
    fn constraints_draw_a_wire_triangle_per_held_degree_of_freedom() {
        let u = 0.1;
        let mut m = bar(4);
        m.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
        let g = Glyphs::from_model(&m);
        assert!((g.model_unit() - m_diag(&m) * UNIT_FRACTION).abs() < 1e-12);
        // Solid nodes have no rotations: three crossed triangles (six lines each) on each of the
        // four nodes, wireframe only, in the constraint colour.
        let (s, _) = drawn_at(&g, &m, &view(GlyphShow::ALL, u));
        assert!(s.tris.is_empty() && s.on_top.is_empty());
        assert_eq!(s.lines.len(), 4 * 3 * 6);
        // A hidden definition draws nothing; hiding another leaves it.
        let (hidden, _) =
            drawn_at(&Glyphs::from_model_hiding(&m, &|d, n| d == Definition::BoundaryCondition && n == "Clamp"), &m, &view(GlyphShow::ALL, u));
        assert!(hidden.lines.is_empty());
        let (other, _) = drawn_at(&Glyphs::from_model_hiding(&m, &|d, _| d == Definition::Load), &m, &view(GlyphShow::ALL, u));
        assert_eq!(other.lines.len(), s.lines.len());
        assert!(s.lines.iter().all(|l| l.color == COLORS.constraint && l.width == LINE_WIDTH));
        // Tip on the node, each triangle along its axis (two crossed bases per axis).
        let b = bases(&s, Vec3::ZERO);
        for axis in [Vec3::X, Vec3::Y, Vec3::Z] {
            assert_eq!(near(&b, -axis * u), 2, "{axis:?}: {b:?}");
        }
        // Hidden with Loads & BCs.
        let s = drawn(&g, &m, GlyphShow { loads: false, interactions: true });
        assert!(s.lines.is_empty());
        // One held translation: one triangle a node, along x.
        let mut m = bar(4);
        let values = [Some(0.0), None, None, None, None, None];
        m.add_boundary_condition("Slide", BoundaryCondition::Displacement { target: Target::Set("root".into()), values }).unwrap();
        let (s, _) = drawn_at(&Glyphs::from_model(&m), &m, &view(GlyphShow::ALL, u));
        assert_eq!(s.lines.len(), 4 * 6);
        assert_eq!(near(&bases(&s, Vec3::ZERO), -Vec3::X * u), 2);
    }

    /// A prescribed motion is drawn as an arrow leaving each node along it (a bar, a shaft and a
    /// head), its held degrees of freedom still as supports (#354); a prescribed
    /// temperature is a thermometer, never a support (#355).
    #[test]
    fn prescribed_motion_and_temperature_have_their_own_symbols() {
        let u = 0.1;
        let mut m = bar(4);
        let values = [Some(0.0), None, Some(-2.5), None, None, None];
        m.add_boundary_condition("Press", BoundaryCondition::Displacement { target: Target::Set("root".into()), values }).unwrap();
        let (s, _) = drawn_at(&Glyphs::from_model(&m), &m, &view(GlyphShow::ALL, u));
        // Per node: one support triangle (x held at 0, six lines), and the motion's bar, shaft and
        // four-line head pointing down from the node.
        assert_eq!(s.lines.len(), 4 * (6 + 1 + 1 + 4), "{}", s.lines.len());
        let down = s.lines.iter().filter(|l| l.a.dist(Vec3::ZERO) < 1e-9 && l.b.dist(Vec3::new(0.0, 0.0, -2.2 * u)) < 1e-9).count();
        assert_eq!(down, 1, "a shaft leaves the node along the motion");
        assert!(s.lines.iter().all(|l| l.color == COLORS.constraint));

        let mut m = bar(4);
        m.add_boundary_condition("Heater", BoundaryCondition::Temperature { target: Target::Set("root".into()), value: 120.0 }).unwrap();
        let (s, _) = drawn_at(&Glyphs::from_model(&m), &m, &view(GlyphShow::ALL, u));
        assert!(s.lines.iter().all(|l| l.color == COLORS.thermal) && s.lines.len() >= 4 * (1 + ARC_SEGMENTS), "{}", s.lines.len());
        assert_eq!(s.points.len(), 4);
    }

    #[test]
    fn rotations_are_double_marks_and_symmetry_shows_its_normal_only() {
        let u = 0.5;
        let beam = |bc: BoundaryCondition| {
            let mut m = Model::default();
            m.set_node(1, [0.0; 3]).unwrap();
            m.set_node(2, [10.0, 0.0, 0.0]).unwrap();
            m.set_element(1, ElementKind::Beam2, vec![1, 2]).unwrap();
            m.add_set("end", Set { nodes: vec![1], ..Set::default() }).unwrap();
            m.add_boundary_condition("BC", bc).unwrap();
            let (s, _) = drawn_at(&Glyphs::from_model(&m), &m, &view(GlyphShow::ALL, u));
            s
        };
        let end = || Target::Set("end".into());
        // Encastre on a beam node: three translations and three rotations.
        let s = beam(BoundaryCondition::Fixed { target: end() });
        assert_eq!(s.lines.len(), 6 * 6);
        // The rotation's triangle sits behind the translation's: its tip at the first's base.
        assert_eq!(near(&bases(&s, -Vec3::X * u), -Vec3::X * (1.8 * u)), 2);
        // A symmetry plane normal to y holds u2, ur1 and ur3, and shows u2 alone.
        let s = beam(BoundaryCondition::Symmetry { target: end(), plane: fieldcraft_model::Symmetry::Y });
        assert_eq!(s.lines.len(), 6);
        assert_eq!(near(&bases(&s, Vec3::ZERO), -Vec3::Y * u), 2);
        // An inclined one shows its normal: one triangle along it.
        let plane = fieldcraft_model::Symmetry::from_normal([1.0, 1.0, 0.0]).unwrap();
        let s = beam(BoundaryCondition::Symmetry { target: end(), plane });
        assert_eq!(s.lines.len(), 6);
        let n = Vec3::new(1.0, 1.0, 0.0).normalized().unwrap();
        assert_eq!(near(&bases(&s, Vec3::ZERO), -n * u), 2);
        // A rotation alone keeps its place behind the node, apart from it.
        let s = beam(BoundaryCondition::Displacement { target: end(), values: [None, None, None, Some(0.0), None, None] });
        assert_eq!(s.lines.len(), 6);
        assert!(s.lines.iter().all(|l| l.a.x <= -u + 1e-9 && l.b.x <= -u + 1e-9));
    }

    /// Supports on a face come at it from outside the part.
    #[test]
    fn supports_on_faces_point_in_from_outside() {
        let u = 0.2;
        let mut m = bar(4);
        // The end face at x = 4, outward +x.
        m.add_set("end", Set { faces: vec![ElementFace { element: 4, face: 4 }], ..Set::default() }).unwrap();
        m.add_boundary_condition(
            "Hold",
            BoundaryCondition::Displacement { target: Target::Set("end".into()), values: [Some(0.0), None, None, None, None, None] },
        )
        .unwrap();
        let (s, _) = drawn_at(&Glyphs::from_model(&m), &m, &view(GlyphShow::ALL, u));
        assert_eq!(s.lines.len(), 4 * 6);
        assert!(s.lines.iter().all(|l| l.a.x >= 4.0 - 1e-9 && l.b.x >= 4.0 - 1e-9), "every triangle outside x = 4");
    }

    /// The same symbol measures the same number of pixels on screen at any zoom, orthographic or
    /// perspective.
    #[test]
    fn symbols_keep_their_size_on_screen_at_any_zoom() {
        let mut m = bar(4);
        m.add_set("one", Set { nodes: vec![1], ..Set::default() }).unwrap();
        m.add_boundary_condition(
            "Down",
            BoundaryCondition::Displacement { target: Target::Set("one".into()), values: [None, None, Some(0.0), None, None, None] },
        )
        .unwrap();
        let g = Glyphs::from_model(&m);
        let (w, h) = (800.0, 600.0);
        for fov in [0.0, 0.6] {
            for distance in [3.0, 30.0, 3000.0] {
                let mut cam = Camera { distance, fov, target: Vec3::ZERO, ..Default::default() };
                cam.set_view(StandardView::Front);
                cam.target = Vec3::ZERO;
                let unit = symbol_unit(cam.length_per_pixel(h).unwrap(), DEFAULT_SIZE).unwrap();
                let (s, _) = drawn_at(&g, &m, &view(GlyphShow::ALL, unit));
                let base = bases(&s, Vec3::ZERO);
                let at = |p: Vec3| cam.to_screen(p, w, h, 10.0 * distance).unwrap();
                let (a, b) = (at(Vec3::ZERO), at(base[0]));
                let px = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
                assert!((px - DEFAULT_SIZE).abs() < 0.05 * DEFAULT_SIZE, "fov {fov} distance {distance}: {px} px");
            }
        }
        // The size setting is held to its range; no camera gives no unit.
        assert_eq!(symbol_unit(1.0, 1e9), Some(SIZE_RANGE.1));
        assert_eq!(symbol_unit(1.0, f64::NAN), Some(DEFAULT_SIZE));
        assert_eq!(symbol_unit(f64::INFINITY, 10.0), None);
        // A view at twice the pixels per point has the same unit as at one (#285).
        let cam = crate::Camera::default();
        assert_eq!(unit_for(&cam, 1600.0, DEFAULT_SIZE, 2.0), unit_for(&cam, 800.0, DEFAULT_SIZE, 1.0));
        assert_eq!(unit_for(&cam, 800.0, DEFAULT_SIZE, f64::NAN), unit_for(&cam, 800.0, DEFAULT_SIZE, 1.0));
        assert_eq!(unit_for(&cam, 0.0, DEFAULT_SIZE, 1.0), None);
        assert_eq!(Camera::default().length_per_pixel(0.0), None);
        // Units round to a quarter of a doubling.
        let (k, q) = quantise_unit(1.1).unwrap();
        assert!(k == 1 && (q - 2f64.powf(0.25)).abs() < 1e-12, "{k} {q}");
        assert_eq!(quantise_unit(1.0).map(|x| x.0), Some(0));
        assert!(quantise_unit(-1.0).is_none() && quantise_unit(f64::NAN).is_none());
        assert!(quantise_unit(f64::MAX).is_some_and(|(_, q)| q.is_finite()));
    }

    /// Symbols spread over a set at a density that follows the screen, never more than the cap.
    #[test]
    fn sampling_follows_the_screen_density_and_is_bounded() {
        let mut m = Model::default();
        let n = 101u32;
        for i in 0..n {
            for j in 0..n {
                m.set_node(1 + i * n + j, [i as f64, j as f64, 0.0]).unwrap();
            }
        }
        m.add_set("all", Set { nodes: (1..=n * n).collect(), ..Set::default() }).unwrap();
        m.add_boundary_condition(
            "Clamp",
            BoundaryCondition::Displacement { target: Target::Set("all".into()), values: [None, None, Some(0.0), None, None, None] },
        )
        .unwrap();
        let g = Glyphs::from_model(&m);
        let count = |u: f64| {
            let (s, d) = drawn_at(&g, &m, &view(GlyphShow::ALL, u));
            (s.lines.len() / 6, d.notes)
        };
        // Zoomed in, the cap; zoomed out, one per grid cell of 2.4 units.
        let (close, notes) = count(0.01);
        assert!(close <= MAX_PER_DEFINITION && close > MAX_PER_DEFINITION * 9 / 10, "{close}");
        assert_eq!(notes, vec![format!("Clamp: {close} of {} nodes", n * n)]);
        for u in [2.0f64, 5.0, 20.0] {
            let cells = (100.0 / (2.4 * u)).floor() + 1.0;
            let (k, _) = count(u);
            assert!(k as f64 <= cells * cells && k > 0, "unit {u}: {k}");
        }
        assert!(count(20.0).0 < count(5.0).0 && count(5.0).0 < count(2.0).0);
    }

    #[test]
    fn colours_come_from_the_view_and_the_selected_definition_is_highlighted() {
        let mut m = bar(4);
        m.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: Target::Set("root".into()) }).unwrap();
        m.add_load("Push", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, -10.0], distribution: Default::default() }).unwrap();
        let g = Glyphs::from_model(&m);
        let v = GlyphView { highlight: Some((Definition::Load, "Push")), ..view(GlyphShow::ALL, 0.1) };
        let (s, _) = drawn_at(&g, &m, &v);
        let of = |c: Rgb| s.lines.iter().filter(|l| l.color == c).collect::<Vec<_>>();
        assert_eq!(of(COLORS.constraint).len(), 4 * 3 * 6);
        assert_eq!(of(COLORS.highlight).len(), 4 * 5, "the push's arrows: a shaft and four head lines");
        assert!(of(COLORS.highlight).iter().all(|l| l.width == HIGHLIGHT_WIDTH));
        assert!(of(COLORS.load).is_empty() && s.tris.is_empty());
        // Another definition of the same name and another kind is not highlighted.
        let v = GlyphView { highlight: Some((Definition::BoundaryCondition, "Push")), ..view(GlyphShow::ALL, 0.1) };
        assert!(drawn_at(&g, &m, &v).0.lines.iter().all(|l| l.color != COLORS.highlight));
    }

    #[test]
    fn forces_end_on_their_nodes_and_pressures_point_onto_their_faces() {
        let mut m = bar(4);
        m.add_load("Push", Load::Force { target: Target::Set("tip".into()), vector: [0.0, 0.0, -10.0], distribution: Default::default() }).unwrap();
        let g = Glyphs::from_model(&m);
        let s = drawn(&g, &m, GlyphShow::ALL);
        // A shaft per tip node coming down onto it, and an open head; set loads are depth-tested.
        let shafts: Vec<&SceneLine> = s.lines.iter().filter(|l| l.a.x == l.b.x && l.a.y == l.b.y).collect();
        assert_eq!(shafts.len(), 4);
        for l in shafts {
            assert!(l.a.z > l.b.z && l.a.x == 4.0 && l.color == COLORS.load, "{l:?}");
        }
        assert_eq!(s.lines.len(), 4 * 5);
        assert!(s.on_top.is_empty());
        // A pressure on the end face (x = 4): one arrow onto it, from outside (x > 4).
        let mut m = bar(4);
        m.add_set("end", Set { faces: vec![ElementFace { element: 4, face: 4 }], ..Set::default() }).unwrap();
        m.add_load("Press", Load::Pressure { set: "end".into(), magnitude: 2.0 }).unwrap();
        let s = drawn(&Glyphs::from_model(&m), &m, GlyphShow::ALL);
        let l = s.lines.first().unwrap();
        assert_eq!(s.lines.len(), 5);
        assert!(l.a.x > 4.0 && (l.b.x - 4.0).abs() < 1e-9, "{l:?}");
        // A suction points away.
        m.replace_load("Press", Load::Pressure { set: "end".into(), magnitude: -2.0 }).unwrap();
        let s = drawn(&Glyphs::from_model(&m), &m, GlyphShow::ALL);
        let l = s.lines.first().unwrap();
        assert!((l.a.x - 4.0).abs() < 1e-9 && l.b.x > l.a.x, "{l:?}");
    }

    /// A load at a reference point or a single node is drawn over the part; its support-sized
    /// arrow is four symbol units long.
    #[test]
    fn concentrated_loads_draw_on_top() {
        let mut m = bar(4);
        // Inside the bar: the part would hide it.
        m.add_reference_point("Mid", ReferencePoint { xyz: [2.0, 0.5, 0.5], label: 0 }).unwrap();
        m.add_load("F", Load::Force { target: Target::ReferencePoint("Mid".into()), vector: [0.0, 0.0, -1.0], distribution: Default::default() })
            .unwrap();
        m.add_set("one", Set { nodes: vec![5], ..Set::default() }).unwrap();
        m.add_load("G", Load::Force { target: Target::Set("one".into()), vector: [1.0, 0.0, 0.0], distribution: Default::default() }).unwrap();
        let g = Glyphs::from_model(&m);
        let (s, _) = drawn_at(&g, &m, &view(GlyphShow { loads: true, interactions: false }, 0.1));
        assert!(s.lines.is_empty());
        assert_eq!(s.on_top.len(), 2 * 5);
        let shaft = s.on_top.iter().find(|l| l.b == Vec3::new(2.0, 0.5, 0.5)).unwrap();
        assert!((shaft.a.z - (0.5 + 0.4)).abs() < 1e-9, "{shaft:?}");
        // The reference point's cross and dot too.
        let (s, _) = drawn_at(&g, &m, &view(GlyphShow { loads: false, interactions: true }, 0.1));
        assert_eq!(s.on_top.len(), 3 + 1);
        assert!(s.lines.is_empty());
    }

    #[test]
    fn moments_curl_right_handed_about_their_vector() {
        let mut m = bar(4);
        let c = Vec3::new(6.0, 0.5, 0.5);
        m.add_reference_point("RP", ReferencePoint { xyz: [c.x, c.y, c.z], label: 0 }).unwrap();
        m.add_load("Twist", Load::Moment { target: Target::ReferencePoint("RP".into()), vector: [3.0, 0.0, 0.0], distribution: Default::default() })
            .unwrap();
        let u = 0.2;
        let (s, _) = drawn_at(&Glyphs::from_model(&m), &m, &view(GlyphShow { loads: true, interactions: false }, u));
        let lines = &s.on_top;
        assert_eq!(lines.len(), 1 + ARC_SEGMENTS + 4);
        // The axis line through the point along x.
        let axis = lines.first().unwrap();
        assert!((axis.a + axis.b) * 0.5 == c && axis.a.y == c.y && axis.a.z == c.z && axis.a.x < axis.b.x, "{axis:?}");
        // The arc lies in the plane x = 6, a radius from the point.
        let radius = 1.3 * u;
        let arc = lines.get(1..1 + ARC_SEGMENTS).unwrap();
        for l in arc {
            for p in [l.a, l.b] {
                assert!((p.x - c.x).abs() < 1e-12 && (p.dist(c) - radius).abs() < 1e-9, "{p:?}");
            }
        }
        // The head at the arc's end points along the right-handed tangent, x × r.
        let end = arc.last().unwrap().b;
        let head = lines.get(1 + ARC_SEGMENTS..).unwrap();
        assert!(head.iter().all(|l| l.a == end));
        let back = head.iter().fold(Vec3::ZERO, |s, l| s + l.b) * 0.25;
        let pointing = (end - back).normalized().unwrap();
        let want = Vec3::X.cross(end - c).normalized().unwrap();
        assert!(pointing.dist(want) < 1e-9, "{pointing:?} vs {want:?}");
    }

    #[test]
    fn reference_points_couplings_contact_and_gravity() {
        let mut m = bar(4);
        m.add_reference_point("RP-1", ReferencePoint { xyz: [6.0, 0.5, 0.5], label: 0 }).unwrap();
        m.add_coupling("Tie", Coupling { reference_point: "RP-1".into(), set: "tip".into(), kind: CouplingKind::Kinematic, dofs: Dofs::all() })
            .unwrap();
        m.add_load("g", Load::Gravity { acceleration: [0.0, 0.0, -9810.0], set: None }).unwrap();
        m.add_set("end", Set { faces: vec![ElementFace { element: 4, face: 4 }], ..Set::default() }).unwrap();
        m.add_set("start", Set { faces: vec![ElementFace { element: 1, face: 6 }], ..Set::default() }).unwrap();
        m.add_contact("Touch", fieldcraft_model::Contact::new(ContactSurfaces::Pair { main: "end".into(), secondary: "start".into() })).unwrap();
        let g = Glyphs::from_model(&m);
        let s = drawn(&g, &m, GlyphShow::ALL);
        assert_eq!(s.labels, vec![SceneLabel { p: Vec3::new(6.0, 0.5, 0.5), text: "RP-1".into(), color: COLORS.label }]);
        assert_eq!(s.lines.iter().filter(|l| l.color == COLORS.interaction && l.width == COUPLING_WIDTH).count(), 4, "a line to each tip node");
        // A square on each contact face.
        assert_eq!(s.lines.iter().filter(|l| l.color == COLORS.contact).count(), 2 * 4);
        // Gravity: one arrow.
        assert_eq!(s.lines.iter().filter(|l| l.color == COLORS.load).count(), 5);
        // Interactions hidden: no point, cross, coupling or contact; the loads stay.
        let s = drawn(&g, &m, GlyphShow { loads: true, interactions: false });
        assert!(s.labels.is_empty() && s.on_top.is_empty() && s.lines.iter().all(|l| l.color == COLORS.load));
        // A results plot moves the symbols with the nodes; the reference point keeps its place
        // when the plot does not know it.
        let mut s = Scene::default();
        g.draw(|l| m.node(l).map(|n| Vec3::from(n.xyz) + Vec3::new(0.0, 0.0, 1.0)), &view(GlyphShow::ALL, g.model_unit()), &mut s);
        assert!(
            s.lines
                .iter()
                .filter(|l| l.width == COUPLING_WIDTH && l.color == COLORS.interaction)
                .all(|l| l.a == Vec3::new(6.0, 0.5, 0.5) && l.b.z >= 1.0)
        );
    }

    #[test]
    fn rigid_walls_back_their_set_and_initial_velocities_leave_it() {
        use fieldcraft_model::{ExplicitDynamic, Friction, InitialVelocity, RigidWall, Step};
        let mut m = bar(4);
        m.add_set("all", Set { nodes: m.nodes().keys().copied().collect(), ..Set::default() }).unwrap();
        let e = ExplicitDynamic {
            initial_velocities: vec![InitialVelocity { set: "all".into(), velocity: [-10.0, 0.0, 0.0] }],
            rigid_walls: vec![RigidWall { set: "root".into(), point: [0.0, 0.0, 0.0], normal: [1.0, 0.0, 0.0], friction: Friction::NONE }],
            ..ExplicitDynamic::new(1.0)
        };
        m.add_step(Step::new("Impact", Procedure::ExplicitDynamic(e)), None).unwrap();
        let g = Glyphs::from_model(&m);
        let s = drawn(&g, &m, GlyphShow { loads: false, interactions: true });
        // The wall: in the plane x = 0, centred behind the root face, hatched toward -x.
        let outline: Vec<&SceneLine> = s.lines.iter().filter(|l| l.width == LINE_WIDTH).collect();
        assert_eq!(outline.len(), 4);
        assert!(outline.iter().all(|l| l.a.x.abs() < 1e-9 && l.b.x.abs() < 1e-9));
        let mid = outline.iter().fold(Vec3::ZERO, |a, l| a + l.a) * 0.25;
        assert!(mid.dist(Vec3::new(0.0, 0.5, 0.5)) < 1e-9, "{mid:?}");
        assert!(s.lines.iter().filter(|l| l.width == COUPLING_WIDTH).all(|l| l.b.x < 0.0));
        // The velocity: one arrow behind the bar's +x end, pointing along -x toward it.
        let s = drawn(&g, &m, GlyphShow { loads: true, interactions: false });
        let shaft = s.lines.first().unwrap();
        assert_eq!(s.lines.len(), 5);
        assert!(shaft.b.x > 4.0 && shaft.b.x < shaft.a.x, "{shaft:?}");
        // Selecting the step highlights both.
        let v = GlyphView { highlight: Some((Definition::Step, "Impact")), ..view(GlyphShow::ALL, g.model_unit()) };
        assert!(drawn_at(&g, &m, &v).0.lines.iter().all(|l| l.color == COLORS.highlight));
    }

    #[test]
    fn colours_mute_toward_the_background() {
        let m = COLORS.muted(Rgb(255, 255, 255), 0.5);
        assert_eq!(m.constraint, Rgb(128, 128, 255));
        assert_eq!(m.highlight, COLORS.highlight);
        assert_eq!(COLORS.muted(Rgb(1, 2, 3), f64::NAN), COLORS);
    }

    /// A set of a quarter million nodes draws a capped sample and says so; hostile coordinates,
    /// vectors and units draw nothing rather than fail.
    #[test]
    fn huge_and_hostile_setups_are_sampled_and_bounded() {
        let mut m = Model::default();
        let n = 500u32;
        for i in 0..n {
            for j in 0..n {
                m.set_node(1 + i * n + j, [i as f64, j as f64, 0.0]).unwrap();
            }
        }
        m.add_set("all", Set { nodes: (1..=n * n).collect(), ..Set::default() }).unwrap();
        m.add_boundary_condition("Clamp", BoundaryCondition::Fixed { target: Target::Set("all".into()) }).unwrap();
        m.add_load("Push", Load::Force { target: Target::Set("all".into()), vector: [0.0, 0.0, -1.0], distribution: Default::default() }).unwrap();
        let g = Glyphs::from_model(&m);
        let (s, d) = drawn_at(&g, &m, &view(GlyphShow::ALL, 1e-3));
        assert!(d.notes.iter().any(|n| n.starts_with("Clamp: ") && n.ends_with(" of 250000 nodes")), "{:?}", d.notes);
        assert!(s.lines.len() <= MAX_PER_DEFINITION * (3 * 6 + 5));
        // Positions that are not finite are left out; a unit that is not a length is the model's.
        let mut s = Scene::default();
        g.draw(|_| Some(Vec3::new(f64::NAN, 0.0, 0.0)), &view(GlyphShow::ALL, f64::NAN), &mut s);
        assert!(s.tris.is_empty() && s.lines.is_empty());
        for unit in [f64::INFINITY, -1.0, 0.0, 1e300, 1e-300] {
            let _ = drawn_at(&g, &m, &view(GlyphShow::ALL, unit));
        }
        // An empty model, and one with a reference point only.
        assert!(Glyphs::from_model(&Model::default()).groups.is_empty());
        let mut m = Model::default();
        m.add_reference_point("RP", ReferencePoint { xyz: [1e9, -1e9, 0.0], label: 0 }).unwrap();
        let g = Glyphs::from_model(&m);
        let mut s = Scene::default();
        g.draw(|_| None, &view(GlyphShow::ALL, 1.0), &mut s);
        assert_eq!(s.labels.len(), 1);
    }

    #[test]
    fn grid_thinning_keeps_one_per_cell_and_respects_the_cap() {
        let items: Vec<(u32, Vec3)> = (0..100).map(|i| (i, Vec3::new(f64::from(i) * 0.1, 0.0, 0.0))).collect();
        assert_eq!(thin(items.iter().copied(), 1.0, 1000).len(), 10);
        assert_eq!(thin(items.iter().copied(), 1e-9, 7).len(), 7);
        assert_eq!(thin(items.iter().copied(), f64::NAN, 1000).len(), 10);
        assert!(thin::<u32>([], 1.0, 10).is_empty());
        assert!(thin([(1u32, Vec3::new(f64::INFINITY, 0.0, 0.0))], 1.0, 10).is_empty());
    }

    fn m_diag(m: &Model) -> f64 {
        let mut b = Aabb3::EMPTY;
        for n in m.nodes().values() {
            b.add(Vec3::from(n.xyz));
        }
        b.diagonal()
    }

    fn cube(lo: f64, hi: f64) -> Aabb3 {
        let mut b = Aabb3::EMPTY;
        b.add(Vec3::new(lo, lo, lo));
        b.add(Vec3::new(hi, hi, hi));
        b
    }

    #[test]
    fn datums_draw_at_the_model_with_their_names() {
        let mut m = Model::default();
        m.add_datum("Mid", fieldcraft_model::Datum::Plane { point: [5.0, 0.0, 0.0], normal: [1.0, 0.0, 0.0] }).unwrap();
        m.add_datum("Ax", fieldcraft_model::Datum::Axis { point: [0.0; 3], direction: [0.0, 0.0, 1.0] }).unwrap();
        m.add_datum("Pt", fieldcraft_model::Datum::Point { point: [1.0, 2.0, 3.0] }).unwrap();
        let bounds = cube(0.0, 10.0);
        let mut scene = Scene::default();
        let drawn = datums(&m, &bounds, &|_| true, None, &COLORS, &mut scene);
        // A rectangle (4 lines, filled by two translucent triangles), the axis line and the
        // point's three arms, each crossing at the point.
        assert_eq!(scene.lines.len(), 8);
        assert_eq!(scene.fills.len(), 2);
        assert!(scene.fills.iter().all(|f| f.alpha == DATUM_FILL_ALPHA && f.tri.normal == Vec3::X));
        let p = Vec3::new(1.0, 2.0, 3.0);
        assert!(scene.lines.iter().rev().take(3).all(|l| ((l.a + l.b) * 0.5).dist(p) < 1e-12));
        let mut names: Vec<&str> = scene.labels.iter().map(|l| l.text.as_str()).collect();
        names.sort_unstable();
        assert_eq!(names, ["Ax", "Mid", "Pt"]);
        // The axis spans the model along z with a 5 % margin: z from -0.5 to 10.5.
        let axis = scene.lines.iter().find(|l| l.width == SHAFT_WIDTH).unwrap();
        assert!((axis.a.z + 0.5).abs() < 1e-9 && (axis.b.z - 10.5).abs() < 1e-9, "{axis:?}");
        // Everything drawn stays within the model's box and its margin.
        let mut reach = bounds;
        reach.add(Vec3::new(-0.5, -0.5, -0.5));
        reach.add(Vec3::new(10.5, 10.5, 10.5));
        assert!(drawn.min.x >= reach.min.x - 1e-9 && drawn.max.z <= reach.max.z + 1e-9, "{drawn:?}");
        // A datum left out draws nothing; the highlighted one draws anyway, on top.
        let mut only = Scene::default();
        datums(&m, &bounds, &|n| n == "Ax", None, &COLORS, &mut only);
        assert_eq!((only.lines.len(), only.labels.len()), (1, 1));
        let mut picked = Scene::default();
        datums(&m, &bounds, &|_| false, Some("Mid"), &COLORS, &mut picked);
        assert_eq!((picked.lines.len(), picked.on_top.len(), picked.fills.len()), (0, 4, 2));
        assert!(picked.on_top.iter().all(|l| l.color == COLORS.highlight));
    }

    /// A datum plane at any offset or angle is drawn as the model's extent within it plus the
    /// margin: every corner lies on the plane and within the model's box grown by the margin.
    #[test]
    fn datum_planes_fit_the_model_at_any_offset_and_angle() {
        let bounds = cube(-20.0, 20.0);
        let size = 40.0;
        let normals = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [1.0, 2.0, 3.0], [0.3, -0.2, 0.9]];
        for nr in normals {
            let n = Vec3::from(nr).normalized().unwrap();
            for offset in [-19.0, -7.5, 0.0, 3.0, 15.0] {
                let p = n * offset;
                let mut m = Model::default();
                m.add_datum("D", fieldcraft_model::Datum::Plane { point: [p.x, p.y, p.z], normal: [n.x, n.y, n.z] }).unwrap();
                let mut scene = Scene::default();
                datums(&m, &bounds, &|_| true, None, &COLORS, &mut scene);
                // Every line lies on the plane.
                assert!(scene.lines.iter().all(|l| ((l.a - p).dot(n)).abs() < 1e-9 && ((l.b - p).dot(n)).abs() < 1e-9), "{nr:?} at {offset}");
                // The outline is the tightest rectangle in the plane's axes around the points
                // where the plane cuts the box's edges, plus the margin on each side.
                let (u, v) = perpendiculars(n);
                let mut cut = Vec::new();
                for a in 0..8usize {
                    for bit in [1usize, 2, 4] {
                        if a & bit == 0 {
                            let at = |k: usize| Vec3::new([-20.0, 20.0][k & 1], [-20.0, 20.0][(k >> 1) & 1], [-20.0, 20.0][(k >> 2) & 1]);
                            let (qa, qb) = (at(a), at(a | bit));
                            let (ha, hb) = ((qa - p).dot(n), (qb - p).dot(n));
                            if ha * hb <= 0.0 && ha != hb {
                                cut.push(qa + (qb - qa) * (ha / (ha - hb)));
                            }
                        }
                    }
                }
                let range = |pts: &mut dyn Iterator<Item = f64>| pts.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), x| (a.min(x), b.max(x)));
                let (cu, cv) = (range(&mut cut.iter().map(|q| (*q - p).dot(u))), range(&mut cut.iter().map(|q| (*q - p).dot(v))));
                let margin = DATUM_MARGIN * (cu.1 - cu.0).max(cv.1 - cv.0);
                let (ru, rv) = (
                    range(&mut scene.lines.iter().take(4).map(|l| (l.a - p).dot(u))),
                    range(&mut scene.lines.iter().take(4).map(|l| (l.a - p).dot(v))),
                );
                let close = |a: (f64, f64), b: (f64, f64)| (a.0 - (b.0 - margin)).abs() < 1e-9 && (a.1 - (b.1 + margin)).abs() < 1e-9;
                assert!(close(ru, cu) && close(rv, cv), "{nr:?} at {offset}: {ru:?} {rv:?} around {cu:?} {cv:?}");
                // Never larger than the model's section can be, plus the margins.
                assert!((ru.1 - ru.0).hypot(rv.1 - rv.0) <= (1.0 + 2.0 * DATUM_MARGIN) * 2f64.sqrt() * bounds.diagonal() + 1e-9);
                if nr.iter().filter(|x| **x != 0.0).count() == 1 {
                    // An axis-aligned plane: the cube's face plus 5 % on each side.
                    assert!((ru.1 - ru.0 - size * 1.1).abs() < 1e-9 && (rv.1 - rv.0 - size * 1.1).abs() < 1e-9);
                }
            }
        }
        // A plane missing the model is drawn as the model's shadow on it.
        let mut m = Model::default();
        m.add_datum("Far", fieldcraft_model::Datum::Plane { point: [0.0, 0.0, 100.0], normal: [0.0, 0.0, 1.0] }).unwrap();
        let mut scene = Scene::default();
        let drawn = datums(&m, &bounds, &|_| true, None, &COLORS, &mut scene);
        assert!((drawn.max.x - 22.0).abs() < 1e-9 && (drawn.min.y + 22.0).abs() < 1e-9, "{drawn:?}");
    }
}
