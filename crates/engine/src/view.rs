//! Headless pictures of a session: the model, its geometry and the results display drawn with
//! the CPU rasterizer of [`fieldcraft_render`], with no window or GPU (the MCP server's
//! screenshots without an app, docs/mcp.md).
//!
//! The picture follows the desktop viewport's conventions: geometry draws as faces with their
//! boundaries while the model has no element faces, then as boundaries only; elements in their
//! kind's colour with their edges; a results view ([`crate::results::ResultsView`]) draws the
//! deformed shape (automatic or chosen scale) with its contour bands, or a fluid step's slice,
//! coloured obstacles and streamlines over the undeformed part ([`crate::flow`]). Animation is not drawn: a
//! picture shows the full deformation (a mode at its positive shape). The analysis setup's
//! symbols ([`fieldcraft_render::glyphs`]) are drawn over the model as the viewport draws them,
//! muted over a results plot (and only those its step used), at their default size on screen
//! once the camera is known.

use fieldcraft_model::ElementKind;
use fieldcraft_model::units::{Quantity, Typed};
use fieldcraft_render::contour::{banded, bands_color};
use fieldcraft_render::glyphs::{self, GlyphColors, GlyphShow, GlyphView, Glyphs};
use fieldcraft_render::{Aabb3, Camera, DisplayMesh, Orientation, Rgb, Scene, SceneLine, ScenePoint, SceneTri, Vec3};
use serde_json::Value;

use crate::results::{Kind, deformation_unit, displacements, field_values, range, shown_scale};
use crate::{EngineError, Result, Session};

/// The colours of a picture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Palette {
    /// Background gradient, top to bottom.
    pub background: (Rgb, Rgb),
    pub hex: Rgb,
    pub tet: Rgb,
    pub wedge: Rgb,
    pub pyramid: Rgb,
    pub shell: Rgb,
    pub line_element: Rgb,
    pub mesh_edge: Rgb,
    /// Streamlines when the legend shows another field than velocity: stands out from the
    /// background.
    pub streamline: Rgb,
    /// Nodes, when they are shown.
    pub node: Rgb,
    pub geometry: Rgb,
    pub geometry_edge: Rgb,
    /// The bodies of a design with several, in turn.
    pub bodies: [Rgb; 8],
    /// Faces of nodes without a value in a contour plot.
    pub no_result: Rgb,
    /// The contour colour scale, low to high.
    pub contour: [Rgb; 5],
    /// Values above a clamped legend's maximum.
    pub contour_over: Rgb,
    /// The elements `mesh.quality.mark` marked.
    pub marked: Rgb,
    /// Faces merged for meshing (`geometry.faces.merge`).
    pub merged: Rgb,
    /// The analysis setup's symbols.
    pub glyphs: GlyphColors,
}

/// Values above a clamped legend's maximum (both themes): a deep purple apart from the scale's red.
pub const CONTOUR_OVER: Rgb = Rgb(0x6e, 0x2c, 0x80);

impl Palette {
    /// The desktop app's light theme (its `Tokens::light().palette()`, a test there checks).
    pub const LIGHT: Palette = Palette {
        background: (Rgb(240, 243, 247), Rgb(203, 210, 221)),
        hex: Rgb(0x78, 0x9c, 0xc0),
        tet: Rgb(0x80, 0xb4, 0x8a),
        wedge: Rgb(0xc0, 0xaa, 0x78),
        pyramid: Rgb(0x8c, 0x80, 0xb0),
        shell: Rgb(0xc0, 0x8c, 0x78),
        line_element: Rgb(0xc8, 0x8c, 0x1e),
        mesh_edge: Rgb(0x1e, 0x23, 0x2a),
        streamline: Rgb(0x28, 0x2e, 0x38),
        node: Rgb(0xb4, 0x3c, 0x8c),
        geometry: Rgb(176, 186, 198),
        geometry_edge: Rgb(40, 44, 52),
        bodies: [
            Rgb(0x8f, 0xb4, 0xd9),
            Rgb(0xe0, 0xa4, 0x6c),
            Rgb(0x8c, 0xc4, 0x8a),
            Rgb(0xd6, 0x8c, 0x96),
            Rgb(0xb0, 0x9c, 0xd8),
            Rgb(0xd8, 0xc8, 0x70),
            Rgb(0x7c, 0xc4, 0xc0),
            Rgb(0xc0, 0xa0, 0x84),
        ],
        no_result: Rgb(0x9a, 0xa0, 0xaa),
        contour: [Rgb(0x2b, 0x4f, 0xd8), Rgb(0x2f, 0xb8, 0xe0), Rgb(0x3c, 0xc8, 0x64), Rgb(0xf0, 0xd2, 0x3c), Rgb(0xe0, 0x40, 0x2e)],
        contour_over: CONTOUR_OVER,
        marked: Rgb(214, 54, 54),
        merged: Rgb(222, 150, 30),
        glyphs: GlyphColors {
            constraint: Rgb(0x00, 0x8a, 0xa8),
            load: Rgb(0xd0, 0x24, 0x24),
            interaction: Rgb(0x9a, 0x30, 0xc0),
            thermal: Rgb(0xe8, 0x4c, 0x00),
            label: Rgb(0x26, 0x2a, 0x32),
            datum: Rgb(0x1e, 0x8c, 0x6e),
            contact: Rgb(0x6a, 0x8a, 0x00),
            highlight: Rgb(0xd6, 0x00, 0x7f),
        },
    };

    /// The colour elements of this kind are drawn in.
    pub fn element(&self, kind: ElementKind) -> Rgb {
        match kind {
            ElementKind::Hex8 | ElementKind::Hex20 => self.hex,
            ElementKind::Tet4 | ElementKind::Tet10 => self.tet,
            ElementKind::Wedge6 => self.wedge,
            ElementKind::Pyramid5 => self.pyramid,
            ElementKind::Tri3 | ElementKind::Tri6 | ElementKind::Quad4 | ElementKind::Quad8 => self.shell,
            ElementKind::Rod2 | ElementKind::Rod3 | ElementKind::Beam2 => self.line_element,
        }
    }
}

/// How the model is drawn.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Style {
    /// Lit faces with the element edges on them.
    #[default]
    ShadedEdges,
    /// Lit faces only.
    Shaded,
    /// Element edges only.
    Wireframe,
}

impl Style {
    pub const ALL: [Style; 3] = [Style::ShadedEdges, Style::Shaded, Style::Wireframe];
    pub fn label(self) -> &'static str {
        match self {
            Style::ShadedEdges => "Shaded with edges",
            Style::Shaded => "Shaded",
            Style::Wireframe => "Wireframe",
        }
    }
}

/// What is shown (the browser's eye toggles).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Show {
    /// The open design and imported surfaces.
    pub geometry: bool,
    /// Elements.
    pub mesh: bool,
    pub nodes: bool,
    /// Boundary conditions and loads.
    pub loads: bool,
    /// Reference points, couplings and rigid walls.
    pub interactions: bool,
}

impl Show {
    /// Geometry, mesh and the analysis setup, without nodes.
    pub const DEFAULT: Show = Show { geometry: true, mesh: true, nodes: false, loads: true, interactions: true };
}

/// Mesh edge width, geometry edge width, line element width and node size, in points.
const EDGE_WIDTH: f64 = 1.0;
const GEOMETRY_EDGE_WIDTH: f64 = 1.5;
const LINE_ELEMENT_WIDTH: f64 = 3.0;
const NODE_SIZE: f64 = 5.0;

/// Coordinates beyond this magnitude in a camera command are refused.
const MAX_TARGET: f64 = 1e12; // unit-factors: a cap on hostile input
/// A closed surface's volume is well away from 0 when above this fraction of its area times
/// its size.
const ENCLOSES: f64 = 1e-6; // unit-factors: relative to the area times the size
/// The smallest sine of the half angle a perspective fit divides by.
const MIN_SINE: f64 = 1e-6; // unit-factors: the sine of an angle

/// The scene for `display` drawn in `style` with `p`'s colours: the one builder of the window's
/// viewport and the headless pictures, so both draw the same model.
///
/// Geometry (the open design, imported surfaces) draws as faces with their boundaries while
/// the model has no element faces; once it is meshed, as its boundaries drawn along the mesh
/// (the element edges on geometry edges), so the mesh is not hidden behind coincident faces and
/// every line lies on the element faces.
pub fn build(display: &DisplayMesh, style: Style, show: Show, p: &Palette) -> Scene {
    let mut s = Scene { background: Some(p.background), ..Default::default() };
    // Geometry faces draw where no element faces cover them.
    let mesh_faces = show.mesh && !display.tris.is_empty();
    if style != Style::Wireframe {
        if show.mesh {
            s.tris = display.tris.iter().map(|d| SceneTri { p: d.p, normal: d.normal, color: p.element(d.kind) }).collect();
        }
        if !mesh_faces && show.geometry {
            s.tris.extend(display.surface_tris.iter().map(|d| SceneTri {
                p: d.p,
                normal: d.normal,
                color: display.body_color(d, p.geometry, &p.bodies),
            }));
        }
    }
    if style != Style::Shaded && show.mesh {
        // Without faces behind them the edges take the element colour, so they read on the
        // background.
        let color = if style == Style::Wireframe { p.hex } else { p.mesh_edge };
        s.lines = display.edges.iter().map(|e| SceneLine { a: e.p[0], b: e.p[1], color, width: EDGE_WIDTH }).collect();
    }
    // Face boundaries in every style: they are the geometry's outline. Over the mesh they are
    // the element edges on geometry edges, so they lie exactly on the element faces; the
    // geometry's own tessellation lies off them by its chord error.
    if show.geometry {
        let color = p.geometry_edge;
        if mesh_faces {
            s.lines.extend(display.edges.iter().filter(|e| e.feature).map(|e| SceneLine { a: e.p[0], b: e.p[1], color, width: GEOMETRY_EDGE_WIDTH }));
        } else {
            s.lines.extend(display.surface_edges.iter().map(|[a, b]| SceneLine { a: *a, b: *b, color, width: GEOMETRY_EDGE_WIDTH }));
        }
    }
    if show.mesh {
        s.lines.extend(display.lines.iter().map(|l| SceneLine { a: l.a, b: l.b, color: p.element(l.kind), width: LINE_ELEMENT_WIDTH }));
    }
    if show.nodes {
        s.points = display.nodes.iter().map(|n| ScenePoint { p: n.p, color: p.node, size: NODE_SIZE }).collect();
    }
    s
}

/// The chordal tolerance bodies are drawn with, in the geometry's millimetres (`b` too): a
/// thousandth of their size, within 1 µm–1 mm.
pub fn display_tolerance(b: &Aabb3) -> f64 {
    const RELATIVE: f64 = 1e-3; // unit-factors: relative to the body's size
    static RANGE: std::sync::OnceLock<(f64, f64)> = std::sync::OnceLock::new();
    let (lo, hi) = *RANGE.get_or_init(|| {
        let mm = |t: Typed| t.express(Quantity::Length, crate::units::geometry()).unwrap_or(f64::NAN);
        (mm(Typed::new(1.0, "µm")), mm(Typed::new(1.0, "mm")))
    });
    let d = b.diagonal();
    // `max` and `min`, not `clamp`: no panic whatever the bounds.
    if d.is_finite() { (d * RELATIVE).max(lo).min(hi) } else { 0.1 }
}

/// The findings of a `geometry.check`, to draw over the design it checked.
#[derive(Clone, Debug, Default)]
pub struct GeometryMarks {
    /// The design checked: the marks are stale once another is open.
    pub design: std::sync::Weak<fieldcraft_geometry::Design>,
    pub marks: Vec<GeometryMark>,
}

/// One finding to draw: its faces, its edge and the point to look at.
#[derive(Clone, Debug, PartialEq)]
pub struct GeometryMark {
    /// The faces by body name and face index.
    pub faces: Vec<(String, u32)>,
    /// A short edge's polyline (empty for the other kinds).
    pub path: Vec<[f64; 3]>,
    pub at: [f64; 3],
    /// A small or sliver face, by body and persistent name: what `geometry.face.absorb`
    /// `from_check` absorbs.
    pub small: Option<(String, String)>,
}

/// The findings of the last `geometry.check`, while the design it checked is still open.
pub fn geometry_marks(s: &Session) -> Option<&GeometryMarks> {
    let marks = s.geometry_marks.as_deref()?;
    let design = s.design.as_ref()?;
    std::sync::Weak::ptr_eq(&marks.design, &std::sync::Arc::downgrade(design)).then_some(marks)
}

/// Marker line width, and the marker dot's diameter, in output pixels.
const MARK_WIDTH: f64 = 4.0;
const MARK_DOT: f64 = 12.0;

/// Draw `marks` over `scene` in `color`: the faces they name recoloured when `faces` (the
/// scene's triangles start with `display`'s surface triangles in order, as they do while no
/// element faces are drawn; `bodies` names the drawn surfaces in order), short edges as thick
/// lines and a dot at each finding's point, both over everything so a feature too small to see
/// at the zoom, or inside the mesh, still shows where it is.
pub fn draw_geometry_marks(scene: &mut Scene, display: &DisplayMesh, bodies: &[&str], marks: &GeometryMarks, faces: bool, color: Rgb) {
    if faces {
        let marked: std::collections::BTreeSet<(u32, u32)> = marks
            .marks
            .iter()
            .flat_map(|m| &m.faces)
            .filter_map(|(body, f)| bodies.iter().position(|b| b == body).and_then(|b| u32::try_from(b).ok()).map(|b| (b, *f)))
            .collect();
        for (t, d) in scene.tris.iter_mut().zip(&display.surface_tris) {
            if marked.contains(&(d.body, d.face)) {
                t.color = color;
            }
        }
    }
    // The marks are on the geometry, in millimetres.
    let point = |p: [f64; 3]| Vec3::from(crate::units::geometry_point(p));
    for m in &marks.marks {
        for w in m.path.windows(2) {
            if let [a, b] = w {
                scene.on_top.push(SceneLine { a: point(*a), b: point(*b), color, width: MARK_WIDTH });
            }
        }
        let at = point(m.at);
        scene.on_top.push(SceneLine { a: at, b: at, color, width: MARK_DOT });
    }
}

/// Draw the merged faces (`geometry.faces.merge`) of the drawn bodies (`bodies`, in the order of
/// their surfaces) in `color`, so each shows as one face in one colour. The scene's triangles
/// start with `display`'s surface triangles in order, as they do while no element faces are drawn.
pub fn draw_merged_faces(scene: &mut Scene, display: &DisplayMesh, s: &Session, bodies: &[&str], color: Rgb) {
    let Some(design) = s.design.as_deref() else { return };
    let merged: std::collections::BTreeSet<(u32, u32)> = bodies
        .iter()
        .enumerate()
        .filter_map(|(i, name)| Some((u32::try_from(i).ok()?, design.body(name)?)))
        .flat_map(|(i, b)| b.merged_faces().into_iter().map(move |(f, _)| (i, f)))
        .collect();
    if merged.is_empty() {
        return;
    }
    for (t, d) in scene.tris.iter_mut().zip(&display.surface_tris) {
        if merged.contains(&(d.body, d.face)) {
            t.color = color;
        }
    }
}

/// What a picture shows besides the image.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Shown {
    /// The solved step drawn, when the results display is on.
    pub step: Option<String>,
    /// The deformation scale factor drawn (0: undeformed).
    pub scale: f64,
    /// The contour field drawn, with its smallest and largest value and their nodes (a fluid
    /// step's: their cells).
    pub contour: Option<(String, (u32, f64), (u32, f64))>,
    /// Bodies that could not be tessellated for display.
    pub failed: Vec<String>,
    /// Definitions whose symbols are a sample (`Clamp: 400 of 12000 nodes`).
    pub sampled: Vec<String>,
}

/// A section through the mesh: the elements whose centroid lies on the camera's side of a plane
/// normal to an axis are not drawn, so the element faces inside show (an element cut: whole elements, not
/// sliced ones). Geometry surfaces are not drawn over a sectioned mesh.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Section {
    /// 0, 1, 2: the plane is normal to X, Y, Z.
    pub axis: usize,
    /// Where the plane crosses the axis; `None` at the middle of the mesh.
    pub at: Option<f64>,
}

impl Section {
    /// `x`, `y`, `z` (the middle of the mesh) or `x=12.5`.
    pub fn parse(s: &str) -> Option<Section> {
        let (axis, at) = match s.split_once('=') {
            Some((a, v)) => (a, Some(v.trim().parse::<f64>().ok().filter(|v| v.is_finite())?)),
            None => (s, None),
        };
        let axis = ["x", "y", "z"].iter().position(|n| n.eq_ignore_ascii_case(axis.trim()))?;
        Some(Section { axis, at })
    }

    /// The elements on the side of the plane `toward` points to (the direction toward the eye).
    fn beyond(&self, m: &fieldcraft_model::Model, toward: Vec3) -> std::collections::BTreeSet<u32> {
        let sign = match [toward.x, toward.y, toward.z].get(self.axis) {
            Some(v) if *v < 0.0 => -1.0,
            _ => 1.0,
        };
        let coord = |label: &u32| m.node(*label).and_then(|n| n.xyz.get(self.axis).copied());
        let centroid = |nodes: &[u32]| {
            let v: Vec<f64> = nodes.iter().filter_map(coord).collect();
            (!v.is_empty()).then(|| v.iter().sum::<f64>() / v.len() as f64)
        };
        let at = self.at.unwrap_or_else(|| {
            let (lo, hi) = m
                .nodes()
                .values()
                .filter_map(|n| n.xyz.get(self.axis).copied())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)));
            if lo <= hi { (lo + hi) * 0.5 } else { 0.0 }
        });
        m.elements().iter().filter(|(_, e)| centroid(&e.nodes).is_some_and(|c| (c - at) * sign > 0.0)).map(|(l, _)| *l).collect()
    }
}

/// The model's display data with the open design's bodies and the imported surfaces.
pub fn display(s: &Session) -> (DisplayMesh, Vec<String>) {
    display_cut(s, None, Vec3::Z)
}

/// [`display`] with a section through the mesh, cut away on the side `toward` points to.
pub fn display_cut(s: &Session, section: Option<Section>, toward: Vec3) -> (DisplayMesh, Vec<String>) {
    let (d, failed, _) = display_parts(s, section, toward);
    (d, failed)
}

/// [`display_cut`], and the names of the bodies drawn in the order of their surfaces
/// ([`fieldcraft_render::SurfaceTri::body`]).
fn display_parts(s: &Session, section: Option<Section>, toward: Vec3) -> (DisplayMesh, Vec<String>, Vec<String>) {
    let mut bodies = Vec::new();
    let mut hidden = s.hidden_elements();
    // The section's position is given in the display units; the mesh is in metres.
    let section = section.map(|c| Section { at: c.at.map(|x| length_si(s, x)), ..c });
    let cut = section.filter(|_| !s.model.elements().is_empty());
    if let Some(c) = cut {
        hidden.extend(c.beyond(&s.model, toward));
    }
    let mut d = DisplayMesh::from_model_hiding(&s.model, &hidden);
    let mut failed = Vec::new();
    if cut.is_some() {
        return (d, failed, bodies);
    }
    for body in s.design.iter().flat_map(|d| &d.bodies).filter(|b| s.shows(fieldcraft_io::ItemKind::Body, &b.name)) {
        let mesh = body.measure().and_then(|m| body.tessellate(display_tolerance(&Aabb3 { min: m.bbox_min.into(), max: m.bbox_max.into() })));
        match mesh {
            Ok(mut m) => {
                // A merged face's B-rep faces were tessellated apart: their copies of the points
                // they share differ by up to the weld distance, and would draw as edges.
                if !body.merged_faces().is_empty() {
                    weld_points(&mut m.positions, m.weld_distance);
                }
                // The geometry is in millimetres, the model in metres.
                let positions: Vec<[f64; 3]> = m.positions.iter().map(|p| crate::units::geometry_point(*p)).collect();
                d.add_surface(&positions, &m.triangles, &m.tri_face);
                bodies.push(body.name.clone());
            }
            Err(e) => failed.push(format!("{}: {e}", body.name)),
        }
    }
    for surface in s.surfaces.iter().filter(|x| s.shows(fieldcraft_io::ItemKind::Surface, &x.name)) {
        let positions: Vec<[f64; 3]> = surface.mesh.positions.iter().map(|p| crate::units::geometry_point(*p)).collect();
        d.add_surface(&positions, &surface.mesh.triangles, &surface.mesh.tri_face);
    }
    (d, failed, bodies)
}

/// Move every point within `distance` of an earlier one (in x order) onto it, so copies of one
/// point are identical.
fn weld_points(positions: &mut [[f64; 3]], distance: f64) {
    let mut order: Vec<usize> = (0..positions.len()).collect();
    order.sort_by(|&a, &b| positions.get(a).map_or(0.0, |p| p[0]).total_cmp(&positions.get(b).map_or(0.0, |p| p[0])));
    let close = |p: [f64; 3], q: [f64; 3]| (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>() <= distance * distance;
    for (k, &i) in order.iter().enumerate() {
        let Some(p) = positions.get(i).copied() else { continue };
        // The earliest point within reach, among those not far behind in x.
        let mut to = p;
        for &j in order.get(..k).unwrap_or_default().iter().rev() {
            let Some(q) = positions.get(j).copied() else { continue };
            if p[0] - q[0] > distance {
                break;
            }
            if close(p, q) {
                to = q;
            }
        }
        if let Some(slot) = positions.get_mut(i) {
            *slot = to;
        }
    }
}

/// The model drawn as the viewport draws it by default (shaded with element edges).
fn model_scene(d: &DisplayMesh, p: &Palette) -> Scene {
    build(d, Style::ShadedEdges, Show::DEFAULT, p)
}

/// The elements `mesh.quality.mark` marked (its set), drawn in the marker colour.
pub fn marked(s: &Session) -> std::collections::BTreeSet<u32> {
    s.model.set(crate::cmd::mesh::MARKED_SET).map(|set| set.elements.iter().copied().collect()).unwrap_or_default()
}

/// Colour the faces of the `marked` elements in `scene`, whose triangles start with the element
/// faces of `display` in order (as [`model_scene`] and the app's viewport build them).
pub fn mark(scene: &mut Scene, display: &DisplayMesh, marked: &std::collections::BTreeSet<u32>, color: Rgb) {
    if marked.is_empty() {
        return;
    }
    for (t, d) in scene.tris.iter_mut().zip(&display.tris) {
        if marked.contains(&d.element) {
            t.color = color;
        }
    }
}

/// The model part of a picture: the model or the results view's plot, the shape the symbols sit
/// on (the deformed one of a plot), the box around it and what it shows.
struct Picture {
    scene: Scene,
    drawn: DisplayMesh,
    bounds: Aabb3,
    shown: Shown,
    glyphs: Glyphs,
    colors: GlyphColors,
    /// The box around the datums drawn: inside the clipping radius, never framed (#295).
    datums: Aabb3,
    /// The bodies drawn, in the order of their surfaces.
    bodies: Vec<String>,
}

impl Picture {
    /// The setup's symbols at symbol unit `unit` ([`glyphs::symbol_unit`]), in a scene of their
    /// own; notes the definitions drawn with a sample. Returns the box around them.
    fn symbols(&mut self, unit: f64, into: &mut Scene) -> Aabb3 {
        let drawn = &self.drawn;
        let view = GlyphView { show: GlyphShow::ALL, colors: &self.colors, unit, highlight: None };
        let d = self.glyphs.draw(|label| drawn.node_position(label), &view, into);
        self.shown.sampled = d.notes;
        d.bounds
    }
}

fn picture(s: &Session, p: &Palette, section: Option<Section>, toward: Vec3) -> Picture {
    let (base, failed, bodies) = display_parts(s, section, toward);
    let mut shown = Shown { failed, ..Shown::default() };
    let plot = s.results_view.as_ref().and_then(|v| s.latest_result(&v.step).map(|(r, _)| (v, r)));
    // Over a results plot, the symbols of what its step used.
    let glyphs = plot.map_or_else(|| s.glyphs(), |(v, _)| s.step_glyphs(&v.step));
    let Some((view, result)) = plot else {
        let mut scene = model_scene(&base, p);
        mark(&mut scene, &base, &marked(s), p.marked);
        // `bodies` names the drawn surfaces in order.
        let drawn: Vec<&str> = bodies.iter().map(String::as_str).collect();
        if base.tris.is_empty() {
            draw_merged_faces(&mut scene, &base, s, &drawn, p.merged);
        }
        if let Some(marks) = geometry_marks(s) {
            draw_geometry_marks(&mut scene, &base, &drawn, marks, base.tris.is_empty(), p.marked);
        }
        let bounds = base.bounds;
        let shows = |n: &str| s.shows(fieldcraft_io::ItemKind::Datum, n);
        let datums = fieldcraft_render::glyphs::datums(&s.model, &base.bounds, &shows, None, &p.glyphs, &mut scene);
        return Picture { scene, drawn: base, bounds, shown, glyphs, colors: p.glyphs, datums, bodies };
    };
    shown.step = Some(view.step.clone());
    let colors = p.glyphs.over_results(p.background.1);
    if matches!(result.kind(), Kind::Fluid(_)) {
        // The part as it is, with the flow drawn through and around it.
        let mut scene = model_scene(&base, p);
        let mut bounds = base.bounds;
        let names: Vec<&str> = bodies.iter().map(String::as_str).collect();
        if let Ok(fp) = crate::flow::plot_with_surfaces(view, result, &s.traced, s, &base, &names) {
            if let Some(c) = &fp.contour {
                let cell = |(k, v): (usize, f64)| (u32::try_from(k).unwrap_or(u32::MAX), v);
                shown.contour = Some((c.field.id().to_string(), cell(c.min), cell(c.max)));
            }
            crate::flow::draw(&fp, &mut scene, &crate::flow::FlowColors { stops: &p.contour, over: p.contour_over, plain: p.streamline });
            bounds.union(fp.bounds());
        }
        return Picture { scene, drawn: base, bounds, shown, glyphs, colors, datums: Aabb3::EMPTY, bodies };
    }
    let pick = view.pick();
    let nodes = result.nodes();
    let u = displacements(result, pick).unwrap_or(&[]);
    let scale = if !view.deformed { 0.0 } else { shown_scale(view, result, pick, base.bounds.diagonal(), s.units()) };
    shown.scale = scale;
    let factor = scale * deformation_unit(result, s.units());
    let drawn = base.displaced(|label| {
        let i = nodes.binary_search(&label).ok()?;
        u.get(i).map(|d| Vec3::from(*d) * factor)
    });
    let mut scene = model_scene(&drawn, p);
    let contour = view.contour.and_then(|field| {
        let values = field_values(&s.model, result, field, pick).ok()?;
        let (min, max) = range(nodes, &values)?;
        Some((field, min, max, values))
    });
    if let Some((field, min, max, values)) = contour {
        shown.contour = Some((field.id().to_string(), min, max));
        let bands = view.bands(min.1, max.1);
        let at = |label: u32| nodes.binary_search(&label).ok().and_then(|i| values.get(i).copied()).unwrap_or(f64::NAN);
        let color = |band: Option<u32>| band.map_or(p.no_result, |k| bands_color(&p.contour, p.contour_over, &bands, k));
        scene.tris.clear();
        for tri in &drawn.tris {
            banded(tri.p, tri.nodes.map(at), tri.normal, &bands, &color, &mut scene.tris);
        }
    }
    let mut bounds = base.bounds;
    bounds.union(drawn.bounds);
    Picture { scene, drawn, bounds, shown, glyphs, colors, datums: Aabb3::EMPTY, bodies }
}

/// The scene of the session (the model, or the results view's plot) with the setup's symbols
/// at the model's own symbol unit, the box around it, and what it shows.
pub fn scene(s: &Session, p: &Palette) -> (Scene, Aabb3, Shown) {
    let mut pic = picture(s, p, None, Vec3::Z);
    let mut symbols = Scene::default();
    let b = pic.symbols(pic.glyphs.model_unit(), &mut symbols);
    let mut scene = pic.scene;
    scene.append(symbols);
    pic.bounds.union(b);
    (scene, pic.bounds, pic.shown)
}

/// A camera orientation from JSON, `{"yaw": …, "pitch": …, "distance"?: …}`: degrees, and the
/// distance from the target in the display units' length (converted where the session draws it) (the control channel's `ui.camera`, the MCP
/// `screenshot` tool's `camera`). Hostile input is refused or clamped
/// ([`Orientation::checked`]).
pub fn orientation(p: &Value) -> Result<Orientation> {
    let bad = |m: &str| EngineError::Params(format!("camera: {m}"));
    let o = p.as_object().ok_or_else(|| bad("an object {yaw, pitch, distance?} in degrees"))?;
    if let Some(k) = o.keys().find(|k| !matches!(k.as_str(), "yaw" | "pitch" | "distance")) {
        return Err(bad(&format!("unknown field `{}` (yaw, pitch, distance)", k.chars().take(40).collect::<String>())));
    }
    let number = |k: &str| o.get(k).filter(|v| !v.is_null()).map(|v| v.as_f64().ok_or_else(|| bad(&format!("`{k}` must be a number"))));
    let yaw = number("yaw").ok_or_else(|| bad("`yaw` (degrees) is required"))??;
    let pitch = number("pitch").ok_or_else(|| bad("`pitch` (degrees) is required"))??;
    let distance = number("distance").transpose()?;
    Orientation::checked(yaw, pitch, distance).ok_or_else(|| bad("`yaw` and `pitch` must be finite and `distance` finite and above 0"))
}

/// What a camera looks at: a point, or the centre of mass or of the bounding box of the bodies
/// named (every body drawn when `None`).
#[derive(Clone, Debug, PartialEq)]
pub enum Target {
    Point(Vec3),
    Mass(Option<Vec<String>>),
    Bbox(Option<Vec<String>>),
}

impl Target {
    /// The default: the centre of mass of the model drawn, so a sweep of yaws orbits it.
    pub const MODEL: Target = Target::Mass(None);

    fn bodies(&self) -> Option<&[String]> {
        match self {
            Target::Point(_) => None,
            Target::Mass(b) | Target::Bbox(b) => b.as_deref(),
        }
    }

    /// The JSON form [`camera`] reads.
    pub fn json(&self) -> Value {
        let center = |what: &str, bodies: &Option<Vec<String>>| match bodies {
            Some(b) => serde_json::json!({"center_of": what, "bodies": b}),
            None => serde_json::json!({"center_of": what}),
        };
        match self {
            Target::Point(p) => serde_json::json!([p.x, p.y, p.z]),
            Target::Mass(b) => center("mass", b),
            Target::Bbox(b) => center("bbox", b),
        }
    }
}

/// A camera aimed as agents give it ([`camera`]): an orientation, what it looks at, whether its
/// distance fits the target's bodies to the frame, and a perspective field of view.
#[derive(Clone, Debug, PartialEq)]
pub struct Aim {
    pub at: Orientation,
    /// `None` keeps the camera's target (the app's `ui.camera`); a picture of its own aims at
    /// [`Target::MODEL`].
    pub target: Option<Target>,
    /// Without a distance: the target's bodies fill the frame from this direction. Otherwise the
    /// distance keeps the whole model in the frame from every direction, so a sweep keeps its
    /// size.
    pub fit: bool,
    /// Vertical field of view in degrees; 0 is orthographic.
    pub fov: Option<f64>,
}

impl Aim {
    /// Largest field of view, in degrees.
    pub const MAX_FOV: f64 = 150.0;
    /// Most bodies a target names.
    pub const MAX_BODIES: usize = 1000;

    /// Looking from `at`, with the target left as it is.
    pub fn new(at: Orientation) -> Self {
        Aim { at, target: None, fit: false, fov: None }
    }

    /// The JSON form [`camera`] reads.
    pub fn json(&self) -> Value {
        let mut v = serde_json::json!({"yaw": self.at.yaw, "pitch": self.at.pitch});
        if let Some(d) = self.at.distance {
            v["distance"] = d.into();
        }
        if let Some(t) = &self.target {
            v["target"] = t.json();
        }
        if self.fit {
            v["fit"] = true.into();
        }
        if let Some(f) = self.fov {
            v["fov"] = f.into();
        }
        v
    }
}

/// A camera from JSON: [`orientation`]'s `{yaw, pitch, distance?}` and `target`: `[x, y, z]` or
/// `{"center_of": "mass" | "bbox", "bodies"?: [names]}`; `fit` (the target's bodies fill the
/// frame; not with a distance) and `fov` (vertical, degrees, 0 orthographic). Hostile input is
/// refused or clamped.
pub fn camera(p: &Value) -> Result<Aim> {
    let bad = |m: &str| EngineError::Params(format!("camera: {m}"));
    let o = p.as_object().ok_or_else(|| bad("an object {yaw, pitch, distance?, target?, fit?, fov?}"))?;
    if let Some(k) = o.keys().find(|k| !matches!(k.as_str(), "yaw" | "pitch" | "distance" | "target" | "fit" | "fov")) {
        return Err(bad(&format!("unknown field `{}` (yaw, pitch, distance, target, fit, fov)", k.chars().take(40).collect::<String>())));
    }
    let given = |k: &str| o.get(k).filter(|v| !v.is_null());
    let mut core = serde_json::Map::new();
    for k in ["yaw", "pitch", "distance"] {
        if let Some(v) = given(k) {
            core.insert(k.into(), v.clone());
        }
    }
    let at = orientation(&Value::Object(core))?;
    let target = given("target").map(target).transpose()?;
    let fit = match given("fit") {
        None => false,
        Some(v) => v.as_bool().ok_or_else(|| bad("`fit` must be true or false"))?,
    };
    if fit && at.distance.is_some() {
        return Err(bad("give `distance` or `fit`, not both"));
    }
    let fov = match given("fov") {
        None => None,
        Some(v) => Some(
            v.as_f64()
                .filter(|f| f.is_finite() && (0.0..=Aim::MAX_FOV).contains(f))
                .ok_or_else(|| bad(&format!("`fov` must be a number of degrees from 0 (orthographic) to {}", Aim::MAX_FOV)))?,
        ),
    };
    Ok(Aim { at, target, fit, fov })
}

fn target(v: &Value) -> Result<Target> {
    let bad = |m: &str| EngineError::Params(format!("camera: target {m}"));
    if let Some(a) = v.as_array() {
        let c: Vec<f64> = a.iter().filter_map(Value::as_f64).filter(|x| x.is_finite() && x.abs() <= MAX_TARGET).collect();
        return match c[..] {
            [x, y, z] if a.len() == 3 => Ok(Target::Point(Vec3::new(x, y, z))),
            _ => Err(bad("must be [x, y, z] (finite numbers)")),
        };
    }
    let o = v.as_object().ok_or_else(|| bad("must be [x, y, z] or {center_of: mass | bbox, bodies?}"))?;
    if let Some(k) = o.keys().find(|k| !matches!(k.as_str(), "center_of" | "bodies")) {
        return Err(bad(&format!("has an unknown field `{}` (center_of, bodies)", k.chars().take(40).collect::<String>())));
    }
    let bodies = match o.get("bodies").filter(|v| !v.is_null()) {
        None => None,
        Some(b) => {
            let a = b
                .as_array()
                .filter(|a| !a.is_empty() && a.len() <= Aim::MAX_BODIES)
                .ok_or_else(|| bad(&format!("bodies must be a list of 1 to {} names", Aim::MAX_BODIES)))?;
            Some(a.iter().map(|n| n.as_str().map(str::to_string).ok_or_else(|| bad("bodies must be names"))).collect::<Result<Vec<_>>>()?)
        }
    };
    match o.get("center_of").and_then(Value::as_str).unwrap_or("mass") {
        "mass" => Ok(Target::Mass(bodies)),
        "bbox" => Ok(Target::Bbox(bodies)),
        other => Err(bad(&format!("center_of `{}`: mass or bbox", other.chars().take(40).collect::<String>()))),
    }
}

/// The centre of mass of the solid a closed triangle surface bounds (the divergence theorem), or
/// the centre of area of an open one; `None` without triangles.
fn centre_of_mass<'a>(tris: impl Iterator<Item = &'a [Vec3; 3]>) -> Option<Vec3> {
    let (mut volume, mut moment, mut area, mut areal) = (0.0, Vec3::ZERO, 0.0, Vec3::ZERO);
    for [a, b, c] in tris {
        let v = a.dot(b.cross(*c)) / 6.0;
        volume += v;
        moment = moment + (*a + *b + *c) * (v / 4.0);
        let s = (*b - *a).cross(*c - *a).len() * 0.5;
        area += s;
        areal = areal + (*a + *b + *c) * (s / 3.0);
    }
    if !(area > 0.0) {
        return None;
    }
    // A closed surface's volume is well away from 0 compared with its area times its size;
    // an open one (a shell) has none worth the name.
    let size = area.sqrt();
    let c = if volume.abs() > ENCLOSES * area * size { moment * (1.0 / volume) } else { areal * (1.0 / area) };
    c.is_finite().then_some(c)
}

impl Picture {
    /// The point `target` names, and the box around what it frames: the bodies it names, or
    /// everything drawn.
    fn aim_point(&self, target: &Target) -> Result<(Vec3, Aabb3)> {
        let surfaces: Vec<&[Vec3; 3]> = match target.bodies() {
            Some(names) => {
                let mut picked = Vec::new();
                for n in names {
                    let i = self.bodies.iter().position(|b| b == n).ok_or_else(|| {
                        EngineError::Params(format!(
                            "camera: no body `{}` is drawn (drawn: {})",
                            n.chars().take(80).collect::<String>(),
                            if self.bodies.is_empty() { "none".into() } else { self.bodies.join(", ") }
                        ))
                    })?;
                    picked.push(u32::try_from(i).unwrap_or(u32::MAX));
                }
                self.drawn.surface_tris.iter().filter(|t| picked.contains(&t.body)).map(|t| &t.p).collect()
            }
            // What is drawn: the elements' faces over the geometry's.
            None if !self.drawn.tris.is_empty() => self.drawn.tris.iter().map(|t| &t.p).collect(),
            None => self.drawn.surface_tris.iter().map(|t| &t.p).collect(),
        };
        let mut b = Aabb3::EMPTY;
        for p in surfaces.iter().flat_map(|t| t.iter()) {
            b.add(*p);
        }
        let framed = if target.bodies().is_some() { b } else { self.bounds };
        let whole = if b.is_empty() { self.bounds } else { b };
        let point = match target {
            Target::Point(p) => *p,
            Target::Bbox(_) => whole.center(),
            Target::Mass(_) => centre_of_mass(surfaces.into_iter()).unwrap_or_else(|| whole.center()),
        };
        Ok((point, framed))
    }
}

/// `cam` aimed as `aim` says with the picture `pic` and a view `aspect` wide per unit of height:
/// turned, its field of view set, then (with a target) moved to look at it from `aim`'s distance,
/// from the distance that fits the target's bodies (`fit`), or from the distance at which the
/// sphere about the target through the farthest corner of everything drawn fits the view.
fn aimed(pic: &Picture, aim: &Aim, mut cam: Camera, aspect: f64) -> Result<Camera> {
    cam.orient(&Orientation { distance: None, ..aim.at });
    if let Some(f) = aim.fov {
        cam.fov = f.to_radians();
    }
    let Some(target) = &aim.target else {
        cam.orient(&aim.at);
        return Ok(cam);
    };
    let (point, framed) = pic.aim_point(target)?;
    if !point.is_finite() {
        return Ok(cam);
    }
    cam.target = point;
    if aim.at.distance.is_some() {
        cam.orient(&aim.at);
    } else if aim.fit {
        cam.fit_about(&framed, aspect);
    } else if !pic.bounds.is_empty() {
        let aspect = if aspect.is_finite() && aspect > fieldcraft_render::ASPECT_RANGE.0 { aspect } else { 1.0 };
        let radius = pic.bounds.corners().iter().map(|c| c.dist(point)).fold(fieldcraft_render::MIN_SIZE, f64::max)
            / (1.0 - 2.0 * fieldcraft_render::FIT_MARGIN);
        // The half height the sphere needs, then the distance that shows it.
        let half = radius / aspect.min(1.0);
        let d = if cam.fov > 0.0 {
            let half_angle = ((cam.fov * 0.5).tan() * aspect.min(1.0)).atan();
            radius / half_angle.sin().max(MIN_SINE)
        } else {
            half / Camera { distance: 1.0, ..cam }.half_height()
        };
        if d.is_finite() {
            cam.distance = d.clamp(Orientation::DISTANCE.0, Orientation::DISTANCE.1);
        }
    }
    Ok(cam)
}

/// The app's camera `cam` aimed as `aim` says at the session `s` drawn in a view `aspect` wide
/// per unit of height (the control channel's `ui.camera`); `cam` as it is but turned when `aim`
/// has no target.
pub fn aim(s: &Session, aim: &Aim, cam: Camera, aspect: f64) -> Result<Camera> {
    let pic = picture(s, &Palette::LIGHT, None, cam.back());
    aimed(&pic, &aim_si(s, aim), cam, aspect)
}

/// A picture ready to draw: the scene (with its background and clipping radius), the setup's
/// symbols as a layer over it, and the camera.
pub struct Framed {
    pub scene: Scene,
    pub symbols: Scene,
    pub cam: Camera,
    pub shown: Shown,
}

/// The session seen from `at` in a `w` × `h` pixel picture: framing everything drawn but the
/// datums (they are drawn at the model's extent, and framing them would only add margin), or from
/// `at`'s distance when it has one. The symbols are sized for their default size on screen
/// ([`glyphs::DEFAULT_SIZE`] pixels) for that camera. The CPU rasterizer ([`render_png`]) and
/// the GPU (`fieldcraft-cli snapshot`) draw the same [`Framed`]. A `section` draws the mesh cut
/// ([`Section`]).
pub fn frame(s: &Session, at: &Orientation, w: usize, h: usize, p: &Palette, section: Option<Section>) -> Framed {
    let at = &orientation_si(s, at);
    let mut cam = Camera::default();
    cam.orient(at);
    let mut pic = picture(s, p, section, cam.back());
    let aspect = w as f64 / h.max(1) as f64;
    let unit = |cam: &Camera| glyphs::unit_for(cam, h as f64, glyphs::DEFAULT_SIZE, 1.0).unwrap_or(f64::NAN);
    // Framed on the model, then again with its symbols at that camera's size, which reach a
    // little beyond it.
    cam.orient(&Orientation { distance: None, ..*at });
    cam.fit_in(&pic.bounds, aspect, fieldcraft_render::FULL_VIEW);
    let mut symbols = Scene::default();
    let mut framed = pic.bounds;
    framed.union(pic.symbols(unit(&cam), &mut symbols));
    cam.fit_in(&framed, aspect, fieldcraft_render::FULL_VIEW);
    cam.orient(at);
    symbols = Scene::default();
    let mut bounds = pic.bounds;
    bounds.union(pic.symbols(unit(&cam), &mut symbols));
    bounds.union(pic.datums);
    pic.scene.radius =
        if bounds.is_empty() { 1.0 } else { (cam.target.dist(bounds.center()) + bounds.diagonal() * 0.5).max(fieldcraft_render::MIN_SIZE) };
    Framed { scene: pic.scene, symbols, cam, shown: pic.shown }
}

/// The session seen through `aim` ([`camera`]) in a `w` × `h` pixel picture, aimed at
/// [`Target::MODEL`] unless it names a target; the symbols at their default size on screen and
/// a `section` as [`frame`] draws them. An error names a body the target asks for that is not
/// drawn.
pub fn frame_at(s: &Session, aim: &Aim, w: usize, h: usize, p: &Palette, section: Option<Section>) -> Result<Framed> {
    let aim = aim_si(s, aim);
    let aim = Aim { target: Some(aim.target.clone().unwrap_or(Target::MODEL)), ..aim };
    let mut cam = Camera::default();
    cam.orient(&aim.at);
    let mut pic = picture(s, p, section, cam.back());
    let cam = aimed(&pic, &aim, Camera { fov: 0.0, ..cam }, w as f64 / h.max(1) as f64)?;
    let mut symbols = Scene::default();
    let mut bounds = pic.bounds;
    bounds.union(pic.symbols(glyphs::unit_for(&cam, h as f64, glyphs::DEFAULT_SIZE, 1.0).unwrap_or(f64::NAN), &mut symbols));
    bounds.union(pic.datums);
    pic.scene.radius =
        if bounds.is_empty() { 1.0 } else { (cam.target.dist(bounds.center()) + bounds.diagonal() * 0.5).max(fieldcraft_render::MIN_SIZE) };
    Ok(Framed { scene: pic.scene, symbols, cam, shown: pic.shown })
}

/// A length given in the session's display units, in metres (NaN when it cannot be one).
fn length_si(s: &Session, x: f64) -> f64 {
    s.units().to_si(fieldcraft_model::units::Quantity::Length, x).unwrap_or(f64::NAN)
}

/// `at` with its distance in metres (within [`Orientation::DISTANCE`]): it is given in the
/// display units.
fn orientation_si(s: &Session, at: &Orientation) -> Orientation {
    let (near, far) = Orientation::DISTANCE;
    Orientation { distance: at.distance.map(|d| length_si(s, d)).filter(|d| d.is_finite()).map(|d| d.clamp(near, far)), ..*at }
}

/// `aim` with its lengths (distance, target point) in metres: they are given in the display
/// units.
fn aim_si(s: &Session, aim: &Aim) -> Aim {
    let target = match &aim.target {
        Some(Target::Point(p)) => Some(Target::Point(Vec3::new(length_si(s, p.x), length_si(s, p.y), length_si(s, p.z)))),
        t => t.clone(),
    };
    Aim { at: orientation_si(s, &aim.at), target, ..aim.clone() }
}

/// The camera `cam` as an [`Aim`] in the session's display units: what `view.save` keeps of a
/// host's camera (the app's viewport, an MCP picture), looking at its own target from its own
/// distance.
pub fn aim_of(cam: &Camera, s: &Session) -> Aim {
    let length = |m: f64| s.units().express(fieldcraft_model::units::Quantity::Length, m).unwrap_or(f64::NAN);
    let at = Orientation {
        yaw: cam.yaw.to_degrees(),
        pitch: cam.pitch.to_degrees(),
        distance: Some(length(cam.distance)).filter(|d| d.is_finite() && *d > 0.0),
    };
    let target = Vec3::new(length(cam.target.x), length(cam.target.y), length(cam.target.z));
    Aim {
        at,
        target: target.is_finite().then_some(Target::Point(target)),
        fit: false,
        fov: Some(cam.fov.to_degrees()).filter(|f| f.is_finite() && *f > 0.0),
    }
}

/// The camera `cam` as agents read it: its fields, with the distance and the target in the
/// session's display units (the camera itself works in metres).
pub fn camera_shown(cam: &Camera, s: &Session) -> Value {
    let mut v = serde_json::json!(cam);
    let length = |m: f64| s.units().express(fieldcraft_model::units::Quantity::Length, m).unwrap_or(f64::NAN);
    if let Some(o) = v.as_object_mut() {
        o.insert("distance".into(), serde_json::json!(length(cam.distance)));
        o.insert("target".into(), serde_json::json!({"x": length(cam.target.x), "y": length(cam.target.y), "z": length(cam.target.z)}));
    }
    v
}

/// A PNG of the session seen from `at`, `w` × `h` pixels ([`frame`]), drawn on the CPU.
pub fn render_png(s: &Session, at: &Orientation, w: usize, h: usize, p: &Palette) -> Result<(Vec<u8>, Shown)> {
    let f = frame(s, at, w, h, p, None);
    let png = fieldcraft_render::render_layers_png(&f.scene, &[&f.symbols], &f.cam, w, h).map_err(|e| EngineError::Params(format!("render: {e}")))?;
    Ok((png, f.shown))
}

#[cfg(test)]
mod tests {
    use fieldcraft_render::StandardView;
    use serde_json::json;

    use super::*;

    #[test]
    fn renders_the_model_and_its_results() {
        let iso = Orientation::of(StandardView::Iso);
        let mut s = Session::default();
        let empty = render_png(&s, &iso, 32, 24, &Palette::LIGHT).unwrap();
        assert!(empty.0.starts_with(&[0x89, b'P', b'N', b'G']));
        s.run_script(&crate::sample::cantilever()).unwrap();
        let (model, shown) = render_png(&s, &iso, 64, 48, &Palette::LIGHT).unwrap();
        assert_eq!(shown, Shown::default());
        // The clamp's triangles and the push's arrows are drawn in their colours, as wireframes.
        let (scene, ..) = super::scene(&s, &Palette::LIGHT);
        assert!(scene.lines.iter().any(|l| l.color == Palette::LIGHT.glyphs.constraint));
        assert!(scene.lines.iter().any(|l| l.color == Palette::LIGHT.glyphs.load));
        assert!(!scene.tris.iter().any(|t| t.color == Palette::LIGHT.glyphs.constraint));
        // And in the picture: pixels of the clamp's colour.
        let (png, _) = render_png(&s, &iso, 320, 240, &Palette::LIGHT).unwrap();
        assert!(count_color(&png, Palette::LIGHT.glyphs.constraint) > 20);
        s.execute("solve.static", &json!({})).unwrap();
        let (plot, shown) = render_png(&s, &iso, 64, 48, &Palette::LIGHT).unwrap();
        assert_ne!(model, plot);
        // A static step first shows its contours on the undeformed shape.
        assert!(shown.scale == 0.0 && shown.step.is_some(), "{shown:?}");
        let (field, min, max) = shown.contour.unwrap();
        assert_eq!(field, "von_mises");
        assert!(max.1 > min.1);
        assert!(render_png(&s, &Orientation::of(StandardView::Top), 0, 48, &Palette::LIGHT).is_err());
        // Over the plot the symbols are muted.
        let (scene, ..) = super::scene(&s, &Palette::LIGHT);
        let muted = Palette::LIGHT.glyphs.over_results(Palette::LIGHT.background.1);
        assert!(scene.lines.iter().any(|l| l.color == muted.constraint));
        assert!(!scene.lines.iter().any(|l| l.color == Palette::LIGHT.glyphs.constraint));
    }

    /// A session of hex8 boxes `[min, max]`, one element each.
    fn boxes(b: &[([f64; 3], [f64; 3])]) -> Session {
        let mut s = Session::default();
        for (e, (lo, hi)) in b.iter().enumerate() {
            let e = e as u32;
            let corner = |i: usize| {
                let pick = |k: usize, bit: bool| if bit { hi[k] } else { lo[k] };
                let (x, y, z) = ([false, true, true, false][i % 4], [false, false, true, true][i % 4], i >= 4);
                [pick(0, x), pick(1, y), pick(2, z)]
            };
            for i in 0..8 {
                s.execute("mesh.node.add", &json!({"label": e * 8 + i as u32 + 1, "xyz": corner(i)})).unwrap();
            }
            s.execute("mesh.element.add", &json!({"label": e + 1, "kind": "hex8", "nodes": (1..=8).map(|i| e * 8 + i).collect::<Vec<_>>()})).unwrap();
        }
        s
    }

    /// RGBA pixels of a PNG.
    fn pixels(png: &[u8]) -> Vec<u8> {
        let mut r = png::Decoder::new(std::io::Cursor::new(png)).read_info().unwrap();
        let mut buf = vec![0; r.output_buffer_size().unwrap()];
        let info = r.next_frame(&mut buf).unwrap();
        buf.truncate(info.buffer_size());
        buf
    }

    /// The centroid of the pixels of `png` that differ from the empty picture `empty`, and how
    /// many there are.
    fn drawn_centroid(png: &[u8], empty: &[u8], w: usize) -> (f64, f64, usize) {
        let (a, b) = (pixels(png), pixels(empty));
        let (mut sx, mut sy, mut n) = (0.0, 0.0, 0usize);
        for (i, (p, q)) in a.chunks(4).zip(b.chunks(4)).enumerate() {
            if p.iter().zip(q).take(3).any(|(x, y)| x.abs_diff(*y) > 8) {
                sx += (i % w) as f64 + 0.5;
                sy += (i / w) as f64 + 0.5;
                n += 1;
            }
        }
        (sx / n.max(1) as f64, sy / n.max(1) as f64, n)
    }

    fn shot(s: &Session, aim: &Aim, w: usize, h: usize) -> Vec<u8> {
        let f = frame_at(s, aim, w, h, &Palette::LIGHT, None).unwrap();
        fieldcraft_render::render_layers_png(&f.scene, &[&f.symbols], &f.cam, w, h).unwrap()
    }

    /// A yaw sweep of a box far from the origin orbits it: the box stays at the picture's
    /// centre and keeps its size, from below and above, orthographic and in perspective (#467).
    #[test]
    fn a_yaw_sweep_orbits_the_model() {
        let s = boxes(&[([200.0, 40.0, 30.0], [260.0, 70.0, 50.0])]);
        let (w, h) = (240, 160);
        let empty = shot(&Session::default(), &camera(&json!({"yaw": 0, "pitch": 0})).unwrap(), w, h);
        for extra in [json!({}), json!({"fov": 40}), json!({"fit": true})] {
            for pitch in [25.0, -30.0] {
                let mut sizes = Vec::new();
                for yaw in (0..360).step_by(30) {
                    let mut p = json!({"yaw": yaw, "pitch": pitch});
                    p.as_object_mut().unwrap().extend(extra.as_object().unwrap().clone());
                    let (x, y, n) = drawn_centroid(&shot(&s, &camera(&p).unwrap(), w, h), &empty, w);
                    assert!(n > 500, "{p}: {n} pixels");
                    // In perspective the near end of the box is drawn larger than the far end, which
                    // moves its silhouette's centroid off the centre of mass by a few percent of
                    // the width; the centre of mass itself projects to the centre (checked below).
                    let tolerance = if extra.get("fov").is_some() { 0.05 * w as f64 } else { 2.0 };
                    let off = (x - w as f64 / 2.0).abs().max((y - h as f64 / 2.0).abs());
                    assert!(off < tolerance, "{p}: centroid ({x:.1}, {y:.1})");
                    sizes.push(n);
                    let cam = frame_at(&s, &camera(&p).unwrap(), w, h, &Palette::LIGHT, None).unwrap().cam;
                    let (cx, cy, _) = cam.to_screen(Vec3::new(0.230, 0.055, 0.040), w as f64, h as f64, 100.0).unwrap();
                    assert!((cx - w as f64 / 2.0).abs() < 1e-6 && (cy - h as f64 / 2.0).abs() < 1e-6, "{p}");
                }
                // Without `fit` the distance is the same from every side; with it the box fills
                // the frame, so it is larger.
                if extra.get("fit").is_none() {
                    let (lo, hi) = (sizes.iter().min().unwrap(), sizes.iter().max().unwrap());
                    assert!(*hi < lo * 3, "{extra} {pitch}: {sizes:?}");
                }
            }
        }
    }

    /// A mode shape drawn at a fixed scale moves a node by the shape shown (mass-normalised in
    /// the display units, as the probe reads it) times the scale, in the display length: as
    /// large as before SI storage, and as large again in other display units' own terms.
    #[test]
    fn a_mode_drawn_at_a_fixed_scale_moves_by_the_shape_shown() {
        let mut s = Session::default();
        s.run_script(&crate::sample::cantilever()).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        s.execute("results.modes", &json!({"mode": 1})).unwrap();
        s.execute("results.deformed", &json!({"show": true, "scale": 5})).unwrap();
        let tip = *s.model.nodes().iter().find(|(_, n)| n.xyz == [0.1, 0.0, 0.0]).unwrap().0;
        for (system, length) in [("mm-t-s-N-MPa", 1e-3), ("m-kg-s-N-Pa", 1.0)] {
            s.execute("model.units", &json!({"system": system})).unwrap();
            let shape = s.execute("solve.probe", &json!({"step": "Modes", "nodes": [tip], "mode": 1})).unwrap()["nodes"][0]["shape"].clone();
            let shape: [f64; 3] = serde_json::from_value(shape).unwrap();
            let pic = picture(&s, &Palette::LIGHT, None, Vec3::Z);
            let moved = pic.drawn.node_position(tip).unwrap() - Vec3::new(0.1, 0.0, 0.0);
            let want = Vec3::from(shape) * (5.0 * length);
            assert!(moved.dist(want) <= 1e-12 * want.dist(Vec3::ZERO), "{system}: {moved:?} vs {want:?}");
        }
    }

    /// The target is the centre of mass by default, the bounding box's centre or a point when
    /// asked; bodies name what it frames, and a body not drawn is an error.
    #[test]
    fn camera_targets() {
        // An L: a long box and a short one on top of its end.
        let s = boxes(&[([0.0, 0.0, 0.0], [90.0, 10.0, 10.0]), ([0.0, 0.0, 10.0], [10.0, 10.0, 40.0])]);
        let pic = picture(&s, &Palette::LIGHT, None, Vec3::Z);
        let (mass, _) = pic.aim_point(&Target::MODEL).unwrap();
        // Volumes 9000 and 3000: x = (45·9000 + 5·3000)/12000 = 35, z = (5·9000 + 25·3000)/12000 = 10.
        // (The scene is in metres.)
        assert!(mass.dist(Vec3::new(0.035, 0.005, 0.010)) < 1e-12, "{mass:?}");
        let (bbox, _) = pic.aim_point(&Target::Bbox(None)).unwrap();
        assert!(bbox.dist(Vec3::new(0.045, 0.005, 0.020)) < 1e-12, "{bbox:?}");
        let (point, _) = pic.aim_point(&Target::Point(Vec3::new(1.0, 2.0, 3.0))).unwrap();
        assert_eq!(point, Vec3::new(1.0, 2.0, 3.0));
        // An open surface (one face of a box) has its centre of area.
        let face = [[Vec3::ZERO, Vec3::X, Vec3::new(1.0, 1.0, 0.0)], [Vec3::ZERO, Vec3::new(1.0, 1.0, 0.0), Vec3::Y]];
        assert!(centre_of_mass(face.iter()).unwrap().dist(Vec3::new(0.5, 0.5, 0.0)) < 1e-12);
        assert_eq!(centre_of_mass([].iter()), None);
        let err = frame_at(
            &s,
            &camera(&json!({"yaw": 0, "pitch": 0, "target": {"center_of": "mass", "bodies": ["Gear"]}})).unwrap(),
            64,
            48,
            &Palette::LIGHT,
            None,
        );
        assert!(err.is_err_and(|e| e.to_string().contains("no body `Gear`")));
        // The app's camera without a target keeps its own and only turns.
        let cam = Camera { target: Vec3::new(7.0, 8.0, 9.0), distance: 321.0, ..Camera::default() };
        let turned = aim(&s, &camera(&json!({"yaw": 90, "pitch": 10})).unwrap(), cam, 1.5).unwrap();
        assert_eq!((turned.target, turned.distance), (cam.target, cam.distance));
        // A point and a distance are given in the display units (mm here) and drawn in metres.
        let given = aim(&s, &camera(&json!({"yaw": 0, "pitch": 0, "target": [35, 5, 10], "distance": 200})).unwrap(), cam, 1.5).unwrap();
        assert!(given.target.dist(mass) < 1e-12 && (given.distance - 0.2).abs() < 1e-15, "{given:?}");
        let aimed = aim(&s, &camera(&json!({"yaw": 90, "pitch": 10, "target": {"center_of": "mass"}, "fit": true})).unwrap(), cam, 1.5).unwrap();
        assert!(aimed.target.dist(mass) < 1e-9 && aimed.distance != cam.distance);
    }

    #[test]
    fn camera_json_is_checked() {
        let a = camera(&json!({"yaw": 10, "pitch": 20, "target": [1, 2, 3], "fov": 35})).unwrap();
        assert_eq!((a.target, a.fov, a.fit), (Some(Target::Point(Vec3::new(1.0, 2.0, 3.0))), Some(35.0), false));
        let a = camera(&json!({"yaw": 0, "pitch": 0, "target": {"center_of": "bbox", "bodies": ["A"]}, "fit": true})).unwrap();
        assert_eq!(camera(&a.json()).unwrap(), a);
        assert_eq!(a.target, Some(Target::Bbox(Some(vec!["A".into()]))));
        assert_eq!(camera(&json!({"yaw": 0, "pitch": 0, "target": {}})).unwrap().target, Some(Target::MODEL));
        for bad in [
            json!({"yaw": 0, "pitch": 0, "target": [1, 2]}),
            json!({"yaw": 0, "pitch": 0, "target": [1, 2, "x"]}),
            json!({"yaw": 0, "pitch": 0, "target": [1, 2, 1e300]}),
            json!({"yaw": 0, "pitch": 0, "target": {"center_of": "gravity"}}),
            json!({"yaw": 0, "pitch": 0, "target": {"center_of": "mass", "bodies": []}}),
            json!({"yaw": 0, "pitch": 0, "target": {"center_of": "mass", "bodies": [3]}}),
            json!({"yaw": 0, "pitch": 0, "target": {"middle": 1}}),
            json!({"yaw": 0, "pitch": 0, "target": "mass"}),
            json!({"yaw": 0, "pitch": 0, "fit": true, "distance": 10}),
            json!({"yaw": 0, "pitch": 0, "fit": 1}),
            json!({"yaw": 0, "pitch": 0, "fov": -1}),
            json!({"yaw": 0, "pitch": 0, "fov": 179}),
            json!({"yaw": 0, "pitch": 0, "fov": "wide"}),
            json!({"yaw": 0, "pitch": 0, "zoom": 2}),
            json!({"pitch": 0}),
        ] {
            assert!(camera(&bad).is_err(), "{bad}");
        }
    }

    /// A fluid step draws the part with its faces coloured by the flow beside them, a coloured
    /// slice through the domain and streamlines (#193: the part vanished and nothing was drawn).
    #[test]
    fn a_fluid_step_draws_its_flow_around_the_part() {
        let mut s = crate::flow::tests::plate_in_a_tank();
        let p = Palette::LIGHT;
        let (base, ..) = display(&s);
        assert!(!base.surface_tris.is_empty());
        let part: Vec<Rgb> = model_scene(&base, &p).tris.iter().map(|t| t.color).collect();
        let (scene, bounds, shown) = super::scene(&s, &p);
        let (field, min, max) = shown.contour.clone().unwrap();
        assert!(field == "velocity" && max.1 > min.1, "{shown:?}");
        assert!(scene.tris.len() > base.surface_tris.len(), "the slice is drawn");
        assert!(!scene.lines.is_empty(), "streamlines are drawn");
        // The tank, in metres.
        let tank = Aabb3 { min: Vec3::new(-0.02, -0.02, -0.02), max: Vec3::new(0.06, 0.05, 0.03) };
        assert!(bounds.min.dist(tank.min) < 1e-12 && bounds.max.dist(tank.max) < 1e-12, "the picture frames the tank: {bounds:?}");
        // The part is coloured by the flow; in its own colours without the surface colours.
        let in_part_colour = |scene: &Scene| scene.tris.iter().filter(|t| part.contains(&t.color)).count();
        assert!(in_part_colour(&scene) < base.surface_tris.len());
        s.execute("results.surface", &json!({"show": false})).unwrap();
        s.execute("results.slice", &json!({"axis": "none"})).unwrap();
        s.execute("results.streamlines", &json!({"show": false})).unwrap();
        let (scene, ..) = super::scene(&s, &p);
        assert_eq!((in_part_colour(&scene), scene.tris.len()), (base.surface_tris.len(), base.surface_tris.len()), "the part alone");
        let (png, _) = render_png(&s, &Orientation::of(StandardView::Iso), 160, 120, &p).unwrap();
        assert!(count_color(&png, p.background.0) < 160 * 120, "the part is in the picture");
    }

    /// Pixels of colour `c` (within a few levels) in a PNG.
    fn count_color(png: &[u8], c: Rgb) -> usize {
        let d = png::Decoder::new(std::io::Cursor::new(png));
        let mut r = d.read_info().unwrap();
        let mut buf = vec![0; r.output_buffer_size().unwrap()];
        let info = r.next_frame(&mut buf).unwrap();
        let close = |a: u8, b: u8| a.abs_diff(b) <= 6;
        buf[..info.buffer_size()].chunks(4).filter(|p| close(p[0], c.0) && close(p[1], c.1) && close(p[2], c.2)).count()
    }

    /// The symbols are as large on screen in a close-up as in a distant view.
    #[test]
    fn symbols_keep_their_size_in_pictures_at_any_distance() {
        let mut s = Session::default();
        s.run_script(&crate::sample::cantilever()).unwrap();
        let at = |d: f64| orientation(&json!({"yaw": -45, "pitch": 30, "distance": d})).unwrap();
        let clamp = |d: f64| count_color(&render_png(&s, &at(d), 400, 300, &Palette::LIGHT).unwrap().0, Palette::LIGHT.glyphs.constraint);
        let (near, far) = (clamp(120.0), clamp(1200.0));
        // Zoomed out ten times, the symbols are as large (a few overlap and thin out).
        assert!(far > 30 && far * 4 > near, "near {near}, far {far}");
    }

    #[test]
    fn camera_orientation_is_checked() {
        let o = orientation(&json!({"yaw": -45, "pitch": 30})).unwrap();
        assert_eq!((o.yaw, o.pitch, o.distance), (-45.0, 30.0, None));
        let o = orientation(&json!({"yaw": 400, "pitch": -400, "distance": 1e300})).unwrap();
        assert_eq!((o.yaw, o.pitch, o.distance), (40.0, -90.0, Some(Orientation::DISTANCE.1)));
        for hostile in [
            json!(null),
            json!([1, 2]),
            json!({}),
            json!({"yaw": 0}),
            json!({"pitch": 0}),
            json!({"yaw": "left", "pitch": 0}),
            json!({"yaw": 0, "pitch": 0, "distance": 0}),
            json!({"yaw": 0, "pitch": 0, "distance": -5}),
            json!({"yaw": 0, "pitch": 0, "distance": "far"}),
            json!({"yaw": 0, "pitch": 0, "roll": 3}),
            json!({"yaw": f64::INFINITY, "pitch": 0}),
        ] {
            assert!(orientation(&hostile).is_err(), "{hostile}");
        }
        // A picture from a chosen orientation and distance differs from the framed iso one.
        let mut s = Session::default();
        s.run_script(&crate::sample::cantilever()).unwrap();
        let iso = render_png(&s, &Orientation::of(StandardView::Iso), 64, 48, &Palette::LIGHT).unwrap().0;
        let side = render_png(&s, &orientation(&json!({"yaw": -80, "pitch": 10, "distance": 2000})).unwrap(), 64, 48, &Palette::LIGHT).unwrap().0;
        assert_ne!(iso, side);
    }
}
