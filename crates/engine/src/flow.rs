//! The display of a fluid step's result: its fields per cell drawn as a coloured slice through the
//! domain, on the obstacles' surfaces beside the fluid, and as streamlines coloured by speed or temperature
//! ([`fieldcraft_render::flow`]), as the results view's [`FlowView`] asks. The headless pictures
//! ([`crate::view`]) and the desktop viewport draw the same [`FlowPlot`]; the part stays drawn
//! (a fluid step does not deform it).
//!
//! A fluid step is first shown with velocity magnitude on a slice through the middle of the
//! domain, normal to the horizontal axis across the flow (the slice holds the flow and the
//! vertical: Z is up), on the obstacles, and with streamlines seeded on the inlet face upstream
//! of the obstacles.

use std::sync::{Arc, Mutex, Weak};

use fieldcraft_lbm::{FluidResult, Frame};
use fieldcraft_render::contour::{Bands, banded, bands_color};
use fieldcraft_render::flow::{self, CellGrid, Slice, Streamline};
use fieldcraft_render::{Aabb3, Rgb, Scene, SceneLine, Vec3};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::results::{Field, Kind, Pick, ResultsView, magnitude};
use crate::{EngineError, Result, StepResult};

pub mod line;
pub mod wall;

/// Streamlines traced when no count is given, and the steps (half a cell each) each may take.
pub const DEFAULT_STREAMLINES: u32 = 64;
pub const DEFAULT_STREAMLINE_STEPS: u32 = 4_000;
/// Most streamlines and most steps of one, as commands accept them.
pub const MAX_STREAMLINES: u32 = flow::MAX_STREAMLINES as u32;
pub const MAX_STREAMLINE_STEPS: u32 = flow::MAX_STREAMLINE_STEPS as u32;
/// A streamline's width, in points.
const STREAMLINE_WIDTH: f64 = 1.5;

/// A coordinate axis: the normal of a slice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    pub const ALL: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

    /// 0, 1, 2.
    pub fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }

    pub fn id(self) -> &'static str {
        match self {
            Axis::X => "x",
            Axis::Y => "y",
            Axis::Z => "z",
        }
    }

    pub fn parse(id: &str) -> Option<Axis> {
        Axis::ALL.into_iter().find(|a| a.id().eq_ignore_ascii_case(id.trim()))
    }
}

/// A slice through the domain: normal to `axis` at `position` along it (`None`: the domain's
/// middle).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SliceView {
    pub axis: Axis,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub position: Option<f64>,
}

/// Streamlines: from `seeds`, or (`None`) `count` seeds on the inlet face upstream of the
/// obstacles; each at most `max_steps` steps of half a cell.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StreamlineView {
    pub count: u32,
    pub max_steps: u32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seeds: Option<Vec<[f64; 3]>>,
}

impl Default for StreamlineView {
    fn default() -> Self {
        StreamlineView { count: DEFAULT_STREAMLINES, max_steps: DEFAULT_STREAMLINE_STEPS, seeds: None }
    }
}

/// What a fluid step's display shows besides the field's legend (the view's contour field:
/// velocity, pressure, temperature or wall heat flux; `None` draws no colours).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FlowView {
    /// The slice, `None` for none.
    pub slice: Option<SliceView>,
    /// The obstacles' surfaces coloured by the field in the fluid beside them.
    pub surface: bool,
    pub streamlines: Option<StreamlineView>,
}

impl FlowView {
    /// A fluid result's first look: a slice through the middle across the flow, the surfaces
    /// coloured, streamlines from the inlet.
    pub fn for_result(r: &FluidResult) -> FlowView {
        let along = r.frames.last().map_or(Axis::X, |f| direction(r, f).0);
        // The slice holds the flow and the vertical (Z is up); for a flow along y, it is normal to x.
        let axis = if along == Axis::Y { Axis::X } else { Axis::Y };
        FlowView { slice: Some(SliceView { axis, position: None }), surface: true, streamlines: Some(StreamlineView::default()) }
    }

    /// What `view` shows of fluid result `r`: its own settings, or the result's first look.
    pub fn of<'a>(view: &'a ResultsView, r: &FluidResult) -> std::borrow::Cow<'a, FlowView> {
        match &view.flow {
            Some(f) => std::borrow::Cow::Borrowed(f),
            None => std::borrow::Cow::Owned(FlowView::for_result(r)),
        }
    }
}

/// The result's grid of cells.
pub fn grid(r: &FluidResult) -> CellGrid {
    CellGrid { origin: r.grid.origin, spacing: r.grid.spacing, cells: r.grid.cells }
}

/// The index of the frame of a fluid result `pick` asks for (`None`: the last).
pub fn frame_index(r: &FluidResult, pick: Pick) -> Result<usize> {
    let last = r.frames.len().checked_sub(1).ok_or_else(|| EngineError::Params(format!("step `{}` has no frames", r.step)))?;
    let k = pick.frame.unwrap_or(last);
    if k > last {
        return Err(EngineError::Params(format!("step `{}` has frames 0 to {last}; there is no frame {k}", r.step)));
    }
    Ok(k)
}

/// The frame of a fluid result `pick` asks for (`None`: the last).
pub fn frame(r: &FluidResult, pick: Pick) -> Result<&Frame> {
    let k = frame_index(r, pick)?;
    r.frames.get(k).ok_or_else(|| EngineError::Params(format!("step `{}` has no frame {k}", r.step)))
}

/// The axis the flow mostly runs along in a frame (by its mean velocity over the fluid cells) and
/// whether it runs toward its high end; along +x when the fluid is at rest.
pub fn direction(r: &FluidResult, f: &Frame) -> (Axis, bool) {
    let mut sum = [0.0f64; 3];
    for (k, v) in f.velocity.iter().enumerate() {
        if !r.solid.get(k).copied().unwrap_or(false) {
            for (s, c) in sum.iter_mut().zip(v) {
                if c.is_finite() {
                    *s += f64::from(*c);
                }
            }
        }
    }
    let axis = Axis::ALL.into_iter().fold(Axis::X, |best, a| if sum[a.index()].abs() > sum[best.index()].abs() { a } else { best });
    (axis, sum[axis.index()] >= 0.0)
}

/// A cell field's value in every cell of a frame: velocity magnitude, pressure or temperature, NaN in the
/// solid cells.
pub fn cell_values(r: &FluidResult, field: Field, f: &Frame) -> Result<Vec<f64>> {
    let solid = |k: usize| r.solid.get(k).copied().unwrap_or(false);
    let values: Vec<f64> = match field {
        Field::Velocity => f.velocity.iter().map(|v| magnitude(v.map(f64::from))).collect(),
        Field::Pressure => f.pressure.iter().map(|p| f64::from(*p)).collect(),
        Field::Temperature if f.temperature.len() == r.grid.len() && !f.temperature.is_empty() => r.temperatures(f),
        Field::Temperature => {
            return Err(EngineError::Params(format!(
                "temperature is not available for fluid step `{}`: this frame has no complete temperature field (isothermal flow)",
                r.step
            )));
        }
        _ => {
            return Err(EngineError::Params(format!(
                "{} is not available for fluid step `{}` (velocity, pressure or temperature)",
                field.label(),
                r.step
            )));
        }
    };
    Ok(values.into_iter().enumerate().map(|(k, v)| if solid(k) || !v.is_finite() { f64::NAN } else { v }).collect())
}

/// The smallest and largest value over the interior fluid cells, with their cells: the
/// open-boundary (inlet and outlet) layers are excluded because a wake reaching a pressure outlet
/// distorts the last layer (#590), which would otherwise report a spurious peak there.
pub fn interior_range(r: &FluidResult, values: &[f64]) -> Option<((usize, f64), (usize, f64))> {
    let inner: Vec<f64> = values.iter().enumerate().map(|(k, v)| if r.in_open_layer(k) { f64::NAN } else { *v }).collect();
    cell_range(&inner).or_else(|| cell_range(values))
}

/// The smallest and largest value, with their cells.
pub fn cell_range(values: &[f64]) -> Option<((usize, f64), (usize, f64))> {
    let mut out: Option<((usize, f64), (usize, f64))> = None;
    for (k, &v) in values.iter().enumerate().filter(|(_, v)| v.is_finite()) {
        out = Some(match out {
            None => ((k, v), (k, v)),
            Some((lo, hi)) => (if v < lo.1 { (k, v) } else { lo }, if v > hi.1 { (k, v) } else { hi }),
        });
    }
    out
}

/// The centre of cell `k` of a result's grid.
pub fn cell_centre(r: &FluidResult, k: usize) -> [f64; 3] {
    let g = grid(r);
    let c = g.centre(g.coords(k));
    [c.x, c.y, c.z]
}

/// Where a slice lies: its position, or the middle of the domain, held inside it.
pub fn slice_position(r: &FluidResult, s: &SliceView) -> f64 {
    let g = grid(r);
    let a = s.axis.index();
    let (lo, hi) = (g.start(a), g.end(a));
    s.position.filter(|p| p.is_finite()).unwrap_or(0.5 * (lo + hi)).clamp(lo.min(hi), hi.max(lo))
}

/// The seeds streamlines start from: the given ones, or `count` on the inlet face (where the
/// flow enters the domain) upstream of the obstacles.
pub fn seeds(r: &FluidResult, f: &Frame, s: &StreamlineView) -> Vec<Vec3> {
    match &s.seeds {
        Some(points) => points.iter().map(|p| Vec3::from(*p)).collect(),
        None => {
            let (axis, toward_high) = direction(r, f);
            flow::face_seeds(&grid(r), &r.solid, axis.index(), !toward_high, s.count as usize)
        }
    }
}

/// The streamlines a view asks for in a frame, and how many seeds they started from.
pub fn trace(r: &FluidResult, f: &Frame, s: &StreamlineView) -> (usize, Vec<Streamline>) {
    let seeds = seeds(r, f, s);
    (seeds.len(), flow::streamlines(&grid(r), &f.velocity, &r.solid, &seeds, s.max_steps as usize))
}

/// Streamlines traced once and shared: `results.streamlines` traces them for its reply, and the
/// plot drawn right after takes them from here instead of tracing them again (at the caps, a
/// trace is 2 000 000 points). One entry: the lines of one result, frame and settings, kept
/// while that result lives.
#[derive(Debug, Default)]
pub struct TraceCache(Mutex<Option<Traced>>);

#[derive(Clone, Debug)]
struct Traced {
    result: Weak<StepResult>,
    frame: usize,
    view: StreamlineView,
    seeds: usize,
    lines: Arc<Vec<Streamline>>,
}

impl Clone for TraceCache {
    fn clone(&self) -> Self {
        TraceCache(Mutex::new(self.lock().clone()))
    }
}

impl TraceCache {
    fn lock(&self) -> std::sync::MutexGuard<'_, Option<Traced>> {
        self.0.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// The streamlines `s` asks for in frame `k` of the fluid result `result`, and how many
    /// seeds they started from: the kept ones when they are those, else traced now and kept.
    pub fn trace(&self, result: &Arc<StepResult>, k: usize, s: &StreamlineView) -> Result<(usize, Arc<Vec<Streamline>>)> {
        let Kind::Fluid(r) = result.kind() else { return Err(EngineError::Params("streamlines are traced in fluid steps".into())) };
        let f = r.frames.get(k).ok_or_else(|| EngineError::Params(format!("step `{}` has no frame {k}", r.step)))?;
        if let Some(t) = self.lock().as_ref()
            && t.frame == k
            && t.view == *s
            && t.result.upgrade().is_some_and(|a| Arc::ptr_eq(&a, result))
        {
            return Ok((t.seeds, t.lines.clone()));
        }
        let (seeds, lines) = trace(r, f, s);
        let lines = Arc::new(lines);
        *self.lock() = Some(Traced { result: Arc::downgrade(result), frame: k, view: s.clone(), seeds, lines: lines.clone() });
        Ok((seeds, lines))
    }
}

/// A cell field ready to colour: its values, bands and extremes (with their cells).
#[derive(Clone, Debug, PartialEq)]
pub struct CellContour {
    pub field: Field,
    pub values: Vec<f64>,
    pub bands: Bands,
    pub min: (usize, f64),
    pub max: (usize, f64),
}

/// A fluid step's display, ready to draw ([`draw`]).
#[derive(Clone, Debug, PartialEq)]
pub struct FlowPlot {
    /// What it shows: the view's settings, or the result's first look.
    pub view: FlowView,
    pub grid: CellGrid,
    /// The field coloured, when there is one.
    pub contour: Option<CellContour>,
    /// The slice (with its normal axis and position), when there is one and a field to colour it.
    pub slice: Option<(Axis, f64, Slice)>,
    pub surface: bool,
    /// Shared with the [`TraceCache`] they came from.
    pub streamlines: Arc<Vec<Streamline>>,
    /// The seeds they were traced from (the given ones, or those the inlet face's grid fits).
    pub seeds: usize,
    /// The bands streamlines are coloured by speed over: the legend's when it shows velocity;
    /// `None` draws them in one colour (the legend shows another field).
    pub speed: Option<Bands>,
    /// Face-averaged heat values paired with the displayed geometry triangles.
    pub wall_values: Vec<([Vec3; 3], Option<f64>)>,
    /// A point on each recorded wall face for the legend's extreme markers.
    pub wall_positions: Vec<Option<Vec3>>,
}

impl FlowPlot {
    /// The box around what is drawn of the flow: the domain, with a slice or streamlines.
    pub fn bounds(&self) -> Aabb3 {
        if self.slice.is_some() || !self.streamlines.is_empty() { self.grid.bounds() } else { Aabb3::EMPTY }
    }
}

/// The display of the fluid result `result` the view asks for, its streamlines from `traced`.
pub fn plot(view: &ResultsView, result: &Arc<StepResult>, traced: &TraceCache) -> Result<FlowPlot> {
    let Kind::Fluid(r) = result.kind() else { return Err(EngineError::Params("the flow display shows fluid steps".into())) };
    let k = frame_index(r, view.pick())?;
    let f = frame(r, view.pick())?;
    let fv = FlowView::of(view, r).into_owned();
    let contour = match view.contour.filter(|c| matches!(c, Field::Velocity | Field::Pressure | Field::Temperature)) {
        Some(field) => {
            let values = cell_values(r, field, f)?;
            interior_range(r, &values).map(|(min, max)| CellContour { field, bands: view.bands(min.1, max.1), min, max, values })
        }
        None => None,
    };
    let g = grid(r);
    let slice = match (&fv.slice, &contour) {
        (Some(s), Some(c)) => {
            let at = slice_position(r, s);
            flow::slice(&g, &c.values, s.axis.index(), at).map(|sl| (s.axis, at, sl))
        }
        _ => None,
    };
    let (seeds, streamlines) = match &fv.streamlines {
        Some(s) => traced.trace(result, k, s)?,
        None => (0, Arc::default()),
    };
    let speed = contour.as_ref().filter(|c| c.field == Field::Velocity).map(|c| c.bands);
    Ok(FlowPlot {
        grid: g,
        contour,
        slice,
        surface: fv.surface,
        streamlines,
        seeds,
        speed,
        view: fv,
        wall_values: Vec::new(),
        wall_positions: Vec::new(),
    })
}

/// Add the obstacle face mapping shared by the desktop and headless views.
pub fn plot_with_surfaces(
    view: &ResultsView,
    result: &Arc<StepResult>,
    traced: &TraceCache,
    session: &crate::Session,
    base: &fieldcraft_render::DisplayMesh,
    bodies: &[&str],
) -> Result<FlowPlot> {
    let mut out = plot(view, result, traced)?;
    if view.contour == Some(Field::Temperature)
        && let Kind::Fluid(r) = result.kind()
    {
        // An obstacle face held at a temperature is drawn at it: sampling the fluid beside it
        // would mix in the cooler cells further out.
        let held = held_temperatures(session, r, base, bodies);
        if held.iter().any(Option::is_some) {
            if let Some(c) = &mut out.contour {
                // A held face beyond the fluid's range is the extreme: the legend, its Max/Min
                // lines and the viewport's marker all name it. Its label continues past the
                // cells, into `wall_positions` (see [`FlowPlot::wall_positions`]).
                let cells = c.values.len();
                for (t, v) in base.surface_tris.iter().zip(&held) {
                    let Some(v) = *v else { continue };
                    let at = || (t.p[0] + t.p[1] + t.p[2]) * (1.0 / 3.0);
                    if v > c.max.1 {
                        c.max = (cells + out.wall_positions.len(), v);
                        out.wall_positions.push(Some(at()));
                    } else if v < c.min.1 {
                        c.min = (cells + out.wall_positions.len(), v);
                        out.wall_positions.push(Some(at()));
                    }
                }
                c.bands = view.bands(c.min.1, c.max.1);
            }
            out.wall_values = base.surface_tris.iter().zip(held).map(|(t, v)| (t.p, v)).collect();
        }
    }
    if view.contour == Some(Field::WallHeatFlux) {
        let Kind::Fluid(r) = result.kind() else { return Ok(out) };
        let faces = wall::faces(session, r, view.pick())?;
        let positions: std::collections::BTreeMap<_, _> = base
            .surface_tris
            .iter()
            .filter_map(|t| bodies.get(t.body as usize).map(|name| ((*name, t.face), (t.p[0] + t.p[1] + t.p[2]) * (1.0 / 3.0))))
            .collect();
        out.wall_positions = faces.iter().map(|f| positions.get(&(f.body.as_str(), f.index)).copied()).collect();
        let values: Vec<f64> = faces.iter().map(|f| f.flux.unwrap_or(f64::NAN)).collect();
        out.contour =
            cell_range(&values).map(|(min, max)| CellContour { field: Field::WallHeatFlux, bands: view.bands(min.1, max.1), min, max, values });
        out.wall_values = base.surface_tris.iter().zip(wall::triangle_values(&faces, base, bodies)).map(|(t, v)| (t.p, v)).collect();
    }
    Ok(out)
}

/// Per displayed surface triangle, the temperature its obstacle face is held at, when the step's
/// domain holds that face at one and the result is of the current model.
fn held_temperatures(session: &crate::Session, r: &FluidResult, base: &fieldcraft_render::DisplayMesh, bodies: &[&str]) -> Vec<Option<f64>> {
    let faces = held_faces(session, r);
    let held: std::collections::BTreeMap<_, _> = faces.iter().map(|h| ((h.body.as_str(), h.index), h.value)).collect();
    base.surface_tris.iter().map(|t| bodies.get(t.body as usize).and_then(|name| held.get(&(*name, t.face)).copied())).collect()
}

/// An obstacle face the step's fluid domain holds at a temperature.
pub struct HeldFace {
    pub body: String,
    pub face: String,
    pub index: u32,
    pub value: f64,
}

/// The faces the step's fluid domain holds at a temperature, when the result is of the current
/// model (a stale result was solved with other conditions).
pub fn held_faces(session: &crate::Session, r: &FluidResult) -> Vec<HeldFace> {
    if session.latest_result(&r.step).is_some_and(|(_, current)| !current) {
        return Vec::new();
    }
    let Some(thermal) = session
        .model
        .step(&r.step)
        .and_then(|st| st.procedure.fluid_domain())
        .and_then(|d| session.model.fluid_domain(d))
        .and_then(|d| d.thermal.as_ref())
    else {
        return Vec::new();
    };
    let design = session.design.as_deref();
    let mut out = Vec::new();
    for sb in &thermal.surfaces {
        let (fieldcraft_model::ThermalCondition::Temperature { value }, fieldcraft_model::Entity::Face(face)) = (&sb.condition, &sb.face.entity)
        else {
            continue;
        };
        let index = design.and_then(|d| d.bodies.iter().find(|b| b.name == sb.face.body)).and_then(|b| b.face_index(face));
        if let Some(index) = index {
            out.push(HeldFace { body: sb.face.body.clone(), face: face.clone(), index, value: *value });
        }
    }
    out
}

/// The colours a flow plot is drawn in.
#[derive(Clone, Copy, Debug)]
pub struct FlowColors<'a> {
    /// The contour scale, low to high, and values above a clamped legend's maximum.
    pub stops: &'a [Rgb],
    pub over: Rgb,
    /// Streamlines when the legend shows neither velocity nor temperature.
    pub plain: Rgb,
}

/// Draw `plot` into `scene`, which holds the model as drawn: its faces beside the fluid
/// coloured by the field (when the surface is shown), then the slice and the streamlines.
pub fn draw(plot: &FlowPlot, scene: &mut Scene, colors: &FlowColors) {
    if let Some(c) = &plot.contour {
        let color = |band: Option<u32>| band.map_or(colors.plain, |k| bands_color(colors.stops, colors.over, &c.bands, k));
        if plot.surface && c.field == Field::WallHeatFlux {
            // Only tagged geometry triangles can receive wall heat. A structural mesh or
            // hidden geometry must not accidentally inherit a face's value by index.
            if scene.tris.len() == plot.wall_values.len() {
                for (tri, (points, value)) in scene.tris.iter_mut().zip(&plot.wall_values) {
                    if tri.p == *points
                        && let Some(v) = value
                    {
                        tri.color = color(c.bands.band(*v));
                    }
                }
            }
        } else if plot.surface {
            let model = std::mem::take(&mut scene.tris);
            // Faces held at a temperature (`wall_values`, matched by position as for wall heat).
            let held = (c.field == Field::Temperature && model.len() == plot.wall_values.len()).then_some(&plot.wall_values);
            for (k, t) in model.into_iter().enumerate() {
                if let Some((points, Some(v))) = held.and_then(|h| h.get(k))
                    && t.p == *points
                {
                    banded(t.p, [*v; 3], t.normal, &c.bands, &color, &mut scene.tris);
                    continue;
                }
                let v = t.p.map(|p| flow::sample(&plot.grid, &c.values, p));
                let (sum, known) = v.iter().filter(|x| !x.is_nan()).fold((0.0, 0u32), |(s, n), x| (s + x, n + 1));
                if known == 0 {
                    scene.tris.push(t);
                    continue;
                }
                let mean = sum / f64::from(known);
                banded(t.p, v.map(|x| if x.is_nan() { mean } else { x }), t.normal, &c.bands, &color, &mut scene.tris);
            }
        }
        if let Some((_, _, s)) = &plot.slice {
            for t in &s.tris {
                banded(t.p, t.v, s.normal, &c.bands, &color, &mut scene.tris);
            }
        }
    }
    for l in plot.streamlines.iter() {
        for (k, w) in l.points.windows(2).enumerate() {
            let (Some(&a), Some(&b)) = (w.first(), w.get(1)) else { continue };
            let speed = 0.5 * (l.speed.get(k).copied().unwrap_or(0.0) + l.speed.get(k + 1).copied().unwrap_or(0.0));
            let temperature = plot.contour.as_ref().filter(|c| c.field == Field::Temperature);
            let color = if let Some(c) = temperature {
                let value = flow::sample(&plot.grid, &c.values, (a + b) * 0.5);
                c.bands.band(value).map(|k| bands_color(colors.stops, colors.over, &c.bands, k))
            } else {
                plot.speed.as_ref().and_then(|b| b.band(speed).map(|k| bands_color(colors.stops, colors.over, b, k)))
            }
            .unwrap_or(colors.plain);
            scene.lines.push(SceneLine { a, b, color, width: STREAMLINE_WIDTH });
        }
    }
}

/// What a fluid step's display shows, for the `results.*` replies: the field's range with the
/// cells (and their centres) of its extremes, the slice and the streamlines' settings.
pub fn report(view: &ResultsView, r: &FluidResult, display: &fieldcraft_model::DisplaySystem) -> Value {
    use crate::results::{Quantity, shown, shown_all, unit_of};
    let mut out = json!({});
    let Ok(f) = frame(r, view.pick()) else { return out };
    out["time"] = json!(shown(Quantity::Time, f.time, display));
    if let Some(field) = view.contour
        && let Ok(values) = cell_values(r, field, f)
    {
        let q = field.quantity();
        let value = |v: f64| q.map_or(v, |q| shown(q, v, display));
        let end = |(k, v): (usize, f64)| json!({"cell": k, "at": shown_all(Quantity::Length, cell_centre(r, k), display), "value": value(v)});
        let range = interior_range(r, &values).map(|(lo, hi)| json!({"min": end(lo), "max": end(hi)}));
        let unit = q.and_then(|q| unit_of(q, display));
        out["contour"] = json!({"field": field, "label": field.label(), "unit": unit, "range": range});
        if view.range_min.is_some() || view.range_max.is_some() {
            out["contour"]["legend"] = json!({"min": view.range_min.map(value), "max": view.range_max.map(value)});
        }
    }
    let fv = FlowView::of(view, r);
    out["slice"] = match fv.slice.as_ref().filter(|_| view.contour != Some(Field::WallHeatFlux)) {
        Some(s) => json!({"axis": s.axis, "position": shown(Quantity::Length, slice_position(r, s), display)}),
        None => Value::Null,
    };
    out["surface"] = json!(fv.surface);
    out["streamlines"] = match &fv.streamlines {
        Some(s) => {
            let mut v = json!({"count": s.count, "max_steps": s.max_steps});
            match &s.seeds {
                Some(points) => v["seeds"] = json!(points.len()),
                None => {
                    let (axis, high) = direction(r, f);
                    v["inlet"] = json!(format!("{}_{}", axis.id(), if high { "min" } else { "max" }));
                }
            }
            v
        }
        None => Value::Null,
    };
    out
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use fieldcraft_lbm::{FluidStats, Grid};

    /// A synthetic fluid result: a 10 × 4 × 4 grid of unit cells, a solid 2 × 2 × 2 block at
    /// x, y, z in [5, 7] × [1, 3] × [1, 3], the fluid moving along +x at the speed of its x and
    /// the pressure falling along x.
    pub(crate) fn channel() -> FluidResult {
        let grid = Grid { origin: [0.0; 3], spacing: 1.0, cells: [10, 4, 4] };
        let n = grid.len();
        let solid: Vec<bool> = (0..n).map(|c| matches!(grid.coords(c), [5 | 6, 1 | 2, 1 | 2])).collect();
        let velocity: Vec<[f32; 3]> = (0..n).map(|c| if solid[c] { [0.0; 3] } else { [grid.coords(c)[0] as f32 + 1.0, 0.0, 0.0] }).collect();
        let pressure: Vec<f32> = (0..n).map(|c| if solid[c] { 0.0 } else { 10.0 - grid.coords(c)[0] as f32 }).collect();
        let at_rest = Frame { time: 0.0, step: 0, velocity: vec![[0.0; 3]; n], pressure: vec![0.0; n], temperature: Vec::new(), mass: 1.0 };
        let last = Frame { time: 1.0, step: 10, velocity, pressure, temperature: Vec::new(), mass: 1.0 };
        let stats = FluidStats {
            cells: n,
            fluid_cells: n - 8,
            steps: 10,
            time_step: 0.1,
            end_time: 1.0,
            tau_plus: 1.0,
            tau_minus: 1.0,
            magic: 0.25,
            collision: fieldcraft_model::Collision::Trt,
            velocity: 1.0,
            lattice_velocity: 0.1,
            mach: 0.1,
            cell_updates: 0,
            backend: "cpu".into(),
            tau_thermal: None,
            eddy_viscosity_ratio: None,
        };
        FluidResult {
            step: "Flow".into(),
            grid,
            solid,
            open_faces: [false; 6],
            frames: vec![at_rest, last],
            history_times: Vec::new(),
            forces: Vec::new(),
            heat: Vec::new(),
            zones: Vec::new(),
            solids: Vec::new(),
            stats,
            warnings: Vec::new(),
            steady: None,
            temperature_reference: 0.0,
        }
    }

    /// A 40 × 30 × 10 mm plate (the geometry sample) in a tank 20 mm larger on every side, of
    /// 5 mm cells, a liquid flowing along x past it: solved, so its result is shown.
    pub(crate) fn plate_in_a_tank() -> crate::Session {
        use fieldcraft_geometry::sample::Plate;
        let (_dir, path) = crate::tests::sample_design_file(&Plate::default());
        let mut s = crate::Session::default();
        s.execute("geometry.open", &json!({"path": path})).unwrap();
        let script = json!([
            {"command": "model.material.create", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-7}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Water", "min": [-20, -20, -20], "max": [60, 50, 30], "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [10, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "steps": 60, "frames": 2}}},
        ]);
        s.run_script(&script).unwrap();
        s.execute("solve.fluid", &json!({})).unwrap();
        s
    }

    /// The plate conducts heat and generates power in the tank's flow, through commands: the solve
    /// reports its temperatures, power and the heat its faces give the water; its faces take no
    /// thermal condition, its material must conduct, and a build without a GPU backend refuses
    /// `backend: "gpu"`.
    #[test]
    fn a_conducting_plate_reports_its_temperatures_and_heat() {
        let mut s = plate_in_a_tank();
        let face = s.design.as_ref().unwrap().bodies[0].faces()[0].name.clone();
        let script = json!([
            {"command": "model.material.edit", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-7}, "thermal": {"conductivity": 400, "specific_heat": 4e9}}},
            {"command": "model.material.create", "params": {"name": "Alu", "density": 2.7e-9, "thermal": {"conductivity": 237, "specific_heat": 9e8}}},
            {"command": "model.material.create", "params": {"name": "Steel", "density": 7.8e-9}},
        ]);
        s.run_script(&script).unwrap();
        let domain = |solids: Value| {
            json!({"name": "Tank", "material": "Water", "min": [-20, -20, -20], "max": [60, 50, 30], "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [10, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}},
                "thermal": {"initial_temperature": 20, "solids": solids}})
        };
        s.execute("model.fluid_domain.edit", &domain(json!([{"body": "Body1", "material": "Alu", "power": 2e5}]))).unwrap();
        // typed in the display system (mm-t-s: mW), stored in SI (W); the reply converts back
        let stored = s.model.fluid_domain("Tank").and_then(|d| d.solids().next()).map(|c| c.power);
        assert!(stored.is_some_and(|p| (p - 200.0).abs() < 1e-9), "{stored:?}");
        let r = s.execute("solve.fluid", &json!({})).unwrap();
        let solid = &r["solids"][0];
        assert_eq!(solid["body"], "Body1", "{r}");
        assert!(solid["cells"].as_u64().unwrap() > 0, "{solid}");
        assert!((solid["power"].as_f64().unwrap() - 2e5).abs() < 1e-6 * 2e5, "{solid}");
        let (t, given) = (&solid["temperature"], solid["heat_to_fluid"]["final"].as_f64().unwrap());
        let [lo, mean, hi] = ["min", "mean", "max"].map(|k| t[k].as_f64().unwrap());
        assert!(20.0 <= lo && lo <= mean && mean <= hi && hi > 20.0, "{solid}");
        assert!(given > 0.0, "the heated plate warms the water: {solid}");
        // no limit: a peak (this transient run's highest over its frames), no verdict
        let peak = solid["peak"].as_f64().unwrap();
        assert!(hi <= peak && solid["limit"].is_null() && solid["ok"].is_null() && r["components_ok"].is_null(), "{solid}");
        // A limit typed in the display system (°C), stored in kelvin, replied in °C; passing and
        // failing either side of the peak, the margin a temperature difference.
        for (limit, ok) in [(peak + 1.0, true), (peak - 1.0, false)] {
            s.execute("model.fluid_domain.edit", &domain(json!([{"body": "Body1", "material": "Alu", "power": 2e5, "limit": limit}]))).unwrap();
            let stored = s.model.fluid_domain("Tank").and_then(|d| d.solids().next()).and_then(|c| c.limit).unwrap();
            assert!((stored - (limit + 273.15)).abs() < 1e-9, "{stored}");
            let r = s.execute("solve.fluid", &json!({})).unwrap();
            let solid = &r["solids"][0];
            assert!((solid["limit"].as_f64().unwrap() - limit).abs() < 1e-9, "{solid}");
            assert!((solid["margin"].as_f64().unwrap() - (limit - solid["peak"].as_f64().unwrap())).abs() < 1e-9, "{solid}");
            assert_eq!((&solid["ok"], &r["components_ok"]), (&json!(ok), &json!(ok)), "{r}");
            assert_eq!((r["units"]["temperature"].as_str(), r["units"]["temperature_difference"].as_str()), (Some("°C"), Some("K")), "{r}");
        }
        // saved and read back: still kelvin (the file is SI)
        let back = fieldcraft_io::read(&s.file_bytes(None).unwrap()).unwrap();
        let kept = back.model.fluid_domain("Tank").and_then(|d| d.solids().next()).and_then(|c| c.limit).unwrap();
        assert!((kept - (peak - 1.0 + 273.15)).abs() < 1e-9, "{kept}");
        // the focused command sets and clears it; a body that is not a conducting solid has none
        s.execute("model.fluid_domain.limit", &json!({"domain": "Tank", "body": "Body1", "limit": "90 °C"})).unwrap();
        assert_eq!(s.model.fluid_domain("Tank").and_then(|d| d.solids().next()).and_then(|c| c.limit), Some(363.15));
        s.execute("model.fluid_domain.limit", &json!({"domain": "Tank", "body": "Body1"})).unwrap();
        assert_eq!(s.model.fluid_domain("Tank").and_then(|d| d.solids().next()).and_then(|c| c.limit), None);
        let e = s.execute("model.fluid_domain.limit", &json!({"domain": "Tank", "body": "Body2", "limit": 90})).unwrap_err().to_string();
        assert!(e.contains("`Body2` is not a conducting solid"), "{e}");
        assert!(s.execute("model.fluid_domain.limit", &json!({"domain": "Tank", "body": "Body1", "limit": "-300 °C"})).is_err());
        // a thermal condition on a conducting face, a material that does not conduct, the GPU
        assert!(
            s.execute(
                "model.fluid_domain.thermal_boundary",
                &json!({"domain": "Tank", "at": {"body": "Body1", "face": face}, "condition": {"type": "temperature", "value": 80}}),
            )
            .is_err()
        );
        assert!(s.execute("model.fluid_domain.edit", &domain(json!([{"body": "Body1", "material": "Steel"}]))).is_err());
        assert!(s.execute("model.fluid_domain.edit", &domain(json!([{"body": "Nope", "material": "Alu"}]))).is_err());
        assert!(s.execute("model.fluid_domain.edit", &domain(json!([{"body": "Body1", "material": "Alu", "power": "nan"}]))).is_err());
        // (this build registers no GPU backend; lbm-gpu's tests run the solids on the device)
        let e = s.execute("solve.fluid", &json!({"backend": "gpu"})).unwrap_err().to_string();
        assert!(e.contains("GPU") && e.contains("cpu"), "{e}");
        // a thermal contact's resistance, or its layer, typed in the display system (mm-t-s:
        // mm²·K/mW, mm, mW/(mm·K)) or with its unit, stored in SI: 1e-5 K·m²/W each time
        let two = |contact: Value| {
            let mut d = domain(json!([{"body": "Body1", "material": "Alu"}, {"body": "Lid", "material": "Alu"}]));
            d["obstacles"] = json!(["Body1", "Lid"]);
            d["thermal"]["contacts"] = json!([contact]);
            d
        };
        for contact in [
            json!({"bodies": ["Body1", "Lid"], "resistance": 0.01}),
            json!({"bodies": ["Body1", "Lid"], "resistance": "10 mm²·K/W"}),
            json!({"bodies": ["Lid", "Body1"], "thickness": 0.05, "conductivity": 5}),
        ] {
            s.execute("model.fluid_domain.edit", &two(contact.clone())).unwrap();
            let r = s.model.fluid_domain("Tank").and_then(|d| d.thermal.as_ref()).and_then(|t| t.contacts.first()).map(|c| c.resistance());
            assert!(r.is_some_and(|r| (r / 1e-5 - 1.0).abs() < 1e-12), "{contact}: {r:?}");
        }
        let e = s.execute("model.fluid_domain.edit", &two(json!({"bodies": ["Body1", "Lid"], "resistance": 0}))).unwrap_err().to_string();
        assert!(e.contains("thermal contact `Body1–Lid`") && e.contains("positive"), "{e}");
    }

    /// A face held at a temperature is drawn at that temperature, and the legend reaches it,
    /// rather than at the mean of the fluid cells around it (which mixes in the cooler ones
    /// further out: an 80 °C bracket looked 50 °C).
    #[test]
    fn a_face_held_at_a_temperature_is_drawn_at_it() {
        let mut s = plate_in_a_tank();
        let face = s.design.as_ref().unwrap().bodies[0].faces()[0].name.clone();
        s.execute(
            "model.material.edit",
            &json!({"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-7}, "thermal": {"conductivity": 400, "specific_heat": 4e9}}),
        )
        .unwrap();
        s.execute(
            "model.fluid_domain.edit",
            &json!({"name": "Tank", "material": "Water", "min": [-20, -20, -20], "max": [60, 50, 30], "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [10, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}},
                "thermal": {"initial_temperature": 20}}),
        )
        .unwrap();
        s.execute(
            "model.fluid_domain.thermal_boundary",
            &json!({"domain": "Tank", "at": {"body": "Body1", "face": face}, "condition": {"type": "temperature", "value": 80}}),
        )
        .unwrap();
        s.execute("solve.fluid", &json!({})).unwrap();
        s.execute("results.surface", &json!({"field": "temperature", "show": true})).unwrap();
        let (result, current) = s.latest_result("Flow").unwrap();
        assert!(current);
        let view = s.results_view.clone().unwrap();
        let base = crate::view::display(&s).0;
        let names: Vec<&str> = s.design.as_ref().unwrap().bodies.iter().map(|b| b.name.as_str()).collect();
        let fp = plot_with_surfaces(&view, result, &s.traced, &s, &base, &names).unwrap();
        let c = fp.contour.as_ref().unwrap();
        // The fluid never quite reaches the wall's 80 °C; the legend does, and its maximum is the
        // held face (a label past the cells, at a point on that face). The plot is in kelvin.
        let hot = fieldcraft_model::DisplaySystem::default().to_si(fieldcraft_model::units::Quantity::Temperature, 80.0).unwrap();
        assert!(c.values.iter().all(|t| t.is_nan() || *t < hot));
        assert_eq!((c.max.1, c.bands.max), (hot, hot), "{:?}", c.bands);
        let at = c.max.0.checked_sub(c.values.len()).and_then(|k| fp.wall_positions.get(k).copied().flatten());
        assert!(at.is_some(), "{:?}", c.max);
        let held: Vec<f64> = fp.wall_values.iter().filter_map(|(_, v)| *v).collect();
        assert!(!held.is_empty() && held.iter().all(|t| *t == hot), "{held:?}");
        // A script reads the same maximum as the app shows: the held face.
        let reply = s.results_report();
        assert_eq!(reply["contour"]["range"]["max"], json!({"body": "Body1", "face": face, "value": 80.0}), "{reply}");
        // The headless picture draws the held face in the top band's colour.
        let p = crate::view::Palette::LIGHT;
        let top = bands_color(&p.contour, p.contour_over, &c.bands, c.bands.band(hot).unwrap());
        let (scene, ..) = crate::view::scene(&s, &p);
        assert!(scene.tris.iter().any(|t| t.color == top));
    }

    #[test]
    fn a_fluid_result_is_first_shown_across_the_flow() {
        let r = channel();
        assert_eq!(direction(&r, &r.frames[1]), (Axis::X, true));
        let v = FlowView::for_result(&r);
        assert_eq!(v.slice, Some(SliceView { axis: Axis::Y, position: None }));
        assert!(v.surface && v.streamlines.is_some());
        assert_eq!(slice_position(&r, v.slice.as_ref().unwrap()), 2.0);
        assert_eq!(slice_position(&r, &SliceView { axis: Axis::Z, position: Some(99.0) }), 4.0);
        assert_eq!(slice_position(&r, &SliceView { axis: Axis::Z, position: Some(f64::NAN) }), 2.0);
        // A flow along -y is sliced normal to x and seeded on the y_max face.
        let mut r = channel();
        for v in &mut r.frames[1].velocity {
            *v = [0.0, -v[0], 0.0];
        }
        assert_eq!(direction(&r, &r.frames[1]), (Axis::Y, false));
        assert_eq!(FlowView::for_result(&r).slice.unwrap().axis, Axis::X);
        let s = seeds(&r, &r.frames[1], &StreamlineView { count: 4, ..StreamlineView::default() });
        assert!(s.iter().all(|p| p.y == 3.5), "{s:?}");
    }

    #[test]
    fn temperature_colours_cells_and_streamlines_and_masks_solids() {
        let mut r = channel();
        for f in &mut r.frames {
            f.temperature = (0..r.grid.len()).map(|k| 20.0 + 2.0 * r.grid.coords(k)[0] as f32).collect();
        }
        let values = cell_values(&r, Field::Temperature, &r.frames[1]).unwrap();
        assert_eq!(values.iter().filter(|v| v.is_nan()).count(), 8);
        assert_eq!(cell_range(&values).map(|(a, b)| (a.1, b.1)), Some((20.0, 38.0)));
        let result = Arc::new(crate::StepResult::Fluid(r));
        let mut view = ResultsView::for_result(&result);
        view.contour = Some(Field::Temperature);
        let p = plot(&view, &result, &TraceCache::default()).unwrap();
        let mut scene = Scene::default();
        let colors = FlowColors { stops: &[Rgb(0, 0, 255), Rgb(255, 0, 0)], over: Rgb(255, 255, 255), plain: Rgb(7, 7, 7) };
        draw(&p, &mut scene, &colors);
        assert!(!scene.tris.is_empty());
        assert!(!scene.lines.is_empty());
        assert!(scene.lines.iter().all(|l| l.color != colors.plain));
        assert!(scene.lines.windows(2).any(|l| l[0].color != l[1].color));
        let cold = channel();
        assert!(cell_values(&cold, Field::Temperature, &cold.frames[1]).unwrap_err().to_string().contains("isothermal"));
    }

    #[test]
    fn cell_fields_leave_the_solid_out() {
        let r = channel();
        let f = &r.frames[1];
        let speed = cell_values(&r, Field::Velocity, f).unwrap();
        assert_eq!(speed.iter().filter(|v| v.is_nan()).count(), 8);
        let ((lo, a), (hi, b)) = cell_range(&speed).unwrap();
        assert_eq!((a, b, r.grid.coords(lo)[0], r.grid.coords(hi)[0]), (1.0, 10.0, 0, 9));
        assert_eq!(cell_centre(&r, hi), [9.5, 0.5, 0.5]);
        let p = cell_values(&r, Field::Pressure, f).unwrap();
        assert_eq!(cell_range(&p).map(|(lo, hi)| (lo.1, hi.1)), Some((1.0, 10.0)));
        assert!(cell_values(&r, Field::VonMises, f).is_err());
        assert_eq!(cell_range(&[f64::NAN]), None);
        assert!(frame(&r, Pick { mode: 1, frame: Some(2) }).is_err());
        assert_eq!(frame(&r, Pick::default()).unwrap().time, 1.0);
    }

    /// A trace is kept for its result, frame and settings: the reply and the plot share one;
    /// another frame, other settings or another result trace anew.
    #[test]
    fn a_trace_is_shared_until_its_result_frame_or_settings_change() {
        let result = Arc::new(crate::StepResult::Fluid(channel()));
        let view = ResultsView::for_result(&result);
        let traced = TraceCache::default();
        let s = StreamlineView::default();
        let (seeds, a) = traced.trace(&result, 1, &s).unwrap();
        assert!(seeds > 0 && !a.is_empty());
        let (_, b) = traced.trace(&result, 1, &s).unwrap();
        assert!(Arc::ptr_eq(&a, &b), "the same trace again");
        let p = plot(&view, &result, &traced).unwrap();
        assert!(Arc::ptr_eq(&a, &p.streamlines), "the plot draws the kept trace");
        let (_, c) = traced.trace(&result, 0, &s).unwrap();
        assert!(!Arc::ptr_eq(&a, &c), "another frame: another trace");
        let (_, d) = traced.trace(&result, 1, &StreamlineView { count: 3, ..s.clone() }).unwrap();
        assert!(!Arc::ptr_eq(&a, &d) && d.len() < a.len(), "other settings: another trace ({} lines, then {})", a.len(), d.len());
        let other = Arc::new(crate::StepResult::Fluid(channel()));
        let (_, e) = traced.trace(&other, 1, &StreamlineView { count: 3, ..s.clone() }).unwrap();
        assert!(!Arc::ptr_eq(&d, &e), "another result");
        assert!(traced.trace(&result, 2, &s).is_err(), "no frame 2");
        let cloned = traced.clone();
        assert!(Arc::ptr_eq(&cloned.trace(&other, 1, &StreamlineView { count: 3, ..s }).unwrap().1, &e), "a clone keeps the trace");
        let solid = Arc::new(crate::StepResult::Fluid(channel()));
        let mut v = view.clone();
        v.step = "Other".into();
        assert!(plot(&v, &solid, &traced).is_ok());
    }

    /// The plot colours the slice, the faces beside the fluid and the streamlines by speed.
    #[test]
    fn the_plot_draws_the_slice_surface_and_streamlines() {
        let r = channel();
        let result = Arc::new(crate::StepResult::Fluid(r.clone()));
        let mut view = ResultsView::for_result(&result);
        assert_eq!(view.contour, Some(Field::Velocity));
        let traced = TraceCache::default();
        let p = plot(&view, &result, &traced).unwrap();
        let (axis, at, slice) = p.slice.as_ref().unwrap();
        assert_eq!((*axis, *at), (Axis::Y, 2.0));
        assert!(!slice.tris.is_empty());
        assert!(!p.streamlines.is_empty() && p.speed.is_some());
        assert_eq!(p.bounds(), p.grid.bounds());
        // A face of the block, at x = 5, and one far off the grid.
        let tri = |x: f64| fieldcraft_render::SceneTri {
            p: [Vec3::new(x, 1.0, 1.0), Vec3::new(x, 3.0, 1.0), Vec3::new(x, 1.0, 3.0)],
            normal: Vec3::X,
            color: Rgb(1, 2, 3),
        };
        let mut scene = Scene { tris: vec![tri(5.0), tri(500.0)], ..Scene::default() };
        let stops = [Rgb(0, 0, 255), Rgb(255, 0, 0)];
        let colors = FlowColors { stops: &stops, over: Rgb(9, 9, 9), plain: Rgb(7, 7, 7) };
        draw(&p, &mut scene, &colors);
        assert!(scene.tris.iter().any(|t| t.color == Rgb(1, 2, 3)), "the face off the grid keeps its colour");
        let coloured = scene.tris.iter().filter(|t| t.color != Rgb(1, 2, 3)).count();
        assert!(coloured > slice.tris.len(), "the slice and the block's face are coloured");
        assert_eq!(scene.lines.len(), p.streamlines.iter().map(|l| l.points.len() - 1).sum::<usize>());
        assert!(scene.lines.iter().all(|l| l.color != colors.plain));
        // Pressure: streamlines in the plain colour; no field: no slice, no colours.
        view.contour = Some(Field::Pressure);
        let p = plot(&view, &result, &traced).unwrap();
        let mut scene = Scene::default();
        draw(&p, &mut scene, &colors);
        assert!(scene.lines.iter().all(|l| l.color == colors.plain));
        view.contour = None;
        let p = plot(&view, &result, &traced).unwrap();
        assert!(p.slice.is_none() && p.contour.is_none() && !p.streamlines.is_empty());
        view.flow = Some(FlowView { slice: None, surface: false, streamlines: None });
        view.contour = Some(Field::Velocity);
        let p = plot(&view, &result, &traced).unwrap();
        assert!(p.slice.is_none() && p.streamlines.is_empty() && p.bounds().is_empty());
    }
}
