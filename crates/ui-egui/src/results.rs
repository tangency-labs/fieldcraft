//! The results display in the viewport: what the engine's [`ResultsView`] asks for, drawn on the
//! model's display data. The deformed shape moves every element by its nodes' displacements
//! (or a mode shape) times the scale; contours split the faces into flat-coloured bands
//! ([`fieldcraft_render::contour`]); animation steps the deformation through a fixed number of
//! frames per cycle (an explicit step's animation runs through its field frames, contours
//! following). A plot of a stale result (the model changed since the solve) is drawn the
//! same way and marked. A fluid step draws the part as it is, with the engine's flow display
//! ([`fieldcraft_engine::flow`]: a slice, the obstacles coloured and streamlines) through it.

use std::sync::Arc;

use fieldcraft_engine::flow::FlowPlot;
use fieldcraft_engine::results::{Field, Kind, Pick, ResultsView, displacements, field_values, range};
use fieldcraft_engine::{Session, StepResult};
use fieldcraft_render::contour::{Bands, banded, bands_color};
use fieldcraft_render::{DisplayMesh, Scene, Vec3};

use crate::Style;
use crate::theme::{Tokens, rgb};

/// Animation frames per cycle, and the cycle's length in seconds.
pub const FRAMES: u64 = 24;
pub const PERIOD: f64 = 2.0;
/// Significant digits of the numbers the results show.
const DIGITS: i32 = 4;
/// A value below this fraction of the largest it is read with is round-off, shown as 0.
const ROUND_OFF: f64 = 1e-9; // unit-factors: relative to the largest value
/// Magnitudes written as plain decimals; others in scientific notation.
const PLAIN: std::ops::Range<f64> = 1e-3..1e6; // unit-factors: magnitudes of a written number

/// A contour plot: its field, bands, extremes and the value at every node.
#[derive(Clone, Debug, PartialEq)]
pub struct Contour {
    pub field: Field,
    pub bands: Bands,
    /// The smallest and largest value and their nodes (a fluid step's: their cells).
    pub min: (u32, f64),
    pub max: (u32, f64),
    /// The value at every node, in the result's node order (empty for a fluid step's cells).
    pub values: Vec<f64>,
}

impl Contour {
    /// The value at a node (NaN when it has none).
    pub fn at(&self, nodes: &[u32], label: u32) -> f64 {
        nodes.binary_search(&label).ok().and_then(|i| self.values.get(i).copied()).unwrap_or(f64::NAN)
    }
}

/// A results plot ready to draw.
#[derive(Clone, Debug)]
pub struct Plot {
    pub view: ResultsView,
    pub result: Arc<StepResult>,
    /// Whether the model changed since the step was solved.
    pub stale: bool,
    /// The deformation scale factor, before animation: the view's or the automatic one, also
    /// while the view is undeformed (what ticking Deformed would draw it with).
    pub scale: f64,
    /// The elements as drawn (deformed at the current animation frame).
    pub display: DisplayMesh,
    pub contour: Option<Arc<Contour>>,
    /// A fluid step's slice, coloured obstacles and streamlines.
    pub flow: Option<Arc<FlowPlot>>,
    /// The mode and field frame drawn.
    pub pick: Pick,
    /// The display units, which the numbers shown with the plot (the contour's too) are in.
    pub units: fieldcraft_model::DisplaySystem,
}

/// The part of a plot that does not change with the animation frame.
#[derive(Clone, Debug, Default)]
pub struct Prepared {
    pub scale: f64,
    pub contour: Option<Arc<Contour>>,
    pub flow: Option<Arc<FlowPlot>>,
}

/// A mode's value as the results show it: its frequency, or a buckling mode's load factor.
pub fn mode_value(buckling: bool, eigenvalue: f64, frequency: f64, units: &fieldcraft_model::DisplaySystem) -> String {
    use fieldcraft_engine::results::{Quantity, shown, unit_of};
    if buckling {
        format!("λ = {}", number(eigenvalue))
    } else {
        with_unit(number(shown(Quantity::Frequency, frequency, units)), unit_of(Quantity::Frequency, units).as_deref())
    }
}

impl Plot {
    /// A frequency or buckling step's mode shown, and its frequency (Hz) or load factor.
    pub fn mode(&self) -> Option<(usize, String)> {
        match self.result.kind() {
            Kind::Frequency(f) => f
                .modes
                .get(self.pick.mode.checked_sub(1)?)
                .map(|m| (self.pick.mode, mode_value(f.buckling.is_some(), m.eigenvalue, m.frequency, &self.units))),
            Kind::Static(_) | Kind::Explicit(_) | Kind::Fluid(_) | Kind::Heat(_) => None,
        }
    }

    /// An explicit or transient heat step's frame shown, its time in the display units (a nonlinear static step's:
    /// its load factor) and the number of frames.
    pub fn frame(&self) -> Option<(usize, f64, usize)> {
        let (n, time): (usize, &dyn Fn(usize) -> Option<f64>) = match self.result.kind() {
            Kind::Explicit(e) => (e.frames.len(), &|k| e.frames.get(k).map(|f| f.time)),
            Kind::Static(r) if r.nonlinear.is_some() => {
                let frames = r.nonlinear.as_ref().map(|h| h.frames.as_slice()).unwrap_or_default();
                (frames.len(), &|k| frames.get(k).map(|f| f.load_factor))
            }
            Kind::Heat(h) if h.times.len() > 1 => (h.times.len(), &|k| h.times.get(k).copied()),
            Kind::Fluid(f) => (f.frames.len(), &|k| f.frames.get(k).map(|f| f.time)),
            _ => return None,
        };
        let k = self.pick.frame.unwrap_or(n.checked_sub(1)?);
        let at = time(k)?;
        let at = if matches!(self.result.kind(), Kind::Static(_)) {
            at
        } else {
            fieldcraft_engine::results::shown(fieldcraft_engine::results::Quantity::Time, at, &self.units)
        };
        Some((k, at, n))
    }

    /// The node's position in the plot, if it is drawn.
    pub fn node_position(&self, label: u32) -> Option<Vec3> {
        self.display.node_position(label)
    }

    /// Where a contour extreme is: its node's position, or a fluid step's cell centre.
    pub fn extreme_position(&self, label: u32) -> Option<Vec3> {
        if self.view.contour == Some(Field::WallHeatFlux) {
            return self.flow.as_ref()?.wall_positions.get(label as usize).copied().flatten();
        }
        match self.result.kind() {
            Kind::Fluid(f) => {
                let label = usize::try_from(label).ok()?;
                // Past the cells: a face held at a temperature that is the field's extreme.
                match label.checked_sub(f.grid.len()) {
                    Some(k) => self.flow.as_ref()?.wall_positions.get(k).copied().flatten(),
                    None => Some(Vec3::from(fieldcraft_engine::flow::cell_centre(f, label))),
                }
            }
            _ => self.node_position(label),
        }
    }

    /// Whether it shows a fluid step (fields per cell, no deformation).
    pub fn is_fluid(&self) -> bool {
        matches!(self.result.kind(), Kind::Fluid(_))
    }
}

/// The result the session's view shows, and whether it is current.
pub fn shown(session: &Session) -> Option<(&ResultsView, &Arc<StepResult>, bool)> {
    let view = session.results_view.as_ref()?;
    let (result, current) = session.latest_result(&view.step)?;
    Some((view, result, current))
}

/// The animation frame shown at time `now` (seconds), when animating. A deformed step swings its
/// shape; a fluid step plays its frames though it does not deform; a step with nothing to animate
/// (Animate on but Deformed off, and not fluid) stays still.
pub fn frame(view: &ResultsView, result: &StepResult, now: f64) -> Option<u64> {
    let plays = view.deformed || matches!(result.kind(), Kind::Fluid(_));
    (view.animate && plays && now.is_finite()).then(|| ((now.max(0.0) / PERIOD * FRAMES as f64).floor() as u64) % FRAMES)
}

/// The deformation multiplier of an animation frame: a mode swings through ±1, a static step
/// grows from 0 to 1 and back; an explicit step's frames are drawn as they are.
pub fn swing(result: &StepResult, frame: Option<u64>) -> f64 {
    let Some(k) = frame else { return 1.0 };
    let a = std::f64::consts::TAU * k as f64 / FRAMES as f64;
    match result.kind() {
        Kind::Frequency(_) => a.sin(),
        // A nonlinear step plays its increments' frames.
        Kind::Static(r) if r.nonlinear.is_some() => 1.0,
        Kind::Static(_) => 0.5 * (1.0 - a.cos()),
        Kind::Explicit(_) | Kind::Fluid(_) | Kind::Heat(_) => 1.0,
    }
}

/// What is drawn at animation `frame`: the view's mode and field frame, an animated explicit
/// step running through its field frames once per cycle.
pub fn pick(view: &ResultsView, result: &StepResult, frame: Option<u64>) -> Pick {
    let mut p = view.pick();
    let frames = match result.kind() {
        Kind::Explicit(e) => e.frames.len(),
        Kind::Fluid(f) => f.frames.len(),
        Kind::Static(r) => r.nonlinear.as_ref().map_or(0, |h| h.frames.len()),
        _ => 0,
    };
    if let (Some(k), Some(last)) = (frame, frames.checked_sub(1)) {
        p.frame = Some((k as usize * frames / FRAMES as usize).min(last));
    }
    p
}

/// The deformation scale factor: the view's, or automatic ([`auto_scale`]), whether or not the
/// view is deformed.
pub fn scale(view: &ResultsView, result: &StepResult, pick: Pick, base: &DisplayMesh, units: &fieldcraft_model::DisplaySystem) -> f64 {
    fieldcraft_engine::results::shown_scale(view, result, pick, base.bounds.diagonal(), units)
}

/// The part of the session's results plot that depends on what is drawn (`pick`), not on the
/// animation's deformation multiplier: the deformation scale and the contour values.
pub fn prepare(session: &Session, base: &DisplayMesh, pick: Pick, bodies: &[&str]) -> Option<Prepared> {
    let (view, result, _) = shown(session)?;
    if matches!(result.kind(), Kind::Fluid(_)) {
        let mut view = view.clone();
        view.frame = pick.frame;
        let flow = fieldcraft_engine::flow::plot_with_surfaces(&view, result, &session.traced, session, base, bodies).ok();
        // The legend's numbers in the display units (the cells are coloured in SI: the same
        // bands, as the conversion keeps the order and the spacing of values).
        let units = session.units();
        let cell = |field: Field, (k, v): (usize, f64)| (u32::try_from(k).unwrap_or(u32::MAX), field.shown(result, v, units));
        let contour = flow.as_ref().and_then(|p| p.contour.as_ref()).map(|c| {
            let bands = Bands { min: c.field.shown(result, c.bands.min, units), max: c.field.shown(result, c.bands.max, units), ..c.bands };
            Arc::new(Contour { field: c.field, bands, min: cell(c.field, c.min), max: cell(c.field, c.max), values: Vec::new() })
        });
        return Some(Prepared { scale: 0.0, contour, flow: flow.map(Arc::new) });
    }
    let units = session.units();
    let contour = view.contour.and_then(|field| {
        let mut values = field_values(&session.model, result, field, pick).ok()?;
        // In the display units, the legend's range too (stored in SI).
        for v in &mut values {
            *v = field.shown(result, *v, units);
        }
        // The range (and the bands) follow what is drawn: a hidden body's nodes have no value.
        for (v, label) in values.iter_mut().zip(result.nodes()) {
            if base.node_position(*label).is_none() {
                *v = f64::NAN;
            }
        }
        let (min, max) = range(result.nodes(), &values)?;
        let shown = |v: Option<f64>| v.map(|v| field.shown(result, v, units));
        let bands = Bands::clamped(min.1, max.1, shown(view.range_min), shown(view.range_max), view.bands);
        Some(Arc::new(Contour { field, bands, min, max, values }))
    });
    Some(Prepared { scale: scale(view, result, pick, base, units), contour, flow: None })
}

/// The plot of the session's results view on the model's display data at animation `frame`,
/// with its [`prepare`]d part.
pub fn plot(session: &Session, base: &DisplayMesh, frame: Option<u64>, prepared: &Prepared) -> Option<Plot> {
    let (view, result, current) = shown(session)?;
    let nodes = result.nodes();
    let pick = pick(view, result, frame);
    // The shown scale multiplies the shape as the display units show it.
    let unit = fieldcraft_engine::results::deformation_unit(result, session.units());
    let factor = if view.deformed { prepared.scale * unit * swing(result, frame) } else { 0.0 };
    let u = displacements(result, pick).unwrap_or(&[]);
    // A fluid step does not move the part: it stays drawn as it is, geometry included.
    let display = if matches!(result.kind(), Kind::Fluid(_)) {
        base.clone()
    } else {
        base.displaced(|label| {
            let i = nodes.binary_search(&label).ok()?;
            u.get(i).map(|d| Vec3::from(*d) * factor)
        })
    };
    Some(Plot {
        view: view.clone(),
        result: result.clone(),
        stale: !current,
        scale: prepared.scale,
        display,
        contour: prepared.contour.clone(),
        flow: prepared.flow.clone(),
        pick,
        units: session.units().clone(),
    })
}

/// The scene of a plot: the deformed elements drawn as the model is ([`scene::build`]), with
/// their faces split into contour bands when a field is shown; a fluid step's part as it is with
/// the flow drawn through it ([`fieldcraft_engine::flow::draw`]). Over a plot, hiding the mesh
/// hides the element edges: the result itself never vanishes with them.
pub fn scene(plot: &Plot, style: Style, show: crate::scene::Show, t: &Tokens) -> Scene {
    let (style, show) = if show.mesh { (style, show) } else { (Style::Shaded, crate::scene::Show { mesh: true, ..show }) };
    let mut s = crate::scene::build(&plot.display, style, show, t);
    if let Some(flow) = &plot.flow {
        let stops = t.contour_stops();
        let colors = fieldcraft_engine::flow::FlowColors { stops: &stops, over: rgb(t.contour_over), plain: rgb(t.streamline) };
        fieldcraft_engine::flow::draw(flow, &mut s, &colors);
        return s;
    }
    let Some(c) = plot.contour.as_ref().filter(|_| style != Style::Wireframe && show.mesh) else { return s };
    let nodes = plot.result.nodes();
    let stops = t.contour_stops();
    let color = |band: Option<u32>| band.map_or(rgb(t.no_result), |k| bands_color(&stops, rgb(t.contour_over), &c.bands, k));
    s.tris.clear();
    for tri in &plot.display.tris {
        banded(tri.p, tri.nodes.map(|l| c.at(nodes, l)), tri.normal, &c.bands, &color, &mut s.tris);
    }
    s
}

/// A number as the results show it (legend, probe, reaction sums, solve summary): four
/// significant digits, thousands separators, scientific notation when very large or small.
pub fn number(v: f64) -> String {
    if !v.is_finite() {
        return "—".into();
    }
    let a = v.abs();
    if a == 0.0 {
        return "0".into();
    }
    if !PLAIN.contains(&a) {
        return format!("{:.*e}", (DIGITS - 1) as usize, v);
    }
    let decimals = (DIGITS - 1 - a.log10().floor() as i32).clamp(0, 6) as usize;
    let text = format!("{v:.decimals$}");
    // Rounding can make it -0.000.
    if text.trim_start_matches('-').chars().all(|c| c == '0' || c == '.') {
        return "0".into();
    }
    let (sign, digits) = text.strip_prefix('-').map_or(("", text.as_str()), |d| ("-", d));
    let (whole, fraction) = digits.split_once('.').map_or((digits, None), |(w, f)| (w, Some(f)));
    let grouped = grouped(whole);
    match fraction {
        Some(f) => format!("{sign}{grouped}.{f}"),
        None => format!("{sign}{grouped}"),
    }
}

/// `v` read among values as large as `scale` (the largest magnitude beside it): round-off below
/// a billionth of it is 0.
pub fn clean(v: f64, scale: f64) -> f64 {
    if scale.is_finite() && v.abs() <= ROUND_OFF * scale.abs() { 0.0 } else { v }
}

/// [`number`] of `v` [`clean`]ed against `scale`.
pub fn number_in(v: f64, scale: f64) -> String {
    number(clean(v, scale))
}

/// A count with thousands separators: `221,850,309`.
pub fn count(n: u64) -> String {
    grouped(&n.to_string())
}

/// Digits with a comma between each group of three from the right.
fn grouped(digits: &str) -> String {
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// A number with its unit: `0.1048 mm` (the number alone without one).
pub fn with_unit(text: String, unit: Option<&str>) -> String {
    match unit {
        Some(u) => format!("{text} {u}"),
        None => text,
    }
}

/// A vector's components, round-off among them (against the largest, or `scale` when larger)
/// shown as 0: `0, 0, -0.1048`.
pub fn components(v: &[f64], scale: f64) -> String {
    let largest = v.iter().map(|c| c.abs()).filter(|c| c.is_finite()).fold(scale.abs(), f64::max);
    v.iter().map(|c| number_in(*c, largest)).collect::<Vec<_>>().join(", ")
}

#[cfg(test)]
pub(crate) mod tests {
    use serde_json::json;

    use fieldcraft_engine::results::{AUTO_FRACTION, magnitude};
    use fieldcraft_render::contour::band_color;

    use super::*;

    fn plot_at(s: &Session, base: &DisplayMesh, frame: Option<u64>) -> Option<Plot> {
        let (view, result, _) = shown(s)?;
        let bodies: Vec<&str> = s.design.iter().flat_map(|d| &d.bodies).map(|b| b.name.as_str()).collect();
        plot(s, base, frame, &prepare(s, base, pick(view, result, frame), &bodies)?)
    }

    /// The sample cantilever, solved (static step `Load`, frequency step `Modes`).
    pub(crate) fn solved_cantilever() -> Session {
        let mut s = Session::default();
        s.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        s.execute("solve.static", &json!({})).unwrap();
        s.execute("solve.frequency", &json!({})).unwrap();
        s
    }

    /// The label of the node at `xyz`, given in millimetres (the model is in metres).
    pub(crate) fn node_at(s: &Session, xyz: [f64; 3]) -> u32 {
        let at = fieldcraft_engine::units::geometry_point(xyz);
        *s.model.nodes().iter().find(|(_, n)| n.xyz == at).unwrap().0
    }

    /// A point given in millimetres, in the scene's metres.
    pub(crate) fn mm(xyz: [f64; 3]) -> Vec3 {
        Vec3::from(fieldcraft_engine::units::geometry_point(xyz))
    }

    #[test]
    fn the_deformed_shape_moves_by_the_scaled_displacement() {
        let mut s = solved_cantilever();
        let base = DisplayMesh::from_model(&s.model);
        s.execute("results.clear", &json!({})).unwrap();
        assert!(plot_at(&s, &base, None).is_none(), "nothing shown after results.clear");
        s.execute("results.deformed", &json!({"scale": 2.0})).unwrap();
        let p = plot_at(&s, &base, None).unwrap();
        assert_eq!(p.scale, 2.0);
        // A tip corner moves twice its displacement.
        let tip = node_at(&s, [100.0, 0.0, 0.0]);
        let Kind::Static(r) = p.result.kind() else { panic!("a static result") };
        let u = r.displacement[r.nodes.binary_search(&tip).unwrap()];
        let at = p.node_position(tip).unwrap();
        assert!((at.z - 2.0 * u[2]).abs() < 1e-15 && u[2] < -0.2e-3, "{at:?} {u:?}");
        // Automatic: the largest displacement is a tenth of the model's size.
        s.execute("results.deformed", &json!({"scale": "auto"})).unwrap();
        let p = plot_at(&s, &base, None).unwrap();
        let largest = r.displacement.iter().map(|v| magnitude(*v)).fold(0.0, f64::max);
        let auto = AUTO_FRACTION * base.bounds.diagonal() / largest;
        assert!((p.scale / auto - 1.0).abs() < 1e-12, "{}", p.scale);
        // Undeformed, the plot still knows the scale ticking Deformed would draw (#325).
        s.execute("results.deformed", &json!({"show": false})).unwrap();
        let p = plot_at(&s, &base, None).unwrap();
        assert_eq!(p.display.tris, base.tris);
        assert!((p.scale / auto - 1.0).abs() < 1e-12, "{}", p.scale);
    }

    #[test]
    fn contours_colour_known_nodes_with_their_band() {
        let mut s = solved_cantilever();
        s.execute("results.contour", &json!({"field": "displacement_z", "bands": 4})).unwrap();
        s.execute("results.deformed", &json!({"show": false})).unwrap();
        let base = DisplayMesh::from_model(&s.model);
        let p = plot_at(&s, &base, None).unwrap();
        let legend = p.contour.clone().unwrap();
        // The tip moves most (lowest band), the clamped root not at all (highest band).
        assert!(legend.min.1 < -0.2 && legend.max.1.abs() < 1e-12, "{legend:?}");
        let t = Tokens::dark();
        let scene = scene(&p, Style::Shaded, crate::scene::Show::DEFAULT, &t);
        let stops = t.contour_stops();
        let (low, high) = (band_color(&stops, 0, 4), band_color(&stops, 3, 4));
        // Faces whose corners are all near the root are in the top band; near the tip, the bottom.
        // (The scene is in metres: the beam runs from x = 0 to 0.1.)
        let root_end = scene.tris.iter().filter(|tri| tri.p.iter().all(|q| q.x < 0.0005)).map(|tri| tri.color).collect::<Vec<_>>();
        assert!(!root_end.is_empty() && root_end.iter().all(|c| *c == high), "{root_end:?}");
        let tip_end = scene.tris.iter().filter(|tri| tri.p.iter().all(|q| q.x > 0.0995)).map(|tri| tri.color).collect::<Vec<_>>();
        assert!(!tip_end.is_empty() && tip_end.iter().all(|c| *c == low), "{tip_end:?}");
        // Every band appears, and the pieces still cover the faces (more triangles than before).
        for k in 0..4 {
            assert!(scene.tris.iter().any(|tri| tri.color == band_color(&stops, k, 4)), "band {k}");
        }
        assert!(scene.tris.len() > base.tris.len());
        // Hiding the mesh hides the element edges, not the contours (#317).
        let hidden = super::scene(&p, Style::ShadedEdges, crate::scene::Show { mesh: false, ..crate::scene::Show::DEFAULT }, &t);
        assert_eq!(hidden.tris, scene.tris);
        assert!(hidden.lines.iter().all(|l| l.color != rgb(t.mesh_edge)), "no element edges");
    }

    /// A 40 × 30 × 10 mm plate in a tank of 5 mm cells, a liquid flowing past it along x, solved.
    fn plate_in_a_tank(thermal: bool) -> Session {
        use fieldcraft_geometry::sample::{Plate, plate_design};
        // A folder of its own per call: tests run in parallel, and one removing a shared folder
        // took the design from under another's `geometry.open`.
        static CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-flow-{}-{n}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("plate.solidcraft");
        std::fs::write(&path, plate_design(&Plate::default()).unwrap()).unwrap();
        let mut s = Session::default();
        s.execute("geometry.open", &json!({"path": path.to_str().unwrap()})).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        let mut script = json!([
            {"command": "model.material.create", "params": {"name": "Water", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-7}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Tank", "material": "Water", "min": [-20, -20, -20], "max": [60, 50, 30], "spacing": 5, "obstacles": ["Body1"],
                "faces": {"x_min": {"type": "velocity", "velocity": [10, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Tank", "steps": 60, "frames": 2}}},
        ]);
        if thermal {
            script[0]["params"]["thermal"] = json!({"conductivity": 1.0, "specific_heat": 1e8});
            let faces: Vec<_> = s.design.as_ref().unwrap().bodies[0]
                .faces()
                .iter()
                .map(|f| {
                    json!({
                        "face": {"body": "Body1", "face": f.name}, "condition": {"type": "temperature", "value": 60.0}
                    })
                })
                .collect();
            script[1]["params"]["thermal"] =
                json!({"initial_temperature": 20.0, "surfaces": faces, "faces": {"x_min": {"type": "temperature", "value": 20.0}}});
        }
        s.run_script(&script).unwrap();
        if thermal {
            s.execute("model.units", &json!({"system": "mm-t-s-N-MPa"})).unwrap();
        }
        s.execute("solve.fluid", &json!({})).unwrap();
        s
    }

    /// After a fluid solve the part stays drawn and the flow is drawn through it: velocity on a
    /// slice and on the part, with its legend, and streamlines (#193).
    #[test]
    fn a_fluid_step_keeps_the_part_and_draws_the_flow() {
        let mut s = plate_in_a_tank(false);
        let (base, _) = fieldcraft_engine::view::display(&s);
        assert!(!base.surface_tris.is_empty());
        let p = plot_at(&s, &base, None).unwrap();
        assert!(p.is_fluid() && p.flow.is_some());
        assert_eq!(p.display.surface_tris, base.surface_tris, "the part stays");
        let legend = p.contour.clone().unwrap();
        assert!(legend.field == Field::Velocity && legend.max.1 > legend.min.1, "{legend:?}");
        let at = p.extreme_position(legend.max.0).unwrap();
        assert!(at.x > -20.0 && at.x < 60.0, "{at:?}");
        assert_eq!(p.frame().map(|(k, _, n)| (k, n)), Some((2, 3)));
        let t = Tokens::light();
        let scene = scene(&p, Style::ShadedEdges, crate::scene::Show::DEFAULT, &t);
        assert!(scene.tris.len() > base.surface_tris.len(), "the slice is drawn");
        assert!(!scene.lines.is_empty(), "streamlines and the part's edges are drawn");
        // Pressure on the slice; the first frame. This plate_in_a_tank fixture sets no initial
        // velocity, so it starts at rest (unlike the flow-past sample's uniform inlet field).
        s.execute("results.slice", &json!({"field": "pressure", "frame": 0})).unwrap();
        let p = plot_at(&s, &base, None).unwrap();
        assert_eq!((p.contour.as_ref().map(|c| c.field), p.pick.frame), (Some(Field::Pressure), Some(0)));
    }

    /// Animating a fluid step plays its frames: the clock maps across the frames once per cycle,
    /// and results.animate turns it on without deforming the part (the flow start-up transient).
    #[test]
    fn a_fluid_step_animates_its_frames() {
        let mut s = plate_in_a_tank(false);
        let (base, _) = fieldcraft_engine::view::display(&s);
        s.execute("results.animate", &json!({})).unwrap();
        let view = s.results_view.as_ref().unwrap();
        let result = s.latest_result(&view.step).unwrap().0;
        assert!(view.animate && !view.deformed, "fluid animates its frames, not a deformation");
        // The animation clock runs even though the fluid step does not deform.
        assert_eq!(frame(view, result, 0.0), Some(0));
        assert_eq!(frame(view, result, PERIOD * 0.5), Some(FRAMES / 2));
        // The clock maps across the three frames: the plate_in_a_tank field, which starts at rest
        // (the fixture sets no initial velocity, unlike the flow-past sample), then the settled flow.
        let n = plot_at(&s, &base, None).unwrap().frame().unwrap().2;
        assert_eq!(n, 3);
        assert_eq!(plot_at(&s, &base, Some(0)).unwrap().pick.frame, Some(0));
        assert_eq!(plot_at(&s, &base, Some(FRAMES - 1)).unwrap().pick.frame, Some(n - 1));
    }

    #[test]
    fn heated_flow_fields_have_units_colours_and_face_markers() {
        let mut s = plate_in_a_tank(true);
        let (base, _) = fieldcraft_engine::view::display(&s);
        s.execute("results.slice", &json!({"field": "temperature"})).unwrap();
        let p = plot_at(&s, &base, None).unwrap();
        assert_eq!(p.contour.as_ref().unwrap().field, Field::Temperature);
        assert_eq!(Field::Temperature.title(&p.result, &p.units), "Temperature [°C]");
        assert!(Field::available(&p.result).any(|f| f == Field::Temperature));
        assert!(p.contour.as_ref().unwrap().max.1 > 20.0);
        let last = p.frame().unwrap().2 - 1;
        s.execute("results.frame", &json!({"frame": 0})).unwrap();
        let at_rest = plot_at(&s, &base, None).unwrap();
        let c = at_rest.contour.as_ref().unwrap();
        let cells = &at_rest.flow.as_ref().unwrap().contour.as_ref().unwrap().values;
        // The cells are coloured from SI (20 °C), the legend in the display units.
        assert!(cells.iter().filter(|t| t.is_finite()).all(|t| *t == 293.15), "the frame selector changes the temperature shown");
        // The faces are held at 60 °C from the start: the legend's maximum is a face, at a point on it.
        assert_eq!(c.max.1, 60.0);
        assert!(at_rest.extreme_position(c.max.0).is_some());
        s.execute("results.frame", &json!({"frame": last})).unwrap();
        s.execute("results.surface", &json!({"field": "wall_heat_flux"})).unwrap();
        let p = plot_at(&s, &base, None).unwrap();
        let c = p.contour.as_ref().unwrap();
        assert_eq!(c.field, Field::WallHeatFlux);
        assert_eq!(c.field.unit(&p.result, &p.units).as_deref(), Some("mW/mm²"));
        assert!(c.min.1 < 0.0, "hot wall gives heat to fluid");
        assert!(p.extreme_position(c.min.0).is_some());
        assert!(p.flow.as_ref().unwrap().slice.is_none());
        for t in [Tokens::light(), Tokens::dark()] {
            let scene = scene(&p, Style::ShadedEdges, crate::scene::Show::DEFAULT, &t);
            assert_eq!(scene.tris.len(), base.surface_tris.len());
            assert!(scene.tris.windows(2).any(|t| t[0].color != t[1].color));
        }
        s.execute("results.frame", &json!({"frame": 0})).unwrap();
        assert!(plot_at(&s, &base, None).unwrap().contour.is_none());
    }

    #[test]
    fn modes_swing_and_static_steps_grow_when_animated() {
        let mut s = solved_cantilever();
        s.execute("results.modes", &json!({"mode": 1})).unwrap();
        s.execute("results.animate", &json!({})).unwrap();
        let base = DisplayMesh::from_model(&s.model);
        let view = s.results_view.clone().unwrap();
        let result = s.latest_result(&view.step).unwrap().0;
        assert_eq!(frame(&view, result, 0.0), Some(0));
        assert_eq!(frame(&view, result, PERIOD * 0.25), Some(FRAMES / 4));
        assert_eq!(frame(&view, result, PERIOD * 1.25), Some(FRAMES / 4));
        // Animate on but Deformed off, and not a fluid step: nothing to animate, so it stays still.
        let mut still = view.clone();
        still.deformed = false;
        assert_eq!(frame(&still, result, PERIOD * 0.25), None);
        assert!((swing(result, Some(FRAMES / 4)) - 1.0).abs() < 1e-12);
        assert!((swing(result, Some(3 * FRAMES / 4)) + 1.0).abs() < 1e-12);
        // At frame 0 a mode is undeformed; a quarter cycle later fully deformed.
        let p0 = plot_at(&s, &base, Some(0)).unwrap();
        let p1 = plot_at(&s, &base, Some(FRAMES / 4)).unwrap();
        assert_eq!(p0.display.tris, base.tris);
        assert_ne!(p1.display.tris, base.tris);
        assert_eq!(p1.mode().map(|(m, _)| m), Some(1));
        s.execute("results.deformed", &json!({"step": "Load"})).unwrap();
        let (_, bend, _) = shown(&s).unwrap();
        assert_eq!(swing(bend, Some(0)), 0.0);
        assert!((swing(bend, Some(FRAMES / 2)) - 1.0).abs() < 1e-12);
        assert_eq!(swing(bend, None), 1.0);
    }

    #[test]
    fn numbers_read_well() {
        assert_eq!(number(0.0), "0");
        assert_eq!(number(1.23456), "1.235");
        assert_eq!(number(-4.024), "-4.024");
        assert_eq!(number(123.456), "123.5");
        assert_eq!(number(1.5e-6), "1.500e-6");
        assert_eq!(number(2.5e7), "2.500e7");
        assert_eq!(number(f64::NAN), "—");
        assert_eq!(number(123_456.7), "123,457");
        assert_eq!(number(-4_321.06), "-4,321");
        assert_eq!(number(-0.000_000_4), "-4.000e-7");
        assert_eq!(count(221_850_309), "221,850,309");
        assert_eq!(count(945), "945");
        // Round-off beside larger values is 0.
        assert_eq!(number_in(-3.16e-11, 1000.0), "0");
        assert_eq!(number_in(4.34e-9, 1e5), "0");
        assert_eq!(number_in(0.5, 1000.0), "0.5000");
        assert_eq!(components(&[9.497e-17, -3.463e-14, -0.2415], 0.0), "0, 0, -0.2415");
        assert_eq!(with_unit(number(0.1048), Some("mm")), "0.1048 mm");
    }
}
