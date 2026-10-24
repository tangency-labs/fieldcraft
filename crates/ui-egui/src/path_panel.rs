//! Path editor, scalar-distance plot and viewport overlay. Numerical sampling is exclusively
//! `results.path` on the existing background-job runner; editor coordinates remain SI.

use egui::{Align2, FontId, Pos2, Rect, RichText, Stroke, pos2, vec2};
use fieldcraft_engine::results::{Field, Kind, Pick, Quantity};
use fieldcraft_engine::{Session, StepResult};
use fieldcraft_render::Vec3;
use serde_json::{Value, json};

use crate::ModelApp;
use crate::results::number;
use crate::theme::Tokens;

const GRAPH_POINTS: usize = 400;
const OVERLAY_POINTS: usize = 256;
const MIN_GRAPH_RELATIVE_RANGE: f64 = 1e-3; // unit-factors: a minimum plotted relative range, not a conversion

#[derive(Debug, Clone)]
pub struct Editor {
    pub step: String,
    pub field: Field,
    pub points: Vec<[f64; 3]>,
    pub samples: usize,
    pub frame: Option<usize>,
    pub mode: Option<usize>,
    pub picking: bool,
    pub selected: Option<usize>,
}

fn fields(result: &StepResult) -> Vec<Field> {
    Field::available(result).filter(|f| *f != Field::ContactStatus).collect()
}

fn stored_frames(result: &StepResult, field: Field) -> usize {
    match result.kind() {
        Kind::Heat(h) if field != Field::HeatFlux => h.frames.len(),
        Kind::Explicit(e) => e.frames.len(),
        Kind::Static(r) if matches!(field, Field::Displacement | Field::DisplacementX | Field::DisplacementY | Field::DisplacementZ) => {
            r.nonlinear.as_ref().map_or(0, |h| h.frames.len())
        }
        _ => 0,
    }
}

impl Editor {
    fn new(step: String, result: &StepResult, contour: Option<Field>, pick: Pick) -> Self {
        let available = fields(result);
        let field = contour.filter(|f| available.contains(f)).or_else(|| available.first().copied()).unwrap_or(Field::Displacement);
        Self {
            step,
            field,
            points: Vec::new(),
            samples: 101,
            frame: pick.frame.filter(|k| *k < stored_frames(result, field)),
            mode: match result.kind() {
                Kind::Frequency(f) => Some(if (1..=f.modes.len()).contains(&pick.mode) { pick.mode } else { 1 }),
                _ => None,
            },
            picking: false,
            selected: None,
        }
    }

    pub fn params(&self, session: &Session) -> Result<Value, String> {
        let points = self
            .points
            .iter()
            .map(|p| {
                Ok([
                    session.units().express(Quantity::Length, p[0]).map_err(|e| e.to_string())?,
                    session.units().express(Quantity::Length, p[1]).map_err(|e| e.to_string())?,
                    session.units().express(Quantity::Length, p[2]).map_err(|e| e.to_string())?,
                ])
            })
            .collect::<Result<Vec<_>, String>>()?;
        let params = json!({"step":self.step,"field":self.field.id(),"points":points,"samples":self.samples,"frame":self.frame,"mode":self.mode});
        // File pickers/jobs may outlive the display setting under which numbers were typed.
        Ok(session.params_with_units("results.path", &params))
    }
}

impl ModelApp {
    pub fn open_path_plot(&mut self) {
        if self.path_editor.is_some() {
            self.follow_stage("results.path");
            return;
        }
        self.restore_path_plot();
    }

    /// Called only after a successful sample/export, or when opening without an existing
    /// draft. Repeated ribbon activation must not discard unsampled edits.
    pub(crate) fn restore_path_plot(&mut self) {
        let shown = self.session.results_view.as_ref();
        let stored = self.session.path_plot.as_ref().filter(|p| p.current(&self.session));
        let step = stored
            .map(|p| p.step.clone())
            .or_else(|| shown.map(|v| v.step.clone()))
            .or_else(|| self.session.solved_steps().next().map(str::to_owned));
        let Some(step) = step else { return };
        let Some(result) = self.session.result(&step).cloned() else {
            self.set_status("Solve the step again before plotting a path", true);
            return;
        };
        if matches!(result.kind(), Kind::Fluid(_)) {
            self.set_status("Path plots currently use solid nodal fields", true);
            return;
        }
        let visible = shown.filter(|v| v.step == step);
        let pick =
            visible.map(|v| crate::results::pick(v, &result, crate::results::frame(v, &result, self.now))).unwrap_or(Pick { mode: 1, frame: None });
        let mut editor = Editor::new(step.clone(), &result, visible.and_then(|v| v.contour), pick);
        if let Some(plot) = stored {
            editor.step = plot.step.clone();
            editor.field = plot.field;
            editor.points = plot.points.clone();
            editor.samples = plot.xyz.len();
            editor.frame = plot.pick.frame;
            editor.mode = plot.mode;
        }
        if self.session.results_view.is_none() {
            let _ = self.run("results.contour", json!({"step":step,"field":editor.field.id()}));
        }
        self.path_editor = Some(editor);
        self.follow_stage("results.path");
    }

    pub fn sample_path(&mut self) -> Result<u64, String> {
        let editor = self.path_editor.as_ref().ok_or("open Path Plot first")?;
        let params = editor.params(&self.session)?;
        self.start_job("results.path", params)
    }

    /// Consume the same viewport click used by the actual pointer handler, without changing
    /// node/face selection or probing a nearest node. A miss simply adds no path vertex.
    pub fn path_pick_at(&mut self, at: Pos2) -> bool {
        if self.dialog.is_some() || !self.path_editor.as_ref().is_some_and(|p| p.picking) {
            return false;
        }
        if self.busy().is_some() {
            return true;
        }
        let Some(rect) = self.viewport.rect.filter(|r| r.contains(at)) else { return false };
        let Some(plot) = self.plot() else { return true };
        let Some(editor) = self.path_editor.as_ref() else { return true };
        if editor.step != plot.view.step || plot.stale {
            return true;
        }
        let (origin, dir) =
            self.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
        let Some(hit) = plot.display.hit_element_face(origin, dir) else { return true };
        let Some(triangle) = plot.display.tris.get(hit.index) else { return true };
        let Some(weights) = barycentric(triangle.p, origin + dir * hit.t) else { return true };
        let Some(point) = fieldcraft_engine::path::reference_surface_point(&self.session.model, triangle.element, triangle.nodes, weights) else {
            return true;
        };
        if let Some(editor) = self.path_editor.as_mut() {
            if editor.points.len() >= fieldcraft_engine::path::MAX_POINTS {
                self.set_status("Path already has 1000 points", true);
                return true;
            }
            editor.points.push(point);
        }
        true
    }
}

/// Scale-relative dominant-plane barycentrics of a genuine ray-triangle hit. Computing on
/// translated, normalized edge vectors avoids squaring world coordinates or a dot-product
/// Gram determinant near cancellation. No closest-node approximation is involved.
fn barycentric(triangle: [Vec3; 3], point: Vec3) -> Option<[f64; 3]> {
    let edge = |p: Vec3| {
        let p = p - triangle[0];
        [p.x, p.y, p.z]
    };
    let (a, b, p) = (edge(triangle[1]), edge(triangle[2]), edge(point));
    let scale = a.iter().chain(&b).map(|v| v.abs()).fold(0.0, f64::max);
    if !(scale > 0.0 && scale.is_finite()) {
        return None;
    }
    let (a, b, p) = (a.map(|v| v / scale), b.map(|v| v / scale), p.map(|v| v / scale));
    let cross = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let k = (0..3).max_by(|a, b| cross[*a].abs().total_cmp(&cross[*b].abs()))?;
    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
    let denominator = cross[k];
    if !(denominator.abs() > f64::MIN_POSITIVE && denominator.is_finite()) {
        return None;
    }
    let u = (p[i] * b[j] - p[j] * b[i]) / denominator;
    let v = (a[i] * p[j] - a[j] * p[i]) / denominator;
    let weights = [1.0 - u - v, u, v];
    weights.iter().all(|v| v.is_finite()).then_some(weights)
}

pub fn controls(app: &mut ModelApp, view: &crate::results::Plot, ui: &mut egui::Ui, t: &Tokens) {
    let Some(mut editor) = app.path_editor.take() else { return };
    let busy = app.busy().is_some();
    let mut close = false;
    ui.horizontal(|ui| {
        ui.label(RichText::new("Path Plot").strong());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            close = ui.small_button("×").on_hover_text("Close path plot").clicked();
        });
    });
    if close {
        if app.loading.as_ref().is_some_and(|l| l.command == "results.path") {
            app.cancel_loading();
        }
        return;
    }
    let mut sample = false;
    let mut export = false;
    ui.add_enabled_ui(!busy, |ui| {
        egui::ComboBox::from_id_salt("fc_path_step").selected_text(&editor.step).width(ui.available_width() - 8.0).show_ui(ui, |ui| {
            for step in app.session.solved_steps() {
                if ui.selectable_label(editor.step == step, step).clicked() {
                    editor.step = step.to_owned();
                    editor.frame = None;
                    editor.mode = app.session.result(step).filter(|r| matches!(r.kind(), Kind::Frequency(_))).map(|_| 1);
                }
            }
        });
        let Some(result) = app.session.result(&editor.step) else {
            ui.label(RichText::new("Result is stale; solve it again").color(t.error));
            return;
        };
        let available = fields(result);
        if !available.contains(&editor.field) {
            editor.field = available.first().copied().unwrap_or(Field::Displacement);
        }
        egui::ComboBox::from_id_salt("fc_path_field").selected_text(editor.field.label()).width(ui.available_width() - 8.0).show_ui(ui, |ui| {
            for field in available {
                ui.selectable_value(&mut editor.field, field, field.label());
            }
        });
        ui.label(RichText::new(format!("Reference coordinates [{}]", app.session.units().symbol(Quantity::Length))).small().color(t.text_dim));
        let mut remove = None;
        egui::ScrollArea::vertical().id_salt("fc_path_rows").max_height(130.0).show(ui, |ui| {
            for (i, point) in editor.points.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 3.0;
                    ui.label(format!("{}", i + 1));
                    // Account for the remove button and all three remaining widget gaps.
                    let width = ((ui.available_width() - 20.0 - 3.0 * ui.spacing().item_spacing.x) / 3.0).max(30.0);
                    for (k, axis) in ["x", "y", "z"].into_iter().enumerate() {
                        let mut value = app.session.units().express(Quantity::Length, point[k]).unwrap_or(f64::NAN);
                        if ui.add_sized([width, 20.0], egui::DragValue::new(&mut value).speed(0.1).prefix(format!("{axis} "))).changed()
                            && let Ok(stored) = app.session.units().to_si(Quantity::Length, value)
                        {
                            point[k] = stored;
                        }
                    }
                    if ui.small_button("×").on_hover_text("Remove point").clicked() {
                        remove = Some(i);
                    }
                });
            }
        });
        if let Some(i) = remove {
            editor.points.remove(i);
        }
        ui.horizontal(|ui| {
            if ui.add_enabled(editor.points.len() < fieldcraft_engine::path::MAX_POINTS, egui::Button::new("Add point")).clicked() {
                editor.points.push(editor.points.last().copied().unwrap_or([0.0; 3]));
            }
            ui.checkbox(&mut editor.picking, "Pick surface points")
                .on_hover_text("Click the visible mesh to add reference-coordinate points, including on deformed and curved element faces.");
        });
        if editor.picking && editor.step != view.view.step {
            ui.label(RichText::new("Show this step before picking its surface").small().color(t.text_dim));
        }
        ui.horizontal(|ui| {
            ui.label("Samples");
            ui.add(egui::DragValue::new(&mut editor.samples).range(2..=fieldcraft_engine::path::MAX_SAMPLES));
        });
        match result.kind() {
            Kind::Frequency(f) => {
                let mode = editor.mode.get_or_insert(1);
                ui.horizontal(|ui| {
                    ui.label("Mode");
                    ui.add(egui::DragValue::new(mode).range(1..=f.modes.len().max(1)));
                });
            }
            _ => {
                let frames = stored_frames(result, editor.field);
                if frames > 0 {
                    let mut chosen = editor.frame.is_some();
                    ui.checkbox(&mut chosen, "Stored frame (otherwise final)");
                    if chosen {
                        let frame = editor.frame.get_or_insert(frames.saturating_sub(1));
                        ui.add(egui::DragValue::new(frame).range(0..=frames.saturating_sub(1)));
                    } else {
                        editor.frame = None;
                    }
                } else {
                    editor.frame = None;
                    ui.label(RichText::new("Final state").small().color(t.text_dim));
                }
            }
        }
        ui.horizontal(|ui| {
            sample = ui.add_enabled(editor.points.len() >= 2, egui::Button::new("Sample path")).clicked();
            export = ui.add_enabled(editor.points.len() >= 2, egui::Button::new("Export CSV…")).clicked();
        });
    });
    if let Some(plot) = app.session.path_plot.as_ref() {
        ui.separator();
        let frame = plot.pick.frame.map_or_else(|| "final state".to_string(), |k| format!("frame {k}"));
        let mode = plot.mode.map(|mode| format!("mode {mode}"));
        let result = plot.result.upgrade();
        let title = result.as_ref().map_or_else(|| plot.field.label().to_owned(), |r| plot.field.title(r, app.session.units()));
        ui.label(RichText::new(format!("{} · {}", plot.step, mode.as_deref().unwrap_or(&frame))).strong());
        ui.label(RichText::new(title).small().color(t.text_dim));
        if !plot.current(&app.session) {
            ui.label(RichText::new("Plotted result is stale; sample again after solving").color(t.error));
        } else if let Some(result) = result {
            let distance: Vec<_> =
                plot.distance.iter().map(|v| fieldcraft_engine::results::shown(Quantity::Length, *v, app.session.units())).collect();
            let values: Vec<_> = plot.values.iter().map(|v| v.map(|v| plot.field.shown(&result, v, app.session.units()))).collect();
            if let Some(selected) = graph(ui, &distance, &values, app.session.units().symbol(Quantity::Length), t) {
                editor.selected = Some(selected);
            }
            let valid = plot.values.iter().filter(|v| v.is_some()).count();
            ui.label(RichText::new(format!("{valid}/{} values · {} gaps", plot.values.len(), plot.values.len() - valid)).small().color(t.text_dim));
        }
    }
    app.path_editor = Some(editor);
    if sample && let Err(e) = app.sample_path() {
        app.set_status(e, true);
    }
    if export {
        let params = app.path_editor.as_ref().and_then(|p| p.params(&app.session).ok());
        if let Some(params) = params {
            app.export_csv_with("results.path", params);
        }
    }
}

/// Keep endpoints and both extrema in their original station order. Any bin with an
/// undefined contributor is a break, so decimation cannot hide a peak or bridge a gap.
fn graph_stations(values: &[Option<f64>], bins: usize) -> Vec<Option<usize>> {
    let stride = values.len().div_ceil(bins.max(1)).max(1);
    let mut selected = Vec::new();
    for start in (0..values.len()).step_by(stride) {
        let end = (start + stride).min(values.len());
        let bin = values.get(start..end).unwrap_or_default();
        if bin.iter().any(|v| v.is_none_or(|v| !v.is_finite())) {
            selected.push(None);
            continue;
        }
        let mut lo = start;
        let mut hi = start;
        for i in start..end {
            if values.get(i).copied().flatten() < values.get(lo).copied().flatten() {
                lo = i;
            }
            if values.get(i).copied().flatten() > values.get(hi).copied().flatten() {
                hi = i;
            }
        }
        let mut candidates = [start, lo, hi, end.saturating_sub(1)];
        candidates.sort_unstable();
        let mut last = None;
        for i in candidates {
            if last != Some(i) {
                selected.push(Some(i));
                last = Some(i);
            }
        }
    }
    selected
}

/// Normalize by the actual field magnitude, so changing its display unit preserves shape.
/// Subtract only normalized bounds, including opposite finite extremes without overflow.
fn graph_range(lo: f64, hi: f64) -> (f64, [f64; 2]) {
    let magnitude = lo.abs().max(hi.abs());
    let magnitude = if magnitude > 0.0 { magnitude } else { 1.0 };
    let (low, high) = (lo / magnitude, hi / magnitude);
    // A nonconstant field uses its actual span. A floor tied to its absolute level
    // would flatten a temperature gradient differently in kelvin and Celsius.
    let span = if high > low { high - low } else { MIN_GRAPH_RELATIVE_RANGE };
    let pad = span * 0.05;
    (magnitude, [low - pad, high + pad])
}

/// Axis limits and hover use all stored stations; drawing is bounded to at most
/// 400 bins/four vertices per fully defined bin, retaining its extrema.
fn graph(ui: &mut egui::Ui, distance: &[f64], values: &[Option<f64>], unit: &str, t: &Tokens) -> Option<usize> {
    let finite: Vec<_> = values.iter().flatten().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        ui.label("No defined values along this path");
        return None;
    }
    let (outer, response) = ui.allocate_exact_size(vec2(ui.available_width(), 145.0), egui::Sense::click());
    let rect = Rect::from_min_max(outer.min + vec2(47.0, 4.0), outer.max - vec2(0.0, 19.0));
    let painter = ui.painter_at(outer);
    painter.rect_filled(rect, 2.0, t.field);
    let lo = finite.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (magnitude, [low, high]) = graph_range(lo, hi);
    let total = distance.last().copied().unwrap_or(1.0).max(f64::MIN_POSITIVE);
    let at = |x: f64, y: f64| {
        pos2(rect.left() + (x / total) as f32 * rect.width(), rect.bottom() - ((y / magnitude - low) / (high - low)) as f32 * rect.height())
    };
    for (y, label) in [(lo, Align2::RIGHT_BOTTOM), (hi, Align2::RIGHT_TOP)] {
        painter.text(pos2(rect.left() - 4.0, at(0.0, y).y), label, number(y), FontId::proportional(10.0), t.text_dim);
    }
    painter.text(pos2(rect.left(), rect.bottom() + 2.0), Align2::LEFT_TOP, "0", FontId::proportional(10.0), t.text_dim);
    painter.text(
        pos2(rect.right(), rect.bottom() + 2.0),
        Align2::RIGHT_TOP,
        format!("s {} {unit}", number(total)),
        FontId::proportional(10.0),
        t.text_dim,
    );
    let mut line = Vec::new();
    let flush = |line: &mut Vec<Pos2>| {
        if let [point] = line.as_slice() {
            painter.circle_filled(*point, 2.0, t.accent);
        }
        if !line.is_empty() {
            painter.add(egui::Shape::line(std::mem::take(line), Stroke::new(1.8, t.accent)));
        }
    };
    for index in graph_stations(values, GRAPH_POINTS) {
        let Some(i) = index else {
            flush(&mut line);
            continue;
        };
        if let (Some(x), Some(Some(y))) = (distance.get(i), values.get(i)) {
            line.push(at(*x, *y));
        }
    }
    flush(&mut line);
    let hover = response.hover_pos().filter(|p| rect.contains(*p))?;
    let ratio = f64::from(((hover.x - rect.left()) / rect.width()).clamp(0.0, 1.0));
    let index = (ratio * distance.len().saturating_sub(1) as f64).round() as usize;
    let x = distance.get(index).copied()?;
    let value = values.get(index).copied().flatten();
    painter.line_segment([pos2(at(x, 0.0).x, rect.top()), pos2(at(x, 0.0).x, rect.bottom())], Stroke::new(1.0, t.text_dim));
    response.on_hover_text(format!("s {} {unit}\n{}", number(x), value.map_or_else(|| "Gap / undefined".to_string(), number)));
    Some(index)
}

pub fn overlay(app: &mut ModelApp, painter: &egui::Painter, rect: Rect, t: &Tokens) {
    let Some(editor) = app.path_editor.as_ref() else { return };
    let selected = editor.selected;
    let Some(path) = app.session.path_plot.clone().filter(|p| p.current(&app.session)) else { return };
    let Some(view) = app.plot().filter(|v| v.view.step == path.step) else { return };
    let frame = crate::results::frame(&view.view, &view.result, app.now);
    let factor = if view.view.deformed {
        view.scale * fieldcraft_engine::results::deformation_unit(&view.result, &view.units) * crate::results::swing(&view.result, frame)
    } else {
        0.0
    };
    let stride = path.xyz.len().div_ceil(OVERLAY_POINTS).max(1);
    let mut previous = None;
    for start in (0..path.xyz.len()).step_by(stride) {
        let end = (start + stride).min(path.xyz.len());
        if path.elements.get(start..end).unwrap_or_default().iter().any(Option::is_none) {
            previous = None;
            continue;
        }
        for i in [start, end.saturating_sub(1)] {
            let at = path
                .displayed_point(&app.session, i, &view.view.step, view.pick, factor)
                .and_then(|p| crate::viewport::to_screen(app, Vec3::from(p)));
            if let Some(at) = at.filter(|p| rect.contains(*p)) {
                if let Some(before) = previous {
                    painter.line_segment([before, at], Stroke::new(2.0, t.accent));
                }
                previous = Some(at);
                if i.is_multiple_of(stride.saturating_mul(8).max(1)) {
                    painter.circle(at, 2.8, t.accent, Stroke::new(1.0, t.overlay_text));
                }
            } else {
                previous = None;
            }
        }
    }
    if let Some(i) = selected
        && let Some(at) = path
            .displayed_point(&app.session, i, &view.view.step, view.pick, factor)
            .and_then(|p| crate::viewport::to_screen(app, Vec3::from(p)))
            .filter(|p| rect.contains(*p))
    {
        painter.circle(at, 5.0, t.warm, Stroke::new(1.5, t.overlay_text));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_model::units::{DisplaySystem, Preset};
    use std::sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    };

    fn app() -> ModelApp {
        let mut session = Session::default();
        session.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        session.execute("solve.static", &json!({"step":"Load"})).unwrap();
        ModelApp::new(session)
    }

    fn draw(app: &mut ModelApp, ctx: &egui::Context) {
        let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, vec2(1280.0, 960.0))), ..Default::default() };
        let mut output = ctx.run_ui(raw, |ui| {
            app.logic(ui.ctx());
            app.ui(ui);
        });
        output.textures_delta.clear();
    }

    fn query(editor: &mut Editor) {
        editor.points = vec![[0.02, 0.005, 0.02], [0.08, 0.005, 0.02]];
        editor.samples = 17;
        editor.field = Field::StressXx;
    }

    #[test]
    fn typed_editor_coordinates_remain_si_across_all_display_presets() {
        let mut session = Session::default();
        for preset in Preset::ALL {
            session.default_units = DisplaySystem::new(preset);
            let editor = Editor {
                step: "Load".into(),
                field: Field::StressXx,
                points: vec![[0.02, 0.005, 0.02], [0.08, 0.005, 0.02]],
                samples: 17,
                frame: None,
                mode: None,
                picking: false,
                selected: None,
            };
            let params = editor.params(&session).unwrap();
            let form = fieldcraft_engine::find_command("results.path").unwrap().form.unwrap();
            let converted = fieldcraft_engine::units::convert_params(form, &params, session.units(), fieldcraft_engine::units::si()).unwrap();
            let points: Vec<[f64; 3]> = serde_json::from_value(converted["points"].clone()).unwrap();
            assert!(points.iter().zip(&editor.points).all(|(a, b)| a.iter().zip(b).all(|(a, b)| (a - b).abs() < 1e-15)));
        }
    }

    #[test]
    fn opening_sampling_and_cancelling_use_the_command_and_background_job() {
        let mut app = app();
        let ctx = egui::Context::default();
        draw(&mut app, &ctx);
        app.run("results.contour", json!({"field":"stress_xx"})).unwrap();
        app.start("results.path").unwrap();
        assert!(app.dialog.is_none());
        assert_eq!(app.path_editor.as_ref().unwrap().field, Field::StressXx);
        query(app.path_editor.as_mut().unwrap());
        app.start("results.path").unwrap();
        assert_eq!(app.path_editor.as_ref().unwrap().points, vec![[0.02, 0.005, 0.02], [0.08, 0.005, 0.02]]);
        assert_eq!(app.path_editor.as_ref().unwrap().samples, 17);
        let journal = app.session.journal.len();
        let model = app.session.model.clone();
        app.sample_path().unwrap();
        assert_eq!(app.loading.as_ref().unwrap().command, "results.path");
        app.wait_loaded();
        assert_eq!(app.session.path_plot.as_ref().unwrap().values.len(), 17);
        assert_eq!(app.session.model, model);
        assert_eq!(app.session.journal.len(), journal);
        let previous = app.session.path_plot.clone().unwrap();
        app.path_editor.as_mut().unwrap().samples = 10_000;
        app.sample_path().unwrap();
        app.cancel_loading();
        assert!(Arc::ptr_eq(app.session.path_plot.as_ref().unwrap(), &previous));
        app.run("solve.frequency", json!({"step":"Modes"})).unwrap();
        app.run("results.modes", json!({"step":"Modes","mode":2})).unwrap();
        app.open_path_plot();
        let editor = app.path_editor.as_ref().unwrap();
        assert_eq!(editor.step, "Load");
        assert_eq!(editor.field, Field::StressXx);
        assert!(editor.mode.is_none());
        assert_eq!(app.session.results_view.as_ref().unwrap().step, "Modes");
        draw(&mut app, &ctx);
    }

    #[test]
    fn fresh_draft_inherits_visible_mode_and_compatible_temperature_frame() {
        let mut app = app();
        app.run("solve.frequency", json!({"step":"Modes"})).unwrap();
        app.run("results.modes", json!({"step":"Modes","mode":2})).unwrap();
        app.open_path_plot();
        assert_eq!(app.path_editor.as_ref().unwrap().mode, Some(2));
        assert!(app.path_editor.as_ref().unwrap().frame.is_none());

        let mut script: Value = serde_json::from_str(include_str!("../../../docs/recipes/path-temperature.json")).unwrap();
        let commands = script.as_array_mut().unwrap();
        commands.retain(|c| c["command"] != "results.path");
        for command in commands {
            match command["command"].as_str().unwrap() {
                "model.material.create" => {
                    command["params"]["density"] = json!("7800 kg/m^3");
                    command["params"]["thermal"]["specific_heat"] = json!("500 J/(kg*K)");
                }
                "model.step.create" => {
                    command["params"]["procedure"]["transient"] = json!({"end_time":2,"increment":1,"initial_temperature":"20 °C"})
                }
                _ => {}
            }
        }
        let mut session = Session::default();
        session.run_script(&script).unwrap();
        // A real stored heat result with the visible view on its first frame. The baseline
        // results.frame command omits heat; setting view state here isolates draft inheritance.
        session.results_view.as_mut().unwrap().frame = Some(0);
        let mut app = ModelApp::new(session);
        app.open_path_plot();
        assert_eq!(app.path_editor.as_ref().unwrap().field, Field::Temperature);
        assert_eq!(app.path_editor.as_ref().unwrap().frame, Some(0));
        app.path_editor = None;
        app.session.results_view.as_mut().unwrap().contour = Some(Field::HeatFlux);
        app.open_path_plot();
        assert_eq!(app.path_editor.as_ref().unwrap().field, Field::HeatFlux);
        assert!(app.path_editor.as_ref().unwrap().frame.is_none(), "heat flux is final-only");
    }

    #[test]
    fn programmatic_path_never_opens_a_picker_and_ui_export_runs_as_a_job() {
        let mut app = app();
        let called = Arc::new(AtomicBool::new(false));
        let observed = called.clone();
        app.services.pick_save = Some(Box::new(move |_, _| {
            observed.store(true, Ordering::Relaxed);
            None
        }));
        let params = json!({"step":"Load","field":"stress_xx","points":[[20,5,20],[80,5,20]],"samples":3});
        app.run("results.path", params.clone()).unwrap();
        assert!(!called.load(Ordering::Relaxed));
        let file = std::env::temp_dir().join(format!("fieldcraft-path-ui-export-{}.csv", std::process::id()));
        let name = file.to_string_lossy().into_owned();
        app.services.pick_save = Some(Box::new(move |_, _| Some(name.clone())));
        app.export_csv_with("results.path", params);
        assert_eq!(app.loading.as_ref().unwrap().command, "results.path");
        app.wait_loaded();
        let csv = std::fs::read_to_string(&file).unwrap();
        assert!(csv.starts_with("distance,x,y,z,value\r\n"));
        assert_eq!(csv.lines().count(), 4);
        std::fs::remove_file(file).unwrap();
    }

    #[test]
    fn genuine_deformed_viewport_hits_add_reference_points_instead_of_nearest_nodes() {
        let mut app = app();
        let ctx = egui::Context::default();
        draw(&mut app, &ctx);
        app.run("results.deformed", json!({"show":true,"scale":5})).unwrap();
        app.start("results.path").unwrap();
        draw(&mut app, &ctx);
        app.path_editor.as_mut().unwrap().picking = true;
        let plot = app.plot().unwrap();
        let rect = app.viewport.rect.unwrap();
        app.open_form("model.set.create", json!({}));
        assert!(!app.path_pick_at(rect.center()), "an open command form keeps its target picking");
        assert!(app.path_editor.as_ref().unwrap().points.is_empty());
        app.dialog = None;
        let mut added = 0;
        for (index, triangle) in plot.display.tris.iter().enumerate().take(500) {
            let middle = (triangle.p[0] + triangle.p[1] + triangle.p[2]) * (1.0 / 3.0);
            let Some(at) = crate::viewport::to_screen(&app, middle).filter(|p| rect.contains(*p)) else { continue };
            let (origin, dir) =
                app.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
            if plot.display.hit_element_face(origin, dir).is_none_or(|h| h.index != index) {
                continue;
            }
            assert!(app.path_pick_at(at));
            if app.path_editor.as_ref().unwrap().points.len() > added {
                added += 1;
            }
            if added == 2 {
                break;
            }
        }
        assert_eq!(added, 2);
        let editor = app.path_editor.as_ref().unwrap();
        for point in &editor.points {
            assert!(
                app.session.model.nodes().values().all(|n| n.xyz.iter().zip(point).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt() > 1e-6),
                "pick was snapped to a node"
            );
        }
        app.sample_path().unwrap();
        app.wait_loaded();
        assert!(app.session.path_plot.as_ref().unwrap().values.iter().all(Option::is_some));
    }

    #[test]
    fn dominant_plane_barycentrics_are_scaled_and_refuse_degenerate_input() {
        for scale in [1e-12, 1.0, 1e12] {
            for triangle in [
                [Vec3::ZERO, Vec3::new(scale, 0.0, 0.0), Vec3::new(0.0, scale, 0.0)],
                [Vec3::ZERO, Vec3::new(0.0, 0.0, scale), Vec3::new(0.0, scale, 0.0)],
                [Vec3::ZERO, Vec3::new(0.0, scale, 0.0), Vec3::new(0.0, 0.0, scale)],
            ] {
                let point = triangle[1] * 0.3 + triangle[2] * 0.5;
                let weights = barycentric(triangle, point).unwrap();
                assert!(weights.iter().zip([0.2, 0.3, 0.5]).all(|(a, b)| (a - b).abs() < 1e-12));
            }
        }
        assert!(barycentric([Vec3::ZERO; 3], Vec3::ZERO).is_none());
        assert!(barycentric([Vec3::ZERO, Vec3::new(f64::INFINITY, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)], Vec3::ZERO).is_none());
    }

    #[test]
    fn graph_decimation_retains_extrema_in_station_order_and_breaks_undefined_bins() {
        let values = [0.0, 100.0, 0.0, 0.0, -100.0, 0.0].map(Some);
        assert_eq!(graph_stations(&values, 2), vec![Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)]);
        for values in [[0.0, -100.0, 100.0, 0.0], [0.0, 100.0, -100.0, 0.0]] {
            assert_eq!(graph_stations(&values.map(Some), 1), vec![Some(0), Some(1), Some(2), Some(3)]);
        }
        assert_eq!(graph_stations(&[Some(0.0), None, Some(100.0), Some(4.0), Some(5.0), Some(4.0)], 2), vec![None, Some(3), Some(4), Some(5)]);
        assert_eq!(graph_stations(&[Some(f64::NAN)], 1), vec![None]);
        let full = vec![Some(1.0); 10_000];
        assert!(graph_stations(&full, GRAPH_POINTS).len() <= 4 * GRAPH_POINTS);
    }

    #[test]
    fn graph_scaling_is_unit_invariant_for_small_values_and_finite_at_extremes() {
        let position = |value: f64, lo: f64, hi: f64| {
            let (magnitude, [low, high]) = graph_range(lo, hi);
            (value / magnitude - low) / (high - low)
        };
        for scale in [1e-12, 1.0, 1e12] {
            for value in [1e-9, 1.5e-9, 2e-9] {
                assert!((position(value * scale, 1e-9 * scale, 2e-9 * scale) - position(value, 1e-9, 2e-9)).abs() < 1e-14);
            }
        }
        assert!((position(0.0, 0.0, 0.0) - 0.5).abs() < 1e-14);
        for value in [-f64::MAX, 0.0, f64::MAX] {
            let at = position(value, -f64::MAX, f64::MAX);
            assert!(at.is_finite() && (0.0..=1.0).contains(&at));
        }
        assert!((position(1e-300, 1e-300, 2e-300) - position(1.0, 1.0, 2.0)).abs() < 1e-14);
        for value in [20.0, 20.005, 20.01] {
            let expected = position(value, 20.0, 20.01);
            for (scale, offset) in [(1.0, 273.15), (1.8, 32.0)] {
                assert!((position(value * scale + offset, 20.0 * scale + offset, 20.01 * scale + offset) - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn coordinate_rows_fit_the_context_width_without_rounding_stored_points() {
        let mut app = app();
        let ctx = egui::Context::default();
        draw(&mut app, &ctx);
        app.run("results.path", json!({"step":"Load","field":"stress_xx","points":[[20,5,20],[80,5,20]],"samples":17})).unwrap();
        let original = app.path_editor.as_ref().unwrap().points.clone();
        for preset in Preset::ALL {
            app.session.execute("model.units", &json!({"system":preset.label()})).unwrap();
            let plot = app.plot().unwrap();
            let tokens = Tokens::get(app.ui.dark);
            let mut sizes = (0.0, 0.0);
            let mut output =
                ctx.run_ui(egui::RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, vec2(640.0, 1000.0))), ..Default::default() }, |ui| {
                    ui.scope_builder(egui::UiBuilder::new().max_rect(Rect::from_min_size(Pos2::ZERO, vec2(294.0, 960.0))), |ui| {
                        let width = ui.available_width();
                        controls(&mut app, &plot, ui, &tokens);
                        sizes = (ui.min_rect().width(), width);
                    });
                });
            output.textures_delta.clear();
            assert!(sizes.0 <= sizes.1 + 0.5, "{preset:?}: content {} exceeds {}", sizes.0, sizes.1);
            assert_eq!(app.path_editor.as_ref().unwrap().points, original);
        }
    }

    #[test]
    fn curved_deformed_surface_picks_round_trip_without_pretty_rounding_in_every_preset() {
        let mut session = Session::default();
        let script: Value = serde_json::from_str(include_str!("../../../docs/recipes/path-curved-picking.json")).unwrap();
        session.run_script(&script).unwrap();
        let mut app = ModelApp::new(session);
        let ctx = egui::Context::default();
        draw(&mut app, &ctx);
        app.run("results.deformed", json!({"show":true,"scale":5})).unwrap();
        app.start("results.path").unwrap();
        draw(&mut app, &ctx);
        app.path_editor.as_mut().unwrap().picking = true;
        let plot = app.plot().unwrap();
        let rect = app.viewport.rect.unwrap();
        for (index, triangle) in plot.display.tris.iter().enumerate() {
            let middle = (triangle.p[0] + triangle.p[1] + triangle.p[2]) * (1.0 / 3.0);
            let Some(at) = crate::viewport::to_screen(&app, middle).filter(|p| rect.contains(*p)) else { continue };
            let (origin, dir) =
                app.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
            if plot.display.hit_element_face(origin, dir).is_none_or(|h| h.index != index) {
                continue;
            }
            assert!(app.path_pick_at(at));
            if app.path_editor.as_ref().unwrap().points.len() == 2 {
                break;
            }
        }
        let original = app.path_editor.as_ref().unwrap().points.clone();
        assert_eq!(original.len(), 2);
        app.path_editor.as_mut().unwrap().samples = 2;
        // Two stations deliberately test the picked curved-face points themselves. A straight
        // chord between surface picks can leave a concave solid and must still report gaps.
        for preset in Preset::ALL {
            app.session.execute("model.units", &json!({"system":preset.label()})).unwrap();
            draw(&mut app, &ctx);
            let params = app.path_editor.as_ref().unwrap().params(&app.session).unwrap();
            app.run("results.path", params).unwrap();
            assert!(app.session.path_plot.as_ref().unwrap().values.iter().all(Option::is_some), "{preset:?}");
            let restored = &app.path_editor.as_ref().unwrap().points;
            assert!(original.iter().zip(restored).all(|(a, b)| a.iter().zip(b).all(|(a, b)| (a - b).abs() < 1e-15)));
        }
    }
}
