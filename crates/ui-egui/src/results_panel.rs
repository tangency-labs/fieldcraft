//! What the viewport shows over a results plot: the legend (field and unit, bands, range,
//! deformation scale, a stale-result warning), markers at the contour's largest and smallest
//! value, and the probed node; and the results controls of the context panel (right), which run
//! the `results.*` commands (a fluid step's: its slice, coloured obstacles, streamlines and
//! frames). Numbers go through one formatter ([`crate::results::number`]) with
//! the units of the model's unit system ([`Quantity`]).

use egui::{Align2, Color32, FontId, Pos2, Rect, RichText, Stroke, pos2, vec2};
use fieldcraft_engine::results::{Field, Kind, Quantity, true_scale, with_unit as titled};
use fieldcraft_explicit::{Energies, ExplicitResult};
use serde_json::{Value, json};

use std::sync::Arc;

use crate::ModelApp;
use crate::results::{Plot, components, count, number, number_in, with_unit};
use crate::theme::Tokens;

/// Width of the panel's drop-down lists.
const COMBO_WIDTH: f32 = 160.0;
const SWATCH: egui::Vec2 = vec2(18.0, 15.0);
/// Height of the energy plot's drawing area (points), the room its axis labels take, and the
/// most points one of its lines draws.
const PLOT_HEIGHT: f32 = 110.0;
const AXIS_LEFT: f32 = 48.0;
const AXIS_BOTTOM: f32 = 16.0;
const MAX_PLOT_POINTS: usize = 400;
/// Gap between a marker and its label, and the label's padding (points).
const LABEL_GAP: f32 = 12.0;
const LABEL_PAD: egui::Vec2 = vec2(8.0, 4.0);

/// Draw the legend, the extremes' markers and the probe when a results plot is shown. A click
/// on the legend's Max or Min line turns the view to that node and probes it.
pub fn show(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect, t: &Tokens) {
    let Some(plot) = app.plot() else { return };
    let painter = ui.painter_at(rect);
    let extremes = legend(&plot, &painter, rect, t);
    if app.ui.result_extremes {
        markers(app, &plot, &painter, rect, t);
    }
    probe(app, &plot, &painter, rect, t);
    let mut focus = None;
    for (i, (line, node)) in extremes.into_iter().enumerate() {
        let tip = if plot.is_fluid() { "Centre this field extreme".into() } else { format!("Show and probe node {node}") };
        let r = ui.interact(line, ui.id().with(("fc_legend_extreme", i)), egui::Sense::click()).on_hover_text(tip);
        if r.hovered() {
            ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
        }
        if r.clicked() {
            focus = Some(node);
        }
    }
    if let Some(node) = focus {
        focus_node(app, &plot, node);
    }
}

/// Turn the view about node `node` of the plot (keeping the zoom) and probe it; a fluid step's
/// extreme is a cell, turned about only.
fn focus_node(app: &mut ModelApp, plot: &Plot, node: u32) {
    if let Some(at) = plot.extreme_position(node).filter(|p| p.is_finite()) {
        app.cam.target = at;
    }
    if plot.is_fluid() {
        return;
    }
    // A failure (a stale result) shows in the status bar.
    let _ = app.probe_node(node);
}

/// What the legend shows and where: its frame, title and text lines.
struct Legend {
    frame: Rect,
    title: String,
    lines: Vec<(String, Color32)>,
}

const LEGEND_LINE: f32 = 16.0;

/// Where the legend of `plot` goes in the viewport `rect` (a fit leaves it free).
pub fn legend_rect(plot: &Plot, painter: &egui::Painter, rect: Rect, t: &Tokens) -> Rect {
    legend_layout(plot, painter, rect, t).frame
}

/// The legend's title: the field with its unit (a mode shape normalised), or the shape drawn.
fn title(plot: &Plot) -> String {
    match (&plot.contour, plot.mode()) {
        (Some(l), _) => l.field.title(&plot.result, plot.units),
        (None, Some(_)) => "Mode shape (normalised)".into(),
        (None, None) if plot.is_fluid() => plot.view.contour.map_or_else(|| "Flow".into(), |f| f.title(&plot.result, plot.units)),
        (None, None) => "Deformed shape".into(),
    }
}

/// The deformation line: its scale, and whether it is automatic or true.
fn deformation(plot: &Plot) -> String {
    if !plot.view.deformed {
        return "Undeformed".into();
    }
    let how = match plot.view.scale {
        _ if plot.scale == 1.0 && true_scale(&plot.result) => " (true scale)",
        None => " (auto)",
        Some(_) => "",
    };
    format!("Deformation × {}{how}", number(plot.scale))
}

/// The contour's Max and Min lines.
fn extreme_lines(plot: &Plot) -> Option<[(String, u32); 2]> {
    let l = plot.contour.as_ref()?;
    let unit = l.field.unit(&plot.result, plot.units);
    let scale = l.max.1.abs().max(l.min.1.abs());
    let at = |node: u32| match plot.extreme_position(node).filter(|_| plot.is_fluid()) {
        Some(p) => format!("at {}", components(&[p.x, p.y, p.z], 0.0)),
        None => format!("at node {node}"),
    };
    let line = |what: &str, (node, v): (u32, f64)| (format!("{what} {} {}", with_unit(number_in(v, scale), unit.as_deref()), at(node)), node);
    Some([line("Max", l.max), line("Min", l.min)])
}

/// The legend's text and frame: top left, as wide as its longest line and as tall as its lines
/// and bands.
fn legend_layout(plot: &Plot, painter: &egui::Painter, rect: Rect, t: &Tokens) -> Legend {
    let mut lines: Vec<(String, Color32)> = Vec::new();
    let title = title(plot);
    let mut step = format!("Step {}", plot.view.step);
    if let Some((mode, value)) = plot.mode() {
        step.push_str(&format!(", mode {mode}: {value}"));
    }
    if let Some((k, at, n)) = plot.frame() {
        let last = n.saturating_sub(1);
        match plot.result.kind() {
            // A nonlinear static step's frames are its increments, at their load factors.
            Kind::Static(_) => step.push_str(&format!(", increment {k} of {last}: load factor {}", number(at))),
            _ => step.push_str(&format!(", frame {k} of {last}: t = {}", with_unit(number(at), Quantity::Time.unit(plot.units).as_deref()))),
        }
    }
    lines.push((step, t.overlay_text));
    if !plot.is_fluid() {
        lines.push((deformation(plot), t.text_dim));
    }
    if let Some(l) = plot.contour.as_ref().filter(|_| plot.view.range_min.is_some() || plot.view.range_max.is_some()) {
        let scale = l.max.1.abs().max(l.min.1.abs());
        lines.push((format!("Legend clamped to {} … {}", number_in(l.bands.min, scale), number_in(l.bands.max, scale)), t.warm));
    }
    if plot.stale {
        lines.push(("Stale: the model changed since this step was solved".into(), t.warm));
    }
    let origin = rect.left_top() + vec2(12.0, 12.0);
    let font = FontId::proportional(12.0);
    let line_h = LEGEND_LINE;
    let bands = plot.contour.as_ref().map_or(0, |l| swatches(&l.bands)) as f32;
    // As wide as its longest line (and the extremes' lines under the swatches).
    let extremes = extreme_lines(plot).into_iter().flatten().map(|(text, _)| text);
    let widest = lines
        .iter()
        .map(|(text, _)| text.clone())
        .chain(extremes)
        .map(|text| painter.layout_no_wrap(text, font.clone(), t.overlay_text).size().x)
        .fold(0.0_f32, f32::max);
    let title_w = painter.layout_no_wrap(title.clone(), FontId::proportional(14.0), t.overlay_text).size().x;
    let width = (widest.max(title_w) + 20.0).max(170.0);
    // Padding and title, the lines, then the swatches with the range lines under them.
    let height = 36.0 + line_h * lines.len() as f32 + if bands > 0.0 { 16.0 + bands * SWATCH.y + 2.0 * line_h } else { 0.0 };
    Legend { frame: Rect::from_min_size(origin, vec2(width, height)), title, lines }
}

/// The legend's swatches: the bands, and the "over" band above a clamped maximum.
fn swatches(b: &fieldcraft_render::contour::Bands) -> u32 {
    b.count + u32::from(b.over)
}

/// The legend, top left under the view buttons. Returns the Max and Min lines' rects with their
/// nodes.
fn legend(plot: &Plot, painter: &egui::Painter, rect: Rect, t: &Tokens) -> Vec<(Rect, u32)> {
    let Legend { frame, title, lines } = legend_layout(plot, painter, rect, t);
    let origin = frame.min;
    let (font, line_h) = (FontId::proportional(12.0), LEGEND_LINE);
    painter.rect_filled(frame, 4.0, t.overlay);
    let mut y = origin.y + 8.0;
    painter.text(pos2(origin.x + 10.0, y), Align2::LEFT_TOP, title, FontId::proportional(14.0), t.overlay_text);
    y += 20.0;
    for (text, color) in &lines {
        painter.text(pos2(origin.x + 10.0, y), Align2::LEFT_TOP, text, font.clone(), *color);
        y += line_h;
    }
    let (Some(l), Some(extremes)) = (&plot.contour, extreme_lines(plot)) else { return Vec::new() };
    y += 6.0;
    let stops = t.contour_stops();
    // Round-off among the levels (a range from -1e-15 to 3) reads as 0.
    let scale = l.max.1.abs().max(l.min.1.abs());
    // High values on top, each band's lower boundary written beside its swatch's bottom edge;
    // over a clamped maximum, the "over" band runs up to the field's true maximum.
    let rows = swatches(&l.bands);
    let over = crate::theme::rgb(t.contour_over);
    for k in (0..rows).rev() {
        let row = (rows - 1 - k) as f32;
        let r = Rect::from_min_size(pos2(origin.x + 10.0, y + row * SWATCH.y), SWATCH);
        let c = fieldcraft_render::contour::bands_color(&stops, over, &l.bands, k);
        painter.rect_filled(r, 0.0, Color32::from_rgb(c.0, c.1, c.2));
        if k + 1 == rows {
            let top = if l.bands.is_over(k) { l.max.1 } else { l.bands.level(k + 1) };
            painter.text(pos2(r.right() + 6.0, r.top()), Align2::LEFT_CENTER, number_in(top, scale), font.clone(), t.overlay_text);
        }
        painter.text(pos2(r.right() + 6.0, r.bottom()), Align2::LEFT_CENTER, number_in(l.bands.level(k), scale), font.clone(), t.overlay_text);
    }
    painter.rect_stroke(
        Rect::from_min_size(pos2(origin.x + 10.0, y), vec2(SWATCH.x, SWATCH.y * rows as f32)),
        0.0,
        Stroke::new(1.0, t.border),
        egui::StrokeKind::Outside,
    );
    y += SWATCH.y * rows as f32 + 10.0;
    let mut out = Vec::new();
    for (i, (text, node)) in extremes.into_iter().enumerate() {
        let r = painter.text(pos2(origin.x + 10.0, y + line_h * i as f32), Align2::LEFT_TOP, text, font.clone(), t.text_dim);
        out.push((r.expand(2.0), node));
    }
    out
}

/// Where a label box of `size` goes beside the screen point `at`, kept inside `bounds`: right of
/// the point, or left of it when it would run off the right edge; held inside top and bottom.
/// When it would cover `avoid` (a label already placed), it goes under that box, or over it
/// when under runs off the bottom.
fn label_box(at: Pos2, size: egui::Vec2, bounds: Rect, avoid: Option<Rect>) -> Rect {
    let right = at + vec2(LABEL_GAP, -size.y / 2.0);
    let x = if right.x + size.x > bounds.right() { at.x - LABEL_GAP - size.x } else { right.x }.max(bounds.left());
    let inside = |y: f32| y.clamp(bounds.top(), (bounds.bottom() - size.y).max(bounds.top()));
    let at_y = |y: f32| Rect::from_min_size(pos2(x, inside(y)), size);
    let beside = at_y(right.y);
    match avoid.filter(|a| a.intersects(beside)) {
        None => beside,
        Some(a) => {
            let under = at_y(a.bottom() + LABEL_PAD.y);
            if under.intersects(a) { at_y(a.top() - LABEL_PAD.y - size.y) } else { under }
        }
    }
}

/// A box of `text` beside the screen point `at` ([`label_box`]); returns where it went.
fn label_beside(painter: &egui::Painter, at: Pos2, text: String, bounds: Rect, avoid: Option<Rect>, t: &Tokens) -> Rect {
    let galley = painter.layout_no_wrap(text, FontId::proportional(12.0), t.overlay_text);
    let box_rect = label_box(at, galley.size() + 2.0 * LABEL_PAD, bounds, avoid);
    painter.rect_filled(box_rect, 4.0, t.overlay);
    painter.galley(box_rect.left_top() + LABEL_PAD, galley, t.overlay_text);
    box_rect
}

/// Markers with their values at the contour's largest and smallest value.
fn markers(app: &ModelApp, plot: &Plot, painter: &egui::Painter, rect: Rect, t: &Tokens) {
    let Some(l) = &plot.contour else { return };
    let unit = l.field.unit(&plot.result, plot.units);
    let scale = l.max.1.abs().max(l.min.1.abs());
    let (high, low) = (t.contour[4], t.contour[0]);
    // The Min label keeps off the Max label when the two cells are close on screen.
    let mut placed = None;
    for (what, (node, v), color) in [("Max", l.max, high), ("Min", l.min, low)] {
        let Some(at) = plot.extreme_position(node).and_then(|w| crate::viewport::to_screen(app, w)) else { continue };
        if !rect.contains(at) {
            continue;
        }
        painter.circle(at, 5.0, color, Stroke::new(2.0, t.overlay_text));
        placed = Some(label_beside(painter, at, format!("{what} {}", with_unit(number_in(v, scale), unit.as_deref())), rect, placed, t));
    }
}

/// `results.contour`'s parameters that change one end of the legend's range (`key` `min` or
/// `max`), keeping the field and the other end.
fn legend_params(plot: &Plot, key: &str, value: Value) -> Value {
    let field = plot.view.contour.map_or("none", Field::id);
    let mut range = json!({"min": plot.view.range_min, "max": plot.view.range_max});
    range[key] = value;
    json!({"step": plot.view.step, "field": field, "range": range})
}

/// The probed node: a marker on it and its values beside it, inside the viewport.
fn probe(app: &ModelApp, plot: &Plot, painter: &egui::Painter, rect: Rect, t: &Tokens) {
    let Some(p) = app.probe.as_ref().filter(|p| p.step == plot.view.step && p.mode == plot.pick.mode && p.frame == plot.pick.frame) else { return };
    let Some(at) = plot.node_position(p.node).and_then(|w| crate::viewport::to_screen(app, w)) else { return };
    painter.circle(at, 5.0, t.accent, Stroke::new(2.0, t.overlay_text));
    label_beside(painter, at, p.lines.join("\n"), rect, None, t);
}

/// The results controls (the context panel's Results section): step, contour field,
/// deformation, modes, animation, what the viewport marks, reaction sums and the probed node.
/// Every control runs a command or a view setting.
pub fn controls(app: &mut ModelApp, plot: &Plot, ui: &mut egui::Ui, t: &Tokens) {
    let mut run: Vec<(&'static str, Value)> = Vec::new();
    let mut export: Option<Value> = None;
    let steps: Vec<(String, bool)> = app.session.result_steps().map(|(n, c)| (n.to_string(), c)).collect();
    egui::Grid::new("fc_results_grid").num_columns(2).spacing(vec2(8.0, 4.0)).show(ui, |ui| {
        ui.label("Step");
        egui::ComboBox::from_id_salt("fc_results_step").selected_text(&plot.view.step).width(COMBO_WIDTH).truncate().show_ui(ui, |ui| {
            for (name, current) in &steps {
                let label = if *current { name.clone() } else { format!("{name} (stale)") };
                if ui.selectable_label(*name == plot.view.step, label).clicked() && *name != plot.view.step {
                    run.push(("results.deformed", json!({"step": name})));
                }
            }
        });
        ui.end_row();
        ui.label("Contours");
        let shown = plot.view.contour.map_or("None", Field::label);
        egui::ComboBox::from_id_salt("fc_results_field").selected_text(shown).width(COMBO_WIDTH).truncate().show_ui(ui, |ui| {
            if ui.selectable_label(plot.view.contour.is_none(), "None").clicked() {
                run.push(("results.contour", json!({"field": "none"})));
            }
            for f in Field::available(&plot.result) {
                if ui.selectable_label(plot.view.contour == Some(f), f.label()).clicked() {
                    run.push(("results.contour", json!({"field": f.id()})));
                }
            }
        });
        ui.end_row();
        // The legend's range: each end follows the field (Auto) or is clamped, e.g. below a
        // singular peak (#191).
        if let Some(l) = &plot.contour {
            for (label, key, set, value) in [("Legend max", "max", plot.view.range_max, l.max.1), ("Legend min", "min", plot.view.range_min, l.min.1)]
            {
                ui.label(label);
                ui.horizontal(|ui| {
                    let mut auto = set.is_none();
                    let tip = format!("Follow the field's {}; untick to clamp the legend", if key == "max" { "maximum" } else { "minimum" });
                    if ui.checkbox(&mut auto, "Auto").on_hover_text(tip).changed() {
                        run.push(("results.contour", legend_params(plot, key, if auto { json!("auto") } else { json!(value) })));
                    }
                    let mut v = set.unwrap_or(value);
                    let speed = (l.max.1 - l.min.1).abs().max(1e-9) * 0.005;
                    let drag = egui::DragValue::new(&mut v)
                        .speed(speed)
                        .custom_formatter(|v, _| number(v))
                        .custom_parser(|s| s.replace(',', "").trim().parse::<f64>().ok());
                    if ui.add_enabled(set.is_some(), drag).changed() {
                        run.push(("results.contour", legend_params(plot, key, json!(v))));
                    }
                });
                ui.end_row();
            }
        }
        if plot.is_fluid() {
            flow_rows(plot, ui, &mut run);
        } else {
            deformation_rows(plot, ui, &mut run);
        }
        ui.label("Show");
        ui.horizontal(|ui| {
            if plot.contour.is_some() {
                ui.checkbox(&mut app.ui.result_extremes, "Max/Min").on_hover_text("Mark the nodes of the largest and smallest value");
            }
            ui.checkbox(&mut app.ui.result_symbols, "Loads & BCs").on_hover_text("The symbols of the step's boundary conditions and loads");
        });
        ui.end_row();
    });
    if let Kind::Explicit(e) = plot.result.kind() {
        ui.separator();
        explicit(ui, e, plot, t, &mut run);
    } else if let Some((k, at, n)) = plot.frame().filter(|_| plot.is_fluid()) {
        // A fluid step: its frames, the slice, the obstacles and the streamlines following.
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Frame");
            let mut frame = k;
            if ui.add(egui::Slider::new(&mut frame, 0..=n.saturating_sub(1))).changed() {
                run.push(("results.frame", json!({"frame": frame})));
            }
            let mut animate = plot.view.animate;
            if ui.checkbox(&mut animate, "Animate").on_hover_text("Play the frames from the start-up transient to the settled flow").changed() {
                run.push(("results.animate", json!({"on": animate})));
            }
        });
        ui.label(RichText::new(format!("t = {}", with_unit(number(at), Quantity::Time.unit(plot.units).as_deref()))).color(t.text_dim));
    } else if let Some((k, factor, n)) = plot.frame().filter(|_| matches!(plot.result.kind(), Kind::Static(_))) {
        // A nonlinear static step: its increments, the deformed shape following the slider.
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Increment");
            let mut frame = k;
            if ui.add(egui::Slider::new(&mut frame, 0..=n.saturating_sub(1))).changed() {
                run.push(("results.frame", json!({"frame": frame})));
            }
        });
        ui.label(RichText::new(format!("load factor {}", number(factor))).color(t.text_dim));
    } else if let Kind::Frequency(f) = plot.result.kind() {
        ui.separator();
        ui.label(RichText::new("Modes").color(t.text_dim));
        egui::ScrollArea::vertical().id_salt("fc_results_modes").max_height(150.0).show(ui, |ui| {
            for (i, m) in f.modes.iter().enumerate() {
                let mode = i + 1;
                if ui
                    .selectable_label(
                        plot.view.mode == mode,
                        format!("{mode}   {}", crate::results::mode_value(f.buckling.is_some(), m.eigenvalue, m.frequency)),
                    )
                    .clicked()
                {
                    run.push(("results.modes", json!({"mode": mode})));
                }
            }
        });
    } else if let Kind::Heat(h) = plot.result.kind() {
        ui.separator();
        heat_flows(ui, &fieldcraft_engine::results::heat_flows(h), t);
    }
    if matches!(plot.result.kind(), Kind::Static(_)) {
        ui.separator();
        ui.horizontal(|ui| {
            if ui.button("Reaction sums").on_hover_text("Reaction force and moment per boundary-condition set (results.reactions)").clicked() {
                run.push(("results.reactions", json!({"step": plot.view.step})));
            }
            if ui.button("Export CSV…").on_hover_text("Write the reaction sums to a CSV file (results.export.csv)").clicked() {
                export = Some(json!({"what": "reactions", "step": plot.view.step}));
            }
        });
        if let Some(r) = app.reactions.as_ref().filter(|r| r["step"] == plot.view.step.as_str()) {
            reactions_table(ui, r, plot, t);
        }
    }
    let mut xy_export = None;
    if app.xy.is_some() {
        ui.separator();
        xy_export = xy_section(app, plot, ui, t, &mut run);
    }
    ui.separator();
    if let Some(p) = probed(app, plot, ui, t) {
        export = Some(p);
    }
    if let Some(p) = export {
        app.export_csv(p);
    }
    if let Some(p) = xy_export {
        app.export_csv_with("results.xy", p);
    }
    for (id, params) in run {
        // The reaction sums land in `app.reactions` (ModelApp::run_form).
        let _ = app.run(id, params);
    }
}

/// The deformation rows of the results grid: Deformed (with its scale) and Animate.
fn deformation_rows(plot: &Plot, ui: &mut egui::Ui, run: &mut Vec<(&'static str, Value)>) {
    ui.label("Deformed");
    ui.horizontal(|ui| {
        let mut deformed = plot.view.deformed;
        if ui.checkbox(&mut deformed, "").changed() {
            run.push(("results.deformed", json!({"show": deformed})));
        }
        let mut auto = plot.view.scale.is_none();
        let tip = if true_scale(&plot.result) {
            "Draw the displacements at true scale (they are large already)"
        } else {
            "Scale the largest displacement to a tenth of the model"
        };
        if ui.checkbox(&mut auto, "Auto").on_hover_text(tip).changed() {
            run.push(("results.deformed", json!({"scale": if auto { json!("auto") } else { json!(plot.scale) }})));
        }
        // The scale ticking Deformed draws with, also while it is off.
        let mut scale = plot.scale;
        let drag = egui::DragValue::new(&mut scale)
            .speed(plot.scale.abs().max(1e-6) * 0.01)
            .range(0.0..=1e12)
            .prefix("× ")
            .custom_formatter(|v, _| number(v))
            .custom_parser(|s| s.replace(',', "").trim().parse::<f64>().ok());
        if ui.add_enabled(plot.view.deformed, drag).changed() {
            run.push(("results.deformed", json!({"scale": scale})));
        }
    });
    ui.end_row();
    ui.label("Animate");
    let mut animate = plot.view.animate;
    if ui.checkbox(&mut animate, "").changed() {
        run.push(("results.animate", json!({"on": animate})));
    }
    ui.end_row();
}

/// A fluid step's rows of the results grid: the slice (its normal, position), the obstacles'
/// colours and the streamlines (their count from the inlet).
fn flow_rows(plot: &Plot, ui: &mut egui::Ui, run: &mut Vec<(&'static str, Value)>) {
    use fieldcraft_engine::flow::{Axis, slice_position};
    let (Kind::Fluid(f), Some(flow)) = (plot.result.kind(), plot.flow.as_ref()) else { return };
    if plot.view.contour == Some(Field::WallHeatFlux) {
        ui.label("Wall heat");
        ui.label("Positive into wall; face average");
        ui.end_row();
        if plot.contour.is_none() {
            ui.label("History");
            ui.label("No wall heat recorded at this frame");
            ui.end_row();
        }
    }
    let view = &flow.view;
    ui.label("Slice");
    ui.horizontal(|ui| {
        if plot.view.contour == Some(Field::WallHeatFlux) {
            ui.label("Hidden for wall heat");
            return;
        }
        let shown = view.slice.as_ref().map_or("None", |s| s.axis.id());
        egui::ComboBox::from_id_salt("fc_results_slice").selected_text(format!("Normal {}", shown.to_uppercase())).width(90.0).show_ui(ui, |ui| {
            if ui.selectable_label(view.slice.is_none(), "None").clicked() {
                run.push(("results.slice", json!({"axis": "none"})));
            }
            for a in Axis::ALL {
                if ui.selectable_label(view.slice.as_ref().is_some_and(|s| s.axis == a), format!("Normal {}", a.id().to_uppercase())).clicked() {
                    run.push(("results.slice", json!({"axis": a.id()})));
                }
            }
        });
        if let Some(s) = &view.slice {
            let (lo, hi) = (flow.grid.start(s.axis.index()), flow.grid.end(s.axis.index()));
            let mut at = slice_position(f, s);
            let drag = egui::DragValue::new(&mut at)
                .speed((hi - lo).abs().max(1e-9) * 0.002)
                .range(lo..=hi)
                .custom_formatter(|v, _| number(v))
                .custom_parser(|s| s.replace(',', "").trim().parse::<f64>().ok());
            if ui.add(drag).on_hover_text("Where the slice crosses its axis").changed() {
                run.push(("results.slice", json!({"position": at})));
            }
        }
    });
    ui.end_row();
    ui.label("Obstacles");
    let mut surface = view.surface;
    if ui.checkbox(&mut surface, "Coloured").on_hover_text("Colour the part by the selected fluid field or wall heat").changed() {
        run.push(("results.surface", json!({"show": surface})));
    }
    ui.end_row();
    ui.label("Streamlines");
    ui.horizontal(|ui| {
        let mut on = view.streamlines.is_some();
        if ui.checkbox(&mut on, "").on_hover_text("Streamlines from the inlet, coloured by speed or temperature when selected").changed() {
            run.push(("results.streamlines", json!({"show": on})));
        }
        if let Some(sl) = &view.streamlines {
            let mut count = sl.count;
            let drag = egui::DragValue::new(&mut count).range(1..=fieldcraft_engine::flow::MAX_STREAMLINES).speed(1.0).suffix(" seeds");
            if ui.add(drag).on_hover_text("Seeds on the inlet face, upstream of the obstacles").changed() {
                run.push(("results.streamlines", json!({"count": count})));
            }
        }
    });
    ui.end_row();
}

/// The probed node's values, one row each, or how to probe one; the export's parameters when
/// its Export CSV… is clicked.
fn probed(app: &ModelApp, plot: &Plot, ui: &mut egui::Ui, t: &Tokens) -> Option<Value> {
    let Some(p) = app.probe.as_ref().filter(|p| p.step == plot.view.step && p.mode == plot.pick.mode && p.frame == plot.pick.frame) else {
        ui.label(RichText::new("Click the model to probe a node").small().color(t.text_dim));
        return None;
    };
    ui.label(RichText::new(format!("Node {}", p.node)).strong());
    egui::Grid::new("fc_results_probe").num_columns(2).spacing(vec2(8.0, 2.0)).show(ui, |ui| {
        for (what, value) in &p.rows {
            ui.label(RichText::new(what).color(t.text_dim));
            ui.add(egui::Label::new(value).extend());
            ui.end_row();
        }
    });
    ui.button("Export CSV…").on_hover_text("Write the node's results to a CSV file (results.export.csv)").clicked().then(|| {
        let mut params = json!({"what": "probe", "step": p.step, "nodes": [p.node]});
        if matches!(plot.result.kind(), Kind::Frequency(_)) {
            params["mode"] = json!(p.mode);
        }
        if let Some(f) = p.frame {
            params["frame"] = json!(f);
        }
        params
    })
}

/// An explicit step: a frame slider (the deformed shape and contours follow it) and the energy
/// histories, plotted against time with a cursor at the frame shown (a click or drag moves it)
/// and read there.
fn explicit(ui: &mut egui::Ui, e: &ExplicitResult, plot: &Plot, t: &Tokens, run: &mut Vec<(&'static str, Value)>) {
    let time_unit = Quantity::Time.unit(plot.units);
    if let Some((k, time, n)) = plot.frame() {
        ui.horizontal(|ui| {
            ui.label("Frame");
            let mut frame = k;
            if ui.add(egui::Slider::new(&mut frame, 0..=n.saturating_sub(1))).changed() {
                run.push(("results.frame", json!({"frame": frame})));
            }
        });
        ui.label(RichText::new(format!("t = {}", with_unit(number(time), time_unit.as_deref()))).color(t.text_dim));
    }
    let h = &e.history;
    if h.time.is_empty() {
        return;
    }
    let series: [(&str, Color32, fn(&Energies) -> f64); 4] = [
        ("Kinetic", t.accent, |e| e.kinetic),
        ("Internal", t.ok, |e| e.internal),
        ("Hourglass", t.warm, |e| e.hourglass),
        ("External work", t.error, |e| e.external),
    ];
    let energy_unit = Quantity::Energy.unit(plot.units);
    ui.label(RichText::new(titled("Energies", energy_unit.as_deref())).color(t.text_dim));
    let lines: Vec<Line> = series
        .iter()
        .map(|(_, color, f)| Line { color: *color, points: h.time.iter().zip(&h.energies).map(|(time, e)| (*time, f(e))).collect() })
        .collect();
    let t_end = h.time.last().copied().unwrap_or(0.0);
    let now = plot.frame().map_or(t_end, |(_, time, _)| time);
    let (response, axes) = line_plot(ui, &lines, "t", time_unit.as_deref(), Some(now), t);
    let scale = axes.lo.abs().max(axes.hi.abs());
    let nearest = |time: f64| h.time.iter().enumerate().min_by(|a, b| (a.1 - time).abs().total_cmp(&(b.1 - time).abs())).map_or(0, |(i, _)| i);
    // A click or drag on the plot shows the frame nearest that time.
    if (response.clicked() || response.dragged())
        && let Some(p) = response.interact_pointer_pos()
    {
        run.push(("results.frame", json!({"time": axes.x_of(p.x)})));
    }
    // Hover reads the curves at the pointer's time.
    if let Some(p) = response.hover_pos().filter(|p| axes.rect.contains(*p))
        && let Some(e) = h.energies.get(nearest(axes.x_of(p.x)))
    {
        let mut text = vec![format!("t = {}", with_unit(number(axes.x_of(p.x)), time_unit.as_deref()))];
        text.extend(series.iter().map(|(name, _, f)| format!("{name} {}", number_in(f(e), scale))));
        response.on_hover_text(text.join("\n"));
    }
    // The values at the frame's time.
    let Some(at_now) = h.energies.get(nearest(now)) else { return };
    let frame = plot.frame().map_or_else(String::new, |(k, ..)| format!(", frame {k}"));
    ui.label(RichText::new(format!("At t = {}{frame}", with_unit(number(now), time_unit.as_deref()))).small().color(t.text_dim));
    egui::Grid::new("fc_results_energies").num_columns(2).spacing(vec2(8.0, 2.0)).show(ui, |ui| {
        for (name, color, f) in &series {
            ui.label(RichText::new(*name).color(*color));
            ui.label(with_unit(number_in(f(at_now), scale), energy_unit.as_deref()));
            ui.end_row();
        }
        ui.label(RichText::new("Balance").color(t.text_dim)).on_hover_text("The energy balance's relative error");
        ui.label(number(at_now.balance));
        ui.end_row();
    });
}

/// A line of a plot: its colour and its (x, y) points (those not finite are left out).
struct Line {
    color: Color32,
    points: Vec<(f64, f64)>,
}

/// Where a plot's values land: its drawing area and the x and value ranges it spans.
struct Axes {
    rect: Rect,
    x0: f64,
    x1: f64,
    lo: f64,
    hi: f64,
}

impl Axes {
    fn at(&self, x: f64, y: f64) -> Pos2 {
        let (w, h) = ((self.x1 - self.x0).max(f64::MIN_POSITIVE), (self.hi - self.lo).max(f64::MIN_POSITIVE));
        pos2(
            self.rect.left() + ((x - self.x0) / w) as f32 * self.rect.width(),
            self.rect.bottom() - 2.0 - ((y - self.lo) / h) as f32 * (self.rect.height() - 4.0),
        )
    }

    /// The x value under the screen x `px` (clamped to the range).
    fn x_of(&self, px: f32) -> f64 {
        self.x0 + f64::from(((px - self.rect.left()) / self.rect.width()).clamp(0.0, 1.0)) * (self.x1 - self.x0)
    }
}

/// A line plot the panel's width: the lines (at most [`MAX_PLOT_POINTS`] points each), the
/// value range left (spanning 0, with a 0 line when it crosses it), the x range under it named
/// `x_name` with `x_unit`, and a cursor at x `cursor`. The response senses clicks and drags.
fn line_plot(ui: &mut egui::Ui, lines: &[Line], x_name: &str, x_unit: Option<&str>, cursor: Option<f64>, t: &Tokens) -> (egui::Response, Axes) {
    let (outer, response) = ui.allocate_exact_size(vec2(ui.available_width(), PLOT_HEIGHT + AXIS_BOTTOM), egui::Sense::click_and_drag());
    let rect = Rect::from_min_max(outer.min + vec2(AXIS_LEFT, 0.0), outer.max - vec2(0.0, AXIS_BOTTOM));
    let painter = ui.painter_at(outer);
    painter.rect_filled(rect, 2.0, t.field);
    let finite = || lines.iter().flat_map(|l| l.points.iter()).filter(|(x, y)| x.is_finite() && y.is_finite());
    let (mut x0, mut x1) = finite().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), (x, _)| (a.min(*x), b.max(*x)));
    if !(x0 <= x1) {
        (x0, x1) = (0.0, 1.0);
    } else if x0 == x1 {
        x1 = x0 + 1.0;
    }
    let hi = finite().map(|p| p.1).fold(0.0, f64::max);
    let lo = finite().map(|p| p.1).fold(0.0, f64::min);
    let hi = if hi > lo { hi } else { lo + 1.0 };
    let axes = Axes { rect, x0, x1, lo, hi };
    let small = FontId::proportional(10.0);
    let scale = hi.abs().max(lo.abs());
    painter.text(pos2(rect.left() - 4.0, rect.top()), Align2::RIGHT_TOP, number_in(hi, scale), small.clone(), t.text_dim);
    painter.text(pos2(rect.left() - 4.0, rect.bottom()), Align2::RIGHT_BOTTOM, number_in(lo, scale), small.clone(), t.text_dim);
    if lo < 0.0 && hi > 0.0 {
        let zero = axes.at(x0, 0.0).y;
        painter.line_segment([pos2(rect.left(), zero), pos2(rect.right(), zero)], Stroke::new(1.0, t.border));
    }
    painter.text(pos2(rect.left(), rect.bottom() + 2.0), Align2::LEFT_TOP, number(x0), small.clone(), t.text_dim);
    let end = format!("{x_name} {}", with_unit(number(x1), x_unit));
    painter.text(pos2(rect.right(), rect.bottom() + 2.0), Align2::RIGHT_TOP, end, small, t.text_dim);
    for line in lines {
        let stride = (line.points.len() / MAX_PLOT_POINTS).max(1);
        let points: Vec<Pos2> =
            line.points.iter().step_by(stride).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| axes.at(*x, *y)).collect();
        // A single point (a steady heat step) shows as a dot.
        if let [p] = points.as_slice() {
            painter.circle_filled(*p, 2.5, line.color);
        }
        painter.add(egui::Shape::line(points, Stroke::new(1.5, line.color)));
    }
    if let Some(x) = cursor.filter(|x| x.is_finite()) {
        let x = axes.at(x, 0.0).x.clamp(rect.left(), rect.right());
        painter.line_segment([pos2(x, rect.top()), pos2(x, rect.bottom())], Stroke::new(1.5, t.text));
    }
    if let Some(p) = response.hover_pos().filter(|p| rect.contains(*p)) {
        painter.line_segment([pos2(p.x, rect.top()), pos2(p.x, rect.bottom())], Stroke::new(1.0, t.border));
    }
    (response, axes)
}

/// The quantities the XY plot offers for a result, as (label, `results.xy`'s `y`); a nodal
/// field's `node` is filled in from the probed node.
pub fn xy_options(r: &fieldcraft_engine::StepResult) -> Vec<(String, Value)> {
    use fieldcraft_engine::cmd::xy;
    let mut out = Vec::new();
    if let Kind::Frequency(f) = r.kind() {
        out.push(if f.buckling.is_some() {
            ("Buckling load factor".to_string(), json!("load_factor"))
        } else {
            ("Frequency".to_string(), json!("frequency"))
        });
    }
    if let Kind::Explicit(e) = r.kind() {
        out.extend(xy::ENERGIES.iter().map(|n| (format!("Energy: {}", n.replace('_', " ")), json!({"energy": n}))));
        out.extend(xy::HISTORIES.iter().map(|n| (n.replace('_', " "), json!({"history": n}))));
        for node in e.history.nodes.iter().take(8) {
            out.extend(
                xy::NODE_QUANTITIES
                    .iter()
                    .map(|q| (format!("History node {node}: {}", q.replace('_', " ")), json!({"history_node": node, "quantity": q}))),
            );
        }
    }
    let linear = matches!(r.kind(), Kind::Static(s) if s.nonlinear.is_none());
    if !linear {
        out.extend(xy::plottable(r).into_iter().map(|f| (format!("{} at a node", f.label()), json!({"field": f.id()}))));
    }
    out
}

/// The XY plot (results.xy): a quantity picker, the series plotted against the step's own axis
/// with the frame (mode) shown marked (a click shows the nearest one), and its export's
/// parameters when Export CSV… is clicked. A nodal field is read at the probed node.
fn xy_section(app: &mut ModelApp, plot: &Plot, ui: &mut egui::Ui, t: &Tokens, run: &mut Vec<(&'static str, Value)>) -> Option<Value> {
    let step = plot.view.step.clone();
    let mut close = false;
    ui.horizontal(|ui| {
        ui.label(RichText::new("XY Plot").strong());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            close = ui.small_button("×").on_hover_text("Close the XY plot").clicked();
        });
    });
    if close {
        app.xy = None;
        return None;
    }
    let options = xy_options(&plot.result);
    if options.is_empty() {
        ui.label(
            RichText::new("A linear static step has one state: plot a nonlinear static, explicit, heat or frequency step").small().color(t.text_dim),
        );
        return None;
    }
    let probed = app.probe.as_ref().filter(|p| p.step == step).map(|p| p.node);
    let xy = app.xy.get_or_insert_with(Default::default);
    if probed.is_some() {
        xy.node = probed;
    }
    let template = xy.y.clone().filter(|y| options.iter().any(|(_, o)| o == y)).unwrap_or_else(|| options[0].1.clone());
    let shown = options.iter().find(|(_, o)| *o == template).map_or("", |(l, _)| l.as_str());
    let mut picked = template.clone();
    egui::ComboBox::from_id_salt("fc_results_xy").selected_text(shown).width(ui.available_width() - 8.0).truncate().show_ui(ui, |ui| {
        for (label, y) in &options {
            if ui.selectable_label(*y == template, label).clicked() {
                picked = y.clone();
            }
        }
    });
    xy.y = Some(picked.clone());
    let y = if picked.get("field").is_some() {
        let Some(node) = xy.node else {
            (xy.key, xy.series) = (None, None);
            ui.label(RichText::new("Click the model to pick the node").small().color(t.text_dim));
            return None;
        };
        let mut y = picked;
        y["node"] = json!(node);
        y
    } else {
        picked
    };
    let key = (step.clone(), Arc::as_ptr(&plot.result) as usize, y.clone());
    if xy.key.as_ref() != Some(&key) {
        xy.key = Some(key);
        xy.series = None;
        run.push(("results.xy", json!({"step": step, "y": y})));
        return None;
    }
    let s = xy.series.clone().filter(|s| s["step"] == step.as_str())?;
    let floats = |v: &Value| -> Vec<f64> { v.as_array().into_iter().flatten().map(|x| x.as_f64().unwrap_or(f64::NAN)).collect() };
    let (xs, ys) = (floats(&s["x"]["values"]), floats(&s["y"]["values"]));
    let frames: Vec<usize> = s["frames"].as_array().into_iter().flatten().filter_map(|f| f.as_u64().and_then(|f| usize::try_from(f).ok())).collect();
    let axis = s["x"]["axis"].as_str().unwrap_or_default();
    // The point shown: the frame's, or the mode's.
    let nearest =
        |x: f64| xs.iter().enumerate().filter(|(_, v)| v.is_finite()).min_by(|a, b| (a.1 - x).abs().total_cmp(&(b.1 - x).abs())).map(|(i, _)| i);
    let current = match plot.frame() {
        Some((k, ..)) if !frames.is_empty() => frames.iter().position(|f| *f == k),
        // A history against time: the sample nearest the frame's time.
        Some((_, time, _)) if axis == "time" => nearest(time),
        _ if axis == "mode" => Some(plot.view.mode.saturating_sub(1)),
        _ => None,
    };
    ui.label(RichText::new(s["y"]["label"].as_str().unwrap_or_default()).small().color(t.text_dim));
    let x_name = match axis {
        "load_factor" => "λ",
        "time" => "t",
        _ => "mode",
    };
    let x_unit = (axis == "time").then(|| Quantity::Time.unit(plot.units)).flatten();
    let line = Line { color: t.accent, points: xs.iter().copied().zip(ys.iter().copied()).collect() };
    let (response, axes) = line_plot(ui, &[line], x_name, x_unit.as_deref(), current.and_then(|i| xs.get(i).copied()), t);
    let x_label = s["x"]["label"].as_str().unwrap_or_default().to_string();
    if response.clicked()
        && let Some(i) = response.interact_pointer_pos().and_then(|p| nearest(axes.x_of(p.x)))
    {
        if let Some(f) = frames.get(i) {
            run.push(("results.frame", json!({"step": step, "frame": f})));
        } else if axis == "time" {
            run.push(("results.frame", json!({"step": step, "time": xs.get(i)})));
        } else if axis == "mode" {
            run.push(("results.modes", json!({"step": step, "mode": i + 1})));
        }
    }
    if let Some(i) = response.hover_pos().filter(|p| axes.rect.contains(*p)).and_then(|p| nearest(axes.x_of(p.x))) {
        let (x, y) = (xs.get(i).copied().unwrap_or(f64::NAN), ys.get(i).copied().unwrap_or(f64::NAN));
        let scale = axes.lo.abs().max(axes.hi.abs());
        response.on_hover_text(format!("{x_label} {}\n{}", number(x), number_in(y, scale)));
    }
    if let Some(i) = current {
        let at = |v: &[f64]| v.get(i).map_or("—".into(), |x| number(*x));
        ui.label(RichText::new(format!("{x_label} {}: {}", at(&xs), at(&ys))).small().color(t.text_dim));
    }
    ui.button("Export CSV…").on_hover_text("Write the series to a CSV file (results.xy)").clicked().then(|| json!({"step": step, "y": y}))
}

impl ModelApp {
    /// XY Plot from the ribbon: the plot opens in the results panel, on the step shown (the
    /// first solved one when none is).
    pub(crate) fn open_xy(&mut self) {
        if self.session.results_view.is_none() {
            let _ = self.run("results.deformed", json!({}));
        }
        self.xy.get_or_insert_with(Default::default);
    }
}

/// A heat step's balance: the heat each heat load carries in (negative: out), the prescribed
/// temperatures' total, and what is left over (about 0 when steady).
fn heat_flows(ui: &mut egui::Ui, flows: &Value, t: &Tokens) {
    ui.label(RichText::new("Heat flow into the body").color(t.text_dim));
    egui::Grid::new("fc_results_heat").num_columns(2).spacing(vec2(8.0, 2.0)).show(ui, |ui| {
        let row = |ui: &mut egui::Ui, name: &str, v: Option<f64>| {
            ui.label(name);
            ui.label(v.map_or("—".into(), number));
            ui.end_row();
        };
        for l in flows["loads"].as_array().into_iter().flatten() {
            row(ui, &format!("{} ({})", l["name"].as_str().unwrap_or("?"), l["type"].as_str().unwrap_or("").replace('_', " ")), l["heat"].as_f64());
        }
        row(ui, "Prescribed temperatures", flows["reactions"].as_f64());
        row(ui, "Imbalance", flows["imbalance"].as_f64());
    });
}

/// The reaction sums: per set and over the model, the force and the moment about the point, and
/// the applied total beside the reacted one (they balance in equilibrium). Round-off below a
/// billionth of the largest force (moment) is 0.
fn reactions_table(ui: &mut egui::Ui, r: &Value, plot: &Plot, t: &Tokens) {
    let vector = |v: &Value| -> Option<[f64; 3]> { serde_json::from_value(v.clone()).ok() };
    let rows: Vec<(String, Option<[f64; 3]>, Option<[f64; 3]>)> = r["sets"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|s| {
            let nodes = s["nodes"].as_u64().unwrap_or(0);
            (format!("{} ({} nodes)", s["set"].as_str().unwrap_or("?"), count(nodes)), vector(&s["force"]), vector(&s["moment"]))
        })
        .chain(std::iter::once(("All nodes".to_string(), vector(&r["total"]["force"]), vector(&r["total"]["moment"]))))
        .collect();
    let applied = vector(&r["applied"]);
    let largest = |vs: &mut dyn Iterator<Item = [f64; 3]>| vs.flatten().map(f64::abs).filter(|x| x.is_finite()).fold(0.0, f64::max);
    let force_scale = largest(&mut rows.iter().filter_map(|r| r.1).chain(applied));
    let moment_scale = largest(&mut rows.iter().filter_map(|r| r.2));
    let (force, moment) = (Quantity::Force.unit(plot.units), Quantity::Moment.unit(plot.units));
    let about = vector(&r["about"]).unwrap_or_default();
    ui.label(RichText::new(format!("Moments about ({})", components(&about, 0.0))).small().color(t.text_dim));
    let cells = |ui: &mut egui::Ui, v: Option<[f64; 3]>, scale: f64| {
        for k in 0..3 {
            ui.add(egui::Label::new(v.and_then(|v| v.get(k).copied()).map_or("—".into(), |x| number_in(x, scale))).extend());
        }
    };
    egui::Grid::new("fc_results_reactions").num_columns(4).spacing(vec2(8.0, 2.0)).striped(true).show(ui, |ui| {
        for h in ["", "X", "Y", "Z"] {
            ui.label(RichText::new(h).small().color(t.text_dim));
        }
        ui.end_row();
        for (name, f, m) in &rows {
            ui.label(RichText::new(name).strong());
            ui.end_row();
            ui.label(RichText::new(titled("F", force.as_deref())).color(t.text_dim));
            cells(ui, *f, force_scale);
            ui.end_row();
            ui.label(RichText::new(titled("M", moment.as_deref())).color(t.text_dim));
            cells(ui, *m, moment_scale);
            ui.end_row();
        }
        if applied.is_some() {
            ui.label(RichText::new(titled("Applied F", force.as_deref())).color(t.text_dim));
            cells(ui, applied, force_scale);
            ui.end_row();
        }
    });
    // Equilibrium: the reactions over the whole model balance the applied loads.
    if let (Some(total), Some(applied)) = (rows.last().and_then(|r| r.1), applied) {
        let residual = (0..3).map(|k| (total[k] + applied[k]).abs()).fold(0.0, f64::max);
        let balanced = residual <= 1e-6 * force_scale.max(f64::MIN_POSITIVE);
        let (text, color) = if balanced {
            ("✔ Reactions balance the applied loads".to_string(), t.ok)
        } else {
            (format!("Reactions and applied loads differ by {}", with_unit(number(residual), force.as_deref())), t.warm)
        };
        ui.label(RichText::new(text).small().color(color));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A label goes right of its point, left when right runs off, inside the top and bottom;
    /// one that would cover a placed label goes under it, or over it at the bottom edge.
    #[test]
    fn labels_keep_inside_the_view_and_off_each_other() {
        let bounds = Rect::from_min_max(pos2(0.0, 0.0), pos2(400.0, 300.0));
        let size = vec2(80.0, 20.0);
        let a = label_box(pos2(100.0, 150.0), size, bounds, None);
        assert_eq!(a, Rect::from_min_size(pos2(100.0 + LABEL_GAP, 140.0), size));
        assert_eq!(label_box(pos2(380.0, 150.0), size, bounds, None).min, pos2(380.0 - LABEL_GAP - 80.0, 140.0));
        assert_eq!(label_box(pos2(100.0, 2.0), size, bounds, None).min.y, 0.0);
        assert_eq!(label_box(pos2(100.0, 299.0), size, bounds, None).max.y, 300.0);
        // Far apart: unmoved. Overlapping: under the placed one, apart by the label padding.
        assert_eq!(label_box(pos2(100.0, 250.0), size, bounds, Some(a)), Rect::from_min_size(pos2(100.0 + LABEL_GAP, 240.0), size));
        let b = label_box(pos2(104.0, 156.0), size, bounds, Some(a));
        assert_eq!(b.min, pos2(104.0 + LABEL_GAP, a.max.y + LABEL_PAD.y));
        assert!(!b.intersects(a));
        // Placed at the bottom edge: the second goes over it.
        let low = label_box(pos2(100.0, 299.0), size, bounds, None);
        let over = label_box(pos2(100.0, 295.0), size, bounds, Some(low));
        assert_eq!(over.max.y, low.min.y - LABEL_PAD.y);
        assert!(!over.intersects(low));
    }
}
