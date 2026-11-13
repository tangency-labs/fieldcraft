//! What the viewport shows over a results plot: the legend (field and unit, bands, range,
//! deformation scale, a stale-result warning), markers at the contour's largest and smallest
//! value, and the probed node; and the results controls of the context panel (right), which run
//! the `results.*` commands (a fluid step's: its slice, coloured obstacles, streamlines and
//! frames). Numbers go through one formatter ([`crate::results::number`]), in the display
//! units ([`Plot::units`]): the result is stored in SI and converted here as it is shown.

use egui::{Align2, Color32, FontId, Pos2, Rect, Stroke, pos2, vec2};
use fieldcraft_engine::results::{ComponentCheck, Field, FluidResult, Kind, Quantity, shown as shown_in, true_scale, unit_of, with_unit as titled};
use fieldcraft_explicit::{Energies, ExplicitResult};
use serde_json::{Value, json};

use std::sync::Arc;

use crate::ModelApp;
use crate::results::{Plot, components, count, number, number_in, with_unit};
use crate::theme::Tokens;
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

const SWATCH: egui::Vec2 = vec2(18.0, 15.0);
/// The smallest span a drag speed is taken of, so a zero range still drags.
const MIN_DRAG_SPAN: f64 = 1e-9; // unit-factors: a floor on a span in any unit
/// The smallest scale factor a drag speed is taken of.
const MIN_DRAG_SCALE: f64 = 1e-6; // unit-factors: a floor on a scale factor
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
    let bar = crate::viewport::selection_bar_rect(ui.ctx());
    let extremes = legend(&plot, &painter, rect, bar, t);
    // The markers and their labels keep off the legend: its Max and Min lines name a node under it.
    let covered = plot.contour.as_ref().map(|_| legend_layout(&plot, &painter, rect, bar, t).frame);
    if app.ui.result_extremes {
        markers(app, &plot, &painter, rect, covered, t);
    }
    probe(app, &plot, &painter, rect, covered, t);
    let mut focus = None;
    for (i, (line, node)) in extremes.into_iter().enumerate() {
        let tip = if plot.is_fluid() { "Centre this field extreme".into() } else { format!("Show and probe node {node}") };
        let r = ui.interact(line, ui.scope_id().with(("fc_legend_extreme", i)), egui::Sense::click()).backing_on_hover_text(tip);
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
    band_height: f32,
}

const LEGEND_LINE: f32 = 16.0;
/// The room the legend keeps from the selection bar (points).
const LEGEND_GAP: f32 = 6.0;

/// Where the legend of `plot` goes in the viewport `rect`, clear of the selection `bar` (a fit
/// leaves it free).
pub fn legend_rect(plot: &Plot, painter: &egui::Painter, rect: Rect, bar: Option<Rect>, t: &Tokens) -> Rect {
    legend_layout(plot, painter, rect, bar, t).frame
}

/// The legend's title: the field with its unit (a mode shape normalised), or the shape drawn.
fn title(plot: &Plot) -> String {
    if plot.unavailable {
        return plot.view.contour.map_or_else(|| "Field unavailable".into(), |f| f.title(&plot.result, &plot.units));
    }
    match (&plot.contour, plot.mode()) {
        (Some(l), _) => l.field.title(&plot.result, &plot.units),
        (None, Some(_)) => "Mode shape (normalised)".into(),
        (None, None) if plot.is_fluid() => plot.view.contour.map_or_else(|| "Flow".into(), |f| f.title(&plot.result, &plot.units)),
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
    let unit = l.field.unit(&plot.result, &plot.units);
    let scale = l.max.1.abs().max(l.min.1.abs());
    let at = |node: u32| match plot.extreme_position(node).filter(|_| plot.is_fluid()) {
        Some(p) => {
            let at = fieldcraft_engine::results::shown_all(Quantity::Length, [p.x, p.y, p.z], &plot.units);
            format!("at {}", with_unit(components(&at, 0.0), unit_of(Quantity::Length, &plot.units).as_deref()))
        }
        None => format!("at node {node}"),
    };
    let line = |what: &str, (node, v): (u32, f64)| (format!("{what} {} {}", with_unit(number_in(v, scale), unit.as_deref()), at(node)), node);
    Some([line("Max", l.max), line("Min", l.min)])
}

/// The legend's frame of `n`, `k` counted from 0 like the slider and `results.frame`: a time
/// step's `frame k of n` (frame 0 the initial field), a nonlinear static step's `increment k of
/// n − 1` at its load factor (increment 0 the unloaded state).
fn frame_text(increments: bool, k: usize, n: usize, at: f64, time_unit: Option<&str>) -> String {
    if increments {
        format!(", increment {k} of {}: load factor {}", n.saturating_sub(1), number(at))
    } else {
        format!(", frame {k} of {n}: t = {}", with_unit(number(at), time_unit))
    }
}

/// The legend's text and frame: top left (below the selection bar when they would meet), as wide
/// as its longest line. All bands fit above navigation; short plots thin their tick labels.
fn legend_layout(plot: &Plot, painter: &egui::Painter, rect: Rect, bar: Option<Rect>, t: &Tokens) -> Legend {
    let mut lines: Vec<(String, Color32)> = Vec::new();
    let title = title(plot);
    let mut step = format!("Step {}", plot.view.step);
    if let Some((mode, value)) = plot.mode() {
        step.push_str(&format!(", mode {mode}: {value}"));
    }
    if let Some((k, at, n)) = plot.frame() {
        let increments = matches!(plot.result.kind(), Kind::Static(_));
        step.push_str(&frame_text(increments, k, n, at, unit_of(Quantity::Time, &plot.units).as_deref()));
    }
    if let Some(a) = plot.amplitude {
        let what = if matches!(plot.result.kind(), Kind::Frequency(_)) { "Mode amplitude" } else { "Proportional load" };
        lines.push((format!("{what} × {} (display interpolation)", number(a)), t.text_dim));
    }
    lines.push((step, t.overlay_text));
    if !plot.is_fluid() {
        lines.push((deformation(plot), t.text_dim));
    }
    if plot.unavailable {
        lines.push(("Field not saved at this increment".into(), t.warm));
        lines.push(("No stress output, or past the step's field cap".into(), t.warm));
    }
    if let Some(l) = plot.contour.as_ref().filter(|_| plot.view.range_min.is_some() || plot.view.range_max.is_some()) {
        let scale = l.max.1.abs().max(l.min.1.abs());
        lines.push((format!("Legend clamped to {} … {}", number_in(l.bands.min, scale), number_in(l.bands.max, scale)), t.warm));
    }
    if plot.stale {
        lines.push(("Stale: the model changed since this step was solved".into(), t.warm));
    }
    let k = solidcraft_theme::Tokens::get();
    let origin = rect.left_top() + vec2(k.space_12, k.space_12);
    let font = FontId::proportional(k.font_body);
    let line_h = LEGEND_LINE;
    let bands = plot.contour.as_ref().map(|l| &l.bands);
    // As wide as its longest line (and the extremes' lines under the swatches).
    let extremes = extreme_lines(plot).into_iter().flatten().map(|(text, _)| text);
    let widest = lines
        .iter()
        .map(|(text, _)| text.clone())
        .chain(extremes)
        .map(|text| painter.layout_no_wrap(text, font.clone(), t.overlay_text).size().x)
        .fold(0.0_f32, f32::max);
    let title_w = painter.layout_no_wrap(title.clone(), FontId::proportional(k.font_heading), t.overlay_text).size().x;
    let width = (widest.max(title_w) + 20.0).max(170.0);
    // Padding and title, the lines, then the swatches with the range lines under them.
    let fixed_height = 36.0 + line_h * lines.len() as f32 + if bands.is_some() { 16.0 + 2.0 * line_h } else { 0.0 };
    let height = fixed_height + bands.map_or(0.0, |b| swatches(b) as f32 * SWATCH.y);
    let mut frame = Rect::from_min_size(origin, vec2(width, height));
    // A viewport too narrow for the legend beside the selection bar puts it below the bar (#559).
    if let Some(b) = bar.filter(|b| b.expand(LEGEND_GAP).intersects(frame)) {
        frame = Rect::from_min_size(pos2(origin.x, b.bottom() + LEGEND_GAP), frame.size());
    }
    // Reserve the full-size title, metadata and clickable extremes before giving the bands
    // their room. Keep the over-limit band full height so its true/clamped maxima stay legible.
    let nav = crate::viewport::nav_bar_rect(rect);
    let bottom = if frame.left() < nav.right() && frame.right() > nav.left() { nav.top() - LEGEND_GAP } else { rect.bottom() - 12.0 };
    let over_height = bands.filter(|b| b.over).map_or(0.0, |_| SWATCH.y);
    let band_height = bands.map_or(SWATCH.y, |b| ((bottom - frame.top() - fixed_height - over_height) / b.count.max(1) as f32).clamp(1.0, SWATCH.y));
    frame.max.y = frame.top() + fixed_height + over_height + bands.map_or(0.0, |b| b.count as f32 * band_height);
    Legend { frame, title, lines, band_height }
}

/// The legend's swatches: the bands, and the "over" band above a clamped maximum.
fn swatches(b: &fieldcraft_render::contour::Bands) -> u32 {
    b.count + u32::from(b.over)
}

/// The legend, top left under the view buttons. Returns the Max and Min lines' rects with their
/// nodes.
fn legend(plot: &Plot, painter: &egui::Painter, rect: Rect, bar: Option<Rect>, t: &Tokens) -> Vec<(Rect, u32)> {
    let Legend { frame, title, lines, band_height } = legend_layout(plot, painter, rect, bar, t);
    let origin = frame.min;
    let k = solidcraft_theme::Tokens::get();
    let (font, line_h) = (FontId::proportional(k.font_body), LEGEND_LINE);
    painter.rect_filled(frame, k.radius_sm, t.overlay);
    let mut y = origin.y + 8.0;
    painter.text(pos2(origin.x + 10.0, y), Align2::LEFT_TOP, title, FontId::proportional(k.font_heading), t.overlay_text);
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
    let swatch_top = y;
    let swatch_bottom = y + l.bands.count as f32 * band_height + if l.bands.over { SWATCH.y } else { 0.0 };
    let mut last_label = y;
    // The values are right-aligned on their widest, so their decimal points line up.
    let highest = if l.bands.over { l.max.1 } else { l.bands.level(l.bands.count) };
    let value_w = (0..=l.bands.count)
        .map(|k| l.bands.level(k))
        .chain([highest])
        .map(|v| painter.layout_no_wrap(number_in(v, scale), font.clone(), t.overlay_text).size().x)
        .fold(0.0_f32, f32::max);
    for k in (0..rows).rev() {
        let h = if l.bands.is_over(k) { SWATCH.y } else { band_height };
        let r = Rect::from_min_size(pos2(origin.x + 10.0, y), vec2(SWATCH.x, h));
        let c = fieldcraft_render::contour::bands_color(&stops, over, &l.bands, k);
        painter.rect_filled(r, 0.0, Color32::from_rgb(c.0, c.1, c.2));
        if k + 1 == rows {
            let top = if l.bands.is_over(k) { l.max.1 } else { l.bands.level(k + 1) };
            painter.text(pos2(r.right() + 6.0 + value_w, r.top()), Align2::RIGHT_CENTER, number_in(top, scale), font.clone(), t.overlay_text);
        }
        // Keep both end values and leave a readable gap between the remaining tick labels.
        // The selected contour bands and their colours never change when the window shrinks.
        if k == 0 || (r.bottom() - last_label >= SWATCH.y && swatch_bottom - r.bottom() >= SWATCH.y) {
            painter.text(
                pos2(r.right() + 6.0 + value_w, r.bottom()),
                Align2::RIGHT_CENTER,
                number_in(l.bands.level(k), scale),
                font.clone(),
                t.overlay_text,
            );
            last_label = r.bottom();
        }
        y = r.bottom();
    }
    painter.rect_stroke(
        Rect::from_min_max(pos2(origin.x + 10.0, swatch_top), pos2(origin.x + 10.0 + SWATCH.x, swatch_bottom)),
        0.0,
        Stroke::new(1.0, t.border),
        egui::StrokeKind::Outside,
    );
    y += 10.0;
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

/// Where the label of a marker at `at` goes ([`label_box`]) clear of the legend and of the label
/// `placed` already: under (or over) the legend when beside it would cover it. None when the
/// point is under the legend, or when no clear place is left.
fn clear_label(at: Pos2, size: egui::Vec2, bounds: Rect, legend: Option<Rect>, placed: Option<Rect>) -> Option<Rect> {
    if legend.is_some_and(|l| l.contains(at)) {
        return None;
    }
    let mut r = label_box(at, size, bounds, placed);
    if let Some(l) = legend.filter(|l| l.intersects(r)) {
        r = label_box(at, size, bounds, Some(l));
    }
    let covers = |other: Option<Rect>| other.is_some_and(|o| o.intersects(r));
    (!covers(legend) && !covers(placed)).then_some(r)
}

/// A box of `text` beside the screen point `at`, clear of the legend and of `placed`
/// ([`clear_label`]); returns where it went, or None when it was left out.
fn label_beside(
    painter: &egui::Painter,
    at: Pos2,
    text: String,
    bounds: Rect,
    legend: Option<Rect>,
    placed: Option<Rect>,
    t: &Tokens,
) -> Option<Rect> {
    let k = solidcraft_theme::Tokens::get();
    let galley = painter.layout_no_wrap(text, FontId::proportional(k.font_body), t.overlay_text);
    let box_rect = clear_label(at, galley.size() + 2.0 * LABEL_PAD, bounds, legend, placed)?;
    painter.rect_filled(box_rect, k.radius_sm, t.overlay);
    painter.galley(box_rect.left_top() + LABEL_PAD, galley, t.overlay_text);
    Some(box_rect)
}

/// Markers with their values at the contour's largest and smallest value.
fn markers(app: &ModelApp, plot: &Plot, painter: &egui::Painter, rect: Rect, legend: Option<Rect>, t: &Tokens) {
    let Some(l) = &plot.contour else { return };
    let unit = l.field.unit(&plot.result, &plot.units);
    let scale = l.max.1.abs().max(l.min.1.abs());
    let (high, low) = (t.contour[4], t.contour[0]);
    // The Min label keeps off the Max label when the two cells are close on screen.
    let mut placed = None;
    for (what, (node, v), color) in [("Max", l.max, high), ("Min", l.min, low)] {
        let Some(at) = plot.extreme_position(node).and_then(|w| crate::viewport::to_screen(app, w)) else { continue };
        if !rect.contains(at) || legend.is_some_and(|l| l.contains(at)) {
            continue;
        }
        painter.circle(at, 5.0, color, Stroke::new(2.0, t.overlay_text));
        let text = format!("{what} {}", with_unit(number_in(v, scale), unit.as_deref()));
        placed = label_beside(painter, at, text, rect, legend, placed, t).or(placed);
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
fn probe(app: &ModelApp, plot: &Plot, painter: &egui::Painter, rect: Rect, legend: Option<Rect>, t: &Tokens) {
    let Some(p) = app.probe.as_ref().filter(|p| p.step == plot.view.step && p.mode == plot.pick.mode && p.frame == plot.pick.frame) else { return };
    let Some(at) = plot.node_position(p.node).and_then(|w| crate::viewport::to_screen(app, w)) else { return };
    if legend.is_some_and(|l| l.contains(at)) {
        return;
    }
    painter.circle(at, 5.0, t.accent, Stroke::new(2.0, t.overlay_text));
    label_beside(painter, at, p.lines.join("\n"), rect, legend, None, t);
}

/// The results controls (the context panel's Results section): step, contour field,
/// deformation, modes, animation, what the viewport marks, reaction sums and the probed node.
/// Every control runs a command or a view setting.
pub fn controls(app: &mut ModelApp, plot: &Plot, ui: &mut egui::Ui, t: &Tokens) {
    let mut run: Vec<(&'static str, Value)> = Vec::new();
    let mut export: Option<Value> = None;
    let steps: Vec<(String, bool)> = app.session.result_steps().map(|(n, c)| (n.to_string(), c)).collect();
    kit::row(ui, "Step", |ui| {
        kit::bounded_dropdown_menu(ui, "fc_results_step", &plot.view.step, ui.available_width(), |ui| {
            for (name, current) in &steps {
                let label = if *current { name.clone() } else { format!("{name} (stale)") };
                if ui.selectable_label(*name == plot.view.step, label).clicked() && *name != plot.view.step {
                    run.push(("results.deformed", json!({"step": name})));
                }
            }
        })
    });
    let shown = plot.view.contour.map_or("None", Field::label);
    kit::row(ui, "Contours", |ui| {
        kit::bounded_dropdown_menu(ui, "fc_results_field", shown, ui.available_width(), |ui| {
            if ui.selectable_label(plot.view.contour.is_none(), "None").clicked() {
                run.push(("results.contour", json!({"field": "none"})));
            }
            for f in Field::available(&plot.result) {
                if ui.selectable_label(plot.view.contour == Some(f), f.label()).clicked() {
                    run.push(("results.contour", json!({"field": f.id()})));
                }
            }
        })
    });
    // The legend's range: each end follows the field (Auto) or is clamped, e.g. below a
    // singular peak (#191).
    if let Some(l) = &plot.contour {
        for (label, key, set, value) in
            [("Legend max", "max", plot.view.range_max, l.bands.max), ("Legend min", "min", plot.view.range_min, l.bands.min)]
        {
            kit::row(ui, label, |ui| {
                let mut auto = set.is_none();
                let tip = format!("Follow the field's {}; untick to clamp the legend", if key == "max" { "maximum" } else { "minimum" });
                if kit::checkbox(ui, "Auto", &mut auto).backing_on_hover_text(tip).changed() {
                    run.push(("results.contour", legend_params(plot, key, if auto { json!("auto") } else { json!(value) })));
                }
                // The clamp is stored in SI; the legend shows the display units.
                let mut v = set.map_or(value, |v| l.field.shown(&plot.result, v, &plot.units));
                let speed = (l.max.1 - l.min.1).abs().max(MIN_DRAG_SPAN) * 0.005;
                if ui.add_enabled_ui(set.is_some(), |ui| column_number_field(ui, &mut v, f64::MIN..=f64::MAX, speed, "")).inner.changed() {
                    run.push(("results.contour", legend_params(plot, key, json!(v))));
                }
            });
        }
    }
    if plot.is_fluid() {
        flow_rows(plot, ui, &mut run);
    } else {
        deformation_rows(plot, ui, &mut run);
    }
    kit::row(ui, "Show", |ui| {
        if plot.contour.is_some() {
            kit::checkbox(ui, "Max/Min", &mut app.ui.result_extremes).backing_on_hover_text("Mark the nodes of the largest and smallest value");
        }
        kit::checkbox(ui, "Loads & BCs", &mut app.ui.result_symbols).backing_on_hover_text("The symbols of the step's boundary conditions and loads");
    });
    if let Kind::Explicit(e) = plot.result.kind() {
        kit::divider(ui, "fc_results_explicit");
        explicit(ui, e, plot, t, &mut run);
    } else if let Some((k, at, n)) = plot.frame().filter(|_| plot.is_fluid()) {
        // A fluid step: its frames, the slice, the obstacles and the streamlines following.
        kit::divider(ui, "fc_results_flow");
        kit::row(ui, "Frame", |ui| {
            let mut frame = k;
            if kit::slider(
                ui,
                &mut frame,
                0..=n.saturating_sub(1),
                kit::Slider { integer: true, format: |v| format!("{v:.0}"), ..Default::default() },
            )
            .changed()
            {
                run.push(("results.frame", json!({"frame": frame})));
            }
            let mut animate = plot.view.animate;
            if kit::checkbox(ui, "Animate", &mut animate)
                .backing_on_hover_text("Play the frames from the start-up transient to the settled flow")
                .changed()
            {
                run.push(("results.animate", json!({"on": animate})));
            }
        });
        kit::status_line(ui, &format!("t = {}", with_unit(number(at), unit_of(Quantity::Time, &plot.units).as_deref())), kit::StatusTone::Dim);
        if let Kind::Fluid(f) = plot.result.kind() {
            component_limits(ui, f, &plot.units);
            port_transport(ui, f, plot.pick, &plot.units);
        }
    } else if let Some((k, factor, n)) = plot.frame().filter(|_| matches!(plot.result.kind(), Kind::Static(_))) {
        // A nonlinear static step: its increments, the deformed shape following the slider.
        kit::divider(ui, "fc_results_increments");
        kit::row(ui, "Increment", |ui| {
            let mut frame = k;
            if kit::slider(
                ui,
                &mut frame,
                0..=n.saturating_sub(1),
                kit::Slider { integer: true, format: |v| format!("{v:.0}"), ..Default::default() },
            )
            .changed()
            {
                run.push(("results.frame", json!({"frame": frame})));
            }
        });
        kit::status_line(ui, &format!("load factor {}", number(factor)), kit::StatusTone::Dim);
    } else if let Kind::Frequency(f) = plot.result.kind() {
        kit::divider(ui, "fc_results_modes");
        kit::text(ui, "Modes", kit::Text::Heading("fc_results_modes"));
        egui::ScrollArea::vertical().id_salt("fc_results_modes").max_height(kit::row_height() * 5.0).show(ui, |ui| {
            // Justified: the selected mode highlights across the panel, not only behind its text.
            ui.with_layout(egui::Layout::top_down_justified(egui::Align::LEFT), |ui| {
                for (i, m) in f.modes.iter().enumerate() {
                    let mode = i + 1;
                    if ui
                        .selectable_label(
                            plot.view.mode == mode,
                            format!("{mode}   {}", crate::results::mode_value(f.buckling.is_some(), m.eigenvalue, m.frequency, &plot.units)),
                        )
                        .clicked()
                    {
                        run.push(("results.modes", json!({"mode": mode})));
                    }
                }
            });
        });
    } else if let Kind::Heat(h) = plot.result.kind() {
        kit::divider(ui, "fc_results_heat");
        heat_flows(ui, &fieldcraft_engine::results::heat_flows(h, &plot.units), &plot.units);
    }
    if matches!(plot.result.kind(), Kind::Static(_)) {
        kit::divider(ui, "fc_results_reactions");
        ui.horizontal(|ui| {
            if kit::button(ui, "Reaction sums", kit::Button::Normal)
                .backing_on_hover_text("Reaction force and moment per boundary-condition set (results.reactions)")
                .clicked()
            {
                run.push(("results.reactions", json!({"step": plot.view.step})));
            }
            if kit::button(ui, "Export CSV…", kit::Button::Normal)
                .backing_on_hover_text("Write the reaction sums to a CSV file (results.export.csv)")
                .clicked()
            {
                export = Some(json!({"what": "reactions", "step": plot.view.step}));
            }
        });
        if let Some(r) = app.reactions.as_ref().filter(|r| r["step"] == plot.view.step.as_str()) {
            reactions_table(ui, r, plot);
        }
    }
    let mut xy_export = None;
    if app.xy.is_some() {
        kit::divider(ui, "fc_results_xy");
        xy_export = xy_section(app, plot, ui, t, &mut run);
    }
    if app.path_editor.is_some() {
        kit::divider(ui, "fc_results_path");
        crate::path_panel::controls(app, plot, ui, t);
    }
    kit::divider(ui, "fc_results_probe");
    if let Some(p) = probed(app, plot, ui) {
        export = Some(p);
    }
    if let Some(p) = export {
        app.export_csv(p);
    }
    if let Some(p) = xy_export {
        app.export_csv_with("results.xy", p);
    }
    kit::divider(ui, "fc_results_views");
    let save_view = views_row(app, ui, &mut run);
    for (id, params) in run {
        // The reaction sums land in `app.reactions` (ModelApp::run_form).
        let _ = app.run(id, params);
    }
    if save_view {
        // The form starts with a free name and the viewport's camera.
        let taken: Vec<&str> = app.session.view.views.iter().map(|v| v.name.as_str()).collect();
        let name = (1..).map(|i| format!("View {i}")).find(|n| !taken.contains(&n.as_str())).unwrap_or_default();
        let camera = fieldcraft_engine::view::aim_of(&app.cam, &app.session).json();
        app.open_form("view.save", json!({"name": name, "camera": camera}));
    }
}

/// The saved views: a list that shows one (`view.restore`) and Save View… (`view.save`'s
/// form, the viewport's camera filled in); whether Save View… was clicked.
fn views_row(app: &ModelApp, ui: &mut egui::Ui, run: &mut Vec<(&'static str, Value)>) -> bool {
    let views = &app.session.view.views;
    let current = app.session.view.current.as_deref();
    kit::row(ui, "Views", |ui| {
        let shown = current.unwrap_or(if views.is_empty() { "None saved" } else { "Choose…" });
        ui.add_enabled_ui(!views.is_empty(), |ui| {
            kit::bounded_dropdown_menu(ui, "fc_results_views", shown, ui.available_width(), |ui| {
                for v in views {
                    if ui.selectable_label(current == Some(v.name.as_str()), &v.name).backing_on_hover_text("Show this view (view.restore)").clicked()
                    {
                        run.push(("view.restore", json!({"name": v.name})));
                    }
                }
            });
        });
    });
    kit::button(ui, "Save View…", kit::Button::Normal)
        .backing_on_hover_text("Keep the camera and this results display as a named view, saved with the model (view.save)")
        .clicked()
}

/// A kit number field as wide as what is left of its row's control column, so the rows' fields
/// and drop-downs end at one right edge.
fn column_number_field(ui: &mut egui::Ui, value: &mut f64, range: std::ops::RangeInclusive<f64>, speed: f64, suffix: &str) -> egui::Response {
    ui.spacing_mut().interact_size.x = ui.available_width();
    kit::number_field(ui, value, range, speed, suffix)
}

/// The deformation rows of the results grid: Deformed (with its scale) and Animate.
fn deformation_rows(plot: &Plot, ui: &mut egui::Ui, run: &mut Vec<(&'static str, Value)>) {
    kit::row(ui, "Deformed", |ui| {
        let mut deformed = plot.view.deformed;
        if kit::checkbox(ui, "Show", &mut deformed).backing_on_hover_text("Draw the deformed shape").changed() {
            run.push(("results.deformed", json!({"show": deformed})));
        }
    });
    // Its own row, so Auto lines up with the legend's.
    kit::row(ui, "Scale", |ui| {
        let mut auto = plot.view.scale.is_none();
        let tip = if true_scale(&plot.result) {
            "Draw the displacements at true scale (they are large already)"
        } else {
            "Scale the largest displacement to a tenth of the model"
        };
        if kit::checkbox(ui, "Auto", &mut auto).backing_on_hover_text(tip).changed() {
            run.push(("results.deformed", json!({"scale": if auto { json!("auto") } else { json!(plot.scale) }})));
        }
        // The scale ticking Deformed draws with, also while it is off.
        let mut scale = plot.scale;
        let speed = plot.scale.abs().max(MIN_DRAG_SCALE) * 0.01;
        let field = ui.add_enabled_ui(plot.view.deformed, |ui| {
            column_number_field(ui, &mut scale, 0.0..=fieldcraft_engine::results::MAX_DEFORMATION_SCALE, speed, "×")
        });
        if field.inner.backing_on_hover_text("The scale the displacements are drawn at").changed() {
            run.push(("results.deformed", json!({"scale": scale})));
        }
    });
    kit::row(ui, "Animate", |ui| {
        let mut animate = plot.view.animate;
        if kit::checkbox(ui, "", &mut animate)
            .backing_on_hover_text(
                "Animate shape and contours. Linear static results interpolate a proportional load; nonlinear steps play their converged increments with their own stresses, explicit steps play saved fields.",
            )
            .changed()
        {
            run.push(("results.animate", json!({"on": animate})));
        }
    });
}

/// A fluid step's rows of the results grid: the slice (its normal, position), the obstacles'
/// colours and the streamlines (their count from the inlet).
fn flow_rows(plot: &Plot, ui: &mut egui::Ui, run: &mut Vec<(&'static str, Value)>) {
    use fieldcraft_engine::flow::{Axis, slice_position};
    let (Kind::Fluid(f), Some(flow)) = (plot.result.kind(), plot.flow.as_ref()) else { return };
    if plot.view.contour == Some(Field::WallHeatFlux) {
        kit::row(ui, "Wall heat", |ui| kit::text(ui, "Positive into wall; face average", kit::Text::Cell));
        if plot.contour.is_none() {
            kit::row(ui, "History", |ui| kit::text(ui, "No wall heat recorded at this frame", kit::Text::Cell));
        }
    }
    let view = &flow.view;
    kit::row(ui, "Slice", |ui| {
        if plot.view.contour == Some(Field::WallHeatFlux) {
            kit::text(ui, "Hidden for wall heat", kit::Text::Cell);
            return;
        }
        let shown = view.slice.as_ref().map_or("None", |s| s.axis.id());
        kit::bounded_dropdown_menu(ui, "fc_results_slice", &format!("Normal {}", shown.to_uppercase()), ui.available_width(), |ui| {
            if ui.selectable_label(view.slice.is_none(), "None").clicked() {
                run.push(("results.slice", json!({"axis": "none"})));
            }
            for a in Axis::ALL {
                if ui.selectable_label(view.slice.as_ref().is_some_and(|s| s.axis == a), format!("Normal {}", a.id().to_uppercase())).clicked() {
                    run.push(("results.slice", json!({"axis": a.id()})));
                }
            }
        });
    });
    if let Some(s) = view.slice.as_ref().filter(|_| plot.view.contour != Some(Field::WallHeatFlux)) {
        kit::row(ui, "Position", |ui| {
            let length = |m: f64| shown_in(Quantity::Length, m, &plot.units);
            let (lo, hi) = (length(flow.grid.start(s.axis.index())), length(flow.grid.end(s.axis.index())));
            let mut at = length(slice_position(f, s));
            let speed = (hi - lo).abs().max(MIN_DRAG_SPAN) * 0.002;
            if column_number_field(ui, &mut at, lo..=hi, speed, "").backing_on_hover_text("Where the slice crosses its axis").changed() {
                run.push(("results.slice", json!({"position": at})));
            }
        });
    }
    kit::row(ui, "Obstacles", |ui| {
        let mut surface = view.surface;
        if kit::checkbox(ui, "Coloured", &mut surface).backing_on_hover_text("Colour the part by the selected fluid field or wall heat").changed() {
            run.push(("results.surface", json!({"show": surface})));
        }
    });
    kit::row(ui, "Streamlines", |ui| {
        let mut on = view.streamlines.is_some();
        if kit::checkbox(ui, "", &mut on)
            .backing_on_hover_text("Streamlines from the inlet, coloured by speed or temperature when selected")
            .changed()
        {
            run.push(("results.streamlines", json!({"show": on})));
        }
        if let Some(sl) = &view.streamlines {
            if sl.seeds.is_some() {
                // Given seeds (or a line of them): their number; `count` does not change them.
                kit::text(ui, &format!("{} seeds", flow.seeds), kit::Text::Cell).backing_on_hover_text("The seeds given to results.streamlines");
            } else {
                let mut count = sl.count;
                let field = kit::whole_number_field(ui, &mut count, 1..=fieldcraft_engine::flow::MAX_STREAMLINES, " seeds");
                if field.backing_on_hover_text("Seeds on the inlet face, upstream of the obstacles").changed() {
                    run.push(("results.streamlines", json!({"count": count})));
                }
                if flow.seeds != sl.count as usize {
                    kit::text(ui, &format!("{} traced", flow.seeds), kit::Text::Inline)
                        .backing_on_hover_text("The seeds the grid across the inlet face fits");
                }
            }
        }
    });
}

/// The probed node's values, one row each, or how to probe one; the export's parameters when
/// its Export CSV… is clicked.
fn probed(app: &ModelApp, plot: &Plot, ui: &mut egui::Ui) -> Option<Value> {
    let Some(p) = app.probe.as_ref().filter(|p| p.step == plot.view.step && p.mode == plot.pick.mode && p.frame == plot.pick.frame) else {
        kit::status_line(ui, "Click the model to probe a node", kit::StatusTone::Dim);
        return None;
    };
    kit::text(ui, &format!("Node {}", p.node), kit::Text::Heading("fc_results_probe"));
    for (what, value) in &p.rows {
        kit::row(ui, what, |ui| kit::copy_cell(ui, value, &format!("{what}: {value}"), false, "Copy the value"));
    }
    kit::button(ui, "Export CSV…", kit::Button::Normal)
        .backing_on_hover_text("Write the node's results to a CSV file (results.export.csv)")
        .clicked()
        .then(|| {
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
    let time_unit = unit_of(Quantity::Time, &plot.units);
    if let Some((k, time, n)) = plot.frame() {
        kit::row(ui, "Frame", |ui| {
            let mut frame = k;
            if kit::slider(
                ui,
                &mut frame,
                0..=n.saturating_sub(1),
                kit::Slider { integer: true, format: |v| format!("{v:.0}"), ..Default::default() },
            )
            .changed()
            {
                run.push(("results.frame", json!({"frame": frame})));
            }
        });
        kit::status_line(ui, &format!("t = {}", with_unit(number(time), time_unit.as_deref())), kit::StatusTone::Dim);
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
    let energy_unit = unit_of(Quantity::Energy, &plot.units);
    kit::text(ui, &titled("Energies", energy_unit.as_deref()), kit::Text::Heading("fc_results_energies"));
    let (time, energy) = (|s: f64| shown_in(Quantity::Time, s, &plot.units), |j: f64| shown_in(Quantity::Energy, j, &plot.units));
    let lines: Vec<Line> = series
        .iter()
        .map(|(_, color, f)| Line { color: *color, points: h.time.iter().zip(&h.energies).map(|(t, e)| (time(*t), energy(f(e)))).collect() })
        .collect();
    let t_end = h.time.last().copied().map_or(0.0, time);
    let now = plot.frame().map_or(t_end, |(_, time, _)| time);
    let (response, axes) = line_plot(ui, &lines, "t", time_unit.as_deref(), Some(now), t);
    let scale = axes.lo.abs().max(axes.hi.abs());
    let nearest =
        |at: f64| h.time.iter().enumerate().min_by(|a, b| (time(*a.1) - at).abs().total_cmp(&(time(*b.1) - at).abs())).map_or(0, |(i, _)| i);
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
        text.extend(series.iter().map(|(name, _, f)| format!("{name} {}", number_in(energy(f(e)), scale))));
        response.backing_on_hover_text(text.join("\n"));
    }
    // The values at the frame's time.
    let Some(at_now) = h.energies.get(nearest(now)) else { return };
    let frame = plot.frame().map_or_else(String::new, |(k, ..)| format!(", frame {k}"));
    kit::status_line(ui, &format!("At t = {}{frame}", with_unit(number(now), time_unit.as_deref())), kit::StatusTone::Dim);
    let k = solidcraft_theme::Tokens::get();
    for (name, color, f) in &series {
        // The series' line as its key beside the value, the name in the label colour: line
        // colours are not text colours (the accent, warm and error fall short of 4.5:1 on a dark
        // panel).
        kit::row(ui, name, |ui| {
            let (key, _) = ui.allocate_exact_size(vec2(k.space_14, k.space_10), egui::Sense::hover());
            ui.painter().line_segment([key.left_center(), key.right_center()], Stroke::new(2.0, *color));
            kit::text(ui, &with_unit(number_in(energy(f(at_now)), scale), energy_unit.as_deref()), kit::Text::Cell);
        });
    }
    kit::row(ui, "Balance", |ui| {
        kit::text(ui, &number(at_now.balance), kit::Text::Cell).backing_on_hover_text("The energy balance's relative error")
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
    let k = solidcraft_theme::Tokens::get();
    painter.rect_filled(rect, k.radius_sm, t.field);
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
    let small = FontId::proportional(k.font_micro);
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

/// The quantities the XY plot offers for a result, as (label, `results.xy`'s `y`, its `x` when
/// not the step's own axis); a nodal field's `node` is filled in from the probed node. A
/// nonlinear static step also offers the reaction on each of its boundary-condition sets,
/// against the load factor and against the set's mean travel (a force–travel curve); all of
/// them (the picker scrolls).
pub fn xy_options(s: &fieldcraft_engine::Session, r: &fieldcraft_engine::StepResult) -> Vec<(String, Value, Option<Value>)> {
    use fieldcraft_engine::cmd::xy;
    let mut out = Vec::new();
    if let Kind::Frequency(f) = r.kind() {
        out.push(if f.buckling.is_some() {
            ("Buckling load factor".to_string(), json!("load_factor"), None)
        } else {
            ("Frequency".to_string(), json!("frequency"), None)
        });
    }
    if let Kind::Explicit(e) = r.kind() {
        out.extend(xy::ENERGIES.iter().map(|n| (format!("Energy: {}", n.replace('_', " ")), json!({"energy": n}), None)));
        out.extend(xy::HISTORIES.iter().map(|n| (n.replace('_', " "), json!({"history": n}), None)));
        for node in e.history.nodes.iter().take(8) {
            out.extend(
                xy::NODE_QUANTITIES
                    .iter()
                    .map(|q| (format!("History node {node}: {}", q.replace('_', " ")), json!({"history_node": node, "quantity": q}), None)),
            );
        }
    }
    let linear = matches!(r.kind(), Kind::Static(st) if st.nonlinear.is_none());
    if !linear {
        out.extend(xy::plottable(r).into_iter().map(|f| (format!("{} at a node", f.label()), json!({"field": f.id()}), None)));
    }
    if matches!(r.kind(), Kind::Static(st) if st.nonlinear.is_some()) {
        let mut sets: Vec<&str> = s
            .model
            .step(r.step())
            .map(|st| st.boundary_conditions.iter().filter_map(|bc| s.model.boundary_condition(bc)?.target().set()).collect())
            .unwrap_or_default();
        sets.sort_unstable();
        sets.dedup();
        for set in sets {
            for c in ["x", "y", "z"] {
                let upper = c.to_uppercase();
                let y = json!({"reaction": set, "component": c});
                out.push((format!("Force–travel {upper}: {set}"), y.clone(), Some(json!({"displacement": set, "component": c}))));
                out.push((format!("Reaction {upper} on {set}"), y, None));
            }
        }
    }
    out
}

/// The XY plot (results.xy): a quantity picker, the series plotted against the step's own axis
/// with the frame (mode) shown marked (a click shows the nearest one), and its export's
/// parameters when Export CSV… is clicked. A nodal field is read at the probed node.
fn xy_section(app: &mut ModelApp, plot: &Plot, ui: &mut egui::Ui, t: &Tokens, run: &mut Vec<(&'static str, Value)>) -> Option<Value> {
    let step = plot.view.step.clone();
    let mut close = false;
    kit::panel::header(ui, "results.xy", "XY Plot", || close = true);
    if close {
        app.xy = None;
        return None;
    }
    let options = xy_options(&app.session, &plot.result);
    if options.is_empty() {
        kit::status_line(ui, "A linear static step has one state: plot a nonlinear static, explicit, heat or frequency step", kit::StatusTone::Dim);
        return None;
    }
    let probed = app.probe.as_ref().filter(|p| p.step == step).map(|p| p.node);
    let xy = app.xy.get_or_insert_with(Default::default);
    if probed.is_some() {
        xy.node = probed;
    }
    let current = (xy.y.clone(), xy.x.clone());
    let template = options
        .iter()
        .find(|(_, y, x)| current.0.as_ref() == Some(y) && current.1 == *x)
        .or(options.first())
        .map(|(_, y, x)| (y.clone(), x.clone()))
        .unwrap_or_default();
    let shown = options.iter().find(|(_, y, x)| (y, x) == (&template.0, &template.1)).map_or("", |(l, ..)| l.as_str());
    let (mut picked, mut against) = template.clone();
    kit::row(ui, "Quantity", |ui| {
        kit::bounded_dropdown_menu(ui, "fc_results_xy", shown, ui.available_width(), |ui| {
            for (label, y, x) in &options {
                if ui.selectable_label((y, x) == (&template.0, &template.1), label).clicked() {
                    (picked, against) = (y.clone(), x.clone());
                }
            }
        })
    });
    (xy.y, xy.x) = (Some(picked.clone()), against.clone());
    let y = if picked.get("field").is_some() {
        let Some(node) = xy.node else {
            (xy.key, xy.series) = (None, None);
            kit::status_line(ui, "Click the model to pick the node", kit::StatusTone::Dim);
            return None;
        };
        let mut y = picked;
        y["node"] = json!(node);
        y
    } else {
        picked
    };
    let mut params = json!({"step": step, "y": y});
    if let Some(x) = against {
        params["x"] = x;
    }
    let key = (step.clone(), Arc::as_ptr(&plot.result) as usize, params.clone());
    if xy.key.as_ref() != Some(&key) {
        xy.key = Some(key);
        xy.series = None;
        run.push(("results.xy", params));
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
    kit::status_line(ui, s["y"]["label"].as_str().unwrap_or_default(), kit::StatusTone::Dim);
    let x_name = match axis {
        "load_factor" => "λ",
        "time" => "t",
        "displacement" => "u",
        _ => "mode",
    };
    let x_unit = match axis {
        "time" => unit_of(Quantity::Time, &plot.units),
        "displacement" => unit_of(Quantity::Length, &plot.units),
        _ => None,
    };
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
        response.backing_on_hover_text(format!("{x_label} {}\n{}", number(x), number_in(y, scale)));
    }
    if let Some(i) = current {
        let at = |v: &[f64]| v.get(i).map_or("—".into(), |x| number(*x));
        kit::status_line(ui, &format!("{x_label} {}: {}", at(&xs), at(&ys)), kit::StatusTone::Dim);
    }
    kit::button(ui, "Export CSV…", kit::Button::Normal)
        .backing_on_hover_text("Write the series to a CSV file (results.xy)")
        .clicked()
        .then_some(params)
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
fn heat_flows(ui: &mut egui::Ui, flows: &Value, units: &fieldcraft_model::DisplaySystem) {
    kit::text(ui, &titled("Heat flow into the body", unit_of(Quantity::Power, units).as_deref()), kit::Text::Heading("fc_results_heat"));
    let row = |ui: &mut egui::Ui, name: &str, v: Option<f64>| {
        kit::row(ui, name, |ui| kit::text(ui, &v.map_or("—".into(), number), kit::Text::Cell));
    };
    for l in flows["loads"].as_array().into_iter().flatten() {
        row(ui, &format!("{} ({})", l["name"].as_str().unwrap_or("?"), l["type"].as_str().unwrap_or("").replace('_', " ")), l["heat"].as_f64());
    }
    row(ui, "Prescribed temperatures", flows["reactions"].as_f64());
    row(ui, "Imbalance", flows["imbalance"].as_f64());
}

/// Directional open-face transport from the displayed fluid frame.
fn port_transport(ui: &mut egui::Ui, f: &FluidResult, pick: fieldcraft_engine::results::Pick, units: &fieldcraft_model::DisplaySystem) {
    let Ok(frame) = fieldcraft_engine::flow::frame_index(f, pick) else { return };
    if !f.open_faces.contains(&true) {
        return;
    }
    kit::section(ui, "fc_inlets_outlets", "Inlets and outlets", false, |ui| {
        let report = fieldcraft_engine::results::fluid_ports(f, frame, units);
        let Some(faces) = report["faces"].as_array() else {
            kit::status_line(ui, report["error"].as_str().unwrap_or("Port fields unavailable"), kit::StatusTone::Error);
            return;
        };
        kit::status_line(ui, "From the displayed frame", kit::StatusTone::Dim);
        let format =
            |v: &Value, quantity: &str| v.as_f64().map(|v| with_unit(number(v), report["units"][quantity].as_str())).unwrap_or_else(|| "—".into());
        for face in faces {
            let name = face["face"].as_str().unwrap_or("?");
            kit::text(ui, name, kit::Text::Heading("fc_port"));
            let rows: Vec<Vec<String>> = [
                ("Volume flow", "volume_flow", "volume_flow"),
                ("Mass flow", "mass_flow", "mass_flow"),
                ("Mixed temperature", "temperature", "temperature"),
                ("Static pressure", "pressure", "pressure"),
                ("Sensible heat", "sensible_heat", "power"),
            ]
            .into_iter()
            .map(|(label, key, quantity)| {
                let amount = |v: &Value| v.as_f64().map(number).unwrap_or_else(|| "—".into());
                vec![titled(label, report["units"][quantity].as_str()), amount(&face["inflow"][key]), amount(&face["outflow"][key])]
            })
            .collect();
            crate::reply::text_table(ui, &format!("fc_port_{name}"), &["", "In", "Out"], &rows);
        }
        kit::row(ui, "Net mass out", |ui| kit::text(ui, &format(&report["net_mass_flow"], "mass_flow"), kit::Text::Cell));
        kit::row(ui, "Net heat out", |ui| kit::text(ui, &format(&report["net_sensible_heat"], "power"), kit::Text::Cell));
        kit::status_line(
            ui,
            "Mass uses reference density. Heat is advective transport, relative to the mixed inflow temperature below; wall heat, conduction and storage are separate.",
            kit::StatusTone::Dim,
        );
        kit::row(ui, "Mixed inflow temperature", |ui| kit::text(ui, &format(&report["heat_reference"], "temperature"), kit::Text::Cell));
    });
}

/// The conducting solids in the order the limits list shows them: those over their limit
/// first, then those under it, each from the smallest margin up, then those without a limit.
fn by_margin(checks: &[ComponentCheck]) -> Vec<&ComponentCheck> {
    let mut rows: Vec<_> = checks.iter().collect();
    rows.sort_by(|a, b| {
        let rank = |c: &ComponentCheck| match c.ok {
            Some(false) => 0,
            Some(true) => 1,
            None => 2,
        };
        rank(a).cmp(&rank(b)).then(a.margin.unwrap_or(f64::INFINITY).total_cmp(&b.margin.unwrap_or(f64::INFINITY)))
    });
    rows
}

/// A fluid step's conducting solids against their temperature limits: the verdict, then per
/// solid its peak, limit and margin, those over their limit first and marked.
fn component_limits(ui: &mut egui::Ui, f: &FluidResult, units: &fieldcraft_model::DisplaySystem) {
    let checks = f.components();
    if checks.iter().all(|c| c.peak.is_none()) {
        return;
    }
    kit::divider(ui, "fc_results_components");
    let unit = unit_of(Quantity::Temperature, units);
    kit::text(ui, &titled("Component temperatures", unit.as_deref()), kit::Text::Heading("fc_results_components"));
    let over = checks.iter().filter(|c| c.ok == Some(false)).count();
    let limited = checks.iter().filter(|c| c.limit.is_some()).count();
    match f.components_ok() {
        Some(true) => {
            kit::check_status_line(
                ui,
                "fc_results_components:verdict",
                &format!("Every component stays under its limit ({limited} checked)"),
                kit::CheckStatus::Ok,
            );
        }
        Some(false) => {
            kit::check_status_line(
                ui,
                "fc_results_components:verdict",
                &format!("{over} of {limited} components over their limit"),
                kit::CheckStatus::Error,
            );
        }
        None if limited == 0 => {
            kit::status_line(ui, "No temperature limits: set one per conducting solid", kit::StatusTone::Dim);
        }
        None => {}
    }
    let temperature = |k: f64| number(shown_in(Quantity::Temperature, k, units));
    let difference = unit_of(Quantity::TemperatureDifference, units);
    let margin = titled("Margin", difference.as_deref());
    let columns = [
        kit::Column::new("", 1.6).flex(),
        kit::Column::new("Peak", 1.0).flex(),
        kit::Column::new("Limit", 1.0).flex(),
        kit::Column::new(&margin, 1.0).flex(),
    ];
    kit::table_custom(ui, "fc_results_components", &columns, "", |ui, table| {
        table.headers(ui, |_, _| false);
        for c in by_margin(&checks) {
            let tip = c.peak.map(|p| {
                let at = with_unit(number(shown_in(Quantity::Time, p.time, units)), unit_of(Quantity::Time, units).as_deref());
                if p.steady { format!("Steady value (t = {at})") } else { format!("Highest over the frames, at t = {at}") }
            });
            let over = c.ok == Some(false);
            table.row(ui, &c.body, false, |row| {
                row.cell(|ui| {
                    // Over its limit: the error mark and the word in the margin say so, not only a colour.
                    let name = match c.ok {
                        Some(ok) => kit::check_status_line(
                            ui,
                            &format!("fc_results_components:{}", c.body),
                            &c.body,
                            if ok { kit::CheckStatus::Ok } else { kit::CheckStatus::Error },
                        ),
                        None => kit::text(ui, &c.body, kit::Text::Cell),
                    };
                    if let Some(tip) = &tip {
                        name.backing_on_hover_text(tip);
                    }
                });
                row.cell(|ui| kit::text(ui, &c.peak.map_or("—".into(), |p| temperature(p.temperature)), kit::Text::Cell));
                row.cell(|ui| kit::text(ui, &c.limit.map_or("—".into(), temperature), kit::Text::Cell));
                let margin = c.margin.map_or("—".into(), |m| number(shown_in(Quantity::TemperatureDifference, m, units)));
                row.cell(|ui| kit::text(ui, &margin, kit::Text::Value("fc_results_components:margin", over)));
            });
        }
    });
}

/// The reaction sums: per set and over the model, the force and the moment about the point, and
/// the applied total beside the reacted one (they balance in equilibrium). Round-off below a
/// billionth of the largest force (moment) is 0.
fn reactions_table(ui: &mut egui::Ui, r: &Value, plot: &Plot) {
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
    let (force, moment) = (unit_of(Quantity::Force, &plot.units), unit_of(Quantity::Moment, &plot.units));
    let about = vector(&r["about"]).unwrap_or_default();
    kit::status_line(ui, &format!("Moments about ({})", components(&about, 0.0)), kit::StatusTone::Dim);
    let cells = |v: Option<[f64; 3]>, scale: f64| -> Vec<String> {
        (0..3).map(|k| v.and_then(|v| v.get(k).copied()).map_or("—".into(), |x| number_in(x, scale))).collect()
    };
    let (f_label, m_label) = (titled("F", force.as_deref()), titled("M", moment.as_deref()));
    for (name, f, m) in &rows {
        kit::text(ui, name, kit::Text::Heading("fc_results_reactions"));
        let table = [[vec![f_label.clone()], cells(*f, force_scale)].concat(), [vec![m_label.clone()], cells(*m, moment_scale)].concat()];
        crate::reply::text_table(ui, &format!("fc_results_reactions_{name}"), &["", "X", "Y", "Z"], &table);
    }
    if applied.is_some() {
        let table = [[vec![titled("Applied F", force.as_deref())], cells(applied, force_scale)].concat()];
        crate::reply::text_table(ui, "fc_results_reactions_applied", &["", "X", "Y", "Z"], &table);
    }
    // Equilibrium: the reactions over the whole model balance the applied loads.
    if let (Some(total), Some(applied)) = (rows.last().and_then(|r| r.1), applied) {
        let residual = (0..3).map(|k| (total[k] + applied[k]).abs()).fold(0.0, f64::max);
        let balanced = residual <= fieldcraft_engine::results::BALANCE_TOLERANCE * force_scale.max(f64::MIN_POSITIVE);
        if balanced {
            kit::check_status_line(ui, "fc_results_reactions:balance", "Reactions balance the applied loads", kit::CheckStatus::Ok);
        } else {
            let text = format!("Reactions and applied loads differ by {}", with_unit(number(residual), force.as_deref()));
            kit::check_status_line(ui, "fc_results_reactions:balance", &text, kit::CheckStatus::Warning);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::results;

    /// #205: nonlinear playback colours each increment by its own stress; a frame without
    /// its own stress (no stress output, or past the step's field cap) says so in the legend.
    #[test]
    fn nonlinear_playback_shows_increment_stress() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("model.section.edit", json!({"name": "Beam", "set": "beam", "material": "Steel", "kind": {"type": "solid"}, "formulation": {}}))
            .unwrap();
        app.run("model.step.create", json!({"name": "Increments", "procedure": {"type": "nonlinear_static", "nonlinear_geometry": false, "initial_increment": 0.25}, "boundary_conditions": ["Clamp"], "loads": ["Push"]})).unwrap();
        app.run("solve.nonlinear", json!({"step": "Increments"})).unwrap();
        app.run("results.contour", json!({"field": "von_mises"})).unwrap();
        let final_max = app.plot().unwrap().contour.as_ref().unwrap().max.1;
        app.run("results.animate", json!({})).unwrap();
        app.now = 0.0;
        let plot = app.plot().unwrap();
        assert!(!plot.unavailable && plot.contour.is_some() && plot.amplitude.is_none());
        assert_eq!(title(&plot), "Von Mises stress [MPa]");
        app.run("results.deformed", json!({"show": false})).unwrap();
        app.now = results::PERIOD * (results::FRAMES - 1) as f64 / results::FRAMES as f64;
        let final_frame = app.plot().unwrap();
        assert_ne!(plot.pick.frame, final_frame.pick.frame, "recorded nonlinear increments keep playing without deformation");
        assert!(final_frame.contour.is_some() && !final_frame.unavailable && final_frame.amplitude.is_none());
        app.run("results.animate", json!({"on": false})).unwrap();
        app.run("results.frame", json!({"frame": 1})).unwrap();
        let first = app.plot().unwrap();
        let max = first.contour.as_ref().unwrap().max.1;
        assert!(!first.unavailable && max > 0.0 && max < 0.5 * final_max, "a quarter of the load: {max} of {final_max}");
        let mut capped = first.as_ref().clone();
        capped.unavailable = true;
        let ctx = egui::Context::default();
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            let legend =
                legend_layout(&capped, &ui.ctx().debug_painter(), Rect::from_min_size(Pos2::ZERO, vec2(1280.0, 800.0)), None, &Tokens::dark());
            assert!(legend.lines.iter().any(|(line, _)| line == "Field not saved at this increment"));
        });
        out.textures_delta.clear();
    }

    /// The clickable extreme rows keep full-sized text inside the legend, even with a stale
    /// warning and the extra over-limit swatch at the native minimum window size (#610).
    #[test]
    fn short_legends_keep_both_extreme_targets_readable() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.contour", json!({"field": "von_mises", "bands": 24, "range": {"min": 0, "max": 10}})).unwrap();
        let mut plot = app.plot().unwrap().as_ref().clone();
        let ctx = egui::Context::default();
        let t = Tokens::get(true);
        let view = Rect::from_min_size(pos2(250.0, 126.0), vec2(401.0, 406.0));
        let bar = Some(Rect::from_min_size(pos2(347.0, 136.0), vec2(207.0, 31.0)));
        for stale in [false, true] {
            plot.stale = stale;
            let mut output = ctx.run_ui(Default::default(), |ui| {
                let painter = ui.painter_at(view);
                let layout = legend_layout(&plot, &painter, view, bar, &t);
                let targets = legend(&plot, &painter, view, bar, &t);
                assert_eq!(targets.len(), 2);
                assert!(layout.band_height < SWATCH.y);
                for (target, node) in targets {
                    assert!(layout.frame.contains_rect(target), "{target:?} outside {:?}", layout.frame);
                    assert!(target.height() >= 12.0, "extremes retain their readable font");
                    assert!(!target.intersects(crate::viewport::nav_bar_rect(view)));
                    focus_node(&mut app, &plot, node);
                    assert_eq!(app.cam.target, plot.extreme_position(node).unwrap());
                    assert!(app.probe.is_some());
                }
            });
            output.textures_delta.clear();
            let mut ticks: Vec<_> = output
                .shapes
                .iter()
                .filter_map(|s| match &s.shape {
                    egui::Shape::Text(text) => text.galley.text().parse::<f64>().ok().map(|value| (text.pos.y, text.galley.size().y, value)),
                    _ => None,
                })
                .collect();
            ticks.sort_by(|a, b| a.0.total_cmp(&b.0));
            assert!(ticks.iter().any(|t| t.2 == 10.0), "the clamped maximum stays labelled");
            assert!(ticks.first().is_some_and(|t| t.2 > 150.0), "the true maximum stays labelled");
            assert!(ticks.last().is_some_and(|t| t.2 == 0.0), "the lower end stays labelled");
            for pair in ticks.windows(2) {
                assert!(pair[0].0 + pair[0].1 <= pair[1].0, "numeric ticks overlap: {pair:?}");
            }
        }
    }

    /// A label goes right of its point, left when right runs off, inside the top and bottom;
    /// one that would cover a placed label goes under it, or over it at the bottom edge.
    /// #590: 65 frames (0 to 64) end at "frame 64 of 65", not "64 of 64"; four increments after
    /// the unloaded state end at "increment 4 of 4".
    #[test]
    fn the_legend_counts_frames_from_zero_of_all_frames() {
        assert!(frame_text(false, 64, 65, 16.0, Some("s")).starts_with(", frame 64 of 65: t = 16"));
        assert!(frame_text(false, 0, 65, 0.0, Some("s")).starts_with(", frame 0 of 65: t = "));
        assert!(frame_text(true, 4, 5, 1.0, None).starts_with(", increment 4 of 4: load factor "));
    }

    /// #672: the components over their limit come first, the worst on top, then those under it
    /// from the smallest margin, then those without a limit.
    #[test]
    fn components_over_their_limit_are_listed_first() {
        let check = |body: &str, margin: Option<f64>| ComponentCheck {
            body: body.into(),
            peak: Some(fieldcraft_engine::results::Peak { temperature: 350.0, time: 1.0, steady: true }),
            limit: margin.map(|_| 360.0),
            margin,
            ok: margin.map(|m| m >= 0.0),
        };
        let checks = [check("Sink", None), check("Vrm", Some(12.0)), check("Die", Some(-3.0)), check("Ram", Some(2.0)), check("Gpu", Some(-8.0))];
        let order: Vec<&str> = by_margin(&checks).iter().map(|c| c.body.as_str()).collect();
        assert_eq!(order, ["Gpu", "Die", "Ram", "Vrm", "Sink"]);
    }

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

    /// A3-3: a marker under the legend is not drawn over it, and a label that would cover the
    /// legend goes under (or over) it instead; with no free place the label is left out.
    #[test]
    fn markers_and_their_labels_keep_off_the_legend() {
        let bounds = Rect::from_min_max(pos2(0.0, 0.0), pos2(300.0, 400.0));
        let legend = Rect::from_min_max(pos2(12.0, 12.0), pos2(200.0, 300.0));
        let size = vec2(80.0, 20.0);
        assert_eq!(clear_label(pos2(100.0, 100.0), size, bounds, Some(legend), None), None, "under the legend");
        assert_eq!(clear_label(pos2(220.0, 100.0), size, bounds, None, None), Some(label_box(pos2(220.0, 100.0), size, bounds, None)));
        // Near the right edge the label goes left, over the legend: it moves under the legend.
        let r = clear_label(pos2(210.0, 150.0), size, bounds, Some(legend), None).expect("a free place");
        assert!(!r.intersects(legend) && r.top() >= legend.bottom(), "{r:?}");
        // A second label keeps off the first and off the legend, or is left out.
        if let Some(second) = clear_label(pos2(212.0, 152.0), size, bounds, Some(legend), Some(r)) {
            assert!(!second.intersects(legend) && !second.intersects(r), "{second:?}");
        }
    }

    /// A3-14: the Results rows share one control column: the drop-downs and the number fields
    /// end at the same right edge, and the deformed shape has its own labelled toggle.
    #[test]
    fn the_results_rows_share_one_control_column() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.contour", json!({"field": "von_mises", "range": {"max": 10}})).unwrap();
        let plot = app.plot().unwrap();
        let ctx = egui::Context::default();
        let t = Tokens::get(true);
        for _ in 0..2 {
            let mut out = ctx.run_ui(Default::default(), |ui| {
                ui.allocate_ui(vec2(320.0, 2000.0), |ui| controls(&mut app, &plot, ui, &t));
            });
            out.textures_delta.clear();
        }
        let right = |row: &str| solidcraft_ui_registry::rect_of(&format!("dialog:widget:{row}")).map(|r| r.right());
        let step = right("Step").expect("the Step drop-down");
        for row in ["Contours", "Legend max", "Legend min", "Scale", "Views"] {
            let edge = right(row).unwrap_or_else(|| panic!("no {row} control"));
            assert!((edge - step).abs() <= 1.0, "{row} ends at {edge}, the Step drop-down at {step}");
        }
        assert!(right("Deformed").is_some(), "the deformed shape has its own toggle");
    }

    /// A3-14: the legend's values are right-aligned, so their decimal points line up.
    #[test]
    fn legend_values_are_right_aligned() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.contour", json!({"field": "von_mises", "bands": 5, "range": {"min": 0, "max": 100}})).unwrap();
        let plot = app.plot().unwrap();
        let ctx = egui::Context::default();
        let t = Tokens::get(true);
        let view = Rect::from_min_size(Pos2::ZERO, vec2(1280.0, 800.0));
        let mut out = ctx.run_ui(Default::default(), |ui| {
            legend(&plot, &ui.painter_at(view), view, None, &t);
        });
        out.textures_delta.clear();
        let rights: Vec<f32> = out
            .shapes
            .iter()
            .filter_map(|s| match &s.shape {
                egui::Shape::Text(text) if text.galley.text().replace(',', "").parse::<f64>().is_ok() => Some(text.pos.x + text.galley.size().x),
                _ => None,
            })
            .collect();
        assert!(rights.len() >= 3, "{rights:?}");
        for r in &rights {
            assert!((r - rights[0]).abs() <= 0.5, "the values end at different x: {rights:?}");
        }
    }

    /// The results controls, Path Plot and command replies draw with the kit's rows, tables,
    /// number fields, dividers, buttons and status marks, and their painted text with the theme's
    /// font sizes: none of egui's raw grids, drag values, separators, buttons, links or
    /// checkboxes, no literal font size and no emoji for a severity (UI audit K4, K5, K7, K8).
    #[test]
    fn results_path_and_replies_draw_with_the_kit() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let raw = [
            concat!("egui::", "Grid"),
            concat!("egui::", "DragValue"),
            concat!("ui.", "separator()"),
            concat!("ui.", "button("),
            concat!(".", "button(\""),
            concat!("ui.", "link("),
            concat!("egui::", "Link"),
            concat!("ui.", "checkbox("),
            concat!("ui.", "selectable_value("),
            concat!("RichText::new(\"", "Path Plot\")"),
            concat!("RichText::new(\"", "XY Plot\")"),
            "\u{26d4}",
            "\u{26a0}",
            "\u{2714}",
        ];
        let mut bad = Vec::new();
        for file in ["results_panel.rs", "reply.rs", "path_panel.rs", "results.rs"] {
            let text = std::fs::read_to_string(dir.join(file)).unwrap();
            // The source before its tests: the tests name what they look for.
            let code = text.split("#[cfg(test)]").next().unwrap_or_default();
            for (n, line) in code.lines().enumerate() {
                for r in raw {
                    if line.contains(r) {
                        bad.push(format!("{file}:{}: {r}", n + 1));
                    }
                }
                // A painted font's size is a theme token, never a number.
                if let Some(at) = line.find(concat!("FontId::", "proportional(")) {
                    let arg = &line[at + concat!("FontId::", "proportional(").len()..];
                    if arg.starts_with(|c: char| c.is_ascii_digit()) {
                        bad.push(format!("{file}:{}: literal font size", n + 1));
                    }
                }
            }
        }
        assert!(bad.is_empty(), "use the kit's rows, tables, fields, dividers, buttons and theme tokens: {bad:#?}");
    }
}
