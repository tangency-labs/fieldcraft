//! The 3D viewport: the model drawn on the GPU (or on the CPU without one), camera navigation,
//! the navigation bar (bottom centre) and the orientation gizmo (top right), as in SolidCraft.
//!
//! Navigation ([`crate::navigation`] presets, chosen in the display menu): by default left drag
//! orbits (or pans or zooms: the navigation bar's modes), right or middle drag (or Shift + left
//! drag) pans, Ctrl + left drag from empty space orbits. An orbit turns about the model point
//! under the press (marked while orbiting), or the centre of what is drawn from empty space. The
//! wheel or a pinch zooms toward the model point under the cursor; a double click or F frames
//! the model. The gizmo's axis tips turn the view to look along them. The bar and the gizmo
//! follow SolidCraft's `viewport.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.
//!
//! Selection ([`crate::selection`]): what a click would pick is highlighted under the pointer; a
//! click selects it (Shift adds, Cmd/Ctrl adds or removes, a click on nothing clears), Escape
//! clears. The selection bar (top centre) chooses what is picked and turns on box selection,
//! where a left drag selects the visible nodes, element faces or elements inside the box. A
//! right click offers a set made from the selection, and on a datum or a partition's cut face
//! that item's menu ([`crate::menu`]), as its browser row does.

use egui::{Color32, Pos2, Rect, Sense, Stroke, Vec2, pos2, vec2};
use fieldcraft_render::{StandardView, Vec3};

use crate::browser::Def;
use crate::gpu::{Frame, GpuScene, ViewportCallback};
use crate::icons::{self, Look};
use crate::menu::{self, Target};
use crate::navigation::{self, NavButton};
use crate::selection::{Filter, Mode};
use crate::theme::Tokens;
use crate::{ModelApp, NavMode, Style};

/// Wheel zoom per point scrolled.
const ZOOM_PER_POINT: f64 = 0.0015;
/// The smallest pinch factor a zoom divides by.
const MIN_PINCH: f64 = 1e-3; // unit-factors: a zoom ratio
/// A zoom factor this close to 1 is no zoom.
const NO_ZOOM: f64 = 1e-9; // unit-factors: a zoom ratio

pub fn show(app: &mut ModelApp, ui: &mut egui::Ui) {
    let rect = ui.available_rect_before_wrap();
    app.viewport.rect = Some(rect);
    app.place("viewport", rect);
    let t = Tokens::get(app.ui.dark);
    app.viewport.overlays = overlays(app, ui, rect, &t);
    app.apply_reframe();
    let response = ui.allocate_rect(rect, Sense::click_and_drag());
    let boxing = app.pick.box_mode && app.pick.filter.boxes();
    node_drag(app, ui, &response, boxing);
    navigate(app, ui, &response, rect, boxing);
    if boxing && app.viewport.drag.is_none() {
        box_drag(app, ui, &response);
    }
    if response.clicked()
        && let Some(at) = response.interact_pointer_pos()
    {
        // With Relax Nodes open a plain click on the mesh places its sphere (#398); Shift or
        // Ctrl picks faces to relax (#471), and so does every click once the form lists faces.
        let plain = mode(ui.input(|i| i.modifiers)) == Mode::Replace;
        let placed =
            app.path_pick_at(at) || app.relaxing() && plain && !app.relax_by_faces() && app.surface_point(at).is_some_and(|p| app.relax_place(p));
        if !placed {
            let _ = app.select_at(at, mode(ui.input(|i| i.modifiers)));
            probe_at(app, at, rect);
        }
    }
    if !app.capturing {
        hover(app, &response);
    }
    context_menu(app, &response);
    let painter = ui.painter_at(rect);
    gradient(&painter, rect, t.viewport_top, t.viewport_bottom);
    if app.viewport.gpu.is_some() {
        paint_gpu(app, ui, &painter, rect);
        labels(app, &painter, &t);
    } else {
        paint_cpu(app, ui, &painter, rect);
    }
    sampled_note(app, &painter, rect, &t);
    hover_label(app, &painter, &t);
    pivot_marker(app, &painter, &t);
    if let Some((a, b)) = app.pick.dragging {
        let r = Rect::from_two_pos(a, b);
        painter.rect(r, 0.0, t.pick_selected.gamma_multiply(0.12), Stroke::new(1.0, t.pick_selected), egui::StrokeKind::Inside);
    }
    if app.display().bounds.is_empty() && app.busy().is_none() {
        start_screen(app, ui, rect, &t);
    } else {
        selection_bar(app, ui, rect, &t);
    }
    crate::results_panel::show(app, ui, rect, &t);
    crate::path_panel::overlay(app, &painter, rect, &t);
    axis_view(app, ui, rect, &t);
    nav_bar(app, ui, rect, &t);
}

/// Where this frame's overlays go over the viewport `rect`: the results legend, the selection
/// bar (where it was drawn last), the orientation gizmo and the navigation bar.
fn overlays(app: &mut ModelApp, ui: &egui::Ui, rect: Rect, t: &Tokens) -> Vec<Rect> {
    let mut out: Vec<Rect> = [axis_view_rect(rect), nav_bar_rect(rect)].into_iter().filter(|r| rect.contains_rect(*r)).collect();
    if let Some(plot) = app.plot() {
        let bar = selection_bar_rect(ui.ctx());
        let legend = crate::results_panel::legend_rect(&plot, &ui.painter_at(rect), rect, bar, t);
        app.place("legend", legend);
        out.push(legend);
    }
    out.extend(selection_bar_rect(ui.ctx()));
    out
}

/// The Start card's largest size; a smaller viewport shrinks it and its samples scroll.
const START_W: f32 = 820.0;
const START_H: f32 = 620.0;
/// The margins the view cube and the navigation bar keep.
const START_RIGHT: f32 = 130.0;
const START_BOTTOM: f32 = 56.0;

/// The empty viewport: open a file or a sample (the samples are ready to solve).
fn start_screen(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect, t: &Tokens) {
    let mut open = None;
    let mut pick = false;
    // Clear of the view cube (top right) and the navigation bar (bottom centre).
    let free = Rect::from_min_max(rect.min + vec2(12.0, 12.0), rect.max - vec2(START_RIGHT, START_BOTTOM));
    let area = Rect::from_center_size(free.center(), vec2(START_W.min(free.width()), START_H.min(free.height())));
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(area).layout(egui::Layout::top_down(egui::Align::Min)));
    egui::Frame::NONE.fill(t.overlay).corner_radius(6.0).inner_margin(egui::Margin::same(14)).show(&mut child, |ui| {
        ui.label(egui::RichText::new("Start").size(16.0).color(t.overlay_text));
        ui.label(egui::RichText::new("Open a model, a SolidCraft design or a tagged mesh — or drop one here.").color(t.text_dim));
        ui.add_space(4.0);
        pick = ui.button("Open…").clicked();
        ui.add_space(8.0);
        ui.label(egui::RichText::new("Samples").color(t.overlay_text));
        // The grid takes what the card has left and scrolls the rest.
        let left = (area.bottom() - ui.cursor().top() - 14.0).max(60.0);
        crate::samples::grid(ui, t, left, &mut open);
    });
    if pick {
        app.open_dialog(None);
    } else if let Some(n) = open {
        app.open_sample(n);
    }
}

/// A click on a results plot probes the node nearest to the face under the pointer.
fn probe_at(app: &mut ModelApp, at: Pos2, rect: Rect) {
    let Some(plot) = app.plot() else { return };
    let (origin, dir) = app.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
    if let Some((node, _)) = plot.display.pick(origin, dir) {
        // A failure (a stale result) shows in the status bar.
        let _ = app.probe_node(node);
    }
}

/// How a click with these modifier keys changes the selection: Cmd/Ctrl adds or removes, Shift
/// adds, a plain click replaces.
fn mode(m: egui::Modifiers) -> Mode {
    if m.command || m.ctrl {
        Mode::Toggle
    } else if m.shift {
        Mode::Add
    } else {
        Mode::Replace
    }
}

/// Box selection: a left drag draws the box; on release, what is visible inside it is selected
/// (with the release's modifier keys).
fn box_drag(app: &mut ModelApp, ui: &egui::Ui, response: &egui::Response) {
    let pointer = response.interact_pointer_pos();
    if response.drag_started_by(egui::PointerButton::Primary) {
        app.pick.dragging = pointer.map(|p| (p, p));
    } else if let (Some((a, _)), Some(p)) = (app.pick.dragging, pointer)
        && response.dragged_by(egui::PointerButton::Primary)
    {
        app.pick.dragging = Some((a, p));
    }
    if response.drag_stopped_by(egui::PointerButton::Primary)
        && let Some((a, b)) = app.pick.dragging.take()
        && a.distance(b) > 3.0
    {
        app.select_box(a, b, mode(ui.input(|i| i.modifiers)));
    }
}

/// The Move Nodes tool: with its form open, a plain left drag that starts on the mesh moves the
/// node under the press in the screen plane ([`crate::node_move`]); release runs the move.
fn node_drag(app: &mut ModelApp, ui: &egui::Ui, response: &egui::Response, boxing: bool) {
    if response.drag_started_by(egui::PointerButton::Primary)
        && app.moving_nodes()
        && !boxing
        && mode(ui.input(|i| i.modifiers)) == Mode::Replace
        && let Some(at) = ui.input(|i| i.pointer.press_origin()).or(response.interact_pointer_pos())
    {
        app.node_move_grab(at);
    }
    // A changed form is prepared once it settles (here, before the overlay is drawn): draw
    // again then.
    if app.moving_nodes() {
        app.node_move_preview();
    }
    if app.node_move_pending() {
        ui.ctx().request_repaint_after(std::time::Duration::from_secs_f64(crate::node_move::DEBOUNCE));
    }
    if app.node_drag.is_none() {
        return;
    }
    // The release is handled on the next frame even when no more input comes.
    ui.ctx().request_repaint();
    if response.dragged_by(egui::PointerButton::Primary)
        && let Some(at) = response.interact_pointer_pos()
    {
        app.node_move_drag(at);
    }
    if response.drag_stopped_by(egui::PointerButton::Primary) || !response.dragged() {
        // A failure shows in the status bar.
        let _ = app.node_move_release();
    }
}

/// What a click would pick under the pointer, worked out again when the pointer, the camera, the
/// scene or the filter change.
fn hover(app: &mut ModelApp, response: &egui::Response) {
    let at = response.hover_pos().filter(|_| !response.dragged() && app.pick.dragging.is_none());
    let key = (at, app.cam, app.scene_key(), app.pick.filter, app.viewport.rect);
    if app.pick.hover_key == Some(key) {
        return;
    }
    app.pick.hover_key = Some(key);
    app.pick.hover_at = at;
    let item = at.and_then(|p| app.pick_at(p));
    app.set_hover(item);
}

/// The hovered item's name beside the pointer.
fn hover_label(app: &ModelApp, painter: &egui::Painter, t: &Tokens) {
    let (Some(item), Some(at)) = (&app.pick.hover, app.pick.hover_at) else { return };
    if !painter.clip_rect().contains(at) {
        return;
    }
    let galley = painter.layout_no_wrap(item.describe(), egui::FontId::proportional(12.0), t.overlay_text);
    let r = Rect::from_min_size(at + vec2(14.0, 14.0), galley.size()).expand2(vec2(5.0, 2.0));
    painter.rect_filled(r, 3.0, t.overlay);
    painter.galley(r.min + vec2(5.0, 2.0), galley, t.overlay_text);
}

/// The right-click menu: a set from the selection, clearing it, and the filter.
/// Most screen distance (points) from a datum's drawn lines at which a right click is on it.
const DATUM_REACH: f32 = 6.0;

/// The datum drawn under window point `at`: the nearest shown datum whose lines pass within
/// [`DATUM_REACH`].
fn datum_at(app: &mut ModelApp, at: Pos2) -> Option<String> {
    if app.plot().is_some() || !app.ui.show_geometry {
        return None;
    }
    let bounds = app.display().bounds;
    let colors = Tokens::get(app.ui.dark).glyph_colors();
    let names: Vec<String> =
        app.session.model.datums().keys().filter(|n| app.session.shows(fieldcraft_engine::ItemKind::Datum, n)).cloned().collect();
    let mut best: Option<(f32, String)> = None;
    for name in names {
        let mut drawn = fieldcraft_render::Scene::default();
        fieldcraft_render::glyphs::datums(&app.session.model, &bounds, &|n| n == name, None, &colors, &mut drawn);
        for l in &drawn.lines {
            let (Some(a), Some(b)) = (to_screen(app, l.a), to_screen(app, l.b)) else { continue };
            let xy = |p: Pos2| [f64::from(p.x), f64::from(p.y)];
            let d = fieldcraft_render::pick::segment_distance(xy(at), xy(a), xy(b)).0 as f32;
            if d <= DATUM_REACH && best.as_ref().is_none_or(|(e, _)| d < *e) {
                best = Some((d, name.clone()));
            }
        }
    }
    best.map(|(_, n)| n)
}

/// The partition whose cut face is under window point `at` (`partition:<name>:a|b`).
fn partition_at(app: &mut ModelApp, at: Pos2) -> Option<String> {
    let rect = app.viewport.rect.filter(|r| r.contains(at))?;
    let proj = app.projector()?;
    let local = [f64::from(at.x - rect.left()), f64::from(at.y - rect.top())];
    let item = app.with_picker(|p| p.pick(&proj, local, Filter::Faces))?;
    let crate::selection::Item::Geometry(g) = item else { return None };
    let name = g.entity.name()?.strip_prefix("partition:")?.split(':').next()?.to_string();
    app.session.model.partitions().iter().any(|p| p.name == name).then_some(name)
}

/// What a right click at `at` offers its item menu for: a datum drawn there, else a partition
/// whose cut face is there.
pub(crate) fn menu_target(app: &mut ModelApp, at: Pos2) -> Option<Target> {
    if let Some(name) = datum_at(app, at) {
        return Some(Target::Def(Def::Datum, name));
    }
    partition_at(app, at).map(|name| Target::Def(Def::Partition, name))
}

fn context_menu(app: &mut ModelApp, response: &egui::Response) {
    if response.secondary_clicked() {
        app.viewport.menu_target = response.interact_pointer_pos().and_then(|at| menu_target(app, at));
    }
    response.context_menu(|ui| {
        if let Some(target) = app.viewport.menu_target.clone() {
            let Target::Def(def, name) = &target else { return };
            ui.label(egui::RichText::new(format!("{} {name}", def.label())).strong());
            if let Some(action) = menu::show(app, ui, &target) {
                menu::apply(app, target, action);
                return;
            }
            ui.separator();
        }
        let has = !app.pick.selection.is_empty();
        if ui
            .add_enabled(has, egui::Button::new("Create set from selection…"))
            .on_hover_text("model.set.create with the selection as members")
            .clicked()
        {
            let _ = app.create_set_from_selection();
            ui.close();
        }
        if ui.add_enabled(has, egui::Button::new("Clear selection")).on_hover_text("Esc").clicked() {
            app.clear_selection();
            ui.close();
        }
        ui.separator();
        for f in Filter::ALL {
            if ui.radio(app.pick.filter == f, format!("Select {}", f.label().to_lowercase())).clicked() {
                app.pick.filter = f;
                ui.close();
            }
        }
    });
}

const SELECTION_BAR: &str = "fc_selection_bar";

/// Where the selection bar was drawn last: egui keeps an area's rectangle after it stops being
/// drawn, so a hidden bar (the Start screen) still answers with where it was.
pub(crate) fn selection_bar_rect(ctx: &egui::Context) -> Option<Rect> {
    ctx.memory(|m| m.area_rect(egui::Id::new(SELECTION_BAR)))
}

/// The selection bar, top centre: what a click picks, box selection, and what is selected.
fn selection_bar(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect, t: &Tokens) {
    let area = egui::Area::new(egui::Id::new(SELECTION_BAR))
        .fixed_pos(pos2(rect.center().x, rect.top() + 10.0))
        .pivot(egui::Align2::CENTER_TOP)
        .constrain_to(rect)
        .order(egui::Order::Middle);
    let shown = area.show(ui.ctx(), |ui| {
        egui::Frame::NONE.fill(t.overlay).stroke(Stroke::new(1.0, t.border)).corner_radius(6.0).inner_margin(egui::Margin::symmetric(8, 4)).show(
            ui,
            |ui| {
                ui.horizontal(|ui| {
                    ui.label(egui::RichText::new("Select").color(t.text_dim));
                    egui::ComboBox::from_id_salt("fc_selection_filter").selected_text(app.pick.filter.label()).width(110.0).show_ui(ui, |ui| {
                        for (i, f) in Filter::ALL.into_iter().enumerate() {
                            // Geometry, then the mesh.
                            if i > 0 && Filter::ALL.get(i - 1).is_some_and(|p| p.geometry() != f.geometry()) {
                                ui.separator();
                            }
                            if ui.selectable_label(app.pick.filter == f, f.label()).clicked() {
                                app.pick.filter = f;
                            }
                        }
                    });
                    let boxes = app.pick.filter.boxes();
                    let tip = if boxes {
                        "Box selection: a left drag selects the visible items inside the box (right or middle drag pans)"
                    } else {
                        "Box selection picks nodes, element faces and elements"
                    };
                    if ui.add_enabled(boxes, egui::Button::selectable(app.pick.box_mode && boxes, "Box")).on_hover_text(tip).clicked() {
                        app.pick.box_mode = !app.pick.box_mode;
                    }
                    if !app.pick.selection.is_empty() {
                        ui.separator();
                        let mut text = app.pick.selection.summary();
                        if app.pick.selection.capped() {
                            text.push_str(" (capped)");
                        }
                        ui.label(egui::RichText::new(text).color(t.pick_selected));
                        if ui.small_button("×").on_hover_text("Clear the selection (Esc)").clicked() {
                            app.clear_selection();
                        }
                    }
                });
            },
        );
    });
    app.place("selection_bar", shown.response.rect);
}

/// Camera navigation. A drag's effect is chosen when it starts ([`navigation::drag_mode`]: the
/// preset's bindings, the bar's mode for the left button); an orbit turns about the model point
/// under the press, or the centre of what is drawn when it starts on empty space, and keeps that
/// point still on screen. The wheel (raw, so a notch zooms in full on the frame it arrives) and a
/// pinch zoom toward the model point under the cursor, or the point at the target's depth.
fn navigate(app: &mut ModelApp, ui: &egui::Ui, response: &egui::Response, rect: Rect, boxing: bool) {
    let (w, h) = (f64::from(rect.width()), f64::from(rect.height()));
    let (mods, press) = ui.input(|i| (i.modifiers, i.pointer.press_origin()));
    for (pressed, button) in [
        (egui::PointerButton::Primary, NavButton::Left),
        (egui::PointerButton::Middle, NavButton::Middle),
        (egui::PointerButton::Secondary, NavButton::Right),
    ] {
        // A Move Nodes drag is no navigation.
        if response.drag_started_by(pressed) && app.node_drag.is_none() {
            let at = press.or(response.interact_pointer_pos());
            let preset = navigation::preset(&app.ui.nav_preset);
            // The pick is needed only for an orbit's pivot and an unbound Ctrl + left drag.
            let picked = preset.mode(button, mods).is_none() && button == NavButton::Left && (mods.ctrl || mods.command);
            let hit = if picked { at.and_then(|a| app.surface_point(a)) } else { None };
            app.viewport.drag = navigation::drag_mode(preset, button, mods, app.nav, boxing, hit.is_some());
            app.viewport.orbit_pivot = None;
            if app.viewport.drag == Some(NavMode::Orbit) {
                let hit = if picked { hit } else { at.and_then(|a| app.surface_point(a)) };
                app.viewport.orbit_pivot = Some(hit.unwrap_or_else(|| app.drawn_centre()));
            }
        }
    }
    let d = response.drag_delta();
    let (dx, dy) = (f64::from(d.x), f64::from(d.y));
    if response.dragged() && (dx != 0.0 || dy != 0.0) {
        match app.viewport.drag {
            Some(NavMode::Orbit) => match app.viewport.orbit_pivot {
                Some(p) => app.cam.orbit_about(p, dx, dy),
                None => app.cam.orbit(dx, dy),
            },
            Some(NavMode::Pan) => app.cam.pan(dx, dy, h),
            Some(NavMode::Zoom) => app.cam.zoom_at((dy * ZOOM_PER_POINT * 4.0 * zoom_sign(app)).exp(), None),
            None => {}
        }
    }
    if !response.dragged() {
        app.viewport.drag = None;
        app.viewport.orbit_pivot = None;
    }
    if response.double_clicked() || (response.hovered() && ui.input(|i| i.key_pressed(egui::Key::F) && i.modifiers.is_none())) {
        app.fit_view();
    }
    if response.hovered() {
        let options = ui.ctx().options(|o| o.input_options);
        let (scroll, pinch, at) = ui.input(|i| (wheel_points(i, &options), i.zoom_delta(), i.pointer.hover_pos()));
        let factor = (-f64::from(scroll) * ZOOM_PER_POINT * zoom_sign(app)).exp() / f64::from(pinch).max(MIN_PINCH);
        // Alt + scroll sizes the Relax Nodes sphere, or the Move Nodes radius, instead of zooming.
        let alt = ui.input(|i| i.modifiers.alt) && scroll != 0.0;
        let step = (f64::from(scroll) * ZOOM_PER_POINT).exp();
        let sized = alt && ((app.relaxing() && app.relax_scale_radius(step)) || (app.moving_nodes() && app.node_move_scale_radius(step)));
        if !sized && (factor - 1.0).abs() > NO_ZOOM {
            let anchor = at.and_then(|p| {
                app.surface_point(p).or_else(|| app.cam.target_plane_point(f64::from(p.x - rect.left()), f64::from(p.y - rect.top()), w, h))
            });
            app.cam.zoom_at(factor, anchor);
        }
    }
}

/// This frame's wheel and trackpad scrolling, vertical, in points, as it arrived: egui's smoothed
/// delta spreads a notch over the next frames and loses its tail when no more frames come (a
/// hidden window). Lines count as egui counts them; scrolling with the zoom modifier is egui's
/// zoom ([`egui::InputState::zoom_delta`]) and left out.
fn wheel_points(i: &egui::InputState, options: &egui::InputOptions) -> f32 {
    i.events
        .iter()
        .filter_map(|e| match e {
            egui::Event::MouseWheel { unit, delta, modifiers, .. } if !modifiers.matches_any(options.zoom_modifier) => Some(match unit {
                egui::MouseWheelUnit::Point => delta.y,
                egui::MouseWheelUnit::Line => delta.y * options.line_scroll_speed,
                egui::MouseWheelUnit::Page => delta.y * i.viewport_rect().height(),
            }),
            _ => None,
        })
        .filter(|d| d.is_finite())
        .sum()
}

/// −1 when the wheel is reversed (the preference), else 1.
fn zoom_sign(app: &ModelApp) -> f64 {
    if app.ui.reverse_zoom { -1.0 } else { 1.0 }
}

/// The orbit's pivot, while orbiting: a small ring with a dot.
fn pivot_marker(app: &ModelApp, painter: &egui::Painter, t: &Tokens) {
    if let Some(at) = app.viewport.orbit_pivot.and_then(|p| to_screen(app, p)) {
        painter.circle_filled(at, 2.5, t.accent);
        painter.circle_stroke(at, 6.0, Stroke::new(1.5, t.accent));
    }
}

/// A vertical gradient behind the model.
fn gradient(painter: &egui::Painter, rect: Rect, top: Color32, bottom: Color32) {
    let mut mesh = egui::Mesh::default();
    mesh.colored_vertex(rect.left_top(), top);
    mesh.colored_vertex(rect.right_top(), top);
    mesh.colored_vertex(rect.right_bottom(), bottom);
    mesh.colored_vertex(rect.left_bottom(), bottom);
    mesh.add_triangle(0, 1, 2);
    mesh.add_triangle(0, 2, 3);
    painter.add(egui::Shape::mesh(mesh));
}

fn paint_gpu(app: &mut ModelApp, ui: &egui::Ui, painter: &egui::Painter, rect: Rect) {
    let ppp = ui.ctx().pixels_per_point();
    // The symbols first: the scene's clipping radius holds them.
    let symbols_key = app.symbols().0;
    if app.viewport.symbols_uploaded != Some(symbols_key)
        && let Ok(mut slot) = app.viewport.symbols_slot.lock()
    {
        *slot = Some(GpuScene::from_scene(&app.symbols.1, ppp));
        app.viewport.symbols_uploaded = Some(symbols_key);
    }
    let key = app.scene().0;
    let radius = app.scene.1.radius;
    if app.viewport.uploaded != Some(key) {
        let g = GpuScene::from_scene(&app.scene.1, ppp);
        if let Ok(mut slot) = app.viewport.slot.lock() {
            *slot = Some(g);
            app.viewport.uploaded = Some(key);
        }
    }
    let cam = app.cam;
    let (overlay_key, overlay) = app.overlay();
    if app.viewport.overlay_uploaded != Some(overlay_key)
        && let Ok(mut slot) = app.viewport.overlay_slot.lock()
    {
        *slot = Some(GpuScene::from_scene(&overlay, ppp));
        app.viewport.overlay_uploaded = Some(overlay_key);
    }
    let cb = ViewportCallback {
        key,
        slot: app.viewport.slot.clone(),
        symbols_key,
        symbols_slot: app.viewport.symbols_slot.clone(),
        overlay_key,
        overlay_slot: app.viewport.overlay_slot.clone(),
        frame: Frame::of(&cam, rect.width() * ppp, rect.height() * ppp, radius),
    };
    painter.add(egui_wgpu::Callback::new_paint_callback(rect, cb));
}

/// Without a GPU: the CPU rasterizer into a texture, redrawn when the scene, camera or size
/// changes.
fn paint_cpu(app: &mut ModelApp, ui: &egui::Ui, painter: &egui::Painter, rect: Rect) {
    let ppp = ui.ctx().pixels_per_point();
    // A window larger than the renderer's limit, or than the largest texture the context takes,
    // is drawn smaller and scaled up.
    let (w, h) = cpu_size((rect.width() * ppp) as usize, (rect.height() * ppp) as usize, ui.ctx().input(|i| i.max_texture_side));
    let size = [w, h];
    let cam = app.cam;
    let (overlay_key, overlay) = app.overlay();
    let symbols_key = app.symbols().0;
    let key = (app.scene().0, symbols_key, overlay_key);
    let fresh = app.viewport.cpu.as_ref().is_some_and(|(_, k, c, s)| *k == key && *c == cam && *s == size);
    // A scene the renderer refused is reported once and not retried until it changes.
    if !fresh && app.viewport.cpu_failed != Some(key.0) {
        let c = match fieldcraft_render::render_layers(&app.scene.1, &[&app.symbols.1, &overlay], &cam, w, h) {
            Ok(c) => c,
            Err(e) => {
                app.viewport.cpu = None;
                app.viewport.cpu_failed = Some(key.0);
                app.set_status(format!("Viewport: {e}"), true);
                return;
            }
        };
        app.viewport.cpu_failed = None;
        let image = egui::ColorImage::from_rgba_unmultiplied([c.w, c.h], &c.rgba);
        let tex = match app.viewport.cpu.take() {
            Some((mut tex, ..)) => {
                tex.set(image, egui::TextureOptions::LINEAR);
                tex
            }
            None => ui.ctx().load_texture("fc_viewport_cpu", image, egui::TextureOptions::LINEAR),
        };
        app.viewport.cpu = Some((tex, key, cam, size));
    }
    if let Some((tex, ..)) = &app.viewport.cpu {
        painter.image(tex.id(), rect, Rect::from_min_max(pos2(0.0, 0.0), pos2(1.0, 1.0)), Color32::WHITE);
    }
}

/// The CPU viewport image's size for a viewport of `w` × `h` pixels: within the renderer's limits
/// and no side above `max_side`, the aspect kept.
fn cpu_size(w: usize, h: usize, max_side: usize) -> (usize, usize) {
    let (w, h) = fieldcraft_render::fit_size(w, h);
    let longest = w.max(h);
    if max_side == 0 || longest <= max_side {
        return (w, h);
    }
    let k = max_side as f64 / longest as f64;
    let fit = |v: usize| ((v as f64 * k).floor() as usize).clamp(1, max_side);
    (fit(w), fit(h))
}

/// The scene's labels (reference point names) beside their points, with the interface's font.
/// The CPU renderer draws them into its image.
fn labels(app: &mut ModelApp, painter: &egui::Painter, t: &Tokens) {
    app.symbols();
    app.scene();
    let labels: Vec<(Vec3, String)> = app.scene.1.labels.iter().chain(&app.symbols.1.labels).map(|l| (l.p, l.text.clone())).collect();
    let offset = fieldcraft_render::LABEL_OFFSET as f32;
    for (p, text) in labels {
        if let Some(at) = to_screen(app, p) {
            let galley = painter.layout_no_wrap(text, egui::FontId::proportional(12.0), t.overlay_text);
            let r = Rect::from_min_size(at + vec2(offset, -offset - galley.size().y), galley.size()).expand2(vec2(4.0, 1.0));
            painter.rect_filled(r, 3.0, t.overlay);
            painter.galley(r.min + vec2(4.0, 1.0), galley, t.overlay_text);
        }
    }
}

/// Bottom left: the definitions whose symbols are a sample of their nodes or faces, wrapped to
/// the room left of the navigation bar (lines are added upwards).
fn sampled_note(app: &ModelApp, painter: &egui::Painter, rect: Rect, t: &Tokens) {
    let show = app.ui.show();
    let notes = &app.glyph_notes();
    if notes.is_empty() || !(show.loads || show.interactions) {
        return;
    }
    let text = format!("Symbols sampled: {}", notes.join("; "));
    let room = sampled_note_width(rect);
    let galley = painter.layout(text, egui::FontId::proportional(11.0), t.text_dim, room);
    painter.galley(rect.left_bottom() + vec2(10.0, -10.0 - galley.size().y), galley, t.text_dim);
}

/// How wide the sampled-symbols note may be: up to the navigation bar, with a gap.
fn sampled_note_width(rect: Rect) -> f32 {
    let bar = nav_bar_rect(rect);
    let right = if rect.contains_rect(bar) { bar.left() - 12.0 } else { rect.right() - 10.0 };
    (right - rect.left() - 10.0).max(80.0)
}

/// The orientation gizmo's arm length and tip radius (points).
const AXIS_LENGTH: f32 = 38.0;
const AXIS_RADIUS: f32 = 11.0;

/// The standard view looking along `d` from outside (the one whose eye side is nearest `d`).
fn view_along(d: Vec3) -> StandardView {
    let mut best = (StandardView::Front, f64::NEG_INFINITY);
    for v in StandardView::ALL.into_iter().filter(|v| *v != StandardView::Iso) {
        let mut c = fieldcraft_render::Camera::default();
        c.set_view(v);
        let score = c.back().dot(d);
        if score > best.1 {
            best = (v, score);
        }
    }
    best.0
}

fn axis_view_rect(rect: Rect) -> Rect {
    Rect::from_center_size(pos2(rect.right() - 68.0, rect.top() + 68.0), Vec2::splat(116.0))
}

/// The orientation gizmo (top right): the six signed axes in the current view, back to front;
/// clicking a tip looks along it, the house returns to the home view.
fn axis_view(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect, t: &Tokens) {
    let area = axis_view_rect(rect);
    if !rect.contains_rect(area) {
        return;
    }
    app.place("axis_view", area);
    let c = area.center();
    let (r, u, b) = app.cam.basis();
    let mut tips: Vec<(&str, Vec3, Pos2, f64, Color32)> = [
        ("X", Vec3::X, t.axis_x),
        ("-X", Vec3::new(-1.0, 0.0, 0.0), t.axis_x),
        ("Y", Vec3::Y, t.axis_y),
        ("-Y", Vec3::new(0.0, -1.0, 0.0), t.axis_y),
        ("Z", Vec3::Z, t.axis_z),
        ("-Z", Vec3::new(0.0, 0.0, -1.0), t.axis_z),
    ]
    .into_iter()
    .map(|(name, d, col)| (name, d, c + vec2(d.dot(r) as f32 * AXIS_LENGTH, -(d.dot(u) as f32) * AXIS_LENGTH), d.dot(b), col))
    .collect();
    tips.sort_by(|a, b2| a.3.total_cmp(&b2.3));
    let hover = ui.input(|i| i.pointer.hover_pos()).filter(|p| area.contains(*p));
    let hovered = hover.and_then(|p| tips.iter().rposition(|(_, _, at, _, _)| at.distance(p) <= AXIS_RADIUS));
    let painter = ui.painter_at(rect);
    if hover.is_some() {
        painter.circle_filled(c, 55.0, t.overlay);
    }
    for (i, (name, _, at, depth, col)) in tips.iter().enumerate() {
        let positive = !name.starts_with('-');
        let ink = if *depth < -0.01 { col.gamma_multiply(0.55) } else { *col };
        painter.line_segment([c, *at], Stroke::new(if positive { 2.0 } else { 1.0 }, ink));
        let hot = hovered == Some(i);
        painter.circle_filled(*at, AXIS_RADIUS, if hot { t.accent_soft } else { t.overlay });
        painter.circle_stroke(*at, AXIS_RADIUS, Stroke::new(if hot { 2.0 } else { 1.0 }, if hot { t.accent } else { ink }));
        painter.text(*at, egui::Align2::CENTER_CENTER, *name, egui::FontId::proportional(11.0), ink);
    }
    let resp = ui.interact(area, ui.id().with("fc_axis_view"), Sense::click());
    if let Some(i) = hovered
        && let Some((name, d, ..)) = tips.get(i)
    {
        let view = view_along(*d);
        let resp = resp.clone().on_hover_text(format!("Look from {name} ({})", view.label()));
        if resp.clicked() {
            app.set_view(view);
        }
    }
    let home = Rect::from_center_size(pos2(c.x - 44.0, c.y - 44.0), Vec2::splat(18.0));
    let hr = ui.interact(home, ui.id().with("fc_axis_home"), Sense::click());
    if hr.hovered() {
        painter.rect_filled(home.expand(2.0), 3.0, t.hover);
    }
    icons::paint(app, ui, "ui.home", home, Look::Normal);
    if hr.on_hover_text("Home view").clicked() {
        app.home_view();
    }
}

/// The navigation bar's buttons: (key, icon, tooltip).
const NAV: [(&str, &str, &str); 7] = [
    ("orbit", "ui.orbit", "Orbit: left drag turns the view"),
    ("pan", "ui.pan", "Pan: left drag moves the view (or right or middle drag)"),
    ("zoom", "ui.zoom", "Zoom: left drag up and down (or the wheel)"),
    ("fit", "ui.fit", "Fit the model (F, or a double click)"),
    ("home", "ui.home", "Home view"),
    ("perspective", "ui.perspective", "Perspective or orthographic"),
    ("style", "ui.style", "Visual style and what is shown"),
];

fn nav_bar_rect(rect: Rect) -> Rect {
    Rect::from_center_size(pos2(rect.center().x, rect.bottom() - 22.0), vec2(NAV.len() as f32 * 30.0 + 10.0, 30.0))
}

/// The navigation bar, bottom centre: navigation modes, fit, home, projection, and the display
/// menu (style, geometry, mesh, nodes).
fn nav_bar(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect, t: &Tokens) {
    let bar = nav_bar_rect(rect);
    if !rect.contains_rect(bar) {
        return;
    }
    app.place("nav_bar", bar);
    ui.painter().rect(bar, 6.0, t.overlay, Stroke::new(1.0, t.border), egui::StrokeKind::Inside);
    let mut style_button = None;
    for (i, (key, icon, tip)) in NAV.iter().enumerate() {
        let br = Rect::from_min_size(pos2(bar.left() + 5.0 + i as f32 * 30.0, bar.top() + 2.0), Vec2::splat(26.0));
        let resp = ui.interact(br, ui.id().with(("fc_nav", *key)), Sense::click());
        let active = match *key {
            "orbit" => app.nav == NavMode::Orbit,
            "pan" => app.nav == NavMode::Pan,
            "zoom" => app.nav == NavMode::Zoom,
            "perspective" => app.ui.perspective,
            "style" => app.ui.menu.as_deref() == Some("view"),
            _ => false,
        };
        if active || resp.hovered() {
            ui.painter().rect_filled(br, 4.0, if active { t.accent_soft } else { t.hover });
        }
        icons::paint(app, ui, icon, br.shrink(4.0), Look::Normal);
        resp.widget_info(|| egui::WidgetInfo::selected(egui::WidgetType::Button, true, active, *tip));
        if resp.on_hover_text(*tip).clicked() {
            match *key {
                "orbit" => app.nav = NavMode::Orbit,
                "pan" => app.nav = NavMode::Pan,
                "zoom" => app.nav = NavMode::Zoom,
                "fit" => app.fit_view(),
                "home" => {
                    app.home_view();
                }
                "perspective" => app.ui.perspective = !app.ui.perspective,
                _ => app.toggle_menu("view"),
            }
        }
        if *key == "style" {
            style_button = Some(br);
        }
    }
    if let Some(br) = style_button {
        crate::frame::settings_menu(app, ui, "view", br, |app, ui| {
            let t = Tokens::get(app.ui.dark);
            crate::frame::menu_heading(ui, &t, "STYLE");
            for s in Style::ALL {
                if ui.radio(app.ui.style == s, s.label()).clicked() {
                    app.ui.style = s;
                }
            }
            ui.separator();
            crate::frame::menu_heading(ui, &t, "SHOW");
            ui.checkbox(&mut app.ui.show_geometry, "Geometry");
            ui.checkbox(&mut app.ui.show_mesh, "Mesh");
            ui.checkbox(&mut app.ui.show_nodes, "Nodes");
            ui.checkbox(&mut app.ui.show_loads, "Loads & BCs");
            ui.checkbox(&mut app.ui.show_interactions, "Interactions");
            let (lo, hi) = fieldcraft_render::glyphs::SIZE_RANGE;
            ui.add(egui::Slider::new(&mut app.ui.symbol_size, lo..=hi).step_by(1.0).suffix(" pt").text("Symbol size"))
                .on_hover_text("The size of boundary condition, load and interaction symbols on screen, at any zoom");
            ui.separator();
            ui.checkbox(&mut app.ui.show_tree, "Browser");
            ui.checkbox(&mut app.ui.show_timeline, "Timeline");
            ui.separator();
            crate::frame::menu_heading(ui, &t, "NAVIGATION");
            let current = navigation::preset(&app.ui.nav_preset);
            egui::ComboBox::from_id_salt("fc_nav_preset")
                .selected_text(current.label)
                .show_ui(ui, |ui| {
                    for p in &navigation::NAV_PRESETS {
                        if ui.selectable_label(p.id == current.id, p.label).on_hover_text(p.help).clicked() {
                            app.ui.nav_preset = p.id.into();
                            app.ui.reverse_zoom = p.reverse_zoom;
                        }
                    }
                })
                .response
                .on_hover_text(current.help);
            ui.checkbox(&mut app.ui.reverse_zoom, "Reverse zoom").on_hover_text("The wheel zooms the other way");
        });
    }
}

/// The screen point (points) of a world point in the viewport drawn last frame.
pub fn to_screen(app: &ModelApp, p: Vec3) -> Option<Pos2> {
    let rect = app.viewport.rect?;
    // The projection picking uses ([`ModelApp::projector`]).
    let (x, y, _) = app.cam.to_screen(p, f64::from(rect.width()), f64::from(rect.height()), app.scene.1.radius)?;
    Some(pos2(rect.left() + x as f32, rect.top() + y as f32))
}

#[cfg(test)]
mod tests {
    use super::{cpu_size, nav_bar_rect, sampled_note_width};

    /// The CPU viewport never asks for a texture larger than the context takes (#237: a Retina
    /// window 2300 pixels wide on a 2048 limit).
    #[test]
    fn the_cpu_image_fits_the_texture_limit() {
        assert_eq!(cpu_size(2300, 1412, 2048), (2048, 1257));
        assert_eq!(cpu_size(1412, 2300, 2048), (1257, 2048));
        assert_eq!(cpu_size(800, 600, 2048), (800, 600));
        assert_eq!(cpu_size(800, 600, 0), (800, 600));
        let (w, h) = cpu_size(100_000, 3, 2048);
        assert!(w <= 2048 && h >= 1);
    }

    /// The sampled-symbols note stops short of the navigation bar (#206).
    #[test]
    fn the_sampled_note_stays_left_of_the_navigation_bar() {
        let view = egui::Rect::from_min_size(egui::pos2(300.0, 80.0), egui::vec2(900.0, 700.0));
        let w = sampled_note_width(view);
        assert!(view.left() + 10.0 + w < nav_bar_rect(view).left());
    }
}
