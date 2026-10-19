//! The application bar: the window's title bar. The app's mark and name, the quick actions (File
//! menu, save, undo, redo, Solve, command search, theme, help), the model's document tab and the
//! version; on Windows and Linux the caption buttons ([`crate::titlebar`]).
//!
//! Layout follows SolidCraft's `toolbar::app_bar` (MIT OR Apache-2.0), see ATTRIBUTION.md; the
//! quick-action buttons are craft-ui `IconButton`s and the document tab a craft-ui `Tab`.

use craft_ui::buttons::{ButtonVisuals, IconButton, IconButtonStyle};
use egui::{Align2, Color32, CornerRadius, FontId, Rect, Sense, Stroke, Vec2, pos2, vec2};
use serde_json::json;

use crate::ModelApp;
use crate::icons::{self, Look};
use crate::theme::Tokens;

/// The bar's height.
pub const HEIGHT: f32 = 34.0;
/// A quick-action button's side, and the room each takes.
const BUTTON: f32 = 26.0;
const STRIDE: f32 = 30.0;

/// The title bar's quick actions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Quick {
    File,
    Save,
    Undo,
    Redo,
    Solve,
    Search,
    Theme,
    Help,
}

impl Quick {
    /// Its region in `ui.inspect`'s layout.
    fn region(self) -> &'static str {
        match self {
            Quick::File => "app_bar:file",
            Quick::Save => "app_bar:save",
            Quick::Undo => "app_bar:undo",
            Quick::Redo => "app_bar:redo",
            Quick::Solve => "app_bar:solve",
            Quick::Search => "app_bar:search",
            Quick::Theme => "app_bar:theme",
            Quick::Help => "app_bar:help",
        }
    }
}

/// The quick actions' look: flat on the dark bar, lighter under the pointer.
fn button_style(t: &Tokens) -> IconButtonStyle {
    let normal = ButtonVisuals { fill: Color32::TRANSPARENT, stroke: Stroke::NONE, icon: t.app_bar_text };
    IconButtonStyle {
        normal,
        hovered: ButtonVisuals { fill: Color32::from_white_alpha(30), ..normal },
        pressed: Some(ButtonVisuals { fill: Color32::from_white_alpha(18), ..normal }),
        selected: ButtonVisuals { fill: Color32::from_white_alpha(40), ..normal },
        corner_radius: CornerRadius::same(4),
        focus_stroke: Stroke::new(1.0, t.accent),
    }
}

/// FieldCraft's mark: the application icon's fringe cube (`assets/app-icon/icon-small.svg`, after
/// `xtask/src/app_icon.rs`) drawn in code without its tile or streamline: each visible face split
/// into small cells coloured by the fringe band at the cell's centre, shaded per face.
fn mark(p: &egui::Painter, c: egui::Pos2) {
    // The cube is 2 units tall on screen (half side 0.5, isometric); draw it 22 points tall.
    const S: f32 = 11.0;
    const N: usize = 8;
    const FRINGE: [(u8, u8, u8); 7] =
        [(0x2c, 0x4f, 0xbf), (0x2f, 0x86, 0xe0), (0x22, 0xbc, 0xd8), (0x3f, 0xc5, 0x84), (0xf0, 0xd1, 0x3a), (0xf2, 0x8a, 0x2b), (0xde, 0x3b, 0x2f)];
    const RADII: [f32; 6] = [1.28, 1.07, 0.87, 0.67, 0.47, 0.27];
    let at = |q: [f32; 3]| c + vec2((q[0] - q[1]) * 0.866 * S, ((q[0] + q[1]) * 0.5 - q[2]) * S);
    // Each face: origin corner, two edge vectors, and how much it is darkened (0 to 1).
    let faces: [([f32; 3], [f32; 3], [f32; 3], f32); 3] = [
        ([-0.5, -0.5, 0.5], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
        ([-0.5, 0.5, -0.5], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.12),
        ([0.5, -0.5, -0.5], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 0.3),
    ];
    let point = |o: [f32; 3], u: [f32; 3], v: [f32; 3], a: f32, b: f32| [0, 1, 2].map(|k| o[k] + u[k] * a + v[k] * b);
    for (o, u, v, shade) in faces {
        for i in 0..N {
            for j in 0..N {
                let (a0, a1, b0, b1) = (i as f32 / N as f32, (i + 1) as f32 / N as f32, j as f32 / N as f32, (j + 1) as f32 / N as f32);
                let m = point(o, u, v, (a0 + a1) / 2.0, (b0 + b1) / 2.0);
                let r = m.iter().map(|x| (0.5 - x) * (0.5 - x)).sum::<f32>().sqrt();
                let band = RADII.iter().filter(|&&rad| r < rad).count();
                let (cr, cg, cb) = FRINGE.get(band).copied().unwrap_or(FRINGE[0]);
                let k = 1.0 - shade;
                let col = Color32::from_rgb((cr as f32 * k) as u8, (cg as f32 * k) as u8, (cb as f32 * k) as u8);
                let quad = vec![at(point(o, u, v, a0, b0)), at(point(o, u, v, a1, b0)), at(point(o, u, v, a1, b1)), at(point(o, u, v, a0, b1))];
                // A hairline in the cell's own colour closes the seams between neighbouring cells.
                p.add(egui::Shape::convex_polygon(quad, col, Stroke::new(0.5, col)));
            }
        }
    }
}

/// Where the title bar's parts go.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BarLayout {
    /// The mark's centre, and where the app's name starts (none when compact).
    pub mark: egui::Pos2,
    pub name: Option<egui::Pos2>,
    /// File, Save, Undo, Redo: after the mark, at the left.
    pub quick: [Rect; 4],
    /// Solve, then Search, Theme and Help: at the right end, before the version.
    pub solve: Rect,
    pub tools: [Rect; 3],
    /// The document tab's room, between the two groups.
    pub tab_room: (f32, f32),
    /// The version's right edge.
    pub version_right: f32,
    /// Too narrow for labels: no name, and Solve without its word.
    pub compact: bool,
}

/// The title bar's layout in `bar`. `traffic_lights`: the macOS window buttons sit over its left
/// end (an integrated title bar out of full screen, where there are none). `captions`: the width
/// of the custom caption buttons at its right end; `version_w`: the version text's width.
pub fn bar_layout(bar: Rect, traffic_lights: bool, captions: f32, version_w: f32) -> BarLayout {
    let left = bar.left() + if traffic_lights { crate::titlebar::TRAFFIC_LIGHTS } else { 0.0 };
    let compact = bar.right() - left < 760.0;
    let cy = bar.center().y;
    let button = |x: f32| Rect::from_center_size(pos2(x + BUTTON / 2.0, cy), Vec2::splat(BUTTON));
    let x0 = left + if compact { 44.0 } else { 132.0 };
    let quick = [0.0, 1.0, 2.0, 3.0].map(|i| button(x0 + i * STRIDE));
    let version_right = bar.right() - 12.0 - captions;
    let tools_end = version_right - version_w - 14.0;
    let tools = [3.0, 2.0, 1.0].map(|i| button(tools_end - i * STRIDE));
    let solve_w = if compact { 30.0 } else { 74.0 };
    let solve = Rect::from_min_size(pos2(tools[0].left() - 10.0 - solve_w, cy - 12.0), vec2(solve_w, 24.0));
    BarLayout {
        mark: pos2(left + 20.0, cy),
        name: (!compact).then(|| pos2(left + 36.0, cy)),
        quick,
        solve,
        tools,
        tab_room: (quick[3].right() + 16.0, solve.left() - 16.0),
        version_right,
        compact,
    }
}

/// Whether the window is in full screen (no traffic lights to make room for).
fn fullscreen(ctx: &egui::Context) -> bool {
    ctx.input(|i| i.viewport().fullscreen.unwrap_or(false))
}

pub fn app_bar(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let style = button_style(&t);
    let mut anchors: Vec<(Quick, Rect)> = Vec::new();
    egui::Panel::top("fc_appbar").exact_size(HEIGHT).show_separator_line(false).frame(egui::Frame::NONE.fill(t.app_bar)).show(ui, |ui| {
        let r = ui.max_rect();
        app.place("app_bar", r);
        if app.custom_titlebar || app.integrated_titlebar {
            crate::titlebar::drag_area(ui, r);
        }
        let captions = if app.custom_titlebar { crate::titlebar::WIDTH } else { 0.0 };
        let version = format!("v{}", crate::VERSION);
        let version_w = ui.painter().layout_no_wrap(version.clone(), FontId::proportional(11.0), t.text_dim).size().x;
        let lay = bar_layout(r, app.integrated_titlebar && !fullscreen(ui.ctx()), captions, version_w);
        mark(ui.painter(), lay.mark);
        if let Some(at) = lay.name {
            ui.painter().text(at, Align2::LEFT_CENTER, "FieldCraft", FontId::proportional(14.0), t.app_bar_text);
        }
        let (undo, redo) = app.session.undo_depth();
        let busy = app.busy().is_some();
        let file_open = app.ui.menu.as_deref() == Some("file");
        let actions: [(Quick, &'static str, String, bool, bool); 4] = [
            (Quick::File, "file.open", "File: open, save, decks, samples".to_string(), true, file_open),
            (Quick::Save, "file.save", format!("Save ({})", shortcut("S")), !busy, false),
            (Quick::Undo, "edit.undo", format!("Undo ({})", shortcut("Z")), undo > 0 && !busy, false),
            (Quick::Redo, "edit.redo", format!("Redo ({})", shortcut("Shift+Z")), redo > 0 && !busy, false),
        ];
        let tools: [(Quick, &'static str, String, bool, bool); 3] = [
            (Quick::Search, "ui.search", format!("Search commands (S or {})", shortcut("K")), true, app.ui.palette),
            (
                Quick::Theme,
                if app.ui.dark { "ui.sun" } else { "ui.moon" },
                if app.ui.dark { "Light theme".into() } else { "Dark theme".into() },
                true,
                false,
            ),
            (Quick::Help, "ui.help", "Help: samples, navigation, about".into(), true, app.ui.menu.as_deref() == Some("help")),
        ];
        let mut clicked = None;
        for ((q, icon, tip, enabled, selected), rect) in actions.into_iter().zip(lay.quick).chain(tools.into_iter().zip(lay.tools)) {
            if quick_button(app, ui, &style, q.region(), icon, rect, &tip, enabled, selected) {
                clicked = Some(q);
            }
            anchors.push((q, rect));
        }
        // Solve: the one prominent action, a filled button with its label.
        let solve = lay.solve;
        // Solve runs `solve.run`: greyed while the model misses what it needs, every missing
        // piece listed in the tooltip.
        let missing = app.unavailable("solve.run");
        let can = missing.is_none() && !busy;
        let resp = ui.interact(solve, ui.id().with("fc_solve"), if can { Sense::click() } else { Sense::hover() });
        let fill = if resp.hovered() && can {
            t.accent.gamma_multiply(1.15)
        } else if can {
            t.accent
        } else {
            Color32::from_white_alpha(24)
        };
        ui.painter().rect_filled(solve, 5.0, fill);
        // A white play mark: the accent fill is the button's colour.
        let c = pos2(solve.left() + 15.0, solve.center().y);
        ui.painter().add(egui::Shape::convex_polygon(
            vec![c + vec2(-3.5, -5.5), c + vec2(5.5, 0.0), c + vec2(-3.5, 5.5)],
            Color32::WHITE,
            Stroke::NONE,
        ));
        if !lay.compact {
            ui.painter().text(
                pos2(solve.left() + 28.0, solve.center().y),
                Align2::LEFT_CENTER,
                "Solve",
                crate::theme::label_font(ui.ctx(), 13.0),
                Color32::WHITE,
            );
        }
        let head = if missing.is_some() { "Solve: not ready yet (or try Help › Samples)" } else { "Solve every step and show the results" };
        let tip = fieldcraft_engine::requirement::Unavailable::tip(head, missing.as_ref());
        resp.widget_info(|| egui::WidgetInfo::labeled(egui::WidgetType::Button, can, "Solve"));
        if resp.on_hover_text(tip).clicked() && can {
            clicked = Some(Quick::Solve);
        }
        anchors.push((Quick::Solve, solve));
        // The model's document tab, joined to the ribbon below it, centred in its room.
        let migrated = app.session.migrated_from.map(|n| format!(" (from format {n})")).unwrap_or_default();
        let label = format!("{}{migrated}{}", app.model_name(), if app.session.is_dirty() { " •" } else { "" });
        let text_w = ui.painter().layout_no_wrap(label.clone(), FontId::proportional(13.0), t.text).size().x;
        let (room_l, room_r) = lay.tab_room;
        let tab_w = (text_w + 40.0).clamp(120.0, 260.0).min((room_r - room_l).max(0.0));
        let tab_x = (r.center().x - tab_w / 2.0).clamp(room_l, (room_r - tab_w).max(room_l));
        if tab_w > 40.0 {
            let tab = Rect::from_min_size(pos2(tab_x, r.top() + 5.0), vec2(tab_w, r.height() - 5.0));
            let mut tip = app.session.path.clone().unwrap_or_else(|| "Not saved yet".into());
            if let Some(n) = app.session.migrated_from {
                tip.push_str(&format!("\nOpened from a format {n} file: the first save keeps that file as a copy beside it"));
                for note in &app.session.migration_notes {
                    tip.push_str(&format!("\n{note}"));
                }
            }
            craft_ui::tabs::Tab::new(egui::Id::new("fc_document_tab"), &label, true).show_at(ui, tab, |ui, resp| {
                ui.painter().rect_filled(tab, CornerRadius { nw: 5, ne: 5, sw: 0, se: 0 }, t.toolbar);
                let clip = tab.shrink2(vec2(10.0, 0.0));
                ui.painter().with_clip_rect(clip).text(
                    pos2(tab.left() + 14.0, tab.center().y),
                    Align2::LEFT_CENTER,
                    &label,
                    FontId::proportional(13.0),
                    t.text,
                );
                let _ = resp.clone().on_hover_text(&tip);
            });
            app.place("document_tab", tab);
        }
        ui.painter().text(
            pos2(lay.version_right, r.center().y),
            Align2::RIGHT_CENTER,
            version,
            FontId::proportional(11.0),
            Color32::from_rgb(160, 168, 182),
        );
        if app.custom_titlebar {
            crate::titlebar::caption_buttons(app, ui, r);
        }
        match clicked {
            Some(Quick::File) => app.toggle_menu("file"),
            Some(Quick::Help) => app.toggle_menu("help"),
            Some(Quick::Save) => {
                app.save();
            }
            Some(Quick::Undo) => {
                let _ = app.run("edit.undo", json!({}));
            }
            Some(Quick::Redo) => {
                let _ = app.run("edit.redo", json!({}));
            }
            Some(Quick::Solve) => app.solve_all(),
            Some(Quick::Search) => {
                app.ui.palette = !app.ui.palette;
                app.ui.menu = None;
            }
            Some(Quick::Theme) => app.ui.dark = !app.ui.dark,
            None => {}
        }
        for (q, rect) in &anchors {
            app.place(q.region(), *rect);
        }
        let anchor = |q: Quick| anchors.iter().find(|(k, _)| *k == q).map_or(r, |(_, a)| *a);
        crate::frame::menu(app, ui, "file", anchor(Quick::File), crate::files::file_menu);
        crate::frame::menu(app, ui, "help", anchor(Quick::Help), help_menu);
    });
}

/// One quick action: its icon on the dark bar, dimmed when it cannot run now. `true` when clicked.
#[allow(clippy::too_many_arguments)]
fn quick_button(
    app: &mut ModelApp,
    ui: &mut egui::Ui,
    style: &IconButtonStyle,
    key: &str,
    icon: &'static str,
    rect: Rect,
    tip: &str,
    enabled: bool,
    selected: bool,
) -> bool {
    let resp = ui
        .scope_builder(egui::UiBuilder::new().max_rect(rect).id_salt(("fc_quick", key)), |ui| {
            ui.add_enabled_ui(enabled, |ui| {
                IconButton::new(tip, rect.size(), style).selected(selected).show(ui, |_, _, _| {}).on_hover_text(tip).on_disabled_hover_text(tip)
            })
            .inner
        })
        .inner;
    let tint = if enabled { Color32::WHITE } else { Color32::from_white_alpha(90) };
    icons::paint_tinted(app, ui, icon, rect.shrink(4.0), Look::AppBar, tint);
    enabled && resp.clicked()
}

/// The shortcut text with this platform's command key.
pub fn shortcut(key: &str) -> String {
    format!("{}+{key}", if cfg!(target_os = "macos") { "Cmd" } else { "Ctrl" })
}

/// Help: the samples (ready to solve), how to navigate, keys, and the version.
fn help_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    ui.set_max_width(340.0);
    crate::frame::menu_heading(ui, &t, "SAMPLES");
    let mut open = None;
    crate::samples::rows(ui, &t, &mut open);
    ui.separator();
    crate::frame::menu_heading(ui, &t, "NAVIGATION");
    ui.label("Left drag orbits (or pans or zooms: the navigation bar), right or middle drag pans, the wheel zooms, F or a double click fits.");
    ui.label(format!("S or {} searches every command; the command line at the bottom runs `id {{json}}`.", shortcut("K")));
    ui.separator();
    let report = ui
        .button("Report a Bug…")
        .on_hover_text("Draft a GitHub issue with the version, the OS and the last error; nothing is sent until you submit it")
        .clicked();
    let (mut check, mut reinstall) = (false, false);
    if app.updates_available() {
        check =
            ui.button("Check for Updates").on_hover_text("Ask GitHub Releases whether a newer FieldCraft is published (app.update.check)").clicked();
        reinstall = ui
            .button("Reinstall FieldCraft…")
            .on_hover_text("Download this version again, verify it and swap it in, to repair a damaged install")
            .clicked();
    }
    let crashes = ui.button("Crash Reports…").on_hover_text("The crash reports saved on this computer: report, copy or delete them").clicked();
    let settings = ui.button("Settings…").on_hover_text("Privacy and updates").clicked();
    ui.label(egui::RichText::new(format!("FieldCraft {}", crate::VERSION)).color(t.text_dim));
    if report {
        app.report_bug(None);
    } else if crashes {
        app.open_crash_reports();
    } else if check {
        app.check_for_updates();
    } else if reinstall {
        app.reinstall();
    } else if settings {
        app.open_settings();
    }
    if let Some(n) = open {
        app.ui.menu = None;
        app.open_sample(n);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Solve and the tools sit at the right end; the left group starts after the traffic lights,
    /// at the very left in full screen (where there are none).
    #[test]
    fn the_title_bar_keeps_solve_right_and_drops_the_traffic_light_inset_in_full_screen() {
        let bar = Rect::from_min_size(egui::Pos2::ZERO, vec2(1400.0, HEIGHT));
        let windowed = bar_layout(bar, true, 0.0, 30.0);
        let full = bar_layout(bar, false, 0.0, 30.0);
        assert_eq!(windowed.mark.x - full.mark.x, crate::titlebar::TRAFFIC_LIGHTS);
        assert!(full.mark.x < 30.0 && full.quick[0].left() < 150.0, "{full:?}");
        for l in [windowed, full] {
            // Left to right: mark, quick actions, (tab), Solve, search, theme, help, version.
            assert!(l.quick.windows(2).all(|w| w[0].right() <= w[1].left()));
            assert!(l.quick[3].right() < l.tab_room.0 && l.tab_room.1 < l.solve.left());
            assert!(l.solve.right() <= l.tools[0].left() && l.tools.windows(2).all(|w| w[0].right() <= w[1].left()));
            assert!(l.tools[2].right() < l.version_right - 30.0 && l.version_right <= 1400.0);
            assert!(l.solve.left() > 1100.0, "Solve is at the right end: {:?}", l.solve);
        }
        // The right group does not move with full screen; custom caption buttons push it left.
        assert_eq!(windowed.solve, full.solve);
        let captions = bar_layout(bar, false, crate::titlebar::WIDTH, 30.0);
        assert!((full.tools[2].right() - captions.tools[2].right() - crate::titlebar::WIDTH).abs() < 0.01);
        // Narrow: compact, still in order.
        let narrow = bar_layout(Rect::from_min_size(egui::Pos2::ZERO, vec2(700.0, HEIGHT)), true, 0.0, 30.0);
        assert!(narrow.compact && narrow.name.is_none() && narrow.quick[3].right() < narrow.solve.left());
    }
}
