//! The application bar: the window's title bar, drawn by the kit's application bar
//! (`solidcraft_ui_kit::app_bar`) as SolidCraft's is: FieldCraft's mark and name, Home (the
//! samples), File and its actions (save, undo, redo), a quick-access button (the theme), the
//! model's document tab, Solve, the command search, Settings, Help, Report a bug and the version.
//! The kit's window chrome (`solidcraft_ui_kit::window_chrome`, [`crate::titlebar`]) makes its
//! empty space move the window and draws the caption buttons on Windows and Linux; macOS keeps
//! its traffic lights over the bar's left end.
//!
//! What stays FieldCraft's: its mark, Solve (the one prominent action, in the place the kit
//! leaves for an offer pill) and the single document tab (drawn in the kit's tab shape; the kit's
//! tab strip closes, adds and reorders documents, which a FieldCraft window does not).

use egui::{Align2, Color32, FontId, Rect, Sense, Stroke, pos2, vec2};
use serde_json::json;
use solidcraft_ui_kit::app_bar::{self, BarAction, BarButton};
use solidcraft_ui_kit::window_chrome;
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _};

use crate::ModelApp;
use crate::icons::{self, Look};
use crate::theme::Tokens;

/// The bar's height.
pub const HEIGHT: f32 = app_bar::HEIGHT;

/// Where the kit's bar publishes its parts in the widget registry, under the names `ui.inspect`'s
/// layout gives them.
const HANDLES: app_bar::Handles<'static> = app_bar::Handles {
    home: &["app_bar:home"],
    home_rect: &["app_bar:home"],
    file: &["app_bar:file"],
    file_rect: &["app_bar:file"],
    button: &["app_bar:"],
    separator: "app_bar:separator:",
    search: "app_bar:search",
    report: "app_bar:report",
    settings: "app_bar:settings",
    version: "app_bar:version",
};

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

/// Whether the window is in full screen (no traffic lights to make room for).
fn fullscreen(ctx: &egui::Context) -> bool {
    ctx.input(|i| i.viewport().fullscreen.unwrap_or(false))
}

/// This version's release notes.
fn release_notes() -> String {
    format!("https://github.com/tangency-labs/fieldcraft/releases/tag/v{}", crate::VERSION)
}

/// Solve's label and icon: the kit's primary icon-label action.
const SOLVE: &str = "Solve";
const SOLVE_ICON: &str = "play";

pub fn app_bar(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    // The title bar is an application-bar surface: flat and full width whatever the layout.
    let host = kit::Container::new(ui.ctx(), kit::Surface::ApplicationBar);
    let panel = host.edge("fc_appbar", kit::Edge::Top).exact_size(HEIGHT).show_separator_line(false);
    host.show_edge(ui, panel, host.edge_frame(t.app_bar), |ui| {
        let r = ui.max_rect();
        app.place("app_bar", r);
        let chrome =
            window_chrome::TitleBar { traffic_lights: app.integrated_titlebar && !fullscreen(ui.ctx()), caption_buttons: app.custom_titlebar };
        let action = chrome.show(ui, r, |ui, content| bar(app, ui, r, content));
        if let Some(action) = action {
            crate::titlebar::run(app, ui.ctx(), action);
        }
    });
}

/// The bar's contents within `content`, the room the window's own buttons leave in `r`.
fn bar(app: &mut ModelApp, ui: &mut egui::Ui, r: Rect, content: Rect) {
    let (undo, redo) = app.session.undo_depth();
    let busy = app.busy().is_some();
    let (save_tip, undo_tip, redo_tip) =
        (format!("Save ({})", shortcut("S")), format!("Undo ({})", shortcut("Z")), format!("Redo ({})", shortcut("Shift+Z")));
    let search_key = shortcut("K");
    let version = format!("v{}", crate::VERSION);
    let url = release_notes();
    let button = |id, icon, tip, enabled, wide_only| BarButton { id, icon, tip, enabled, narrow: true, wide_only };
    let chrome = |id, icon, tip| BarButton { id, icon, tip, enabled: true, narrow: false, wide_only: false };
    let data = app_bar::AppBar {
        mark: &mark,
        name: "FieldCraft",
        home: chrome("home", "home", "Samples: models ready to solve"),
        file: "File",
        file_open: app.ui.menu.as_deref() == Some("file"),
        actions: vec![
            button("save", "save", &save_tip, !busy, false),
            button("undo", "undo", &undo_tip, undo > 0 && !busy, true),
            button("redo", "redo", &redo_tip, redo > 0 && !busy, true),
        ],
        quick: 1,
        help: chrome("help", "help", "Help: navigation, keys, updates, about"),
        settings: chrome("settings", "settings", "Settings: privacy and updates"),
        report: app_bar::Report {
            label: "Report a bug",
            tip: "Report a bug: draft a GitHub issue with the version, the OS and the last error; nothing is sent until you submit it",
            icon_only: false,
        },
        search: true,
        search_shortcut: Some(&search_key),
        assistant: false,
        version: Some(app_bar::Version { text: &version, url: &url }),
        handles: HANDLES,
    };
    let out = app_bar::show(ui, &data, r, content);
    // The kit's buttons, where it published them this frame (`ui.inspect`'s layout).
    for region in
        ["app_bar:home", "app_bar:file", "app_bar:save", "app_bar:undo", "app_bar:redo", "app_bar:help", "app_bar:settings", "app_bar:report"]
    {
        if let Some(rect) = solidcraft_ui_registry::rect_of(region) {
            app.place(region, rect);
        }
    }
    for action in &out.actions {
        match action {
            BarAction::Home => app.toggle_menu("home"),
            BarAction::File { .. } => {
                app.ui.palette = false;
                app.toggle_menu("file");
            }
            BarAction::Button { id } => match id.as_str() {
                "save" => {
                    app.save();
                }
                "undo" => {
                    let _ = app.run("edit.undo", json!({}));
                }
                "redo" => {
                    let _ = app.run("edit.redo", json!({}));
                }
                "settings" => app.open_settings(),
                _ => {}
            },
            BarAction::Help { .. } => app.toggle_menu("help"),
            BarAction::Report => app.report_bug(None),
            BarAction::Version => ui.ctx().open_url(egui::OpenUrl::new_tab(&url)),
        }
    }
    let cy = r.center().y;
    // The quick-access button: the theme.
    if let Some((first, _)) = out.quick {
        let rect = Rect::from_center_size(pos2(first, cy), vec2(32.0, 32.0));
        let (icon, tip) = if app.ui.dark { ("ui.sun", "Light theme") } else { ("ui.moon", "Dark theme") };
        if quick_button(ui, "theme", icon, rect, tip) {
            app.ui.dark = !app.ui.dark;
        }
        app.place("app_bar:theme", rect);
    }
    // The closed command search: a click opens the command palette.
    if let Some(search) = out.search {
        search_field(app, ui, search, &search_key);
        app.place("app_bar:search", search);
    }
    // Solve: right-aligned where the kit leaves room for an offer pill, before the right cluster.
    let k = solidcraft_theme::Tokens::get();
    let solve_w = kit::icon_label_width(ui, if out.offer_labeled { SOLVE } else { "" }, Some(SOLVE_ICON), true);
    let solve = Rect::from_center_size(pos2(out.offer.right() - solve_w / 2.0, cy), vec2(solve_w, k.control_height));
    let tabs_end = if out.offer.width() >= solve_w {
        solve_button(app, ui, solve, out.offer_labeled);
        app.place("app_bar:solve", solve);
        solve.left() - 8.0
    } else {
        out.controls_left - 8.0
    };
    // The open search field covers the bar from its right edge leftward, as SolidCraft's does:
    // the document tab ends before it (and is left out when no room remains).
    let tabs_end = match out.search.filter(|_| app.ui.palette) {
        Some(field) => tabs_end.min(crate::palette::open_field_at(ui.ctx().content_rect(), field).x - 8.0),
        None => tabs_end,
    };
    document_tab(app, ui, r, out.tabs_start, tabs_end);
    let home_at = solidcraft_ui_registry::rect_of("app_bar:home").unwrap_or(out.file);
    let help_at = solidcraft_ui_registry::rect_of("app_bar:help").unwrap_or(r);
    crate::frame::menu(app, ui, "file", out.file, crate::files::file_menu);
    crate::frame::wide_menu(app, ui, "home", home_at, samples_menu);
    crate::frame::wide_menu(app, ui, "help", help_at, help_menu);
}

/// The theme's quick-access button, drawn as the kit's bar draws its icon buttons. `true` when
/// clicked.
fn quick_button(ui: &mut egui::Ui, key: &str, icon: &'static str, rect: Rect, tip: &str) -> bool {
    let k = solidcraft_theme::Tokens::get();
    let resp = ui.interact(rect, ui.scope_id().with(("fc_quick", key)), Sense::click());
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, tip));
    if resp.hovered() || resp.has_focus() {
        ui.painter().rect_filled(rect.shrink(3.0), k.radius_sm, k.app_bar_hover_button);
    }
    icons::paint(ui, icon, Rect::from_center_size(rect.center(), vec2(18.0, 18.0)), Look::AppBar);
    resp.backing_on_hover_text(tip).clicked()
}

/// The command search's closed field (the kit's), opening the command palette when clicked.
fn search_field(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect, key: &str) {
    // While the search is open its own field is drawn here, over the bar (`palette::popup`): the
    // kit's one search editor, drawn once.
    if app.ui.palette {
        return;
    }
    let geometry = solidcraft_ui_kit::SearchFieldGeometry {
        height: rect.height(),
        radius: solidcraft_ui_kit::CommandSearch::field_radius(0.0, false),
        progress: 0.0,
    };
    let mut field = ui.new_child(egui::UiBuilder::new().max_rect(rect).id_salt("fc_app_bar_search"));
    let mut empty = String::new();
    solidcraft_ui_kit::CommandSearch::field(&mut field, &mut empty, false, false, false, geometry, Some(key), 0);
    let resp = ui.interact(rect, ui.scope_id().with("fc_app_bar_search_open"), Sense::click());
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, "Search commands"));
    if resp.backing_on_hover_text(format!("Search commands (S or {key})")).clicked() {
        app.ui.palette = !app.ui.palette;
        app.ui.menu = None;
    }
}

/// Solve: the one prominent action, the kit's primary icon-label button (with its label when
/// `labeled`). It runs `solve.run`: greyed (the kit's disabled look) while the model misses what
/// it needs, every missing piece listed in the tooltip.
fn solve_button(app: &mut ModelApp, ui: &mut egui::Ui, solve: Rect, labeled: bool) {
    let missing = app.unavailable("solve.run");
    let can = missing.is_none() && app.busy().is_none();
    let head = if missing.is_some() { "Solve: not ready yet (or try Home › Samples)" } else { "Solve every step and show the results" };
    let tip = fieldcraft_engine::requirement::Unavailable::tip(head, missing.as_ref());
    let mut slot = ui.new_child(egui::UiBuilder::new().max_rect(solve).id_salt("fc_solve").layout(egui::Layout::left_to_right(egui::Align::Center)));
    let resp = kit::icon_label_action(&mut slot, if labeled { SOLVE } else { "" }, Some(SOLVE_ICON), can, true);
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, can, SOLVE));
    let resp = if can { resp.backing_on_hover_text(tip) } else { resp.backing_on_disabled_hover_text(tip) };
    if resp.clicked() && can {
        app.solve_all();
    }
}

/// The model's document tab, centred in the room from `x0` to `x1`, in the kit's tab shape.
fn document_tab(app: &mut ModelApp, ui: &mut egui::Ui, r: Rect, x0: f32, x1: f32) {
    let k = solidcraft_theme::Tokens::get();
    let migrated = app.session.migrated_from.map(|n| format!(" (from format {n})")).unwrap_or_default();
    let label = format!("{}{migrated}{}", app.model_name(), if app.session.is_dirty() { " •" } else { "" });
    let font = FontId::proportional(k.font_button);
    let text_w = crate::elide::width(ui, &label, &font);
    let tab_w = (text_w + 40.0).clamp(120.0, 260.0).min((x1 - x0).max(0.0));
    if tab_w <= 40.0 {
        return;
    }
    let tab_x = (r.center().x - tab_w / 2.0).clamp(x0, (x1 - tab_w).max(x0));
    let tab = Rect::from_min_size(pos2(tab_x, r.center().y - 13.0), vec2(tab_w, 26.0));
    let mut place = app.session.path.clone().unwrap_or_else(|| "Not saved yet".into());
    if let Some(n) = app.session.migrated_from {
        place.push_str(&format!("\nOpened from a format {n} file: the first save keeps that file as a copy beside it"));
        for note in &app.session.migration_notes {
            place.push_str(&format!("\n{note}"));
        }
    }
    let resp = ui.interact(tab, egui::Id::unique("fc_document_tab"), Sense::hover());
    resp.widget_info(|| egui::WidgetInfo::selected(egui::Role::Tab, true, true, &label));
    app_bar::tab_shape(ui.painter(), &k, tab, true, resp.hovered());
    let shown = crate::elide::middle_in(ui, &label, &font, tab_w - 28.0);
    let cut = shown != label.as_str();
    ui.painter().text(pos2(tab.left() + 14.0, tab.center().y), Align2::LEFT_CENTER, shown, font, k.app_bar_text);
    let _ = resp.backing_on_hover_text(app_bar::tab_tip(&label, &place, cut));
    app.place("document_tab", tab);
}

/// The shortcut text with this platform's command key.
pub fn shortcut(key: &str) -> String {
    format!("{}+{key}", if cfg!(target_os = "macos") { "Cmd" } else { "Ctrl" })
}

/// Home: the samples, ready to solve.
fn samples_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    crate::frame::menu_heading(ui, "Samples");
    let mut open = None;
    crate::samples::rows(ui, &t, &mut open);
    if let Some(n) = open {
        app.ui.menu = None;
        app.open_sample(n);
    }
}

/// Help: the samples (ready to solve), how to navigate, keys, updates, crash reports and the version.
fn help_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    samples_menu(app, ui);
    crate::frame::menu_separator(ui);
    crate::frame::menu_heading(ui, "Navigation");
    let note = |ui: &mut egui::Ui, text: &str| {
        ui.horizontal(|ui| {
            ui.add_space(solidcraft_theme::Tokens::get().space_8);
            crate::files::note(ui, text);
        });
    };
    note(ui, "Left drag orbits (or pans or zooms: the navigation bar), right or middle drag pans, the wheel zooms, F or a double click fits.");
    note(ui, &format!("S or {} searches every command; the command line at the bottom runs `id {{json}}`.", shortcut("K")));
    crate::frame::menu_separator(ui);
    let row =
        |ui: &mut egui::Ui, id: &str, label: &str, tip: &str| crate::frame::menu_row(ui, id, label, None, true).backing_on_hover_text(tip).clicked();
    let report = row(
        ui,
        "help.report_bug",
        "Report a Bug…",
        "Draft a GitHub issue with the version, the OS and the last error; nothing is sent until you submit it",
    );
    let (mut check, mut reinstall) = (false, false);
    if app.updates_available() {
        check = row(ui, "help.check_updates", "Check for Updates", "Ask GitHub Releases whether a newer FieldCraft is published (app.update.check)");
        reinstall =
            row(ui, "help.reinstall", "Reinstall FieldCraft…", "Download this version again, verify it and swap it in, to repair a damaged install");
    }
    let crashes = row(ui, "help.crash_reports", "Crash Reports…", "The crash reports saved on this computer: report, copy or delete them");
    let settings = row(ui, "help.settings", "Settings…", "Privacy and updates");
    crate::frame::menu_separator(ui);
    note(ui, &format!("FieldCraft {}", crate::VERSION));
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
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;

    use super::*;

    /// The bar's regions after a frame of a window `width` wide.
    fn regions(width: f32, integrated: bool) -> std::collections::BTreeMap<&'static str, Rect> {
        regions_with(width, integrated, false)
    }

    /// [`regions`], with the command search open or not.
    fn regions_with(width: f32, integrated: bool, search: bool) -> std::collections::BTreeMap<&'static str, Rect> {
        let mut app = ModelApp::new(Session::default());
        app.integrated_titlebar = integrated;
        app.ui.palette = search;
        let ctx = egui::Context::default();
        for _ in 0..2 {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(width, 600.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app_bar(&mut app, ui));
            out.textures_delta.clear();
        }
        app.layout
    }

    /// The open command search covers the bar from its field leftward; the document tab ends
    /// before it, so no part of the model's name shows beside the field.
    #[test]
    fn the_open_search_leaves_no_part_of_the_document_tab_beside_it() {
        for width in [1440.0, 1100.0] {
            let open = regions_with(width, false, true);
            let field = open.get("app_bar:search").copied().unwrap_or(Rect::NOTHING);
            let left = crate::palette::open_field_at(Rect::from_min_size(egui::Pos2::ZERO, vec2(width, 600.0)), field).x;
            if let Some(tab) = open.get("document_tab") {
                assert!(tab.right() <= left, "{width}: the tab {tab:?} runs under the open field from {left}");
            }
            assert!(regions(width, false).contains_key("document_tab"), "{width}: the tab shows while the search is closed");
        }
    }

    /// Left to right as in SolidCraft: Home, File, its actions, the theme, the document tab, Solve,
    /// the search, Settings, Help and Report a bug; the left group starts after the traffic
    /// lights, and Solve and the right cluster stay put.
    #[test]
    fn the_title_bar_is_laid_out_as_solidcrafts() {
        let at = |m: &std::collections::BTreeMap<&'static str, Rect>, k: &str| m.get(k).copied().unwrap_or_else(|| panic!("{k} not drawn"));
        let windowed = regions(1400.0, true);
        let full = regions(1400.0, false);
        assert!((at(&windowed, "app_bar:home").left() - at(&full, "app_bar:home").left() - window_chrome::TRAFFIC_LIGHTS).abs() < 0.5);
        assert_eq!(at(&windowed, "app_bar:solve"), at(&full, "app_bar:solve"));
        for m in [windowed, full] {
            let order = [
                "app_bar:home",
                "app_bar:file",
                "app_bar:save",
                "app_bar:undo",
                "app_bar:redo",
                "app_bar:theme",
                "document_tab",
                "app_bar:solve",
                "app_bar:search",
                "app_bar:settings",
                "app_bar:help",
                "app_bar:report",
            ];
            for pair in order.windows(2) {
                let (a, b) = (at(&m, pair[0]), at(&m, pair[1]));
                assert!(a.right() <= b.left() + 0.5, "{} {a:?} before {} {b:?}", pair[0], pair[1]);
            }
            assert!(at(&m, "app_bar:report").right() < 1400.0 && at(&m, "app_bar:solve").left() > 700.0);
        }
    }

    /// Solve is the kit's primary icon-label action, at the kit's control height (it was a
    /// hand-painted 24 pt rect with its own fill, mark and label font).
    #[test]
    fn solve_is_the_kits_primary_action() {
        let m = regions(1400.0, false);
        let solve = m.get("app_bar:solve").copied().unwrap_or_else(|| panic!("Solve not drawn"));
        let k = solidcraft_theme::Tokens::get();
        assert!((solve.height() - k.control_height).abs() < 0.5, "{solve:?}");
    }
}
