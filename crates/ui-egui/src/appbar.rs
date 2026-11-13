//! The application bar: the window's title bar, drawn by the kit's application bar
//! (`solidcraft_ui_kit::app_bar`) as SolidCraft's is: FieldCraft's mark and name, Home (the
//! samples), File and its actions (save, undo, redo), a quick-access button (the theme), the
//! model's document tab, Solve, Settings, Help, Report a bug and the version. The command search
//! is at the end of the ribbon's tab row, SolidCraft's default place for it ([`crate::toolbar`]).
//! The kit's window chrome (`solidcraft_ui_kit::window_chrome`, [`crate::titlebar`]) makes its
//! empty space move the window and draws the caption buttons on Windows and Linux; macOS keeps
//! its traffic lights over the bar's left end.
//!
//! What stays FieldCraft's: its mark and Solve (the one prominent action, in the place the kit
//! leaves for an offer pill). The document tab is the kit's tab strip holding the window's one
//! model, where SolidCraft's tabs start.

use egui::{Color32, FontId, Rect, Sense, Stroke, pos2, vec2};
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
    search: "ribbon:search",
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
        search: false,
        search_shortcut: None,
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
    // An update on offer: SolidCraft's Update pill at the offer span's right end, before the
    // right cluster; Solve goes left of it.
    let k = solidcraft_theme::Tokens::get();
    let mut offer = out.offer;
    if let Some(pill) = update_pill(app, ui, offer, out.offer_labeled) {
        offer.max.x = pill.left() - 8.0;
        app.place("app_bar:update", pill);
    }
    // Solve: right-aligned where the kit leaves room for an offer pill, before the right cluster.
    let solve_w = kit::icon_label_width(ui, if out.offer_labeled { SOLVE } else { "" }, Some(SOLVE_ICON), true);
    let solve = Rect::from_center_size(pos2(offer.right() - solve_w / 2.0, cy), vec2(solve_w, k.control_height));
    let tabs_end = if offer.width() >= solve_w {
        solve_button(app, ui, solve, out.offer_labeled);
        app.place("app_bar:solve", solve);
        solve.left() - 8.0
    } else {
        out.controls_left - 8.0
    };
    document_tab(app, ui, r, out.tabs_start, tabs_end);
    let home_at = solidcraft_ui_registry::rect_of("app_bar:home").unwrap_or(out.file);
    let help_at = solidcraft_ui_registry::rect_of("app_bar:help").unwrap_or(r);
    crate::frame::menu(app, ui, "file", out.file, crate::files::file_menu);
    crate::frame::wide_menu(app, ui, "home", home_at, samples_menu);
    crate::frame::wide_menu(app, ui, "help", help_at, help_menu);
    if let Some(pill) = app.layout.get("app_bar:update").copied() {
        crate::frame::menu(app, ui, "update", pill, crate::updates::offer_menu);
    }
}

/// The Update pill's height and the room it leaves the document tab before it shortens its
/// label (SolidCraft's `PILL_HEIGHT`, `update::TAB_ROOM`).
const PILL_HEIGHT: f32 = 22.0;
const PILL_TAB_ROOM: f32 = 300.0;

/// SolidCraft's Update pill (`update::button`, `toolbar::offer_pill`): an accent pill at the right
/// end of `span` reading "Update to X", or "Update" when the document tab would be squeezed, or
/// its icon alone in a narrow bar. A click opens the update menu ([`crate::updates::offer_menu`]).
/// None when nothing is on offer or there is no room.
fn update_pill(app: &mut ModelApp, ui: &mut egui::Ui, span: Rect, labeled: bool) -> Option<Rect> {
    use crate::updates::Offer;
    let offer = crate::updates::offer(app)?;
    let (full, tip) = match &offer {
        Offer::Release(version) => (format!("Update to {version}"), format!("FieldCraft {version} is available")),
        Offer::Busy(busy) => ("Update".to_string(), busy.clone()),
        Offer::Reload => ("Reload to update".to_string(), "A new version of FieldCraft is available".to_string()),
    };
    let k = solidcraft_theme::Tokens::get();
    let font = FontId::proportional(k.font_label);
    let (icon_x, text_x, right_pad) = (16.0, 29.0, 11.0);
    let labels = [full.as_str(), "Update"];
    let label = labels.iter().enumerate().filter(|_| labeled).find_map(|(i, text)| {
        let width = text_x + ui.painter().layout_no_wrap(text.to_string(), font.clone(), k.app_bar_text).size().x + right_pad;
        let room = if i + 1 < labels.len() { span.width() - PILL_TAB_ROOM } else { span.width() };
        (width <= room).then_some((*text, width))
    });
    let width = label.map_or(PILL_HEIGHT + 4.0, |(_, w)| w);
    if width > span.width() {
        return None;
    }
    let rect = Rect::from_min_size(pos2(span.right() - width, span.center().y - PILL_HEIGHT / 2.0), vec2(width, PILL_HEIGHT));
    let resp = ui.interact(rect, ui.scope_id().with("fc_update_pill"), Sense::click());
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, label.map_or("Update", |(t, _)| t)));
    let hot = resp.hovered() || resp.has_focus();
    let radius = PILL_HEIGHT / 2.0 * f32::from(k.radius_xl) / 10.0;
    ui.painter().rect_filled(rect, radius, if hot { k.app_bar_hover_button } else { k.accent });
    let glyph_x = if label.is_some() { rect.left() + icon_x } else { rect.center().x };
    icons::paint(ui, crate::updates::UPDATE_ICON, Rect::from_center_size(pos2(glyph_x, rect.center().y), vec2(16.0, 16.0)), Look::AppBar);
    if let Some((text, _)) = label {
        ui.painter().text(pos2(rect.left() + text_x, rect.center().y), egui::Align2::LEFT_CENTER, text, font, k.app_bar_text);
    }
    if resp.backing_on_hover_text(tip).clicked() {
        app.toggle_menu("update");
    }
    Some(rect)
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

/// The model's document tab: the kit's tab strip (SolidCraft's `documents::tabs`) with the one
/// tab a FieldCraft window holds, from the documents separator at `x0` up to `x1`. Its × closes
/// the model and its + starts a new one: both are File › New (`file.new`), which asks first
/// when that would discard unsaved changes.
fn document_tab(app: &mut ModelApp, ui: &mut egui::Ui, r: Rect, x0: f32, x1: f32) {
    let migrated = app.session.migrated_from.map(|n| format!(" (from format {n})")).unwrap_or_default();
    let name = format!("{}{migrated}", app.model_name());
    let dirty = app.session.is_dirty();
    if x1 - x0 < TAB_STRIP_ROOM + TAB_MIN {
        return;
    }
    let mut place = app.session.path.clone().unwrap_or_else(|| "Not saved yet".into());
    if let Some(n) = app.session.migrated_from {
        place.push_str(&format!("\nOpened from a format {n} file: the first save keeps that file as a copy beside it"));
        for note in &app.session.migration_notes {
            place.push_str(&format!("\n{note}"));
        }
    }
    let close_key = shortcut("W");
    let menu = [app_bar::TabMenuRow { id: "close", label: "Close", key_on_selected: &close_key, needs_others: false }];
    let mut strip = app_bar::TabStrip {
        tabs: vec![app_bar::DocTab { name: &name, list_name: &name, dirty, place: &place }],
        selected: 0,
        shown_selected: Some(0),
        title_bar: true,
        double_click_renames: false,
        rename: None,
        menu: &menu,
        close_name: "Close the model",
        close_tip: "Close the model: a new, empty one takes its place (File › New)",
        more_tip: "Every model",
        new_tip: "New model (File › New)",
        id: egui::Id::unique("fc_document_tabs"),
        handles: TAB_HANDLES,
    };
    // The strip's hit regions stay clear of what is right of it.
    let room = Rect::from_min_max(pos2(x0, r.top()), pos2(x1, r.bottom()));
    let mut tabs_ui = ui.new_child(egui::UiBuilder::new().max_rect(room));
    tabs_ui.set_clip_rect(room.intersect(ui.clip_rect()));
    let actions = app_bar::tabs(&mut tabs_ui, &mut strip, r, x0, x1);
    if let Some(tab) = solidcraft_ui_registry::rect_of("document_tab:0") {
        app.place("document_tab", tab);
    }
    for action in actions {
        match action {
            app_bar::TabAction::Close(_) | app_bar::TabAction::New => app.new_model(),
            app_bar::TabAction::Menu { id, .. } if id == "close" => app.new_model(),
            app_bar::TabAction::Maximize => window_chrome::toggle_maximize(ui.ctx()),
            _ => {}
        }
    }
}

/// Room the kit's tab strip keeps beside its tabs (its + button), and its narrowest tab.
const TAB_STRIP_ROOM: f32 = 34.0;
const TAB_MIN: f32 = 72.0;

/// Where the tab strip publishes its parts in the widget registry.
const TAB_HANDLES: app_bar::TabHandles<'static> = app_bar::TabHandles {
    tab: "document_tab:",
    close: "document_tab_close:",
    more: "document_tabs_more",
    new: "document_tab_new",
    rename: "document_tab_rename",
    menu_row: "document_tab_menu:",
    more_row: "document_tabs_more:",
};

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
    // The samples are Home's (and File › Samples…): listed here too, they pushed Report a Bug,
    // Crash Reports, Settings and the version below a laptop window's edge.
    crate::frame::menu_heading(ui, "Navigation");
    let note = |ui: &mut egui::Ui, text: &str| {
        ui.horizontal(|ui| {
            ui.add_space(solidcraft_theme::Tokens::get().space_8);
            crate::files::note(ui, text);
        });
    };
    note(ui, "Left drag orbits (or pans or zooms: the navigation bar), right or middle drag pans, the wheel zooms, F or a double click fits.");
    note(ui, &format!("S or {} searches every command; `id {{json}}` typed there runs it.", shortcut("K")));
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
    // Everything the menu drew, down to the version (control channel, tests).
    app.place("menu:help", ui.min_rect());
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
        let mut app = ModelApp::new(Session::default());
        app.integrated_titlebar = integrated;
        let ctx = egui::Context::default();
        for _ in 0..2 {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(width, 600.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app_bar(&mut app, ui));
            out.textures_delta.clear();
        }
        app.layout
    }

    /// Left to right as in SolidCraft: Home, File, its actions, the theme, the document tab, Solve,
    /// Settings, Help and Report a bug; the left group starts after the traffic
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

    /// An update on offer is SolidCraft's Update pill before the right cluster, Solve left of it;
    /// with none on offer there is no pill.
    #[test]
    fn an_update_on_offer_is_a_pill_before_settings() {
        let draw = |reload: bool| {
            let mut app = ModelApp::new(Session::default());
            if reload {
                app.reload = Some(Box::new(|| {}));
            }
            let ctx = egui::Context::default();
            for _ in 0..2 {
                let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1400.0, 600.0))), ..Default::default() };
                let mut out = ctx.run_ui(raw, |ui| app_bar(&mut app, ui));
                out.textures_delta.clear();
            }
            app.layout
        };
        assert!(!draw(false).contains_key("app_bar:update"));
        let m = draw(true);
        let at = |k: &str| m.get(k).copied().unwrap_or_else(|| panic!("{k} not drawn"));
        let (solve, update, settings) = (at("app_bar:solve"), at("app_bar:update"), at("app_bar:settings"));
        assert!(solve.right() < update.left() && update.right() < settings.left(), "{solve:?} {update:?} {settings:?}");
    }

    /// The document tab is the kit's tab strip, as SolidCraft's: it starts just after the
    /// documents separator, right of the quick-access buttons, instead of centred on the bar, and
    /// its + follows it.
    #[test]
    fn the_document_tab_starts_the_tab_strip_after_the_actions() {
        let m = regions(1400.0, false);
        let at = |k: &str| m.get(k).copied().unwrap_or_else(|| panic!("{k} not drawn"));
        let (theme, tab) = (at("app_bar:theme"), at("document_tab"));
        assert!(tab.left() > theme.right() && tab.left() - theme.right() < 24.0, "theme {theme:?}, tab {tab:?}");
        assert!(tab.center().x < 700.0, "left-aligned, not centred: {tab:?}");
        let new = solidcraft_ui_registry::handle_at("document_tab_new").unwrap_or_else(|| panic!("no + button"));
        assert!(new.x > tab.right() && new.x < tab.right() + 24.0, "+ at {new:?} after {tab:?}");
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

    /// The Help menu fits a laptop window (1280 × 800): every row, down to the version, can be
    /// reached.
    #[test]
    fn the_help_menu_fits_the_window() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        let screen = Rect::from_min_size(egui::Pos2::ZERO, vec2(1280.0, 800.0));
        for _ in 0..3 {
            app.ui.menu = Some("help".into());
            let raw = egui::RawInput { screen_rect: Some(screen), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app_bar(&mut app, ui));
            out.textures_delta.clear();
        }
        let menu = app.layout.get("menu:help").copied().unwrap_or_else(|| panic!("Help menu not drawn"));
        assert!(menu.bottom() <= screen.bottom(), "{menu:?}");
    }
}
