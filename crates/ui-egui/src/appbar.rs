//! The application bar: the window's title bar, drawn by the kit's application bar
//! (`solidcraft_ui_kit::app_bar`) as SolidCraft's is: FieldCraft's mark and name, Home (the
//! start page), File and its actions (save, undo, redo), a quick-access button (the theme), the
//! model's document tab, the community pills (Donate, Discord: [`crate::links`]), Settings, Help,
//! Report a bug, the version, and Solve. The command search
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
    pill: crate::links::PILL_PREFIX,
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
/// App-owned item identifiers shared between [`bar`] and [`fill_overflow_rows`].
const QUICK_THEME_ID: &str = "quick:0";
const OFFERED_ID: &str = "offered";
const UPDATE_ID: &str = "update";
/// Smallest width at which Solve still renders cleanly as a compact icon button on the bar.
const COMPACT_SOLVE_MIN_W: f32 = 28.0;
const COMPACT_SOLVE_MAX_W: f32 = 32.0;

pub fn app_bar(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    // The title bar is an application-bar surface: flat and full width whatever the layout.
    let host = kit::Container::new(ui.ctx(), kit::Surface::ApplicationBar);
    let panel = host.edge("fc_appbar", kit::Edge::Top).exact_size(HEIGHT).show_separator_line(false);
    host.show_edge(ui, panel, host.edge_frame(t.app_bar), |ui| {
        let r = ui.max_rect();
        app.place("app_bar", r);
        let chrome = window_chrome::TitleBar {
            traffic_lights: app.integrated_titlebar && !fullscreen(ui.ctx()),
            caption_buttons: app.custom_titlebar,
            hide_caption_buttons: false,
        };
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
        items: None,
        quick_ids: &[],
        mark: &mark,
        name: "FieldCraft",
        leading: None,
        home: chrome("home", "home", "Start: new, open, recent files and samples"),
        file: "File",
        file_open: app.ui.menu.as_deref() == Some("file"),
        actions: vec![
            button("save", "save", &save_tip, !busy, false),
            button("undo", "undo", &undo_tip, undo > 0 && !busy, true),
            button("redo", "redo", &redo_tip, redo > 0 && !busy, true),
        ],
        quick: 1,
        help: chrome("help", "help", "Help: navigation, keys, updates, about"),
        settings: chrome("settings", "settings", "Settings: units, the window, privacy and updates"),
        report: app_bar::Report {
            label: "Report a bug",
            tip: "Report a bug: draft a GitHub issue with the version, the OS and the last error; nothing is sent until you submit it",
            icon_only: false,
        },
        pills: crate::links::pills(&app.prefs, true),
        search: false,
        search_shortcut: None,
        assistant: false,
        version: Some(app_bar::Version { text: &version, url: &url }),
        handles: HANDLES,
    };
    let solve_full = kit::icon_label_width(ui, SOLVE, Some(SOLVE_ICON), true);
    let update_full = update_pill_label(app).map(|(full, _)| update_pill_width(ui, &full));
    let offered = (solve_full + update_full.map_or(0.0, |w| w + 8.0)).ceil();
    let mut out = app_bar::show_with_offers(ui, &data, r, content, app_bar::Offers { update: Some(offered), managed: None });
    // The kit's buttons, where it published them this frame (`ui.inspect`'s layout).
    for region in
        ["app_bar:home", "app_bar:file", "app_bar:save", "app_bar:undo", "app_bar:redo", "app_bar:help", "app_bar:settings", "app_bar:report"]
    {
        if let Some(rect) = solidcraft_ui_registry::rect_of(region) {
            app.place(region, rect);
        }
    }
    for link in crate::links::LINKS {
        for region in [link.region, link.hide_region] {
            if let Some(rect) = solidcraft_ui_registry::rect_of(region) {
                app.place(region, rect);
            }
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
    // Reserve Solve's slot before placing the Update pill so Update never starves Solve.
    let k = solidcraft_theme::Tokens::get();
    let mut offer = out.offer;
    let solve_w = kit::icon_label_width(ui, if out.offer_labeled { SOLVE } else { "" }, Some(SOLVE_ICON), true);
    let solve_reserve = if out.offer_labeled { solve_w } else { COMPACT_SOLVE_MIN_W };
    let update_span = if offer.width() >= solve_reserve + 8.0 {
        Rect::from_min_max(pos2(offer.left() + solve_reserve + 8.0, offer.top()), offer.right_bottom())
    } else {
        Rect::NOTHING
    };
    let update_drawn = if let Some(pill) = update_pill(app, ui, update_span, out.offer_labeled) {
        offer.max.x = pill.left() - 8.0;
        app.place("app_bar:update", pill);
        true
    } else {
        false
    };
    let (drawn_left, solve_drawn) = if offer.width() >= solve_w {
        let solve = Rect::from_center_size(pos2(offer.right() - solve_w / 2.0, cy), vec2(solve_w, k.control_height));
        solve_button(app, ui, solve, out.offer_labeled);
        app.place("app_bar:solve", solve);
        (solve.left(), true)
    } else if !out.offer_labeled && offer.width() >= COMPACT_SOLVE_MIN_W {
        // The kit's "offered" slot steps down to 32 px before overflowing; `icon_label_width`
        // for a primary button without text is 52 px, so paint a compact icon button instead of
        // letting Solve vanish from both the bar and the More menu.
        let w = offer.width().min(COMPACT_SOLVE_MAX_W);
        let solve = Rect::from_center_size(pos2(offer.right() - w / 2.0, cy), vec2(w, k.control_height));
        solve_compact_button(app, ui, solve);
        app.place("app_bar:solve", solve);
        (solve.left(), true)
    } else {
        (offer.right(), false)
    };
    fill_overflow_rows(app, &mut out.overflow_items, !solve_drawn, !update_drawn);
    if let Some(action) = app_bar::overflow_menu(ui.ctx(), &out.overflow_items, out.overflow_button) {
        out.actions.push(action);
    }
    for action in &out.actions {
        match action {
            BarAction::Home => app.open_start(),
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
                other => log::warn!("unhandled app bar button action: {other}"),
            },
            BarAction::Help { .. } => app.toggle_menu("help"),
            BarAction::Report => app.report_bug(None),
            BarAction::Version => ui.ctx().open_url(egui::OpenUrl::new_tab(&url)),
            BarAction::Pill { id } => {
                if let Some(link) = crate::links::by_id(id) {
                    app.open_link(ui.ctx(), link.url);
                }
            }
            BarAction::HidePill { id } => app.hide_link(id),
            BarAction::App { id } => match id.as_str() {
                QUICK_THEME_ID => app.ui.dark = !app.ui.dark,
                OFFERED_ID => {
                    if app.unavailable("solve.run").is_none() && app.busy().is_none() {
                        app.solve_all();
                    }
                }
                UPDATE_ID => app.toggle_menu("update"),
                other => log::warn!("unhandled app bar app action: {other}"),
            },
        }
    }
    // The document tab ends before the right cluster, or before Solve if it ever reaches left of it.
    let tabs_end = drawn_left.min(out.controls_left) - 8.0;
    document_tab(app, ui, r, out.tabs_start, tabs_end);
    let help_at = solidcraft_ui_registry::rect_of("app_bar:help").unwrap_or(r);
    crate::frame::menu(app, ui, "file", out.file, crate::files::file_menu);
    crate::frame::wide_menu(app, ui, "help", help_at, help_menu);
    if let Some(pill) = app.layout.get("app_bar:update").copied() {
        crate::frame::menu(app, ui, "update", pill, crate::updates::offer_menu);
    }
}

/// Fill labels, icons and availability for app-owned rows in the `More` overflow menu, ensuring
/// Solve and Update appear there whenever they could not be painted on the bar.
fn fill_overflow_rows(app: &ModelApp, rows: &mut Vec<app_bar::OverflowItem>, append_solve: bool, append_update: bool) {
    let solve_enabled = app.unavailable("solve.run").is_none() && app.busy().is_none();
    let update_label = update_pill_label(app);
    for row in rows.iter_mut() {
        match row.id.as_str() {
            QUICK_THEME_ID => {
                let (icon, label) = if app.ui.dark { ("ui.sun", "Light theme") } else { ("ui.moon", "Dark theme") };
                row.label = label.to_string();
                row.icon = icon.to_string();
                row.enabled = true;
            }
            OFFERED_ID => {
                row.label = SOLVE.to_string();
                row.icon = SOLVE_ICON.to_string();
                row.enabled = solve_enabled;
            }
            UPDATE_ID => {
                if let Some((full, _)) = &update_label {
                    row.label = full.clone();
                } else {
                    row.label = "Update".to_string();
                }
                row.icon = crate::updates::UPDATE_ICON.to_string();
                row.enabled = true;
            }
            other => {
                if !other.starts_with("action:") && !other.starts_with("pill:") {
                    log::warn!("unhandled app bar overflow item: {other}");
                }
            }
        }
    }
    if append_solve && !rows.iter().any(|r| r.id == OFFERED_ID) {
        rows.push(app_bar::OverflowItem {
            id: OFFERED_ID.to_string(),
            label: SOLVE.to_string(),
            icon: SOLVE_ICON.to_string(),
            enabled: solve_enabled,
            action: BarAction::App { id: OFFERED_ID.to_string() },
        });
    }
    if append_update
        && let Some((full, _)) = update_label
        && !rows.iter().any(|r| r.id == UPDATE_ID)
    {
        rows.push(app_bar::OverflowItem {
            id: UPDATE_ID.to_string(),
            label: full,
            icon: crate::updates::UPDATE_ICON.to_string(),
            enabled: true,
            action: BarAction::App { id: UPDATE_ID.to_string() },
        });
    }
}

/// The Update pill's height (SolidCraft's `PILL_HEIGHT`).
const PILL_HEIGHT: f32 = 22.0;

fn update_pill_label(app: &ModelApp) -> Option<(String, String)> {
    use crate::updates::Offer;
    let offer = crate::updates::offer(app)?;
    Some(match &offer {
        Offer::Release(version) => (format!("Update to {version}"), format!("FieldCraft {version} is available")),
        Offer::Busy(busy) => ("Update".to_string(), busy.clone()),
        Offer::Reload => ("Reload to update".to_string(), "A new version of FieldCraft is available".to_string()),
    })
}

fn update_pill_width(ui: &egui::Ui, text: &str) -> f32 {
    let k = solidcraft_theme::Tokens::get();
    let font = FontId::proportional(k.font_label);
    let (text_x, right_pad) = (29.0, 11.0);
    text_x + ui.painter().layout_no_wrap(text.to_string(), font, k.app_bar_text).size().x + right_pad
}

/// SolidCraft's Update pill (`update::button`, `toolbar::offer_pill`): an accent pill at the right
/// end of `span` reading "Update to X", or "Update", or its icon alone in a narrow bar. A click
/// opens the update menu ([`crate::updates::offer_menu`]). None when nothing is on offer or there
/// is no room.
fn update_pill(app: &mut ModelApp, ui: &mut egui::Ui, span: Rect, labeled: bool) -> Option<Rect> {
    let (full, tip) = update_pill_label(app)?;
    let k = solidcraft_theme::Tokens::get();
    let font = FontId::proportional(k.font_label);
    let (icon_x, text_x) = (16.0, 29.0);
    let labels = [full.as_str(), "Update"];
    let label = labels.iter().filter(|_| labeled).find_map(|text| {
        let width = update_pill_width(ui, text);
        (width <= span.width()).then_some((*text, width))
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

/// Compact primary icon-only Solve button for the kit's 32 px narrow-bar offer slot.
fn solve_compact_button(app: &mut ModelApp, ui: &mut egui::Ui, solve: Rect) {
    let missing = app.unavailable("solve.run");
    let can = missing.is_none() && app.busy().is_none();
    let head = if missing.is_some() { "Solve: not ready yet (or try Home › Samples)" } else { "Solve every step and show the results" };
    let tip = fieldcraft_engine::requirement::Unavailable::tip(head, missing.as_ref());
    let k = solidcraft_theme::Tokens::get();
    let resp = ui
        .add_enabled_ui(can, |ui| {
            let resp = ui.interact(solve, ui.scope_id().with("fc_solve_compact"), Sense::click());
            let fill = if resp.hovered() && can { k.app_bar_hover_button } else { k.form_primary };
            ui.painter().rect_filled(solve, k.control_radius, fill);
            let icon_rect = Rect::from_center_size(solve.center(), vec2(k.space_14, k.space_14));
            kit::icons::paint_ink(ui.painter(), icon_rect, SOLVE_ICON, k.on_accent);
            resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, can, SOLVE));
            resp
        })
        .inner;
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
    let menu = [app_bar::TabMenuRow::new("close", "Close", &close_key)];
    let mut strip = app_bar::TabStrip {
        tabs: vec![app_bar::DocTab { name: &name, list_name: &name, dirty, place: &place, security: None, views: &[], collapsed: false }],
        selected: 0,
        shown_selected: (!app.start_shown()).then_some(0),
        pop_out: false,
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
        views_as_tabs: false,
        switcher_always: false,
        clear_view_tabs: false,
        view_text: app_bar::ViewText {
            close_view: "",
            close_design: "",
            new_view: "",
            more: "",
            switch_view: "",
            collapse_views: "",
            expand_views: "",
        },
        view_menu: &[],
        view_right_row: "",
        view_rename_row: "",
        view_rename: None,
        plus_menu: &[],
        layout_menu: &[],
        view_layout: None,
        layout_current: "",
        layout_tip: "",
        layout_fill: "",
        layout_fill_heading: "",
        switcher_menu: &[],
        tear_off_below: 0.0,
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
            app_bar::TabAction::Select(_) => app.close_start(),
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
    view: "document_view:",
    view_close: "document_view_close:",
    view_new: "document_view_new",
    view_more: "document_view_more",
    switcher: "document_view_switcher",
    view_collapse: "document_view_collapse",
    file_drop: "document_file_drop",
};

/// The shortcut text with this platform's command key.
pub fn shortcut(key: &str) -> String {
    format!("{}+{key}", if cfg!(target_os = "macos") { "Cmd" } else { "Ctrl" })
}

/// Help: how to navigate, keys, updates, crash reports and the version.
fn help_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    // The samples are on the start page (Home) and in File › Samples…, not here: listed here
    // too, they pushed Report a Bug, Crash Reports, Settings and the version below a laptop
    // window's edge.
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
        check = row(ui, "help.check_updates", "Check for Updates", "Ask GitHub Releases whether a newer FieldCraft is published");
        reinstall =
            row(ui, "help.reinstall", "Reinstall FieldCraft…", "Download this version again, verify it and swap it in, to repair a damaged install");
    }
    let crashes = row(ui, "help.crash_reports", "Crash Reports…", "The crash reports saved on this computer: report, copy or delete them");
    let settings = row(ui, "help.settings", "Settings…", "Units, the window, privacy and updates");
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

    /// The bar's Settings opens Settings as the kit's Preferences window (its page list and
    /// footer published under `prefs:`), and no other window.
    #[test]
    fn the_bar_settings_opens_the_kits_preferences_window() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        let run = |app: &mut ModelApp, events: Vec<egui::Event>| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1400.0, 900.0))), events, ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| {
                app.logic(ui.ctx());
                app.ui(ui);
            });
            out.textures_delta.clear();
        };
        for _ in 0..2 {
            run(&mut app, Vec::new());
        }
        assert!(solidcraft_ui_registry::rect_of("prefs:window").is_none(), "closed at launch");
        let at = app.layout.get("app_bar:settings").copied().expect("the bar draws Settings").center();
        let press = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: egui::Modifiers::NONE };
        run(&mut app, vec![egui::Event::PointerMoved(at), press(true)]);
        run(&mut app, vec![press(false)]);
        for _ in 0..2 {
            run(&mut app, Vec::new());
        }
        for name in ["prefs:window", "prefs:page:General", "prefs:page:Privacy", "prefs:page:Updates", "prefs:ok", "prefs:cancel", "prefs:apply"] {
            assert!(solidcraft_ui_registry::rect_of(name).is_some() || solidcraft_ui_registry::handle_at(name).is_some(), "{name} is drawn");
        }
        assert_eq!(app.prompt, None, "no other window opens");
    }

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

    /// Left to right as in SolidCraft: Home, File, its actions, the theme, the document tab,
    /// Settings, Help, Report a bug, and Solve in the reserved offer span; the left group starts
    /// after the traffic lights, and Solve and the right cluster stay put.
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
                "app_bar:settings",
                "app_bar:help",
                "app_bar:report",
                "app_bar:solve",
            ];
            for pair in order.windows(2) {
                let (a, b) = (at(&m, pair[0]), at(&m, pair[1]));
                assert!(a.right() <= b.left() + 0.5, "{} {a:?} before {} {b:?}", pair[0], pair[1]);
            }
            assert!(at(&m, "app_bar:solve").right() < 1400.0 && at(&m, "app_bar:solve").left() > 700.0);
        }
    }

    /// An update on offer is SolidCraft's Update pill at the right end of the offer span after
    /// Settings/Report and Solve; with none on offer there is no pill.
    #[test]
    fn an_update_on_offer_is_a_pill_in_the_offer_span() {
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
        let (settings, solve, update) = (at("app_bar:settings"), at("app_bar:solve"), at("app_bar:update"));
        assert!(settings.right() < solve.left() && solve.right() < update.left(), "{settings:?} {solve:?} {update:?}");
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

    /// One frame of the bar in a 1400 × 600 window with `events`, after the control channel's
    /// synthetic input as the desktop app feeds it; the layout and the kit's registry hold only
    /// this frame's regions, as [`ModelApp::ui`] leaves them.
    fn frame(app: &mut ModelApp, ctx: &egui::Context, events: Vec<egui::Event>) -> egui::FullOutput {
        let mut raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1400.0, 600.0))), events, ..Default::default() };
        app.raw_input_hook(ctx, &mut raw);
        app.layout.clear();
        solidcraft_ui_registry::clear_handles();
        let mut out = ctx.run_ui(raw, |ui| app_bar(app, ui));
        out.textures_delta.clear();
        out
    }

    /// A click of the primary button at `at`, as the window system sends it.
    fn click(app: &mut ModelApp, ctx: &egui::Context, at: egui::Pos2) -> egui::FullOutput {
        frame(app, ctx, vec![egui::Event::PointerMoved(at)]);
        let button = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
        frame(app, ctx, vec![egui::Event::PointerMoved(at), button(true), button(false)])
    }

    fn region(app: &ModelApp, k: &str) -> Rect {
        app.layout.get(k).copied().unwrap_or_else(|| panic!("{k} not drawn"))
    }

    /// Preferences kept in memory, as a host's store keeps them on disk.
    fn memory_store(saved: &std::rc::Rc<std::cell::RefCell<Option<Vec<u8>>>>) -> crate::files::PrefsStore {
        let (read, write) = (saved.clone(), saved.clone());
        crate::files::PrefsStore {
            load: Box::new(move || read.borrow().clone()),
            save: Box::new(move |bytes| {
                *write.borrow_mut() = Some(bytes.to_vec());
                Ok(())
            }),
        }
    }

    /// Donate, then Discord, sit after the document tab and before Report a bug and Solve, each a
    /// capsule in its service's own colour (the kit's `app_bar::brand`), and `ui.inspect`'s layout
    /// names them.
    #[test]
    fn the_community_pills_sit_before_report_a_bug_in_their_own_colours() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        frame(&mut app, &ctx, Vec::new());
        let out = frame(&mut app, &ctx, Vec::new());
        let (tab, donate, discord) = (region(&app, "document_tab"), region(&app, "app_bar:donate"), region(&app, "app_bar:discord"));
        let (solve, report) = (region(&app, "app_bar:solve"), region(&app, "app_bar:report"));
        for (a, b) in [(tab, donate), (donate, discord), (discord, report), (report, solve)] {
            assert!(a.right() <= b.left(), "{a:?} before {b:?}: tab {tab:?} donate {donate:?} discord {discord:?} report {report:?} solve {solve:?}");
        }
        for (pill, fill) in [(donate, app_bar::brand::KOFI), (discord, app_bar::brand::DISCORD)] {
            let painted = out.shapes.iter().any(|s| matches!(&s.shape, egui::Shape::Rect(r) if r.fill == fill && r.rect.contains_rect(pill)));
            assert!(painted, "{pill:?} is not filled with {fill:?}");
        }
        let (tx, _rx) = std::sync::mpsc::channel();
        let crate::control::Outcome::Done(reply) = crate::control::handle(&mut app, &ctx, &crate::ControlRequest::new("ui.inspect", json!({}), tx))
        else {
            panic!("ui.inspect replies at once")
        };
        for k in ["app_bar:donate", "app_bar:discord"] {
            assert!(reply["result"]["layout"][k].is_array(), "{k} not in ui.inspect's layout");
        }
    }

    /// The × on a pill hides that pill alone, and the window remembers it in its preferences: the
    /// next window does not show it either.
    #[test]
    fn hiding_a_pill_hides_only_it_and_is_remembered() {
        let saved = std::rc::Rc::new(std::cell::RefCell::new(None));
        let services = || crate::files::Services { prefs: Some(memory_store(&saved)), ..Default::default() };
        let mut app = ModelApp::new(Session::default()).with_services(services());
        let ctx = egui::Context::default();
        frame(&mut app, &ctx, Vec::new());
        frame(&mut app, &ctx, Vec::new());
        let hide = region(&app, "app_bar:donate:hide").center();
        click(&mut app, &ctx, hide);
        frame(&mut app, &ctx, Vec::new());
        assert!(!app.layout.contains_key("app_bar:donate") && app.layout.contains_key("app_bar:discord"));
        let prefs = serde_json::to_value(app.prefs()).unwrap_or_default();
        assert_eq!((prefs["show_donate_button"].as_bool(), prefs["show_discord_button"].as_bool()), (Some(false), Some(true)), "{prefs}");
        let mut next = ModelApp::new(Session::default()).with_services(services());
        let ctx = egui::Context::default();
        frame(&mut next, &ctx, Vec::new());
        frame(&mut next, &ctx, Vec::new());
        assert!(!next.layout.contains_key("app_bar:donate") && next.layout.contains_key("app_bar:discord"));
    }

    /// A click on a pill opens its page through the host's URL opener; the same click sent over
    /// the control channel opens nothing.
    #[test]
    fn a_click_opens_the_page_and_a_control_channel_click_does_not() {
        let opened = std::rc::Rc::new(std::cell::RefCell::new(Vec::<String>::new()));
        let log = opened.clone();
        let services = crate::files::Services {
            open_url: Some(Box::new(move |url| {
                log.borrow_mut().push(url.to_string());
                Ok(())
            })),
            ..Default::default()
        };
        let mut app = ModelApp::new(Session::default()).with_services(services);
        let ctx = egui::Context::default();
        frame(&mut app, &ctx, Vec::new());
        frame(&mut app, &ctx, Vec::new());
        for k in ["app_bar:donate", "app_bar:discord"] {
            let at = region(&app, k).center();
            click(&mut app, &ctx, at);
        }
        assert_eq!(*opened.borrow(), ["https://ko-fi.com/jyviko/tiers", "https://discord.gg/FtU4NQNVDH"]);
        for k in ["app_bar:donate", "app_bar:discord"] {
            let at = region(&app, k).center();
            let (tx, _rx) = std::sync::mpsc::channel();
            crate::control::handle(&mut app, &ctx, &crate::ControlRequest::new("ui.click", json!({"x": at.x, "y": at.y}), tx));
            while !app.synthetic.is_empty() {
                let out = frame(&mut app, &ctx, Vec::new());
                assert!(!out.platform_output.commands.iter().any(|c| matches!(c, egui::OutputCommand::OpenUrl(_))), "{k}");
            }
            frame(&mut app, &ctx, Vec::new());
        }
        assert_eq!(opened.borrow().len(), 2, "{:?}", opened.borrow());
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

    /// In a narrow bar that hides items into `More`, clicking `More` opens the overflow menu and
    /// choosing an app-owned control (such as the theme button) dispatches its action.
    #[test]
    fn narrow_bar_overflow_menu_runs_hidden_controls() {
        let mut app = ModelApp::new(Session::default());
        let start_dark = app.ui.dark;
        let ctx = egui::Context::default();
        let draw = |app: &mut ModelApp, events: Vec<egui::Event>| {
            let mut raw =
                egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(420.0, 600.0))), events, ..Default::default() };
            app.raw_input_hook(&ctx, &mut raw);
            app.layout.clear();
            solidcraft_ui_registry::clear_handles();
            let mut out = ctx.run_ui(raw, |ui| app_bar(app, ui));
            out.textures_delta.clear();
        };
        let click_at = |app: &mut ModelApp, at: egui::Pos2| {
            draw(app, vec![egui::Event::PointerMoved(at)]);
            let button =
                |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
            draw(app, vec![egui::Event::PointerMoved(at), button(true), button(false)]);
        };
        draw(&mut app, Vec::new());
        draw(&mut app, Vec::new());
        let more = solidcraft_ui_registry::handle_at("appbar:overflow").unwrap_or_else(|| panic!("no appbar:overflow handle"));
        click_at(&mut app, more);
        let theme_row =
            solidcraft_ui_registry::handle_at("menu:appbar:overflow:quick:0").unwrap_or_else(|| panic!("no menu:appbar:overflow:quick:0 handle"));
        click_at(&mut app, theme_row);
        assert_ne!(app.ui.dark, start_dark, "clicking overflowed theme row should toggle dark mode");
    }

    /// Across narrow widths (down to where the kit steps the offer slot to 32 px or moves items
    /// into More), Solve never silently disappears: it is always either drawn on the bar or
    /// present in the More overflow menu.
    #[test]
    fn narrow_bar_always_exposes_solve_on_bar_or_in_overflow() {
        for has_update in [false, true] {
            for width in [360.0_f32, 420.0, 480.0, 540.0, 600.0, 680.0, 760.0] {
                let mut app = ModelApp::new(Session::default());
                if has_update {
                    app.reload = Some(Box::new(|| {}));
                }
                let ctx = egui::Context::default();
                let draw = |app: &mut ModelApp, events: Vec<egui::Event>| {
                    let mut raw =
                        egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(width, 600.0))), events, ..Default::default() };
                    app.raw_input_hook(&ctx, &mut raw);
                    app.layout.clear();
                    solidcraft_ui_registry::clear_handles();
                    let mut out = ctx.run_ui(raw, |ui| app_bar(app, ui));
                    out.textures_delta.clear();
                };
                draw(&mut app, Vec::new());
                draw(&mut app, Vec::new());
                if let Some(solve) = app.layout.get("app_bar:solve").copied() {
                    assert!(solve.width() >= COMPACT_SOLVE_MIN_W, "width={width} has_update={has_update}: drawn Solve too narrow: {solve:?}");
                } else {
                    let more = solidcraft_ui_registry::handle_at("appbar:overflow")
                        .unwrap_or_else(|| panic!("width={width} has_update={has_update}: neither Solve nor More overflow button is drawn"));
                    draw(&mut app, vec![egui::Event::PointerMoved(more)]);
                    let button = |pressed| egui::Event::PointerButton {
                        pos: more,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: Default::default(),
                    };
                    draw(&mut app, vec![egui::Event::PointerMoved(more), button(true), button(false)]);
                    assert!(
                        solidcraft_ui_registry::handle_at("menu:appbar:overflow:offered").is_some(),
                        "width={width} has_update={has_update}: Solve missing from both bar and More overflow menu"
                    );
                }
            }
        }
    }
}
