//! The ribbon: the workspace switcher, the workspace's tabs, and the current tab's panels, all
//! from the capability catalogue (`fieldcraft_engine::catalogue`).
//!
//! Each panel shows its ready entries as icon buttons (up to [`MAX_BUTTONS`]; a panel with none
//! ready shows its first entry, greyed with its milestone) over its label with a caret. The label
//! opens the panel's drop-down, which lists every entry of the panel: ready ones with their
//! related commands (edit, delete…) under them, coming-soon ones greyed with their milestone.
//! When the window is too narrow, panels collapse from the right into one button that opens the
//! same drop-down, then into the » menu at the end. Nothing is dropped.
//!
//! A ready entry starts its command ([`ModelApp::start`]); a coming-soon entry does nothing and
//! its tooltip says when it comes. Structure follows SolidCraft's `crates/ui-egui/src/toolbar.rs`
//! (MIT OR Apache-2.0), see ATTRIBUTION.md; the tab row is the kit's (`solidcraft_ui_kit::ribbon`).

use egui::{Align2, Color32, FontId, Rect, Sense, Stroke, Vec2, pos2, vec2};
use fieldcraft_engine::catalogue::{self, Capability, Tab, Workspace};
use fieldcraft_engine::requirement::Unavailable;
use serde::Serialize;
use solidcraft_theme::TabAccent;
use solidcraft_ui_kit::ribbon;

use crate::ModelApp;
use crate::icons::{self, Look};
use crate::theme::{Tokens, label_font};

/// The ribbon's height: the tab row, then the panels.
pub const HEIGHT: f32 = TABS_HEIGHT + 66.0;
/// The tab row's height (the kit's).
const TABS_HEIGHT: f32 = ribbon::TAB_ROW_HEIGHT;
/// Most buttons a panel shows; the rest are in its drop-down.
pub const MAX_BUTTONS: usize = 4;
/// A button's side and the room it takes.
const BUTTON: f32 = 38.0;
const STRIDE: f32 = 40.0;
/// Where the tabs and the panels start, right of the workspace switcher.
const PANELS_X: f32 = 112.0;
/// A collapsed panel's width, and the gap after each panel (its separator).
pub const COLLAPSED_WIDTH: f32 = 72.0;
const PANEL_GAP: f32 = 2.0;
/// The » menu's width.
const MORE_WIDTH: f32 = 40.0;
/// A drop-down row's height.
const ROW_H: f32 = 24.0;

/// How an entry is shown on a toolbar that may be too narrow for every panel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Fit {
    /// A button of its panel.
    Full,
    /// Its panel is one button that opens the drop-down.
    Collapsed,
    /// Its panel is in the » menu.
    More,
    /// In its panel's drop-down (a row while it is open).
    Menu,
}

/// A panel's width with `buttons` buttons.
pub fn panel_width(buttons: usize) -> f32 {
    (buttons.max(1) as f32 * STRIDE).max(64.0) + 8.0
}

/// Fit panels of these full widths into `room`: all in full when they fit; otherwise panels
/// collapse into their drop-down from the right, then (all collapsed and still too wide) the last
/// ones go into the » menu. Nothing is dropped.
pub fn fit(widths: &[f32], room: f32) -> Vec<Fit> {
    let mut fits = vec![Fit::Full; widths.len()];
    let total = |fits: &[Fit]| -> f32 {
        let more = if fits.contains(&Fit::More) { MORE_WIDTH } else { 0.0 };
        more + widths
            .iter()
            .zip(fits)
            .map(|(w, f)| match f {
                Fit::Full | Fit::Menu => w + PANEL_GAP,
                Fit::Collapsed => COLLAPSED_WIDTH.min(*w) + PANEL_GAP,
                Fit::More => 0.0,
            })
            .sum::<f32>()
    };
    for stage in [Fit::Collapsed, Fit::More] {
        for i in (0..widths.len()).rev() {
            if total(&fits) <= room {
                return fits;
            }
            if let Some(f) = fits.get_mut(i) {
                *f = stage;
            }
        }
    }
    fits
}

/// The buttons a panel shows: its ready entries (at most [`MAX_BUTTONS`]), or its first entry
/// when none is ready yet.
pub fn buttons(entries: &'static [Capability]) -> Vec<&'static Capability> {
    let ready: Vec<_> = entries.iter().filter(|c| c.command().is_some()).take(MAX_BUTTONS).collect();
    if ready.is_empty() { entries.first().into_iter().collect() } else { ready }
}

/// The key of a panel's drop-down in [`crate::UiState::menu`]: `<tab>/<panel>`.
pub fn panel_key(tab: Tab, panel: &str) -> String {
    format!("{}/{panel}", tab.id())
}

/// A small drop-down caret.
pub fn caret(p: &egui::Painter, c: egui::Pos2, col: Color32) {
    p.add(egui::Shape::convex_polygon(vec![pos2(c.x - 3.5, c.y - 2.0), pos2(c.x + 3.5, c.y - 2.0), pos2(c.x, c.y + 2.5)], col, Stroke::NONE));
}

/// A milestone badge ("M5") with its top right corner at `corner`.
fn badge(p: &egui::Painter, t: &Tokens, corner: egui::Pos2, m: &str) {
    let text = p.layout_no_wrap(m.to_string(), FontId::proportional(9.0), t.text_dim);
    let r = Rect::from_min_size(pos2(corner.x - text.size().x - 4.0, corner.y), text.size() + vec2(4.0, 0.0));
    p.rect_filled(r, 3.0, t.panel_header);
    p.rect_stroke(r, 3.0, Stroke::new(1.0, t.border), egui::StrokeKind::Inside);
    p.galley(r.min + vec2(2.0, 0.0), text, t.text_dim);
}

pub fn toolbar(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let mut start: Option<&'static str> = None;
    egui::Panel::top("fc_toolbar").exact_size(HEIGHT).show_separator_line(false).frame(egui::Frame::NONE.fill(t.toolbar)).show(ui, |ui| {
        let r = ui.max_rect();
        app.place("ribbon", r);
        app.toolbar_rects.clear();
        app.panel_rects.clear();
        let painter = ui.painter().clone();
        workspace_switcher(app, ui, &t, r);
        // Tabs: the kit's tab row, in the main accent.
        let tabs = app.ui.workspace.tabs();
        let accent = TabAccent::current(ui.ctx());
        let row = ribbon::Tabs {
            tabs: tabs.iter().map(|tab| ribbon::Tab { id: tab.id(), label: tab.label(), accent, contextual: false, aliases: &[] }).collect(),
            active: app.ui.tab.id(),
            add: None,
            menu_open_from: &|_| false,
            search_width: 0.0,
            buttons: Vec::new(),
            floating: false,
            underline: None,
            handles: ribbon::TabHandles { tab: "tab:", tab_rect: "tab_rect:", add: "tab_add", action: "tab_action:" },
        };
        let row_rect = Rect::from_min_max(pos2(r.left() + PANELS_X, r.top()), pos2(r.right(), r.top() + TABS_HEIGHT));
        for action in ribbon::tabs(ui, &row, row_rect).actions {
            if let ribbon::TabAction::Select { index, .. } = action
                && let Some(&tab) = tabs.get(index)
            {
                app.ui.tab = tab;
                app.ui.menu = None;
            }
        }
        // Panels.
        let tab = app.ui.tab;
        let top = r.top() + TABS_HEIGHT;
        let panels: Vec<&'static [Capability]> = catalogue::panels(tab).collect();
        let shown: Vec<Vec<&'static Capability>> = panels.iter().map(|p| buttons(p)).collect();
        // A panel is as wide as its buttons, and wide enough for its label's longest word (plus
        // the caret beside it).
        let label_word = |entries: &[Capability]| {
            let label = entries.first().map_or(String::new(), |c| c.panel.to_uppercase());
            label.split_whitespace().map(|w| painter.layout_no_wrap(w.to_string(), label_font(ui.ctx(), 10.5), t.text).size().x).fold(0.0, f32::max)
        };
        let widths: Vec<f32> = shown.iter().zip(&panels).map(|(b, p)| panel_width(b.len()).max(label_word(p) + 26.0)).collect();
        let fits = fit(&widths, r.right() - 8.0 - (r.left() + PANELS_X));
        let mut px = r.left() + PANELS_X;
        let mut overflow: Vec<&'static [Capability]> = Vec::new();
        for ((entries, shown), (fit, full)) in panels.iter().zip(&shown).zip(fits.iter().zip(&widths)) {
            let Some(panel) = entries.first().map(|c| c.panel) else { continue };
            if *fit == Fit::More {
                for c in entries.iter() {
                    app.toolbar_rects.push((c.id, Rect::NOTHING, Fit::More));
                }
                overflow.push(entries);
                continue;
            }
            let key = panel_key(tab, panel);
            let collapsed = *fit == Fit::Collapsed;
            let width = if collapsed { COLLAPSED_WIDTH.min(*full) } else { *full };
            if collapsed {
                // One button: the panel's first button's icon over a caret; it opens the drop-down.
                let br = Rect::from_center_size(pos2(px + width / 2.0, top + 21.0), Vec2::splat(BUTTON));
                let resp = ui.interact(br, ui.scope_id().with(("fc_collapsed", &key)), Sense::click());
                if resp.hovered() {
                    painter.rect_filled(br, 4.0, t.hover);
                }
                if let Some(c) = shown.first() {
                    icons::paint(
                        app,
                        ui,
                        c.id,
                        Rect::from_center_size(br.center() - vec2(0.0, 3.0), Vec2::splat(26.0)),
                        Look::of(c.command().is_some() && app.unavailable(c.id).is_none()),
                    );
                }
                caret(&painter, pos2(br.center().x, br.bottom() - 3.0), t.text_dim);
                for c in entries.iter() {
                    app.toolbar_rects.push((c.id, br, Fit::Collapsed));
                }
                if resp.on_hover_text(format!("{panel} ({})", entries.len())).clicked() {
                    app.toggle_menu(&key);
                }
            } else {
                let inset = ((width - 8.0 - (shown.len() as f32 * STRIDE - 2.0)) / 2.0).max(0.0).floor();
                for (i, c) in shown.iter().enumerate() {
                    let br = Rect::from_min_size(pos2(px + 4.0 + inset + i as f32 * STRIDE, top + 2.0), Vec2::splat(BUTTON));
                    if let Some(id) = button(app, ui, &t, c, br) {
                        start = Some(id);
                    }
                }
                for c in entries.iter().filter(|c| !shown.iter().any(|s| s.id == c.id)) {
                    app.toolbar_rects.push((c.id, Rect::NOTHING, Fit::Menu));
                }
            }
            // The panel's label with its caret opens the drop-down.
            // A label too wide wraps between words onto a second line; a single word wider than
            // the panel stays whole.
            let label = panel.to_uppercase();
            let widest_word = label
                .split_whitespace()
                .map(|w| painter.layout_no_wrap(w.to_string(), label_font(ui.ctx(), 10.5), t.text).size().x)
                .fold(0.0, f32::max);
            let mut job = egui::text::LayoutJob::simple(label, label_font(ui.ctx(), 10.5), t.text, (width - 14.0).max(widest_word + 1.0));
            job.wrap.overflow_character = None;
            job.wrap.max_rows = 2;
            job.halign = egui::Align::Center;
            let g = painter.layout_job(job);
            let gw = g.size().x;
            // One line centred under the buttons; two lines start higher to stay inside the ribbon.
            let lh = (g.size().y + 2.0).max(18.0);
            let lr = if lh > 20.0 {
                Rect::from_min_size(pos2(px, top + 40.0), vec2(width, lh))
            } else {
                Rect::from_center_size(pos2(px + width / 2.0, top + 51.0), vec2(width, lh))
            };
            let lresp = ui.interact(lr, ui.scope_id().with(("fc_panel", &key)), Sense::click());
            let open = app.ui.menu.as_deref() == Some(key.as_str());
            if lresp.hovered() || open {
                painter.rect_filled(lr, 3.0, t.hover);
            }
            let gx = lr.center().x - 4.5;
            painter.galley(pos2(gx, lr.center().y - g.size().y / 2.0), g, t.text);
            caret(&painter, pos2(gx + gw / 2.0 + 6.0, lr.center().y), t.text);
            app.panel_rects.push((key.clone(), lr));
            if lresp.on_hover_text(format!("All {panel} commands")).clicked() {
                app.toggle_menu(&key);
            }
            if let Some(rect) = crate::frame::menu(app, ui, &key, lr, |app, ui| dropdown(app, ui, entries)) {
                app.place("toolbar_dropdown", rect);
            }
            px += width;
            painter.line_segment([pos2(px, top + 4.0), pos2(px, top + 58.0)], Stroke::new(1.0, t.border));
            px += PANEL_GAP;
        }
        if !overflow.is_empty() {
            let mr = Rect::from_min_size(pos2(px.min(r.right() - MORE_WIDTH), top + 4.0), vec2(MORE_WIDTH - 2.0, 52.0));
            let resp = ui.interact(mr, ui.scope_id().with("fc_more"), Sense::click());
            if resp.hovered() || app.ui.menu.as_deref() == Some("more") {
                painter.rect_filled(mr, 4.0, t.hover);
            }
            painter.text(mr.center(), Align2::CENTER_CENTER, "»", FontId::proportional(18.0), t.text);
            if resp.on_hover_text("More toolbar panels").clicked() {
                app.toggle_menu("more");
            }
            crate::frame::menu(app, ui, "more", mr, |app, ui| {
                for entries in &overflow {
                    let panel = entries.first().map_or("", |c| c.panel);
                    ui.menu_button(panel, |ui| dropdown(app, ui, entries));
                }
            });
        }
        painter.line_segment([pos2(r.left(), r.bottom() - 0.5), pos2(r.right(), r.bottom() - 0.5)], Stroke::new(1.0, t.border));
    });
    if let Some(id) = start {
        app.ui.menu = None;
        // Errors already show in the status bar.
        let _ = app.start(id);
    }
}

/// The workspace switcher: a box spanning the ribbon's height with the workspace's icon and name;
/// it opens the list of workspaces.
fn workspace_switcher(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, r: Rect) {
    let ws = Rect::from_min_max(pos2(r.left() + 6.0, r.top() + 4.0), pos2(r.left() + 104.0, r.bottom() - 5.0));
    let resp = ui.interact(ws, ui.scope_id().with("fc_workspace"), Sense::click());
    let open = app.ui.menu.as_deref() == Some("workspace");
    ui.painter().rect(ws, 4.0, if resp.hovered() || open { t.hover } else { t.field }, Stroke::new(1.0, t.border), egui::StrokeKind::Inside);
    let w = app.ui.workspace;
    icons::paint(app, ui, w.icon(), Rect::from_center_size(pos2(ws.center().x, ws.top() + 24.0), Vec2::splat(26.0)), Look::Normal);
    let p = ui.painter();
    let label = w.label().to_uppercase();
    let font = label_font(ui.ctx(), 12.0);
    let label_w = p.layout_no_wrap(label.clone(), font.clone(), t.text).size().x;
    p.text(pos2(ws.center().x - 5.0, ws.top() + 52.0), Align2::CENTER_CENTER, label, font, t.text);
    caret(p, pos2(ws.center().x - 5.0 + label_w / 2.0 + 7.0, ws.top() + 52.0), t.text);
    p.text(pos2(ws.center().x, ws.top() + 69.0), Align2::CENTER_CENTER, "workspace", FontId::proportional(10.0), t.text_dim);
    app.place("workspace_switcher", ws);
    if resp.on_hover_text("Change workspace").clicked() {
        app.toggle_menu("workspace");
    }
    crate::frame::menu(app, ui, "workspace", ws, |app, ui| {
        let t = Tokens::get(app.ui.dark);
        for w in Workspace::ALL {
            let (rect, resp) = ui.allocate_exact_size(vec2(360.0, 40.0), Sense::click());
            if resp.hovered() || w == app.ui.workspace {
                ui.painter().rect_filled(rect, 3.0, if w == app.ui.workspace { t.accent_soft } else { t.hover });
            }
            icons::paint(app, ui, w.icon(), Rect::from_center_size(pos2(rect.left() + 20.0, rect.center().y), Vec2::splat(24.0)), Look::Normal);
            let p = ui.painter();
            p.text(pos2(rect.left() + 40.0, rect.top() + 12.0), Align2::LEFT_CENTER, w.label().to_uppercase(), label_font(ui.ctx(), 12.5), t.text);
            p.text(pos2(rect.left() + 40.0, rect.top() + 28.0), Align2::LEFT_CENTER, w.about(), FontId::proportional(11.0), t.text_dim);
            resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, w.label()));
            if resp.clicked() {
                app.ui.select_workspace(w);
                app.ui.menu = None;
            }
        }
    });
}

/// One entry's button; returns the id to start when it was clicked.
fn button(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, c: &'static Capability, br: Rect) -> Option<&'static str> {
    let missing = app.unavailable(c.id);
    let enabled = c.command().is_some() && missing.is_none();
    let resp = ui.interact(br, ui.scope_id().with(("fc_button", c.id)), if enabled { Sense::click() } else { Sense::hover() });
    app.toolbar_rects.push((c.id, br, Fit::Full));
    let running = app.dialog.as_ref().is_some_and(|d| c.command() == Some(d.command));
    if running {
        ui.painter().rect_filled(br, 4.0, t.border.gamma_multiply(0.7));
    } else if resp.hovered() && enabled {
        ui.painter().rect_filled(br, 4.0, t.hover);
    }
    icons::paint(app, ui, c.id, br.shrink(6.0), Look::of(enabled));
    if let Some(m) = c.milestone() {
        badge(ui.painter(), t, br.right_top() + vec2(4.0, -2.0), m);
    }
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, enabled, c.label));
    let resp = resp.on_hover_text(Unavailable::tip(&c.tooltip(), missing.as_ref()));
    resp.clicked().then_some(c.id)
}

/// The rows of a panel's drop-down: every entry with its related commands under it, coming-soon
/// entries greyed with their milestone. A click starts the row's command.
pub fn dropdown(app: &mut ModelApp, ui: &mut egui::Ui, entries: &'static [Capability]) {
    let t = Tokens::get(app.ui.dark);
    let font = FontId::proportional(13.0);
    let width = entries
        .iter()
        .flat_map(|c| std::iter::once(c.label).chain(c.more.iter().map(|id| crate::command_label(id))))
        .map(|l| ui.painter().layout_no_wrap(l.to_string(), font.clone(), t.text).size().x)
        .fold(0.0_f32, f32::max)
        .max(140.0)
        + 96.0;
    let mut start = None;
    for c in entries {
        let missing = app.unavailable(c.id);
        let enabled = c.command().is_some() && missing.is_none();
        let (rect, resp) = ui.allocate_exact_size(vec2(width, ROW_H), if enabled { Sense::click() } else { Sense::hover() });
        app.toolbar_rects.push((c.id, rect, Fit::Menu));
        if resp.hovered() && enabled {
            ui.painter().rect_filled(rect, 3.0, t.hover);
        }
        icons::paint(app, ui, c.id, Rect::from_center_size(pos2(rect.left() + 15.0, rect.center().y), Vec2::splat(18.0)), Look::of(enabled));
        let p = ui.painter();
        p.text(pos2(rect.left() + 32.0, rect.center().y), Align2::LEFT_CENTER, c.label, font.clone(), if enabled { t.text } else { t.text_dim });
        if let Some(m) = c.milestone() {
            p.text(pos2(rect.right() - 8.0, rect.center().y), Align2::RIGHT_CENTER, format!("soon · {m}"), FontId::proportional(11.0), t.text_dim);
        }
        resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, enabled, c.label));
        if resp.on_hover_text(Unavailable::tip(&c.tooltip(), missing.as_ref())).clicked() && enabled {
            start = Some(c.id);
        }
        for id in c.more {
            let missing = app.unavailable(id);
            let enabled = missing.is_none();
            let (rect, resp) = ui.allocate_exact_size(vec2(width, ROW_H - 2.0), if enabled { Sense::click() } else { Sense::hover() });
            if resp.hovered() && enabled {
                ui.painter().rect_filled(rect, 3.0, t.hover);
            }
            let label = crate::command_label(id);
            let p = ui.painter();
            p.text(
                pos2(rect.left() + 32.0, rect.center().y),
                Align2::LEFT_CENTER,
                label,
                FontId::proportional(12.5),
                if enabled { t.text } else { t.text_dim },
            );
            p.line_segment([pos2(rect.left() + 15.0, rect.top()), pos2(rect.left() + 15.0, rect.center().y)], Stroke::new(1.0, t.border));
            p.line_segment([pos2(rect.left() + 15.0, rect.center().y), pos2(rect.left() + 25.0, rect.center().y)], Stroke::new(1.0, t.border));
            resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, enabled, label));
            let tip = Unavailable::tip(id, missing.as_ref());
            if resp.on_hover_text(tip).clicked() && enabled {
                start = Some(id);
            }
        }
    }
    if let Some(id) = start {
        app.ui.menu = None;
        ui.close();
        let _ = app.start(id);
    }
}

/// The toolbar as the control channel reports it: the workspaces and tabs, the current ones, the
/// panels' labels (screen points, last frame), the open menu, and every entry with its status
/// (`enabled`, and `unavailable`: what its command is missing, `Unavailable::to_json`) and where
/// it is: a button (`full`), a collapsed panel, the » menu, or its panel's drop-down
/// (`menu`, with the row's rect while the drop-down is open).
pub fn describe(app: &ModelApp) -> serde_json::Value {
    let rect = |r: &Rect| {
        if r.is_finite() && r.is_positive() { serde_json::json!([r.left(), r.top(), r.width(), r.height()]) } else { serde_json::Value::Null }
    };
    let entries: Vec<serde_json::Value> = catalogue::CATALOGUE
        .iter()
        .map(|c| {
            let shown = app.toolbar_rects.iter().rev().find(|(id, ..)| *id == c.id);
            let missing = app.unavailable(c.id);
            serde_json::json!({
                "id": c.id,
                "label": c.label,
                "tab": c.tab,
                "panel": c.panel,
                "command": c.command(),
                "more": c.more,
                "covers": c.covers,
                "coming_soon": c.milestone(),
                "enabled": c.command().is_some() && missing.is_none(),
                "unavailable": missing.map(|u| u.to_json()),
                "rect": shown.map(|(_, r, _)| rect(r)),
                "fit": shown.map(|(.., f)| f),
            })
        })
        .collect();
    let panels: Vec<serde_json::Value> = app.panel_rects.iter().map(|(k, r)| serde_json::json!({"menu": k, "rect": rect(r)})).collect();
    serde_json::json!({
        "workspace": app.ui.workspace,
        "workspaces": Workspace::ALL,
        "tab": app.ui.tab,
        "tabs": app.ui.workspace.tabs(),
        "menu": app.ui.menu,
        "panels": panels,
        "entries": entries,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panels_collapse_from_the_right_then_go_to_the_more_menu() {
        let widths = [200.0, 300.0, 100.0];
        assert_eq!(fit(&widths, 700.0), vec![Fit::Full; 3]);
        assert_eq!(fit(&widths, 500.0), vec![Fit::Full, Fit::Collapsed, Fit::Collapsed]);
        assert_eq!(fit(&widths, 230.0), vec![Fit::Collapsed; 3]);
        let tight = fit(&widths, 120.0);
        assert_eq!(tight.first(), Some(&Fit::Collapsed));
        assert!(tight.contains(&Fit::More), "nothing is dropped: {tight:?}");
        assert!(fit(&[], 0.0).is_empty());
    }

    #[test]
    fn panels_show_their_ready_entries_or_their_first_one() {
        for tab in Tab::ALL {
            for p in catalogue::panels(tab) {
                let b = buttons(p);
                assert!(!b.is_empty() && b.len() <= MAX_BUTTONS, "{tab:?}");
                if p.iter().any(|c| c.command().is_some()) {
                    assert!(b.iter().all(|c| c.command().is_some()));
                }
            }
        }
        assert_eq!(panel_key(Tab::Loads, "Loads"), "loads/Loads");
    }
}
