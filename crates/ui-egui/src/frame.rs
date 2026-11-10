//! The one place windows and panel headings are made, so every prompt and panel shares the same
//! width limits and heading bands.
//!
//! Adapted from SolidCraft's `crates/ui-egui/src/frame.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.

use egui::{Align2, Ui, WidgetText, vec2};

/// Floating windows: prompts and panels, and those holding a table.
pub const WINDOW_MAX_W: f32 = 460.0;
pub const WINDOW_WIDE_MAX_W: f32 = 760.0;

/// How wide a window may get.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Width {
    Normal,
    Wide,
}

/// A floating window, centred on the screen, shrink-wrapped to its content up to its width limit.
pub fn window<'a>(ctx: &egui::Context, title: impl Into<WidgetText>, width: Width) -> egui::Window<'a> {
    let max = match width {
        Width::Normal => WINDOW_MAX_W,
        Width::Wide => WINDOW_WIDE_MAX_W,
    };
    let c = ctx.content_rect().center();
    egui::Window::new(title).collapsible(false).resizable(false).pivot(Align2::CENTER_CENTER).default_pos(c - vec2(0.0, 60.0)).max_width(max)
}

/// Heading levels inside a panel: its title, and the sections below it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    /// The panel's title row (`Tokens::heading_bg`).
    Title,
    /// A section heading (`Tokens::subheading_bg`).
    Section,
}

/// Horizontal inner margin of a heading band, on each side.
pub const HEADING_MARGIN: f32 = 4.0;

/// A heading row drawn on a band as wide as the panel; `add` lays out the row left to right.
pub fn heading<R>(ui: &mut Ui, t: &crate::theme::Tokens, level: Level, add: impl FnOnce(&mut Ui) -> R) -> R {
    let fill = match level {
        Level::Title => t.heading_bg,
        Level::Section => t.subheading_bg,
    };
    let inner = (ui.available_width() - 2.0 * HEADING_MARGIN).max(0.0);
    egui::Frame::new()
        .fill(fill)
        .corner_radius(3)
        .inner_margin(egui::Margin::symmetric(HEADING_MARGIN as i8, 2))
        .show(ui, |ui| {
            ui.set_width(inner);
            ui.horizontal(add).inner
        })
        .inner
}

/// A menu below `anchor`, open while [`crate::UiState::menu`] is `key` (so `ui.set {menu}` and
/// offscreen captures open it too). A click on one of its rows, a click outside or Escape closes
/// it. Returns where it was drawn.
pub fn menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    popup(app, ui, key, anchor, egui::PopupCloseBehavior::CloseOnClick, add)
}

/// As [`menu`], but clicks inside keep it open (settings with checkboxes).
pub fn settings_menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    popup(app, ui, key, anchor, egui::PopupCloseBehavior::CloseOnClickOutside, add)
}

fn popup<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    close: egui::PopupCloseBehavior,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    let mut open = app.ui.menu.as_deref() == Some(key);
    if !open {
        return None;
    }
    let id = egui::Id::unique(("fc_menu", key));
    let shown =
        egui::Popup::new(id, ui.ctx().clone(), anchor, ui.layer_id()).kind(egui::PopupKind::Menu).open_bool(&mut open).close_behavior(close).show(
            |ui| {
                ui.set_min_width(200.0);
                add(app, ui)
            },
        );
    if (!open || ui.input(|i| i.key_pressed(egui::Key::Escape))) && app.ui.menu.as_deref() == Some(key) {
        app.ui.menu = None;
    }
    shown.map(|r| r.response.rect)
}

/// A menu row's width and height.
pub const MENU_W: f32 = 240.0;
const MENU_ROW_H: f32 = 24.0;

/// A menu row: its label, and its shortcut (or a note) at the right; dimmed and inert when
/// disabled. Painted like the toolbar's drop-down rows.
pub fn menu_row(ui: &mut Ui, t: &crate::theme::Tokens, label: &str, right: Option<&str>, enabled: bool) -> egui::Response {
    let (rect, resp) = ui.allocate_exact_size(vec2(MENU_W, MENU_ROW_H), if enabled { egui::Sense::click() } else { egui::Sense::hover() });
    if resp.hovered() && enabled {
        ui.painter().rect_filled(rect, 3.0, t.hover);
    }
    let font = egui::FontId::proportional(13.0);
    let clip = rect.shrink2(vec2(8.0, 0.0));
    ui.painter().with_clip_rect(clip).text(
        egui::pos2(rect.left() + 10.0, rect.center().y),
        Align2::LEFT_CENTER,
        label,
        font.clone(),
        if enabled { t.text } else { t.text_dim },
    );
    if let Some(r) = right {
        ui.painter().text(egui::pos2(rect.right() - 10.0, rect.center().y), Align2::RIGHT_CENTER, r, egui::FontId::proportional(12.0), t.text_dim);
    }
    resp.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, enabled, label));
    resp
}

/// A menu section's small heading.
pub fn menu_heading(ui: &mut Ui, t: &crate::theme::Tokens, text: &str) {
    ui.add_space(2.0);
    ui.label(egui::RichText::new(text).font(crate::theme::label_font(ui.ctx(), 10.5)).color(t.text_dim));
}

#[cfg(test)]
mod tests {
    /// Windows are made in this file only, so none can skip the width limits.
    #[test]
    fn windows_come_from_the_frame_helper() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut bad = Vec::new();
        for e in std::fs::read_dir(&dir).unwrap() {
            let p = e.unwrap().path();
            if p.extension().is_some_and(|x| x == "rs") && p.file_name().is_some_and(|n| n != "frame.rs") {
                let text = std::fs::read_to_string(&p).unwrap();
                if text.contains(concat!("Window", "::new(")) {
                    bad.push(p.display().to_string());
                }
            }
        }
        assert!(bad.is_empty(), "use crate::frame::window instead of egui::Window::new in {bad:?}");
    }
}
