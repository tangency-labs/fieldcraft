//! The one place windows, menus and panel headings are made, so every prompt and panel shares the
//! same frames: the UI kit's windows, confirmations and menus (`solidcraft_ui_kit`), and the
//! heading bands of the context panel.

use egui::Ui;
use solidcraft_ui_kit as kit;

/// How wide a window is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Width {
    /// Prompts and small panels.
    Normal,
    /// Windows holding a table or a list.
    Wide,
}

impl Width {
    /// The window's width in the kit's form cards.
    fn cards(self) -> f32 {
        match self {
            Width::Normal => 1.75,
            Width::Wide => 2.9,
        }
    }
}

/// A system window (the kit's): `icon` (an icon id) and `title` in its header with a close button,
/// centred on the screen, `width` wide. Footer actions are set with [`kit::Window::actions`].
pub fn window<'a>(id: egui::Id, icon: &'a str, title: &'a str, width: Width) -> kit::Window<'a> {
    kit::Window::new(id, icon, title).cards(width.cards())
}

/// Where a window or confirmation `id` was drawn last frame.
pub fn drawn(ctx: &egui::Context, id: egui::Id) -> Option<egui::Rect> {
    ctx.memory(|m| m.area_rect(id))
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
/// offscreen captures open it too): the kit's popup and list surface. A click on one of its rows,
/// a click outside or Escape closes it. Returns where it was drawn.
pub fn menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    popup(app, ui, key, anchor, (MENU_W, egui::PopupCloseBehavior::CloseOnClick), add)
}

/// As [`menu`], [`WIDE_MENU_W`] wide, for a menu with notes and cards (Help).
pub fn wide_menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    popup(app, ui, key, anchor, (WIDE_MENU_W, egui::PopupCloseBehavior::CloseOnClick), add)
}

/// As [`menu`], but clicks inside keep it open (settings with checkboxes).
pub fn settings_menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    popup(app, ui, key, anchor, (MENU_W, egui::PopupCloseBehavior::CloseOnClickOutside), add)
}

fn popup<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    (width, close): (f32, egui::PopupCloseBehavior),
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    let mut open = app.ui.menu.as_deref() == Some(key);
    if !open {
        return None;
    }
    let id = egui::Id::unique(("fc_menu", key));
    // egui's popup places it and closes it; the kit draws its surface and rows.
    let shown = egui::Popup::new(id, ui.ctx().clone(), anchor, ui.layer_id())
        .kind(egui::PopupKind::Menu)
        .open_bool(&mut open)
        .close_behavior(close)
        .frame(egui::Frame::NONE)
        .show(|ui| kit::Menu::show(ui, width, |ui| add(app, ui)));
    if (!open || ui.input(|i| i.key_pressed(egui::Key::Escape))) && app.ui.menu.as_deref() == Some(key) {
        app.ui.menu = None;
    }
    shown.map(|r| r.response.rect)
}

/// The menu a right click on `response` opens: the kit's popup, [`MENU_W`] wide, for `add`'s rows.
pub fn context_menu<R>(response: &egui::Response, add: impl FnOnce(&mut Ui) -> R) -> Option<R> {
    kit::Menu::context(response, |ui| {
        let t = solidcraft_theme::Tokens::get();
        ui.set_width((MENU_W - t.space_6 * 2.0 - 2.0).max(0.0));
        ui.spacing_mut().item_spacing.y = 0.0;
        add(ui)
    })
    .map(|r| r.inner)
}

/// A menu's width.
pub const MENU_W: f32 = kit::LIST_W;
/// A wide menu's width.
pub const WIDE_MENU_W: f32 = 360.0;

/// A menu row: its label, and its shortcut (or a note) at the right; dimmed and inert when
/// disabled. The kit's menu row; `id` names it in the widget registry (`menu:<id>`).
pub fn menu_row(ui: &mut Ui, id: &str, label: &str, right: Option<&str>, enabled: bool) -> egui::Response {
    let mut item = kit::MenuItem::new(id, label, "");
    item.shortcut = right.unwrap_or("");
    item.enabled = enabled;
    item.show(ui)
}

/// A menu section's small heading.
pub fn menu_heading(ui: &mut Ui, text: &str) {
    kit::Menu::heading(ui, text);
}

/// A rule between a menu's groups.
pub fn menu_separator(ui: &mut Ui) {
    kit::Menu::separator(ui);
}

#[cfg(test)]
mod tests {
    /// Windows and popups are made by the kit only, so none can skip its frames.
    #[test]
    fn windows_come_from_the_kit() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut bad = Vec::new();
        for e in std::fs::read_dir(&dir).unwrap() {
            let p = e.unwrap().path();
            if p.extension().is_some_and(|x| x == "rs") {
                let text = std::fs::read_to_string(&p).unwrap();
                for raw in [concat!("egui::Window", "::new("), concat!("Popup", "::new("), concat!("egui::Modal", "::new(")] {
                    if text.contains(raw) && !(raw.starts_with("Popup") && p.file_name().is_some_and(|n| n == "frame.rs")) {
                        bad.push(format!("{}: {raw}", p.display()));
                    }
                }
            }
        }
        assert!(bad.is_empty(), "use crate::frame (the kit's windows and menus) instead of {bad:?}");
    }
}
