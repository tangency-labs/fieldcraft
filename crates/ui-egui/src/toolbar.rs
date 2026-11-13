//! The ribbon: the workspace switcher, the workspace's tabs, and the current tab's panels, all
//! from the capability catalogue (`fieldcraft_engine::catalogue`), drawn by the kit's ribbon
//! (`solidcraft_ui_kit::ribbon`) as SolidCraft's is.
//!
//! Each panel shows its ready entries as labelled buttons (up to [`MAX_BUTTONS`]; a panel with
//! none ready shows its first entry, greyed with its milestone) over its label. The label opens
//! the panel's drop-down, which lists every entry of the panel: ready ones with their related
//! commands (edit, delete…) under them, coming-soon ones greyed with their milestone. When the
//! window is too narrow the kit's fitting narrows the buttons, hands them to the drop-down,
//! collapses panels into one button and moves the last ones into its » menu. Nothing is dropped.
//!
//! A ready entry starts its command ([`ModelApp::start`]); a coming-soon entry does nothing and
//! its tooltip says when it comes. The milestone badge is FieldCraft's; the drop-downs hang from
//! their panel labels as SolidCraft's do (`crate::frame::surface_popup`), without its pinning.

use egui::{FontId, Rect, pos2, vec2};
use fieldcraft_engine::catalogue::{self, Capability, Tab, Workspace};
use fieldcraft_engine::requirement::Unavailable;
use serde::Serialize;
use solidcraft_theme::TabAccent;
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _, app_bar, ribbon};

use crate::ModelApp;
use crate::theme::Tokens;

/// The ribbon card's height: the tab row, then the panels with their buttons' labels (the kit's;
/// SolidCraft's `RIBBON_HEIGHT`).
pub const CARD_HEIGHT: f32 = TABS_HEIGHT + ribbon::GROUP_HEIGHT + 1.0;
/// The band the cards take under the app bar: the kit's shell gap above them, then the cards.
/// What is below them (the browser card) starts its own gap under them.
pub const HEIGHT: f32 = kit::SHELL_GAP + CARD_HEIGHT;
/// The tab row's height (the kit's).
const TABS_HEIGHT: f32 = ribbon::TAB_ROW_HEIGHT;
/// Most buttons a panel shows; the rest are in its drop-down.
pub const MAX_BUTTONS: usize = 4;
/// The workspace switcher card's width (SolidCraft's `WORKSPACE_WIDTH`).
const WORKSPACE_WIDTH: f32 = 88.0;
/// Room a group's label needs beside its text: the drop-down arrow and the group's insets.
const LABEL_ROOM: f32 = 30.0;

/// Where the ribbon publishes its parts in the widget registry.
const HANDLES: ribbon::Handles<'static> = ribbon::Handles {
    button: &["ribbon:"],
    button_rect: "ribbon_rect:",
    split: "ribbon_split:",
    cell: &["ribbon_cell:"],
    collapsed: "ribbon_collapsed:",
    section: "ribbon_section:",
    section_target: "ribbon_section_rect:",
    group: "ribbon_group:",
    running: "ribbon_running:",
    reveal: "ribbon_reveal",
    more: "ribbon_more",
    more_row: "ribbon_more:",
    collapsed_count: "ribbon_collapsed_count",
    more_count: "ribbon_more_count",
    trimmed_count: "ribbon_trimmed_count",
    drag: "ribbon_drag",
    drag_remove: "ribbon_drag_remove",
};

/// Where an entry is on a toolbar that may be too narrow for every panel (`ui.toolbar`'s `fit`).
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

/// The kit's ribbon icon's side (it is centred across the button).
const ICON_SIDE: f32 = 26.0;

/// Where a milestone badge `size` goes on the button `whole`: in its top-right corner, moved right
/// of the icon when the button is too narrow for both (it never covers the icon).
fn badge_rect(whole: Rect, size: egui::Vec2) -> Rect {
    let k = solidcraft_theme::Tokens::get();
    let icon_right = whole.center().x + ICON_SIDE / 2.0;
    let left = (whole.right() - k.space_2 - size.x).max(icon_right + k.space_2);
    Rect::from_min_size(pos2(left, whole.top() + k.space_2), size)
}

/// A milestone badge ("M5") on the button `whole` ([`badge_rect`]).
fn badge(p: &egui::Painter, t: &Tokens, whole: Rect, m: &str) {
    let k = solidcraft_theme::Tokens::get();
    let text = p.layout_no_wrap(m.to_string(), FontId::proportional(k.font_badge), t.text_dim);
    let pad = k.space_2;
    let r = badge_rect(whole, text.size() + vec2(2.0 * pad, 0.0));
    p.rect_filled(r, k.radius_sm, t.panel_header);
    p.rect_stroke(r, k.radius_sm, k.primary_border(t.border), egui::StrokeKind::Inside);
    p.galley(r.min + vec2(pad, 0.0), text, t.text_dim);
}

/// Widest a button grows to show its label whole: twice the kit's widest.
const WHOLE_LABEL_MAX: f32 = 2.0 * ribbon::BUTTON_MAX_WIDTH;

/// Whether `label` shows whole (two lines at most, nothing cut) on a button `width` wide.
fn label_whole(painter: &egui::Painter, label: &str, width: f32) -> bool {
    let shown: String = ribbon::button_lines(painter, label, width - ribbon::BUTTON_LABEL_PADDING).concat();
    shown.chars().filter(|c| !c.is_whitespace()).eq(label.chars().filter(|c| !c.is_whitespace()))
}

/// A button's width: the kit's, widened only as far as its label needs to show whole on its two
/// lines (FieldCraft's labels run longer than the kit's widest button). The ribbon still narrows
/// buttons, and cuts their labels, when the window has no room for them.
fn button_width(painter: &egui::Painter, label: &str) -> f32 {
    let mut width = ribbon::button_width(painter, label);
    while width < WHOLE_LABEL_MAX && !label_whole(painter, label, width) {
        width = (width + 2.0).min(WHOLE_LABEL_MAX);
    }
    width
}

/// One entry as the ribbon shows it: can it run now, is its form open, and its hover text.
struct State {
    enabled: bool,
    running: bool,
    tip: String,
}

fn state(app: &ModelApp, c: &Capability) -> State {
    let missing = app.unavailable(c.id);
    State {
        enabled: c.command().is_some() && missing.is_none(),
        running: app.dialog.as_ref().is_some_and(|d| c.command() == Some(d.command)),
        tip: Unavailable::tip(&c.tooltip(), missing.as_ref()),
    }
}

pub fn toolbar(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let mut start: Option<&'static str> = None;
    // The workspace switcher and the ribbon are two chrome cards over the viewport, as
    // SolidCraft's floating shell lays them out (`toolbar::toolbar`, `ui_kit::Cards`): `ui`'s
    // available space is the window below the application bar, which the viewport fills.
    let host = kit::Container::new(ui.ctx(), kit::Surface::Chrome);
    let available = ui.available_rect_before_wrap();
    let cards = kit::Cards::new(available, kit::LayoutMode::current(ui.ctx()), CARD_HEIGHT, WORKSPACE_WIDTH);
    // The band the cards take, from the application bar down to the cards' bottom.
    let band = Rect::from_min_max(available.min, pos2(available.right(), cards.ribbon.bottom()));
    let frame = host.frame(ui.ctx(), t.toolbar, kit::frame::Side::Left);
    app.place("ribbon", band);
    app.place("ribbon:card", cards.ribbon);
    app.shell.cards = vec![cards.workspace, cards.ribbon];
    app.toolbar_rects.clear();
    app.panel_rects.clear();
    // One area per card, each exactly its card (an area takes the clicks over all of its rect):
    // the gaps around the cards are the view's, as SolidCraft shields only its cards from the
    // view's input (`shell::covered_rects`).
    host.area(egui::Id::unique("fc_toolbar_workspace")).fixed_pos(cards.workspace.min).show(ui.ctx(), |ui| {
        ui.set_min_size(cards.workspace.size());
        host.paint_background(ui.painter(), frame, cards.workspace);
        workspace_switcher(app, ui, cards.workspace);
    });
    host.area(egui::Id::unique("fc_toolbar")).fixed_pos(cards.ribbon.min).show(ui.ctx(), |ui| {
        ui.set_min_size(cards.ribbon.size());
        host.paint_background(ui.painter(), frame, cards.ribbon);
        let r = cards.ribbon;
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
            floating: host.floating(),
            underline: None,
            handles: ribbon::TabHandles { tab: "tab:", tab_rect: "tab_rect:", add: "tab_add", action: "tab_action:" },
        };
        let row_rect = Rect::from_min_max(r.left_top(), pos2(r.right(), r.top() + TABS_HEIGHT));
        for action in ribbon::tabs(ui, &row, row_rect).actions {
            if let ribbon::TabAction::Select { index, .. } = action
                && let Some(&tab) = tabs.get(index)
            {
                app.ui.tab = tab;
                app.ui.menu = None;
            }
        }
        // Panels: the kit's ribbon row, one group per panel.
        let tab = app.ui.tab;
        let top = r.top() + TABS_HEIGHT;
        let panels: Vec<&'static [Capability]> = catalogue::panels(tab).collect();
        let shown: Vec<Vec<&'static Capability>> = panels.iter().map(|p| buttons(p)).collect();
        let keys: Vec<String> = panels.iter().map(|p| panel_key(tab, p.first().map_or("", |c| c.panel))).collect();
        let states: Vec<Vec<State>> = panels.iter().map(|p| p.iter().map(|c| state(app, c)).collect()).collect();
        let painter = ui.painter().clone();
        let groups = panels
            .iter()
            .zip(&shown)
            .zip(keys.iter().zip(&states))
            .map(|((entries, shown), (key, states))| {
                let panel = entries.first().map_or("", |c| c.panel);
                let of = |c: &Capability| entries.iter().position(|e| e.id == c.id).and_then(|i| states.get(i));
                // A group is as wide as its buttons; its label (and its arrow) must fit under them,
                // so a short group's buttons widen to its label ("Boundary Conditions" over one
                // button would run into the next group's label).
                let widths: Vec<f32> = shown.iter().map(|c| button_width(&painter, c.label)).collect();
                let k = solidcraft_theme::Tokens::get();
                let label = painter.layout_no_wrap(panel.to_string(), kit::type_face::semibold(k.font_caption), k.text_dim).size().x;
                let need = label + LABEL_ROOM;
                let row = widths.iter().sum::<f32>() + widths.len().saturating_sub(1) as f32 * ribbon::BUTTON_GAP + 8.0;
                let extra = if widths.is_empty() { 0.0 } else { (need - row).max(0.0) / widths.len() as f32 };
                ribbon::Group {
                    key,
                    name: panel,
                    label: panel,
                    tip: panel,
                    buttons: shown
                        .iter()
                        .zip(&widths)
                        .map(|(c, width)| ribbon::Button {
                            id: c.id,
                            run: c.id,
                            label: c.label,
                            name: c.label,
                            icon: c.id,
                            size: ribbon::ButtonSize { width: width + extra, min_width: ribbon::BUTTON_MIN_WIDTH + extra },
                            enabled: of(c).is_some_and(|s| s.enabled),
                            pressed: of(c).is_some_and(|s| s.running),
                            exit: false,
                            split: None,
                            app_paints: false,
                            covers: &[],
                        })
                        .collect(),
                    cells: Vec::new(),
                    commands: entries.iter().map(|c| c.id).collect(),
                    icon: shown.first().map(|c| c.id),
                    enabled: true,
                    open: app.ui.menu.as_deref() == Some(key.as_str()),
                    running: None,
                    more: entries
                        .iter()
                        .zip(states)
                        .map(|(c, s)| ribbon::MoreRow {
                            id: c.id,
                            label: c.label,
                            icon: c.id,
                            shortcut: "",
                            enabled: s.enabled,
                            checked: None,
                            hint: Some(s.tip.as_str()),
                        })
                        .collect(),
                }
            })
            .collect();
        let row = ribbon::Row { groups, labels: true, accent, exit: false, reveal: None, more_tip: "More toolbar panels", handles: HANDLES };
        let out = ribbon::show(ui, &row, Rect::from_min_max(pos2(r.left(), top), r.max), &mut |_, _| false);
        // FieldCraft's part of each button: its hover text and milestone badge.
        for d in &out.buttons {
            let Some(c) = shown.get(d.group).and_then(|s| s.get(d.index)) else { continue };
            let s = panels.get(d.group).and_then(|p| p.iter().position(|e| e.id == c.id)).and_then(|i| states.get(d.group)?.get(i));
            if let Some(m) = c.milestone() {
                badge(ui.painter(), &t, d.whole, m);
            }
            if let Some(s) = s {
                let _ = d.response.clone().backing_on_hover_text(&s.tip);
            }
            app.toolbar_rects.push((c.id, d.whole, Fit::Full));
        }
        let mut anchors: Vec<(usize, Rect)> = Vec::new();
        for g in &out.groups {
            let (Some(entries), Some(key)) = (panels.get(g.group), keys.get(g.group)) else { continue };
            app.panel_rects.push((key.clone(), g.label));
            let (rect, fit) = if g.collapsed { (g.rect, Fit::Collapsed) } else { (Rect::NOTHING, Fit::Menu) };
            let rest: Vec<&'static str> = entries.iter().map(|c| c.id).filter(|id| !app.toolbar_rects.iter().any(|(on, ..)| on == id)).collect();
            app.toolbar_rects.extend(rest.into_iter().map(|id| (id, rect, fit)));
            anchors.push((g.group, if g.collapsed { g.rect } else { g.label }));
        }
        for (gi, entries) in panels.iter().enumerate() {
            if !out.groups.iter().any(|g| g.group == gi) {
                for c in entries.iter() {
                    app.toolbar_rects.push((c.id, Rect::NOTHING, Fit::More));
                }
            }
        }
        for action in &out.actions {
            match action {
                ribbon::Action::Run { id, .. } => start = panels.iter().flat_map(|p| p.iter()).find(|c| c.id == id).map(|c| c.id),
                ribbon::Action::Section { group } | ribbon::Action::Open { group } => {
                    if let Some(key) = keys.get(*group) {
                        app.toggle_menu(key);
                    }
                }
                ribbon::Action::Menu { .. } | ribbon::Action::Drop { .. } => {}
            }
        }
        for (group, anchor) in anchors {
            let (Some(entries), Some(key)) = (panels.get(group), keys.get(group)) else { continue };
            // SolidCraft's panel drop-down: the kit's menu surface hanging from the panel's label.
            let shown = crate::frame::surface_popup(app, ui, key, anchor, |app, ui| {
                kit::Menu::show_frame(ui, kit::Menu::frame(), |ui| dropdown(app, ui, entries))
            });
            if let Some(shown) = shown {
                app.place("toolbar_dropdown", shown.response.rect);
            }
        }
        for card in [cards.workspace, cards.ribbon] {
            kit::Container::paint_inset(ui.painter(), frame.widget_rect(card), frame);
        }
    });
    if let Some(id) = start {
        app.ui.menu = None;
        // Errors already show in the status bar.
        let _ = app.start(id);
    }
}

/// The message ids of a workspace's line and tab list in the kit's workspace menu; their text is
/// FieldCraft's ([`Workspace::about`] and its tabs, `crate::theme`'s catalogue).
pub fn workspace_ids(w: Workspace) -> (&'static str, &'static str) {
    match w {
        Workspace::Prepare => ("workspace.prepare.about", "workspace.prepare.tabs"),
        Workspace::Model => ("workspace.model.about", "workspace.model.tabs"),
        Workspace::Analysis => ("workspace.analysis.about", "workspace.analysis.tabs"),
        Workspace::Results => ("workspace.results.about", "workspace.results.tabs"),
    }
}

/// The workspace switcher: the kit's switcher on its card `ws`, left of the ribbon's; it opens the
/// kit's workspace menu.
fn workspace_switcher(app: &mut ModelApp, ui: &mut egui::Ui, ws: Rect) {
    let floating = kit::Container::new(ui.ctx(), kit::Surface::Chrome).floating();
    let w = app.ui.workspace;
    let switcher = app_bar::Switcher {
        name: w.label(),
        caption: "Workspace",
        icon: w.icon(),
        tip: "Change workspace",
        open: app.ui.menu.as_deref() == Some("workspace"),
        floating,
        labels: true,
        folded: false,
        stage: false,
        accent: solidcraft_theme::Tokens::get().accent,
        handle_rect: "workspace_switcher",
        handle_target: "ribbon:workspace",
    };
    app.place("workspace_switcher", ws);
    if app_bar::switcher(ui, &switcher, ws) {
        app.toggle_menu("workspace");
    }
    let current = app.ui.workspace;
    let items: Vec<kit::Item> = Workspace::ALL
        .into_iter()
        .map(|w| {
            let mut item = kit::Item::action(w.label(), w.label(), w.icon());
            if w == current {
                item.shortcut = "current".into();
            }
            item
        })
        .collect();
    let style = |item: &kit::Item| {
        let w = Workspace::ALL.into_iter().find(|w| w.label() == item.id)?;
        let (description, tabs) = workspace_ids(w);
        Some(kit::WorkspaceStyle { slug: w.label(), description, tabs: Some(tabs), ramp: None })
    };
    let height = ui.ctx().content_rect().bottom() - ws.bottom();
    let picked = crate::frame::surface_popup(app, ui, "workspace", ws, |_, ui| {
        kit::workspace_menu(
            ui,
            &kit::WorkspaceMenu { items: &items, style: &style, footer: solidcraft_ui_text::tr("ui.workspace_footer") },
            None,
            height,
        )
    });
    if let Some(Some(item)) = picked.map(|p| p.inner)
        && let Some(w) = Workspace::ALL.into_iter().find(|w| w.label() == item.id)
    {
        app.ui.select_workspace(w);
        app.ui.menu = None;
    }
}

/// One row of a panel's drop-down.
struct Row {
    id: &'static str,
    label: &'static str,
    /// A kit icon name ("" for none).
    icon: &'static str,
    /// The milestone of a coming-soon entry, where SolidCraft shows the shortcut.
    note: String,
    enabled: bool,
    /// Its hover text: what it does, or what its command is missing.
    tip: String,
    /// A catalogue entry (not one of its related commands).
    entry: bool,
}

/// The rows of a panel's drop-down: every entry with its related commands (edit, delete…) under
/// it without an icon, coming-soon entries greyed with their milestone.
fn rows(app: &ModelApp, entries: &'static [Capability]) -> Vec<Row> {
    let mut rows = Vec::new();
    for c in entries {
        let missing = app.unavailable(c.id);
        rows.push(Row {
            id: c.id,
            label: c.label,
            icon: c.id,
            note: c.milestone().map(|m| format!("soon · {m}")).unwrap_or_default(),
            enabled: c.command().is_some() && missing.is_none(),
            tip: Unavailable::tip(&c.tooltip(), missing.as_ref()),
            entry: true,
        });
        for id in c.more {
            let missing = app.unavailable(id);
            rows.push(Row {
                id,
                label: crate::command_label(id),
                icon: "",
                note: String::new(),
                enabled: missing.is_none(),
                tip: Unavailable::tip(id, missing.as_ref()),
                entry: false,
            });
        }
    }
    rows
}

/// Where the drop-downs' rows keep their keyboard focus ([`ribbon::list_keys`]); FieldCraft's
/// ribbon has no pinning, so the pin's parts are unused.
fn menu_style(ctx: &egui::Context) -> ribbon::MenuStyle<'static> {
    ribbon::MenuStyle {
        pin_tip: "",
        unpin_tip: "",
        accent: TabAccent::current(ctx),
        focus: egui::Id::unique("fc_dropdown_row"),
        entry: &["dropdown:"],
        pin: "dropdown_pin:",
        hot: "dropdown_hot:",
        reveal: "ribbon:reveal",
    }
}

/// A drop-down row as SolidCraft's ([`ribbon::menu_row`]) without its pin and drag: the kit's
/// menu row [`ribbon::ROW_H`] high and `width` wide, focused by [`ribbon::list_keys`], its hover
/// text beside the menu ([`ribbon::menu_tip`]).
fn menu_row(ui: &mut egui::Ui, row: &Row, width: f32, style: &ribbon::MenuStyle) -> egui::Response {
    let (rect, _) = ui.allocate_exact_size(vec2(width, ribbon::ROW_H), egui::Sense::hover());
    let resp = ui.interact(rect, style.focus_id(row.id), if row.enabled { egui::Sense::click() } else { egui::Sense::hover() });
    if resp.has_focus() && !ui.clip_rect().contains_rect(rect) {
        resp.scroll_to_me_animation(None, egui::style::ScrollAnimation::none());
    }
    let mut item = kit::MenuItem::new(row.id, row.label, row.icon);
    item.enabled = row.enabled;
    item.shortcut = &row.note;
    item.highlighted = row.enabled && (resp.hovered() || resp.has_focus());
    item.paint(ui, rect, &resp, 0.0);
    // The menu's rect from the last pass: the tip sits beside it, never over its rows.
    let menu = egui::containers::AreaState::load(ui.ctx(), ui.layer_id().id).map_or(rect, |s| s.rect());
    ribbon::menu_tip(ui, menu, &resp, |ui| {
        ui.label(&row.tip);
    });
    resp
}

/// A panel's drop-down, as SolidCraft's: its [`rows`] at the kit's menu width, Up and Down through
/// the enabled ones, scrolling when the window is too short. A click (or Enter) starts the row's
/// command and closes it.
pub fn dropdown(app: &mut ModelApp, ui: &mut egui::Ui, entries: &'static [Capability]) {
    let rows = rows(app, entries);
    let style = menu_style(ui.ctx());
    ui.spacing_mut().item_spacing.y = 0.0;
    let width = ribbon::menu_width(ui, rows.iter().map(|r| (r.label, r.note.as_str())));
    ui.set_width(width);
    let enabled: Vec<&str> = rows.iter().filter(|r| r.enabled).map(|r| r.id).collect();
    ribbon::list_keys(ui, &enabled, &style);
    let room = (ui.ctx().content_rect().bottom() - ui.cursor().top() - 64.0).max(ribbon::ROW_H * 3.0);
    let mut start = None;
    egui::ScrollArea::vertical().id_salt("fc_panel_menu_rows").max_height(room).min_scrolled_height(room).show(ui, |ui| {
        for row in &rows {
            let resp = menu_row(ui, row, width, &style);
            if row.entry {
                app.toolbar_rects.push((row.id, resp.rect, Fit::Menu));
            }
            if resp.clicked() && row.enabled {
                start = Some(row.id);
            }
        }
    });
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

    /// Every ribbon label shows whole on its button while the toolbar has room (A3-2).
    #[test]
    fn ribbon_labels_show_whole() {
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let mut cut = Vec::new();
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            let painter = ui.painter();
            for c in catalogue::CATALOGUE.iter() {
                if !label_whole(painter, c.label, button_width(painter, c.label)) {
                    cut.push(c.label);
                }
            }
            assert!(!label_whole(painter, "Fluid Measurements", ribbon::BUTTON_MAX_WIDTH), "the kit's widest cuts it");
        });
        out.textures_delta.clear();
        assert!(cut.is_empty(), "cut labels: {cut:?}");
    }

    /// A panel's drop-down is SolidCraft's: the kit's menu hanging from the panel's label, its rows
    /// at the kit's menu width, Down focusing its first enabled row.
    #[test]
    fn a_drop_down_hangs_from_its_panel_label_with_ribbon_rows() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.ui.select_tab(Tab::Mesh);
        let key = panel_key(Tab::Mesh, "Tet");
        app.ui.menu = Some(key.clone());
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let screen = Rect::from_min_size(pos2(0.0, 0.0), vec2(1400.0, 900.0));
        let run = |app: &mut ModelApp, events: Vec<egui::Event>| {
            let raw = egui::RawInput { screen_rect: Some(screen), events, ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| toolbar(app, ui));
            out.textures_delta.clear();
        };
        run(&mut app, Vec::new());
        run(&mut app, Vec::new());
        let menu = app.layout.get("toolbar_dropdown").copied().unwrap_or_else(|| panic!("drop-down not drawn"));
        let label = app.panel_rects.iter().find(|(k, _)| *k == key).map(|(_, r)| *r).unwrap_or_else(|| panic!("no Tet label"));
        // Within the pixel the area is rounded to.
        assert!((menu.left_top() - label.left_bottom()).length() < 0.5, "menu {menu:?}, label {label:?}");
        let entries = catalogue::panels(Tab::Mesh).find(|p| p.first().is_some_and(|c| c.panel == "Tet")).unwrap_or_default();
        let rows = rows(&app, entries);
        let mut width = 0.0;
        let _ = ctx.run_ui(egui::RawInput::default(), |ui| width = ribbon::menu_width(ui, rows.iter().map(|r| (r.label, r.note.as_str()))));
        let margin = kit::Menu::frame().total_margin().sum().x;
        assert!((menu.width() - width - margin).abs() < 1.0, "menu {} wide, rows {width} + {margin}", menu.width());
        for c in entries.iter() {
            let row = app.toolbar_rects.iter().rev().find(|(id, _, f)| *id == c.id && *f == Fit::Menu).map(|(_, r, _)| *r);
            assert!(row.is_some_and(|r| (r.height() - ribbon::ROW_H).abs() < 0.5 && menu.contains_rect(r)), "{}: {row:?}", c.id);
        }
        let first = rows.iter().find(|r| r.enabled).map(|r| r.id).unwrap_or_else(|| panic!("no enabled row"));
        let down =
            egui::Event::Key { key: egui::Key::ArrowDown, physical_key: None, pressed: true, repeat: false, modifiers: egui::Modifiers::default() };
        run(&mut app, vec![down]);
        run(&mut app, Vec::new());
        assert_eq!(ctx.memory(|m| m.focused()), Some(menu_style(&ctx).focus_id(first)));
        assert_eq!(app.ui.menu.as_deref(), Some(key.as_str()));
    }

    /// The workspace switcher and the ribbon are SolidCraft's two cards (`ui_kit::Cards`): inset
    /// from the window's edges and from the app bar by the kit's shell gap, the switcher
    /// SolidCraft's width and the ribbon a gap right of it, both the ribbon's height. The band they
    /// take ends at the cards, so what is below starts a gap under them; the viewport runs under
    /// them, so they are recorded as the shell's cards.
    #[test]
    fn the_switcher_and_the_ribbon_are_inset_cards() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let screen = Rect::from_min_size(pos2(0.0, 0.0), vec2(1400.0, 900.0));
        for _ in 0..2 {
            let raw = egui::RawInput { screen_rect: Some(screen), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| toolbar(&mut app, ui));
            out.textures_delta.clear();
        }
        let at = |k: &str| app.layout.get(k).copied().unwrap_or_else(|| panic!("{k} not drawn"));
        let (band, ws, card) = (at("ribbon"), at("workspace_switcher"), at("ribbon:card"));
        let gap = kit::SHELL_GAP;
        let near = |a: f32, b: f32| (a - b).abs() < 0.5;
        assert!(near(band.top(), 0.0) && near(band.height(), HEIGHT) && near(band.bottom(), card.bottom()), "{band:?} {card:?}");
        assert!(near(ws.left(), gap) && near(ws.top(), gap) && near(ws.width(), WORKSPACE_WIDTH), "{ws:?}");
        assert!(near(card.left(), ws.right() + gap) && near(card.top(), gap) && near(card.right(), screen.right() - gap), "{card:?}");
        assert!(near(card.height(), CARD_HEIGHT) && near(ws.bottom(), card.bottom()), "{ws:?} {card:?}");
        assert_eq!(app.shell.cards, vec![ws, card], "the viewport keeps its overlays clear of both cards");
        let buttons: Vec<Rect> = app.toolbar_rects.iter().filter(|(.., f)| *f == Fit::Full).map(|(_, r, _)| *r).collect();
        assert!(!buttons.is_empty() && buttons.iter().all(|r| card.contains_rect(*r)), "{buttons:?} outside {card:?}");
    }

    /// A milestone badge keeps clear of its button's icon, even on the narrowest button (A3-15).
    #[test]
    fn badges_keep_off_the_icon() {
        // The kit centres the icon this far below the button's top.
        const ICON_CENTRE: f32 = 16.0;
        for width in [ribbon::BUTTON_MIN_WIDTH, 70.0, ribbon::BUTTON_MAX_WIDTH, 120.0] {
            let whole = Rect::from_min_size(pos2(100.0, 40.0), vec2(width, ribbon::BUTTON_HEIGHT));
            let icon = Rect::from_center_size(pos2(whole.center().x, whole.top() + ICON_CENTRE), vec2(ICON_SIDE, ICON_SIDE));
            let badge = badge_rect(whole, vec2(20.0, 12.0));
            assert!(!badge.intersects(icon), "{width}: badge {badge:?}, icon {icon:?}");
        }
    }
}
