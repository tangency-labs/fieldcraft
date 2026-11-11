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
//! its tooltip says when it comes. The milestone badge and the drop-downs are FieldCraft's.

use egui::{FontId, Rect, Stroke, pos2, vec2};
use fieldcraft_engine::catalogue::{self, Capability, Tab, Workspace};
use fieldcraft_engine::requirement::Unavailable;
use serde::Serialize;
use solidcraft_theme::TabAccent;
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _, app_bar, ribbon};

use crate::ModelApp;
use crate::theme::Tokens;

/// The ribbon's height: the tab row, then the panels with their buttons' labels (the kit's).
pub const HEIGHT: f32 = TABS_HEIGHT + ribbon::GROUP_HEIGHT + 1.0;
/// The tab row's height (the kit's).
const TABS_HEIGHT: f32 = ribbon::TAB_ROW_HEIGHT;
/// Most buttons a panel shows; the rest are in its drop-down.
pub const MAX_BUTTONS: usize = 4;
/// Where the tabs and the panels start, right of the workspace switcher.
const PANELS_X: f32 = 112.0;
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

/// A milestone badge ("M5") with its top right corner at `corner`.
fn badge(p: &egui::Painter, t: &Tokens, corner: egui::Pos2, m: &str) {
    let text = p.layout_no_wrap(m.to_string(), FontId::proportional(9.0), t.text_dim);
    let r = Rect::from_min_size(pos2(corner.x - text.size().x - 4.0, corner.y), text.size() + vec2(4.0, 0.0));
    p.rect_filled(r, 3.0, t.panel_header);
    p.rect_stroke(r, 3.0, Stroke::new(1.0, t.border), egui::StrokeKind::Inside);
    p.galley(r.min + vec2(2.0, 0.0), text, t.text_dim);
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
    egui::Panel::top("fc_toolbar").exact_size(HEIGHT).show_separator_line(false).frame(egui::Frame::NONE.fill(t.toolbar)).show(ui, |ui| {
        let r = ui.max_rect();
        app.place("ribbon", r);
        app.toolbar_rects.clear();
        app.panel_rects.clear();
        workspace_switcher(app, ui, r);
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
                let widths: Vec<f32> = shown.iter().map(|c| ribbon::button_width(&painter, c.label)).collect();
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
        let out = ribbon::show(ui, &row, Rect::from_min_max(pos2(r.left() + PANELS_X, top), r.max), &mut |_, _| false);
        // FieldCraft's part of each button: its hover text and milestone badge.
        for d in &out.buttons {
            let Some(c) = shown.get(d.group).and_then(|s| s.get(d.index)) else { continue };
            let s = panels.get(d.group).and_then(|p| p.iter().position(|e| e.id == c.id)).and_then(|i| states.get(d.group)?.get(i));
            if let Some(m) = c.milestone() {
                badge(ui.painter(), &t, d.whole.right_top() + vec2(-2.0, 2.0), m);
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
            if let Some(rect) = crate::frame::wide_menu(app, ui, key, anchor, |app, ui| dropdown(app, ui, entries)) {
                app.place("toolbar_dropdown", rect);
            }
        }
        painter.line_segment([pos2(r.left(), r.bottom() - 0.5), pos2(r.right(), r.bottom() - 0.5)], Stroke::new(1.0, t.border));
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

/// The workspace switcher: the kit's switcher at the ribbon's left end, spanning its height; it
/// opens the kit's workspace menu.
fn workspace_switcher(app: &mut ModelApp, ui: &mut egui::Ui, r: Rect) {
    let ws = Rect::from_min_max(r.left_top(), pos2(r.left() + PANELS_X - 4.0, r.bottom()));
    let w = app.ui.workspace;
    let switcher = app_bar::Switcher {
        name: w.label(),
        caption: "Workspace",
        icon: w.icon(),
        tip: "Change workspace",
        open: app.ui.menu.as_deref() == Some("workspace"),
        floating: false,
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
    if let Some(Some(item)) = picked
        && let Some(w) = Workspace::ALL.into_iter().find(|w| w.label() == item.id)
    {
        app.ui.select_workspace(w);
        app.ui.menu = None;
    }
}

/// The rows of a panel's drop-down: every entry with its related commands under it, coming-soon
/// entries greyed with their milestone. A click starts the row's command.
pub fn dropdown(app: &mut ModelApp, ui: &mut egui::Ui, entries: &'static [Capability]) {
    let mut start = None;
    for c in entries {
        let missing = app.unavailable(c.id);
        let enabled = c.command().is_some() && missing.is_none();
        let soon = c.milestone().map(|m| format!("soon · {m}"));
        let mut row = solidcraft_ui_kit::MenuItem::new(c.id, c.label, c.id);
        row.enabled = enabled;
        row.shortcut = soon.as_deref().unwrap_or("");
        let resp = row.show(ui);
        app.toolbar_rects.push((c.id, resp.rect, Fit::Menu));
        if resp.on_hover_text(Unavailable::tip(&c.tooltip(), missing.as_ref())).clicked() {
            start = Some(c.id);
        }
        // Its related commands (edit, delete…) follow it, without an icon.
        for id in c.more {
            let missing = app.unavailable(id);
            let mut row = solidcraft_ui_kit::MenuItem::new(id, crate::command_label(id), "");
            row.enabled = missing.is_none();
            let tip = Unavailable::tip(id, missing.as_ref());
            if row.show(ui).on_hover_text(tip).clicked() {
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
