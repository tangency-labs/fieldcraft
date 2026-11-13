//! The timeline (bottom): the command journal in order, as SolidCraft shows its features. Runs of
//! the same command fold into one chip with their count, which a click opens to show each command
//! (as SolidCraft folds a long run of one component's features); the marker stands after the
//! current state, and undone commands follow it greyed, ready to redo. Clicking a command (or the
//! transport buttons) undoes or redoes up to it with `edit.undo` / `edit.redo`, as dragging the
//! marker does; commands older than the undo history cannot be reached and say so.
//!
//! The bar is the kit's History bar (`solidcraft_ui_kit::history_bar`) and its run chip the kit's
//! `HistoryGroupChip`, as SolidCraft draws its own; FieldCraft keeps its journal and undo behind
//! them.

use std::collections::BTreeSet;

use egui::{Align2, Color32, FontId, Rect, Stroke, pos2, vec2};
use serde_json::{Value, json};
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _, history_bar};

use crate::ModelApp;
use crate::theme::Tokens;

/// The timeline's height.
pub const HEIGHT: f32 = 40.0;

/// A run of the same command in the journal (`start..end`), or among the undone commands
/// (positions past the journal's end).
#[derive(Clone, Debug, PartialEq)]
pub struct Item {
    pub command: String,
    /// Journal positions: after this item the journal holds `end` commands.
    pub start: usize,
    pub end: usize,
    /// Each command's parameters, in order (the tooltips).
    pub params: Vec<Value>,
}

impl Item {
    pub fn count(&self) -> usize {
        self.end - self.start
    }

    /// The command at journal position `at` of this run, as an item of its own.
    fn one(&self, at: usize) -> Item {
        let params = at.checked_sub(self.start).and_then(|j| self.params.get(j)).cloned().into_iter().collect();
        Item { command: self.command.clone(), start: at, end: at + 1, params }
    }
}

/// The journal and the undone commands as runs of the same command, rebuilt only when the journal
/// or the undo history changed.
fn cached_items(app: &mut ModelApp) -> (std::sync::Arc<Vec<Item>>, usize) {
    let key = (app.session.revision(), app.session.journal.len(), app.session.undo_depth());
    match &app.timeline {
        Some((k, list)) if *k == key => (list.clone(), key.1),
        _ => {
            let (list, done) = items(app);
            let list = std::sync::Arc::new(list);
            app.timeline = Some((key, list.clone()));
            (list, done)
        }
    }
}

/// The journal and the undone commands as runs of the same command.
pub fn items(app: &ModelApp) -> (Vec<Item>, usize) {
    let done = app.session.journal.len();
    let all = app.session.journal.iter().chain(app.session.redo_entries());
    let mut out: Vec<Item> = Vec::new();
    for (i, (command, params)) in all.enumerate() {
        // A run never crosses the current position.
        match out.last_mut() {
            Some(last) if last.command == *command && last.end == i && i != done => {
                last.end = i + 1;
                last.params.push(params.clone());
            }
            _ => out.push(Item { command: command.clone(), start: i, end: i + 1, params: vec![params.clone()] }),
        }
    }
    (out, done)
}

impl ModelApp {
    /// Undo or redo until the journal holds `target` commands (as far as the history reaches):
    /// one `edit.undo` or `edit.redo` with a `count`.
    pub fn roll_to(&mut self, target: usize) {
        let (undo, redo) = self.session.undo_depth();
        let done = self.session.journal.len();
        let target = target.clamp(done.saturating_sub(undo), done + redo);
        let step = if target < done { "edit.undo" } else { "edit.redo" };
        if target != done {
            // An error shows in the status bar.
            let _ = self.run(step, json!({"count": done.abs_diff(target)}));
        }
    }
}

/// What an item's commands did, in words: "Delete reference point RP1", "Add Node × 3", "Open
/// Sample cantilever". Never the parameters' JSON.
fn describe(item: &Item) -> String {
    let name = item.params.first().and_then(|p| p.get("name")).and_then(Value::as_str);
    let def = crate::browser::Def::ALL
        .into_iter()
        .find_map(|kind| ["create", "edit", "delete"].into_iter().find(|verb| kind.command(verb) == item.command).map(|verb| (kind, verb)));
    let text = match (def, name) {
        (Some((kind, verb)), Some(name)) => {
            let did = match verb {
                "create" => "Add",
                "edit" => "Edit",
                _ => "Delete",
            };
            format!("{did} {} {name}", kind.label().to_lowercase())
        }
        (_, Some(name)) => format!("{} {name}", crate::command_label(&item.command)),
        _ => crate::command_label(&item.command).to_string(),
    };
    if item.count() > 1 { format!("{text} × {}", item.count()) } else { text }
}

/// What a command does to the item its icon shows, when it is not making it: deleting or changing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Badge {
    Delete,
    Edit,
}

fn badge(command: &str) -> Option<Badge> {
    let verb = command.rsplit('.').next()?;
    match verb {
        "delete" | "clear" => Some(Badge::Delete),
        "edit" => Some(Badge::Edit),
        _ => None,
    }
}

/// The badge at an item's top-left corner: a minus on the error colour for a deletion, a pencil
/// stroke on the accent for an edit.
fn paint_badge(p: &egui::Painter, item: Rect, b: Badge, t: &Tokens) {
    let c = item.left_top() + vec2(5.0, 5.0);
    let (fill, ink) = if b == Badge::Delete { (t.error, Color32::WHITE) } else { (t.accent, t.on_accent) };
    let ink = Stroke::new(1.4, ink);
    p.circle(c, 5.0, fill, Stroke::new(1.0, t.timeline_item));
    match b {
        Badge::Delete => {
            p.line_segment([c - vec2(2.5, 0.0), c + vec2(2.5, 0.0)], ink);
        }
        Badge::Edit => {
            p.line_segment([c + vec2(-2.2, 2.2), c + vec2(2.2, -2.2)], ink);
        }
    }
}

/// Text for an item's tooltip: what it did, and whether a click reaches it.
fn tip(item: &Item, reachable: bool) -> String {
    let mut s = describe(item);
    s.push_str(if reachable { "\nClick to roll the model to here" } else { "\nOlder than the undo history" });
    s
}

/// Where the History bar publishes its grip, marker, steps and controls.
const HANDLES: history_bar::Handles<'static> =
    history_bar::Handles { grip: "timeline:grip", marker: "timeline:marker", step: "timeline:step:", controls: "timeline:" };

/// What the strip shows, in order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Entry {
    /// One command: the item at `item`, the command at journal position `at` of it.
    Step { item: usize, at: usize },
    /// A run of one command (the item at `item`) as one chip; open, its commands follow one by
    /// one, in a tray up to its [`Entry::End`].
    Run { item: usize, open: bool },
    /// The end of an open run's tray.
    End,
}

/// The strip's entries for `list`: a single command is a step; a run of one command is a chip,
/// followed by each of its commands while it is open (`open` holds the first journal position of
/// each open run).
fn entries(list: &[Item], open: &BTreeSet<usize>) -> Vec<Entry> {
    let mut out = Vec::with_capacity(list.len());
    for (item, run) in list.iter().enumerate() {
        if run.count() == 1 {
            out.push(Entry::Step { item, at: run.start });
            continue;
        }
        let is_open = open.contains(&run.start);
        out.push(Entry::Run { item, open: is_open });
        if is_open {
            out.extend((run.start..run.end).map(|at| Entry::Step { item, at }));
            out.push(Entry::End);
        }
    }
    out
}

/// The journal's items as the kit's History strip lays them out. Its slots are journal
/// positions: a step's id is its command's position, and a run's chip covers its commands'.
struct Steps<'a> {
    list: &'a [Item],
    entries: &'a [Entry],
    done: usize,
    oldest: usize,
    t: &'a Tokens,
    accent: solidcraft_theme::TabAccent,
    /// The run (its first position) whose chip was clicked: it opens or folds.
    toggle: Option<usize>,
}

impl Steps<'_> {
    fn step(&self, strip: &mut history_bar::Strip<'_, '_>, item: &Item) {
        let reachable = item.end >= self.oldest;
        let label = describe(item);
        let tip = || tip(item, reachable);
        let step = history_bar::Step {
            id: item.start as u64,
            slot: item.start,
            label: &label,
            icon: icon(item),
            icon_accent: None,
            state: history_bar::StepState {
                rolled_back: item.start >= self.done,
                suppressed: !reachable,
                // The floating card dots the last command done.
                last_done: item.end == self.done,
                ..Default::default()
            },
            error: false,
            warning: false,
            ink: None,
            tip: &tip,
        };
        let chip = strip.step(&step);
        if strip.ui().is_rect_visible(chip)
            && let Some(b) = badge(&item.command)
        {
            paint_badge(strip.painter(), chip, b, self.t);
        }
    }

    /// A run's chip: the kit's group chip (its icon, name and count); a click opens or folds it.
    fn run(&mut self, strip: &mut history_bar::Strip<'_, '_>, item: &Item, open: bool) {
        let (done, oldest, t, accent) = (self.done, self.oldest, self.t, self.accent);
        let mut clicked = false;
        strip.chip(item.start, item.end, open, 0, |ui, p, at| {
            let chip = kit::HistoryGroupChip {
                icon: icon(item),
                name: crate::command_label(&item.command),
                count: item.count(),
                open,
                selected: false,
                hovered: false,
                dim: item.start >= done,
                error: false,
                warning: false,
                accent,
            };
            let rect = Rect::from_min_size(at, chip.size(p));
            let response = ui.interact(rect, ui.scope_id().with(("fc_timeline_run", item.start)), egui::Sense::click());
            solidcraft_ui_registry::publish_handle(&format!("timeline:run:{}", item.start), rect.center());
            if ui.is_rect_visible(rect) {
                chip.paint(p, rect, response.hover_pos());
                if let Some(b) = badge(&item.command) {
                    paint_badge(p, rect, b, t);
                }
            }
            let what = if open { "click to fold them" } else { "click to show each one" };
            let mut text = format!("{}: {} commands, {what}", describe(item), item.count());
            if item.end < oldest {
                text.push_str("\nOlder than the undo history");
            }
            clicked = response.on_hover_cursor(egui::CursorIcon::PointingHand).backing_on_hover_text(text).clicked();
            rect.width()
        });
        if clicked {
            self.toggle = Some(item.start);
        }
    }
}

impl history_bar::Contents for Steps<'_> {
    fn entries(&mut self, strip: &mut history_bar::Strip<'_, '_>) {
        for entry in self.entries {
            match *entry {
                Entry::Step { item, at } => {
                    if let Some(run) = self.list.get(item) {
                        self.step(strip, &run.one(at));
                    }
                }
                Entry::Run { item, open } => {
                    if let Some(run) = self.list.get(item) {
                        self.run(strip, run, open);
                    }
                }
                Entry::End => strip.end_tray(),
            }
        }
    }
}

/// An item's icon: its command's catalogue entry's.
fn icon(item: &Item) -> &'static str {
    fieldcraft_engine::catalogue::entry_of(&item.command).map_or("ui.views", |c| c.id)
}

/// The History bar (bottom), as SolidCraft's: in the floating layout a kit Panel card along the
/// bottom of the viewport, over it (`Container::floating_card(..).bottom_strip_rect`), recorded in
/// `app.shell.history`; in the compact layout a strip docked along the window's bottom edge.
pub fn timeline(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let host = kit::Container::new(ui.ctx(), kit::Surface::Panel);
    let (roll, toggle) = if host.floating() {
        let rect = floating_rect(ui.available_rect_before_wrap());
        // A window too short for the card leaves History out rather than drawing into a
        // negative-size rectangle.
        if !rect.is_finite() || rect.width() < 1.0 || rect.height() < HEIGHT - 1.0 {
            return;
        }
        app.shell.history = Some(rect);
        app.place("timeline", rect);
        let frame = host.frame(ui.ctx(), t.timeline, kit::frame::Side::Left);
        host.area(egui::Id::unique("fc_timeline_card"))
            .fixed_pos(rect.min)
            .show(ui.ctx(), |ui| {
                // The card takes its size in the area (and keeps the pointer off the view under it).
                let (card, _) = ui.allocate_exact_size(rect.size(), egui::Sense::hover());
                kit::Container::paint_panel_frame(ui.painter(), frame, card);
                let mut child = ui.new_child(egui::UiBuilder::new().id_salt("fc_timeline").max_rect(card));
                child.set_clip_rect(card);
                contents(app, &mut child, &t, true)
            })
            .inner
    } else {
        let panel = host.edge("fc_timeline", kit::Edge::Bottom).exact_size(HEIGHT).show_separator_line(false);
        let shown = host.show_edge(ui, panel, host.edge_frame(t.timeline), |ui| {
            let r = ui.max_rect();
            app.place("timeline", r);
            ui.painter().line_segment([r.left_top(), r.right_top()], Stroke::new(1.0, t.border));
            contents(app, ui, &t, false)
        });
        shown.inner
    };
    if let Some(start) = toggle
        && !app.timeline_open.remove(&start)
    {
        app.timeline_open.insert(start);
    }
    if let Some(target) = roll {
        app.roll_to(target);
    }
}

/// Where the floating History card goes in `available` (the viewport below the application
/// bar): the kit's shared bottom strip, a shell gap from the edges.
pub(crate) fn floating_rect(available: Rect) -> Rect {
    kit::Container::floating_card(kit::Surface::Panel).bottom_strip_rect(available, HEIGHT)
}

/// The kit's header (grip, title and the History controls) and strip (a chip per command or run,
/// the marker after the current state) in `ui`'s rect. Returns where to roll the model and the
/// run to open or fold. Its controls, a click on a command and a dragged marker undo or redo
/// through [`ModelApp::roll_to`]; the journal cannot be reordered.
fn contents(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, floating: bool) -> (Option<usize>, Option<usize>) {
    let mut roll: Option<usize> = None;
    let mut toggle: Option<usize> = None;
    let r = ui.max_rect();
    let (undo, redo) = app.session.undo_depth();
    let (list, done) = cached_items(app);
    let oldest = done.saturating_sub(undo);
    let accent = solidcraft_theme::TabAccent::current(ui.ctx());
    let bar = history_bar::Bar {
        title: "History",
        marker_tip: "The model is here: drag to roll it back or forward",
        play: false,
        playing: false,
        play_labels: kit::PlayLabels { play: "Play", stop: "Stop" },
        accent,
        drag_accent: accent,
        // Slots are journal positions: the marker stands after the commands done.
        marker: done,
        steps: list.last().map_or(0, |i| i.end),
        floating,
        // The journal is a record of what was done, in order: its steps never move.
        reorder: false,
        handles: HANDLES,
    };
    let head = history_bar::header(ui, &bar, r);
    let mut actions = head.actions;
    if list.is_empty() {
        let k = solidcraft_theme::Tokens::get();
        ui.painter().text(
            pos2(head.x, r.center().y),
            Align2::LEFT_CENTER,
            "Every command that changes the model appears here",
            FontId::proportional(k.font_label),
            k.text_dim,
        );
    } else {
        let strip = Rect::from_min_max(pos2(head.x, r.top()), r.right_bottom());
        // A run that is no longer one (undone, another model) forgets it was open.
        app.timeline_open.retain(|start| list.iter().any(|i| i.start == *start && i.count() > 1));
        let shown = entries(&list, &app.timeline_open);
        let mut steps = Steps { list: &list, entries: &shown, done, oldest, t, accent, toggle: None };
        actions.extend(history_bar::strip(ui, &bar, r, strip, "fc_timeline_strip", &mut steps).actions);
        toggle = steps.toggle;
    }
    for action in actions {
        roll = match action {
            history_bar::Action::Seek(history_bar::Seek::First) => Some(oldest),
            history_bar::Action::Seek(history_bar::Seek::Back) => Some(done.saturating_sub(1)),
            history_bar::Action::Seek(history_bar::Seek::Forward) => Some(done + 1),
            history_bar::Action::Seek(history_bar::Seek::Last) => Some(done + redo),
            history_bar::Action::Seek(history_bar::Seek::Slot { slot, .. }) => Some(slot),
            // A command's chip rolls the model to just after it, when the history reaches it.
            history_bar::Action::Select { step, .. } => usize::try_from(step).ok().map(|at| at + 1).filter(|end| *end >= oldest),
            _ => roll,
        };
    }
    (roll, toggle)
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;

    use super::*;

    #[test]
    fn runs_of_a_command_are_one_item_and_rolling_undoes_and_redoes() {
        let mut app = ModelApp::new(Session::default());
        for l in 1..=3 {
            app.run("mesh.node.add", json!({"label": l, "xyz": [l, 0, 0]})).unwrap();
        }
        app.run(
            "model.material.create",
            json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}, "density": 7.85e-9}),
        )
        .unwrap();
        let (list, done) = items(&app);
        assert_eq!(done, 4);
        assert_eq!(list.iter().map(|i| (i.command.as_str(), i.count())).collect::<Vec<_>>(), [("mesh.node.add", 3), ("model.material.create", 1)]);
        app.roll_to(2);
        assert_eq!(app.session.journal.len(), 2);
        assert!(app.session.model.node(3).is_none());
        // The undone commands follow the marker, split from the done ones.
        let (list, done) = items(&app);
        assert_eq!(done, 2);
        assert_eq!(list.iter().map(|i| (i.start, i.end)).collect::<Vec<_>>(), [(0, 2), (2, 3), (3, 4)]);
        app.roll_to(4);
        assert!(app.session.model.materials().contains_key("Steel"));
        app.roll_to(100);
        assert_eq!(app.session.journal.len(), 4, "clamped to what can be redone");
    }

    /// A deletion or an edit reads apart from the creation it undoes (#359): its own badge, and
    /// a tooltip in words.
    #[test]
    fn deletions_and_edits_read_apart_from_creations() {
        let item = |command: &str, params: Value| Item { command: command.into(), start: 0, end: 1, params: vec![params] };
        let create = item("model.reference_point.create", json!({"name": "RP1", "xyz": [0, 0, 0]}));
        let delete = item("model.reference_point.delete", json!({"name": "RP1"}));
        assert_eq!((badge(&create.command), badge(&delete.command)), (None, Some(Badge::Delete)));
        assert_eq!(badge("model.load.edit"), Some(Badge::Edit));
        assert_eq!(describe(&create), "Add reference point RP1");
        assert_eq!(describe(&delete), "Delete reference point RP1");
        assert!(!tip(&delete, true).contains('{'));
        let nodes = Item { end: 3, ..item("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})) };
        assert_eq!(describe(&nodes), "Add Node × 3");
    }

    /// A run of one command is the kit's group chip covering its commands' slots, and opened it
    /// lists each command after it, in a tray, as SolidCraft folds a run of features.
    #[test]
    fn runs_fold_into_a_group_chip_that_opens_to_each_command() {
        let mut app = ModelApp::new(Session::default());
        for l in 1..=3 {
            app.run("mesh.node.add", json!({"label": l, "xyz": [l, 0, 0]})).unwrap();
        }
        app.run(
            "model.material.create",
            json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}, "density": 7.85e-9}),
        )
        .unwrap();
        let (list, _) = items(&app);
        let closed = entries(&list, &BTreeSet::new());
        assert_eq!(closed, [Entry::Run { item: 0, open: false }, Entry::Step { item: 1, at: 3 }]);
        let open = entries(&list, &BTreeSet::from([0]));
        assert_eq!(
            open,
            [
                Entry::Run { item: 0, open: true },
                Entry::Step { item: 0, at: 0 },
                Entry::Step { item: 0, at: 1 },
                Entry::Step { item: 0, at: 2 },
                Entry::End,
                Entry::Step { item: 1, at: 3 },
            ]
        );
        // Each command of an open run reads as itself.
        let second = list[0].one(1);
        assert_eq!((second.start, second.end, second.params[0]["label"].as_i64()), (1, 2, Some(2)));
        assert_eq!(describe(&second), "Add Node");
    }

    /// A click on a run's chip opens it to each command, and a second click folds it.
    #[test]
    fn a_click_on_a_run_opens_and_folds_it() {
        let mut app = ModelApp::new(Session::default());
        for l in 1..=3 {
            app.run("mesh.node.add", json!({"label": l, "xyz": [l, 0, 0]})).unwrap();
        }
        let ctx = egui::Context::default();
        let size = vec2(1400.0, 900.0);
        let draw = |app: &mut ModelApp, events: Vec<egui::Event>| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, size)), events, ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        };
        let click = |app: &mut ModelApp, at: egui::Pos2| {
            let press =
                |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
            draw(app, vec![egui::Event::PointerMoved(at)]);
            draw(app, vec![press(true)]);
            draw(app, vec![press(false)]);
            draw(app, Vec::new());
        };
        for _ in 0..3 {
            draw(&mut app, Vec::new());
        }
        let chip = solidcraft_ui_registry::handle_at("timeline:run:0").expect("the run's chip");
        click(&mut app, chip);
        assert_eq!(app.timeline_open, BTreeSet::from([0]));
        assert_eq!(solidcraft_ui_registry::handles("timeline:step:Add Node").len(), 3, "each command shown");
        click(&mut app, chip);
        assert!(app.timeline_open.is_empty());
        assert_eq!(app.session.journal.len(), 3, "opening a run rolls nothing back");
    }

    /// History is an inset card over the viewport, a shell gap from the window's edges, as
    /// SolidCraft's; the view runs under it, and its navigation bar sits above it.
    #[test]
    fn history_is_a_card_over_the_viewport() {
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let ctx = egui::Context::default();
        let size = vec2(1400.0, 900.0);
        for _ in 0..3 {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, size)), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        }
        let card = app.shell.history.expect("the History card");
        let gap = kit::SHELL_GAP;
        assert_eq!(app.layout.get("timeline"), Some(&card));
        assert!((card.height() - HEIGHT).abs() < 0.5, "{card:?}");
        assert!((card.left() - gap).abs() < 0.5 && (card.right() - (size.x - gap)).abs() < 0.5, "{card:?}");
        let view = app.layout.get("viewport").copied().unwrap_or(Rect::NOTHING);
        assert!(view.contains_rect(card), "the view runs under the card: {view:?} {card:?}");
        assert!((card.bottom() - (view.bottom() - gap)).abs() < 0.5, "{card:?} {view:?}");
        let content = crate::viewport::beside_card(&app, view);
        assert!(content.bottom() <= card.top() - gap + 0.5, "{content:?}");
        assert!(crate::viewport::nav_bar_rect(content).bottom() < card.top(), "the navigation bar sits above History");
        // The card's area is the card, frame after frame: the view above it keeps the pointer.
        let card_layer = egui::LayerId::new(egui::Order::Background, egui::Id::unique("fc_timeline_card"));
        assert_eq!(ctx.layer_id_at(card.center()), Some(card_layer));
        assert_ne!(ctx.layer_id_at(card.center_top() - vec2(0.0, 20.0)), Some(card_layer));
    }
}
