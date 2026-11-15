//! The timeline (bottom): the command journal in order, as SolidCraft shows its features. Runs of
//! the same command fold into one chip with their count, which a click opens to show each command
//! (as SolidCraft folds a long run of one component's features); the marker stands after the
//! current state, and undone commands follow it greyed, ready to redo. A run the marker stands
//! inside is drawn open, so the marker shows between its commands. Clicking a command (or the
//! transport buttons) undoes or redoes up to it with `edit.undo` / `edit.redo`, as dragging the
//! marker does; commands older than the undo history cannot be reached and say so.
//!
//! The strip is the kit's History strip (`solidcraft_ui_kit::history_bar`) and its run chip the
//! kit's `HistoryGroupChip`, as SolidCraft draws its own; FieldCraft keeps its journal and undo
//! behind them. Its header is a row of the kit's panel header parts ([`header`]).

use std::collections::BTreeSet;

use egui::{Align2, Color32, FontId, Rect, Stroke, pos2, vec2};
use serde_json::{Value, json};
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _, history_bar};

use crate::ModelApp;
use crate::shell::HistoryDock;
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

    /// Whether the model stands inside this run, `done` commands done: some of its commands are
    /// done and the rest undone.
    fn holds_marker(&self, done: usize) -> bool {
        self.start < done && done < self.end
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

/// The journal and the undone commands as runs of the same command. A run spans the marker when
/// the model stands inside it (the strip draws the marker between its commands), so rolling
/// never regroups the journal: an open run stays open, its commands where they were.
pub fn items(app: &ModelApp) -> (Vec<Item>, usize) {
    let done = app.session.journal.len();
    let all = app.session.journal.iter().chain(app.session.redo_entries());
    let mut out: Vec<Item> = Vec::new();
    for (i, (command, params)) in all.enumerate() {
        match out.last_mut() {
            Some(last) if last.command == *command && last.end == i => {
                last.end = i + 1;
                last.params.push(params.clone());
            }
            _ => out.push(Item { command: command.clone(), start: i, end: i + 1, params: vec![params.clone()] }),
        }
    }
    (out, done)
}

/// Why a command older than the undo history cannot be reached: opening a model starts the
/// history empty (its journal came with the file), and the oldest steps go past the undo memory.
const BEYOND_HISTORY: &str = "older than the undo history (it came with the opened model, or is past the undo memory), so it cannot be undone";

/// Where the History bar asks the model to go.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Roll {
    /// Until the journal holds this many commands (a click on a command, the dragged marker).
    To(usize),
    /// Undo this many commands, as many as the history holds (Back, First).
    Undo(usize),
    /// Redo this many commands, as many as the history holds (Forward, Last).
    Redo(usize),
}

impl ModelApp {
    /// Undo or redo until the journal holds `target` commands: one `edit.undo` or `edit.redo`
    /// with a `count`. A target the history does not reach rolls nothing and says why in the
    /// status bar.
    pub fn roll_to(&mut self, target: usize) {
        let (undo, redo) = self.session.undo_depth();
        let done = self.session.journal.len();
        if target < done.saturating_sub(undo) {
            self.set_status(format!("That command is {BEYOND_HISTORY}"), false);
        } else if target > done + redo {
            self.set_status("Nothing to redo", false);
        } else if target != done {
            let step = if target < done { "edit.undo" } else { "edit.redo" };
            // An error shows in the status bar.
            let _ = self.run(step, json!({"count": done.abs_diff(target)}));
        }
    }

    /// Undo or redo as `roll` asks; with nothing to undo or redo, says why instead.
    fn roll(&mut self, roll: Roll) {
        let (undo, redo) = self.session.undo_depth();
        let done = self.session.journal.len();
        match roll {
            Roll::To(target) => self.roll_to(target),
            Roll::Undo(_) if undo == 0 && done > 0 => self.set_status(format!("Nothing to undo: every command shown is {BEYOND_HISTORY}"), false),
            Roll::Undo(_) if undo == 0 => self.set_status("Nothing to undo", false),
            Roll::Redo(_) if redo == 0 => self.set_status("Nothing to redo", false),
            Roll::Undo(n) => self.roll_to(done.saturating_sub(n.min(undo))),
            Roll::Redo(n) => self.roll_to(done + n.min(redo)),
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
    if reachable {
        s.push_str("\nClick to roll the model to here");
    } else {
        s.push_str(&format!("\nThis command is {BEYOND_HISTORY}"));
    }
    s
}

/// Where the History strip publishes its marker and steps.
const HANDLES: history_bar::Handles<'static> = history_bar::Handles { marker: "timeline:marker", step: "timeline:step:" };
/// Where the header publishes its grip.
const GRIP: &str = "timeline:grip";
/// Where the header publishes each History control (this followed by the control's id).
const CONTROLS: &str = "timeline:";

/// The header drawn: where the strip may start, the grip's response, the History control clicked.
struct Header {
    x: f32,
    grip: egui::Response,
    clicked: Option<kit::HistoryTransport>,
}

/// The header along the left of `rect`, made of the kit's panel header parts: its grip
/// (`panel_grip`), the strip heading when it leaves room for the controls, and the four History
/// controls (no Play: the journal is not played back). The grip and the heading together are the
/// handle that drags the card to the other dock ([`follow_grip`]).
fn header(ui: &mut egui::Ui, rect: Rect) -> Header {
    let t = solidcraft_theme::Tokens::get();
    let inner = Rect::from_min_max(pos2(rect.left() + t.space_10, rect.top()), rect.right_bottom());
    let mut row = ui.new_child(egui::UiBuilder::new().max_rect(inner).layout(egui::Layout::left_to_right(egui::Align::Center)));
    row.spacing_mut().item_spacing.x = 0.0;
    let mut handle = kit::panel_grip(&mut row).rect;
    row.add_space(t.space_8);
    let heading = kit::StripHeading::new(row.painter(), "History");
    if heading.width() + kit::history_transport_width() <= row.available_width() {
        handle = handle.union(heading.show(&mut row).rect);
    }
    let grip = ui.interact(handle, ui.scope_id().with(GRIP), egui::Sense::drag());
    solidcraft_ui_registry::publish_handle(GRIP, handle.center());
    let mut clicked = None;
    for (k, control) in kit::HistoryTransport::ALL.into_iter().enumerate() {
        if k > 0 {
            row.add_space(t.space_2);
        }
        let (_, bounds) = row.allocate_space(vec2(t.space_24, t.space_24));
        // Play's labels go unused: it is not among the four.
        if control.show_as(&mut row, bounds, false, CONTROLS, kit::PlayLabels { play: "", stop: "" }).clicked() {
            clicked = Some(control);
        }
    }
    Header { x: row.min_rect().right() + t.space_16, grip, clicked }
}

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
            handle: None,
        };
        let chip = strip.step(&step);
        if strip.ui().is_rect_visible(chip)
            && let Some(b) = badge(&item.command)
        {
            paint_badge(strip.painter(), chip, b, self.t);
        }
    }

    /// A run's chip: the kit's group chip (its icon, name and count); a click opens or folds it,
    /// unless the marker stands inside it, which keeps it open.
    fn run(&mut self, strip: &mut history_bar::Strip<'_, '_>, item: &Item, open: bool) {
        let (done, oldest, t, accent) = (self.done, self.oldest, self.t, self.accent);
        let held = item.holds_marker(done);
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
            let what = match (held, open) {
                (true, _) => "shown one by one while the model stands among them",
                (false, true) => "click to fold them",
                (false, false) => "click to show each one",
            };
            let mut text = format!("{}: {} commands, {what}", describe(item), item.count());
            if item.end < oldest {
                text.push_str(&format!("\nThese commands are {BEYOND_HISTORY}"));
            }
            clicked = response.on_hover_cursor(egui::CursorIcon::PointingHand).backing_on_hover_text(text).clicked();
            rect.width()
        });
        if clicked && !held {
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

/// The History bar, as SolidCraft's: in the floating layout a kit Panel card over the viewport
/// at its dock, along the bottom (`Container::floating_card(..).bottom_strip_rect`) or a shell gap
/// under the top cards (`app.ui.history_dock`, saved in the preferences), recorded in
/// `app.shell.history`; in the compact layout a strip docked along the window's bottom edge. The
/// floating card's grip and title drag it to the other dock ([`follow_grip`]).
pub fn timeline(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let host = kit::Container::new(ui.ctx(), kit::Surface::Panel);
    let (roll, toggle) = if host.floating() {
        let available = ui.available_rect_before_wrap();
        let bottom = floating_rect(available);
        // A window too short for the card leaves History out rather than drawing into a
        // negative-size rectangle.
        if !bottom.is_finite() || bottom.width() < 1.0 || bottom.height() < HEIGHT - 1.0 {
            return;
        }
        let top = top_dock(bottom, &app.shell.cards, available);
        let rect = match app.ui.history_dock {
            HistoryDock::Top => top,
            HistoryDock::Bottom => bottom,
        };
        app.shell.history = Some(rect);
        app.place("timeline", rect);
        let frame = host.frame(ui.ctx(), t.timeline, kit::frame::Side::Left);
        let (roll, toggle, dock) = host
            .area(egui::Id::unique("fc_timeline_card"))
            .fixed_pos(rect.min)
            .show(ui.ctx(), |ui| {
                // The card takes its size in the area (and keeps the pointer off the view under it).
                let (card, _) = ui.allocate_exact_size(rect.size(), egui::Sense::hover());
                kit::Container::paint_panel_frame(ui.painter(), frame, card);
                let mut child = ui.new_child(egui::UiBuilder::new().id_salt("fc_timeline").max_rect(card));
                child.set_clip_rect(card);
                let (roll, toggle, grip) = contents(app, &mut child, &t, true);
                (roll, toggle, follow_grip(ui, &grip, rect, top, bottom))
            })
            .inner;
        if let Some(dock) = dock {
            app.ui.history_dock = dock;
        }
        (roll, toggle)
    } else {
        let panel = host.edge("fc_timeline", kit::Edge::Bottom).exact_size(HEIGHT).show_separator_line(false);
        let shown = host.show_edge(ui, panel, host.edge_frame(t.timeline), |ui| {
            let r = ui.max_rect();
            app.place("timeline", r);
            ui.painter().line_segment([r.left_top(), r.right_top()], Stroke::new(1.0, t.border));
            let (roll, toggle, _) = contents(app, ui, &t, false);
            (roll, toggle)
        });
        shown.inner
    };
    if let Some(start) = toggle
        && !app.timeline_open.remove(&start)
    {
        app.timeline_open.insert(start);
    }
    if let Some(roll) = roll {
        app.roll(roll);
    }
}

/// Where the floating History card goes in `available` (the viewport below the application
/// bar): the kit's shared bottom strip, a shell gap from the edges.
pub(crate) fn floating_rect(available: Rect) -> Rect {
    kit::Container::floating_card(kit::Surface::Panel).bottom_strip_rect(available, HEIGHT)
}

/// The top dock: the card at the `bottom` dock moved up to a shell gap under the top `cards` (the
/// top of `available` without any), never below the bottom dock in a short window.
fn top_dock(bottom: Rect, cards: &[Rect], available: Rect) -> Rect {
    let under_cards = cards.iter().fold(available.top(), |top, card| top.max(card.bottom())) + kit::SHELL_GAP;
    Rect::from_min_size(pos2(bottom.left(), under_cards.min(bottom.top())), bottom.size())
}

/// The grip's drag between the docks `top` and `bottom` of the card drawn `at` one of them: while
/// dragged, the card's graphics follow the pointer up and down between the docks (egui's
/// drag-and-drop preview, `transform_layer_shapes`: nothing is remembered); released, the dock
/// nearer the pointer is returned. Called once the card is drawn.
fn follow_grip(ui: &egui::Ui, grip: &egui::Response, at: Rect, top: Rect, bottom: Rect) -> Option<HistoryDock> {
    if !(grip.dragged() || grip.drag_stopped()) {
        if grip.hovered() {
            ui.ctx().set_cursor_icon(egui::CursorIcon::Grab);
        }
        return None;
    }
    ui.ctx().set_cursor_icon(egui::CursorIcon::Grabbing);
    // The grip sits on the card's middle line, so the pointer's height is the card's middle.
    let y = ui.input(|i| i.pointer.latest_pos()).map_or(at.center().y, |p| p.y);
    if grip.drag_stopped() {
        let nearer_top = (y - top.center().y).abs() < (y - bottom.center().y).abs();
        return Some(if nearer_top { HistoryDock::Top } else { HistoryDock::Bottom });
    }
    let dy = kit::within(y - at.center().y, top.top() - at.top(), bottom.top() - at.top());
    ui.ctx().transform_layer_shapes(ui.layer_id(), egui::emath::TSTransform::from_translation(vec2(0.0, dy)));
    None
}

/// The [`header`] (grip, title and the History controls) and the kit's strip (a chip per command
/// or run, the marker after the current state) in `ui`'s rect. Returns where to roll the model and
/// the run to open or fold, and the grip's response (the grip and the title). Its controls, a click
/// on a command and a dragged marker undo or redo through [`ModelApp::roll`]; the journal cannot be
/// reordered.
fn contents(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, floating: bool) -> (Option<Roll>, Option<usize>, egui::Response) {
    let mut toggle: Option<usize> = None;
    let r = ui.max_rect();
    let (undo, _) = app.session.undo_depth();
    let (list, done) = cached_items(app);
    let oldest = done.saturating_sub(undo);
    let accent = solidcraft_theme::TabAccent::current(ui.ctx());
    let bar = history_bar::Bar {
        marker_tip: "The model is here: drag to roll it back or forward",
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
    let head = header(ui, r);
    let mut roll = head.clicked.and_then(|control| match control {
        // First and Last go as far as the history holds.
        kit::HistoryTransport::First => Some(Roll::Undo(usize::MAX)),
        kit::HistoryTransport::Back => Some(Roll::Undo(1)),
        kit::HistoryTransport::Forward => Some(Roll::Redo(1)),
        kit::HistoryTransport::Last => Some(Roll::Redo(usize::MAX)),
        // Not shown: the journal is not played back.
        kit::HistoryTransport::Play => None,
    });
    let mut actions = Vec::new();
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
        // A run the marker stands inside is drawn open, so the marker shows between its commands
        // and its undone ones dimmed; it is not remembered as opened.
        let open: BTreeSet<usize> = app.timeline_open.iter().copied().chain(list.iter().filter(|i| i.holds_marker(done)).map(|i| i.start)).collect();
        let shown = entries(&list, &open);
        let mut steps = Steps { list: &list, entries: &shown, done, oldest, t, accent, toggle: None };
        actions.extend(history_bar::strip(ui, &bar, r, strip, "fc_timeline_strip", &mut steps).actions);
        toggle = steps.toggle;
    }
    for action in actions {
        roll = match action {
            history_bar::Action::Seek { slot, .. } => Some(Roll::To(slot)),
            // A command's chip rolls the model to just after it (or says why it cannot).
            history_bar::Action::Select { step, .. } => usize::try_from(step).ok().map(|at| Roll::To(at + 1)),
            _ => roll,
        };
    }
    (roll, toggle, head.grip)
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
        // The marker stands inside the run of nodes, which stays one item.
        let (list, done) = items(&app);
        assert_eq!(done, 2);
        assert_eq!(list.iter().map(|i| (i.start, i.end)).collect::<Vec<_>>(), [(0, 3), (3, 4)]);
        app.roll_to(4);
        assert!(app.session.model.materials().contains_key("Steel"));
        // From below the tip, a redo past what was undone goes as far as the history holds.
        app.roll_to(1);
        app.roll(Roll::Redo(100));
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

    /// Draws `app` in a window of `size` with these input events, as a frame of the app does.
    fn frame(ctx: &egui::Context, app: &mut ModelApp, size: egui::Vec2, events: Vec<egui::Event>) {
        let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, size)), events, ..Default::default() };
        let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
        out.textures_delta.clear();
    }

    fn press(at: egui::Pos2, pressed: bool) -> egui::Event {
        egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() }
    }

    /// A click as a hand makes it: move there, press, release, one frame each.
    fn click_at(ctx: &egui::Context, app: &mut ModelApp, size: egui::Vec2, at: egui::Pos2) {
        frame(ctx, app, size, vec![egui::Event::PointerMoved(at)]);
        frame(ctx, app, size, vec![press(at, true)]);
        frame(ctx, app, size, vec![press(at, false)]);
        frame(ctx, app, size, Vec::new());
    }

    /// Four reference points, one run of one command: rolling to any of its commands from the
    /// open run leaves the run open, so each next command is still where the hand reaches for
    /// it. (Shown in the app on 2026-11-15, shots H1-03/H1-04: the run split at the marker into
    /// a new, folded run, so every click after the first found nothing.)
    #[test]
    fn rolling_within_an_open_run_keeps_each_command_clickable() {
        let mut app = ModelApp::new(Session::default());
        for k in 0..4 {
            app.run("model.reference_point.create", json!({"name": format!("RP{k}"), "xyz": [0, 0, k]})).unwrap();
        }
        let ctx = egui::Context::default();
        let size = vec2(1400.0, 900.0);
        for _ in 0..3 {
            frame(&ctx, &mut app, size, Vec::new());
        }
        let run = solidcraft_ui_registry::handle_at("timeline:run:0").expect("the run's chip");
        click_at(&ctx, &mut app, size, run);
        assert_eq!(app.timeline_open, BTreeSet::from([0]));
        for k in [0, 3, 1, 2, 0, 3] {
            let name = format!("timeline:step:Add reference point RP{k}");
            let Some(at) = solidcraft_ui_registry::handle_at(&name) else { panic!("{name} is shown after rolling to {}", app.session.journal.len()) };
            click_at(&ctx, &mut app, size, at);
            assert_eq!(app.session.journal.len(), k + 1, "a click on RP{k} rolls the model to just after it");
            assert_eq!(app.timeline_open, BTreeSet::from([0]), "the run stays open");
        }
    }

    /// A sample opened from the start page comes with its journal but no undo history (opening
    /// starts it empty): no History control and no command older than the history rolls the
    /// model, and each click says why instead of doing nothing. A command made after it rolls.
    /// (Shown in the app on 2026-11-15, W1-22 and F2-01: 13/13 clicks did nothing, with no word.)
    #[test]
    fn clicks_beyond_the_undo_history_roll_nothing_and_say_why() {
        let mut app = ModelApp::new(Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        let done = app.session.journal.len();
        assert!(done > 1, "the sample has a journal");
        assert_eq!(app.session.undo_depth(), (0, 0), "opening starts the undo history empty");
        let ctx = egui::Context::default();
        let size = vec2(1400.0, 900.0);
        for _ in 0..3 {
            frame(&ctx, &mut app, size, Vec::new());
        }
        // A single command older than the history (the sample's are each one of a kind or runs).
        let (list, _) = items(&app);
        let old = list.iter().find(|i| i.count() == 1 && i.end < done).expect("a single command before the last");
        let step = format!("{}{}", HANDLES.step, describe(old));
        let controls = ["first", "back", "forward", "last"].map(|c| format!("{}{c}", CONTROLS));
        for name in controls.iter().chain([&step]) {
            app.status = None;
            let at = solidcraft_ui_registry::handle_at(name).unwrap_or_else(|| panic!("{name} is drawn"));
            click_at(&ctx, &mut app, size, at);
            assert_eq!(app.session.journal.len(), done, "{name} rolls nothing");
            let said = app.status.as_ref().map(|s| s.text.clone()).unwrap_or_default();
            let why = if name.ends_with("forward") || name.ends_with("last") { "Nothing to redo" } else { "undo history" };
            assert!(said.contains(why), "{name} says why it did nothing: {said:?}");
            assert!(!said.contains('{'), "in words: {said:?}");
        }
        // A command made now is in the history: Back undoes it and Forward redoes it.
        app.run("model.reference_point.create", json!({"name": "RP1", "xyz": [0, 0, 0]})).unwrap();
        frame(&ctx, &mut app, size, Vec::new());
        let back = solidcraft_ui_registry::handle_at(&controls[1]).expect("Back");
        click_at(&ctx, &mut app, size, back);
        assert_eq!(app.session.journal.len(), done);
        let forward = solidcraft_ui_registry::handle_at(&controls[2]).expect("Forward");
        click_at(&ctx, &mut app, size, forward);
        assert_eq!(app.session.journal.len(), done + 1);
    }

    /// The marker inside a run marks it apart from its done and undone commands, but the run
    /// stays one item: rolling never regroups the journal.
    #[test]
    fn a_run_stays_one_item_across_the_marker() {
        let mut app = ModelApp::new(Session::default());
        for l in 1..=3 {
            app.run("mesh.node.add", json!({"label": l, "xyz": [l, 0, 0]})).unwrap();
        }
        app.roll_to(1);
        let (list, done) = items(&app);
        assert_eq!(done, 1);
        assert_eq!(list.iter().map(|i| (i.start, i.end)).collect::<Vec<_>>(), [(0, 3)]);
    }

    /// The colours `shapes` paint inside `within`.
    fn colours_in(shapes: &[egui::epaint::ClippedShape], within: Rect) -> Vec<Color32> {
        fn walk(shape: &egui::Shape, within: Rect, out: &mut Vec<Color32>) {
            if let egui::Shape::Vec(v) = shape {
                v.iter().for_each(|s| walk(s, within, out));
                return;
            }
            if !within.contains_rect(shape.visual_bounding_rect()) {
                return;
            }
            match shape {
                egui::Shape::Rect(r) => out.extend([r.fill, r.stroke.color]),
                egui::Shape::Circle(c) => out.extend([c.fill, c.stroke.color]),
                egui::Shape::Path(p) => {
                    out.push(p.fill);
                    if let egui::epaint::ColorMode::Solid(c) = p.stroke.color {
                        out.push(c);
                    }
                }
                egui::Shape::LineSegment { stroke, .. } => out.push(stroke.color),
                egui::Shape::Text(t) => {
                    out.push(t.fallback_color);
                    out.extend(t.galley.job.sections.iter().map(|s| s.format.color));
                }
                egui::Shape::Mesh(m) => out.extend(m.vertices.iter().map(|v| v.color)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        shapes.iter().for_each(|c| walk(&c.shape, within, &mut out));
        out
    }

    /// Back from the tip of a folded run of four reference points: the run is drawn open, as a
    /// tray, the marker between its third and fourth commands and the fourth dimmed, undone; at
    /// the tip it is folded again. Nothing is remembered as opened (review of pr/ux-right-cards
    /// 8da16a22 and asks/lead-2: a folded run spanning the marker looked the same after a step
    /// back).
    #[test]
    fn a_folded_run_holding_the_marker_is_drawn_open() {
        let mut app = ModelApp::new(Session::default());
        for k in 0..4 {
            app.run("model.reference_point.create", json!({"name": format!("RP{k}"), "xyz": [0, 0, k]})).unwrap();
        }
        let ctx = egui::Context::default();
        let size = vec2(1400.0, 900.0);
        for _ in 0..3 {
            frame(&ctx, &mut app, size, Vec::new());
        }
        let steps = || (0..4).map(|k| solidcraft_ui_registry::handle_at(&format!("{}Add reference point RP{k}", HANDLES.step))).collect::<Vec<_>>();
        assert!(steps().iter().all(Option::is_none), "at the tip the run is folded");
        let back = solidcraft_ui_registry::handle_at(&format!("{}back", CONTROLS)).expect("Back");
        click_at(&ctx, &mut app, size, back);
        assert_eq!(app.session.journal.len(), 3);
        let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, size)), ..Default::default() };
        let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
        out.textures_delta.clear();
        let Some(at) = steps().into_iter().collect::<Option<Vec<_>>>() else { panic!("each command is drawn: {:?}", steps()) };
        let marker = solidcraft_ui_registry::rect_of(HANDLES.marker).expect("the marker");
        assert!(at[2].x < marker.center().x && marker.center().x < at[3].x, "the marker between RP2 and RP3: {marker:?} {at:?}");
        // A step's square chip around its handle: the strip's pitch less the gap to the next
        // chip's edge, so only that step's own shapes lie inside.
        let pitch = at[1].x - at[0].x;
        let chip = |p: egui::Pos2| Rect::from_center_size(p, egui::Vec2::splat(pitch - 2.0));
        let colours = |k: usize| colours_in(&out.shapes, chip(at[k]));
        assert!(!colours(1).is_empty(), "the steps are drawn");
        assert_eq!(colours(0), colours(1), "the done commands are drawn alike");
        assert_ne!(colours(3), colours(1), "the undone command is dimmed");
        assert!(app.timeline_open.is_empty(), "nothing is remembered as opened");
        let forward = solidcraft_ui_registry::handle_at(&format!("{}forward", CONTROLS)).expect("Forward");
        click_at(&ctx, &mut app, size, forward);
        assert_eq!(app.session.journal.len(), 4);
        assert!(steps().iter().all(Option::is_none), "back at the tip the run is folded");
    }

    /// Drags the grip by `dy` points, a frame per step, and releases it.
    fn drag_grip(ctx: &egui::Context, app: &mut ModelApp, size: egui::Vec2, dy: f32) {
        let grip = solidcraft_ui_registry::handle_at(GRIP).expect("the grip");
        let to = grip + vec2(0.0, dy);
        frame(ctx, app, size, vec![egui::Event::PointerMoved(grip)]);
        frame(ctx, app, size, vec![press(grip, true)]);
        for k in 1..=10 {
            frame(ctx, app, size, vec![egui::Event::PointerMoved(grip + (to - grip) * (k as f32 / 10.0))]);
        }
        frame(ctx, app, size, vec![press(to, false)]);
        frame(ctx, app, size, Vec::new());
    }

    /// The grip drags History to the other dock, as the kit's grip is meant to ("drags the bar
    /// to another dock"): dropped nearer the top it docks under the top cards, nearer the bottom
    /// it stays along the bottom. The dock is the app's (`UiState::history_dock`, saved in the
    /// preferences), so a new egui context, as `ui.capture` runs, draws it there too.
    #[test]
    fn the_grip_docks_history_at_the_top_or_the_bottom() {
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let ctx = egui::Context::default();
        let size = vec2(1400.0, 900.0);
        for _ in 0..3 {
            frame(&ctx, &mut app, size, Vec::new());
        }
        let bottom = app.shell.history.expect("the History card");
        assert_eq!(app.ui.history_dock, HistoryDock::Bottom);
        // The grip heads the card: inside it, on its middle line, before the History controls.
        let grip = solidcraft_ui_registry::handle_at(GRIP).expect("the grip");
        let first = solidcraft_ui_registry::handle_at(&format!("{CONTROLS}first")).expect("First");
        assert!(bottom.contains(grip) && (grip.y - bottom.center().y).abs() < 0.5 && grip.x < first.x, "{grip:?} {first:?} {bottom:?}");
        // Dropped nearer the bottom: it stays there.
        drag_grip(&ctx, &mut app, size, -250.0);
        assert_eq!((app.ui.history_dock, app.shell.history), (HistoryDock::Bottom, Some(bottom)));
        // Dropped nearer the top: it docks a shell gap under the top cards, as wide as before.
        drag_grip(&ctx, &mut app, size, -650.0);
        assert_eq!(app.ui.history_dock, HistoryDock::Top);
        let top = app.shell.history.expect("the History card");
        let below = app.shell.cards.iter().map(|c| c.bottom()).fold(f32::NEG_INFINITY, f32::max) + kit::SHELL_GAP;
        assert!((top.top() - below).abs() < 0.5, "under the top cards: {top:?} {:?}", app.shell.cards);
        assert_eq!((top.left(), top.size()), (bottom.left(), bottom.size()));
        assert_eq!(app.session.journal.len(), 1, "moving the card rolls nothing");
        // A new context (ui.capture's) draws it at the top as well.
        let fresh = egui::Context::default();
        for _ in 0..3 {
            frame(&fresh, &mut app, size, Vec::new());
        }
        assert_eq!(app.shell.history, Some(top));
        // A shorter window keeps it at the top, inside the view.
        let small = vec2(900.0, 500.0);
        for _ in 0..3 {
            frame(&ctx, &mut app, small, Vec::new());
        }
        let card = app.shell.history.expect("the History card");
        let view = app.layout.get("viewport").copied().unwrap_or(Rect::NOTHING);
        assert!(view.shrink(kit::SHELL_GAP - 0.5).contains_rect(card), "{view:?} {card:?}");
        // And back down to the bottom.
        for _ in 0..3 {
            frame(&ctx, &mut app, size, Vec::new());
        }
        drag_grip(&ctx, &mut app, size, 650.0);
        assert_eq!((app.ui.history_dock, app.shell.history), (HistoryDock::Bottom, Some(bottom)));
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
        assert!(crate::navigation::bar_rect(view, content).bottom() < card.top(), "the navigation bar sits above History");
        // The card's area is the card, frame after frame: the view above it keeps the pointer.
        let card_layer = egui::LayerId::new(egui::Order::Background, egui::Id::unique("fc_timeline_card"));
        assert_eq!(ctx.layer_id_at(card.center()), Some(card_layer));
        assert_ne!(ctx.layer_id_at(card.center_top() - vec2(0.0, 20.0)), Some(card_layer));
    }
}
