//! The timeline (bottom): the command journal in order, as SolidCraft shows its features. Runs of
//! the same command make one item with a count; the marker stands after the current state, and
//! undone commands follow it greyed, ready to redo. Clicking an item (or the transport buttons)
//! undoes or redoes up to it with `edit.undo` / `edit.redo`, as dragging the marker does; commands
//! older than the undo history cannot be reached and say so.
//!
//! The bar is the kit's History bar (`solidcraft_ui_kit::history_bar`), as SolidCraft draws its
//! own; FieldCraft keeps its journal and undo behind it.

use egui::{Align2, Color32, FontId, Rect, Stroke, pos2, vec2};
use serde_json::{Value, json};
use solidcraft_ui_kit::{self as kit, history_bar};

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
    /// The first command's parameters (the tooltip).
    pub params: Value,
}

impl Item {
    pub fn count(&self) -> usize {
        self.end - self.start
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
            Some(last) if last.command == *command && last.end == i && i != done => last.end = i + 1,
            _ => out.push(Item { command: command.clone(), start: i, end: i + 1, params: params.clone() }),
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
    let name = item.params.get("name").and_then(Value::as_str);
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

/// The journal's items as the kit's History strip lays them out: one step each, the badge of a
/// deletion or an edit and a run's count drawn over its chip.
struct Steps<'a> {
    list: &'a [Item],
    done: usize,
    oldest: usize,
    t: &'a Tokens,
}

impl history_bar::Contents for Steps<'_> {
    fn entries(&mut self, strip: &mut history_bar::Strip<'_, '_>) {
        for (slot, item) in self.list.iter().enumerate() {
            let reachable = item.end >= self.oldest;
            let label = describe(item);
            let tip = || tip(item, reachable);
            let step = history_bar::Step {
                id: slot as u64,
                slot,
                label: &label,
                icon: fieldcraft_engine::catalogue::entry_of(&item.command).map_or("ui.views", |c| c.id),
                icon_accent: None,
                state: history_bar::StepState { rolled_back: item.start >= self.done, suppressed: !reachable, ..Default::default() },
                error: false,
                warning: false,
                ink: None,
                tip: &tip,
            };
            let chip = strip.step(&step);
            if !strip.ui().is_rect_visible(chip) {
                continue;
            }
            let p = strip.painter().clone();
            if let Some(b) = badge(&item.command) {
                paint_badge(&p, chip, b, self.t);
            }
            if item.count() > 1 {
                let text = p.layout_no_wrap(format!("{}", item.count()), FontId::proportional(9.0), self.t.text);
                let b = Rect::from_min_size(chip.right_bottom() - text.size() - vec2(3.0, 1.0), text.size() + vec2(3.0, 1.0));
                p.rect_filled(b, 3.0, self.t.panel_header);
                p.galley(b.min + vec2(1.5, 0.5), text, self.t.text);
            }
        }
    }
}

/// The History bar (bottom): the kit's header (grip, title and the History controls) and strip
/// (a chip per item, the marker after the current state). Its controls, a click on an item and a
/// dragged marker undo or redo through [`ModelApp::roll_to`]; the journal cannot be reordered.
pub fn timeline(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let mut roll: Option<usize> = None;
    egui::Panel::bottom("fc_timeline").exact_size(HEIGHT).show_separator_line(false).frame(egui::Frame::NONE.fill(t.timeline)).show(ui, |ui| {
        let r = ui.max_rect();
        app.place("timeline", r);
        ui.painter().line_segment([r.left_top(), r.right_top()], Stroke::new(1.0, t.border));
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
            marker: list.iter().position(|i| i.start >= done).unwrap_or(list.len()),
            steps: list.len(),
            floating: false,
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
            let mut steps = Steps { list: &list, done, oldest, t: &t };
            actions.extend(history_bar::strip(ui, &bar, r, strip, "fc_timeline_strip", &mut steps).actions);
        }
        let start = |slot: usize| list.get(slot).map_or(done + redo, |i| i.start);
        for action in actions {
            roll = match action {
                history_bar::Action::Seek(history_bar::Seek::First) => Some(oldest),
                history_bar::Action::Seek(history_bar::Seek::Back) => Some(done.saturating_sub(1)),
                history_bar::Action::Seek(history_bar::Seek::Forward) => Some(done + 1),
                history_bar::Action::Seek(history_bar::Seek::Last) => Some(done + redo),
                history_bar::Action::Seek(history_bar::Seek::Slot { slot, .. }) => Some(start(slot)),
                history_bar::Action::Select { step, .. } => list.get(step as usize).filter(|i| i.end >= oldest).map(|i| i.end),
                _ => roll,
            };
        }
    });
    if let Some(target) = roll {
        app.roll_to(target);
    }
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
        let item = |command: &str, params: Value| Item { command: command.into(), start: 0, end: 1, params };
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
}
