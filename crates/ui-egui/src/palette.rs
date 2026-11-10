//! The command line (bottom) and the command search (S or Cmd/Ctrl+K).
//!
//! The command line runs `id {json}` (`mesh.node.add {"label": 1, "xyz": [0, 0, 0]}`), starts a
//! command given by its id alone, or starts the best match for a name. Beside it, the status and
//! the model's size. The search lists every catalogue entry and related command by name, id or
//! tip; coming-soon entries are listed greyed with their milestone, and commands missing what
//! they need greyed with what is missing. Adapted from SolidCraft's
//! `crates/ui-egui/src/palette.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.

use egui::{FontId, RichText};
use fieldcraft_engine::catalogue::{self, CATALOGUE};
use serde_json::Value;

use crate::ModelApp;

use crate::theme::Tokens;
use solidcraft_ui_kit as kit;

/// The command line's height.
pub const HEIGHT: f32 = 28.0;
/// Most results the search lists.
const MAX_RESULTS: usize = 60;

/// One search result.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// What starting it starts: a catalogue entry or a command.
    pub id: &'static str,
    pub label: &'static str,
    /// Where it is on the toolbar: `Tab › Panel`.
    pub place: String,
    /// The icon (its catalogue entry's).
    pub icon: &'static str,
    /// The milestone that brings it, while it is coming soon.
    pub soon: Option<&'static str>,
}

/// Every entry and related command matching `query` (by label, id, panel or tip), best first:
/// labels starting with it, then ids starting with it, then the rest, ready before coming soon.
pub fn matches(query: &str) -> Vec<Hit> {
    let q = query.trim().to_lowercase();
    let mut hits: Vec<(u8, Hit)> = Vec::new();
    let rank = |label: &str, id: &str| -> u8 {
        if label.to_lowercase().starts_with(&q) {
            0
        } else if id.starts_with(&q) {
            1
        } else {
            2
        }
    };
    for c in CATALOGUE {
        let place = format!("{} › {}", c.tab.label(), c.panel);
        let text = format!("{} {} {} {}", c.label, c.id, c.panel, c.tip).to_lowercase();
        if q.is_empty() || text.contains(&q) {
            let soon = c.milestone();
            hits.push((
                rank(c.label, c.id) + if soon.is_some() { 3 } else { 0 },
                Hit { id: c.id, label: c.label, place: place.clone(), icon: c.id, soon },
            ));
        }
        for m in c.more.iter().chain(c.covers) {
            let label = crate::command_label(m);
            if q.is_empty() || format!("{label} {m}").to_lowercase().contains(&q) {
                hits.push((rank(label, m), Hit { id: m, label, place: place.clone(), icon: c.id, soon: None }));
            }
        }
    }
    hits.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.label.cmp(b.1.label)));
    hits.dedup_by(|a, b| a.1.id == b.1.id);
    hits.into_iter().map(|(_, h)| h).collect()
}

/// Run a line of command-line input: `id {json}` runs the command with those parameters, `id`
/// starts it, anything else starts the best match.
pub fn submit(app: &mut ModelApp, line: &str) {
    app.with_origin(fieldcraft_engine::usage::Origin::Palette, |app| submit_line(app, line));
}

/// A line split into its first word and the rest (`id {json}`).
fn split_line(line: &str) -> (&str, &str) {
    let line = line.trim();
    match line.find(|c: char| c.is_whitespace() || c == '{') {
        Some(i) => (line.get(..i).unwrap_or(line), line.get(i..).unwrap_or("").trim()),
        None => (line, ""),
    }
}

fn is_command(id: &str) -> bool {
    fieldcraft_engine::find_command(id).is_some() || catalogue::find(id).is_some()
}

/// The command a line typed into the search names with its parameters (`id {json}`), which the
/// search offers to run as typed instead of searching for it.
pub fn typed_command(line: &str) -> Option<&str> {
    let (id, rest) = split_line(line);
    (is_command(id) && !rest.is_empty()).then_some(id)
}

fn submit_line(app: &mut ModelApp, line: &str) {
    let line = line.trim();
    if line.is_empty() {
        return;
    }
    let (id, rest) = split_line(line);
    if is_command(id) {
        if rest.is_empty() {
            let _ = app.start(id);
        } else {
            match serde_json::from_str::<Value>(rest) {
                Ok(p) => {
                    let _ = app.run(id, p);
                }
                Err(e) => app.set_status(format!("the parameters must be JSON: {e}"), true),
            }
        }
        return;
    }
    match matches(line).into_iter().find(|h| h.soon.is_none()) {
        Some(h) => {
            let _ = app.start(h.id);
        }
        None => app.set_status(format!("no command matches `{line}`"), true),
    }
}

/// The command line's width: 40 % of the bar, within these.
const MIN_LINE_W: f32 = 150.0;
const MAX_LINE_W: f32 = 440.0;

/// The command line's placeholders, longest first: the longest that fits whole is shown.
const PLACEHOLDERS: [&str; 5] = [
    r#"Type a command, e.g. mesh.node.add {"label": 1, "xyz": [0, 0, 0]} — or press S to search"#,
    "Type a command, e.g. mesh.tet.generate — or press S to search",
    "Type a command — or press S to search",
    "Command, or S to search",
    "Command",
];

/// The longest [`PLACEHOLDERS`] that fits a command line `width` wide, so it is never cut.
fn placeholder(ui: &egui::Ui, width: f32) -> &'static str {
    let font = egui::TextStyle::Body.resolve(ui.style());
    // A text field's default margin is 4 points a side.
    let room = width - 8.0;
    let shortest = PLACEHOLDERS[PLACEHOLDERS.len() - 1];
    PLACEHOLDERS.iter().copied().find(|p| crate::elide::width(ui, p, &font) <= room).unwrap_or(shortest)
}

/// The command line, the status and the model's size, above the timeline.
pub fn command_line(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    egui::Panel::bottom("fc_command_line")
        .exact_size(HEIGHT)
        .frame(egui::Frame::NONE.fill(t.panel).inner_margin(egui::Margin::symmetric(8, 3)))
        .show(ui, |ui| {
            app.place("command_line", ui.max_rect());
            ui.horizontal_centered(|ui| {
                ui.label(RichText::new("›").strong().color(t.accent_ink));
                // The rest of the bar keeps room for the status and the model's size.
                let width = (ui.available_width() * 0.4).clamp(MIN_LINE_W, MAX_LINE_W);
                let hint = if app.command_line.is_empty() { placeholder(ui, width) } else { "" };
                let r = ui.add(
                    egui::TextEdit::singleline(&mut app.command_line)
                        .id(egui::Id::unique("fc_command_line_text"))
                        .hint_text(hint)
                        .desired_width(width),
                );
                if r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    let line = std::mem::take(&mut app.command_line);
                    submit(app, &line);
                }
                let (nodes, elements) = (app.session.model.nodes().len(), app.session.model.elements().len());
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(RichText::new(app.gpu_info.as_deref().unwrap_or("CPU renderer")).size(11.5).color(t.text_dim));
                    ui.separator();
                    ui.label(RichText::new(if app.ui.perspective { "Perspective" } else { "Orthographic" }).size(11.5).color(t.text_dim));
                    ui.separator();
                    ui.label(RichText::new(format!("{nodes} nodes  {elements} elements")).size(11.5).color(t.text_dim));
                    ui.separator();
                    crate::updates::notice(app, ui);
                    let (text, color) = match &app.status {
                        Some(s) => (s.text.as_str(), if s.error { t.error } else { t.text }),
                        None => ("Ready", t.text_dim),
                    };
                    // Elided in the middle (both the subject and the outcome stay), whole on hover.
                    let font = FontId::proportional(11.5);
                    let shown = crate::elide::middle_in(ui, text, &font, ui.available_width());
                    let r = ui.add(egui::Label::new(RichText::new(shown.as_ref()).font(font).color(color)));
                    if shown != text {
                        r.on_hover_text(text);
                    }
                });
            });
        });
}

/// The command search, centred near the top of the window.
pub fn popup(app: &mut ModelApp, ctx: &egui::Context) {
    if !app.ui.palette {
        return;
    }
    let mut picked: Option<&'static str> = None;
    let mut run_typed = false;
    let id = egui::Id::unique("fc_search");
    let closed = crate::frame::window(id, "ui.search", "Search commands", crate::frame::Width::Normal)
        .anchor(egui::Align2::CENTER_TOP, egui::vec2(0.0, 140.0))
        .show_body(ctx, |ui| {
            let r = kit::search_field(ui, egui::Id::unique("fc_search_text"), &mut app.palette_text, "command name, id or what it does");
            let entered = r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            r.request_focus();
            // Searched again only when the query changes.
            if app.search.as_ref().is_none_or(|(q, _)| *q != app.palette_text) {
                app.search = Some((app.palette_text.clone(), std::sync::Arc::new(matches(&app.palette_text))));
            }
            let list = app.search.as_ref().map(|(_, l)| l.clone()).unwrap_or_default();
            // `id {json}` typed here (say, typing began with no field focused) runs as typed.
            let typed = typed_command(&app.palette_text).is_some();
            if entered {
                if typed {
                    run_typed = true;
                } else {
                    picked = list.iter().find(|h| h.soon.is_none() && app.unavailable(h.id).is_none()).map(|h| h.id);
                }
            }
            egui::ScrollArea::vertical().id_salt("fc_search_results").max_height(320.0).auto_shrink([false, true]).show(ui, |ui| {
                ui.spacing_mut().item_spacing.y = 0.0;
                if typed {
                    let text = format!("Run {}", app.palette_text.trim());
                    let mut row = kit::MenuItem::new("search:run_typed", &text, "");
                    row.shortcut = "Enter";
                    row.highlighted = true;
                    if row.show(ui).on_hover_text("Run the command with these parameters").clicked() {
                        run_typed = true;
                    }
                }
                for h in list.iter().take(MAX_RESULTS) {
                    // Coming soon, or missing what its command needs: listed, greyed, with why.
                    let missing = app.unavailable(h.id);
                    let enabled = h.soon.is_none() && missing.is_none();
                    let right = match (h.soon, &missing) {
                        (Some(m), _) => format!("{} · soon ({m})", h.place),
                        (None, Some(u)) => u.missing.first().map_or_else(|| h.place.clone(), |m| format!("needs {}", m.0.what())),
                        (None, None) => h.place.clone(),
                    };
                    let mut row = kit::MenuItem::new(h.id, h.label, h.icon);
                    row.shortcut = &right;
                    row.enabled = enabled;
                    let tip = fieldcraft_engine::requirement::Unavailable::tip(h.id, missing.as_ref());
                    if row.show(ui).on_hover_text(tip).clicked() {
                        picked = Some(h.id);
                    }
                }
            });
        });
    if let Some(r) = crate::frame::drawn(ctx, id) {
        app.place("palette", r);
    }
    if run_typed {
        let line = std::mem::take(&mut app.palette_text);
        app.ui.palette = false;
        submit(app, &line);
    } else if let Some(id) = picked {
        app.palette_text.clear();
        app.ui.palette = false;
        let _ = app.with_origin(fieldcraft_engine::usage::Origin::Palette, |app| app.start(id));
    } else if closed || ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
        app.ui.palette = false;
        app.palette_text.clear();
    }
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;

    use super::*;

    #[test]
    fn search_finds_entries_and_related_commands_ready_first() {
        let hits = matches("material");
        assert_eq!(hits.first().map(|h| h.id), Some("model.material.create"));
        assert!(hits.iter().any(|h| h.id == "model.material.delete"), "related commands are found");
        let soon = hits.iter().position(|h| h.soon.is_some()).unwrap();
        assert!(hits.iter().skip(soon).all(|h| h.soon.is_some()), "coming soon last");
        assert_eq!(matches("").len(), {
            let mut ids: Vec<&str> = CATALOGUE.iter().flat_map(|c| std::iter::once(c.id).chain(c.more.iter().chain(c.covers).copied())).collect();
            ids.sort();
            ids.dedup();
            ids.len()
        });
        assert!(matches("no such thing at all").is_empty());
    }

    #[test]
    fn the_command_line_runs_ids_with_json_and_starts_names() {
        let mut app = ModelApp::new(Session::default());
        submit(&mut app, r#"mesh.node.add {"label": 4, "xyz": [1, 2, 3]}"#);
        assert!(app.session.model.node(4).is_some());
        submit(&mut app, "mesh.node.add");
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("mesh.node.add"));
        app.dialog = None;
        submit(&mut app, "material");
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.material.create"));
        submit(&mut app, "mesh.node.add {not json");
        assert!(app.status.as_ref().is_some_and(|s| s.error));
        submit(&mut app, "zzzz qqqq");
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.contains("no command")));
    }

    #[test]
    fn the_search_offers_to_run_a_typed_command() {
        assert_eq!(typed_command(r#"results.contour {"field":"von_mises"}"#), Some("results.contour"));
        assert_eq!(typed_command(r#"  mesh.node.add{"label": 1}"#), Some("mesh.node.add"));
        assert_eq!(typed_command("mesh.node.add"), None, "an id alone is searched for");
        assert_eq!(typed_command("material steel"), None);
    }
}
