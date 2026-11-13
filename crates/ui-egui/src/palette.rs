//! The command line (bottom) and the command search (S or Cmd/Ctrl+K).
//!
//! The command line runs `id {json}` (`mesh.node.add {"label": 1, "xyz": [0, 0, 0]}`), starts a
//! command given by its id alone, or starts the best match for a name. Beside it, the status and
//! the model's size. The search lists every catalogue entry and related command by name, id or
//! tip; coming-soon entries are listed greyed with their milestone, and commands missing what
//! they need greyed with what is missing. Adapted from SolidCraft's
//! `crates/ui-egui/src/palette.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.

use egui::{FontId, Rect, RichText};
use fieldcraft_engine::catalogue::{self, CATALOGUE};
use serde_json::Value;
use solidcraft_ui_kit::BackingTooltipExt as _;

use crate::ModelApp;

use crate::theme::Tokens;
use solidcraft_ui_kit as kit;

/// The command line's height: the kit's field (28 points) and a 3 point margin above and below.
pub const HEIGHT: f32 = 34.0;
/// Most results the search lists.
const MAX_RESULTS: usize = 60;
/// The longest query the search field takes, in characters.
const MAX_QUERY: usize = 4096;

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

/// The command line's width: 40 % of the bar while it is in use; idle (empty, not focused) it
/// gives the status message the room it needs to show whole, down to [`MIN_LINE_W`].
fn line_width(app: &ModelApp, ui: &egui::Ui) -> f32 {
    let full = (ui.available_width() * 0.4).clamp(MIN_LINE_W, MAX_LINE_W);
    let focused = ui.memory(|m| m.has_focus(egui::Id::unique("fc_command_line_text")));
    if focused || !app.command_line.is_empty() {
        return full;
    }
    let k = solidcraft_theme::Tokens::get();
    let font = FontId::proportional(k.font_body);
    let (nodes, elements) = (app.session.model.nodes().len(), app.session.model.elements().len());
    let items = [
        app.gpu_info.as_deref().unwrap_or("CPU renderer").to_string(),
        (if app.ui.perspective { "Perspective" } else { "Orthographic" }).to_string(),
        format!("{nodes} nodes  {elements} elements"),
        app.status.as_ref().map_or("Ready", |s| s.text.as_str()).to_string(),
    ];
    let gap = ui.spacing().item_spacing.x;
    // Each item, and a rule with its gaps between them.
    let needed: f32 = items.iter().map(|text| crate::elide::width(ui, text, &font) + 2.0 * gap + 1.0).sum::<f32>() + gap;
    (ui.available_width() - needed).clamp(MIN_LINE_W, full)
}

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
    // The kit's field: monospace text, its inset on each side.
    let t = solidcraft_theme::Tokens::get();
    let font = egui::FontId::monospace(t.font_label);
    let room = width - 2.0 * t.form_inset;
    let shortest = PLACEHOLDERS[PLACEHOLDERS.len() - 1];
    PLACEHOLDERS.iter().copied().find(|p| crate::elide::width(ui, p, &font) <= room).unwrap_or(shortest)
}

/// A rule between two of the status bar's items, as SolidCraft's status bar divides its own.
pub(crate) fn rule(ui: &mut egui::Ui) {
    let k = solidcraft_theme::Tokens::get();
    let (r, _) = ui.allocate_exact_size(egui::vec2(1.0, ui.available_height()), egui::Sense::hover());
    k.secondary_vline(ui.painter(), r.center().x, (r.top() + k.space_6)..=(r.bottom() - k.space_6), k.border);
}

/// The command line, the status and the model's size, above the timeline. The status items
/// follow SolidCraft's status bar (its body text, a hairline along the top, rules between the
/// items); the kit has no status bar of its own, and the command line is FieldCraft's.
pub fn command_line(app: &mut ModelApp, ui: &mut egui::Ui) {
    let t = Tokens::get(app.ui.dark);
    let k = solidcraft_theme::Tokens::get();
    // The status bar is an application-bar surface along the window's bottom edge, as
    // SolidCraft's (docs/ui-container-audit.md).
    let host = kit::Container::new(ui.ctx(), kit::Surface::ApplicationBar);
    let panel = host.edge("fc_command_line", kit::Edge::Bottom).exact_size(HEIGHT);
    let margin = egui::vec2(k.space_8, k.space_3);
    let frame = host.edge_frame(t.panel).inner_margin(egui::Margin::symmetric(margin.x as i8, margin.y as i8));
    host.show_edge(ui, panel, frame, |ui| {
        app.place("command_line", ui.max_rect());
        let bar = ui.max_rect().expand2(margin);
        ui.painter().hline(bar.x_range(), bar.top() + 0.5, k.primary_border(k.border));
        ui.horizontal_centered(|ui| {
            ui.label(RichText::new("›").strong().color(t.accent_ink));
            // The rest of the bar keeps room for the status and the model's size.
            let width = line_width(app, ui);
            let hint = if app.command_line.is_empty() { placeholder(ui, width) } else { "" };
            let r = ui
                .allocate_ui(egui::vec2(width, kit::row_height()), |ui| {
                    kit::cell_editor(ui, egui::Id::unique("fc_command_line_text"), &mut app.command_line, hint, false)
                })
                .inner;
            app.place("command_line:field", r.rect);
            if r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                let line = std::mem::take(&mut app.command_line);
                submit(app, &line);
            }
            let (nodes, elements) = (app.session.model.nodes().len(), app.session.model.elements().len());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let item = |text: &str| RichText::new(text).size(k.font_body).color(t.text_dim);
                ui.label(item(app.gpu_info.as_deref().unwrap_or("CPU renderer")));
                rule(ui);
                ui.label(item(if app.ui.perspective { "Perspective" } else { "Orthographic" }));
                rule(ui);
                ui.label(item(&format!("{nodes} nodes  {elements} elements")));
                rule(ui);
                crate::updates::notice(app, ui);
                let (text, color) = match &app.status {
                    Some(s) => (s.text.as_str(), if s.error { t.error } else { t.text }),
                    None => ("Ready", t.text_dim),
                };
                // Elided in the middle (both the subject and the outcome stay), whole on hover.
                let font = FontId::proportional(k.font_body);
                let shown = crate::elide::middle_in(ui, text, &font, ui.available_width());
                let r = ui.add(egui::Label::new(RichText::new(shown.as_ref()).font(font).color(color)));
                if shown != text {
                    r.backing_on_hover_text(text);
                }
            });
        });
    });
}

/// Where the open search field's top left goes over the title bar's closed `field` in `screen`:
/// its right edge kept, inside the window. The title bar keeps its document tab clear of it.
pub fn open_field_at(screen: Rect, field: Rect) -> egui::Pos2 {
    let k = solidcraft_theme::Tokens::get();
    let width = kit::CommandSearch::width();
    egui::pos2((field.right() - width).clamp(screen.left() + k.space_8, (screen.right() - width - k.space_8).max(screen.left())), field.top())
}

/// The command search: the kit's search field opened where the title bar's closed field is, with
/// the kit's results sheet under it, as SolidCraft's command search opens.
pub fn popup(app: &mut ModelApp, ctx: &egui::Context) {
    if !app.ui.palette {
        return;
    }
    let mut picked: Option<&'static str> = None;
    let mut run_typed = false;
    let k = solidcraft_theme::Tokens::get();
    let width = kit::CommandSearch::width();
    // Opened over the title bar's field, its right edge kept; centred near the top without one.
    let screen = ctx.content_rect();
    let field = app
        .layout
        .get("app_bar:search")
        .copied()
        .unwrap_or_else(|| Rect::from_center_size(egui::pos2(screen.center().x, 60.0), egui::vec2(width, 26.0)));
    let at = open_field_at(screen, field);
    let key = crate::appbar::shortcut("K");
    // An anchored popup host over the title bar's field.
    let host = kit::Container::new(ctx, kit::Surface::Popup);
    let shown = host.area(egui::Id::unique("fc_search")).order(egui::Order::Foreground).fixed_pos(at).show(ctx, |ui| {
        ui.set_width(width);
        ui.spacing_mut().item_spacing.y = 0.0;
        let geometry = kit::SearchFieldGeometry { height: field.height(), radius: kit::CommandSearch::field_radius(1.0, false), progress: 1.0 };
        let r = kit::CommandSearch::field(ui, &mut app.palette_text, false, true, true, geometry, Some(&key), MAX_QUERY);
        let entered = r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
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
        ui.add_space(k.space_4);
        let query = app.palette_text.clone();
        kit::CommandSearch::sheet(ui, width, 1.0, typed, |ui| {
            egui::ScrollArea::vertical().id_salt("fc_search_results").max_height(360.0).auto_shrink([false, true]).show(ui, |ui| {
                if typed {
                    let text = format!("Run {}", query.trim());
                    let mut row = kit::MenuItem::new("search:run_typed", &text, "");
                    row.shortcut = "Enter";
                    row.highlighted = true;
                    if row.search_row(ui, "", "Runs the command with these parameters", None).clicked() {
                        run_typed = true;
                    }
                }
                kit::CommandSearch::heading(ui, "search.commands", list.len().min(MAX_RESULTS));
                // Enter starts the first command that can run: lit as the keyboard's row.
                let first = list.iter().take(MAX_RESULTS).position(|h| h.soon.is_none() && app.unavailable(h.id).is_none());
                for (i, h) in list.iter().take(MAX_RESULTS).enumerate() {
                    // Coming soon, or missing what its command needs: listed, greyed, with why.
                    let missing = app.unavailable(h.id);
                    let enabled = h.soon.is_none() && missing.is_none();
                    let path = match (h.soon, &missing) {
                        (Some(m), _) => format!("{} · soon ({m})", h.place),
                        (None, Some(u)) => u.missing.first().map_or_else(|| h.place.clone(), |m| format!("needs {}", m.0.what())),
                        (None, None) => h.place.clone(),
                    };
                    let mut row = kit::MenuItem::new(h.id, h.label, h.icon);
                    row.enabled = enabled;
                    row.highlighted = !typed && first == Some(i);
                    let tip = fieldcraft_engine::requirement::Unavailable::tip(h.id, missing.as_ref());
                    if row.search_row(ui, &query, &path, None).backing_on_hover_text(tip).clicked() {
                        picked = Some(h.id);
                    }
                }
                if list.is_empty() && !typed {
                    kit::status_line(ui, &solidcraft_ui_text::format("search.no_match", &[("query", query.trim())]), kit::StatusTone::Dim);
                }
            });
        });
    });
    app.place("palette", shown.response.rect);
    let clicked_outside =
        ctx.input(|i| i.pointer.any_click() && i.pointer.interact_pos().is_some_and(|p| !shown.response.rect.contains(p) && !field.contains(p)));
    if run_typed {
        let line = std::mem::take(&mut app.palette_text);
        app.ui.palette = false;
        submit(app, &line);
    } else if let Some(id) = picked {
        app.palette_text.clear();
        app.ui.palette = false;
        let _ = app.with_origin(fieldcraft_engine::usage::Origin::Palette, |app| app.start(id));
    } else if clicked_outside || ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
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

    /// Idle (empty, not focused) the command line gives a long status message its room; typed
    /// in, it takes its usual share of the bar (A3-17).
    #[test]
    fn an_idle_command_line_yields_to_the_status() {
        let field = |typed: &str| {
            let mut app = ModelApp::new(fieldcraft_engine::Session::default());
            app.set_status("Model Info: 315 nodes, 1200 elements, 2 steps and 4 materials (details in the panel)", false);
            app.command_line = typed.to_string();
            let ctx = egui::Context::default();
            for _ in 0..2 {
                let raw =
                    egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1024.0, 700.0))), ..Default::default() };
                let mut out = ctx.run_ui(raw, |ui| command_line(&mut app, ui));
                out.textures_delta.clear();
            }
            app.layout.get("command_line:field").copied().unwrap_or(Rect::NOTHING)
        };
        let (idle, typing) = (field(""), field("mesh.node.add"));
        assert!(idle.width() < typing.width(), "idle {idle:?}, typing {typing:?}");
        assert!(idle.width() >= MIN_LINE_W - 0.5, "{idle:?}");
    }
}
