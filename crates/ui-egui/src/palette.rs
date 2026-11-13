//! The command search (S or Cmd/Ctrl+K), where commands are typed, as SolidCraft's.
//!
//! The search lists every catalogue entry and related command by name, id or tip; coming-soon
//! entries are listed greyed with their milestone, and commands missing what they need greyed
//! with what is missing. Enter starts the first that can run. A command typed with its parameters
//! (`mesh.node.add {"label": 1, "xyz": [0, 0, 0]}`) is offered to run as typed. Adapted from
//! SolidCraft's `crates/ui-egui/src/palette.rs` (MIT OR Apache-2.0), see ATTRIBUTION.md.

use egui::Rect;
use fieldcraft_engine::catalogue::{self, CATALOGUE};
use serde_json::Value;
use solidcraft_ui_kit::BackingTooltipExt as _;

use crate::ModelApp;

use solidcraft_ui_kit as kit;

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

/// Run a typed line: `id {json}` runs the command with those parameters, `id` starts it, anything
/// else starts the best match.
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

/// Where the open search field's top left goes over the ribbon's closed `field` in `screen`:
/// its right edge kept, inside the window.
pub fn open_field_at(screen: Rect, field: Rect) -> egui::Pos2 {
    let k = solidcraft_theme::Tokens::get();
    let width = kit::CommandSearch::width();
    egui::pos2((field.right() - width).clamp(screen.left() + k.space_8, (screen.right() - width - k.space_8).max(screen.left())), field.top())
}

/// The command search: the kit's search field opened where the ribbon's closed field is, with
/// the kit's results sheet under it, as SolidCraft's command search opens.
pub fn popup(app: &mut ModelApp, ctx: &egui::Context) {
    if !app.ui.palette {
        return;
    }
    let mut picked: Option<&'static str> = None;
    let mut run_typed = false;
    let k = solidcraft_theme::Tokens::get();
    let width = kit::CommandSearch::width();
    // Opened over the ribbon's field, its right edge kept; centred near the top without one.
    let screen = ctx.content_rect();
    let field = app
        .layout
        .get("ribbon:search")
        .copied()
        .unwrap_or_else(|| Rect::from_center_size(egui::pos2(screen.center().x, 60.0), egui::vec2(width, 26.0)));
    let at = open_field_at(screen, field);
    let key = crate::appbar::shortcut("K");
    // An anchored popup host over the ribbon's field.
    let host = kit::Container::new(ctx, kit::Surface::Popup);
    let shown = host.area(egui::Id::unique("fc_search")).order(egui::Order::Foreground).fixed_pos(at).show(ctx, |ui| {
        ui.set_width(width);
        ui.spacing_mut().item_spacing.y = 0.0;
        let geometry = kit::SearchFieldGeometry { height: field.height(), radius: kit::CommandSearch::field_radius(1.0, false), progress: 1.0 };
        kit::CommandSearch::field(ui, &mut app.palette_text, false, true, true, geometry, Some(&key), MAX_QUERY);
        // Enter while the search is open, as SolidCraft's search reads it: the field holds the focus
        // (it asks for it every frame), so its own lost-focus never reports the key.
        let entered = ui.input(|i| i.key_pressed(egui::Key::Enter));
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
    fn a_typed_line_runs_ids_with_json_and_starts_names() {
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

    /// Commands are typed in the search, as SolidCraft's: Enter runs a command typed with its
    /// parameters, and starts the first match for a name.
    #[test]
    fn enter_in_the_search_runs_a_typed_command_or_starts_the_first_match() {
        let enter = |app: &mut ModelApp, line: &str| {
            let ctx = egui::Context::default();
            app.ui.palette = true;
            app.palette_text = line.into();
            let key = egui::Event::Key { key: egui::Key::Enter, physical_key: None, pressed: true, repeat: false, modifiers: Default::default() };
            for events in [vec![], vec![], vec![key]] {
                let raw = egui::RawInput {
                    screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0))),
                    events,
                    ..Default::default()
                };
                let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
                out.textures_delta.clear();
            }
        };
        let mut app = ModelApp::new(Session::default());
        enter(&mut app, r#"mesh.node.add {"label": 4, "xyz": [1, 2, 3]}"#);
        assert!(app.session.model.node(4).is_some() && !app.ui.palette);
        enter(&mut app, "material");
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.material.create"));
    }
}
