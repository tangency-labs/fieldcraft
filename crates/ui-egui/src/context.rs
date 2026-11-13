//! The context panel (right): shown only while there is context, as in the other Craft apps.
//! From the top: the open command's form, the properties of the definition selected in the
//! browser (with Edit and Delete), the results controls while results are shown, and last the
//! reply of a report (Check Model, Model Info), so what the person is working on stays on top. Each section has
//! the kit's panel header (its icon, title and close button); everything acts through commands.

use egui::RichText;
use serde_json::{Value, json};

use crate::browser::Def;
use crate::theme::Tokens;
use crate::{ModelApp, selection};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

/// The panel's width.
pub const WIDTH: f32 = 310.0;

/// Whether the panel has anything to show.
pub fn has_context(app: &mut ModelApp) -> bool {
    app.dialog.is_some()
        || app.refused.is_some()
        || app.selected.is_some()
        || app.selected_seed.is_some()
        || app.selected_mesh
        || app.reply.is_some()
        || app.plot().is_some()
}

pub fn panel(app: &mut ModelApp, ui: &mut egui::Ui) {
    // A form's row hovered last frame no longer is once the form is gone.
    if app.dialog.is_none() && app.pick.form_hover {
        app.pick.form_hover = false;
        app.set_hover(None);
    }
    // A Feature Graph row hovered last frame no longer is once the graph is not shown.
    if app.reply.as_ref().is_none_or(|r| r.command != "geometry.features") {
        app.set_feature_hover(Vec::new());
    }
    if !has_context(app) {
        return;
    }
    let t = Tokens::get(app.ui.dark);
    let k = solidcraft_theme::Tokens::get();
    // The kit's workspace panel along the right edge.
    let host = kit::Container::new(ui.ctx(), kit::Surface::Panel);
    let edge = host.edge("fc_context", kit::Edge::Right).exact_size(WIDTH).resizable(false);
    let frame = host.edge_frame(t.dialog_bg).inner_margin(egui::Margin::same(k.space_8 as i8));
    host.show_edge(ui, edge, frame, |ui| {
        app.place("context", ui.max_rect());
        // The form stays on top with its OK and Cancel in view; the cards below it scroll.
        if app.dialog.is_some() {
            form(app, ui, &t, ui.available_rect_before_wrap());
            ui.add_space(k.space_10);
        }
        egui::ScrollArea::vertical().id_salt("fc_context_scroll").auto_shrink([false, false]).show(ui, |ui| {
            ui.set_width(inner_width());
            if app.refused.is_some() {
                refused(app, ui, &t);
                ui.add_space(k.space_10);
            }
            if let Some((kind, name)) = app.selected.clone() {
                properties(app, ui, kind, &name);
                ui.add_space(k.space_10);
            }
            if let Some(edge) = app.selected_seed.clone() {
                seed_card(app, ui, edge);
                ui.add_space(k.space_10);
            }
            if app.selected_mesh {
                mesh_card(app, ui, &t);
                ui.add_space(k.space_10);
            }
            if let Some(plot) = app.plot() {
                let mut close = false;
                kit::panel::header(ui, "results.deformed", "Results", || close = true);
                ui.add_space(k.space_4);
                crate::results_panel::controls(app, &plot, ui, &t);
                if close {
                    let _ = app.run("results.clear", json!({}));
                }
                ui.add_space(k.space_10);
            }
            if app.reply.is_some() {
                crate::reply::show(app, ui, &t);
            }
        });
    });
}

/// The width inside the panel's margins.
pub(crate) fn inner_width() -> f32 {
    WIDTH - 2.0 * solidcraft_theme::Tokens::get().space_8
}

/// A command that cannot run yet: what it misses, each with a button that starts the command
/// that provides it (as the browser's ghost rows do).
fn refused(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens) {
    let Some(u) = app.refused.clone() else { return };
    let mut close = false;
    kit::panel::header(ui, u.command, &format!("{}: not yet", crate::command_label(u.command)), || close = true);
    ui.add_space(solidcraft_theme::Tokens::get().space_4);
    ui.label(RichText::new("It needs:").color(t.text_dim));
    let mut start = None;
    for fieldcraft_engine::requirement::Missing(r) in &u.missing {
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new(format!("• {}", r.what())).color(t.text));
            if kit::button(ui, format!("{}…", r.place()), kit::Button::Normal).backing_on_hover_text(r.provider()).clicked() {
                start = Some(r.provider());
            }
        });
    }
    if close {
        app.refused = None;
    } else if let Some(command) = start {
        let _ = app.start(command);
    }
}

/// JSON shorter than this is shown pretty-printed in the JSON view.
const MAX_PRETTY: usize = 100_000;

/// The open command's form, the kit's feature card on top of the panel (`column`, the room left
/// in it): its header; typed fields drawn from its schema (or its JSON text, the Advanced view),
/// what is missing or wrong as the kit's status lines, and the command's parameter
/// documentation; then OK and Cancel pinned under the fields (or in the header, by the kit's
/// Dialog buttons preference). OK runs it; while a job runs (this form's mesher, a solve) OK
/// waits for it and the card shows it as its compute status.
/// A command that applies to part of the model (a set, a boundary condition, a load, a section, a
/// coupling) has a Target row that takes the viewport's selection.
fn form(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, column: egui::Rect) {
    let selected = (!app.pick.selection.is_empty()).then(|| app.pick.selection.summary());
    let running = app.busy().map(str::to_string);
    let relax_note = app.relax_note().or_else(|| app.node_move_note());
    let fluid_note = app.fluid_note();
    let units = app.session.units().clone();
    // Escape closes an open menu or the command search, not the form under it.
    let escape_elsewhere = app.ui.menu.is_some() || app.ui.palette;
    // The nodes, elements and element faces picked in the viewport, for Use selection.
    let mut picked = crate::form::Picked::default();
    for item in app.pick.selection.items() {
        match item {
            selection::Item::Node(n) => picked.nodes.push(u64::from(*n)),
            selection::Item::Element(e) => picked.elements.push(u64::from(*e)),
            selection::Item::Face(f) => picked.faces.push((u64::from(f.element), u64::from(f.face))),
            _ => {}
        }
    }
    let session = &app.session;
    // The set Use selection will make on OK is offered with the others.
    let pending = app.dialog.as_ref().and_then(|d| d.pending_set.as_ref()).and_then(|s| s.get("name")).and_then(Value::as_str).map(str::to_string);
    let Some(d) = app.dialog.as_mut() else { return };
    let Some(spec) = fieldcraft_engine::find_command(d.command) else {
        app.dialog = None;
        return;
    };
    // The form's widgets are its own, keyed by this opening: the keyboard focus and typed text of
    // the form it replaces (egui drops focus from a widget no longer drawn) never land here.
    let form_id = egui::Id::unique(("fc_form", spec.id, d.opening));
    let params = serde_json::from_str::<Value>(&d.text).unwrap_or(Value::Null);
    // What is missing or wrong, by the fields' labels: OK waits for it (in either view).
    let problems = match (spec.form, params.as_object()) {
        (Some(fields), Some(obj)) => crate::form::messages(fields, obj),
        _ => Vec::new(),
    };
    let missing = problems.iter().map(|(_, m)| m.as_str()).collect::<Vec<_>>().join("\n");
    let (mut invalid, mut create, mut hover_body, mut use_selection) = (0, None, None, false);
    let card = kit::Panel::new(spec.label, form_id.with("card"), column.right_top(), column, kit::frame::Side::Right, egui::Frame::NONE);
    let heading = kit::FeatureOptions {
        icon: spec.id,
        title: spec.label,
        drag_key: "",
        moved: egui::Vec2::ZERO,
        docked: column.right_top(),
        apply_label: None,
    };
    let (closed, action) = card.feature(ui, heading, false, missing.is_empty(), |ui| {
        kit::text(ui, spec.id, kit::Text::Note("form:command"));
        if let Some(target) = selection::target(spec.id, &params) {
            kit::section_heading(ui, "Target");
            ui.horizontal(|ui| {
                let tip = match target {
                    selection::Target::Members => "The selection becomes the set's members (geometry by name, the mesh by label)",
                    selection::Target::Set | selection::Target::TargetSet => {
                        "Names a set of the selection here: an existing set with the same members, or a new one made on OK (model.set.create)"
                    }
                    selection::Target::Edges => "The selected design edges (Select: Edges) become the edges seeded",
                    selection::Target::Position => {
                        "Puts the point at what is selected: a vertex or node, the centre of a hole (its faces or rim), or the centroid of faces and edges"
                    }
                    selection::Target::DatumPoint => {
                        "Puts the point on what is selected: a vertex, an edge's midpoint, a face's centre, or the centre of a hole or of several items"
                    }
                    selection::Target::Bodies => "The selected bodies (Select: Bodies), or the bodies of the selected faces, become the cells",
                    selection::Target::ExtrudeFaces => "Pick the source face and use it, then pick the target surface and use it",
                    selection::Target::Faces => {
                        "The selected design faces (Shift-click to pick, Select: Faces) become the form's faces, in place of a sphere"
                    }
                };
                let button = match target {
                    selection::Target::ExtrudeFaces if selection::extrude_key(&params) == "source" => "Use as source",
                    selection::Target::ExtrudeFaces => "Use as target surface",
                    _ => "Use selection",
                };
                use_selection = kit::button(ui, button, kit::Button::Enabled(selected.is_some())).backing_on_hover_text(tip).clicked();
                match &selected {
                    Some(text) => ui.add(egui::Label::new(RichText::new(text).color(t.pick_selected)).wrap()),
                    None => ui.add(egui::Label::new(RichText::new("Pick in the viewport, or choose below").color(t.text_dim)).wrap()),
                };
            });
        }
        kit::section_heading(ui, "Parameters");
        let mut advanced = d.advanced;
        if spec.form.is_some() {
            kit::checkbox(ui, "Edit as JSON", &mut advanced).backing_on_hover_text("Advanced: edit the parameters as JSON");
        }
        // The typed fields edit the JSON text's object; text that is not an object stays in the
        // JSON view until it is.
        let object = params.as_object().cloned();
        if advanced != d.advanced && !advanced && object.is_none() {
            d.error = Some("fix the JSON before going back to the fields".into());
            advanced = true;
        }
        // Pretty-printed once, when the JSON view opens.
        if advanced && !d.advanced && params.is_object() && d.text.len() < MAX_PRETTY {
            d.text = serde_json::to_string_pretty(&params).unwrap_or_else(|_| d.text.clone());
        }
        d.advanced = advanced;
        match (spec.form, object) {
            (Some(fields), Some(mut obj)) if !d.advanced => {
                let refused = d.refused.as_ref().map(|(k, p)| (k.as_str(), p.as_str()));
                let cx = crate::form::Context {
                    units: &units,
                    names: &|c| {
                        let mut names = session.names(c);
                        if c == fieldcraft_engine::schema::Collection::Set
                            && let Some(p) = &pending
                            && !names.contains(p)
                        {
                            names.push(p.clone());
                        }
                        names
                    },
                    salt: form_id,
                    touched: &d.touched,
                    attempted: d.attempted,
                    refused,
                    problems: &problems,
                    picked: picked.clone(),
                };
                let out = crate::form::show(ui, t, &cx, fields, &mut obj);
                invalid = out.invalid;
                create = out.create;
                hover_body = out.hover_body;
                if out.changed {
                    d.text = Value::Object(obj).to_string();
                    d.error = None;
                    d.refused = None;
                }
                d.touched.extend(out.touched);
            }
            _ => {
                kit::text_area(ui, "Parameters", form_id.with("json"), &mut d.text, 6, usize::MAX);
            }
        }
        if let Some(e) = &d.error {
            kit::status_line(ui, e, kit::StatusTone::Error);
        } else if d.attempted && !missing.is_empty() {
            kit::status_line(ui, &missing, kit::StatusTone::Error);
        }
        if let Some(note) = &relax_note {
            kit::status_line(ui, note, kit::StatusTone::Dim);
        }
        if let Some((note, problem)) = &fluid_note {
            kit::status_line(ui, note, kit::StatusTone::Dim);
            if let Some(problem) = problem {
                kit::status_line(ui, problem, kit::StatusTone::Error);
            }
        }
        kit::section(ui, "fc_form_takes", "What it takes", false, |ui| {
            kit::text(ui, spec.params, kit::Text::Note("form:takes"));
        });
        // Read after the fields are drawn: a field holding something that is not a value, or
        // a running job, keeps OK greyed (the kit's `applying`).
        kit::ComputeStatus {
            text: running.as_ref().map(|job| format!("{job} is running")),
            busy: running.is_some(),
            applying: running.is_some() || invalid > 0,
        }
    });
    let ready = invalid == 0 && missing.is_empty() && running.is_none();
    // The greyed OK says why it waits; a click on it shows every problem, beside its field.
    let mut attempted = false;
    if !ready && let Some(ok) = solidcraft_ui_registry::rect_of("dialog:ok") {
        let not_ready = match &running {
            Some(job) => format!("{job} is running: wait for it or cancel it"),
            None if invalid > 0 => format!("{missing}\nA field holds something that is not a value").trim().to_string(),
            None => missing.clone(),
        };
        let waits = ui.interact(ok, form_id.with("ok_waits"), egui::Sense::click());
        attempted = running.is_none() && waits.clicked();
        waits.backing_on_hover_text(not_ready);
    }
    if attempted && let Some(d) = app.dialog.as_mut() {
        d.attempted = true;
    }
    let escape = ui.input(|i| i.key_pressed(egui::Key::Escape));
    let cancel = closed || (action == kit::Action::Cancel && !(escape && escape_elsewhere));
    // A body's row under the pointer is highlighted in the viewport, as a hover there is.
    match hover_body.and_then(|b| fieldcraft_model::GeometryRef::body(b).ok()) {
        Some(g) => {
            app.pick.form_hover = true;
            app.set_hover(Some(selection::Item::Geometry(g)));
        }
        None if app.pick.form_hover => {
            app.pick.form_hover = false;
            app.set_hover(None);
        }
        None => {}
    }
    if cancel {
        app.close_dialog();
    } else if let Some((collection, path)) = create {
        app.create_from_form(collection, path);
    } else if action == kit::Action::Ok && ready {
        app.submit_dialog();
    } else if use_selection
        && let Err(e) = app.use_selection_in_dialog()
        && let Some(d) = app.dialog.as_mut()
    {
        d.error = Some(e);
    }
}

/// The selected definition's fields as (label, text) in its form's words and units, read again only
/// when the model changes.
fn property_rows(app: &mut ModelApp, kind: Def, name: &str) -> Option<std::sync::Arc<Vec<(String, String)>>> {
    let key = (app.session.revision(), app.session.units().label(), kind, name.to_string());
    if let Some((k, rows)) = &app.properties
        && *k == key
    {
        return Some(rows.clone());
    }
    let rows = crate::describe::rows(&app.session.model, kind, name, app.session.units())?;
    let rows = std::sync::Arc::new(rows);
    app.properties = Some((key, rows.clone()));
    Some(rows)
}

/// The selected definition's fields, with Edit (its form, filled in) and Delete.
fn properties(app: &mut ModelApp, ui: &mut egui::Ui, kind: Def, name: &str) {
    let Some(rows) = property_rows(app, kind, name) else { return };
    let mut deselect = false;
    kit::panel::header(ui, kind.icon(), &format!("{}: {name}", kind.label()), || deselect = true);
    ui.add_space(solidcraft_theme::Tokens::get().space_4);
    property_grid(ui, "fc_properties", &rows);
    ui.add_space(solidcraft_theme::Tokens::get().space_6);
    let (mut edit, mut delete) = (false, false);
    let editing = kind.edit_command(&app.session.model, name);
    ui.horizontal(|ui| {
        if let Some(command) = &editing {
            edit = kit::button(ui, "Edit…", kit::Button::Normal).backing_on_hover_text(command).clicked();
        }
        delete = kit::button(ui, "Delete…", kit::Button::Normal).backing_on_hover_text(format!("{} (Delete key)", kind.command("delete"))).clicked();
    });
    if deselect {
        app.select_definition(None);
    } else if edit && let Some((command, params)) = kind.edit(&app.session.model, name, app.session.units()) {
        app.open_form(&command, params);
    } else if delete {
        app.ask_delete(kind, name);
    }
}

/// The selected seed: its edge or face, divisions, bias and constraint, with Edit (the Seed Edges
/// or Seed Faces form, filled in) and Delete (`mesh.seed.clear`).
fn seed_card(app: &mut ModelApp, ui: &mut egui::Ui, edge: fieldcraft_model::GeometryRef) {
    let Some(seed) = app.seed(&edge) else { return };
    let rows = crate::describe::seed_rows(seed, app.session.units());
    let mut deselect = false;
    kit::panel::header(ui, "mesh.seed.edge", "Seed", || deselect = true);
    ui.add_space(solidcraft_theme::Tokens::get().space_4);
    property_grid(ui, "fc_seed_properties", &rows);
    ui.add_space(solidcraft_theme::Tokens::get().space_6);
    let (mut edit, mut delete) = (false, false);
    ui.horizontal(|ui| {
        edit = kit::button(ui, "Edit…", kit::Button::Normal).backing_on_hover_text(crate::menu::seed_command(&edge)).clicked();
        delete = kit::button(ui, "Delete", kit::Button::Normal).backing_on_hover_text("mesh.seed.clear (undo brings it back)").clicked();
    });
    let target = crate::menu::Target::Seed(edge);
    if deselect {
        app.selected_seed = None;
    } else if edit {
        crate::menu::apply(app, target, crate::menu::Action::Edit);
    } else if delete {
        crate::menu::apply(app, target, crate::menu::Action::Delete);
    }
}

/// The mesh: its size, then each command that made it with its parameters (in its form's
/// words), with Edit… opening that command's form filled in (#343).
fn mesh_card(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens) {
    let k = solidcraft_theme::Tokens::get();
    let key = (app.session.revision(), app.session.units().label());
    let made = match &app.mesh_made {
        Some((k, made)) if *k == key => made.clone(),
        _ => {
            let units = app.session.units();
            let made: Vec<(String, Value, Vec<(String, String)>)> = app
                .mesh_history()
                .into_iter()
                .map(|(c, p)| {
                    let rows = p.as_object().map(|o| crate::describe::command_rows(&c, o, units)).unwrap_or_default();
                    (c, p, rows)
                })
                .collect();
            let made = std::sync::Arc::new(made);
            app.mesh_made = Some((key, made.clone()));
            made
        }
    };
    let mut deselect = false;
    kit::panel::header(ui, "mesh.tet.generate", "Mesh", || deselect = true);
    ui.add_space(k.space_4);
    let (nodes, elements) = (app.session.model.nodes().len(), app.session.model.elements().len());
    let size = [("Nodes".to_string(), crate::describe::quantity(nodes as f64)), ("Elements".to_string(), crate::describe::quantity(elements as f64))];
    property_grid(ui, "fc_mesh_size", &size);
    let mut edit = None;
    if made.is_empty() && elements > 0 {
        ui.add_space(k.space_4);
        ui.label(RichText::new("Not made by a mesher in this model (imported or built by hand)").size(k.font_caption).color(t.text_dim));
    }
    for (i, (command, params, rows)) in made.iter().enumerate() {
        ui.add_space(k.space_8);
        ui.horizontal(|ui| {
            ui.label(RichText::new(crate::command_label(command)).strong());
            ui.label(RichText::new(command.as_str()).size(k.font_small).color(t.text_dim));
        });
        property_grid(ui, &format!("fc_mesh_made_{i}"), rows);
        if kit::button(ui, "Edit…", kit::Button::Normal).backing_on_hover_text(format!("{command}: its form with these parameters")).clicked() {
            edit = Some((command.clone(), params.clone()));
        }
    }
    if deselect {
        app.selected_mesh = false;
    } else if let Some((command, params)) = edit {
        app.open_form(&command, params);
    }
}

/// Rows of (label, text) as a two-column grid.
/// A card's fields as the kit's form rows (name, then its value), as SolidCraft's cards show
/// theirs; a value too long for its row shows whole on hover.
fn property_grid(ui: &mut egui::Ui, id: &str, rows: &[(String, String)]) {
    ui.push_id(id, |ui| {
        for (name, text) in rows {
            kit::row(ui, name, |ui| kit::text(ui, text, kit::Text::Cell).backing_on_hover_text(text.as_str()));
        }
    });
}

/// A field's value as short text: numbers and strings as they are, long lists as their length
/// and first items, other lists and objects as compact JSON (long ones shortened).
pub(crate) fn value_text(v: &Value) -> String {
    const MAX: usize = 160;
    /// Lists longer than this show their length and first items.
    const LIST: usize = 12;
    let s = match v {
        Value::String(s) => s.clone(),
        Value::Null => "—".into(),
        Value::Array(a) if a.len() > LIST => {
            let first: Vec<String> = a.iter().take(LIST).map(Value::to_string).collect();
            format!("{} items: {}, …", a.len(), first.join(", "))
        }
        other => other.to_string(),
    };
    if s.chars().count() > MAX { format!("{}…", s.chars().take(MAX).collect::<String>()) } else { s }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One frame of a 1280 × 800 window with `events`.
    fn frame(app: &mut ModelApp, ctx: &egui::Context, events: Vec<egui::Event>) {
        let raw = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0))),
            events,
            ..Default::default()
        };
        let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
        out.textures_delta.clear();
    }

    /// The command form is the kit's feature card: its header, its OK and Cancel in the footer
    /// pinned under the fields, inside the context panel. A click on the greyed OK shows what is
    /// missing as the kit's status line; Cancel closes the form.
    #[test]
    fn the_form_is_the_kits_feature_card() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.start("model.material.create").unwrap();
        frame(&mut app, &ctx, vec![]);
        frame(&mut app, &ctx, vec![]);
        let panel = app.layout.get("context").copied().expect("the context panel is drawn");
        let footer = solidcraft_ui_registry::rect_of("dialog:footer").expect("the kit's footer");
        let ok = solidcraft_ui_registry::rect_of("dialog:ok").expect("the kit's OK");
        let cancel = solidcraft_ui_registry::rect_of("dialog:cancel").expect("the kit's Cancel");
        assert!(panel.contains_rect(footer) && footer.contains_rect(ok) && footer.contains_rect(cancel), "{panel:?} {footer:?} {ok:?}");
        assert_eq!(solidcraft_ui_registry::count_of("dialog:ok:enabled"), Some(0.0), "nothing named yet: OK waits");
        assert!(solidcraft_ui_registry::rects("dialog:status:").is_empty(), "problems wait for an edit or a try");
        let click = |at: egui::Pos2, app: &mut ModelApp| {
            let button =
                |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
            frame(app, &ctx, vec![egui::Event::PointerMoved(at)]);
            frame(app, &ctx, vec![button(true)]);
            frame(app, &ctx, vec![button(false)]);
            frame(app, &ctx, vec![]);
        };
        click(ok.center(), &mut app);
        assert!(app.dialog.as_ref().is_some_and(|d| d.attempted), "the click on the greyed OK was a try");
        assert!(!solidcraft_ui_registry::rects("dialog:status:").is_empty(), "what is missing, as the kit's status line");
        let cancel = solidcraft_ui_registry::rect_of("dialog:cancel").expect("the kit's Cancel");
        click(cancel.center(), &mut app);
        assert!(app.dialog.is_none(), "Cancel closes the form");
    }

    #[test]
    fn long_values_are_shortened() {
        assert_eq!(value_text(&json!("Steel")), "Steel");
        assert_eq!(value_text(&json!([1, 2])), "[1,2]");
        assert_eq!(value_text(&Value::Null), "—");
        let long = json!((0..500).collect::<Vec<_>>());
        assert!(value_text(&long).starts_with("500 items: 0, 1, 2"));
        let text = json!({"a": "x".repeat(400)});
        assert!(value_text(&text).ends_with('…') && value_text(&text).chars().count() == 161);
    }
}
