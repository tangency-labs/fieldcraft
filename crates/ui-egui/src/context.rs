//! The cards on the right: shown only while there is context, as in the other Craft apps. Each is
//! a fixed floating kit card over the viewport, as the browser is on the left, with the kit's
//! header (its icon, title and close button) and its body scrolling inside. At most two show,
//! stacked by the kit (`Panel::stack_first_bounds`, then `stack_next_bounds`; sharing the column
//! among more is solidcraft#17464): first the newest of the open command's form (the kit's feature
//! card, with what a command misses to run), the reply of a report (Check Model, Model Info) and
//! the properties of the definition, seed or mesh selected (with Edit and Delete), the one it
//! displaced closed ([`opened`]); second the results controls while results are shown, or a
//! report run while a form is open (a form is never displaced). Cards do not move until the kit's
//! panels do (solidcraft#17363). Everything acts through commands.

use egui::{Rect, RichText, Vec2, pos2, vec2};
use serde_json::{Value, json};

use crate::browser::Def;
use crate::theme::Tokens;
use crate::{ModelApp, selection};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

/// The cards' width.
pub const WIDTH: f32 = 310.0;

/// The form card's region in `ui.inspect`'s layout.
const FORM: &str = "context:form";

/// Whether there is a card to show.
pub fn has_context(app: &mut ModelApp) -> bool {
    app.dialog.is_some()
        || app.refused.is_some()
        || app.selected.is_some()
        || app.selected_seed.is_some()
        || app.selected_mesh
        || app.reply.is_some()
        || app.plot().is_some()
}

/// A card the right side may show: the open command's form, which draws itself, or another.
enum Card {
    Form,
    Info(Info),
}

/// A card other than the form, on the kit's floating card with its header and a body.
enum Info {
    /// A command that cannot run yet.
    Refused,
    Reply,
    /// The definition selected in the browser, by its rows.
    Properties(Def, String, std::sync::Arc<Vec<(String, String)>>),
    Seed(fieldcraft_model::GeometryRef),
    Mesh,
    Results(std::sync::Arc<crate::results::Plot>),
}

impl Card {
    /// Its region in `ui.inspect`'s layout, and the id of its card.
    fn key(&self) -> &'static str {
        match self {
            Card::Form => FORM,
            Card::Info(info) => info.key(),
        }
    }
}

impl Info {
    fn key(&self) -> &'static str {
        match self {
            Info::Refused => "context:refused",
            Info::Reply => "context:reply",
            Info::Properties(..) => "context:properties",
            Info::Seed(_) => "context:seed",
            Info::Mesh => "context:mesh",
            Info::Results(_) => "context:results",
        }
    }

    /// Its header's icon and title.
    fn heading(&self, app: &ModelApp) -> (String, String) {
        match self {
            Info::Refused => {
                app.refused.as_ref().map(|u| (u.command.to_string(), format!("{}: not yet", crate::command_label(u.command)))).unwrap_or_default()
            }
            Info::Reply => app.reply.as_ref().map(|r| (r.command.clone(), r.label.clone())).unwrap_or_default(),
            Info::Properties(kind, name, _) => (kind.icon().to_string(), format!("{}: {name}", kind.label())),
            Info::Seed(_) => ("mesh.seed.edge".into(), "Seed".into()),
            Info::Mesh => ("mesh.tet.generate".into(), "Mesh".into()),
            Info::Results(_) => ("results.deformed".into(), "Results".into()),
        }
    }
}

/// What opened a card for the first place.
pub(crate) enum Opened {
    /// A command's form, or the card saying what a command misses to run.
    Form,
    /// A command's reply.
    Reply,
    /// A definition, seed or mesh selected in the browser.
    Selection,
}

/// The newest card leads, so the card it displaces closes, not hidden behind it (review2 of
/// pr/ux-right-cards): a form closes the reply; a reply closes the "not yet" card; a selection
/// closes both. A form is never displaced, so its typed fields stay: a reply takes the other place
/// beside it ([`cards`]), and a selection feeds its Target row. Closing a card never deselects: a
/// form over a selection takes the first place and the selection stays (the item a Rename… or
/// Edit… form acts on). A card's own row selecting what it names (a report's Show) opens nothing,
/// so that card stays ([`mark_selected`]).
pub(crate) fn opened(app: &mut ModelApp, what: Opened) {
    match what {
        Opened::Form => app.reply = None,
        Opened::Reply => app.refused = None,
        Opened::Selection => {
            app.reply = None;
            app.refused = None;
        }
    }
}

/// Select definition `selected` (or with `None` deselect it) and nothing else: no card opens or
/// closes. A browser click does this after [`opened`]; a card's own row does only this, so the
/// card stays open with its other rows.
pub(crate) fn mark_selected(app: &mut ModelApp, selected: Option<(Def, String)>) {
    if selected.is_some() {
        app.selected_seed = None;
        app.selected_mesh = false;
    }
    app.selected = selected;
}

/// The cards shown, top down: the first place's card (what a command misses, the form, the reply,
/// else the selection's properties: [`opened`] closed what a newer card displaced, and a selection
/// under a form or a reply stays selected), then the reply when the form holds the first place,
/// else Results while results are shown, else the next of those. The kit stacks two
/// (solidcraft#17464).
fn cards(app: &mut ModelApp) -> Vec<Card> {
    let mut ranked = Vec::new();
    if app.refused.is_some() {
        ranked.push(Card::Info(Info::Refused));
    }
    if app.dialog.is_some() {
        ranked.push(Card::Form);
    }
    if app.reply.is_some() {
        ranked.push(Card::Info(Info::Reply));
    }
    if let Some((kind, name)) = app.selected.clone()
        && let Some(rows) = property_rows(app, kind, &name)
    {
        ranked.push(Card::Info(Info::Properties(kind, name, rows)));
    }
    if let Some(edge) = app.selected_seed.clone().filter(|e| app.seed(e).is_some()) {
        ranked.push(Card::Info(Info::Seed(edge)));
    }
    if app.selected_mesh {
        ranked.push(Card::Info(Info::Mesh));
    }
    let mut ranked = ranked.into_iter();
    let first = ranked.next();
    // A reply under the form is newer than it (a form opening closes the reply).
    let reply_beside_form = matches!(first, Some(Card::Form)) && app.reply.is_some();
    let mut out: Vec<Card> = first.into_iter().collect();
    match app.plot() {
        Some(plot) if !reply_beside_form => out.push(Card::Info(Info::Results(plot))),
        _ => out.extend(ranked.next()),
    }
    out
}

/// The cards on the right, floating over `ui`'s rect (the viewport) at its right edge, stacked by
/// the kit.
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
    let ctx = ui.ctx().clone();
    let gap = kit::Container::new(&ctx, kit::Surface::Panel).gap();
    // The view's room between the top cards and History; the cards keep a shell gap from its
    // top and left, and the kit's feature card its own from its right and bottom.
    let space = app.shell.panel_space(ui.available_rect_before_wrap());
    let area = Rect::from_min_max(space.min + vec2(gap, gap), space.max.max(space.min + vec2(gap, gap)));
    // The column's right edge where the kit's feature card docks, so every card lines up with
    // the form.
    let right = kit::Panel::new("", egui::Id::NULL, pos2(space.right(), area.top()), space, kit::frame::Side::Right, egui::Frame::NONE)
        .width(WIDTH)
        .bounded_anchor(&ctx)
        .x;
    let shown = cards(app);
    let mut drawn: Vec<Rect> = Vec::new();
    for card in &shown {
        let bounds = match drawn.last() {
            Some(previous) => kit::Panel::stack_next_bounds(area, *previous),
            None if shown.len() > 1 => kit::Panel::stack_first_bounds(area),
            None => area,
        };
        let anchor = pos2(right, bounds.top());
        let rect = match card {
            Card::Form => form(app, ui, &t, bounds, anchor),
            Card::Info(info) => Some(info_card(app, &ctx, &t, info, bounds, anchor)),
        };
        if let Some(rect) = rect {
            app.place(card.key(), rect);
            drawn.push(rect);
        }
    }
    // The viewport's overlays and framing keep clear of them.
    app.shell.right = drawn.clone();
    if let Some(all) = drawn.into_iter().reduce(|a, b| a.union(b)) {
        app.place("context", all);
    }
}

/// The width inside a card's body margins, where a form's fields are laid out (the tests draw
/// fields at it).
#[cfg(test)]
pub(crate) fn inner_width() -> f32 {
    WIDTH - 2.0 * solidcraft_theme::Tokens::get().space_12
}

/// A card other than the form, on the kit's floating card as the form is (so their edges line
/// up), its top right corner at `anchor` and held to `bounds`; its header stays and its body
/// scrolls. Returns where it was drawn.
fn info_card(app: &mut ModelApp, ctx: &egui::Context, t: &Tokens, card: &Info, bounds: Rect, anchor: egui::Pos2) -> Rect {
    let k = solidcraft_theme::Tokens::get();
    let key = card.key();
    let side = kit::frame::Side::Right;
    let frame = kit::Panel::feature_frame(ctx, kit::LayoutMode::current(ctx), side);
    let panel = kit::Panel::new(key, egui::Id::unique(("fc_card", key)), anchor, bounds, side, frame).width(WIDTH);
    // The room the kit's window gives the card under its anchor, inside its frame: egui sizes a
    // window from its last frame, so the body is told its room or it never grows back.
    let gap = kit::Container::new(ctx, kit::Surface::Feature).gap();
    let room = (bounds.bottom() - gap - panel.bounded_anchor(ctx).y - panel.positioned_frame(ctx).total_margin().sum().y).max(0.0);
    let (icon, title) = card.heading(app);
    let mut closed = false;
    let shown = panel.window(ctx).show(ctx, |ui| {
        ui.set_max_height(room);
        ui.spacing_mut().item_spacing.y = 0.0;
        kit::panel::header(ui, &icon, &title, || closed = true);
        egui::ScrollArea::vertical().id_salt(key).min_scrolled_height(0.0).auto_shrink([false, true]).show(ui, |ui| {
            egui::Frame::NONE.inner_margin(egui::Margin::same(k.space_12 as i8)).show(ui, |ui| {
                ui.spacing_mut().item_spacing.y = k.form_gap_y;
                body(app, ui, t, card);
            });
        });
        kit::Container::finish_content(ui, panel.positioned_frame(ui.ctx()));
    });
    if closed {
        close(app, card);
    }
    shown.map_or(Rect::NOTHING, |r| r.response.rect)
}

/// A card's body.
fn body(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, card: &Info) {
    match card {
        Info::Refused => refused(app, ui, t),
        Info::Properties(kind, name, rows) => properties(app, ui, *kind, name, rows),
        Info::Seed(edge) => seed_card(app, ui, edge.clone()),
        Info::Mesh => mesh_card(app, ui, t),
        Info::Results(plot) => crate::results_panel::controls(app, plot, ui, t),
        Info::Reply => crate::reply::show(app, ui, t),
    }
}

/// What a card's close button does.
fn close(app: &mut ModelApp, card: &Info) {
    match card {
        Info::Refused => app.refused = None,
        Info::Properties(..) => app.select_definition(None),
        Info::Seed(_) => app.selected_seed = None,
        Info::Mesh => app.selected_mesh = false,
        // Errors show in the status bar.
        Info::Results(_) => {
            let _ = app.run("results.clear", json!({}));
        }
        Info::Reply => app.reply = None,
    }
}

/// A command that cannot run yet: what it misses, each with a button that starts the command
/// that provides it (as the browser's ghost rows do).
fn refused(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens) {
    let Some(u) = app.refused.clone() else { return };
    ui.label(RichText::new("It needs:").color(t.text_dim));
    let mut start = None;
    for fieldcraft_engine::requirement::Missing(r) in &u.missing {
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new(format!("• {}", r.what())).color(t.text));
            if kit::button(ui, format!("{}…", r.place()), kit::Button::Normal).backing_on_hover_text(crate::command_label(r.provider())).clicked() {
                start = Some(r.provider());
            }
        });
    }
    if let Some(command) = start {
        let _ = app.start(command);
    }
}

/// The open command's form, the kit's feature card with its top right corner at `anchor`, held to
/// `bounds`: its header; typed fields drawn from its schema, each with its help on hover, and what
/// is missing or wrong as the kit's status lines; then OK and Cancel pinned under the fields (or
/// in the header, by the kit's Dialog buttons preference). Never JSON, never the command's id: the
/// parameters' JSON is for scripts, MCP and the control channel. OK runs it; while a job runs (this form's mesher, a solve) OK
/// waits for it and the card shows it as its compute status.
/// A command that applies to part of the model (a set, a boundary condition, a load, a section, a
/// coupling) has a Target row that takes the viewport's selection. Returns where it was drawn.
fn form(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens, bounds: Rect, anchor: egui::Pos2) -> Option<Rect> {
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
    let can_browse = app.has_save_picker();
    let mut browse = false;
    let session = &app.session;
    let selection = &app.pick.selection;
    // The set Use selection will make on OK is offered with the others.
    let pending = app.dialog.as_ref().and_then(|d| d.pending_set.as_ref()).and_then(|s| s.get("name")).and_then(Value::as_str).map(str::to_string);
    let d = app.dialog.as_mut()?;
    let Some(spec) = fieldcraft_engine::find_command(d.command) else {
        app.dialog = None;
        return None;
    };
    // The form's widgets are its own, keyed by this opening: the keyboard focus and typed text of
    // the form it replaces (egui drops focus from a widget no longer drawn) never land here.
    let form_id = egui::Id::unique(("fc_form", spec.id, d.opening));
    let params = serde_json::from_str::<Value>(&d.text).unwrap_or(Value::Null);
    // What is missing or wrong, by the fields' labels: OK waits for it.
    let problems = match (spec.form, params.as_object()) {
        (Some(fields), Some(obj)) => crate::form::messages(fields, obj),
        _ => Vec::new(),
    };
    let missing = problems.iter().map(|(_, m)| m.as_str()).collect::<Vec<_>>().join("\n");
    let (mut invalid, mut create, mut hover_body, mut use_selection) = (0, None, None, false);
    let ctx = ui.ctx().clone();
    let side = kit::frame::Side::Right;
    let frame = kit::Panel::feature_frame(&ctx, kit::LayoutMode::current(&ctx), side);
    // The card keeps its id across openings; its fields are keyed by `form_id`.
    let card = kit::Panel::new(spec.label, egui::Id::unique(FORM), anchor, bounds, side, frame).width(WIDTH);
    let heading = kit::FeatureOptions { icon: spec.id, title: spec.label, drag_key: "", moved: Vec2::ZERO, docked: anchor, apply_label: None };
    let mut out = (false, kit::Action::None);
    let mut attempted = false;
    let shown = card.window(&ctx).show(&ctx, |ui| {
    out = card.feature(ui, heading, false, missing.is_empty(), |ui| {
        if let Some(target) = selection::target(spec.id, &params) {
            kit::section_heading(ui, "Target");
            ui.vertical(|ui| {
                let tip = match target {
                    selection::Target::Members => "The selection becomes the set's members (geometry by name, the mesh by label)",
                    selection::Target::Set | selection::Target::TargetSet => {
                        "Names a set of the selection here: an existing set with the same members, or a new one made on OK"
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
        // The typed fields edit the parameters' object (every command has a form, and the form
        // and the control channel only ever store an object).
        match (spec.form, params.as_object().cloned()) {
            (Some(fields), Some(mut obj)) => {
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
                    selection: Some(selection),
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
                if spec.id == "mesh.export.vol" && can_browse {
                    browse = kit::button(ui, "Browse…", kit::Button::Normal).backing_on_hover_text("Choose the output mesh file").clicked();
                }
            }
            _ => {
                kit::status_line(ui, "These parameters do not fit this form: cancel it and start the command again", kit::StatusTone::Error);
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
        // Read after the fields are drawn: a field holding something that is not a value, or
        // a running job, keeps OK greyed (the kit's `applying`).
        kit::ComputeStatus {
            text: running.as_ref().map(|job| format!("{job} is running")),
            busy: running.is_some(),
            applying: running.is_some() || invalid > 0,
        }
    });
    // The greyed OK says why it waits; a click on it shows every problem, beside its field.
    let ready = invalid == 0 && missing.is_empty() && running.is_none();
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
    });
    let (closed, action) = out;
    let ready = invalid == 0 && missing.is_empty() && running.is_none();
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
    } else if browse {
        app.browse_vol_path();
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
    shown.map(|r| r.response.rect)
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

/// The selected definition's fields (`rows`), with Edit (its form, filled in) and Delete.
fn properties(app: &mut ModelApp, ui: &mut egui::Ui, kind: Def, name: &str, rows: &[(String, String)]) {
    property_grid(ui, "fc_properties", rows);
    ui.add_space(solidcraft_theme::Tokens::get().space_6);
    let (mut edit, mut delete) = (false, false);
    let editing = kind.edit_command(&app.session.model, name);
    ui.horizontal(|ui| {
        if let Some(command) = &editing {
            edit = kit::button(ui, "Edit…", kit::Button::Normal)
                .backing_on_hover_text(format!("{}, filled in", crate::command_label(command)))
                .clicked();
        }
        delete = kit::button(ui, "Delete…", kit::Button::Normal)
            .backing_on_hover_text(format!("Delete this {} (Delete key)", kind.label().to_lowercase()))
            .clicked();
    });
    if edit && let Some((command, params)) = kind.edit(&app.session.model, name, app.session.units()) {
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
    property_grid(ui, "fc_seed_properties", &rows);
    ui.add_space(solidcraft_theme::Tokens::get().space_6);
    let (mut edit, mut delete) = (false, false);
    ui.horizontal(|ui| {
        edit = kit::button(ui, "Edit…", kit::Button::Normal)
            .backing_on_hover_text(format!("{}, filled in", crate::command_label(crate::menu::seed_command(&edge))))
            .clicked();
        delete = kit::button(ui, "Delete", kit::Button::Normal).backing_on_hover_text("Clear this seed (undo brings it back)").clicked();
    });
    let target = crate::menu::Target::Seed(edge);
    if edit {
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
        ui.label(RichText::new(crate::command_label(command)).strong());
        property_grid(ui, &format!("fc_mesh_made_{i}"), rows);
        if kit::button(ui, "Edit…", kit::Button::Normal)
            .backing_on_hover_text(format!("{}: its form with these parameters", crate::command_label(command)))
            .clicked()
        {
            edit = Some((command.clone(), params.clone()));
        }
    }
    if let Some((command, params)) = edit {
        app.open_form(&command, params);
    }
}

/// The most of a card's width its name column takes; a longer name shows whole on hover.
const NAME_SHARE: f32 = 0.45;

/// A card's fields as the kit's form rows (name, then its value), as SolidCraft's cards show
/// theirs, with every value in one column after the longest name; a value too long for its row
/// shows whole on hover.
fn property_grid(ui: &mut egui::Ui, id: &str, rows: &[(String, String)]) {
    let k = solidcraft_theme::Tokens::get();
    let font = egui::FontId::proportional(k.font_label);
    let longest = rows.iter().map(|(name, _)| crate::elide::width(ui, name, &font)).fold(0.0, f32::max);
    let column = (longest.ceil() + k.space_4).min(ui.available_width() * NAME_SHARE);
    ui.push_id(id, |ui| {
        kit::settings_form(ui, column, |ui| {
            for (name, text) in rows {
                kit::row(ui, name, |ui| kit::text(ui, text, kit::Text::Cell).backing_on_hover_text(text.as_str()));
            }
        });
    });
}

/// A field's value as short words, never JSON: numbers as [`crate::describe::quantity`], strings
/// as they are, lists as their items (a list of lists as its rows in parentheses; a list of
/// objects, or one nested too deep, as its count), objects as `key value` pairs; long text is
/// shortened.
pub(crate) fn value_text(v: &Value) -> String {
    const MAX: usize = 160;
    let s = words(v, 0);
    // An object's text starts with its first key: a sentence.
    let s = match v {
        Value::Object(_) => {
            let mut c = s.chars();
            c.next().map_or_else(String::new, |f| f.to_uppercase().chain(c).collect())
        }
        _ => s,
    };
    if s.chars().count() > MAX { format!("{}…", s.chars().take(MAX).collect::<String>()) } else { s }
}

/// Lists longer than this show their length and first items.
const LIST: usize = 12;
/// Nesting worded at most; a deeper list or object is given by its count.
const DEPTH: usize = 2;

fn words(v: &Value, depth: usize) -> String {
    let count = |n: usize, word: &str| if n == 1 { format!("1 {word}") } else { format!("{n} {word}s") };
    // A nested list or object is bracketed by parentheses, so its commas stay its own.
    let nested = |x: &Value| match x {
        Value::Array(a) if !a.is_empty() && depth + 1 < DEPTH && !a.iter().any(Value::is_object) => format!("({})", words(x, depth + 1)),
        Value::Object(o) if !o.is_empty() && depth + 1 < DEPTH => format!("({})", words(x, depth + 1)),
        _ => words(x, depth + 1),
    };
    match v {
        Value::Null => "—".into(),
        Value::Bool(b) => (if *b { "yes" } else { "no" }).into(),
        Value::Number(n) => n.as_f64().map_or_else(|| n.to_string(), crate::describe::quantity),
        Value::String(s) => s.clone(),
        Value::Array(a) if a.is_empty() => "none".into(),
        Value::Array(a) if depth >= DEPTH || a.iter().any(Value::is_object) => count(a.len(), "item"),
        Value::Array(a) => {
            let first: Vec<String> = a.iter().take(LIST).map(nested).collect();
            if a.len() > LIST { format!("{} items: {}, …", a.len(), first.join(", ")) } else { first.join(", ") }
        }
        Value::Object(o) if o.is_empty() => "none".into(),
        Value::Object(o) if depth >= DEPTH => count(o.len(), "field"),
        // Pairs are split by semicolons, so only an object value needs its parentheses.
        Value::Object(o) => o
            .iter()
            .map(|(k, x)| format!("{} {}", k.replace('_', " "), if x.is_object() { nested(x) } else { words(x, depth + 1) }))
            .collect::<Vec<_>>()
            .join("; "),
    }
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

    /// A press at `from`, a drag to `to` and the release, a frame each.
    fn drag(app: &mut ModelApp, ctx: &egui::Context, from: egui::Pos2, to: egui::Pos2) {
        let button = |pos, pressed| egui::Event::PointerButton { pos, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
        frame(app, ctx, vec![egui::Event::PointerMoved(from)]);
        frame(app, ctx, vec![button(from, true)]);
        for i in 1..=4 {
            frame(app, ctx, vec![egui::Event::PointerMoved(from + (to - from) * i as f32 / 4.0)]);
        }
        frame(app, ctx, vec![button(to, false)]);
        frame(app, ctx, vec![]);
    }

    /// The right side is floating kit cards, as the browser on the left (R1): with a form and
    /// Results open, each card is inside the viewport a shell gap from its right edge, under the
    /// ribbon and above History, they do not overlap, and the viewport reaches the window's right
    /// edge. The Results rows start where the form's fields do.
    #[test]
    fn the_right_side_is_floating_cards() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(crate::results::tests::solved_cantilever());
        app.run("results.deformed", json!({"step": "Load"})).unwrap();
        app.start("model.material.create").unwrap();
        for _ in 0..3 {
            frame(&mut app, &ctx, vec![]);
        }
        let region = |app: &ModelApp, r: &str| app.layout.get(r).copied().unwrap_or_else(|| panic!("{r} is drawn: {:?}", app.layout.keys()));
        let view = region(&app, "viewport");
        assert!((view.right() - 1280.0).abs() < 1.0, "the viewport reaches the window's right edge: {view:?}");
        let gap = kit::Container::new(&ctx, kit::Surface::Panel).gap();
        let (ribbon, history) = (region(&app, "ribbon"), region(&app, "timeline"));
        let cards = ["context:form", "context:results"].map(|r| region(&app, r));
        for (i, c) in cards.iter().enumerate() {
            assert!(
                view.contains_rect(*c) && (c.right() - cards[0].right()).abs() < 0.5 && view.right() - c.right() >= gap - 0.5,
                "card {i}: {c:?} in {view:?}"
            );
            assert!(c.top() >= ribbon.bottom() + gap - 0.5 && c.bottom() <= history.top() - gap + 0.5, "card {i}: {c:?}, {ribbon:?} {history:?}");
        }
        assert!(cards[1].top() > cards[0].bottom(), "Results under the form: {cards:?}");
        // The Results rows start where the form's fields do (shot 06: 17 pt off).
        let form =
            solidcraft_ui_registry::rect_of("dialog:body:viewport").expect("the form's body").left() + solidcraft_theme::Tokens::get().space_12;
        let step = solidcraft_ui_registry::rect_of("dialog:label:Step").expect("the Results card's Step row").left();
        assert!((step - form).abs() < 0.5, "{step} vs the form's {form}");
    }

    /// The right cards' regions as `ui.inspect` names them, top down.
    const CARDS: [&str; 7] =
        ["context:refused", "context:form", "context:reply", "context:properties", "context:seed", "context:mesh", "context:results"];

    /// The right cards drawn last frame, top down.
    fn drawn(app: &ModelApp) -> Vec<(&'static str, Rect)> {
        let mut out: Vec<(&'static str, Rect)> = CARDS.iter().filter_map(|k| app.layout.get(*k).map(|r| (*k, *r))).collect();
        out.sort_by(|a, b| a.1.top().total_cmp(&b.1.top()));
        out
    }

    /// The cards lie in the room between the ribbon and History, one under the other, each with
    /// its whole header.
    fn in_the_room(app: &ModelApp, cards: &[(&str, Rect)]) {
        let (ribbon, history) = (app.layout["ribbon"], app.layout["timeline"]);
        for (key, card) in cards {
            assert!(card.top() >= ribbon.bottom() - 0.5 && card.bottom() <= history.top() + 0.5, "{key} {card:?} between {ribbon:?} and {history:?}");
            assert!(card.height() >= kit::header_height(), "{key}'s header is cut: {card:?}");
        }
        for pair in cards.windows(2) {
            assert!(pair[1].1.top() >= pair[0].1.bottom(), "they overlap: {cards:?}");
        }
    }

    /// A model with everything a right card shows: a definition and the mesh selected, results
    /// up and a reply.
    fn everything() -> ModelApp {
        let mut app = ModelApp::new(crate::results::tests::solved_cantilever());
        app.run("results.deformed", json!({"step": "Load"})).expect("results");
        app.select_definition(Some((Def::Material, "Steel".into())));
        app.selected_mesh = true;
        let value = json!({"nodes": 3, "elements": 4});
        app.reply = Some(crate::reply::Reply { command: "model.info".into(), label: "Model Info".into(), params: json!({}), value });
        app
    }

    /// `frames` frames of a window of `size`.
    fn frames_at(app: &mut ModelApp, ctx: &egui::Context, size: Vec2, frames: usize) {
        for _ in 0..frames {
            let raw = egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size)), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        }
    }

    /// With every kind of card asked for, two show (until the kit shares a column among more,
    /// solidcraft#17464): the active one, the form, then Results; in a tall window and a short
    /// one, both above History with their whole headers, neither on the other (review of
    /// pr/ux-right-cards: five cards piled on one another).
    #[test]
    fn two_cards_at_most_the_active_one_then_results() {
        for height in [800.0, 420.0] {
            let ctx = egui::Context::default();
            let mut app = everything();
            app.start("model.load.create").expect("New Load");
            frames_at(&mut app, &ctx, vec2(1280.0, height), 3);
            let cards = drawn(&app);
            let keys: Vec<&str> = cards.iter().map(|c| c.0).collect();
            assert_eq!(keys, ["context:form", "context:results"], "at {height} high");
            in_the_room(&app, &cards);
        }
    }

    /// A card shows its body in the room it is given, also after a window too short for it or
    /// another card taking the first place (F1: Results under a newly opened form showed only
    /// its header).
    #[test]
    fn a_card_shows_its_body_after_the_stack_changes() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(crate::results::tests::solved_cantilever());
        app.run("results.deformed", json!({"step": "Load"})).expect("results");
        frames_at(&mut app, &ctx, vec2(1280.0, 420.0), 3);
        frames_at(&mut app, &ctx, vec2(1280.0, 800.0), 3);
        app.start("model.load.create").expect("New Load");
        frames_at(&mut app, &ctx, vec2(1280.0, 800.0), 3);
        let results = app.layout["context:results"];
        let step = solidcraft_ui_registry::rect_of("dialog:row:Step").expect("the Results card's Step row");
        assert!(results.contains_rect(step), "Step {step:?} shows in Results {results:?}");
    }

    /// Without a form the reply is the active card, and without a reply the selection's
    /// properties; Results stays second.
    #[test]
    fn without_a_form_the_reply_then_the_properties_lead() {
        let ctx = egui::Context::default();
        let mut app = everything();
        frames_at(&mut app, &ctx, vec2(1280.0, 800.0), 3);
        let keys: Vec<&str> = drawn(&app).iter().map(|c| c.0).collect();
        assert_eq!(keys, ["context:reply", "context:results"]);
        app.reply = None;
        frames_at(&mut app, &ctx, vec2(1280.0, 800.0), 3);
        let keys: Vec<&str> = drawn(&app).iter().map(|c| c.0).collect();
        assert_eq!(keys, ["context:properties", "context:results"]);
        in_the_room(&app, &drawn(&app));
    }

    /// The keys of the cards drawn after a few frames of the tests' window.
    fn shown(app: &mut ModelApp, ctx: &egui::Context) -> Vec<&'static str> {
        frames_at(app, ctx, vec2(1280.0, 800.0), 3);
        drawn(app).iter().map(|c| c.0).collect()
    }

    /// The solved cantilever with its results up.
    fn with_results() -> ModelApp {
        let mut app = ModelApp::new(crate::results::tests::solved_cantilever());
        app.run("results.deformed", json!({"step": "Load"})).expect("results");
        app
    }

    fn steel() -> Option<(Def, String)> {
        Some((Def::Material, "Steel".into()))
    }

    /// The newest card takes the first place and the one it displaces closes (review2 of
    /// pr/ux-right-cards): Model Info's reply, then Steel selected in the browser, shows Steel's
    /// properties, and the reply is gone, not hidden.
    #[test]
    fn a_selection_after_a_reply_takes_its_place() {
        let ctx = egui::Context::default();
        let mut app = with_results();
        app.run("model.info", json!({})).expect("Model Info");
        assert_eq!(shown(&mut app, &ctx), ["context:reply", "context:results"]);
        app.select_definition(steel());
        assert_eq!(shown(&mut app, &ctx), ["context:properties", "context:results"]);
        assert!(app.reply.is_none(), "the reply is closed");
    }

    /// A selection made while a form is open feeds its Target row: the form keeps its place.
    #[test]
    fn a_selection_leaves_an_open_form() {
        let ctx = egui::Context::default();
        let mut app = with_results();
        app.start("model.load.create").expect("New Load");
        assert_eq!(shown(&mut app, &ctx), ["context:form", "context:results"]);
        app.select_definition(steel());
        assert_eq!(shown(&mut app, &ctx), ["context:form", "context:results"]);
        assert!(app.dialog.is_some() && app.selected == steel(), "the form stays and Steel is selected");
    }

    /// A form opened over a selection's properties takes their place; the definition stays
    /// selected (closing a card never deselects).
    #[test]
    fn a_form_after_a_selection_takes_its_place() {
        let ctx = egui::Context::default();
        let mut app = with_results();
        app.select_definition(steel());
        assert_eq!(shown(&mut app, &ctx), ["context:properties", "context:results"]);
        app.start("model.load.create").expect("New Load");
        assert_eq!(shown(&mut app, &ctx), ["context:form", "context:results"]);
        assert_eq!(app.selected, steel(), "the definition stays selected");
    }

    /// A click as a hand makes it at `at`: move there, press, release, a frame each.
    fn click(app: &mut ModelApp, ctx: &egui::Context, at: egui::Pos2) {
        let button = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
        frame(app, ctx, vec![egui::Event::PointerMoved(at)]);
        frame(app, ctx, vec![button(true)]);
        frame(app, ctx, vec![button(false)]);
        frame(app, ctx, vec![]);
    }

    /// A report arriving while a form is open takes the other place: the form is never
    /// displaced, and what was typed in it stays (review of pr/ux-right-cards 8da16a22: Check
    /// Model cancelled New Load and its typed name).
    #[test]
    fn a_reply_leaves_an_open_form_and_what_was_typed_in_it() {
        let ctx = egui::Context::default();
        let mut app = with_results();
        app.start("model.load.create").expect("New Load");
        assert_eq!(shown(&mut app, &ctx), ["context:form", "context:results"]);
        // The keyboard on its Name field, as a click there puts it (the form keys its widgets by
        // its opening and each field's path).
        let opening = app.dialog.as_ref().map(|d| d.opening).expect("the form");
        let name = egui::Id::unique(("fc_form", egui::Id::unique(("fc_form", "model.load.create", opening)), "name"));
        ctx.memory_mut(|m| m.request_focus(name));
        frame(&mut app, &ctx, vec![]);
        frame(&mut app, &ctx, vec![egui::Event::Text("Push".into())]);
        frame(&mut app, &ctx, vec![]);
        let typed = |app: &ModelApp| app.dialog.as_ref().map(|d| (d.command, d.opening, d.text.contains("Push")));
        let before = typed(&app);
        assert!(before.is_some_and(|(_, _, pushed)| pushed), "typed into the form: {:?}", app.dialog.as_ref().map(|d| &d.text));
        for report in ["model.check", "model.info"] {
            app.run(report, json!({})).expect("the report");
            assert_eq!(shown(&mut app, &ctx), ["context:form", "context:reply"], "{report}");
            assert_eq!(typed(&app), before, "{report}: the form stays, with what was typed in it");
        }
    }

    #[test]
    fn vol_file_field_accepts_a_pasted_path() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.open_form("mesh.export.vol", json!({"path": ""}));
        shown(&mut app, &ctx);
        let opening = app.dialog.as_ref().unwrap().opening;
        let id = egui::Id::unique(("fc_form", egui::Id::unique(("fc_form", "mesh.export.vol", opening)), "path"));
        ctx.memory_mut(|m| m.request_focus(id));
        frame(&mut app, &ctx, vec![]);
        frame(&mut app, &ctx, vec![egui::Event::Paste("/mesh output/coil.vol".into())]);
        let params: Value = serde_json::from_str(&app.dialog.as_ref().unwrap().text).unwrap();
        assert_eq!(params["path"], "/mesh output/coil.vol");
    }

    /// Show on a row of Check Model selects what it names and the report stays open, its other
    /// rows still there to show (review of pr/ux-right-cards 8da16a22: Show closed the report).
    #[test]
    fn show_on_a_reports_row_keeps_the_report_open() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        for name in ["A", "B"] {
            app.run("model.set.create", json!({"name": name, "nodes": [1]})).expect("a set");
        }
        app.run("model.check", json!({})).expect("Check Model");
        assert_eq!(shown(&mut app, &ctx), ["context:reply"]);
        let shows = || {
            let mut shows: Vec<Rect> =
                solidcraft_ui_registry::rects("dialog:button:").into_iter().filter(|(k, _)| k.ends_with(":Show")).map(|(_, r)| r).collect();
            shows.sort_by(|a, b| a.top().total_cmp(&b.top()));
            shows
        };
        let first = shows();
        assert!(first.len() >= 2, "a Show for each problem of the sets: {first:?}");
        click(&mut app, &ctx, first[0].center());
        assert_eq!(app.selected, Some((Def::Set, "A".to_string())), "Show selects the set its row names");
        assert!(app.reply.as_ref().is_some_and(|r| r.command == "model.check"), "the report stays open");
        assert_eq!(shown(&mut app, &ctx), ["context:reply", "context:properties"]);
        // Every row is still laid out in the report's body (it scrolls in its share of the
        // column), the one shown in view.
        let report = app.layout["context:reply"];
        let rest = shows();
        assert_eq!(rest.len(), first.len(), "every row stays in the report");
        assert!(rest.first().is_some_and(|r| report.contains_rect(*r)), "the row shown stays in view: {rest:?} in {report:?}");
    }

    /// The newest card leads: a selection after a "not yet" card takes its place, and the "not
    /// yet" card closes (review of pr/ux-right-cards 8da16a22: it kept the first place).
    #[test]
    fn a_selection_after_a_refused_card_takes_its_place() {
        let ctx = egui::Context::default();
        let mut app = with_results();
        let refused = fieldcraft_engine::find_command("results.xy").and_then(|spec| fieldcraft_engine::Session::default().availability(spec).err());
        app.refuse(refused.expect("XY Plot cannot run on an empty model"));
        assert_eq!(shown(&mut app, &ctx), ["context:refused", "context:results"]);
        app.select_definition(steel());
        assert_eq!(shown(&mut app, &ctx), ["context:properties", "context:results"]);
        assert!(app.refused.is_none(), "the \"not yet\" card is closed");
    }

    /// A press-drag on a card's header moves nothing (the cards are fixed, as the browser is,
    /// until the kit's movable panels, solidcraft#17363), and another context, as `ui.capture`
    /// draws in, puts it in the same place.
    #[test]
    fn a_card_stays_where_it_is() {
        let ctx = egui::Context::default();
        let mut app = everything();
        app.reply = None;
        frames_at(&mut app, &ctx, vec2(1280.0, 800.0), 3);
        let before = drawn(&app);
        let props = app.layout["context:properties"];
        let grip = props.center_top() + egui::vec2(-40.0, kit::header_height() / 2.0);
        drag(&mut app, &ctx, grip, grip + egui::vec2(-300.0, 120.0));
        assert_eq!(drawn(&app), before, "a drag of the header moves nothing");
        assert!(app.selected.is_some(), "a drag of the header is no click on its close");
        let fresh = egui::Context::default();
        frames_at(&mut app, &fresh, vec2(1280.0, 800.0), 3);
        assert_eq!(drawn(&app), before, "another context draws them in the same place");
    }

    /// New Load beside Results: its header, OK and Cancel show whole in its card, and its fields
    /// scroll inside its body, every one reachable (W1 "Right cards").
    #[test]
    fn new_load_beside_results_keeps_its_header_and_buttons() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(crate::results::tests::solved_cantilever());
        app.run("results.deformed", json!({"step": "Load"})).expect("results");
        app.start("model.load.create").expect("New Load");
        frames_at(&mut app, &ctx, vec2(1280.0, 800.0), 3);
        let form = app.layout["context:form"];
        let history = app.layout["timeline"];
        for part in ["dialog:header", "dialog:ok", "dialog:cancel", "dialog:body:viewport"] {
            let rect = solidcraft_ui_registry::rects(part).first().map(|r| r.1).unwrap_or_else(|| panic!("{part} is drawn"));
            assert!(form.contains_rect(rect) && rect.bottom() <= history.top(), "{part} {rect:?} in {form:?}");
        }
        let body = solidcraft_ui_registry::rect_of("dialog:body:viewport").expect("the form's body");
        assert!(body.height() >= solidcraft_theme::Tokens::get().control_height, "the fields have room to scroll in: {body:?}");
    }

    /// The Target row's note sits under its button on a line of its own, never wrapped beside
    /// it (A3-10).
    #[test]
    fn the_target_note_is_under_its_button() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.start("model.set.create").unwrap();
        frame(&mut app, &ctx, vec![]);
        let raw = egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0))), ..Default::default() };
        let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
        out.textures_delta.clear();
        let rect = |text: &str| {
            out.shapes.iter().find_map(|c| match &c.shape {
                egui::Shape::Text(t) if t.galley.text() == text => Some(t.galley.rect.translate(t.pos.to_vec2())),
                _ => None,
            })
        };
        let button = rect("Use selection").expect("the button");
        let note = rect("Pick in the viewport, or choose below").expect("the note, on one line");
        assert!(note.top() >= button.bottom(), "{note:?} beside {button:?}");
    }

    /// A card's values start at one column, after its longest name (A3-5).
    #[test]
    fn card_values_share_one_column() {
        let ctx = egui::Context::default();
        let rows: Vec<(String, String)> =
            [("Elasticity", "isotropic"), ("Young's modulus", "210000 MPa"), ("Poisson's ratio", "0.3"), ("Density", "7.85e-9 t/mm³")]
                .iter()
                .map(|(a, b)| (a.to_string(), b.to_string()))
                .collect();
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            ui.set_max_width(inner_width());
            solidcraft_ui_registry::clear_handles();
            property_grid(ui, "t", &rows);
        });
        out.textures_delta.clear();
        let lefts: Vec<f32> =
            rows.iter().filter_map(|(k, _)| solidcraft_ui_registry::rect_of(&format!("dialog:control:{k}"))).map(|r| r.left()).collect();
        assert_eq!(lefts.len(), rows.len());
        assert!(lefts.iter().all(|l| (l - lefts[0]).abs() < 0.5), "{lefts:?}");
    }

    /// Every text `shapes` draws inside `within`.
    fn texts_in(shapes: &[egui::epaint::ClippedShape], within: egui::Rect) -> Vec<String> {
        fn walk(shape: &egui::Shape, within: egui::Rect, out: &mut Vec<String>) {
            match shape {
                egui::Shape::Text(t) if within.intersects(t.galley.rect.translate(t.pos.to_vec2())) => out.push(t.galley.text().to_string()),
                egui::Shape::Vec(v) => v.iter().for_each(|s| walk(s, within, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        shapes.iter().for_each(|c| walk(&c.shape, within, &mut out));
        out
    }

    /// Opens each `(command, parameters)` form on the meshed cantilever, one over the last as a
    /// click on another entry does, and lists every text drawn in its card that is JSON or a
    /// command id: "Edit as JSON", raw parameter documentation, an id under the title, a value
    /// shown as JSON text (owner, 2026-11-15).
    fn json_or_ids_shown(forms: impl IntoIterator<Item = (&'static fieldcraft_engine::CommandSpec, serde_json::Map<String, Value>)>) -> Vec<String> {
        let mut session = fieldcraft_engine::Session::default();
        session.run_script(&fieldcraft_engine::sample::cantilever()).expect("the cantilever sample");
        let ids: Vec<&str> = fieldcraft_engine::commands().iter().map(|c| c.id).collect();
        let mut seen = Vec::new();
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(session);
        for (n, (spec, filled)) in forms.into_iter().enumerate() {
            app.open_form(spec.id, Value::Object(filled));
            frame(&mut app, &ctx, vec![]);
            let raw =
                egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
            let Some(panel) = app.layout.get("context").copied() else {
                seen.push(format!("{} (form {n}): no form drawn", spec.id));
                continue;
            };
            for text in texts_in(&out.shapes, panel) {
                if text.contains(['{', '[']) || text.contains("JSON") || ids.iter().any(|id| text.contains(id)) {
                    seen.push(format!("{} (form {n}): {text:?}", spec.id));
                }
            }
        }
        seen
    }

    /// The forms holding every shape of field ([`crate::form::tests::representatives`]: each
    /// kind, both sides of every name-or, every kind of option and variant), each with a value
    /// in every field, show no JSON text and no command id in their card as a window shows it.
    /// An empty field draws nothing, so it would hide a value shown as JSON text. A field
    /// scrolled out of the card is drawn by `form::tests::every_field_has_a_widget`, which
    /// checks every example of every form's fields for JSON text.
    #[test]
    fn no_form_shows_json_or_a_command_id() {
        let seen = json_or_ids_shown(crate::form::tests::representatives());
        assert!(seen.is_empty(), "{} texts:\n{}", seen.len(), seen.join("\n"));
    }

    /// Every example of every command's form (every variant and option of each) shows no JSON
    /// text and no command id: the whole walk of [`no_form_shows_json_or_a_command_id`].
    #[test]
    #[ignore = "heavy: cargo xtask test --heavy (CI heavy job)"]
    fn no_example_of_any_form_shows_json_or_a_command_id() {
        let all = fieldcraft_engine::commands()
            .iter()
            .filter_map(|c| c.form.filter(|f| !f.is_empty()).map(|f| (c, f)))
            .flat_map(|(c, f)| crate::form::tests::examples(f).into_iter().map(move |e| (c, e)));
        let seen = json_or_ids_shown(all);
        assert!(seen.is_empty(), "{} texts:\n{}", seen.len(), seen.join("\n"));
    }

    /// Lists and objects read as words, never as JSON.
    #[test]
    fn values_are_worded() {
        assert_eq!(value_text(&json!([1, 2])), "1, 2");
        assert_eq!(value_text(&json!(true)), "yes");
        assert_eq!(value_text(&json!([])), "none");
        assert_eq!(value_text(&json!({"yaw": 30, "pitch_angle": -20})), "Pitch angle −20; yaw 30");
        assert_eq!(value_text(&json!([[0, 0, 1], [1, 0, 0]])), "(0, 0, 1), (1, 0, 0)");
        assert_eq!(value_text(&json!([{"a": 1}, {"a": 2}])), "2 items");
        assert_eq!(value_text(&json!({"to": {"x": 1}, "names": ["A", "B"]})), "Names A, B; to (x 1)");
        for v in [json!({"a": {"b": {"c": {"d": [1, {"e": null}]}}}}), json!([[[1]]]), json!({"s": "x"})] {
            let text = value_text(&v);
            assert!(!text.contains(['{', '}', '[', ']', '"']), "{v} → {text}");
        }
    }

    #[test]
    fn long_values_are_shortened() {
        assert_eq!(value_text(&json!("Steel")), "Steel");
        assert_eq!(value_text(&Value::Null), "—");
        let long = json!((0..500).collect::<Vec<_>>());
        assert!(value_text(&long).starts_with("500 items: 0, 1, 2"));
        let text = json!({"a": "x".repeat(400)});
        assert!(value_text(&text).ends_with('…') && value_text(&text).chars().count() == 161);
    }
}
