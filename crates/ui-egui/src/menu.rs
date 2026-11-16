//! The item menu: what a right click on a browser row offers, and on a datum or a partition's cut
//! face in the viewport. Edit… opens the item's edit form, Hide/Show and Isolate change the view
//! (`view.*`), Active in step switches a boundary condition or load on or off per step
//! (`model.step.activate`), Rename… opens `model.rename` (every reference follows the new
//! name), Apply… applies a mesh plan or relaxation field (its apply form, filled in), Restore
//! shows a saved view (`view.restore`), and Delete… asks first (a seed's and a rigid wall's Delete
//! runs at once: undo brings it back) ([`crate::files::Prompt::Delete`]): it lists what uses the
//! item, and refuses, or for a datum deletes its partitions with it.
//! Every entry runs the same command as the ribbon, the forms and scripts.

use fieldcraft_engine::{ItemKind, ViewItem};
use serde::Serialize;
use serde_json::{Value, json};

use crate::ModelApp;
use crate::browser::Def;
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

/// What the menu is for.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Target {
    /// A body of the open design.
    Body(String),
    /// An imported surface.
    Surface(String),
    /// A model definition, datum or partition.
    Def(Def, String),
    /// The seed on a design edge or face.
    Seed(fieldcraft_model::GeometryRef),
    /// A saved view (`view.save`).
    View(String),
    /// A saved render (`results.render.save`).
    Render(String),
    /// A rigid wall of an explicit dynamics step: the step, and the wall's index in its walls.
    RigidWall { step: String, index: usize },
}

impl Target {
    /// The viewport item it is, when the viewport draws it on its own.
    pub fn item(&self) -> Option<ViewItem> {
        let (kind, name) = match self {
            Target::Body(n) => (ItemKind::Body, n),
            Target::Surface(n) => (ItemKind::Surface, n),
            Target::Def(def, n) => (def.item_kind()?, n),
            Target::Seed(_) | Target::View(_) | Target::Render(_) | Target::RigidWall { .. } => return None,
        };
        Some(ViewItem { kind, name: name.clone() })
    }
}

/// A menu entry chosen.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Action {
    Edit,
    /// Rename a definition (`model.rename`'s form, filled in).
    Rename,
    /// Suppress a datum or partition, or bring it back.
    Suppress(bool),
    Hide,
    Show,
    Isolate,
    ShowAll,
    /// Switch a boundary condition or load on or off in a step.
    Active {
        step: String,
        active: bool,
    },
    /// Show a saved view.
    Restore,
    /// Apply a mesh plan or a relaxation field to the mesh (its apply form, filled in).
    Apply,
    Delete,
}

/// Whether `target` is active in each step (boundary conditions and loads; empty for the rest).
fn activity(app: &ModelApp, target: &Target) -> Vec<(String, bool)> {
    let steps = app.session.model.steps();
    match target {
        Target::Def(Def::BoundaryCondition, n) => steps.iter().map(|s| (s.name.clone(), s.boundary_conditions.contains(n))).collect(),
        Target::Def(Def::Load, n) => steps.iter().map(|s| (s.name.clone(), s.loads.contains(n))).collect(),
        _ => Vec::new(),
    }
}

/// The hover text of `action` on `target`, in words: what the entry does, never a command id.
/// `None` for an entry the menu never shows for `target`.
pub fn tip(target: &Target, action: &Action) -> Option<String> {
    let noun = match target {
        Target::Body(_) => "body".to_string(),
        Target::Surface(_) => "surface".to_string(),
        Target::Def(def, _) => def.label().to_lowercase(),
        Target::Seed(_) => "seed".to_string(),
        Target::View(_) => Saved::View.label().to_lowercase(),
        Target::Render(_) => Saved::Render.label().to_lowercase(),
        Target::RigidWall { .. } => "rigid wall".to_string(),
    };
    let tip = match (action, target) {
        (Action::Edit, Target::Body(_) | Target::Surface(_) | Target::View(_)) => return None,
        (Action::Edit, Target::Render(_)) => format!("Change this {noun} in its form: saving replaces it"),
        (Action::Edit, _) => format!("Change this {noun} in its form"),
        (Action::Apply, Target::Def(..)) => format!("Apply this {noun} to the mesh, in its form"),
        (Action::Rename, Target::Def(..)) => format!("Rename this {noun}: every definition that refers to it follows"),
        (Action::Suppress(true), Target::Def(..)) => format!("Keep this {noun}, but leave it out of the model"),
        (Action::Suppress(false), Target::Def(..)) => format!("Bring this {noun} back into the model"),
        (Action::Restore, Target::View(_)) => "Show the model as this view saved it".to_string(),
        (Action::Hide, _) => format!("Hide this {noun} in the viewport"),
        (Action::Show, _) => format!("Show this {noun} in the viewport"),
        (Action::Isolate, _) => format!("Show only this {noun}"),
        (Action::ShowAll, _) => "Show everything that is hidden".to_string(),
        (Action::Active { step, .. }, Target::Def(..)) => format!("Apply this {noun} in step {step}, or leave it out of that step"),
        (Action::Delete, Target::Body(_) | Target::Surface(_)) => return None,
        (Action::Delete, Target::Def(..)) => format!("Delete this {noun}: asks first, listing what uses it (Delete key)"),
        (Action::Delete, Target::View(_) | Target::Render(_)) => format!("Delete this {noun}: asks first"),
        (Action::Delete, Target::Seed(_) | Target::RigidWall { .. }) => format!("Delete this {noun} (undo brings it back)"),
        (Action::Apply | Action::Rename | Action::Suppress(_) | Action::Restore | Action::Active { .. }, _) => return None,
    };
    Some(tip)
}

/// The hover text of `action` on `target` when the menu shows it disabled, in words; `None`
/// when the menu never disables it. Edit… is disabled for a kind that cannot be edited, and
/// for a rigid wall its step no longer has.
pub fn disabled_tip(target: &Target, action: &Action) -> Option<&'static str> {
    match (action, target) {
        (Action::Edit, Target::Def(..)) => Some("Not editable: delete it and make it again"),
        (Action::Edit, Target::RigidWall { .. }) => Some("Its step no longer has this rigid wall"),
        _ => None,
    }
}

/// Draw the menu's entries for `target`; the entry clicked.
pub fn show(app: &ModelApp, ui: &mut egui::Ui, target: &Target) -> Option<Action> {
    let mut chosen = None;
    let mut entry = |ui: &mut egui::Ui, enabled: bool, label: &str, action: Action| {
        let mut row = kit::MenuItem::new(label, label, "");
        row.enabled = enabled;
        let mut response = row.show(ui);
        if enabled {
            if let Some(hover) = tip(target, &action) {
                response = response.backing_on_hover_text(hover);
            }
        } else if let Some(hover) = disabled_tip(target, &action) {
            response = response.backing_on_disabled_hover_text(hover);
        }
        if response.clicked() {
            chosen = Some(action);
        }
    };
    if let Target::Def(def, name) = target {
        entry(ui, def.edit_command(&app.session.model, name).is_some(), "Edit…", Action::Edit);
        if def.apply(&app.session.model, name).is_some() {
            entry(ui, true, "Apply…", Action::Apply);
        }
        if def.rename_kind().is_some() {
            entry(ui, true, "Rename…", Action::Rename);
        }
        if let Some((_, suppressed)) = suppression(app, *def, name) {
            entry(ui, true, if suppressed { "Unsuppress" } else { "Suppress" }, Action::Suppress(!suppressed));
        }
    }
    match target {
        Target::Seed(_) | Target::Render(_) => entry(ui, true, "Edit…", Action::Edit),
        Target::View(_) => entry(ui, true, "Restore", Action::Restore),
        Target::RigidWall { step, index } => entry(ui, rigid_wall_params(app, step, *index).is_some(), "Edit…", Action::Edit),
        Target::Def(..) | Target::Body(_) | Target::Surface(_) => {}
    }
    if let Some(item) = target.item() {
        if app.session.shows(item.kind, &item.name) {
            entry(ui, true, "Hide", Action::Hide);
        } else {
            entry(ui, true, "Show", Action::Show);
        }
        entry(ui, true, "Isolate", Action::Isolate);
    }
    if !app.session.view.is_empty() {
        entry(ui, true, "Show All", Action::ShowAll);
    }
    let steps = activity(app, target);
    if !steps.is_empty() {
        let mut item = kit::MenuItem::new("active_in_step", "Active in step", "model.step.activate");
        item.submenu = true;
        kit::Menu::item_button(ui, &item, |ui| {
            kit::scope(ui, |ui| {
                for (step, on) in steps {
                    let mut checked = on;
                    let action = Action::Active { step: step.clone(), active: !on };
                    let mut response = kit::checkbox(ui, &step, &mut checked);
                    if let Some(hover) = tip(target, &action) {
                        response = response.backing_on_hover_text(hover);
                    }
                    if response.changed() {
                        chosen = Some(Action::Active { step, active: checked });
                        ui.close();
                    }
                }
            });
        });
    }
    if let Some(hover) = tip(target, &Action::Delete) {
        kit::Menu::separator(ui);
        let label = if matches!(target, Target::Seed(_) | Target::RigidWall { .. }) { "Delete" } else { "Delete…" };
        let mut row = kit::MenuItem::new("delete", label, "");
        row.danger = true;
        if row.show(ui).backing_on_hover_text(hover).clicked() {
            chosen = Some(Action::Delete);
        }
    }
    if chosen.is_some() {
        ui.close();
    }
    chosen
}

/// Do what the menu entry chosen asks.
pub fn apply(app: &mut ModelApp, target: Target, action: Action) {
    let items = |t: &Target| json!({"items": t.item().into_iter().collect::<Vec<_>>()});
    // Errors show in the status bar.
    let _ = match action {
        Action::Edit => {
            if let Target::Def(def, name) = &target
                && let Some((command, params)) = def.edit(&app.session.model, name, app.session.units())
            {
                app.select_definition(Some((*def, name.clone())));
                app.open_form(&command, params);
            }
            if let Target::Seed(on) = &target
                && let Some(params) = app.seed(on).and_then(|s| seed_params(s, app.session.units()))
            {
                let command = seed_command(on);
                app.select_seed(on.clone());
                app.open_form(command, params);
            }
            if let Target::Render(name) = &target
                && let Some(params) = render_params(app, name)
            {
                app.open_form("results.render.save", params);
            }
            if let Target::RigidWall { step, index } = &target
                && let Some(params) = rigid_wall_params(app, step, *index)
            {
                app.open_form("interaction.rigid_wall.edit", params);
            }
            Ok(Value::Null)
        }
        Action::Apply => {
            if let Target::Def(def, name) = &target
                && let Some((command, params)) = def.apply(&app.session.model, name)
            {
                app.select_definition(Some((*def, name.clone())));
                app.open_form(command, params);
            }
            Ok(Value::Null)
        }
        Action::Restore => match &target {
            Target::View(name) => app.run("view.restore", json!({"name": name})),
            _ => Ok(Value::Null),
        },
        Action::Rename => {
            if let Target::Def(def, name) = &target
                && let Some(kind) = def.rename_kind()
            {
                app.select_definition(Some((*def, name.clone())));
                app.open_form("model.rename", json!({"kind": kind, "from": name, "to": name}));
            }
            Ok(Value::Null)
        }
        Action::Suppress(on) => match &target {
            Target::Def(def, name) => match suppression(app, *def, name) {
                Some((command, _)) => app.run(command, json!({"name": name, "suppressed": on})),
                None => Ok(Value::Null),
            },
            _ => Ok(Value::Null),
        },
        Action::Hide => app.run("view.hide", items(&target)),
        Action::Show => app.run("view.show", items(&target)),
        Action::Isolate => app.run("view.isolate", items(&target)),
        Action::ShowAll => app.run("view.show_all", json!({})),
        Action::Active { step, active } => match target {
            Target::Def(def @ (Def::BoundaryCondition | Def::Load), name) => {
                let key = if def == Def::Load { "loads" } else { "boundary_conditions" };
                app.run("model.step.activate", json!({"step": step, key: [name], "active": active}))
            }
            _ => Ok(Value::Null),
        },
        Action::Delete => match target {
            Target::Def(def, name) => {
                app.ask_delete(def, &name);
                Ok(Value::Null)
            }
            Target::Seed(on) => app.run("mesh.seed.clear", json!({seed_key(&on): [on]})),
            Target::View(name) => {
                app.ask_delete_saved(Saved::View, &name);
                Ok(Value::Null)
            }
            Target::Render(name) => {
                app.ask_delete_saved(Saved::Render, &name);
                Ok(Value::Null)
            }
            Target::RigidWall { step, index } => app.run("interaction.rigid_wall.delete", json!({"step": step, "index": index})),
            Target::Body(_) | Target::Surface(_) => Ok(Value::Null),
        },
    };
}

/// The command that suppresses `name` of `kind` (a datum or a partition), and whether it is
/// suppressed now.
fn suppression(app: &ModelApp, kind: Def, name: &str) -> Option<(&'static str, bool)> {
    let model = &app.session.model;
    match kind {
        Def::Datum => Some(("geometry.datum.suppress", model.datum_suppressed(name))),
        Def::Partition => Some(("geometry.partition.suppress", model.partitions().iter().find(|p| p.name == name)?.suppressed)),
        _ => None,
    }
}

/// The command that seeds what `on` names: `mesh.seed.face` for a face, else `mesh.seed.edge`.
pub fn seed_command(on: &fieldcraft_model::GeometryRef) -> &'static str {
    if matches!(on.entity, fieldcraft_model::Entity::Face(_)) { "mesh.seed.face" } else { "mesh.seed.edge" }
}

/// The parameter [`seed_command`] and `mesh.seed.clear` name `on` in: `faces` or `edges`.
fn seed_key(on: &fieldcraft_model::GeometryRef) -> String {
    format!("{}s", on.entity.kind())
}

/// The parameters of [`seed_command`] that make `seed` as it is, in `units`: its edit form,
/// filled in.
pub fn seed_params(seed: &fieldcraft_model::Seed, units: &fieldcraft_model::DisplaySystem) -> Option<Value> {
    let mut v = serde_json::to_value(seed).ok()?;
    let o = v.as_object_mut()?;
    let on = o.remove(seed.on.entity.kind())?;
    o.insert(seed_key(&seed.on), json!([on]));
    fieldcraft_engine::units::shown_params(seed_command(&seed.on), &v, units).ok()
}

/// `interaction.rigid_wall.edit`'s parameters that keep wall `index` of explicit step `step` as it
/// is, in the display units (its form, filled in); `None` when it is gone.
fn rigid_wall_params(app: &ModelApp, step: &str, index: usize) -> Option<Value> {
    let fieldcraft_model::Procedure::ExplicitDynamic(e) = &app.session.model.step(step)?.procedure else { return None };
    let mut v = serde_json::to_value(e.rigid_walls.get(index)?).ok()?;
    let o = v.as_object_mut()?;
    o.insert("step".into(), json!(step));
    o.insert("index".into(), json!(index));
    fieldcraft_engine::units::shown_params("interaction.rigid_wall.edit", &v, app.session.units()).ok()
}

/// `results.render.save`'s parameters that make render `name` as it is (its form, filled in).
fn render_params(app: &ModelApp, name: &str) -> Option<Value> {
    let v = serde_json::to_value(app.session.view.render(name)?).ok()?;
    fieldcraft_engine::units::shown_params("results.render.save", &v, app.session.units()).ok()
}

/// What the session keeps beside the model's definitions (view state, never undone).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Saved {
    View,
    Render,
}

impl Saved {
    pub fn label(self) -> &'static str {
        match self {
            Saved::View => "Saved view",
            Saved::Render => "Render",
        }
    }

    /// The command that deletes one.
    pub fn command(self) -> &'static str {
        match self {
            Saved::View => "view.delete",
            Saved::Render => "results.render.delete",
        }
    }
}

/// What a deletion deletes: a model definition, or a saved view or render.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(untagged)]
pub enum Deletable {
    Def(Def),
    Saved(Saved),
}

impl Deletable {
    fn label(self) -> &'static str {
        match self {
            Deletable::Def(d) => d.label(),
            Deletable::Saved(s) => s.label(),
        }
    }

    /// The command that deletes `name`, and its parameters.
    fn command(self, name: &str) -> (String, Value) {
        match self {
            Deletable::Def(d) => (d.command("delete"), d.delete_params(name)),
            Deletable::Saved(s) => (s.command().to_string(), json!({"name": name})),
        }
    }
}

/// The question a deletion asks ([`crate::files::Prompt::Delete`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DeleteAsk {
    pub kind: Deletable,
    pub name: String,
    /// What uses it.
    pub users: Vec<String>,
    /// For a boundary condition or load in use: the steps it is active in, where deactivating it
    /// lets it go.
    pub steps: Vec<String>,
}

impl DeleteAsk {
    /// Whether deleting goes ahead with its users: a datum's partitions are deleted with it; any
    /// other item in use is refused.
    pub fn cascades(&self) -> bool {
        self.kind == Deletable::Def(Def::Datum)
    }

    /// Whether it can be deleted at all.
    pub fn allowed(&self) -> bool {
        self.users.is_empty() || self.cascades()
    }

    /// The prompt's way to the fix when it is refused only for the steps it is active in:
    /// "Deactivate in step Load" ("… in 2 steps").
    pub fn deactivate_label(&self) -> Option<String> {
        if self.allowed() || self.steps.len() != self.users.len() {
            return None;
        }
        match self.steps.as_slice() {
            [one] => Some(format!("Deactivate in step {one}")),
            many => Some(format!("Deactivate in {} steps", many.len())),
        }
    }

    /// The question, in words.
    pub fn question(&self) -> String {
        let what = format!("{} `{}`", self.kind.label().to_lowercase(), self.name);
        match (self.users.is_empty(), self.cascades()) {
            (true, _) => format!("Delete {what}?"),
            (false, true) => format!("Delete {what} and what uses it: {}?", self.users.join(", ")),
            (false, false) => format!("{what} is used by {}: change or delete those first.", self.users.join(", ")),
        }
    }
}

impl ModelApp {
    /// Ask before deleting definition `name` of `kind` (the item menu's Delete…, the Delete key).
    pub fn ask_delete(&mut self, kind: Def, name: &str) {
        if !kind.exists(&self.session.model, name) {
            return;
        }
        let m = &self.session.model;
        let users = kind.users(m, name);
        let active = |s: &fieldcraft_model::Step| match kind {
            Def::BoundaryCondition => s.boundary_conditions.contains(name),
            Def::Load => s.loads.contains(name),
            _ => false,
        };
        let steps = m.steps().iter().filter(|s| active(s)).map(|s| s.name.clone()).collect();
        self.prompt = Some(crate::files::Prompt::Delete(DeleteAsk { kind: Deletable::Def(kind), name: name.to_string(), users, steps }));
    }

    /// Ask before deleting saved view or render `name` (the browser's item menu): view state is
    /// not undone. A view a render starts from is refused (delete the render first).
    pub fn ask_delete_saved(&mut self, kind: Saved, name: &str) {
        let view = &self.session.view;
        let users: Vec<String> = match kind {
            Saved::View if view.saved(name).is_some() => {
                view.renders.iter().filter(|r| r.view == name).map(|r| format!("render `{}`", r.name)).collect()
            }
            Saved::Render if view.render(name).is_some() => Vec::new(),
            Saved::View | Saved::Render => return,
        };
        self.prompt =
            Some(crate::files::Prompt::Delete(DeleteAsk { kind: Deletable::Saved(kind), name: name.to_string(), users, steps: Vec::new() }));
    }

    /// Answer [`crate::files::Prompt::Delete`] for a boundary condition or load in use with
    /// [`DeleteAsk::deactivate_label`]: switch it off in each step that uses it
    /// (`model.step.activate`), then ask again, now with Delete.
    pub fn answer_deactivate(&mut self) {
        let Some(crate::files::Prompt::Delete(ask)) = self.prompt.clone() else { return };
        let (Some(_), Deletable::Def(kind)) = (ask.deactivate_label(), ask.kind) else { return };
        let key = if kind == Def::Load { "loads" } else { "boundary_conditions" };
        self.prompt = None;
        for step in &ask.steps {
            if self.run("model.step.activate", json!({"step": step, key: [ask.name], "active": false})).is_err() {
                return;
            }
        }
        self.ask_delete(kind, &ask.name);
    }

    /// Answer [`crate::files::Prompt::Delete`]: delete (with its partitions, for a datum in use)
    /// or leave it.
    pub fn answer_delete(&mut self, delete: bool) {
        let Some(crate::files::Prompt::Delete(ask)) = self.prompt.clone() else { return };
        self.prompt = None;
        if !delete || !ask.allowed() {
            return;
        }
        let (command, mut params) = ask.kind.command(&ask.name);
        if ask.cascades() && !ask.users.is_empty() {
            params["cascade"] = json!(true);
        }
        if self.run(&command, params).is_ok()
            && let Deletable::Def(kind) = ask.kind
            && self.selected.as_ref() == Some(&(kind, ask.name))
        {
            self.selected = None;
        }
    }
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;
    use fieldcraft_render::StandardView;

    use super::*;
    use crate::viewport;

    /// Every hover text the item menu and the viewport menu show says what the entry does in
    /// words: none names a command id (owner's rule: no command ids as UI text).
    #[test]
    fn hover_texts_name_no_command_id() {
        let ids: Vec<&str> = fieldcraft_engine::commands().iter().map(|c| c.id).collect();
        let on = |entity| fieldcraft_model::GeometryRef { body: "Plate".into(), entity };
        let targets = Def::ALL.iter().map(|d| Target::Def(*d, "A".into())).chain([
            Target::Body("Plate".into()),
            Target::Surface("Skin".into()),
            Target::Seed(on(fieldcraft_model::Entity::Edge("E1".into()))),
            Target::Seed(on(fieldcraft_model::Entity::Face("F1".into()))),
            Target::View("Front".into()),
            Target::Render("Stress".into()),
            Target::RigidWall { step: "Crash".into(), index: 0 },
        ]);
        let actions = [
            Action::Edit,
            Action::Rename,
            Action::Suppress(true),
            Action::Suppress(false),
            Action::Hide,
            Action::Show,
            Action::Isolate,
            Action::ShowAll,
            Action::Active { step: "Load".into(), active: true },
            Action::Restore,
            Action::Apply,
            Action::Delete,
        ];
        let mut texts = vec![viewport::CREATE_SET_TIP.to_string()];
        for target in targets {
            texts.extend(actions.iter().filter_map(|a| tip(&target, a)));
            texts.extend(actions.iter().filter_map(|a| disabled_tip(&target, a)).map(str::to_string));
        }
        let named: Vec<(&str, &String)> = texts.iter().flat_map(|t| ids.iter().filter(|id| t.contains(**id)).map(move |id| (*id, t))).collect();
        assert!(named.is_empty(), "{named:?}");
        // Each Def and every target that can be deleted has a delete entry.
        assert!(Def::ALL.iter().all(|d| tip(&Target::Def(*d, "A".into()), &Action::Delete).is_some()));
    }

    /// The window a hover test draws, tall enough for the whole results panel.
    const SCREEN: egui::Vec2 = egui::vec2(1400.0, 1000.0);

    /// One frame of the whole window with `events`, with AccessKit on: every hover text the
    /// frame drew, as the descriptions the kit gives each widget's accessibility node, and the
    /// centre of each node by its label (what a test clicks). A hover the kit does not describe
    /// is not seen: one shown only while disabled, or on a widget with no node (`kit::text`).
    fn hover_frame(app: &mut ModelApp, ctx: &egui::Context, events: Vec<egui::Event>) -> (Vec<String>, Vec<(String, egui::Pos2)>) {
        ctx.enable_accesskit();
        let raw = egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, SCREEN)), events, ..Default::default() };
        let mut out = ctx.run_ui(raw, |ui| {
            app.logic(ui.ctx());
            app.ui(ui);
        });
        out.textures_delta.clear();
        let nodes = out.platform_output.accesskit_update.iter().flat_map(|u| u.nodes.iter().map(|(_, n)| n));
        let (mut hovers, mut labels) = (Vec::new(), Vec::new());
        for node in nodes {
            hovers.extend(node.description().map(str::to_string));
            if let (Some(label), Some(b)) = (node.label().or(node.value()), node.bounds()) {
                let ppp = f64::from(ctx.pixels_per_point());
                labels.push((label.to_string(), egui::pos2(((b.x0 + b.x1) / 2.0 / ppp) as f32, ((b.y0 + b.y1) / 2.0 / ppp) as f32)));
            }
        }
        (hovers, labels)
    }

    /// A click with `button` at `at`: the press, then the release a frame later.
    fn click(at: egui::Pos2, button: egui::PointerButton) -> [Vec<egui::Event>; 2] {
        let press = |pressed| egui::Event::PointerButton { pos: at, button, pressed, modifiers: egui::Modifiers::NONE };
        [vec![egui::Event::PointerMoved(at), press(true)], vec![press(false)]]
    }

    /// Draws `app` a few frames, then for each of `clicks` (a layout region, a label, a button)
    /// clicks with the button on the widget with that label in that region and draws on; every
    /// hover text drawn.
    fn hovers_drawn(app: &mut ModelApp, ctx: &egui::Context, clicks: &[(&str, &str, egui::PointerButton)]) -> Vec<String> {
        let mut hovers = Vec::new();
        let mut labels = Vec::new();
        for _ in 0..3 {
            let (h, l) = hover_frame(app, ctx, Vec::new());
            hovers.extend(h);
            labels = l;
        }
        for (region, label, button) in clicks {
            let area = app.layout.get(*region).copied().unwrap_or_else(|| panic!("{region} not drawn"));
            let at = labels
                .iter()
                .find(|(l, at)| l == label && area.contains(*at))
                .map(|(_, at)| *at)
                .unwrap_or_else(|| panic!("{label} not drawn in {region}: {labels:?}"));
            for events in click(at, *button).into_iter().chain([Vec::new(), Vec::new()]) {
                let (h, l) = hover_frame(app, ctx, events);
                hovers.extend(h);
                labels = l;
            }
        }
        hovers
    }

    /// Every hover text the window draws says what its widget does in words: none names a
    /// command id (owner's rule: no command ids as UI text). Read from drawn frames (each
    /// widget's hover text is its accessibility description), across the browser, its item
    /// menu, the results panel with reactions, a probed node, the XY plot and
    /// streamlines, and the Help menu, so a new hover that names an id fails here.
    #[test]
    fn no_drawn_hover_names_a_command_id() {
        let ids: Vec<&str> = fieldcraft_engine::commands().iter().map(|c| c.id).collect();
        let mut hovers = Vec::new();
        // A solved static model in a copy that updates itself: reactions, a probed node, a saved
        // view, the Help menu, and a boundary condition's item menu.
        let ctx = egui::Context::default();
        let (_tx, events) = crate::inbox::channel();
        let mut app = ModelApp::new(Session::default()).with_updates(crate::updates::UpdateHost {
            unsupported: None,
            automatic: false,
            start: Box::new(|_| {}),
            events,
            now: Box::new(|| 0),
            open_file: Box::new(|_| Ok(())),
        });
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("solve.run", json!({})).unwrap();
        app.run("view.save", json!({"name": "Front"})).unwrap();
        app.probe_node(1).unwrap();
        app.ui.menu = Some("help".into());
        hovers.extend(hovers_drawn(
            &mut app,
            &ctx,
            &[("browser", "Loads & BCs", egui::PointerButton::Primary), ("browser", "Clamp", egui::PointerButton::Secondary)],
        ));
        // An XY plot of an explicit step.
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(Session::default());
        app.open_sample("impact");
        app.wait_loaded();
        app.run("solve.run", json!({})).unwrap();
        app.start("results.xy").unwrap();
        hovers.extend(hovers_drawn(&mut app, &ctx, &[]));
        // Streamlines from given seeds through a duct's flow.
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(Session::default());
        app.run_script(&json!([
            {"command": "model.material.create", "params": {"name": "Oil", "density": 1.0, "fluid": {"type": "newtonian", "viscosity": 0.1}}},
            {"command": "model.fluid_domain.create", "params": {"name": "Duct", "material": "Oil", "min": [0, 0, 0], "max": [2, 1, 0.125], "spacing": 0.125,
                "faces": {"x_min": {"type": "velocity", "velocity": [0.1, 0, 0]}, "x_max": {"type": "pressure", "pressure": 0}}}},
            {"command": "model.step.create", "params": {"name": "Flow", "procedure": {"type": "fluid", "domain": "Duct", "steps": 20, "frames": 2}}},
        ]))
        .unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.streamlines", json!({"seeds": [[0.5, 0.5, 0.0625]]})).unwrap();
        hovers.extend(hovers_drawn(&mut app, &ctx, &[]));
        // A form whose reference field offers to create what it needs (New Section's Material
        // with no material yet): the Create button's hover.
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(Session::default());
        app.open_form("model.section.create", json!({}));
        hovers.extend(hovers_drawn(&mut app, &ctx, &[]));

        let named: Vec<(&str, &String)> = hovers.iter().flat_map(|t| ids.iter().filter(|id| t.contains(**id)).map(move |id| (*id, t))).collect();
        assert!(named.is_empty(), "{named:?}");
        // Nor is one a bare message id (a kit message the catalogue lacks reads as its id).
        let id_like = |t: &str| t.contains('.') && t.split('.').all(|w| !w.is_empty() && w.chars().all(|c| c.is_ascii_lowercase() || c == '_'));
        let bare: Vec<&String> = hovers.iter().filter(|t| id_like(t)).collect();
        assert!(bare.is_empty(), "{bare:?}");
        // Each scene drew what it is here for.
        for reached in [
            "Write the reaction sums to a CSV file",
            "Write the node's results to a CSV file",
            "Ask GitHub Releases whether a newer FieldCraft is published",
            "Rename this boundary condition: every definition that refers to it follows",
            "Write the series to a CSV file",
            "Streamlines from the inlet",
            "New Material: Properties › Material",
        ] {
            assert!(hovers.iter().any(|h| h.starts_with(reached)), "{reached} not drawn: {hovers:?}");
        }
    }

    /// The plate design open in an app with a datum plane `Mid` cutting it (partition `Half`)
    /// and a datum point `P`, framed in an 800 × 600 viewport.
    fn plate() -> (ModelApp, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-menu-{}-{:?}", std::process::id(), std::thread::current().id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        app.run("geometry.datum", json!({"name": "Mid", "by": {"type": "offset", "plane": "yz", "offset": 15}})).unwrap();
        app.run("geometry.partition", json!({"name": "Half", "datum": "Mid"})).unwrap();
        // The partition hides its plane (#295); shown again here to be picked.
        app.run("view.show", json!({"items": [{"kind": "datum", "name": "Mid"}]})).unwrap();
        app.run("geometry.datum.point", json!({"name": "P", "by": {"type": "coordinates", "point": [0, 0, 0]}})).unwrap();
        app.viewport.rect = Some(egui::Rect::from_min_size(egui::pos2(100.0, 50.0), egui::vec2(800.0, 600.0)));
        app.set_view(StandardView::Iso);
        app.fit_view();
        (app, dir)
    }

    /// The Partition form highlights the bodies its plane would cut: a plane through the plate
    /// cuts it, one beyond it nothing, and listed bodies are those listed (#347).
    #[test]
    fn the_partition_form_names_the_bodies_it_would_cut() {
        let (mut app, dir) = plate();
        app.run("geometry.datum", json!({"name": "Far", "by": {"type": "offset", "plane": "yz", "offset": 1000}})).unwrap();
        app.run("geometry.datum", json!({"name": "Cut", "by": {"type": "offset", "plane": "xy", "offset": 1}})).unwrap();
        assert!(app.partition_preview().is_empty(), "no form, nothing");
        app.open_form("geometry.partition", json!({"name": "Q", "datum": "Cut"}));
        let cut = app.partition_preview();
        assert_eq!(cut.len(), 2, "both cells of Half: {cut:?}");
        app.open_form("geometry.partition", json!({"name": "Q", "datum": "Far"}));
        assert!(app.partition_preview().is_empty());
        let one = cut.first().cloned().unwrap();
        app.open_form("geometry.partition", json!({"name": "Q", "datum": "Cut", "bodies": [one.clone()]}));
        assert_eq!(app.partition_preview(), vec![one]);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Rename… opens `model.rename` filled in; running it renames the datum, its partition
    /// follows and the browser keeps it selected under the new name (#290).
    #[test]
    fn rename_from_the_menu_follows_the_selection() {
        let (mut app, dir) = plate();
        apply(&mut app, Target::Def(Def::Datum, "Mid".into()), Action::Rename);
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.rename"));
        let params: Value = serde_json::from_str(&app.dialog.as_ref().unwrap().text).unwrap();
        assert_eq!((params["kind"].as_str(), params["from"].as_str()), (Some("datum"), Some("Mid")));
        app.dialog = None;
        app.run("model.rename", json!({"kind": "datum", "from": "Mid", "to": "Centre"})).unwrap();
        assert_eq!(app.selected, Some((Def::Datum, "Centre".to_string())));
        assert_eq!(app.session.model.partitions()[0].datum, "Centre");
        assert!(Def::Partition.rename_kind().is_none());
        let _ = std::fs::remove_dir_all(dir);
    }

    /// Suppress and Unsuppress on a partition's and a datum's menu run their suppress commands:
    /// the cells merge back and come again; a suppressed datum is not drawn (#291).
    #[test]
    fn suppress_from_the_menu_merges_the_cells_back() {
        let (mut app, dir) = plate();
        let bodies = |app: &ModelApp| app.session.design.as_ref().map_or(0, |d| d.bodies.len());
        assert_eq!(bodies(&app), 2);
        apply(&mut app, Target::Def(Def::Partition, "Half".into()), Action::Suppress(true));
        assert_eq!(bodies(&app), 1, "{:?}", app.status);
        assert_eq!(suppression(&app, Def::Partition, "Half"), Some(("geometry.partition.suppress", true)));
        apply(&mut app, Target::Def(Def::Partition, "Half".into()), Action::Suppress(false));
        assert_eq!(bodies(&app), 2);
        apply(&mut app, Target::Def(Def::Datum, "Mid".into()), Action::Suppress(true));
        assert_eq!(bodies(&app), 1);
        let drawn = |app: &mut ModelApp| {
            let b = app.display().bounds;
            let mut s = fieldcraft_render::Scene::default();
            fieldcraft_render::glyphs::datums(&app.session.model, &b, &|_| true, None, &crate::theme::Tokens::dark().glyph_colors(), &mut s);
            s.labels.into_iter().map(|l| l.text).collect::<Vec<_>>()
        };
        assert_eq!(drawn(&mut app), ["P"], "the suppressed plane is not drawn");
        assert_eq!(suppression(&app, Def::Set, "x"), None);
        let _ = std::fs::remove_dir_all(dir);
    }

    /// A screen point on datum `name`'s drawn lines.
    fn on_datum(app: &mut ModelApp, name: &str) -> egui::Pos2 {
        let bounds = app.display().bounds;
        let mut s = fieldcraft_render::Scene::default();
        fieldcraft_render::glyphs::datums(&app.session.model, &bounds, &|n| n == name, None, &crate::theme::Tokens::dark().glyph_colors(), &mut s);
        let l = s.lines[0];
        viewport::to_screen(app, (l.a + l.b) * 0.5).unwrap()
    }

    /// The menu's entries run the commands the ribbon, forms and scripts run; hiding is view
    /// state (nothing journalled) and a hidden item is not picked.
    #[test]
    fn menu_entries_run_the_same_commands_and_hidden_items_are_not_picked() {
        let (mut app, dir) = plate();
        let mid = Target::Def(Def::Datum, "Mid".into());
        let at = on_datum(&mut app, "Mid");
        assert_eq!(viewport::menu_target(&mut app, at), Some(mid.clone()));
        let journal = app.session.journal.len();
        apply(&mut app, mid.clone(), Action::Hide);
        assert!(!app.session.shows(ItemKind::Datum, "Mid") && app.session.is_dirty());
        assert_eq!(app.session.journal.len(), journal);
        assert_ne!(viewport::menu_target(&mut app, at), Some(mid.clone()), "a hidden datum is not picked");
        apply(&mut app, mid.clone(), Action::Show);
        assert!(app.session.shows(ItemKind::Datum, "Mid"));

        // A hidden body is neither drawn nor picked.
        let cell = app.session.design.as_ref().unwrap().bodies[0].name.clone();
        let centre = app.session.design.as_ref().unwrap().bodies[0].measure().unwrap().centroid;
        let top = viewport::to_screen(&app, crate::results::tests::mm([centre[0], centre[1], 10.0])).unwrap();
        assert!(app.pick_at(top).is_some());
        apply(&mut app, Target::Body(cell.clone()), Action::Isolate);
        assert!(app.session.shows(ItemKind::Body, &cell) && !app.session.shows(ItemKind::Datum, "P"));
        apply(&mut app, Target::Body(cell.clone()), Action::Hide);
        app.display();
        assert!(app.pick.names.iter().all(|n| n.name != cell));
        assert!(app.pick_at(top).is_none_or(|i| !i.describe().contains(cell.as_str())));
        apply(&mut app, Target::Body(cell), Action::ShowAll);
        assert!(app.session.view.is_empty());

        // Edit… opens the datum's edit form filled in; running it is the edit command.
        apply(&mut app, Target::Def(Def::Datum, "P".into()), Action::Edit);
        let d = app.dialog.as_mut().unwrap();
        assert_eq!(d.command, "geometry.datum.point.edit");
        assert!(d.text.contains("[0.0,0.0,0.0]"), "{}", d.text);
        d.text = json!({"name": "P", "by": {"type": "coordinates", "point": [1, 0, 0]}}).to_string();
        app.run_dialog().unwrap();
        assert_eq!(app.session.journal.last().map(|(c, _)| c.as_str()), Some("geometry.datum.point.edit"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A seed is selected, described, edited and deleted from the browser like a definition
    /// (#352).
    #[test]
    fn seeds_are_selected_edited_and_deleted() {
        let (mut app, dir) = plate();
        let bodies = app.run("geometry.bodies", json!({})).unwrap();
        let edge = json!({"body": bodies[0]["name"], "edge": bodies[0]["edges"][0]["name"]});
        app.run("mesh.seed.edge", json!({"edges": [edge], "count": 6, "bias": {"ratio": 4, "toward": "end"}})).unwrap();
        let seed = app.session.model.seeds()[0].clone();
        app.select_seed(seed.on.clone());
        assert!(app.selected.is_none() && crate::context::has_context(&mut app));
        assert_eq!(app.pick.selection.len(), 1, "its edge is picked");
        let rows = crate::describe::seed_rows(&seed, app.session.units());
        let text: Vec<String> = rows.iter().map(|(k, v)| format!("{k}: {v}")).collect();
        assert!(text.iter().any(|r| r.contains('6')) && text.iter().any(|r| r.contains('4')), "{text:?}");
        assert!(!text.iter().any(|r| r.contains('{')), "{text:?}");
        apply(&mut app, Target::Seed(seed.on.clone()), Action::Edit);
        let d = app.dialog.as_ref().unwrap();
        assert_eq!(d.command, "mesh.seed.edge");
        let params: Value = serde_json::from_str(&d.text).unwrap();
        assert_eq!((params["count"].clone(), params["edges"][0]["body"].clone()), (json!(6), bodies[0]["name"].clone()), "{params}");
        app.close_dialog();
        apply(&mut app, Target::Seed(seed.on.clone()), Action::Delete);
        assert!(app.session.model.seeds().is_empty());
        let ctx = egui::Context::default();
        let mut out = ctx.run_ui(Default::default(), |ui| app.logic(ui.ctx()));
        out.textures_delta.clear();
        assert!(app.selected_seed.is_none(), "a deleted seed is no longer selected");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A face seed is described, edited with Seed Faces and deleted like an edge seed.
    #[test]
    fn face_seeds_are_edited_and_deleted_with_their_own_command() {
        let (mut app, dir) = plate();
        let bodies = app.run("geometry.bodies", json!({})).unwrap();
        let face = json!({"body": bodies[0]["name"], "face": bodies[0]["faces"][0]["name"]});
        app.run("mesh.seed.face", json!({"faces": [face], "size": 2.5, "bias": {"ratio": 3, "toward": "both"}})).unwrap();
        let seed = app.session.model.seeds()[0].clone();
        assert!(seed.on_face());
        let rows = crate::describe::seed_rows(&seed, app.session.units());
        assert_eq!(rows[0].0, "Face", "{rows:?}");
        apply(&mut app, Target::Seed(seed.on.clone()), Action::Edit);
        let d = app.dialog.as_ref().unwrap();
        assert_eq!(d.command, "mesh.seed.face");
        let params: Value = serde_json::from_str(&d.text).unwrap();
        assert_eq!((params["faces"][0].clone(), params["bias"]["toward"].clone()), (face, json!("both")), "{params}");
        app.close_dialog();
        apply(&mut app, Target::Seed(seed.on.clone()), Action::Delete);
        assert!(app.session.model.seeds().is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A load in use is not deleted, but its prompt offers the way there: deactivate it in the
    /// steps that use it, after which the same prompt deletes it (A4-10).
    #[test]
    fn a_load_in_use_is_deactivated_then_deleted() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let (load, step) = ("Push", "Load");
        app.ask_delete(Def::Load, load);
        let Some(crate::files::Prompt::Delete(ask)) = app.prompt.clone() else { panic!("{:?}", app.prompt) };
        assert!(!ask.allowed() && ask.steps == [step], "{ask:?}");
        assert_eq!(ask.deactivate_label().as_deref(), Some("Deactivate in step Load"));
        app.answer_deactivate();
        assert!(app.session.model.steps().iter().all(|s| !s.loads.contains(load)));
        let Some(crate::files::Prompt::Delete(ask)) = app.prompt.clone() else { panic!("{:?}", app.prompt) };
        assert!(ask.allowed() && ask.steps.is_empty(), "{ask:?}");
        app.answer_delete(true);
        assert!(!app.session.model.loads().contains_key(load));
    }

    /// Delete… on every item the browser lists runs its kind's delete command once confirmed,
    /// and the item is gone: point masses, springs, size controls (mesh.size.clear), face
    /// merges and absorbs, splits, stitches, defeatures; a saved view and a render (view state,
    /// asked first too); a rigid wall at once, by its index (#A1-01–03, #A1-05, #A1-08).
    #[test]
    fn delete_on_each_kind_runs_its_delete_command() {
        let mut app = crate::browser::tests::every_kind();
        // The render first: the view it starts from is refused while it stands.
        apply(&mut app, Target::View("Front".into()), Action::Delete);
        let Some(crate::files::Prompt::Delete(ask)) = app.prompt.clone() else { panic!("{:?}", app.prompt) };
        assert!(!ask.allowed() && ask.users == ["render `Clip`"], "{ask:?}");
        app.answer_delete(true);
        assert!(app.session.view.saved("Front").is_some());
        let new = [
            Def::Mass,
            Def::Spring,
            Def::SizeControl,
            Def::FaceMerge,
            Def::FaceSplit,
            Def::Stitch,
            Def::Defeature,
            Def::InitialField,
            Def::MeshPlan,
            Def::Relaxation,
        ];
        let mut targets: Vec<Target> =
            crate::browser::tests::listed(&app).into_iter().filter(|t| !matches!(t, Target::Def(k, _) if !new.contains(k))).collect();
        // The render before the view it starts from.
        targets.sort_by_key(|t| matches!(t, Target::View(_)));
        assert_eq!(targets.len(), 14, "{targets:?}");
        for target in targets {
            let before = app.session.journal.len();
            let deleted = match &target {
                Target::Def(kind, name) => {
                    apply(&mut app, target.clone(), Action::Delete);
                    app.answer_delete(true);
                    let journal = app.session.journal.last().map(|(c, _)| c.clone());
                    assert_eq!(journal, Some(kind.command("delete")), "{target:?}: {:?}", app.status);
                    !kind.exists(&app.session.model, name)
                }
                Target::View(_) | Target::Render(_) => {
                    apply(&mut app, target.clone(), Action::Delete);
                    app.answer_delete(true);
                    assert_eq!(app.session.journal.len(), before, "view state is not journalled");
                    !crate::browser::tests::listed(&app).contains(&target)
                }
                Target::RigidWall { .. } => {
                    apply(&mut app, target.clone(), Action::Delete);
                    assert_eq!(app.session.journal.last().map(|(c, _)| c.as_str()), Some("interaction.rigid_wall.delete"));
                    !crate::browser::tests::listed(&app).contains(&target)
                }
                Target::Body(_) | Target::Surface(_) | Target::Seed(_) => continue,
            };
            assert!(deleted, "{target:?} is gone: {:?}", app.status);
        }
    }

    /// Edit… opens each editable kind's form filled in, and running it unchanged changes
    /// nothing; merges, splits and stitches are not editable (made again instead); Restore shows
    /// a saved view; a rigid wall opens interaction.rigid_wall.edit (#A1-05, #A1-08).
    #[test]
    fn edit_and_restore_on_the_new_kinds() {
        let mut app = crate::browser::tests::every_kind();
        for (kind, name, command) in
            [(Def::Mass, "motor", "model.mass.edit"), (Def::Spring, "mount", "model.spring.edit"), (Def::SizeControl, "Fine", "mesh.size.set")]
        {
            let before = definition_of(&app, kind, name);
            apply(&mut app, Target::Def(kind, name.into()), Action::Edit);
            assert_eq!(app.dialog.as_ref().map(|d| d.command), Some(command), "{kind:?}");
            app.run_dialog().unwrap();
            assert_eq!(definition_of(&app, kind, name), before, "{kind:?}");
        }
        apply(&mut app, Target::Def(Def::Defeature, "Small".into()), Action::Edit);
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("geometry.defeature"));
        app.close_dialog();
        // A mesh plan and a relaxation field are applied, not edited; an initial field is mapped
        // again under its name.
        for (kind, name, command, key) in
            [(Def::MeshPlan, "Plan-Base", "mesh.plan.apply", "plan"), (Def::Relaxation, "Relax1", "mesh.node.relax.apply", "field")]
        {
            apply(&mut app, Target::Def(kind, name.into()), Action::Apply);
            let d = app.dialog.as_ref().unwrap();
            let params: Value = serde_json::from_str(&d.text).unwrap();
            assert_eq!((d.command, params[key].as_str()), (command, Some(name)), "{kind:?}");
            app.close_dialog();
        }
        for (kind, name) in [
            (Def::FaceMerge, "Fillet"),
            (Def::FaceSplit, "Band"),
            (Def::Stitch, "Gaps"),
            (Def::InitialField, "Warm"),
            (Def::MeshPlan, "Plan-Base"),
            (Def::Relaxation, "Relax1"),
        ] {
            assert!(kind.edit_command(&app.session.model, name).is_none() && kind.rename_kind().is_none(), "{kind:?}");
        }
        assert_eq!((Def::Mass.rename_kind(), Def::Spring.rename_kind()), (Some("point_mass"), Some("spring")));
        apply(&mut app, Target::Render("Clip".into()), Action::Edit);
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("results.render.save"));
        app.close_dialog();
        // A rigid wall opens its own edit, filled in; running it unchanged changes nothing.
        let wall = |app: &ModelApp| match &app.session.model.step("Impact").map(|s| &s.procedure) {
            Some(fieldcraft_model::Procedure::ExplicitDynamic(e)) => e.rigid_walls.first().cloned(),
            _ => None,
        };
        let before = wall(&app);
        apply(&mut app, Target::RigidWall { step: "Impact".into(), index: 0 }, Action::Edit);
        let d = app.dialog.as_ref().unwrap();
        let params: Value = serde_json::from_str(&d.text).unwrap();
        assert_eq!(
            (d.command, params["step"].as_str(), params["index"].as_u64()),
            ("interaction.rigid_wall.edit", Some("Impact"), Some(0)),
            "{params}"
        );
        app.run_dialog().unwrap();
        assert!(app.dialog.is_none() && before.is_some() && wall(&app) == before);
        app.session.view.current = None;
        apply(&mut app, Target::View("Front".into()), Action::Restore);
        assert_eq!(app.session.view.current.as_deref(), Some("Front"), "{:?}", app.status);
    }

    fn definition_of(app: &ModelApp, kind: Def, name: &str) -> Option<Value> {
        crate::browser::definition(&app.session.model, kind, name)
    }

    /// Delete asks first: a datum a partition uses names it and deletes both on confirmation;
    /// a definition in use elsewhere is refused.
    #[test]
    fn deleting_lists_what_uses_it_and_cascades_or_refuses() {
        let (mut app, dir) = plate();
        apply(&mut app, Target::Def(Def::Datum, "Mid".into()), Action::Delete);
        let Some(crate::files::Prompt::Delete(ask)) = app.prompt.clone() else { panic!("{:?}", app.prompt) };
        assert_eq!(ask.users, ["partition `Half`"]);
        assert!(ask.allowed() && ask.question().contains("partition `Half`"));
        // Cancel changes nothing; Delete deletes the datum and its partition in one step.
        app.answer_delete(false);
        assert!(app.prompt.is_none() && app.session.model.datum("Mid").is_some());
        app.selected = Some((Def::Datum, "Mid".into()));
        app.ask_delete(Def::Datum, "Mid");
        app.answer_delete(true);
        assert_eq!(app.session.journal.last(), Some(&("geometry.datum.delete".to_string(), json!({"name": "Mid", "cascade": true}))));
        assert!(app.session.model.partitions().is_empty() && app.selected.is_none());

        // A material a section uses is refused: no Delete.
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.ask_delete(Def::Material, "Steel");
        let Some(crate::files::Prompt::Delete(ask)) = app.prompt.clone() else { panic!() };
        assert!(!ask.allowed() && !ask.users.is_empty());
        let journal = app.session.journal.len();
        app.answer_delete(true);
        assert!(app.session.model.materials().contains_key("Steel"));
        assert_eq!(app.session.journal.len(), journal);
        // Active in step switches a boundary condition per step.
        let bc = app.session.model.boundary_conditions().keys().next().unwrap().clone();
        let step = app.session.model.steps()[0].name.clone();
        apply(&mut app, Target::Def(Def::BoundaryCondition, bc.clone()), Action::Active { step: step.clone(), active: false });
        assert!(!app.session.model.step(&step).unwrap().boundary_conditions.contains(&bc));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
