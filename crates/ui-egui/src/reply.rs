//! What a command run from the UI answered, made readable: a one-line summary for the status bar
//! ([`summary`]) and, when the reply holds more than that line, the whole reply in the context
//! panel ([`show`]): Check Model as a list of its diagnostics with a way to each fix, the Step
//! Table as a grid of steps, any other reply as labelled rows and tables. Programmatic callers
//! (the control channel's `engine.execute`, MCP, scripts) still get the reply's JSON.

use egui::{RichText, vec2};
use serde_json::{Map, Value};

use crate::browser::{self, Def};
use crate::frame::{Level, heading};
use crate::results::{components, count, number, with_unit};
use crate::theme::{Tokens, label_font};
use crate::{ModelApp, Tab};
use fieldcraft_engine::results::Quantity;
use fieldcraft_model::DisplaySystem;

/// The reply shown in the context panel.
#[derive(Clone, Debug, PartialEq)]
pub struct Reply {
    pub command: String,
    pub label: String,
    pub params: Value,
    pub value: Value,
}

/// The solvers whose reply is a step's summary ([`solve_summary`]).
const SOLVES: [&str; 6] = ["solve.static", "solve.nonlinear", "solve.heat", "solve.frequency", "solve.buckle", "solve.explicit"];

/// Replies that only read the model, refreshed when it changes; any other reply closes then.
pub const LIVE: [&str; 3] = ["model.check", "model.activity", "geometry.defeature.propose"];

/// The status line's text (after the command's label) is cut to this many characters.
const MAX_LINE: usize = 120;
/// Table rows shown before "… n more".
const MAX_ROWS: usize = 200;
/// Table columns shown at most.
const MAX_COLUMNS: usize = 8;
/// Nesting shown at most; deeper values are compact text.
const MAX_DEPTH: usize = 4;

/// The unit of `q` in `units`; `None` without a system or for a plain number.
fn unit(q: Quantity, units: Option<&DisplaySystem>) -> Option<String> {
    units.and_then(|d| fieldcraft_engine::results::unit_of(q, d))
}

/// One line saying what `command` answered: the reply's own `summary`, a sentence for the
/// commands with a well-known reply, otherwise its top-level numbers and words. Never JSON.
pub fn summary(command: &str, v: &Value) -> String {
    summary_in(command, v, None)
}

/// [`summary`] with the numbers' units in `units`.
pub fn summary_in(command: &str, v: &Value, units: Option<&DisplaySystem>) -> String {
    if let Some(s) = v.get("summary").and_then(Value::as_str) {
        return one_line(s);
    }
    if SOLVES.contains(&command) {
        return one_line(&solve_line(v, units));
    }
    let count = |key: &str| v.get(key).map_or(0, |x| x.as_array().map(Vec::len).or_else(|| x.as_object().map(Map::len)).unwrap_or(0));
    let line = match command {
        "model.check" => {
            let (e, w) = (v["errors"].as_u64().unwrap_or(0), v["warnings"].as_u64().unwrap_or(0));
            if e + w == 0 { "ready to solve, no issues".to_string() } else { format!("{}, {}", plural(e, "error"), plural(w, "warning")) }
        }
        "model.activity" => format!(
            "{}, {}, {}",
            plural(count("steps") as u64, "step"),
            plural(count("boundary_conditions") as u64, "boundary condition"),
            plural(count("loads") as u64, "load")
        ),
        "geometry.features" => {
            let bodies = v["bodies"].as_array().map(Vec::as_slice).unwrap_or_default();
            let n: usize = bodies.iter().map(|b| b["features"].as_array().map_or(0, Vec::len)).sum();
            format!(
                "{} building {}",
                plural(n as u64, "feature"),
                if bodies.len() == 1 { "1 body".to_string() } else { format!("{} bodies", bodies.len()) }
            )
        }
        "geometry.bodies" => match v.as_array().map(Vec::len) {
            Some(0) => "no bodies (open a design first)".to_string(),
            Some(n) => plural(n as u64, "body").replace("bodys", "bodies"),
            None => scalars(v),
        },
        "model.info" => {
            let elements: u64 = v["elements"].as_object().into_iter().flatten().filter_map(|(_, n)| n.as_u64()).sum();
            format!(
                "{}, {}, units {}",
                plural(v["nodes"].as_u64().unwrap_or(0), "node"),
                plural(elements, "element"),
                v["display_units"].as_str().unwrap_or("—")
            )
        }
        "interaction.contact.find" => {
            format!("{} among {}", plural(count("touching") as u64, "touching pair"), plural(v["parts"].as_u64().unwrap_or(0), "part"))
        }
        "results.reactions" => match serde_json::from_value::<Vec<f64>>(v["total"]["force"].clone()) {
            Ok(f) => format!("total reaction force ({})", with_unit(components(&f, 0.0), unit(Quantity::Force, units).as_deref())),
            Err(_) => scalars(v),
        },
        _ => scalars(v),
    };
    if line.is_empty() { "done".into() } else { one_line(&line) }
}

/// The status line after `command` changed the model or the display, written for the person
/// rather than from the command's label: "Showing Load deformed ×2", "Mode 2: 1637 Hz", "Added
/// load Push". Commands without a sentence of their own say "<label> done".
pub fn done(command: &str, label: &str, params: &Value, v: &Value) -> String {
    let view = v.get("view").filter(|x| x.is_object());
    let step = view.and_then(|x| x["step"].as_str()).unwrap_or("the step");
    let line = match command {
        "results.deformed" if view.is_some_and(|x| x["deformed"] == false) => Some(format!("Showing {step} undeformed")),
        "results.deformed" => {
            let scale = view.and_then(|x| x["scale"].as_f64()).map(|f| format!(" ×{}", crate::describe::quantity(f))).unwrap_or_default();
            Some(format!("Showing {step} deformed{scale}"))
        }
        "results.contour" => Some(match v.pointer("/contour/label").and_then(Value::as_str) {
            Some(field) => format!("Showing {field} on {step}"),
            None => format!("Contours off on {step}"),
        }),
        "results.modes" => {
            let mode = view.and_then(|x| x["mode"].as_u64()).unwrap_or(1);
            let hz = v
                .get("modes")
                .and_then(Value::as_array)
                .and_then(|m| m.iter().find(|m| m["mode"].as_u64() == Some(mode)))
                .and_then(|m| m["frequency"].as_f64());
            Some(match hz {
                Some(hz) => format!("Mode {mode}: {} Hz", crate::describe::quantity(hz)),
                None => format!("Showing mode {mode} of {step}"),
            })
        }
        "results.animate" => Some(if view.is_some_and(|x| x["animate"] == true) { "Animating" } else { "Animation stopped" }.to_string()),
        "results.frame" => Some(match view.and_then(|x| x["frame"].as_u64()) {
            Some(f) => format!("Showing {step} frame {f}"),
            None => format!("Showing {step}, last frame"),
        }),
        "results.slice" => Some(match (v["slice"]["axis"].as_str(), v["slice"]["position"].as_f64()) {
            (Some(axis), Some(at)) => {
                let field = v.pointer("/contour/label").and_then(Value::as_str).unwrap_or("The flow");
                format!("{field} on a slice at {axis} = {}", crate::describe::quantity(at))
            }
            _ => format!("Slice off on {step}"),
        }),
        "results.surface" => Some(if v["surface"] == true { "Obstacles coloured by the flow" } else { "Obstacles in their own colours" }.to_string()),
        "results.streamlines" => Some(match v.pointer("/traced/streamlines").and_then(Value::as_u64) {
            Some(n) => format!("Traced {}", plural(n, "streamline")),
            None => "Streamlines off".to_string(),
        }),
        "results.clear" => Some("Back to the model".to_string()),
        "results.xy" => {
            let points = v.get("points").and_then(Value::as_u64).unwrap_or(0);
            let what = v.pointer("/y/label").and_then(Value::as_str).unwrap_or("the series");
            Some(format!("Plotted {what}: {}", plural(points, "point")))
        }
        "file.sample" => {
            let label = v.get("sample").and_then(Value::as_str).and_then(fieldcraft_engine::sample::find).map(|s| s.label);
            label.map(|l| format!("Opened the {l} sample"))
        }
        "solve.run" => v.get("steps").and_then(Value::as_array).map(|s| format!("Solved {}", plural(s.len() as u64, "step"))),
        c if c.starts_with("solve.") => Some(format!("Solved {step}", step = params["step"].as_str().unwrap_or("the step"))),
        _ => None,
    };
    if let Some(line) = line {
        return line;
    }
    // A definition added, changed or deleted, by its kind and name.
    let name = params.get("name").and_then(Value::as_str);
    let def = Def::ALL
        .into_iter()
        .find_map(|kind| ["create", "edit", "delete"].into_iter().find(|verb| kind.command(verb) == command).map(|verb| (kind, verb)));
    match (def, name) {
        (Some((kind, verb)), Some(name)) => {
            let did = match verb {
                "create" => "Added",
                "edit" => "Changed",
                _ => "Deleted",
            };
            format!("{did} {} {name}", kind.label().to_lowercase())
        }
        _ => format!("{label} done"),
    }
}

/// Whether the reply says more than its summary line, so it is worth the panel.
pub fn is_rich(v: &Value) -> bool {
    match v {
        Value::Object(m) => m.iter().any(|(k, x)| k != "summary" && !x.is_null()),
        Value::Array(a) => !a.is_empty(),
        _ => false,
    }
}

fn plural(n: u64, word: &str) -> String {
    if n == 1 { format!("1 {word}") } else { format!("{n} {word}s") }
}

/// A scalar as short text.
fn scalar(v: &Value) -> Option<String> {
    match v {
        Value::Bool(b) => Some(if *b { "yes" } else { "no" }.into()),
        Value::Number(n) => Some(
            n.as_u64().map(|u| u.to_string()).or_else(|| n.as_i64().map(|i| i.to_string())).unwrap_or_else(|| number(n.as_f64().unwrap_or(f64::NAN))),
        ),
        Value::String(s) => Some(s.clone()),
        _ => None,
    }
}

/// The top-level scalars as `key value` pairs (lists as their length).
fn scalars(v: &Value) -> String {
    match v {
        Value::Object(m) => {
            let parts: Vec<String> = m
                .iter()
                .filter_map(|(k, x)| match x {
                    Value::Array(a) => Some(format!("{} {}", k.replace('_', " "), a.len())),
                    _ => scalar(x).map(|s| format!("{} {s}", k.replace('_', " "))),
                })
                .collect();
            parts.join(", ")
        }
        Value::Array(a) => plural(a.len() as u64, "item"),
        other => scalar(other).unwrap_or_default(),
    }
}

/// `s` on one line of at most [`MAX_LINE`] characters.
fn one_line(s: &str) -> String {
    let s = s.split_whitespace().collect::<Vec<_>>().join(" ");
    if s.chars().count() > MAX_LINE { format!("{}…", s.chars().take(MAX_LINE).collect::<String>()) } else { s }
}

/// What a click in the panel asks for, done after drawing.
enum Act {
    Close,
    Select(Def, String),
    Start(&'static str),
    Tab(Tab),
    Rerun,
    /// Pick a design entity in the viewport (a body, face or edge from the Bodies list).
    Pick(fieldcraft_model::GeometryRef),
    /// Open or close a Feature Graph row's details: (body, feature id).
    Open((String, u64)),
    /// Select these faces in the viewport (a Feature Graph row's).
    Faces(Vec<crate::selection::Item>),
    /// Select elements in the viewport and frame them (Mesh Quality's worst and invalid ones).
    Elements(Vec<u32>),
    /// Open Relax Nodes round an element (Mesh Quality's worst).
    Relax(u32),
    /// Switch a feature on (kept in the mesh) or off (defeatured by the named defeature).
    Defeature {
        name: String,
        feature: u64,
        on: bool,
    },
    /// Defeature the proposed features with the named defeature.
    Propose {
        name: String,
        features: Vec<u64>,
    },
}

/// The reply's section of the context panel.
pub fn show(app: &mut ModelApp, ui: &mut egui::Ui, t: &Tokens) {
    let Some(reply) = app.reply.clone() else { return };
    let mut acts = Vec::new();
    heading(ui, t, Level::Title, |ui| {
        ui.label(RichText::new(reply.label.to_uppercase()).font(label_font(ui.ctx(), 12.0)));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.small_button("×").on_hover_text("Close").clicked() {
                acts.push(Act::Close);
            }
            if ui.small_button("Run again").on_hover_text(reply.command.as_str()).clicked() {
                acts.push(Act::Rerun);
            }
        });
    });
    ui.add_space(4.0);
    ui.add(egui::Label::new(RichText::new(summary_in(&reply.command, &reply.value, Some(app.session.units()))).color(t.text_dim)).wrap());
    ui.add_space(4.0);
    let units = Some(app.session.units());
    match reply.command.as_str() {
        "deck.read" => {
            ui.label(
                RichText::new("Model data only. Analysis history was not imported; review unread keywords and warnings before solving.")
                    .color(t.text_dim),
            );
            value(ui, t, egui::Id::new("fc_deck_read"), &reply.value, 0);
        }
        "model.check" => diagnostics(app, ui, t, &reply.value, &mut acts),
        c if SOLVES.contains(&c) => solve_summary(ui, t, &reply.value, units),
        "model.activity" => activity(ui, t, &reply.value),
        "geometry.bodies" => bodies(ui, t, &reply.value, app.session.units(), &mut acts),
        "mesh.quality" => quality(ui, t, &reply.value, &mut acts),
        "geometry.defeature.propose" => defeature_list(ui, t, &reply.value, app.session.units(), &mut acts),
        "geometry.features" => {
            let hovered = feature_graph(ui, t, &reply.value, &app.features_open, &mut acts);
            app.set_feature_hover(hovered);
        }
        _ => value(ui, t, egui::Id::new("fc_reply"), &reply.value, 0),
    }
    for a in acts {
        match a {
            Act::Close => app.reply = None,
            Act::Select(kind, name) => app.selected = Some((kind, name)),
            Act::Start(id) => {
                let _ = app.start(id);
            }
            Act::Tab(tab) => app.ui.select_tab(tab),
            Act::Rerun if reply.command == "deck.read" => {
                app.open_form("deck.read", reply.params.clone());
                app.submit_dialog();
            }
            Act::Rerun => {
                let _ = app.run(&reply.command, reply.params.clone());
            }
            Act::Pick(g) => app.select([crate::selection::Item::Geometry(g)], crate::selection::Mode::Replace),
            Act::Open(key) => {
                if !app.features_open.remove(&key) {
                    app.features_open.insert(key);
                }
            }
            Act::Faces(items) => app.select(items, crate::selection::Mode::Replace),
            Act::Elements(labels) => app.show_elements(&labels),
            Act::Relax(element) => {
                app.show_elements(&[element]);
                app.open_form("mesh.node.relax", serde_json::json!({"around_element": element}));
            }
            Act::Defeature { name, feature, on } => {
                let mut features = defeatured(app, &name);
                features.retain(|f| *f != feature);
                if !on {
                    features.push(feature);
                }
                let _ = set_defeature(app, &name, features);
            }
            Act::Propose { name, features } => {
                let mut all = defeatured(app, &name);
                for f in features {
                    if !all.contains(&f) {
                        all.push(f);
                    }
                }
                let _ = set_defeature(app, &name, all);
            }
        }
    }
}

/// The features the defeature `name` takes out (none when the model has no such defeature).
fn defeatured(app: &ModelApp, name: &str) -> Vec<u64> {
    app.session.model.defeatures().iter().find(|d| d.name == name).map(|d| d.features.clone()).unwrap_or_default()
}

/// Give the defeature `name` these features (none: delete it).
fn set_defeature(app: &mut ModelApp, name: &str, features: Vec<u64>) -> Result<Value, String> {
    match features.is_empty() {
        true if defeatured(app, name).is_empty() => Ok(Value::Null),
        true => app.run("geometry.defeature.delete", serde_json::json!({"name": name})),
        false => app.run("geometry.defeature", serde_json::json!({"name": name, "features": features})),
    }
}

/// Propose Defeaturing: each candidate with a switch (on: in the mesh; off: defeatured), its
/// kind, size and the defeature that takes it out, and a button that defeatures the proposal.
/// The symbol of geometry quantity `q` in `units` (a reply's numbers are in them).
fn geometry_unit(units: &DisplaySystem, q: fieldcraft_model::units::Quantity) -> &str {
    units.symbol(q)
}

fn defeature_list(ui: &mut egui::Ui, t: &Tokens, v: &Value, units: &DisplaySystem, acts: &mut Vec<Act>) {
    use fieldcraft_model::units::Quantity as Q;
    let name = v["name"].as_str().unwrap_or("Defeature").to_string();
    let body = v["body"].as_str().unwrap_or_default();
    let candidates = v["candidates"].as_array().map(Vec::as_slice).unwrap_or_default();
    let proposal: Vec<u64> = v["proposal"].as_array().into_iter().flatten().filter_map(Value::as_u64).collect();
    let pending: Vec<u64> = candidates.iter().filter(|c| c["proposed"] == true && c["suppressed"] != true).filter_map(|c| c["id"].as_u64()).collect();
    if !proposal.is_empty() {
        let button = ui.add_enabled(!pending.is_empty(), egui::Button::new(format!("Defeature the {} proposed", pending.len())));
        if button.on_hover_text(format!("Take them out with defeature {name}")).clicked() {
            acts.push(Act::Propose { name: name.clone(), features: pending });
        }
        ui.add_space(4.0);
    }
    egui::Grid::new("fc_defeature").num_columns(4).spacing(vec2(10.0, 2.0)).striped(true).show(ui, |ui| {
        for h in ["In mesh", "Feature", "Kind", "Size"] {
            ui.label(RichText::new(h).color(t.text_dim));
        }
        ui.end_row();
        for c in candidates.iter().take(MAX_ENTITIES) {
            let Some(id) = c["id"].as_u64() else { continue };
            let mut on = c["suppressed"] != true;
            let by = c["by"].as_str().unwrap_or(&name).to_string();
            let hover = if on { format!("Switch off: defeature {name} takes it out") } else { format!("Switch on: defeature {by} keeps it") };
            if ui.checkbox(&mut on, "").on_hover_text(hover).changed() {
                acts.push(Act::Defeature { name: if on { by } else { name.clone() }, feature: id, on });
            }
            let label = format!("F{id} {}", c["name"].as_str().unwrap_or_default());
            let text = if c["proposed"] == true { RichText::new(label).color(t.warm) } else { RichText::new(label) };
            match c["faces"].as_array().and_then(|f| f.first()).and_then(Value::as_str).filter(|_| on) {
                Some(face) => {
                    if ui.link(text).on_hover_text("Pick its face in the viewport").clicked() {
                        let entity = fieldcraft_model::Entity::Face(face.to_string());
                        acts.push(Act::Pick(fieldcraft_model::GeometryRef { body: body.to_string(), entity }));
                    }
                }
                None => {
                    ui.label(text);
                }
            }
            let method = match c["method"].as_str() {
                Some("virtual") => " (virtual)",
                _ => "",
            };
            ui.label(format!("{}{method}", c["category"].as_str().unwrap_or_default()));
            ui.label(format!("{} {}", crate::describe::quantity(sig(c["size"].as_f64().unwrap_or(0.0))), geometry_unit(units, Q::Length)));
            ui.end_row();
        }
    });
    more_than(ui, t, candidates.len(), MAX_ENTITIES);
}

/// The faces of `body` named `faces`, as selection items.
fn face_items(body: &str, faces: &[String]) -> Vec<crate::selection::Item> {
    faces
        .iter()
        .take(MAX_ENTITIES)
        .filter_map(|f| fieldcraft_model::GeometryRef::new(body, fieldcraft_model::Entity::Face(f.clone())).ok())
        .map(crate::selection::Item::Geometry)
        .collect()
}

/// The Feature Graph: per body, its timeline's features in order, each with its glyph, name and
/// what it does in words; a click selects the faces it made and opens its details. Returns the
/// faces of the row under the pointer, to highlight.
fn feature_graph(
    ui: &mut egui::Ui,
    t: &Tokens,
    v: &Value,
    opened: &std::collections::BTreeSet<(String, u64)>,
    acts: &mut Vec<Act>,
) -> Vec<crate::selection::Item> {
    use crate::feature_graph::{describe, paint_glyph};
    let mut hovered = Vec::new();
    let bodies = v["bodies"].as_array().map(Vec::as_slice).unwrap_or_default();
    for (bi, b) in bodies.iter().enumerate() {
        let body = b["body"].as_str().unwrap_or_default();
        let nodes = b["features"].as_array().map(Vec::as_slice).unwrap_or_default();
        let head = RichText::new(format!("{body}  ·  {}", plural(nodes.len() as u64, "feature"))).strong();
        egui::CollapsingHeader::new(head).id_salt(("fc_features_body", bi)).default_open(true).show(ui, |ui| {
            for (fi, node) in nodes.iter().take(MAX_ENTITIES).enumerate() {
                let f = describe(node, v["units"]["length"].as_str().unwrap_or_default());
                let id = ui.id().with(("fc_feature", bi, fi));
                let key = (body.to_string(), f.id);
                let open = opened.contains(&key);
                let bg = ui.painter().add(egui::Shape::Noop);
                let inner = egui::Frame::NONE.inner_margin(egui::Margin::symmetric(4, 4)).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.horizontal_top(|ui| {
                        let (rect, _) = ui.allocate_exact_size(vec2(18.0, 18.0), egui::Sense::hover());
                        let accent = if f.error.is_some() { t.error } else { t.accent };
                        paint_glyph(ui.painter(), rect, f.kind, t.icon, accent);
                        ui.vertical(|ui| {
                            ui.horizontal(|ui| {
                                ui.add(egui::Label::new(RichText::new(&f.name).strong().color(t.text)).truncate());
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    ui.label(RichText::new(format!("F{}", f.id)).small().color(t.text_dim));
                                });
                            });
                            ui.add(egui::Label::new(RichText::new(&f.summary).color(t.text)).wrap());
                            if let Some(e) = &f.error {
                                ui.add(egui::Label::new(RichText::new(e).color(t.error)).wrap());
                            }
                            if open {
                                ui.add_space(2.0);
                                ui.spacing_mut().item_spacing.y = 0.0;
                                ui.spacing_mut().interact_size.y = 0.0;
                                for (k, val) in &f.details {
                                    ui.horizontal_wrapped(|ui| {
                                        ui.spacing_mut().item_spacing.x = 6.0;
                                        ui.label(RichText::new(k).small().color(t.text_dim));
                                        ui.label(RichText::new(val).small());
                                    });
                                }
                            }
                        });
                    });
                });
                let rect = inner.response.rect;
                let resp = ui.interact(rect, id.with("row"), egui::Sense::click());
                let hover_text = if f.faces.is_empty() { "Show its details" } else { "Select the faces it made; show its details" };
                let resp = resp.on_hover_text(hover_text);
                if resp.hovered() {
                    ui.painter().set(bg, egui::Shape::rect_filled(rect, 3.0, t.hover));
                    hovered = face_items(body, &f.faces);
                } else if open {
                    ui.painter().set(bg, egui::Shape::rect_filled(rect, 3.0, t.subheading_bg));
                }
                if resp.clicked() {
                    acts.push(Act::Open(key));
                    if !f.faces.is_empty() {
                        acts.push(Act::Faces(face_items(body, &f.faces)));
                    }
                }
            }
            more_than(ui, t, nodes.len(), MAX_ENTITIES);
        });
    }
    hovered
}

/// The definition a diagnostic's subject names (``section `web` ``), when the model has it.
fn subject_definition(app: &ModelApp, subject: &str) -> Option<(Def, String)> {
    let (kind, rest) = subject.split_once(" `")?;
    let name = rest.strip_suffix('`')?;
    let kind: Def = serde_json::from_value(Value::String(kind.replace(' ', "_"))).ok()?;
    browser::definition(&app.session.model, kind, name).map(|_| (kind, name.to_string()))
}

/// Where a whole-model diagnostic is fixed: a button's text, tooltip and what it does.
fn model_fix(message: &str) -> Option<(&'static str, Act)> {
    if message.contains("mesh has no elements") {
        Some(("Mesh…", Act::Tab(Tab::Mesh)))
    } else if message.contains("no unit system") {
        Some(("Units…", Act::Start("model.units")))
    } else if message.contains("no analysis step") {
        Some(("Step…", Act::Start("model.step.create")))
    } else if message.contains("have no section") {
        Some(("Section…", Act::Start("model.section.create")))
    } else {
        None
    }
}

/// Check Model: each diagnostic with its severity, the item it is about and its message, and a
/// button to where it is fixed.
fn diagnostics(app: &ModelApp, ui: &mut egui::Ui, t: &Tokens, v: &Value, acts: &mut Vec<Act>) {
    let list = v["diagnostics"].as_array().map(Vec::as_slice).unwrap_or_default();
    if list.is_empty() {
        ui.label(RichText::new("✔ The model is ready to solve.").color(t.ok));
        return;
    }
    for (i, d) in list.iter().enumerate() {
        let error = d["severity"] == "error";
        let subject = d["subject"].as_str().unwrap_or_default();
        let message = d["message"].as_str().unwrap_or_default();
        ui.horizontal_top(|ui| {
            let (mark, colour, tip) = if error { ("⛔", t.error, "Error: blocks solving") } else { ("⚠", t.warm, "Warning") };
            ui.label(RichText::new(mark).color(colour)).on_hover_text(tip);
            ui.vertical(|ui| {
                ui.set_width(ui.available_width() - 70.0);
                ui.label(RichText::new(subject).strong());
                ui.add(egui::Label::new(message).wrap());
            });
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Min), |ui| {
                if let Some((kind, name)) = subject_definition(app, subject) {
                    if ui.small_button("Show").on_hover_text(format!("Select {subject} to see and edit it")).clicked() {
                        acts.push(Act::Select(kind, name));
                    }
                } else if subject == "model"
                    && let Some((text, act)) = model_fix(message)
                    && ui.small_button(text).clicked()
                {
                    acts.push(act);
                }
            });
        });
        if i + 1 < list.len() {
            ui.separator();
        }
    }
}

/// Faces or edges listed under a body at most.
const MAX_ENTITIES: usize = 500;

/// Geometry › Bodies: one entry per body, its name first, then volume, area, face and edge counts
/// and its size; a click on the name picks the body, and its faces and edges unfold, by
/// persistent name (edges with their lengths), each picked by a click: the names the seed and set
/// forms ask for. Design geometry is in millimetres.
/// Mesh Quality's metrics: key, label, and whether a low value is the bad end.
const METRICS: [(&str, &str, bool); 4] = [
    ("scaled_jacobian", "Scaled Jacobian", true),
    ("jacobian_ratio", "Jacobian ratio", true),
    ("shape", "Shape", true),
    ("edge_ratio", "Edge ratio", false),
];

/// Elements a Select button selects at most.
const MAX_SELECT: usize = 10_000;

/// `mesh.quality`'s reply: the counts, then each metric's range, histogram and worst element (a
/// click selects and frames it), buttons selecting the invalid and undetermined elements, and
/// Relax around worst element (#398).
fn quality(ui: &mut egui::Ui, t: &Tokens, v: &Value, acts: &mut Vec<Act>) {
    let n = |k: &str| v.get(k).and_then(Value::as_u64).unwrap_or(0);
    let labels = |k: &str| -> Vec<u32> {
        v.pointer(&format!("/{k}/labels"))
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|l| l.as_u64()?.try_into().ok())
            .take(MAX_SELECT)
            .collect()
    };
    let (invalid, undetermined) =
        (v.pointer("/invalid/count").and_then(Value::as_u64).unwrap_or(0), v.pointer("/undetermined/count").and_then(Value::as_u64).unwrap_or(0));
    egui::Grid::new("fc_quality_counts").num_columns(2).spacing(vec2(10.0, 2.0)).show(ui, |ui| {
        let mut row = |k: &str, text: RichText| {
            ui.label(RichText::new(k).color(t.text_dim));
            ui.label(text);
            ui.end_row();
        };
        row("Elements", RichText::new(count(n("measured"))));
        row("Valid", RichText::new(count(n("valid"))));
        row("Invalid", if invalid > 0 { RichText::new(count(invalid)).color(t.error) } else { RichText::new("0") });
        if undetermined > 0 {
            row("Undetermined", RichText::new(count(undetermined)).color(t.warm));
        }
        if let Some(vol) = v.get("volume").and_then(Value::as_f64).filter(|x| *x > 0.0) {
            row("Volume", RichText::new(number(vol)));
        }
        // A surface mesh's area; a solid mesh has none to show.
        if let Some(area) = v.get("area").and_then(Value::as_f64).filter(|x| *x > 0.0) {
            row("Area", RichText::new(number(area)));
        }
    });
    for (key, what) in [("invalid", "invalid"), ("undetermined", "undetermined")] {
        let l = labels(key);
        if !l.is_empty()
            && ui
                .button(format!("Select the {} {what} elements", count(l.len() as u64)))
                .on_hover_text("Select and frame them in the viewport")
                .clicked()
        {
            acts.push(Act::Elements(l));
        }
    }
    if let Some(worst) = v.pointer("/metrics/scaled_jacobian/worst").and_then(Value::as_u64).and_then(|w| u32::try_from(w).ok())
        && ui
            .button("Relax around worst element")
            .on_hover_text("Open Relax Nodes with its sphere round the element of smallest scaled Jacobian; preview, then OK")
            .clicked()
    {
        acts.push(Act::Relax(worst));
    }
    for (key, label, low_bad) in METRICS {
        let Some(m) = v.pointer(&format!("/metrics/{key}")).filter(|m| m.is_object()) else { continue };
        let f = |k: &str| m.get(k).and_then(Value::as_f64);
        ui.add_space(6.0);
        ui.label(RichText::new(label).strong());
        ui.label(
            RichText::new(format!(
                "min {} · mean {} · max {}",
                f("min").map_or("—".into(), number),
                f("mean").map_or("—".into(), number),
                f("max").map_or("—".into(), number)
            ))
            .color(t.text_dim),
        );
        let counts: Vec<u64> = m.pointer("/histogram/counts").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_u64).collect();
        let range: Vec<f64> = m.pointer("/histogram/range").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_f64).collect();
        if let ([lo, hi], false) = (range.as_slice(), counts.is_empty()) {
            histogram(ui, t, &counts, (*lo, *hi), low_bad);
        }
        if let Some(worst) = m.get("worst").and_then(Value::as_u64).and_then(|w| u32::try_from(w).ok()) {
            let value = if low_bad { f("min") } else { f("max") };
            let text = format!("Worst: element {worst}{}", value.map(|x| format!(" ({})", number(x))).unwrap_or_default());
            if ui.link(text).on_hover_text("Select and frame it in the viewport").clicked() {
                acts.push(Act::Elements(vec![worst]));
            }
        }
    }
    if let Some(skipped) = v.get("skipped").and_then(Value::as_object).filter(|s| !s.is_empty()) {
        ui.add_space(6.0);
        let text: Vec<String> = skipped.iter().map(|(k, n)| format!("{k}: {}", n.as_u64().map_or_else(|| n.to_string(), count))).collect();
        ui.label(RichText::new(format!("Not measured: {}", text.join(", "))).color(t.text_dim));
    }
}

/// A histogram of `counts` over `range`, its bars shaded from the bad end (warm) to the good
/// one (accent); the range's ends under it.
fn histogram(ui: &mut egui::Ui, t: &Tokens, counts: &[u64], range: (f64, f64), low_bad: bool) {
    let width = ui.available_width().clamp(120.0, 320.0);
    let (rect, resp) = ui.allocate_exact_size(vec2(width, 48.0), egui::Sense::hover());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 2.0, t.panel_header);
    let top = counts.iter().copied().max().unwrap_or(0).max(1) as f32;
    let w = rect.width() / counts.len() as f32;
    let bins = counts.len().max(1) as f32;
    for (i, c) in counts.iter().enumerate() {
        let h = (*c as f32 / top).sqrt() * (rect.height() - 2.0);
        let x = rect.left() + i as f32 * w;
        let good = if low_bad { (i as f32 + 0.5) / bins } else { 1.0 - (i as f32 + 0.5) / bins };
        let color = if good < 0.3 { t.warm } else { t.accent };
        if *c > 0 {
            p.rect_filled(
                egui::Rect::from_min_max(egui::pos2(x + 1.0, rect.bottom() - h.max(1.0)), egui::pos2(x + w - 1.0, rect.bottom())),
                0.0,
                color,
            );
        }
    }
    let step = (range.1 - range.0) / bins as f64;
    resp.on_hover_ui(|ui| {
        for (i, c) in counts.iter().enumerate() {
            let a = range.0 + step * i as f64;
            ui.label(format!("{} – {}: {}", number(a), number(a + step), count(*c)));
        }
    });
    ui.horizontal(|ui| {
        ui.label(RichText::new(number(range.0)).small().color(t.text_dim));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| ui.label(RichText::new(number(range.1)).small().color(t.text_dim)));
    });
}

fn bodies(ui: &mut egui::Ui, t: &Tokens, v: &Value, units: &DisplaySystem, acts: &mut Vec<Act>) {
    use fieldcraft_model::units::Quantity as Q;
    let (length, area, volume) = (geometry_unit(units, Q::Length), geometry_unit(units, Q::Area), geometry_unit(units, Q::Volume));
    let q = crate::describe::quantity;
    let pick = |body: &str, entity: fieldcraft_model::Entity| fieldcraft_model::GeometryRef { body: body.to_string(), entity };
    for (i, b) in v.as_array().into_iter().flatten().enumerate() {
        let name = b["name"].as_str().unwrap_or_default();
        let faces = b["faces"].as_array().map(Vec::as_slice).unwrap_or_default();
        let edges = b["edges"].as_array().map(Vec::as_slice).unwrap_or_default();
        let corner = |k: usize| -> Option<[f64; 3]> { serde_json::from_value(b["bbox"][k].clone()).ok() };
        let size = corner(0)
            .zip(corner(1))
            .map(|(lo, hi)| format!("{} × {} × {} {length}", q(sig(hi[0] - lo[0])), q(sig(hi[1] - lo[1])), q(sig(hi[2] - lo[2]))));
        let id = egui::Id::new(("fc_reply_body", i));
        if ui.link(RichText::new(name).strong()).on_hover_text("Pick this body in the viewport").clicked() {
            acts.push(Act::Pick(pick(name, fieldcraft_model::Entity::Body)));
        }
        egui::Grid::new(id.with("grid")).num_columns(2).spacing(vec2(10.0, 2.0)).show(ui, |ui| {
            let mut row = |k: &str, text: String| {
                ui.label(RichText::new(k).color(t.text_dim));
                ui.label(text);
                ui.end_row();
            };
            if let Some(v) = b["volume"].as_f64() {
                row("Volume", format!("{} {volume}", q(sig(v))));
            }
            if let Some(a) = b["area"].as_f64() {
                row("Area", format!("{} {area}", q(sig(a))));
            }
            row("Faces, edges", format!("{}, {}", faces.len(), edges.len()));
            if let Some(size) = size {
                row("Size", size);
            }
        });
        body_list(ui, t, id.with("faces"), "Faces", faces, &|f| (f["name"].as_str().unwrap_or_default().to_string(), None), name, length, acts);
        body_list(
            ui,
            t,
            id.with("edges"),
            "Edges",
            edges,
            &|e| (e["name"].as_str().unwrap_or_default().to_string(), e["length"].as_f64()),
            name,
            length,
            acts,
        );
        ui.separator();
    }
}

/// A body's faces or edges, folded: each name picks it in the viewport (edges with their length).
#[allow(clippy::too_many_arguments)]
fn body_list(
    ui: &mut egui::Ui,
    t: &Tokens,
    id: egui::Id,
    title: &str,
    items: &[Value],
    read: &dyn Fn(&Value) -> (String, Option<f64>),
    body: &str,
    length_unit: &str,
    acts: &mut Vec<Act>,
) {
    if items.is_empty() {
        return;
    }
    egui::CollapsingHeader::new(RichText::new(format!("{title} ({})", items.len())).size(12.0)).id_salt(id).show(ui, |ui| {
        for item in items.iter().take(MAX_ENTITIES) {
            let (name, length) = read(item);
            ui.horizontal(|ui| {
                if ui.link(&name).on_hover_text("Pick it in the viewport").clicked() {
                    let entity =
                        if title == "Edges" { fieldcraft_model::Entity::Edge(name.clone()) } else { fieldcraft_model::Entity::Face(name.clone()) };
                    acts.push(Act::Pick(fieldcraft_model::GeometryRef { body: body.to_string(), entity }));
                }
                if let Some(l) = length {
                    ui.label(RichText::new(format!("{} {length_unit}", crate::describe::quantity(sig(l)))).color(t.text_dim));
                }
            });
        }
        more_than(ui, t, items.len(), MAX_ENTITIES);
    });
}

/// `v` to 3 significant figures.
fn sig(v: f64) -> f64 {
    if v == 0.0 || !v.is_finite() {
        return v;
    }
    let scale = 10f64.powi(2 - v.abs().log10().floor() as i32);
    (v * scale).round() / scale
}

fn more_than(ui: &mut egui::Ui, t: &Tokens, n: usize, max: usize) {
    if n > max {
        ui.label(RichText::new(format!("… {} more", n - max)).small().color(t.text_dim));
    }
}

/// The Step Table: each boundary condition and load against the steps, ticked where active.
fn activity(ui: &mut egui::Ui, t: &Tokens, v: &Value) {
    let steps: Vec<&str> = v["steps"].as_array().into_iter().flatten().filter_map(|s| s["name"].as_str()).collect();
    if steps.is_empty() {
        ui.label(RichText::new("No analysis steps yet: Steps › Step creates one.").color(t.text_dim));
        return;
    }
    egui::ScrollArea::horizontal().id_salt("fc_reply_activity").show(ui, |ui| {
        egui::Grid::new("fc_reply_activity_grid").striped(true).spacing(vec2(10.0, 3.0)).show(ui, |ui| {
            ui.label("");
            for (s, step) in steps.iter().zip(v["steps"].as_array().into_iter().flatten()) {
                ui.label(RichText::new(*s).strong()).on_hover_text(step["procedure"].as_str().unwrap_or_default());
            }
            ui.end_row();
            for (key, kind) in [("boundary_conditions", "BC"), ("loads", "Load")] {
                for (name, active) in v[key].as_object().into_iter().flatten() {
                    ui.label(name.as_str()).on_hover_text(kind);
                    for s in &steps {
                        let on = active.as_array().is_some_and(|a| a.iter().any(|x| x == *s));
                        ui.label(if on { RichText::new("✔").color(t.ok) } else { RichText::new("—").color(t.text_dim) });
                    }
                    ui.end_row();
                }
            }
        });
    });
}

/// A solved step's peaks in one line: `max |U| 0.05063 mm, max Mises 1,368 MPa, 115,665 DOFs`.
fn solve_line(v: &Value, units: Option<&DisplaySystem>) -> String {
    let mut parts = Vec::new();
    let peak = |key: &str| v.get(key).and_then(|p| p.get("value")).and_then(Value::as_f64);
    for (key, what, q) in [
        ("max_displacement", "max |U|", Quantity::Length),
        ("max_von_mises", "max Mises", Quantity::Pressure),
        ("max_temperature", "max T", Quantity::Temperature),
    ] {
        if let Some(x) = peak(key) {
            parts.push(format!("{what} {}", with_unit(number(x), unit(q, units).as_deref())));
        }
    }
    if let Some(m) = v["modes"].as_array().and_then(|m| m.first()) {
        if let Some(f) = m["frequency"].as_f64() {
            parts.push(format!("mode 1 at {} Hz", number(f)));
        } else if let Some(l) = m["load_factor"].as_f64() {
            parts.push(format!("first load factor {}", number(l)));
        }
    }
    if let Some(n) = v["dofs"].as_u64() {
        parts.push(format!("{} DOFs", count(n)));
    }
    parts.join(", ")
}

/// What a solve summary's row shows: a number of a quantity, a node's peak, a vector, a count,
/// or text.
enum Row {
    Number(Quantity),
    Peak(Quantity),
    Vector(Quantity),
    Count,
    Text,
}

/// A solve's summary: the step's peaks and totals with their units first, round-off as 0, rows
/// that do not apply left out (a steady step's end time), counts with thousands separators;
/// lists (modes, contacts) as tables, and the solver's internals folded under Details.
fn solve_summary(ui: &mut egui::Ui, t: &Tokens, v: &Value, units: Option<&DisplaySystem>) {
    let Some(m) = v.as_object() else { return value(ui, t, egui::Id::new("fc_reply_solve"), v, 0) };
    let steady = m.get("frames").and_then(Value::as_u64).is_some_and(|n| n <= 1);
    let rows: [(&str, &str, Row); 18] = [
        ("step", "Step", Row::Text),
        ("max_displacement", "Max |U|", Row::Peak(Quantity::Length)),
        ("max_von_mises", "Max Mises", Row::Peak(Quantity::Pressure)),
        ("max_temperature", "Max T", Row::Peak(Quantity::Temperature)),
        ("min_temperature", "Min T", Row::Peak(Quantity::Temperature)),
        ("max_plastic_strain", "Max plastic strain", Row::Peak(Quantity::Dimensionless)),
        ("applied_total", "Applied force", Row::Vector(Quantity::Force)),
        ("reaction_total", "Reaction force", Row::Vector(Quantity::Force)),
        ("heat_supplied", "Heat flow through", Row::Number(Quantity::Power)),
        ("total_mass", "Total mass", Row::Number(Quantity::Mass)),
        ("load_factor", "Load factor", Row::Number(Quantity::Dimensionless)),
        ("time", "Time reached", Row::Number(Quantity::Time)),
        ("end_time", "End time", Row::Number(Quantity::Time)),
        ("completed", "Completed", Row::Text),
        ("energy_balance", "Energy balance error", Row::Number(Quantity::Dimensionless)),
        ("dofs", "Degrees of freedom", Row::Count),
        ("elements", "Elements", Row::Count),
        ("frames", "Frames", Row::Count),
    ];
    let mut shown: Vec<&str> = Vec::new();
    egui::Grid::new("fc_reply_solve").num_columns(2).spacing(vec2(10.0, 3.0)).show(ui, |ui| {
        for (key, label, row) in &rows {
            let Some(x) = m.get(*key).filter(|x| !x.is_null()) else { continue };
            shown.push(key);
            // A steady step has no time, and one frame.
            if steady && matches!(*key, "end_time" | "frames") {
                continue;
            }
            let text = match row {
                Row::Number(q) => x.as_f64().map(|f| with_unit(number(f), unit(*q, units).as_deref())),
                Row::Peak(q) => x["value"].as_f64().map(|f| {
                    let at = x["node"].as_u64().map_or_else(String::new, |n| format!(" at node {n}"));
                    format!("{}{at}", with_unit(number(f), unit(*q, units).as_deref()))
                }),
                Row::Vector(q) => {
                    serde_json::from_value::<Vec<f64>>(x.clone()).ok().map(|c| with_unit(components(&c, 0.0), unit(*q, units).as_deref()))
                }
                Row::Count => x.as_u64().map(count),
                Row::Text => inline(x),
            };
            if let Some(text) = text {
                ui.label(RichText::new(*label).color(t.text_dim));
                ui.add(egui::Label::new(text).wrap());
                ui.end_row();
            }
        }
    });
    // Heat removed is the supplied flow with the other sign: one row says it.
    shown.push("heat_removed");
    if let Some(modes) = m.get("modes").and_then(Value::as_array).filter(|a| !a.is_empty()) {
        shown.push("modes");
        egui::CollapsingHeader::new(RichText::new(format!("Modes ({})", modes.len())).size(12.0))
            .id_salt("fc_reply_solve_modes")
            .default_open(true)
            .show(ui, |ui| {
                egui::Grid::new("fc_reply_solve_modes_grid").striped(true).spacing(vec2(10.0, 2.0)).show(ui, |ui| {
                    let buckling = modes.iter().any(|x| x.get("load_factor").is_some());
                    ui.label(RichText::new("Mode").small().color(t.text_dim));
                    ui.label(RichText::new(if buckling { "Load factor" } else { "Frequency [Hz]" }).small().color(t.text_dim));
                    ui.end_row();
                    for x in modes.iter().take(MAX_ROWS) {
                        ui.label(x["mode"].as_u64().map_or_else(String::new, |n| n.to_string()));
                        ui.label(x.get(if buckling { "load_factor" } else { "frequency" }).and_then(Value::as_f64).map_or("—".into(), number));
                        ui.end_row();
                    }
                });
                more(ui, t, modes.len());
            });
    }
    if let Some(w) = m.get("warnings").and_then(Value::as_array) {
        shown.push("warnings");
        for x in w {
            ui.add(egui::Label::new(RichText::new(format!("⚠ {}", inline(x).unwrap_or_else(|| crate::context::value_text(x)))).color(t.warm)).wrap());
        }
    }
    // The rest: the solver's internals (equations, factor entries, Lanczos vectors…) and lists.
    let rest: serde_json::Map<String, Value> = m
        .iter()
        .filter(|(k, x)| !shown.contains(&k.as_str()) && !x.is_null() && x.as_array().is_none_or(|a| !a.is_empty()))
        .map(|(k, x)| (k.clone(), x.clone()))
        .collect();
    if !rest.is_empty() {
        egui::CollapsingHeader::new(RichText::new("Details").size(12.0)).id_salt("fc_reply_solve_details").default_open(false).show(ui, |ui| {
            let rest = Value::Object(rest.into_iter().map(|(k, x)| (k, x.as_u64().map_or(x.clone(), |n| Value::String(count(n))))).collect());
            value(ui, t, egui::Id::new("fc_reply_solve_details"), &rest, 1);
        });
    }
}

/// Any reply: objects as labelled rows, lists of objects as tables, nested parts folded.
fn value(ui: &mut egui::Ui, t: &Tokens, id: egui::Id, v: &Value, depth: usize) {
    match v {
        Value::Object(m) => {
            egui::Grid::new(id).num_columns(2).spacing(vec2(10.0, 3.0)).show(ui, |ui| {
                for (k, x) in m.iter().filter(|(k, x)| k.as_str() != "summary" && inline(x).is_some()) {
                    ui.label(RichText::new(k.replace('_', " ")).color(t.text_dim));
                    ui.add(egui::Label::new(inline(x).unwrap_or_default()).wrap());
                    ui.end_row();
                }
            });
            for (k, x) in m.iter().filter(|(_, x)| inline(x).is_none()) {
                let title = match x {
                    Value::Array(a) => format!("{} ({})", k.replace('_', " "), a.len()),
                    _ => k.replace('_', " "),
                };
                egui::CollapsingHeader::new(RichText::new(title).size(12.0)).id_salt(id.with(k)).default_open(depth == 0).show(ui, |ui| {
                    value(ui, t, id.with(k), x, depth + 1);
                });
            }
        }
        Value::Array(a) if depth < MAX_DEPTH && a.iter().all(Value::is_object) => table(ui, t, id, a),
        Value::Array(a) if depth < MAX_DEPTH => {
            for (i, x) in a.iter().take(MAX_ROWS).enumerate() {
                value(ui, t, id.with(i), x, depth + 1);
            }
            more(ui, t, a.len());
        }
        other => {
            ui.add(egui::Label::new(inline(other).unwrap_or_else(|| crate::context::value_text(other))).wrap());
        }
    }
}

/// A value short enough for one cell: scalars, and lists of up to four scalars.
fn inline(v: &Value) -> Option<String> {
    match v {
        Value::Null => Some("—".into()),
        Value::Array(a) if a.len() <= 4 && a.iter().all(|x| !x.is_object() && !x.is_array()) => {
            Some(a.iter().filter_map(scalar).collect::<Vec<_>>().join(", "))
        }
        Value::Array(_) | Value::Object(_) => None,
        other => scalar(other),
    }
}

/// A list of objects as a table: the columns are their keys (first seen first).
fn table(ui: &mut egui::Ui, t: &Tokens, id: egui::Id, rows: &[Value]) {
    let mut columns: Vec<&str> = Vec::new();
    for r in rows.iter().take(MAX_ROWS) {
        for k in r.as_object().into_iter().flat_map(Map::keys) {
            if !columns.contains(&k.as_str()) && columns.len() < MAX_COLUMNS {
                columns.push(k);
            }
        }
    }
    egui::ScrollArea::horizontal().id_salt(id.with("scroll")).show(ui, |ui| {
        egui::Grid::new(id).striped(true).spacing(vec2(8.0, 2.0)).show(ui, |ui| {
            for c in &columns {
                ui.label(RichText::new(c.replace('_', " ")).small().color(t.text_dim));
            }
            ui.end_row();
            for r in rows.iter().take(MAX_ROWS) {
                for c in &columns {
                    ui.add(
                        egui::Label::new(r.get(*c).map_or_else(String::new, |x| inline(x).unwrap_or_else(|| crate::context::value_text(x)))).extend(),
                    );
                }
                ui.end_row();
            }
        });
    });
    more(ui, t, rows.len());
}

fn more(ui: &mut egui::Ui, t: &Tokens, n: usize) {
    if n > MAX_ROWS {
        ui.label(RichText::new(format!("… {} more", n - MAX_ROWS)).small().color(t.text_dim));
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    /// Propose Defeaturing stays in the panel and follows the model: switching the plate's hole
    /// off fills it in the viewport's bodies and the list shows it off; on again brings it back.
    #[test]
    fn the_feature_graph_reads_a_plate_in_words_and_selects_a_features_faces() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-features-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        app.run("geometry.features", json!({})).unwrap();
        let reply = app.reply.clone().expect("the graph is shown");
        assert!(summary("geometry.features", &reply.value).contains("building 1 body"));
        let body = reply.value["bodies"][0]["body"].as_str().unwrap().to_string();
        let features: Vec<_> =
            reply.value["bodies"][0]["features"].as_array().unwrap().iter().map(|f| crate::feature_graph::describe(f, "mm")).collect();
        assert!(!features.is_empty());
        for f in &features {
            assert!(!f.summary.contains('{') && !f.summary.is_empty(), "{f:?}");
        }
        let last = features.iter().rev().find(|f| !f.faces.is_empty()).expect("a feature made faces");
        let items = face_items(&body, &last.faces);
        assert!(!items.is_empty());
        app.select(items.clone(), crate::selection::Mode::Replace);
        assert_eq!(app.pick.selection.len(), items.len());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_defeature_list_switches_a_hole_off_and_on() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-defeature-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        let body = app.session.design.as_ref().unwrap().bodies[0].name.clone();
        app.run("geometry.defeature.propose", json!({"body": body, "below": {"hole_diameter": 12}})).unwrap();
        let hole = app.reply.as_ref().unwrap().value["proposal"][0].as_u64().unwrap();
        let faces = |app: &ModelApp| app.session.design.as_ref().unwrap().bodies[0].faces().len();
        set_defeature(&mut app, "Defeature", vec![hole]).unwrap();
        let reply = app.reply.as_ref().expect("the list stays open");
        assert_eq!((faces(&app), reply.value["candidates"][0]["suppressed"].as_bool()), (6, Some(true)));
        assert_eq!(defeatured(&app, "Defeature"), [hole]);
        set_defeature(&mut app, "Defeature", Vec::new()).unwrap();
        assert_eq!((faces(&app), app.reply.as_ref().unwrap().value["candidates"][0]["suppressed"].as_bool()), (8, Some(false)));
        assert!(app.session.model.defeatures().is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A live reply runs again in the units it was asked in: a defeature proposal below half an
    /// inch still lists the plate's hole after the display switches to millimetres.
    #[test]
    fn a_live_reply_keeps_its_units_when_the_display_changes() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-live-units-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        app.run("model.units", json!({"system": "in-lbf-s-psi"})).unwrap();
        let body = app.session.design.as_ref().unwrap().bodies[0].name.clone();
        app.run("geometry.defeature.propose", json!({"body": body, "below": {"hole_diameter": 0.5}})).unwrap();
        let proposal = |app: &ModelApp| app.reply.as_ref().map(|r| r.value["proposal"].clone());
        let before = proposal(&app);
        assert_eq!(before.as_ref().and_then(|p| p.as_array()).map(Vec::len), Some(1), "{before:?}");
        app.run("model.units", json!({"system": "mm-t-s-N-MPa"})).unwrap();
        assert_eq!(proposal(&app), before);
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Mesh Quality's reply draws as counts, histograms and worst elements, and a worst element
    /// is selected and framed (#113).
    #[test]
    fn mesh_quality_shows_its_metrics_and_frames_the_worst_element() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("mesh.quality", json!({})).unwrap();
        assert_eq!(app.reply.as_ref().map(|r| r.command.as_str()), Some("mesh.quality"));
        let t = Tokens::get(false);
        let ctx = egui::Context::default();
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| show(&mut app, ui, &t));
        out.textures_delta.clear();
        let worst = app.reply.as_ref().and_then(|r| r.value.pointer("/metrics/scaled_jacobian/worst")?.as_u64()).unwrap() as u32;
        let before = app.cam;
        app.show_elements(&[worst]);
        assert_eq!(app.pick.selection.summary(), "1 element");
        assert_ne!(app.cam, before, "framed");
    }

    #[test]
    fn summaries_are_sentences() {
        let check = json!({"ready": false, "errors": 2, "warnings": 1, "diagnostics": []});
        assert_eq!(summary("model.check", &check), "2 errors, 1 warning");
        assert_eq!(summary("model.check", &json!({"errors": 0, "warnings": 0})), "ready to solve, no issues");
        let activity = json!({"boundary_conditions": {"Clamp": ["Load"]}, "loads": {}, "steps": [{"name": "Load"}]});
        assert_eq!(summary("model.activity", &activity), "1 step, 1 boundary condition, 0 loads");
        assert_eq!(summary("geometry.bodies", &json!([])), "no bodies (open a design first)");
        assert_eq!(summary("geometry.bodies", &json!([{}, {}])), "2 bodies");
        assert_eq!(summary("mesh.info", &json!({"dofs": 945, "modes": [1, 2], "x": {"y": 1}})), "dofs 945, modes 2");
        let solved = json!({"dofs": 115665, "max_displacement": {"node": 17, "value": 0.0506333}, "max_von_mises": {"node": 3, "value": 1368.07}});
        assert_eq!(
            summary_in("solve.static", &solved, Some(&DisplaySystem::new(fieldcraft_model::units::Preset::MmTonneS))),
            "max |U| 0.05063 mm, max Mises 1,368 MPa, 115,665 DOFs"
        );
        assert_eq!(summary("any", &json!({"summary": "a\nb"})), "a b");
        assert_eq!(summary("any", &json!({})), "done");
        assert!(summary("any", &json!({"s": "x".repeat(500)})).chars().count() <= MAX_LINE + 1);
    }

    #[test]
    fn model_level_diagnostics_lead_to_their_fix() {
        let mut s = fieldcraft_engine::Session::default();
        let check = s.execute("model.check", &json!({})).unwrap();
        let messages: Vec<&str> = check["diagnostics"].as_array().unwrap().iter().filter_map(|d| d["message"].as_str()).collect();
        assert!(messages.len() >= 2, "{check}");
        for m in messages {
            assert!(model_fix(m).is_some(), "no fix for `{m}`");
        }
    }
}
