//! A command's form drawn from its parameter schema ([`fieldcraft_engine::schema`]): typed fields
//! with units, drop-downs of the model's items, groups and one-of variants. The form edits the
//! command's JSON parameters in place; it never builds anything the command does not take.

use std::collections::{BTreeMap, BTreeSet};

use egui::{Id, RichText, Ui};
use fieldcraft_engine::schema::{self, Collection, Field, Kind, MAX_DEPTH, MAX_ROWS, join};
use fieldcraft_model::DisplaySystem;
use serde_json::{Map, Value, json};

use crate::fields::{self, NumberEdit, Style};
use crate::theme::Tokens;

/// The label column: as wide as the form's longest label, at least this…
const LABEL_MIN: f32 = 72.0;
/// …and at most this share of the form's width; longer labels wrap.
const LABEL_SHARE: f32 = 0.45;
/// The narrowest field.
const FIELD_MIN: f32 = 48.0;
/// Indentation of a group's fields.
const INDENT: f32 = 10.0;

/// What the form reads beside the parameters.
pub struct Context<'a> {
    /// The display system: the numbers' units.
    pub units: &'a DisplaySystem,
    /// The names of a collection's items.
    pub names: &'a dyn Fn(Collection) -> Vec<String>,
    /// A salt for the widget ids: the open form (its command and the opening), so no widget
    /// state is shared with another form.
    pub salt: Id,
    /// The fields the user has edited: their problems show beside them; the others wait.
    pub touched: &'a BTreeSet<String>,
    /// OK was tried: every problem shows.
    pub attempted: bool,
    /// The command refused a parameter: its path and the problem as the form phrases it, shown
    /// beside the field.
    pub refused: Option<(&'a str, &'a str)>,
    /// What is missing or wrong in the parameters ([`messages`]).
    pub problems: &'a [(String, String)],
    /// The nodes and the elements picked in the viewport, for a label list's Use selection.
    pub picked: Picked,
}

/// The mesh labels picked in the viewport.
#[derive(Clone, Debug, Default)]
pub struct Picked {
    pub nodes: Vec<u64>,
    pub elements: Vec<u64>,
    /// Element faces: (element, local face).
    pub faces: Vec<(u64, u64)>,
}

/// What a frame of the form did.
#[derive(Debug, Default)]
pub struct Outcome {
    pub changed: bool,
    /// Fields whose typed text is not a value yet (OK waits for them).
    pub invalid: usize,
    /// A reference field's "Create …" was clicked: the collection to add to and the field's path.
    pub create: Option<(Collection, String)>,
    /// The paths of the fields edited this frame.
    pub touched: Vec<String>,
    /// The body whose row in a list of bodies is under the pointer (the viewport highlights it).
    pub hover_body: Option<String>,
}

struct Form<'a> {
    cx: &'a Context<'a>,
    style: Style,
    /// Each field's problem, by path, phrased with its label ([`messages`]).
    problems: BTreeMap<String, String>,
    names: BTreeMap<&'static str, Vec<String>>,
    out: Outcome,
}

/// The form's look; `label_width` is the label column ([`label_column`]).
pub fn style(t: &Tokens, label_width: f32) -> Style {
    Style { label_width, dim: t.text_dim, error: t.error, required: t.accent }
}

/// The label column of `fields` in a form `width` wide: its longest label on one line with its
/// marks, within [`LABEL_MIN`] and [`LABEL_SHARE`] of the width. Labels longer than that wrap; the
/// column widens when one would need more than [`fields::LABEL_LINES`] lines, as long as a field
/// keeps [`FIELD_MIN`].
pub fn label_column(ui: &Ui, fields: &[Field], units: &DisplaySystem, width: f32) -> f32 {
    let share = (width * LABEL_SHARE).max(LABEL_MIN);
    let most = (width - FIELD_MIN - ui.spacing().item_spacing.x).max(share);
    let mut column = LABEL_MIN;
    each_label(fields, units, 0, &mut |label, required, help| {
        let one_line = fields::label_width(ui, label, required, help);
        let need = if one_line <= share { one_line } else { wrapped_width(ui, label, required, help, share, one_line.min(most)) };
        column = column.max(need);
    });
    column.min(most)
}

/// The narrowest column in `lo..=hi` that shows `label` whole on [`fields::LABEL_LINES`] lines.
fn wrapped_width(ui: &Ui, label: &str, required: bool, help: bool, mut lo: f32, mut hi: f32) -> f32 {
    let marks = fields::marks_width(required, help);
    let fits = |column: f32| !fields::label_galley(ui, label, egui::Color32::PLACEHOLDER, column - marks).elided;
    if fits(lo) {
        return lo;
    }
    // Bisection to a point.
    for _ in 0..12 {
        if hi - lo < 1.0 {
            break;
        }
        let mid = (lo + hi) / 2.0;
        if fits(mid) { hi = mid } else { lo = mid }
    }
    hi.ceil()
}

/// Narrowest box of a vector component on a shared line: wide enough for `-12000.5`.
const VECTOR_MIN: f32 = 64.0;

/// A definition's name field (`name`, `new_name`), whose text must be a deck-safe name.
fn is_name(f: &Field) -> bool {
    matches!(f.kind, Kind::Text) && matches!(f.key, "name" | "new_name")
}

/// A field's width: the rest of its row.
fn fill(ui: &Ui) -> f32 {
    ui.available_width().max(FIELD_MIN)
}

/// Each row label of `fields` (nested too): its text, whether it is required and has help.
fn each_label(fields: &[Field], units: &DisplaySystem, depth: usize, out: &mut dyn FnMut(&str, bool, bool)) {
    if depth > MAX_DEPTH {
        return;
    }
    for f in fields {
        let help = !f.help_text(units).is_empty();
        match f.kind {
            Kind::Variants { variants, .. } => {
                // A flattened variant's row is the field's own; a keyed one's is "Type".
                if f.key.is_empty() {
                    out(f.label, f.required, help);
                } else {
                    out("Type", true, help);
                }
                for v in variants {
                    each_label(v.fields, units, depth + 1, out);
                }
            }
            Kind::Group(inner) | Kind::Rows { fields: inner, .. } => each_label(inner, units, depth + 1, out),
            _ => out(f.label, f.required, help),
        }
    }
}

/// Draw `fields` editing `params` (an object).
pub fn show(ui: &mut Ui, t: &Tokens, cx: &Context, fields: &'static [Field], params: &mut Map<String, Value>) -> Outcome {
    let problems = cx.problems.iter().cloned().collect();
    let label_width = label_column(ui, fields, cx.units, ui.available_width());
    let mut f = Form { cx, style: style(t, label_width), problems, names: BTreeMap::new(), out: Outcome::default() };
    ui.spacing_mut().item_spacing.y = 3.0;
    f.fields(ui, fields, fields, params, "", 0);
    f.out
}

impl Form<'_> {
    /// The problem with the field at `path`, under it: the command's refusal at once, the
    /// form's own once the field was edited or OK was tried.
    fn problem_note(&self, ui: &mut Ui, path: &str) {
        let message = match self.cx.refused {
            Some((key, message)) if key == path => Some(message.to_string()),
            _ => self.problems.get(path).filter(|_| self.cx.attempted || self.cx.touched.contains(path)).cloned(),
        };
        if let Some(message) = message {
            fields::note(ui, &self.style, &message, self.style.error);
        }
    }

    fn touch(&mut self, path: &str) {
        self.out.changed = true;
        self.out.touched.push(path.to_string());
    }

    fn id(&self, path: &str) -> Id {
        Id::new(("fc_form", self.cx.salt, path))
    }

    /// A value the fields cannot show: it stays as it is, editable in the JSON view.
    fn json_only(&self, ui: &mut Ui) {
        ui.label(RichText::new("set in JSON").small().color(self.style.dim)).on_hover_text("Edit it in the JSON view");
    }

    fn names(&mut self, c: Collection) -> &[String] {
        self.names.entry(c.label()).or_insert_with(|| (self.cx.names)(c))
    }

    /// The fields of one object. `common` are the fields beside a flattened variant (kept when it
    /// switches).
    fn fields(&mut self, ui: &mut Ui, fields: &[Field], common: &[Field], obj: &mut Map<String, Value>, prefix: &str, depth: usize) {
        if depth > MAX_DEPTH {
            return;
        }
        for f in fields {
            if let Kind::Variants { tag, variants } = f.kind
                && f.key.is_empty()
            {
                let current = obj.get(tag).and_then(Value::as_str).map(str::to_string);
                let path = join(prefix, tag);
                if let Some(t) = self.variant_choice(ui, f, variants, current.as_deref(), &path) {
                    *obj = schema::switch_variant(f, common, obj, &t);
                    self.touch(&path);
                }
                self.problem_note(ui, &path);
                if let Some(v) = obj.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                    self.fields(ui, v.fields, common, obj, prefix, depth + 1);
                }
                continue;
            }
            self.field(ui, f, obj, prefix, depth);
        }
    }

    /// A Type row: the variant's drop-down. Returns the tag picked.
    fn variant_choice(&mut self, ui: &mut Ui, f: &Field, variants: &[schema::Variant], current: Option<&str>, path: &str) -> Option<String> {
        let selected = variants.iter().find(|v| Some(v.tag) == current).map_or(current.unwrap_or("—"), |v| v.label);
        let labels: Vec<&str> = variants.iter().map(|v| v.label).collect();
        let id = self.id(path);
        let help = f.help_text(self.cx.units);
        let picked = fields::row(ui, &self.style, f.label, f.required, &help, |ui| fields::combo(ui, id, selected, &labels, fill(ui)));
        picked.and_then(|i| variants.get(i)).map(|v| v.tag.to_string()).filter(|t| Some(t.as_str()) != current)
    }

    fn field(&mut self, ui: &mut Ui, f: &Field, obj: &mut Map<String, Value>, prefix: &str, depth: usize) {
        let path = join(prefix, f.key);
        match f.kind {
            Kind::Group(_) | Kind::Variants { .. } | Kind::Rows { .. } => self.section(ui, f, obj, &path, depth),
            Kind::NameOr { name, object, fields } => self.name_or(ui, f, name, object, fields, obj, &path, depth),
            _ => {
                let value = obj.get(f.key).cloned();
                let style = self.style;
                let help = f.help_text(self.cx.units);
                let new = fields::row(ui, &style, f.label, f.required, &help, |ui| self.scalar(ui, f, value.as_ref(), &path));
                if let Some(new) = new {
                    match new {
                        Value::Null => obj.remove(f.key),
                        v => obj.insert(f.key.into(), v),
                    };
                    self.touch(&path);
                }
                if is_name(f) && ui.data(|d| d.get_temp::<bool>(self.id(&path).with("spaced"))).unwrap_or(false) {
                    fields::note(ui, &self.style, "Spaces become _ (names go into solver decks).", self.style.dim);
                }
                if let Some(note) = f.note(obj.get(f.key).and_then(Value::as_str)) {
                    fields::note(ui, &self.style, note, self.style.dim);
                }
                if let Kind::Either(options) = f.kind {
                    self.either_value(ui, f, options, obj, &path);
                }
            }
        }
        self.problem_note(ui, &path);
    }

    /// A one-row field; returns its new value (`Null` removes it).
    fn scalar(&mut self, ui: &mut Ui, f: &Field, value: Option<&Value>, path: &str) -> Option<Value> {
        let id = self.id(path);
        let s = self.style;
        match f.kind {
            Kind::Number { .. } => {
                let unit = f.unit_label(self.cx.units);
                let width = (fill(ui) - fields::unit_width(ui, unit.as_deref())).max(FIELD_MIN);
                let out = self.number(ui, id, value, f.hint, false, width);
                fields::unit(ui, &s, unit.as_deref());
                out
            }
            Kind::Integer { .. } => self.number(ui, id, value, f.hint, true, fill(ui)),
            Kind::Text | Kind::Path { .. } => {
                let mut text = value.and_then(Value::as_str).unwrap_or_default().to_string();
                let changed = fields::text(ui, &s, id, &mut text, f.hint, fill(ui));
                // Definition names go into solver decks, which take no spaces: a typed space
                // becomes `_`, and the field says so.
                if changed && is_name(f) && text.contains(' ') {
                    text = text.replace(' ', "_");
                    ui.data_mut(|d| d.insert_temp(id.with("spaced"), true));
                }
                changed.then(|| if text.is_empty() && !f.required { Value::Null } else { json!(text) })
            }
            Kind::Bool => {
                // An unset key runs with the field's default, so the box shows that.
                let mut on = value.and_then(Value::as_bool).unwrap_or(f.default == Some("true"));
                let r = ui.checkbox(&mut on, "");
                if !f.hint.is_empty() {
                    r.clone().on_hover_text(f.hint);
                }
                r.changed().then(|| json!(on))
            }
            Kind::Choice(options) => {
                let current = value.and_then(Value::as_str);
                let selected = options
                    .iter()
                    .find(|(v, _)| Some(*v) == current)
                    .map_or(current.unwrap_or(if f.hint.is_empty() { "—" } else { f.hint }), |(_, l)| *l);
                let mut labels: Vec<&str> = options.iter().map(|(_, l)| *l).collect();
                let optional = !f.required;
                if optional {
                    labels.push("(default)");
                }
                fields::combo(ui, id, selected, &labels, fill(ui)).map(|i| options.get(i).map_or(Value::Null, |(v, _)| json!(v)))
            }
            Kind::Vector { labels, nullable, .. } => self.vector(ui, id, value, f, labels, nullable, None),
            Kind::Tuple { labels, whole, .. } => self.tuple(ui, id, value, labels, whole),
            Kind::Ref(c) => self.reference(ui, id, path, value, c, f.hint, f.required),
            Kind::RefList(c) => self.reference_list(ui, path, value, c),
            Kind::RefChoices(c, options) => self.reference_choices(ui, id, value, c, options),
            Kind::RefText(c) => self.reference_text(ui, id, value, c, f.hint),
            Kind::Either(options) => {
                let key = match value {
                    Some(Value::String(s)) => Some(s.as_str()),
                    Some(Value::Object(o)) => o.keys().next().map(String::as_str),
                    _ => None,
                };
                let selected = options.iter().find(|o| Some(o.key) == key).map_or(if f.hint.is_empty() { "—" } else { f.hint }, |o| o.label);
                let mut labels: Vec<&str> = options.iter().map(|o| o.label).collect();
                if !f.required {
                    labels.push("(none)");
                }
                fields::combo(ui, id, selected, &labels, fill(ui)).and_then(|i| match options.get(i) {
                    Some(o) if Some(o.key) == key => None,
                    Some(o) if matches!(o.kind, Kind::Unit) => Some(json!(o.key)),
                    Some(o) => Some(json!({o.key: schema::prefill(&[*o], &json!({})).get(o.key).cloned().unwrap_or(Value::Null)})),
                    None => Some(Value::Null),
                })
            }
            Kind::Json => self.json(ui, id, value, f.hint),
            Kind::Labels(of) => self.labels(ui, id, value, of, f.hint),
            Kind::Sides(side) => self.sides(ui, id, value, side, f.hint),
            Kind::Points { columns, max, .. } => self.points(ui, id, value, f, columns, max),
            Kind::Multi(options) => self.multi(ui, value, options),
            Kind::NumberOrAuto { .. } => self.number_or_auto(ui, id, value, f),
            Kind::PerAxis { labels, .. } => self.per_axis(ui, id, value, f, labels),
            Kind::Unit | Kind::Group(_) | Kind::Variants { .. } | Kind::Rows { .. } | Kind::NameOr { .. } => None,
        }
    }

    /// The chosen option's value of an [`Kind::Either`] field, on its own row.
    fn either_value(&mut self, ui: &mut Ui, f: &Field, options: &[Field], obj: &mut Map<String, Value>, path: &str) {
        let Some(Value::Object(inner)) = obj.get(f.key).cloned() else { return };
        let Some((key, value)) = inner.into_iter().next() else { return };
        let Some(option) = options.iter().find(|o| o.key == key) else { return };
        let sub = join(path, &key);
        let style = self.style;
        let help = option.help_text(self.cx.units);
        let new = fields::row(ui, &style, option.label, false, &help, |ui| self.scalar(ui, option, Some(&value), &sub));
        if let Some(new) = new {
            obj.insert(f.key.into(), json!({key: new}));
            self.touch(path);
        }
    }

    fn number(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, hint: &str, integer: bool, width: f32) -> Option<Value> {
        let current = value.and_then(Value::as_f64);
        if value.is_some_and(|v| !v.is_null() && current.is_none()) {
            self.json_only(ui);
            return None;
        }
        match fields::number(ui, &self.style, id, current, hint, integer, width) {
            NumberEdit::Set(Some(x)) if integer => Some(json!(x as i64)),
            NumberEdit::Set(Some(x)) => Some(json!(x)),
            NumberEdit::Set(None) => Some(Value::Null),
            NumberEdit::Invalid => {
                self.out.invalid += 1;
                None
            }
            NumberEdit::Unchanged => None,
        }
    }

    /// A vector: one number per component, each after its caption (X, Y, Z or the field's own).
    /// The components of a vector; with `cells`, all on one line at that width and without
    /// captions (a table row under its column headers).
    #[allow(clippy::too_many_arguments)]
    fn vector(
        &mut self,
        ui: &mut Ui,
        id: Id,
        value: Option<&Value>,
        f: &Field,
        labels: &[&str],
        nullable: bool,
        cells: Option<f32>,
    ) -> Option<Value> {
        let current: Vec<Value> = match value {
            Some(Value::Array(a)) if a.len() == labels.len() => a.clone(),
            Some(v) if !v.is_null() => {
                self.json_only(ui);
                return None;
            }
            // Empty slots show their placeholder (x, y, z), not a 0 nobody typed.
            _ => vec![Value::Null; labels.len()],
        };
        let mut next = current.clone();
        let mut changed = false;
        let unit = f.unit_label(self.cx.units);
        // Up to three components share a line when each still shows a value such as -12000
        // whole; otherwise each takes a line of its own.
        let gap = 3.0;
        let small = egui::TextStyle::Small.resolve(ui.style());
        let caption_w = labels.iter().map(|l| crate::elide::width(ui, &caption(l), &small)).fold(0.0, f32::max) + 2.0 * gap;
        let room = fill(ui) - fields::unit_width(ui, unit.as_deref());
        let three = labels.len().clamp(1, 3);
        let (per_line, width) = match cells {
            Some(w) => (labels.len().max(1), w),
            None => {
                let per_line = if (room / three as f32 - caption_w) >= VECTOR_MIN { three } else { 1 };
                (per_line, (room / per_line as f32 - caption_w - 1.0).max(34.0))
            }
        };
        ui.vertical(|ui| {
            for chunk in (0..labels.len()).collect::<Vec<_>>().chunks(per_line) {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = gap;
                    for &i in chunk {
                        let hint = labels.get(i).copied().unwrap_or("");
                        if cells.is_none() {
                            fields::caption(ui, &self.style, &caption(hint));
                        }
                        let slot = current.get(i).cloned().unwrap_or(Value::Null);
                        let shown = (!slot.is_null()).then_some(&slot);
                        if let Some(v) = self.number(ui, id.with(i), shown, hint, false, width) {
                            if let Some(n) = next.get_mut(i) {
                                *n = v;
                            }
                            changed = true;
                        }
                    }
                    if chunk.first() == Some(&0) || per_line == 1 {
                        fields::unit(ui, &self.style, unit.as_deref());
                    }
                });
            }
        });
        // A vector left empty is not written until edited (a required one shows as missing).
        // The slots left empty of a vector that takes no blanks are 0.
        changed.then(|| Value::Array(next.into_iter().map(|x| if x.is_null() && !nullable { json!(0.0) } else { x }).collect()))
    }

    /// A number per label on one line, captioned; trailing empty slots are left out (a 2D
    /// block's divisions are the first two of three).
    fn tuple(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, labels: &[&str], whole: bool) -> Option<Value> {
        let mut slots: Vec<Value> = match value {
            None | Some(Value::Null) => Vec::new(),
            Some(Value::Array(a)) if a.len() <= labels.len() => a.clone(),
            Some(_) => {
                self.json_only(ui);
                return None;
            }
        };
        slots.resize(labels.len(), Value::Null);
        let gap = 3.0;
        let small = egui::TextStyle::Small.resolve(ui.style());
        let caption_w = labels.iter().map(|l| crate::elide::width(ui, &caption(l), &small)).fold(0.0, f32::max) + 2.0 * gap;
        // All on one line when each still shows a value whole; otherwise one per line.
        let room = fill(ui);
        let n = labels.len().max(1);
        let per_line = if (room - gap * (n as f32 - 1.0)) / n as f32 - caption_w >= VECTOR_MIN { n } else { 1 };
        let width = ((room - gap * (per_line as f32 - 1.0)) / per_line as f32 - caption_w - 1.0).max(1.0);
        let mut changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.spacing_mut().item_spacing.x = gap;
            for (i, label) in labels.iter().enumerate() {
                if per_line == 1 && i > 0 {
                    ui.end_row();
                }
                fields::caption(ui, &self.style, &caption(label));
                let slot = slots.get(i).cloned().unwrap_or(Value::Null);
                let shown = (!slot.is_null()).then_some(&slot);
                if let Some(v) = self.number(ui, id.with(i), shown, label, whole, width) {
                    if let Some(s) = slots.get_mut(i) {
                        *s = v;
                    }
                    changed = true;
                }
            }
        });
        if !changed {
            return None;
        }
        while slots.last().is_some_and(Value::is_null) {
            slots.pop();
        }
        Some(if slots.is_empty() { Value::Null } else { Value::Array(slots) })
    }

    /// A drop-down of the collection's items; never free text. With none, the drop-down is
    /// greyed and a button opens the command that makes one ([`Outcome::create`]).
    #[allow(clippy::too_many_arguments)]
    fn reference(&mut self, ui: &mut Ui, id: Id, path: &str, value: Option<&Value>, c: Collection, hint: &str, required: bool) -> Option<Value> {
        let current = value.and_then(Value::as_str).filter(|n| !n.is_empty()).map(str::to_string);
        let names = self.names(c).to_vec();
        if names.is_empty() {
            ui.vertical(|ui| {
                let none = format!("No {} yet", plural(c.label()));
                ui.add_enabled_ui(false, |ui| fields::combo(ui, id, current.as_deref().unwrap_or(&none), &[], fill(ui)));
                self.create_button(ui, c, path);
            });
            return None;
        }
        let empty = if hint.is_empty() { "—" } else { hint };
        let selected = current.clone().unwrap_or_else(|| empty.to_string());
        let mut labels: Vec<&str> = names.iter().map(String::as_str).collect();
        if !required {
            labels.push("(default)");
        }
        fields::combo(ui, id, &selected, &labels, fill(ui)).map(|i| names.get(i).map_or(Value::Null, |n| json!(n)))
    }

    fn reference_list(&mut self, ui: &mut Ui, path: &str, value: Option<&Value>, c: Collection) -> Option<Value> {
        let chosen: Vec<String> =
            value.and_then(Value::as_array).map(|a| a.iter().filter_map(Value::as_str).map(str::to_string).collect()).unwrap_or_default();
        let names = self.names(c).to_vec();
        let shown = names.get(..names.len().min(MAX_ROWS)).unwrap_or_default();
        let mut changed = None;
        let list = |picked: Vec<String>| if picked.is_empty() { Value::Null } else { json!(picked) };
        ui.vertical(|ui| {
            if names.is_empty() {
                ui.label(RichText::new(format!("No {} yet", plural(c.label()))).color(self.style.dim));
                self.create_button(ui, c, path);
            }
            // All, none and the others, for long lists.
            if shown.len() > 3 {
                ui.horizontal(|ui| {
                    if ui.small_button("All").clicked() {
                        changed = Some(list(shown.to_vec()));
                    }
                    if ui.small_button("None").clicked() {
                        changed = Some(Value::Null);
                    }
                    if ui.small_button("Invert").on_hover_text("Tick the unticked ones and untick the others").clicked() {
                        changed = Some(list(shown.iter().filter(|n| !chosen.contains(n)).cloned().collect()));
                    }
                });
            }
            for n in shown {
                let mut on = chosen.contains(n);
                let r = ui.checkbox(&mut on, n);
                if r.hovered() && c == Collection::Body {
                    self.out.hover_body = Some(n.clone());
                }
                if r.changed() {
                    let mut next: Vec<String> = names.iter().filter(|m| if *m == n { on } else { chosen.contains(m) }).cloned().collect();
                    next.dedup();
                    changed = Some(list(next));
                }
            }
        });
        changed
    }

    /// A drop-down per item of `c` (its row hovered highlights a body), "(default)" leaving it
    /// out of the object.
    fn reference_choices(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, c: Collection, options: &[(&str, &str)]) -> Option<Value> {
        let current = value.and_then(Value::as_object).cloned().unwrap_or_default();
        let names = self.names(c).to_vec();
        let mut labels: Vec<&str> = options.iter().map(|(_, l)| *l).collect();
        labels.push("(default)");
        let mut changed = None;
        ui.vertical(|ui| {
            if names.is_empty() {
                ui.label(RichText::new(format!("No {} yet", plural(c.label()))).color(self.style.dim));
            }
            for (i, n) in names.iter().take(MAX_ROWS).enumerate() {
                let chosen = current.get(n).and_then(Value::as_str);
                let selected = options.iter().find(|(v, _)| Some(*v) == chosen).map_or("(default)", |(_, l)| *l);
                let row = ui.horizontal(|ui| {
                    let width = (ui.available_width() * 0.5).max(FIELD_MIN);
                    ui.add_sized([width, ui.spacing().interact_size.y], egui::Label::new(n.as_str()).truncate());
                    fields::combo(ui, id.with(i), selected, &labels, fill(ui))
                });
                if row.response.contains_pointer() && c == Collection::Body {
                    self.out.hover_body = Some(n.clone());
                }
                if let Some(k) = row.inner {
                    let mut next = current.clone();
                    match options.get(k) {
                        Some((v, _)) => next.insert(n.clone(), json!(v)),
                        None => next.remove(n),
                    };
                    changed = Some(if next.is_empty() { Value::Null } else { Value::Object(next) });
                }
            }
        });
        changed
    }

    /// A text per item of `c`, each after the item's name; an empty text leaves the item out.
    fn reference_text(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, c: Collection, hint: &str) -> Option<Value> {
        let current = value.and_then(Value::as_object).cloned().unwrap_or_default();
        let names = self.names(c).to_vec();
        let mut changed = None;
        ui.vertical(|ui| {
            if names.is_empty() {
                ui.label(RichText::new(format!("No {} yet", plural(c.label()))).color(self.style.dim));
            }
            for (i, n) in names.iter().take(MAX_ROWS).enumerate() {
                let mut text = current.get(n).and_then(Value::as_str).unwrap_or_default().to_string();
                let row = ui.horizontal(|ui| {
                    let width = (ui.available_width() * 0.5).max(FIELD_MIN);
                    ui.add_sized([width, ui.spacing().interact_size.y], egui::Label::new(n.as_str()).truncate());
                    fields::text(ui, &self.style, id.with(i), &mut text, hint, fill(ui))
                });
                if row.response.contains_pointer() && c == Collection::Body {
                    self.out.hover_body = Some(n.clone());
                }
                if row.inner {
                    let mut next = current.clone();
                    if text.is_empty() {
                        next.remove(n);
                    } else {
                        next.insert(n.clone(), json!(text));
                    }
                    changed = Some(if next.is_empty() { Value::Null } else { Value::Object(next) });
                }
            }
        });
        changed
    }

    /// A [`Kind::NameOr`] field: a choice between the name and the object, then the chosen one's
    /// field or fields. Switching starts the other form empty.
    #[allow(clippy::too_many_arguments)]
    fn name_or(
        &mut self,
        ui: &mut Ui,
        f: &Field,
        name: &'static Field,
        object: &str,
        fields: &'static [Field],
        obj: &mut Map<String, Value>,
        path: &str,
        depth: usize,
    ) {
        let is_object = matches!(obj.get(f.key), Some(Value::Object(_)));
        let id = self.id(path);
        let help = f.help_text(self.cx.units);
        let style = self.style;
        let labels = [name.label, object];
        let selected = if is_object { object } else { name.label };
        if let Some(i) = fields::row(ui, &style, f.label, f.required, &help, |ui| fields::combo(ui, id.with("form"), selected, &labels, fill(ui)))
            && (i == 1) != is_object
        {
            if i == 1 {
                obj.insert(f.key.into(), schema::prefill(fields, &json!({})));
            } else {
                obj.remove(f.key);
            }
            self.touch(path);
        }
        if is_object {
            ui.horizontal(|ui| {
                ui.add_space(INDENT);
                ui.vertical(|ui| {
                    if let Some(Value::Object(child)) = obj.get_mut(f.key) {
                        self.fields(ui, fields, fields, child, path, depth + 1);
                    }
                });
            });
        } else {
            let value = obj.get(f.key).cloned();
            let help = name.help_text(self.cx.units);
            let new = fields::row(ui, &style, name.label, false, &help, |ui| self.scalar(ui, name, value.as_ref(), path));
            if let Some(new) = new {
                match new {
                    Value::Null => obj.remove(f.key),
                    v => obj.insert(f.key.into(), v),
                };
                self.touch(path);
            }
        }
    }

    /// "Create material…": opens the command that makes an item of `c`, then returns here.
    fn create_button(&mut self, ui: &mut Ui, c: Collection, path: &str) {
        let r = c.requirement();
        let text = match r.provider().ends_with(".create") {
            true => format!("Create {}…", c.label()),
            false => format!("{}…", fieldcraft_engine::catalogue::entry_of(r.provider()).map_or(r.provider(), |e| e.label)),
        };
        if ui.button(text).on_hover_text(format!("{} ({})", r.place(), r.provider())).clicked() {
            self.out.create = Some((c, path.to_string()));
        }
    }

    /// A list of node or element labels typed as numbers and ranges (`1-5, 8`), kept while it
    /// does not parse; Use selection takes the labels picked in the viewport.
    fn labels(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, of: schema::LabelOf, hint: &str) -> Option<Value> {
        let current: Option<Vec<u64>> = match value {
            None | Some(Value::Null) => Some(Vec::new()),
            Some(v) => v.as_array().and_then(|a| a.iter().map(Value::as_u64).collect()),
        };
        let Some(current) = current else {
            self.json_only(ui);
            return None;
        };
        let shown = fields::format_labels(&current);
        let mut text = ui.data(|d| d.get_temp::<(String, String)>(id)).filter(|(_, over)| *over == shown).map_or_else(|| shown.clone(), |(t, _)| t);
        let picked = match of {
            schema::LabelOf::Node => self.cx.picked.nodes.clone(),
            schema::LabelOf::Element => self.cx.picked.elements.clone(),
        };
        let mut out = None;
        ui.vertical(|ui| {
            let parsed = fields::parse_labels(&text, schema::MAX_LABELS);
            let mut edit =
                egui::TextEdit::singleline(&mut text).id(id).desired_width(fill(ui)).hint_text(RichText::new(hint).italics().color(self.style.dim));
            if parsed.is_none() {
                edit = edit.text_color(self.style.error);
                self.out.invalid += 1;
            }
            if ui.add(edit).changed()
                && let Some(l) = fields::parse_labels(&text, schema::MAX_LABELS)
            {
                out = Some(if l.is_empty() { Value::Null } else { json!(l) });
            }
            let label = format!("Use selection ({} {}{})", picked.len(), of.label(), if picked.len() == 1 { "" } else { "s" });
            let tip = format!("The {}s picked in the viewport (Select: {}s)", of.label(), of.label());
            if ui
                .add_enabled(!picked.is_empty(), egui::Button::new(RichText::new(label).small()))
                .on_hover_text(&tip)
                .on_disabled_hover_text(&tip)
                .clicked()
            {
                text = fields::format_labels(&picked);
                out = Some(json!(picked));
            }
        });
        let over = out.as_ref().map_or(shown, |v| {
            v.as_array().map_or_else(String::new, |a| fields::format_labels(&a.iter().filter_map(Value::as_u64).collect::<Vec<_>>()))
        });
        ui.data_mut(|d| d.insert_temp(id, (text, over)));
        out
    }

    /// Element faces or edges typed as `element:side` ranges ([`fields::parse_sides`]); for faces,
    /// Use selection takes the faces picked in the viewport.
    fn sides(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, side: schema::Side, hint: &str) -> Option<Value> {
        let key = side.key();
        let current: Option<Vec<(u64, u64)>> = match value {
            None | Some(Value::Null) => Some(Vec::new()),
            Some(v) => v.as_array().and_then(|a| a.iter().map(|x| Some((x.get("element")?.as_u64()?, x.get(key)?.as_u64()?))).collect()),
        };
        let Some(current) = current else {
            self.json_only(ui);
            return None;
        };
        let to_value =
            |l: &[(u64, u64)]| if l.is_empty() { Value::Null } else { Value::Array(l.iter().map(|(e, s)| json!({"element": e, key: s})).collect()) };
        let shown = fields::format_sides(&current);
        let mut text = ui.data(|d| d.get_temp::<(String, String)>(id)).filter(|(_, over)| *over == shown).map_or_else(|| shown.clone(), |(t, _)| t);
        let mut out = None;
        let mut over = shown.clone();
        ui.vertical(|ui| {
            let parsed = fields::parse_sides(&text, schema::MAX_LABELS);
            let mut edit =
                egui::TextEdit::singleline(&mut text).id(id).desired_width(fill(ui)).hint_text(RichText::new(hint).italics().color(self.style.dim));
            if parsed.is_none() {
                edit = edit.text_color(self.style.error);
                self.out.invalid += 1;
            }
            if ui.add(edit).changed()
                && let Some(l) = fields::parse_sides(&text, schema::MAX_LABELS)
            {
                over = fields::format_sides(&l);
                out = Some(to_value(&l));
            }
            if side == schema::Side::Face {
                let picked = self.cx.picked.faces.clone();
                let label = format!("Use selection ({} face{})", picked.len(), if picked.len() == 1 { "" } else { "s" });
                let tip = "The element faces picked in the viewport (Select: Element faces)";
                if ui
                    .add_enabled(!picked.is_empty(), egui::Button::new(RichText::new(label).small()))
                    .on_hover_text(tip)
                    .on_disabled_hover_text(tip)
                    .clicked()
                {
                    text = fields::format_sides(&picked);
                    over = text.clone();
                    out = Some(to_value(&picked));
                }
            }
        });
        ui.data_mut(|d| d.insert_temp(id, (text, over)));
        out
    }

    /// A list of points: a row of x, y, z each, removable, and Add point.
    fn points(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, f: &Field, columns: &[&str], max: usize) -> Option<Value> {
        let current: Vec<Value> = match value {
            None | Some(Value::Null) => Vec::new(),
            Some(Value::Array(a)) if a.iter().all(Value::is_array) => a.clone(),
            Some(_) => {
                self.json_only(ui);
                return None;
            }
        };
        if current.len() > MAX_ROWS {
            self.json_only(ui);
            return None;
        }
        let mut next = current.clone();
        let (mut changed, mut remove) = (false, None);
        // One line per row: the columns' names head the table once, the cells share the width
        // left by the remove button and the unit.
        let gap = 3.0;
        let remove_w = ui.spacing().interact_size.y + gap;
        let unit = f.unit_label(self.cx.units);
        let room = fill(ui) - remove_w - fields::unit_width(ui, unit.as_deref());
        let n = columns.len().max(1) as f32;
        let cell = ((room - gap * (n - 1.0)) / n - 1.0).max(34.0);
        ui.vertical(|ui| {
            if !current.is_empty() {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = gap;
                    ui.add_space(remove_w - gap);
                    for c in columns {
                        ui.allocate_ui_with_layout(egui::vec2(cell, 0.0), egui::Layout::left_to_right(egui::Align::Center), |ui| {
                            ui.set_width(cell);
                            ui.add(egui::Label::new(RichText::new(caption(c)).small().color(self.style.dim)).truncate());
                        });
                    }
                });
            }
            for (i, p) in current.iter().enumerate() {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = gap;
                    if ui.small_button("×").on_hover_text("Remove this row").clicked() {
                        remove = Some(i);
                    }
                    if let Some(v) = self.vector(ui, id.with(i), Some(p), f, columns, false, Some(cell))
                        && let Some(slot) = next.get_mut(i)
                    {
                        *slot = v;
                        changed = true;
                    }
                });
            }
            if ui
                .add_enabled(
                    current.len() < max.min(MAX_ROWS),
                    egui::Button::new(RichText::new(if columns == schema::XYZ { "Add point" } else { "Add row" }).small()),
                )
                .clicked()
            {
                next.push(Value::Array(vec![json!(0.0); columns.len()]));
                changed = true;
            }
        });
        if let Some(i) = remove {
            next.remove(i);
            changed = true;
        }
        changed.then_some(if next.is_empty() { Value::Null } else { Value::Array(next) })
    }

    /// Any of `options`: a checkbox each; none ticked leaves the key out.
    fn multi(&mut self, ui: &mut Ui, value: Option<&Value>, options: &[(&str, &str)]) -> Option<Value> {
        let chosen: Vec<Value> = value.and_then(Value::as_array).cloned().unwrap_or_default();
        let mut out = None;
        ui.horizontal_wrapped(|ui| {
            for (o, label) in options {
                let v = schema::multi_value(o);
                let mut on = chosen.contains(&v);
                if ui.checkbox(&mut on, *label).changed() {
                    // In the options' order.
                    let next: Vec<Value> =
                        options.iter().map(|(x, _)| schema::multi_value(x)).filter(|x| if *x == v { on } else { chosen.contains(x) }).collect();
                    out = Some(if next.is_empty() { Value::Null } else { Value::Array(next) });
                }
            }
        });
        out
    }

    /// A number, or Auto: the checkbox, then the number when it is off.
    fn number_or_auto(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, f: &Field) -> Option<Value> {
        let auto = match value {
            None | Some(Value::Null) => f.default.is_some_and(|d| d.contains("auto")),
            Some(Value::String(s)) => s == "auto",
            Some(_) => false,
        };
        let mut on = auto;
        let mut out = None;
        ui.horizontal(|ui| {
            if ui.checkbox(&mut on, "Auto").changed() {
                out = Some(if on { json!("auto") } else { json!(1.0) });
            }
            if !auto {
                let unit = f.unit_label(self.cx.units);
                let width = (fill(ui) - fields::unit_width(ui, unit.as_deref())).max(FIELD_MIN);
                if let Some(v) = self.number(ui, id, value, f.hint, false, width) {
                    out = Some(v);
                }
                fields::unit(ui, &self.style, unit.as_deref());
            }
        });
        out
    }

    /// One number, or with Per axis ticked one per axis: ticking it copies the number to every
    /// axis, unticking keeps the first axis's value.
    fn per_axis(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, f: &Field, labels: &[&str]) -> Option<Value> {
        let split = matches!(value, Some(Value::Array(_)));
        let mut on = split;
        let mut out = None;
        ui.vertical(|ui| {
            if ui.checkbox(&mut on, "Per axis").on_hover_text("A different value along each axis").changed() {
                out = Some(match (on, value) {
                    (true, Some(v)) if v.is_number() => Value::Array(vec![v.clone(); labels.len()]),
                    (true, _) => Value::Null,
                    (false, Some(Value::Array(a))) => a.first().cloned().unwrap_or(Value::Null),
                    (false, _) => Value::Null,
                });
            }
            if split {
                if let Some(v) = self.vector(ui, id.with("axes"), value, f, labels, false, None) {
                    out = Some(v);
                }
            } else {
                ui.horizontal(|ui| {
                    let unit = f.unit_label(self.cx.units);
                    let width = (fill(ui) - fields::unit_width(ui, unit.as_deref())).max(FIELD_MIN);
                    if let Some(v) = self.number(ui, id, value, f.hint, false, width) {
                        out = Some(v);
                    }
                    fields::unit(ui, &self.style, unit.as_deref());
                });
            }
        });
        out
    }

    /// A JSON value typed as text: kept while it does not parse.
    /// Empty, it shows `hint` (the expected format) greyed.
    fn json(&mut self, ui: &mut Ui, id: Id, value: Option<&Value>, hint: &str) -> Option<Value> {
        let shown = value.map(Value::to_string).unwrap_or_default();
        let mut text = ui.data(|d| d.get_temp::<(String, String)>(id)).filter(|(_, over)| *over == shown).map_or_else(|| shown.clone(), |(t, _)| t);
        let bad = !text.trim().is_empty() && serde_json::from_str::<Value>(&text).is_err();
        let mut edit = egui::TextEdit::singleline(&mut text)
            .id(id)
            .code_editor()
            .desired_width(fill(ui))
            .hint_text(RichText::new(hint).italics().color(self.style.dim));
        if bad {
            edit = edit.text_color(self.style.error);
            self.out.invalid += 1;
        }
        let r = ui.add(edit);
        let mut out = None;
        if r.changed() {
            if text.trim().is_empty() {
                out = Some(Value::Null);
            } else if let Ok(v) = serde_json::from_str::<Value>(&text) {
                out = Some(v);
            }
        }
        let over = out.as_ref().map_or(shown, |v| if v.is_null() { String::new() } else { v.to_string() });
        ui.data_mut(|d| d.insert_temp(id, (text, over)));
        out
    }

    /// A group, a keyed variant or a list of rows: a sub-heading (with a checkbox when optional),
    /// then its fields, indented.
    fn section(&mut self, ui: &mut Ui, f: &Field, obj: &mut Map<String, Value>, path: &str, depth: usize) {
        let present = obj.get(f.key).is_some_and(|v| !v.is_null());
        let mut on = present;
        ui.add_space(2.0);
        let help = f.help_text(self.cx.units);
        let heading = if f.required { ui.label(RichText::new(f.label).strong()) } else { ui.checkbox(&mut on, RichText::new(f.label).strong()) };
        if !help.is_empty() {
            heading.on_hover_text(help);
        }
        if !f.hint.is_empty() && !present {
            fields::note(ui, &Style { label_width: INDENT, ..self.style }, f.hint, self.style.dim);
        }
        if on != present {
            if on {
                let filled = schema::prefill(&[Field { required: true, ..*f }], &json!({}));
                let start = filled.get(f.key).cloned().unwrap_or_else(|| if matches!(f.kind, Kind::Rows { .. }) { json!([]) } else { json!({}) });
                obj.insert(f.key.into(), start);
            } else {
                obj.remove(f.key);
            }
            self.touch(path);
        }
        if !on {
            return;
        }
        ui.horizontal(|ui| {
            ui.add_space(INDENT);
            ui.vertical(|ui| match f.kind {
                Kind::Group(fields) => {
                    if let Some(Value::Object(child)) = obj.get_mut(f.key) {
                        self.fields(ui, fields, fields, child, path, depth + 1);
                    }
                }
                Kind::Variants { tag, variants } => {
                    let Some(Value::Object(child)) = obj.get_mut(f.key) else {
                        self.json_only(ui);
                        return;
                    };
                    let current = child.get(tag).and_then(Value::as_str).map(str::to_string);
                    let row = Field { label: "Type", required: true, ..*f };
                    let tag_path = join(path, tag);
                    if let Some(t) = self.variant_choice(ui, &row, variants, current.as_deref(), &tag_path) {
                        *child = schema::switch_variant(f, &[], child, &t);
                        self.touch(&tag_path);
                    }
                    self.problem_note(ui, &tag_path);
                    if let Some(v) = child.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                        self.fields(ui, v.fields, v.fields, child, path, depth + 1);
                    }
                }
                Kind::Rows { fields, max } => self.rows(ui, fields, max, obj, f.key, path, depth),
                _ => {}
            });
        });
    }

    #[allow(clippy::too_many_arguments)]
    fn rows(&mut self, ui: &mut Ui, fields: &[Field], max: usize, obj: &mut Map<String, Value>, key: &str, path: &str, depth: usize) {
        let Some(Value::Array(list)) = obj.get_mut(key) else {
            self.json_only(ui);
            return;
        };
        let mut remove = None;
        for (i, row) in list.iter_mut().enumerate().take(max) {
            ui.horizontal(|ui| {
                ui.label(RichText::new(format!("{}", i + 1)).small().color(self.style.dim));
                if ui.small_button("×").on_hover_text("Remove").clicked() {
                    remove = Some(i);
                }
            });
            if let Value::Object(child) = row {
                self.fields(ui, fields, fields, child, &format!("{path}/{i}"), depth + 1);
            }
        }
        if let Some(i) = remove {
            list.remove(i);
            self.touch(path);
        }
        if list.len() < max && ui.small_button("+ Add").clicked() {
            list.push(schema::prefill(fields, &json!({})));
            self.touch(path);
        }
    }
}

/// A vector component's caption: axis names upper-cased ("x" → "X"), others as they are ("E1").
fn caption(component: &str) -> String {
    if component.len() <= 3 && component.chars().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit()) {
        component.to_uppercase()
    } else {
        component.to_string()
    }
}

/// Whether any of `fields` (nested too) is a number whose unit follows the unit system.
/// What is missing or wrong in `params`, in the form's order: (field path, message), each
/// phrased with the field's label as the form shows it ([`schema::Problem::for_form`]).
pub fn messages(fields: &'static [Field], params: &Map<String, Value>) -> Vec<(String, String)> {
    schema::object_problems(fields, params)
        .into_iter()
        .map(|p| {
            let m = p.for_form();
            (p.path, m)
        })
        .collect()
}

/// A collection's label in the plural: "materials", "point masses", "bodies".
fn plural(label: &str) -> String {
    if let Some(stem) = label.strip_suffix('y') {
        format!("{stem}ies")
    } else if label.ends_with('s') {
        format!("{label}es")
    } else {
        format!("{label}s")
    }
}

/// Fill each empty required reference that has one obvious candidate, `pick(collection)` (the
/// browser's selection, or a collection's only item). The form starts with it, so the fields and
/// the JSON view agree. Returns whether anything was filled.
pub fn preselect(fields: &[Field], obj: &mut Map<String, Value>, pick: &dyn Fn(Collection) -> Option<String>) -> bool {
    preselect_in(fields, obj, pick, 0)
}

fn preselect_in(fields: &[Field], obj: &mut Map<String, Value>, pick: &dyn Fn(Collection) -> Option<String>, depth: usize) -> bool {
    if depth > MAX_DEPTH {
        return false;
    }
    let mut changed = false;
    for f in fields {
        match f.kind {
            Kind::Variants { tag, .. } if f.key.is_empty() => {
                if let Some(v) = obj.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t)) {
                    changed |= preselect_in(v.fields, obj, pick, depth + 1);
                }
            }
            Kind::Group(inner) => {
                if let Some(Value::Object(child)) = obj.get_mut(f.key) {
                    changed |= preselect_in(inner, child, pick, depth + 1);
                }
            }
            Kind::Variants { tag, .. } => {
                if let Some(Value::Object(child)) = obj.get_mut(f.key)
                    && let Some(v) = child.get(tag).and_then(Value::as_str).and_then(|t| f.kind.variant(t))
                {
                    changed |= preselect_in(v.fields, child, pick, depth + 1);
                }
            }
            Kind::Ref(c) if f.required && obj.get(f.key).is_none_or(Value::is_null) => {
                if let Some(name) = pick(c) {
                    obj.insert(f.key.into(), json!(name));
                    changed = true;
                }
            }
            _ => {}
        }
    }
    changed
}

#[cfg(test)]
mod tests {
    use egui::Color32;

    use super::*;

    /// Every command's form draws, with its defaults and with hostile parameters, with and without
    /// items to reference, without panicking, and drawing alone changes nothing.
    #[test]
    fn every_form_draws() {
        let ctx = egui::Context::default();
        let one = |c: Collection| vec![format!("{}1", c.label())];
        let none = |_: Collection| Vec::new();
        let hostile =
            [json!({}), json!({"type": 7, "name": [], "elastic": "x", "target": {"set": 1, "reference_point": 2}, "values": [1], "geometry": {}})];
        for (c, names) in fieldcraft_engine::commands().iter().flat_map(|c| [(c, &one as &dyn Fn(Collection) -> Vec<String>), (c, &none)]) {
            let Some(form) = c.form else { continue };
            for start in hostile.iter().chain([&schema::prefill(form, &json!({}))]) {
                let mut params = start.as_object().cloned().unwrap_or_default();
                let before = params.clone();
                let mut outcome = Outcome::default();
                let width = crate::context::WIDTH - 16.0;
                let mut used = 0.0;
                let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                    let touched = BTreeSet::new();
                    let cx = Context {
                        units: &mm(),
                        names,
                        salt: Id::new(c.id),
                        touched: &touched,
                        attempted: true,
                        refused: None,
                        problems: &messages(form, &params),
                        picked: Picked::default(),
                    };
                    ui.vertical(|ui| {
                        ui.set_max_width(width);
                        outcome = show(ui, &Tokens::get(false), &cx, form, &mut params);
                        used = ui.min_rect().width();
                    });
                });
                // Wider would push the context panel over the viewport.
                assert!(used <= width + 0.5, "{}: {used}", c.id);
                out.textures_delta.clear();
                assert!(!outcome.changed, "{}", c.id);
                assert_eq!(params, before, "{}", c.id);
            }
        }
    }

    /// Every form's labels show whole, wrapped at most onto two lines, with their marks inside
    /// the label column, at the context panel's width and at a narrow one, in a 1280×800 window.
    #[test]
    fn no_label_or_mark_is_clipped() {
        let ctx = egui::Context::default();
        let input =
            || egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0))), ..Default::default() };
        let one = |c: Collection| vec![format!("{} with a long name", c.label())];
        let mut checked = 0;
        for width in [crate::context::WIDTH - 16.0, 200.0] {
            for c in fieldcraft_engine::commands() {
                let Some(form) = c.form else { continue };
                let mut params = schema::prefill(form, &json!({})).as_object().cloned().unwrap_or_default();
                let mut labels = BTreeSet::new();
                each_label(form, &mm(), 0, &mut |l, _, _| {
                    labels.insert(l.to_string());
                });
                let (mut left, mut column, mut cut) = (0.0, 0.0, Vec::new());
                let mut out = ctx.run_ui(input(), |ui| {
                    let touched = BTreeSet::new();
                    let cx = Context {
                        units: &mm(),
                        names: &one,
                        salt: Id::new(c.id),
                        touched: &touched,
                        attempted: false,
                        refused: None,
                        problems: &[],
                        picked: Picked::default(),
                    };
                    ui.vertical(|ui| {
                        ui.set_max_width(width);
                        left = ui.min_rect().left();
                        column = label_column(ui, form, cx.units, width);
                        // Every label (also of sections not open), laid out as its row lays it out.
                        each_label(form, cx.units, 0, &mut |label, required, help| {
                            let wrap = column - fields::marks_width(required, help);
                            if fields::label_galley(ui, label, Color32::WHITE, wrap).elided {
                                cut.push(label.to_string());
                            }
                        });
                        show(ui, &Tokens::get(false), &cx, form, &mut params);
                    });
                });
                out.textures_delta.clear();
                assert!(cut.is_empty(), "{} at {width}: {cut:?} do not fit {column}", c.id);
                assert!(column <= width - FIELD_MIN, "{}: {column}", c.id);
                // A row's column starts at its label (nested rows are indented); its marks follow.
                let mut row_left = left;
                for shape in &out.shapes {
                    let egui::Shape::Text(t) = &shape.shape else { continue };
                    let text = t.galley.text();
                    let rect = t.galley.rect.translate(t.pos.to_vec2());
                    let end = row_left + column + 0.5;
                    if labels.contains(text) && rect.left() < left + width / 2.0 {
                        row_left = rect.left();
                        assert!(!t.galley.elided, "{} at {width}: `{text}` elided", c.id);
                        assert!(rect.right() <= rect.left() + column + 0.5, "{} at {width}: `{text}` past the column", c.id);
                        checked += 1;
                    } else if (text == "*" || text == "i") && rect.left() < end {
                        assert!(rect.right() <= end, "{} at {width}: a mark past the column: {rect:?}", c.id);
                    }
                }
            }
        }
        assert!(checked > 100, "{checked}");
    }

    /// Long labels share the column: the Plasticity form's "Hardening modulus" shows whole.
    #[test]
    fn hardening_modulus_shows_whole() {
        let text = drawn("model.material.plasticity", &mm(), json!({"plastic": {"type": "bilinear"}}));
        assert!(text.lines().any(|l| l == "Hardening modulus"), "{text}");
    }

    /// In every form, each required field left empty is a problem, so OK waits, and the message
    /// names the field by its label, never by a JSON key.
    #[test]
    fn every_required_field_left_empty_is_named_by_its_label() {
        let mut checked = 0;
        for c in fieldcraft_engine::commands() {
            let Some(form) = c.form else { continue };
            let full = schema::prefill(form, &json!({}));
            let Some(full) = full.as_object() else { continue };
            for (path, message) in messages(form, full) {
                assert!(!message.contains('`'), "{}: {path}: {message}", c.id);
            }
            for f in form.iter().filter(|f| f.required) {
                let key = match f.kind {
                    Kind::Variants { tag, .. } if f.key.is_empty() => tag,
                    _ => f.key,
                };
                let mut params = full.clone();
                params.remove(key);
                let found: BTreeMap<String, String> = messages(form, &params).into_iter().collect();
                let message = found.get(key).unwrap_or_else(|| panic!("{}: `{key}` empty is no problem: {found:?}", c.id));
                assert!(message.starts_with(&format!("{}: ", f.label)), "{}: {message}", c.id);
                for m in found.values() {
                    assert!(!m.contains('`'), "{}: {m}", c.id);
                }
                checked += 1;
            }
        }
        assert!(checked > 50, "{checked}");
    }

    #[test]
    fn messages_read_as_the_form_shows_them() {
        let plasticity = fieldcraft_engine::find_command("model.material.plasticity").and_then(|c| c.form).unwrap();
        let m = messages(plasticity, &Map::new());
        assert_eq!(m.first(), Some(&("name".to_string(), "Material: choose a material".to_string())), "{m:?}");
        let load = fieldcraft_engine::find_command("model.load.create").and_then(|c| c.form).unwrap();
        let m: Vec<String> = messages(load, &json!({"type": "force"}).as_object().cloned().unwrap_or_default()).into_iter().map(|(_, m)| m).collect();
        assert_eq!(m, ["Name: enter a name", "Target: choose a set or a reference point", "Force: required"]);
    }

    /// A problem waits until its field is edited or OK is tried; a refused parameter shows at once.
    #[test]
    fn problems_show_once_touched_or_attempted() {
        let form = fieldcraft_engine::find_command("model.material.plasticity").and_then(|c| c.form).unwrap();
        let ctx = egui::Context::default();
        let draw = |touched: &BTreeSet<String>, attempted: bool, refused: Option<(&str, &str)>| {
            let mut shown = String::new();
            let mut params = Map::new();
            let names = |_: Collection| Vec::new();
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                let cx = Context {
                    units: &mm(),
                    names: &names,
                    salt: Id::new("t"),
                    touched,
                    attempted,
                    refused,
                    problems: &messages(form, &params),
                    picked: Picked::default(),
                };
                show(ui, &Tokens::get(false), &cx, form, &mut params);
            });
            out.textures_delta.clear();
            // The text drawn this frame.
            for clipped in &out.shapes {
                if let egui::Shape::Text(t) = &clipped.shape {
                    shown.push_str(t.galley.text());
                    shown.push('\n');
                }
            }
            shown
        };
        let none = BTreeSet::new();
        assert!(!draw(&none, false, None).contains("choose a material"));
        assert!(draw(&none, true, None).contains("Material: choose a material"));
        let touched: BTreeSet<String> = ["name".to_string()].into();
        assert!(draw(&touched, false, None).contains("Material: choose a material"));
        let refused = draw(&none, false, Some(("name", "Material: names no material")));
        assert!(refused.contains("Material: names no material") && !refused.contains("choose a material"), "{refused}");
    }

    /// The text a form draws, with `units`.
    fn drawn(id: &str, units: &DisplaySystem, params: Value) -> String {
        let form = fieldcraft_engine::find_command(id).and_then(|c| c.form).unwrap();
        let ctx = egui::Context::default();
        let mut params = params.as_object().cloned().unwrap_or_default();
        let names = |_: Collection| Vec::new();
        let touched = BTreeSet::new();
        let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
            let cx = Context {
                units,
                names: &names,
                salt: Id::new("t"),
                touched: &touched,
                attempted: false,
                refused: None,
                problems: &[],
                picked: Picked { nodes: vec![4, 5], ..Picked::default() },
            };
            show(ui, &Tokens::get(false), &cx, form, &mut params);
        });
        out.textures_delta.clear();
        out.shapes
            .iter()
            .filter_map(|c| if let egui::Shape::Text(t) = &c.shape { Some(t.galley.text().to_string()) } else { None })
            .collect::<Vec<_>>()
            .join("\n")
    }

    fn mm() -> DisplaySystem {
        DisplaySystem::new(fieldcraft_model::units::Preset::MmTonneS)
    }

    /// A number shows its unit in the display system; a vector's components are captioned.
    #[test]
    fn units_follow_the_display_system_and_vector_captions_show() {
        let mpa = drawn("model.material.plasticity", &mm(), json!({"plastic": {"type": "bilinear"}}));
        assert!(mpa.lines().any(|l| l == "MPa"), "{mpa}");
        let si = DisplaySystem::new(fieldcraft_model::units::Preset::MKgS);
        let pa = drawn("model.material.plasticity", &si, json!({"plastic": {"type": "bilinear"}}));
        assert!(pa.lines().any(|l| l == "Pa") && !pa.contains("MPa"), "{pa}");
        let load = drawn("model.load.create", &si, json!({"type": "force"}));
        let lines: Vec<&str> = load.lines().collect();
        assert!(["X", "Y", "Z", "N"].iter().all(|w| lines.contains(w)), "{load}");
    }

    /// Label lists, point lists, multi-choices and number-or-auto are typed fields (#169): a
    /// label list shows its ranges and offers the picked nodes, points show their coordinates,
    /// output choices are checkboxes, the scale an Auto box.
    #[test]
    fn typed_lists_choices_and_auto_show() {
        let probe = drawn("solve.probe", &mm(), json!({"nodes": [1, 2, 3, 7], "at": [[1.5, 0, 0]]}));
        assert!(
            probe.contains("1-3, 7") && probe.contains("Use selection (2 nodes)") && probe.contains("1.5") && probe.contains("Add point"),
            "{probe}"
        );
        assert!(!probe.contains("set in JSON"), "{probe}");
        let step = drawn("model.step.create", &mm(), json!({"name": "S", "procedure": {"type": "static"}, "output": ["stress"]}));
        assert!(step.contains("Reaction force") && step.contains("Heat flux"), "{step}");
        let deformed = drawn("results.deformed", &mm(), json!({"scale": "auto"}));
        assert!(deformed.contains("Auto"), "{deformed}");
    }
}
