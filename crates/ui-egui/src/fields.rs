//! Form field widgets: a labelled row, a number field that keeps what is being typed until it
//! is a number, a drop-down, a text field. The controls are the UI kit's (its fields, drop-down
//! and help icon, at its row height); what is FieldCraft's own is the row's label column: a long
//! label wraps inside it onto two lines, after it the required mark, so a form shows every label
//! whole and says which fields it needs.
//!
//! [`row`] is a deliberate fork of `kit::row`, which truncates its label to one line and has no
//! required mark. Moving both into the kit is requested upstream (SolidCraft); once `kit::row`
//! wraps and marks required fields, this row goes and the forms use the kit's.

use egui::{Color32, Id, RichText, Ui, vec2};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;

/// Widths and colours of a form.
#[derive(Clone, Copy, Debug)]
pub struct Style {
    /// The label column.
    pub label_width: f32,
    /// The room every row keeps after its field for a unit.
    pub unit_gutter: f32,
    pub dim: Color32,
    pub error: Color32,
    /// The mark of a required field.
    pub required: Color32,
}

/// The room the marks after a label take: the required `*` and the kit's help icon.
pub fn marks_width(required: bool, help: bool) -> f32 {
    MARK_GAP + if required { STAR_W } else { 0.0 } + if help { HELP_W } else { 0.0 }
}

/// Between a label and its marks.
const MARK_GAP: f32 = 2.0;
/// The required mark.
const STAR_W: f32 = 9.0;
/// The help icon and the space before it.
const HELP_W: f32 = 17.0;
/// A label longer than its column wraps onto this many lines (labels are short phrases).
pub const LABEL_LINES: usize = 2;

/// The font of a row's label.
pub fn label_font(ui: &Ui) -> egui::FontId {
    egui::TextStyle::Body.resolve(ui.style())
}

/// The width of `label` on one line, with its marks: what a label column needs to show it whole.
pub fn label_width(ui: &Ui, label: &str, required: bool, help: bool) -> f32 {
    crate::elide::width(ui, label, &label_font(ui)).ceil() + marks_width(required, help)
}

/// `label` laid out in a column `wrap` wide, on at most [`LABEL_LINES`] lines (elided past them).
pub fn label_galley(ui: &Ui, label: &str, color: Color32, wrap: f32) -> std::sync::Arc<egui::Galley> {
    let mut job = egui::text::LayoutJob::simple(label.to_string(), label_font(ui), color, wrap.max(0.0));
    job.wrap.max_rows = LABEL_LINES;
    ui.painter().layout_job(job)
}

/// A row: `label` (with a mark when `required`, and the kit's help icon when there is `help`) in
/// the label column, in the kit's label colour, then `add`, a kit row high. A label wider than the
/// column wraps onto [`LABEL_LINES`] lines; the marks always follow it inside the column.
/// Hovering the label or the icon shows the label and `help`.
pub fn row<R>(ui: &mut Ui, s: &Style, label: &str, required: bool, help: &str, add: impl FnOnce(&mut Ui) -> R) -> R {
    ui.horizontal(|ui| {
        let ink = solidcraft_theme::Tokens::get().dialog_label;
        let galley = label_galley(ui, label, ink, s.label_width - marks_width(required, !help.is_empty()));
        let (rect, response) = ui.allocate_exact_size(vec2(s.label_width, galley.size().y.max(kit::row_height())), egui::Sense::hover());
        let text_rect = egui::Align2::LEFT_CENTER.anchor_size(rect.left_center(), galley.size());
        let tip = if help.is_empty() { label.to_string() } else { format!("{label}\n{help}") };
        response.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Label, true, label));
        response.backing_on_hover_text(&tip);
        // The marks sit after the label's last line.
        let last = galley.rows.last().map_or(egui::Rect::ZERO, |r| r.rect());
        let mark_y = text_rect.top() + last.center().y;
        let mut x = text_rect.left() + last.right() + MARK_GAP;
        ui.painter().galley(text_rect.min, galley, ink);
        if required {
            let star = egui::Rect::from_min_size(egui::pos2(x, mark_y - 7.0), vec2(STAR_W, 14.0));
            ui.painter().text(egui::pos2(x, mark_y), egui::Align2::LEFT_CENTER, "*", label_font(ui), s.required);
            ui.interact(star, ui.scope_id().with("required"), egui::Sense::hover()).backing_on_hover_text("Required");
            x += STAR_W;
        }
        if !help.is_empty() {
            let mark = egui::Rect::from_center_size(egui::pos2(x + HELP_W / 2.0, mark_y), vec2(HELP_W - MARK_GAP, HELP_W - MARK_GAP));
            // In the label column already allocated: a child of its own, so the row's field keeps
            // its place.
            let mut mark_ui =
                ui.new_child(egui::UiBuilder::new().max_rect(mark).layout(egui::Layout::centered_and_justified(egui::Direction::LeftToRight)));
            kit::help_icon(&mut mark_ui, &tip);
        }
        add(ui)
    })
    .inner
}

/// The width of `unit` after a field ([`unit`]), with the space before it; 0 without one.
pub fn unit_width(ui: &Ui, unit: Option<&str>) -> f32 {
    unit.map_or(0.0, |u| crate::elide::width(ui, u, &egui::TextStyle::Small.resolve(ui.style())).ceil() + ui.spacing().item_spacing.x)
}

/// A message under a field, aligned with the fields.
pub fn note(ui: &mut Ui, s: &Style, text: &str, color: Color32) {
    ui.horizontal(|ui| {
        ui.add_space(s.label_width + ui.spacing().item_spacing.x);
        ui.add(egui::Label::new(RichText::new(text).small().color(color)).wrap());
    });
}

/// What a number field did this frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NumberEdit {
    Unchanged,
    /// A new value; `None` when the field was emptied.
    Set(Option<f64>),
    /// What is typed is not a number (yet); the value is unchanged.
    Invalid,
}

/// The typed text of a field and the value it was typed over: kept while the value stays the
/// same, so a half-typed `1e-` survives frames and an outside change (undo, the JSON view)
/// replaces it.
#[derive(Clone, Debug, Default)]
struct Draft {
    text: String,
    over: Option<f64>,
}

/// A number field: shows `value` (or the typed text), greys `hint` when empty, and reports a new
/// value only when the text is a finite number (or empty). `integer` refuses fractions. `id` is the
/// field's own (unique in the form, and differing between forms): it keys both the typed text and
/// the widget, so focus and typing never move to another form's field at the same place.
pub fn number(ui: &mut Ui, id: Id, value: Option<f64>, hint: &str, integer: bool, width: f32) -> NumberEdit {
    let shown = |v: Option<f64>| v.map(format_number).unwrap_or_default();
    let mut draft = ui.data(|d| d.get_temp::<Draft>(id)).filter(|d| d.over == value).unwrap_or_else(|| Draft { text: shown(value), over: value });
    let parse = |t: &str| -> Option<Option<f64>> {
        let t = t.trim();
        if t.is_empty() {
            return Some(None);
        }
        let x = t.parse::<f64>().ok().filter(|x| x.is_finite())?;
        (!integer || x.fract() == 0.0).then_some(Some(x))
    };
    let bad = parse(&draft.text).is_none();
    // The kit's editor, keyed by `id`, in the error colours while the text is not a number.
    let r = sized(ui, width, |ui| kit::cell_editor(ui, id, &mut draft.text, hint, bad));
    let out = if r.changed() {
        match parse(&draft.text) {
            Some(v) => {
                draft.over = v;
                NumberEdit::Set(v)
            }
            None => NumberEdit::Invalid,
        }
    } else if bad {
        NumberEdit::Invalid
    } else {
        NumberEdit::Unchanged
    };
    ui.data_mut(|d| d.insert_temp(id, draft));
    out
}

/// A number as typed: integers without a fraction, others in their shortest round-trip form.
pub fn format_number(x: f64) -> String {
    if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        serde_json::Number::from_f64(x).map_or_else(|| x.to_string(), |n| n.to_string())
    }
}

/// `add` in a box `width` wide and a kit row high, left to right.
fn sized<R>(ui: &mut Ui, width: f32, add: impl FnOnce(&mut Ui) -> R) -> R {
    ui.allocate_ui_with_layout(vec2(width.max(0.0), kit::row_height()), egui::Layout::left_to_right(egui::Align::Center), add).inner
}

/// A drop-down of `options` (labels), the kit's; returns the index picked this frame.
pub fn combo(ui: &mut Ui, id: Id, selected: &str, options: &[&str], width: f32) -> Option<usize> {
    let mut picked = None;
    kit::bounded_dropdown_menu(ui, id, selected, width, |ui| {
        for (i, o) in options.iter().enumerate() {
            if ui.selectable_label(*o == selected, *o).clicked() {
                picked = Some(i);
            }
        }
    });
    picked
}

/// A one-line text field (the kit's editor) keyed by `id` (see [`number`]), in the error colours
/// while `bad`; its response.
pub fn checked_text(ui: &mut Ui, id: Id, value: &mut String, hint: &str, bad: bool, width: f32) -> egui::Response {
    sized(ui, width, |ui| kit::cell_editor(ui, id, value, hint, bad))
}

/// A one-line text field (the kit's editor) with a greyed hint, keyed by `id` (see [`number`]); true
/// when it changed.
pub fn text(ui: &mut Ui, id: Id, value: &mut String, hint: &str, width: f32) -> bool {
    checked_text(ui, id, value, hint, false, width).changed()
}

/// Labels as typed: runs of consecutive labels as `a-b`, the rest by commas ("1-5, 8, 10-12").
pub fn format_labels(labels: &[u64]) -> String {
    let mut out: Vec<String> = Vec::new();
    let mut i = 0;
    while let Some(&a) = labels.get(i) {
        let mut j = i;
        while labels.get(j + 1).is_some_and(|&n| Some(n) == labels.get(j).and_then(|m| m.checked_add(1))) {
            j += 1;
        }
        let b = labels.get(j).copied().unwrap_or(a);
        out.push(match j - i {
            0 => a.to_string(),
            1 => format!("{a}, {b}"),
            _ => format!("{a}-{b}"),
        });
        i = j + 1;
    }
    out.join(", ")
}

/// Labels typed as numbers and `a-b` ranges, separated by commas or spaces, in the order typed;
/// `None` when something is not a label from 1, or when there would be more than `max`.
pub fn parse_labels(text: &str, max: usize) -> Option<Vec<u64>> {
    let mut out = Vec::new();
    for part in text.split(|c: char| c == ',' || c.is_whitespace()).filter(|p| !p.is_empty()) {
        let (a, b) = match part.split_once('-') {
            Some((a, b)) => (a.parse::<u64>().ok()?, b.parse::<u64>().ok()?),
            None => {
                let n = part.parse::<u64>().ok()?;
                (n, n)
            }
        };
        if a == 0 || b < a || out.len() as u64 + (b - a) >= max as u64 {
            return None;
        }
        out.extend(a..=b);
    }
    Some(out)
}

/// Element faces or edges as typed: `element:side`, runs of consecutive elements with the same
/// side as `a-b:side` ("1-5:3, 8:1").
pub fn format_sides(sides: &[(u64, u64)]) -> String {
    let mut out: Vec<String> = Vec::new();
    let mut i = 0;
    while let Some(&(a, side)) = sides.get(i) {
        let mut j = i;
        while sides.get(j + 1).is_some_and(|&(n, s)| s == side && Some(n) == sides.get(j).and_then(|(m, _)| m.checked_add(1))) {
            j += 1;
        }
        let b = sides.get(j).map_or(a, |(n, _)| *n);
        out.push(match j - i {
            0 => format!("{a}:{side}"),
            1 => format!("{a}:{side}, {b}:{side}"),
            _ => format!("{a}-{b}:{side}"),
        });
        i = j + 1;
    }
    out.join(", ")
}

/// Element faces or edges typed as `element:side` or `a-b:side`, separated by commas or spaces;
/// `None` when an element is not a label from 1, a side not 1–255, or there would be more than
/// `max`.
pub fn parse_sides(text: &str, max: usize) -> Option<Vec<(u64, u64)>> {
    let mut out = Vec::new();
    for part in text.split(|c: char| c == ',' || c.is_whitespace()).filter(|p| !p.is_empty()) {
        let (elements, side) = part.split_once(':')?;
        let side = side.parse::<u64>().ok().filter(|s| (1..=u64::from(u8::MAX)).contains(s))?;
        let labels = parse_labels(elements, max.saturating_sub(out.len()))?;
        if labels.is_empty() {
            return None;
        }
        out.extend(labels.into_iter().map(|e| (e, side)));
    }
    Some(out)
}

/// A short label before a field ("X"), dim.
pub fn caption(ui: &mut Ui, s: &Style, text: &str) {
    ui.label(RichText::new(text).small().color(s.dim));
}

/// A unit after a field.
pub fn unit(ui: &mut Ui, s: &Style, unit: Option<&str>) {
    if let Some(u) = unit {
        ui.label(RichText::new(u).small().color(s.dim));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_print_and_parse_as_ranges() {
        assert_eq!(format_labels(&[1, 2, 3, 4, 5, 8, 10, 11]), "1-5, 8, 10, 11");
        assert_eq!(format_labels(&[]), "");
        assert_eq!(format_labels(&[7, 3, 4, 5]), "7, 3-5", "the order is kept");
        assert_eq!(parse_labels("1-5, 8 10,11", 100), Some(vec![1, 2, 3, 4, 5, 8, 10, 11]));
        assert_eq!(parse_labels("", 100), Some(vec![]));
        for bad in ["0", "x", "5-3", "1-", "-4", "1-1000"] {
            assert_eq!(parse_labels(bad, 100), None, "{bad}");
        }
        assert_eq!(parse_labels(&format_labels(&[3, 4, 5, 9]), 100), Some(vec![3, 4, 5, 9]));
    }

    #[test]
    fn sides_print_and_parse_as_ranges() {
        let s = [(1, 3), (2, 3), (3, 3), (5, 1), (6, 2), (7, 2)];
        assert_eq!(format_sides(&s), "1-3:3, 5:1, 6:2, 7:2");
        assert_eq!(parse_sides(&format_sides(&s), 100).as_deref(), Some(&s[..]));
        assert_eq!(parse_sides("4:1 9-10:6", 100), Some(vec![(4, 1), (9, 6), (10, 6)]));
        assert_eq!(parse_sides("", 100), Some(vec![]));
        for bad in ["4", "4:", ":3", "0:1", "4:0", "4:256", "1-200:1", "x:1"] {
            assert_eq!(parse_sides(bad, 100), None, "{bad}");
        }
    }

    #[test]
    fn numbers_print_as_typed() {
        assert_eq!(format_number(210000.0), "210000");
        assert_eq!(format_number(0.3), "0.3");
        assert_eq!(format_number(7.85e-9), "7.85e-9");
        assert_eq!(format_number(-2.0), "-2");
    }
}
