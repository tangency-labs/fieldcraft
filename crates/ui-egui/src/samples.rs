//! The samples as a grid of cards (the Start card, File › Samples) and as rows (Help › Samples):
//! every sample of `fieldcraft_engine::sample::SAMPLES`, the flagship ones first, each with its
//! thumbnail where one exists. Thumbnails are the gallery captures (docs/media/examples) cropped
//! and scaled to 320 px, embedded at build time so they work on the web too.

use egui::{Rect, Sense, TextureHandle, vec2};
use fieldcraft_engine::sample::{SAMPLES, Sample};

use crate::theme::Tokens;

/// The order the samples are offered in: these first, then any other in the registry's order.
const FLAGSHIP: &[&str] = &["crank-arm", "drop-test", "rubber-seal", "heated-bracket", "bracket-hex", "flow-past", "cantilever", "impact", "parts"];

/// The embedded thumbnails (PNG, 320 px wide), by sample name.
const THUMBNAILS: &[(&str, &[u8])] = &[
    ("crank-arm", include_bytes!("../thumbnails/crank-arm.png")),
    ("crank-mirror", include_bytes!("../thumbnails/crank-mirror.png")),
    ("drop-test", include_bytes!("../thumbnails/drop-test.png")),
    ("rubber-seal", include_bytes!("../thumbnails/rubber-seal.png")),
    ("heated-bracket", include_bytes!("../thumbnails/heated-bracket.png")),
    ("bracket-hex", include_bytes!("../thumbnails/bracket-hex.png")),
    ("flow-past", include_bytes!("../thumbnails/flow-past.png")),
    ("heated-flow", include_bytes!("../thumbnails/heated-flow.png")),
    ("cooling-heat-sink", include_bytes!("../thumbnails/cooling-heat-sink.png")),
    ("cooling-enclosure", include_bytes!("../thumbnails/cooling-enclosure.png")),
    ("cooling-aio", include_bytes!("../thumbnails/cooling-aio.png")),
    ("cantilever", include_bytes!("../thumbnails/cantilever.png")),
    ("impact", include_bytes!("../thumbnails/impact.png")),
    ("parts", include_bytes!("../thumbnails/parts.png")),
    ("ibeam", include_bytes!("../thumbnails/ibeam.png")),
];

/// The thumbnails' aspect ratio (height over width).
const THUMB_ASPECT: f32 = 213.0 / 320.0;
/// The narrowest card; the grid fits as many columns as this allows.
const MIN_CARD_W: f32 = 180.0;
const MAX_COLUMNS: usize = 4;
const GAP: f32 = 10.0;

/// Every sample, flagship first.
pub fn ordered() -> Vec<&'static Sample> {
    let rank = |s: &Sample| FLAGSHIP.iter().position(|n| *n == s.name).unwrap_or(FLAGSHIP.len());
    let mut out: Vec<&'static Sample> = SAMPLES.iter().collect();
    // Stable: samples outside the flagship list keep the registry's order.
    out.sort_by_key(|s| rank(s));
    out
}

/// The id of a sample's card in the grid (tests and the control channel read its response).
pub fn card_id(name: &str) -> egui::Id {
    egui::Id::new(("fc_sample_card", name))
}

/// A PNG decoded to an egui image; `None` when it does not decode.
fn decode(bytes: &[u8]) -> Option<egui::ColorImage> {
    let mut decoder = png::Decoder::new(std::io::Cursor::new(bytes));
    decoder.set_transformations(png::Transformations::normalize_to_color8());
    let mut reader = decoder.read_info().ok()?;
    let mut buf = vec![0; reader.output_buffer_size()?];
    let info = reader.next_frame(&mut buf).ok()?;
    let size = [info.width as usize, info.height as usize];
    let pixels = buf.get(..info.buffer_size())?;
    match info.color_type {
        png::ColorType::Rgba => Some(egui::ColorImage::from_rgba_unmultiplied(size, pixels)),
        png::ColorType::Rgb => Some(egui::ColorImage::from_rgb(size, pixels)),
        _ => None,
    }
}

/// The thumbnail of sample `name`, uploaded once; `None` when it has none.
fn thumbnail(ctx: &egui::Context, name: &'static str) -> Option<TextureHandle> {
    let key = egui::Id::new(("fc_sample_thumbnail", name));
    if let Some(cached) = ctx.data(|d| d.get_temp::<Option<TextureHandle>>(key)) {
        return cached;
    }
    let texture = THUMBNAILS.iter().find(|(n, _)| *n == name).and_then(|(_, bytes)| {
        let image = decode(bytes);
        if image.is_none() {
            log::warn!("sample thumbnail {name}: does not decode");
        }
        image.map(|i| ctx.load_texture(format!("sample:{name}"), i, egui::TextureOptions::LINEAR))
    });
    ctx.data_mut(|d| d.insert_temp(key, texture.clone()));
    texture
}

/// The samples as a grid of cards in a vertical scroll area with a solid (always drawn)
/// scrollbar; `open` gets the one clicked.
pub fn grid(ui: &mut egui::Ui, t: &Tokens, max_height: f32, open: &mut Option<&'static str>) {
    ui.spacing_mut().scroll = egui::style::ScrollStyle::solid();
    egui::ScrollArea::vertical().id_salt("fc_sample_grid").max_height(max_height).auto_shrink([false, true]).show(ui, |ui| {
        let width = ui.available_width();
        let columns = (((width + GAP) / (MIN_CARD_W + GAP)).floor() as usize).clamp(1, MAX_COLUMNS);
        let card_w = ((width - GAP * (columns as f32 - 1.0)) / columns as f32).max(1.0);
        for row in ordered().chunks(columns) {
            ui.horizontal_top(|ui| {
                ui.spacing_mut().item_spacing.x = GAP;
                for sample in row {
                    if card(ui, t, sample, card_w) {
                        *open = Some(sample.name);
                    }
                }
            });
            ui.add_space(GAP);
        }
    });
}

/// One sample's card: thumbnail, name and a one-line description (in full on hover). True when
/// clicked.
fn card(ui: &mut egui::Ui, t: &Tokens, sample: &'static Sample, w: f32) -> bool {
    let thumb_h = (w * THUMB_ASPECT).round();
    let text_h = 40.0;
    let (rect, _) = ui.allocate_exact_size(vec2(w, thumb_h + text_h), Sense::hover());
    let response = ui.interact(rect, card_id(sample.name), Sense::click()).on_hover_text(sample.about);
    let painter = ui.painter_at(rect);
    let fill = if response.hovered() { t.hover } else { t.panel };
    painter.rect_filled(rect, 6.0, fill);
    let image_rect = Rect::from_min_size(rect.min, vec2(w, thumb_h)).shrink(4.0);
    match thumbnail(ui.ctx(), sample.name) {
        Some(texture) => {
            egui::Image::new(&texture).corner_radius(4.0).paint_at(ui, image_rect);
        }
        None => {
            // No capture yet: the name's initial on a plain tile.
            painter.rect_filled(image_rect, 4.0, t.accent_soft);
            let initial = sample.label.chars().next().map(String::from).unwrap_or_default();
            painter.text(image_rect.center(), egui::Align2::CENTER_CENTER, initial, egui::FontId::proportional(thumb_h * 0.4), t.accent);
        }
    }
    let stroke = if response.hovered() { t.accent } else { t.border };
    painter.rect_stroke(rect, 6.0, egui::Stroke::new(1.0, stroke), egui::StrokeKind::Inside);
    let text_rect = Rect::from_min_max(egui::pos2(rect.left() + 6.0, rect.top() + thumb_h + 2.0), rect.max - vec2(6.0, 2.0));
    let mut text_ui = ui.new_child(egui::UiBuilder::new().max_rect(text_rect).layout(egui::Layout::top_down(egui::Align::Min)));
    text_ui.add(egui::Label::new(egui::RichText::new(sample.label).strong().color(t.overlay_text)).truncate().selectable(false));
    text_ui.add(egui::Label::new(egui::RichText::new(sample.about).small().color(t.text_dim)).truncate().selectable(false));
    response.clicked()
}

/// The samples as rows (name, what it shows, Open) for a narrow menu; `open` gets the one clicked.
pub fn rows(ui: &mut egui::Ui, t: &Tokens, open: &mut Option<&'static str>) {
    let width = ui.available_width().clamp(280.0, 420.0);
    ui.spacing_mut().scroll = egui::style::ScrollStyle::solid();
    egui::ScrollArea::vertical().id_salt("fc_sample_rows").auto_shrink([false, true]).show(ui, |ui| {
        for sample in ordered() {
            ui.horizontal(|ui| {
                ui.vertical(|ui| {
                    ui.set_width(width - 70.0);
                    ui.label(egui::RichText::new(sample.label).strong());
                    ui.add(egui::Label::new(egui::RichText::new(sample.about).small().color(t.text_dim)).wrap());
                });
                if ui.button("Open").clicked() {
                    *open = Some(sample.name);
                }
            });
            ui.add_space(4.0);
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every registered sample is offered exactly once, flagship first.
    #[test]
    fn every_sample_is_offered_once_flagship_first() {
        let names: Vec<&str> = ordered().iter().map(|s| s.name).collect();
        assert_eq!(names.len(), SAMPLES.len());
        for s in SAMPLES {
            assert_eq!(names.iter().filter(|n| **n == s.name).count(), 1, "{}", s.name);
        }
        assert_eq!(names.first(), Some(&"crank-arm"));
    }

    /// Every flagship name and every thumbnail names a registered sample, every sample has a
    /// thumbnail (#294), and every thumbnail decodes at 320 px wide.
    #[test]
    fn thumbnails_and_flagship_names_are_registered_samples() {
        for name in FLAGSHIP.iter().chain(THUMBNAILS.iter().map(|(n, _)| n)) {
            assert!(fieldcraft_engine::sample::find(name).is_some(), "{name}");
        }
        for s in SAMPLES {
            assert!(THUMBNAILS.iter().any(|(n, _)| *n == s.name), "{} has no thumbnail", s.name);
        }
        for (name, bytes) in THUMBNAILS {
            let image = decode(bytes).unwrap_or_else(|| panic!("{name} does not decode"));
            assert_eq!(image.size[0], 320, "{name}");
        }
    }
}
