//! The samples as the kit's sample cards (the Start card, File › Samples) and as rows (Help › Samples):
//! every sample of `fieldcraft_engine::sample::SAMPLES`, the flagship ones first, each with its
//! thumbnail where one exists. Thumbnails are the gallery captures (docs/media/examples) cropped
//! and scaled to 320 px, embedded at build time so they work on the web too.

use egui::TextureHandle;
use fieldcraft_engine::sample::{SAMPLES, Sample};
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _};

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

/// Every sample, flagship first.
pub fn ordered() -> Vec<&'static Sample> {
    let rank = |s: &Sample| FLAGSHIP.iter().position(|n| *n == s.name).unwrap_or(FLAGSHIP.len());
    let mut out: Vec<&'static Sample> = SAMPLES.iter().collect();
    // Stable: samples outside the flagship list keep the registry's order.
    out.sort_by_key(|s| rank(s));
    out
}

/// Where a sample's card was drawn in the grid, kept in the context's data under this id as its
/// rect and visible rect (tests and the control channel read them).
pub fn card_id(name: &str) -> egui::Id {
    egui::Id::unique(("fc_sample_card", name))
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
    let key = egui::Id::unique(("fc_sample_thumbnail", name));
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

/// The samples as the kit's sample cards (as SolidCraft's start page shows its own); `open` gets
/// the one clicked. With `max_height`, in a vertical scroll area that high with a solid (always
/// drawn) scrollbar; without, laid out in full (the caller scrolls). Each card's rect and visible
/// rect are kept under [`card_id`] (tests and the control channel read them).
pub fn grid(ui: &mut egui::Ui, max_height: Option<f32>, open: &mut Option<&'static str>) {
    match max_height {
        Some(height) => {
            ui.spacing_mut().scroll = egui::style::ScrollStyle::solid();
            egui::ScrollArea::vertical().id_salt("fc_sample_grid").max_height(height).auto_shrink([false, true]).show(ui, |ui| cards(ui, open));
        }
        None => cards(ui, open),
    }
}

fn cards(ui: &mut egui::Ui, open: &mut Option<&'static str>) {
    let samples = ordered();
    let textures: Vec<Option<TextureHandle>> = samples.iter().map(|s| thumbnail(ui.ctx(), s.name)).collect();
    let cards: Vec<kit::SampleCard> = samples
        .iter()
        .zip(&textures)
        .map(|(s, texture)| kit::SampleCard {
            title: s.label,
            description: s.about,
            icon: kit::Icon::Kit("ui.views"),
            thumbnail: texture.as_ref().map(|t| (t.id(), t.size_vec2())),
        })
        .collect();
    for (sample, response) in samples.iter().zip(kit::FileView::sample_grid(ui, "fc_samples", &cards)) {
        ui.ctx().data_mut(|d| d.insert_temp(card_id(sample.name), (response.rect, response.interact_rect)));
        if response.backing_on_hover_text(sample.about).clicked() {
            *open = Some(sample.name);
        }
    }
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
                if solidcraft_ui_kit::button(ui, "Open", solidcraft_ui_kit::Button::Small).clicked() {
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
