//! Toolbar icons: original artwork drawn with VectorCraft (`assets/icons`, `cargo xtask icons`),
//! embedded at build time and rasterised with resvg.
//!
//! Icons are drawn in four palette roles (line, accent, soft fill, warm), each in one fixed colour.
//! Before rasterising, every role colour is replaced by the theme token it stands for (`icon`,
//! `accent`, `icon_fill`, `warm`, as SolidCraft paints its icons), so the same SVG serves the dark
//! and light themes and the greyed "coming soon" state. Icons named `ui.*` belong to the window
//! itself (navigation bar, title bar), the others to catalogue entries. Textures are cached per
//! icon, theme, state and pixel size.

use std::collections::HashMap;

use egui::{Color32, TextureHandle};
use resvg::{tiny_skia, usvg};

use crate::theme::Tokens;

include!(concat!(env!("OUT_DIR"), "/icons.rs"));

/// The palette roles' colours in the SVG files. Keep in step with `PALETTE` in
/// `xtask/src/icons.rs`.
pub const PALETTE: [&str; 4] = ["#dde2ea", "#4f9ce8", "#2d4e70", "#f0b44a"];

/// Largest icon side rasterised (pixels).
const MAX_PX: u32 = 256;

/// The icon's SVG text.
pub fn svg(id: &str) -> Option<&'static str> {
    // Fatigue shares the existing original VectorCraft signal/curve artwork.
    // These are semantic aliases, not external assets or a generic missing-icon fallback.
    let id = match id {
        "solve.fatigue.history" | "solve.fatigue.transient" => "model.amplitude",
        "solve.fatigue.strain_cycles" => "results.xy",
        other => other,
    };
    ICONS.binary_search_by(|(k, _)| (*k).cmp(id)).ok().and_then(|i| ICONS.get(i)).map(|(_, s)| *s)
}

/// The colours the palette roles (line, accent, soft, warm) take in a theme, greyed when the
/// entry is not available yet.
pub fn colours(t: &Tokens, enabled: bool) -> [Color32; 4] {
    if enabled { [t.icon, t.accent, t.icon_fill, t.warm] } else { [t.border, t.border, t.panel_header, t.border] }
}

/// The SVG with each palette colour replaced by its role's colour.
pub fn recolour(svg: &str, colours: [Color32; 4]) -> String {
    let mut out = svg.to_string();
    for (from, to) in PALETTE.iter().zip(colours) {
        let hex = format!("#{:02x}{:02x}{:02x}", to.r(), to.g(), to.b());
        out = out.replace(&format!("\"{from}\""), &format!("\"{hex}\""));
    }
    out
}

/// Rasterise an SVG to a `px` × `px` image.
pub fn rasterise(svg: &str, px: u32) -> Result<egui::ColorImage, String> {
    let px = px.clamp(1, MAX_PX);
    let tree = usvg::Tree::from_str(svg, &usvg::Options::default()).map_err(|e| e.to_string())?;
    let mut pixmap = tiny_skia::Pixmap::new(px, px).ok_or("empty icon")?;
    let size = tree.size();
    let scale = px as f32 / size.width().max(size.height()).max(f32::EPSILON);
    resvg::render(&tree, tiny_skia::Transform::from_scale(scale, scale), &mut pixmap.as_mut());
    Ok(egui::ColorImage::from_rgba_premultiplied([px as usize, px as usize], pixmap.data()))
}

/// How an icon is drawn.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Look {
    /// In the theme's colours.
    Normal,
    /// Greyed: coming soon or unavailable.
    Greyed,
    /// On the title bar, which is dark in both themes: the dark theme's colours.
    AppBar,
}

impl Look {
    pub fn of(enabled: bool) -> Look {
        if enabled { Look::Normal } else { Look::Greyed }
    }

    fn colours(self, dark: bool) -> [Color32; 4] {
        match self {
            Look::Normal => colours(&Tokens::get(dark), true),
            Look::Greyed => colours(&Tokens::get(dark), false),
            Look::AppBar => colours(&Tokens::dark(), true),
        }
    }
}

/// Icon textures by (id, look, pixel size) for one theme; `None` when the icon is missing or
/// does not rasterise, so it is not tried again every frame. The window draws icons at several
/// sizes at once (ribbon, browser, menus), so each size keeps its own texture; a change of theme
/// drops them all.
#[derive(Default)]
pub struct IconCache {
    dark: Option<bool>,
    textures: HashMap<(&'static str, Look, u32), Option<TextureHandle>>,
}

impl IconCache {
    pub fn get(&mut self, ctx: &egui::Context, id: &'static str, dark: bool, look: Look, px: u32) -> Option<TextureHandle> {
        if self.dark != Some(dark) {
            self.textures.clear();
            self.dark = Some(dark);
        }
        self.textures
            .entry((id, look, px))
            .or_insert_with(|| {
                let svg = svg(id)?;
                let image = rasterise(&recolour(svg, look.colours(dark)), px).map_err(|e| log::warn!("icon {id}: {e}")).ok()?;
                Some(ctx.load_texture(format!("icon:{id}:{dark}:{look:?}:{px}"), image, egui::TextureOptions::LINEAR))
            })
            .clone()
    }
}

/// Paint icon `id` into `rect` (a square); a missing icon draws a thin box so nothing is blank.
pub fn paint(app: &mut crate::ModelApp, ui: &egui::Ui, id: &'static str, rect: egui::Rect, look: Look) {
    paint_tinted(app, ui, id, rect, look, Color32::WHITE);
}

/// [`paint`] multiplied by `tint` (a translucent white dims it).
pub fn paint_tinted(app: &mut crate::ModelApp, ui: &egui::Ui, id: &'static str, rect: egui::Rect, look: Look, tint: Color32) {
    let ctx = ui.ctx();
    let px = (rect.width().max(rect.height()) * ctx.pixels_per_point()).round().clamp(1.0, MAX_PX as f32) as u32;
    let p = ui.painter();
    match app.icons.get(ctx, id, app.ui.dark, look, px) {
        Some(tex) => {
            p.image(tex.id(), rect, egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)), tint);
        }
        None => {
            let t = Tokens::get(app.ui.dark);
            p.rect_stroke(rect.shrink(3.0), 3.0, egui::Stroke::new(1.0, t.text_dim), egui::StrokeKind::Inside);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_engine::catalogue::CATALOGUE;

    #[test]
    fn every_catalogue_entry_has_an_icon_and_every_icon_an_entry() {
        for c in CATALOGUE {
            assert!(svg(c.id).is_some(), "no icon assets/icons/{}.svg", c.id);
        }
        for (id, _) in ICONS {
            assert!(id.starts_with("ui.") || CATALOGUE.iter().any(|c| c.id == *id), "assets/icons/{id}.svg belongs to no catalogue entry");
        }
        assert!(ICONS.windows(2).all(|w| w[0].0 < w[1].0), "sorted for the binary search");
    }

    /// Icons use only the palette, so recolouring reaches every colour, and they rasterise.
    #[test]
    fn icons_use_only_palette_colours_and_rasterise() {
        for (id, text) in ICONS {
            for part in text.split('"').skip(1).step_by(2) {
                if part.starts_with('#') {
                    assert!(PALETTE.contains(&part), "{id}: colour {part} is not in the palette");
                }
            }
            let image = rasterise(text, 24).unwrap();
            assert!(image.pixels.iter().any(|p| p.a() > 0), "{id} draws nothing");
        }
    }

    /// Each role takes its own token, so recolouring never merges two roles into one colour.
    #[test]
    fn every_role_takes_a_distinct_theme_colour() {
        for t in [Tokens::dark(), Tokens::light()] {
            let c = colours(&t, true);
            for i in 0..c.len() {
                for j in i + 1..c.len() {
                    assert_ne!(c[i], c[j], "roles {i} and {j}");
                }
            }
        }
    }

    #[test]
    fn recolouring_follows_the_theme_and_greys_planned_entries() {
        let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M0 0 L24 24" stroke="#4f9ce8" fill="#2d4e70"/></svg>"##;
        let light = Tokens::light();
        let out = recolour(svg, colours(&light, true));
        let hex = |c: Color32| format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b());
        assert!(out.contains(&format!("stroke=\"{}\"", hex(light.accent))) && out.contains(&format!("fill=\"{}\"", hex(light.icon_fill))), "{out}");
        let grey = recolour(svg, colours(&light, false));
        let dim = light.border;
        assert!(grey.contains(&format!("stroke=\"#{:02x}{:02x}{:02x}\"", dim.r(), dim.g(), dim.b())), "{grey}");
        assert!(rasterise("not svg", 24).is_err());
        assert_eq!(rasterise(svg, 100_000).unwrap().size, [MAX_PX as usize; 2]);
    }
}
