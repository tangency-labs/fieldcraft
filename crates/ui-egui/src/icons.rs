//! Toolbar icons: original artwork drawn with VectorCraft (`assets/icons`, `cargo xtask icons`),
//! embedded at build time and rasterised with resvg.
//!
//! Icons are drawn in four palette roles (line, accent, soft fill, warm), each in one fixed colour.
//! Before rasterising, every role colour is replaced by the theme token it stands for (`icon`,
//! `accent`, `icon_fill`, `warm`, as SolidCraft paints its icons), so the same SVG serves the dark
//! and light themes and the greyed "coming soon" state. Icons named `ui.*` belong to the window
//! itself (navigation bar, title bar), the others to catalogue entries. Textures are cached per
//! icon, theme, state and pixel size.
//!
//! The same artwork is the UI kit's icon set ([`KIT_SET`]), in the same roles and colours: the
//! kit's windows, confirmations and menus draw from it.

use std::collections::HashMap;

use egui::{Color32, TextureHandle};
use resvg::{tiny_skia, usvg};
use solidcraft_ui_kit::icons::{IconData, IconSet};

use crate::theme::Tokens;

include!(concat!(env!("OUT_DIR"), "/icons.rs"));

/// The palette roles' colours in the SVG files. Keep in step with `PALETTE` in
/// `xtask/src/icons.rs`.
pub const PALETTE: [&str; 4] = ["#dde2ea", "#4f9ce8", "#2d4e70", "#f0b44a"];

/// Largest icon side rasterised (pixels).
const MAX_PX: u32 = 256;

/// The icon's SVG text.
pub fn svg(id: &str) -> Option<&'static str> {
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

/// The palette roles as the UI kit names them, in [`PALETTE`]'s order, then the kit's guide role
/// (it draws dashed line strokes in it): line, accent (`active-3`: the kit keeps an icon's active
/// strokes in their colour when dashed), soft fill, warm (`active-2`, likewise).
const KIT_ROLES: [&str; 5] = ["line", "active-3", "passive-1", "active-2", "guide"];
/// The guide role's key colour: no artwork uses it; the kit draws dashed line strokes in it.
const KIT_GUIDE: &str = "#dde2eb";
const KIT_SOURCE: [&str; 5] = [PALETTE[0], PALETTE[1], PALETTE[2], PALETTE[3], KIT_GUIDE];

const fn kit_roles(t: Tokens) -> [[u8; 3]; 5] {
    const fn rgb(c: Color32) -> [u8; 3] {
        [c.r(), c.g(), c.b()]
    }
    // The guide is the line colour, so a dashed line looks as it does in the toolbar.
    [rgb(t.icon), rgb(t.accent), rgb(t.icon_fill), rgb(t.warm), rgb(t.icon)]
}

static KIT_LIGHT: [[u8; 3]; 5] = kit_roles(Tokens::light());
static KIT_DARK: [[u8; 3]; 5] = kit_roles(Tokens::dark());

/// The names the kit's own components draw (close buttons, disclosures, add buttons…) and the
/// FieldCraft icon each one is.
const KIT_ALIASES: &[(&str, &str, &str)] = &[
    ("close", "ui.close", ""),
    ("chevron", "ui.chevron", ""),
    ("chevron_right", "ui.chevron_right", ""),
    ("plus", "ui.plus", ""),
    ("home", "ui.home", ""),
    ("help", "ui.help", ""),
    ("zoom", "ui.zoom", ""),
    ("fit", "ui.fit", ""),
    ("eye", "ui.eye", ""),
    ("folder", "file.open", ""),
    ("bug_report", "diagnostics.bug_report", ""),
];

/// FieldCraft's icons as the UI kit's icon set, in the theme's colours: what the kit's windows,
/// menus and confirmations draw. [`install_kit`] makes it the kit's set.
pub static KIT_SET: IconSet = IconSet::new(IconData {
    role_names: &KIT_ROLES,
    role_source_hex: &KIT_SOURCE,
    palette: [&KIT_LIGHT, &KIT_DARK],
    accent_names: &[],
    accents: &[],
    accent_variants: &[],
    masters: KIT_MASTERS,
    aliases: KIT_ALIASES,
    icon_accents: &[],
    danger_accent: "",
});

/// Give the UI kit FieldCraft's icons (once per process; again does nothing).
pub fn install_kit() {
    solidcraft_ui_kit::icons::install(&KIT_SET);
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

    /// Every icon is a master of the kit's set that the kit reads (none falls back to its
    /// placeholder), and every name the kit's components draw resolves to one of them.
    #[test]
    fn the_kit_reads_every_icon_and_finds_every_name_it_draws() {
        assert_eq!(KIT_MASTERS.len(), ICONS.len());
        for (id, family, svg) in KIT_MASTERS {
            assert_eq!(*family, "", "{id}");
            if let Err(e) = KIT_SET.parse(svg) {
                panic!("{id}: the kit cannot read it: {e}");
            }
        }
        for (name, master, _) in KIT_ALIASES {
            assert!(svg(master).is_some(), "{name} names the missing icon {master}");
            let r = KIT_SET.resolve(name);
            assert!(!r.placeholder && r.master == *master, "{name} resolves to {r:?}");
        }
    }

    /// The kit draws an icon across its square as the toolbar does: the masters are scaled from
    /// the artwork's 24 units to the kit's 64.
    #[test]
    fn the_kit_draws_icons_at_their_full_size() {
        let look = solidcraft_ui_kit::icons::Look::themed(true, None, solidcraft_theme::ColorVision::Off);
        let doc = KIT_SET.doc_for("ui.close", 64).expect("the close icon");
        let pixmap = solidcraft_ui_kit::icons::render(doc, look, 64).expect("rendered");
        let inked: Vec<usize> = pixmap.data().chunks(4).enumerate().filter(|(_, p)| p[3] > 0).map(|(i, _)| i % 64).collect();
        let (lo, hi) = (inked.iter().min().copied().unwrap_or(64), inked.iter().max().copied().unwrap_or(0));
        // The cross runs from 6.5 to 17.5 of 24: about 17 to 47 of 64.
        assert!(lo < 20 && hi > 44, "the cross spans columns {lo}..{hi} of 64");
    }

    /// The kit draws each role in the token the toolbar draws it in, in both themes.
    #[test]
    fn the_kit_draws_icons_in_the_toolbar_colours() {
        for (dark, t) in [(false, Tokens::light()), (true, Tokens::dark())] {
            let look = solidcraft_ui_kit::icons::Look::themed(dark, None, solidcraft_theme::ColorVision::Off);
            for (role, colour) in KIT_ROLES.into_iter().zip(colours(&t, true).into_iter().chain([t.icon])) {
                assert_eq!(KIT_SET.role_colour(role, look), colour, "{role} dark={dark}");
            }
        }
    }
}
