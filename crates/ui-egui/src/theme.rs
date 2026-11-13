//! Colours and widget style. The interface draws with SolidCraft's UI kit, whose controls take
//! their colours from `solidcraft_theme::Tokens`: FieldCraft gives the kit SolidCraft's dark and
//! light chrome with FieldCraft's own teal accent ([`kit`]), so the apps look alike side by side
//! and are told apart, and installs them with the kit's theme state ([`apply`]). [`Tokens`] are
//! FieldCraft's: the chrome read from the kit's tokens, and what is FieldCraft's own (the
//! viewport, elements, results, setup symbols and selection). Every colour the UI uses comes from
//! one or the other.
//!
//! The token structure was adapted from SolidCraft's `crates/ui-egui/src/theme.rs` (MIT OR
//! Apache-2.0), see ATTRIBUTION.md.

use egui::Color32;
use fieldcraft_engine::view::Palette;
use fieldcraft_model::ElementKind;
use fieldcraft_render::Rgb;
use fieldcraft_render::glyphs::GlyphColors;
use solidcraft_theme as kit_theme;
use solidcraft_ui_kit::icons as kit_icons;
use solidcraft_ui_kit::type_face;

/// The accent in the light theme: teal #0F7F78 (OKLab hue 188°; white on it 4.85:1).
const LIGHT_ACCENT: Color32 = Color32::from_rgb(0x0f, 0x7f, 0x78);
/// Text in the accent on light grounds: the accent darkened (hue kept), #0B6A64, so it reads at
/// 4.5:1 on every light ground, the panel header (#ECEEF2) being the closest.
const LIGHT_ACCENT_INK: Color32 = Color32::from_rgb(0x0b, 0x6a, 0x64);
/// The accent in the dark theme: teal #0E837A, #14958B darkened (hue kept) so white
/// on it reaches 4.5:1 (#14958B gives 3.68:1).
const DARK_ACCENT: Color32 = Color32::from_rgb(0x0e, 0x83, 0x7a);
/// Text in the accent on dark grounds: a light teal, #5CD6CB, at least 7:1 on every dark ground.
const DARK_ACCENT_INK: Color32 = Color32::from_rgb(0x5c, 0xd6, 0xcb);
/// How far a hovered accent fill moves toward black: darker keeps white on it legible.
const ACCENT_HOVER: f32 = 0.12;
/// How far the selected document tab moves from the title bar toward white: the dark theme's
/// step, (53 − 28) / (255 − 28). The bar is dark in both themes, and the light theme's own tab
/// colour is almost the bar's.
const TAB_ON_BAR: f32 = 0.11;

/// `a` moved the fraction `t` of the way to `b`, channel by channel: `Color32::lerp_to_gamma` for
/// opaque colours, as a `const fn` so the themes stay constants.
const fn mix(a: Color32, b: Color32, t: f32) -> Color32 {
    const fn channel(a: u8, b: u8, t: f32) -> u8 {
        let v = a as f32 + (b as f32 - a as f32) * t;
        // In range for t in 0..=1; `as` saturates otherwise.
        (v + 0.5) as u8
    }
    Color32::from_rgb(channel(a.r(), b.r(), t), channel(a.g(), b.g(), t), channel(a.b(), b.b(), t))
}

/// FieldCraft's accent and the colours made from it, on SolidCraft's chrome of one theme.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Accent {
    /// The primary action, the current tab, focus.
    fill: Color32,
    /// Text in the accent on a panel or form.
    ink: Color32,
    /// Text and marks on an accent fill.
    on: Color32,
    /// An accent fill under the pointer.
    hover: Color32,
    /// The accent faintly over the panel: selected rows and pressed controls.
    soft: Color32,
    /// A hovered row or control: the accent more faintly still.
    row_hover: Color32,
}

const fn accent(dark: bool) -> Accent {
    let chrome = if dark { kit_theme::Tokens::dark() } else { kit_theme::Tokens::light() };
    let fill = if dark { DARK_ACCENT } else { LIGHT_ACCENT };
    Accent {
        fill,
        ink: if dark { DARK_ACCENT_INK } else { LIGHT_ACCENT_INK },
        on: Color32::WHITE,
        hover: mix(fill, Color32::BLACK, ACCENT_HOVER),
        soft: if dark { mix(chrome.panel, fill, 0.30) } else { mix(chrome.panel, fill, 0.18) },
        row_hover: if dark { mix(chrome.panel_header, fill, 0.10) } else { mix(chrome.panel, fill, 0.09) },
    }
}

/// The tokens the kit draws with: SolidCraft's chrome for the theme, with FieldCraft's accent in
/// every slot that carries an accent (the primary action, selection, focus, links, the start page
/// and the accents SolidCraft gives its own workspaces), so no SolidCraft blue is left.
pub const fn kit(dark: bool) -> kit_theme::Tokens {
    let mut k = if dark { kit_theme::Tokens::dark() } else { kit_theme::Tokens::light() };
    let a = accent(dark);
    k.accent = a.fill;
    k.accent_soft = a.soft;
    k.hover = a.row_hover;
    k.on_accent = a.on;
    k.pref_on_accent = a.on;
    k.selected_ink = a.ink;
    k.selected_bg = a.soft;
    k.form_primary = a.ink;
    k.form_primary_fill = a.fill;
    k.figure_accent = a.ink;
    k.file_primary = a.soft;
    k.file_primary_ink = a.ink;
    k.file_primary_border = a.fill;
    k.file_rail_hover = a.hover;
    // SolidCraft's workspaces each have an accent; FieldCraft's tabs all take its own.
    k.sketch_accent = a.fill;
    k.sketch_ink = a.ink;
    k.sketch_soft = a.soft;
    k.surface_accent = a.fill;
    k.surface_fill = a.fill;
    k.surface_ink = a.ink;
    k.surface_soft = a.soft;
    k.mesh_accent = a.fill;
    k.mesh_fill = a.fill;
    k.mesh_ink = a.ink;
    k.mesh_soft = a.soft;
    k.form_accent = a.fill;
    k.form_fill = a.fill;
    k.form_tab_ink = a.ink;
    k.form_soft = a.soft;
    k.sheet_fill = a.fill;
    k.sheet_ink = a.ink;
    k.sheet_soft = a.soft;
    k.render_accent = a.fill;
    k.render_fill = a.fill;
    k.render_ink = a.ink;
    k.render_soft = a.soft;
    if !dark {
        k.app_bar_tab_active = mix(k.app_bar, Color32::WHITE, TAB_ON_BAR);
    }
    k
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tokens {
    /// The title bar (dark in both themes).
    pub app_bar: Color32,
    /// The ribbon.
    pub toolbar: Color32,
    pub panel: Color32,
    pub panel_header: Color32,
    /// Command forms and the context panel.
    pub dialog_bg: Color32,
    /// Band behind a panel's heading and, lighter, behind its section headings.
    pub heading_bg: Color32,
    pub subheading_bg: Color32,
    pub border: Color32,
    pub field: Color32,
    pub text: Color32,
    pub text_dim: Color32,
    /// FieldCraft's accent (teal): the primary action, the current tab, focus and links. Teal
    /// keeps green, amber and red free for ok, warning and error.
    pub accent: Color32,
    /// Text in the accent on a panel or form (links, the next step, required marks): the accent
    /// where it reads at 4.5:1, lighter where it does not (the dark theme).
    pub accent_ink: Color32,
    /// Text and marks drawn on an accent fill (at least 4.5:1 against it, WCAG 2.x).
    pub on_accent: Color32,
    /// An accent fill under the pointer.
    pub accent_hover: Color32,
    /// The accent faintly over the panel: selected rows and pressed controls.
    pub accent_soft: Color32,
    /// A hovered row or control: the accent more faintly still.
    pub hover: Color32,
    pub error: Color32,
    pub ok: Color32,
    /// Warnings, loads and heat in icons.
    pub warm: Color32,
    /// Icon line work.
    pub icon: Color32,
    /// The journal timeline: its strip, its items and the current-position marker.
    pub timeline: Color32,
    pub timeline_item: Color32,
    pub timeline_marker: Color32,
    // ---- viewport ----
    /// Background gradient, top and bottom (equal for a flat fill).
    pub viewport_top: Color32,
    pub viewport_bottom: Color32,
    /// Chips drawn over the viewport (navigation bar, legend, start screen).
    pub overlay: Color32,
    pub overlay_text: Color32,
    pub axis_x: Color32,
    pub axis_y: Color32,
    pub axis_z: Color32,
    // ---- elements ----
    pub hex: Color32,
    pub tet: Color32,
    pub wedge: Color32,
    pub pyramid: Color32,
    pub shell: Color32,
    pub line_element: Color32,
    /// Element edges (the mesh wireframe).
    pub mesh_edge: Color32,
    /// Streamlines when the legend shows another field than velocity: stands out from the
    /// viewport background.
    pub streamline: Color32,
    pub node: Color32,
    // ---- geometry ----
    /// Faces of the open design and imported surfaces (SolidCraft's body colour).
    pub geometry: Color32,
    /// Their face boundaries.
    pub geometry_edge: Color32,
    /// The bodies of a design with several (the cells of a partitioned part), in turn.
    pub bodies: [Color32; 8],
    // ---- results ----
    /// The contour colour scale, low to high (the same in both themes).
    pub contour: [Color32; 5],
    /// Faces with no result at a node.
    pub no_result: Color32,
    /// Values above a clamped legend's maximum (the same in both themes).
    pub contour_over: Color32,
    // ---- analysis setup symbols ----
    /// Boundary conditions and rigid walls.
    pub glyph_constraint: Color32,
    /// Forces, moments, pressures, body loads and initial velocities.
    pub glyph_load: Color32,
    /// Reference points and couplings.
    pub glyph_interaction: Color32,
    /// Contact and tie surface symbols.
    pub glyph_contact: Color32,
    /// Temperatures, heat fluxes, convection and heat sources.
    pub glyph_thermal: Color32,
    /// Datum planes and axes.
    pub glyph_datum: Color32,
    // ---- selection ----
    /// What is selected in the viewport, and what a click would pick under the pointer.
    pub pick_selected: Color32,
    pub pick_hover: Color32,
}

/// The contour colour scale: blue, cyan, green, yellow, red.
const CONTOUR: [Color32; 5] = [
    Color32::from_rgb(0x2b, 0x4f, 0xd8),
    Color32::from_rgb(0x2f, 0xb8, 0xe0),
    Color32::from_rgb(0x3c, 0xc8, 0x64),
    Color32::from_rgb(0xf0, 0xd2, 0x3c),
    Color32::from_rgb(0xe0, 0x40, 0x2e),
];

impl Tokens {
    /// Both themes, evaluated at compile time.
    const LIGHT: Tokens = Tokens::light();
    const DARK: Tokens = Tokens::dark();

    pub fn get(dark: bool) -> Tokens {
        if dark { Tokens::DARK } else { Tokens::LIGHT }
    }

    /// These tokens with the chrome of a theme: the kit's ([`kit`]), with FieldCraft's accent.
    const fn with_chrome(mut self, dark: bool) -> Tokens {
        let (k, a) = (kit(dark), accent(dark));
        self.app_bar = k.app_bar;
        self.toolbar = k.toolbar;
        self.panel = k.panel;
        self.panel_header = k.panel_header;
        self.dialog_bg = k.dialog_body;
        self.border = k.border;
        self.field = k.control_fill;
        self.text = k.text;
        self.text_dim = k.text_dim;
        self.accent = a.fill;
        self.accent_ink = a.ink;
        self.on_accent = a.on;
        self.accent_hover = a.hover;
        self.accent_soft = a.soft;
        self.hover = a.row_hover;
        self.error = k.error;
        self.warm = k.warning;
        self.icon = k.icon;
        self.timeline = k.timeline;
        self.timeline_item = k.timeline_item;
        self.timeline_marker = k.timeline_marker;
        self.viewport_top = k.viewport_top;
        self.viewport_bottom = k.viewport_bottom;
        self.overlay = k.overlay;
        self.overlay_text = k.text;
        self
    }

    /// The light theme: the kit's light chrome and FieldCraft's own colours.
    pub const fn light() -> Tokens {
        // The chrome is filled in by `with_chrome`.
        let chrome = Color32::PLACEHOLDER;
        Tokens {
            app_bar: chrome,
            toolbar: chrome,
            panel: chrome,
            panel_header: chrome,
            dialog_bg: chrome,
            heading_bg: Color32::from_rgb(218, 223, 231),
            subheading_bg: Color32::from_rgb(234, 237, 242),
            border: chrome,
            field: chrome,
            text: chrome,
            text_dim: chrome,
            accent: chrome,
            accent_ink: chrome,
            on_accent: chrome,
            accent_hover: chrome,
            accent_soft: chrome,
            hover: chrome,
            error: chrome,
            // Darker than the kit's green so the teal accent stays apart from it for tritan
            // vision (ΔE00 13 rather than 5).
            ok: Color32::from_rgb(0x23, 0x5f, 0x2e),
            warm: chrome,
            icon: chrome,
            timeline: chrome,
            timeline_item: chrome,
            timeline_marker: chrome,
            viewport_top: chrome,
            viewport_bottom: chrome,
            overlay: chrome,
            overlay_text: chrome,
            axis_x: Color32::from_rgb(0xc8, 0x3c, 0x32),
            axis_y: Color32::from_rgb(0x3c, 0x96, 0x32),
            axis_z: Color32::from_rgb(0x2c, 0x64, 0xc8),
            hex: Color32::from_rgb(0x78, 0x9c, 0xc0),
            tet: Color32::from_rgb(0x80, 0xb4, 0x8a),
            wedge: Color32::from_rgb(0xc0, 0xaa, 0x78),
            pyramid: Color32::from_rgb(0x8c, 0x80, 0xb0),
            shell: Color32::from_rgb(0xc0, 0x8c, 0x78),
            line_element: Color32::from_rgb(0xc8, 0x8c, 0x1e),
            mesh_edge: Color32::from_rgb(0x1e, 0x23, 0x2a),
            streamline: Color32::from_rgb(0x28, 0x2e, 0x38),
            node: Color32::from_rgb(0xb4, 0x3c, 0x8c),
            geometry: Color32::from_rgb(176, 186, 198),
            geometry_edge: Color32::from_rgb(40, 44, 52),
            bodies: [
                Color32::from_rgb(0x8f, 0xb4, 0xd9),
                Color32::from_rgb(0xe0, 0xa4, 0x6c),
                Color32::from_rgb(0x8c, 0xc4, 0x8a),
                Color32::from_rgb(0xd6, 0x8c, 0x96),
                Color32::from_rgb(0xb0, 0x9c, 0xd8),
                Color32::from_rgb(0xd8, 0xc8, 0x70),
                Color32::from_rgb(0x7c, 0xc4, 0xc0),
                Color32::from_rgb(0xc0, 0xa0, 0x84),
            ],
            contour: CONTOUR,
            no_result: Color32::from_rgb(0x9a, 0xa0, 0xaa),
            contour_over: Color32::from_rgb(0x6e, 0x2c, 0x80),
            glyph_constraint: Color32::from_rgb(0x00, 0x8a, 0xa8),
            glyph_load: Color32::from_rgb(0xd0, 0x24, 0x24),
            glyph_interaction: Color32::from_rgb(0x9a, 0x30, 0xc0),
            glyph_contact: Color32::from_rgb(0x6a, 0x8a, 0x00),
            glyph_thermal: Color32::from_rgb(0xe8, 0x4c, 0x00),
            glyph_datum: Color32::from_rgb(0x1e, 0x8c, 0x6e),
            pick_selected: Color32::from_rgb(0xd6, 0x00, 0x7f),
            pick_hover: Color32::from_rgb(0xf0, 0x7c, 0xc0),
        }
        .with_chrome(false)
    }

    /// The dark theme: the kit's dark chrome (the viewport #2e3440, the chrome a little darker,
    /// light text) and FieldCraft's own colours.
    pub const fn dark() -> Tokens {
        let mut t = Tokens::light().with_chrome(true);
        t.heading_bg = Color32::from_rgb(26, 29, 35);
        t.subheading_bg = Color32::from_rgb(31, 35, 42);
        t.ok = Color32::from_rgb(0x8c, 0xc8, 0x8c);
        t.axis_x = Color32::from_rgb(0xe0, 0x5a, 0x50);
        t.axis_y = Color32::from_rgb(0x6c, 0xc0, 0x5a);
        t.axis_z = Color32::from_rgb(0x52, 0x8c, 0xe6);
        t.hex = Color32::from_rgb(0x82, 0xa6, 0xc8);
        t.tet = Color32::from_rgb(0x8c, 0xc0, 0x96);
        t.wedge = Color32::from_rgb(0xc8, 0xb4, 0x82);
        t.pyramid = Color32::from_rgb(0x8c, 0x82, 0xb4);
        t.shell = Color32::from_rgb(0xc8, 0x96, 0x82);
        t.line_element = Color32::from_rgb(0xf0, 0xc0, 0x50);
        t.mesh_edge = Color32::from_rgb(0x14, 0x18, 0x1e);
        t.streamline = Color32::from_rgb(0xd8, 0xde, 0xe8);
        t.node = Color32::from_rgb(0xf2, 0xe6, 0x6e);
        t.geometry = Color32::from_rgb(152, 156, 162);
        t.geometry_edge = Color32::from_rgb(16, 16, 18);
        t.bodies = [
            Color32::from_rgb(0x6c, 0x92, 0xb8),
            Color32::from_rgb(0xbc, 0x84, 0x50),
            Color32::from_rgb(0x6c, 0xa2, 0x6a),
            Color32::from_rgb(0xb4, 0x6e, 0x78),
            Color32::from_rgb(0x90, 0x7e, 0xb6),
            Color32::from_rgb(0xb6, 0xa6, 0x54),
            Color32::from_rgb(0x5c, 0xa2, 0x9e),
            Color32::from_rgb(0xa0, 0x82, 0x66),
        ];
        t.no_result = Color32::from_rgb(0x70, 0x76, 0x80);
        t.glyph_constraint = Color32::from_rgb(0x3c, 0xdc, 0xf0);
        t.glyph_load = Color32::from_rgb(0xff, 0x5a, 0x4c);
        t.glyph_interaction = Color32::from_rgb(0xd2, 0x82, 0xff);
        t.glyph_contact = Color32::from_rgb(0xa0, 0xe6, 0x28);
        t.glyph_thermal = Color32::from_rgb(0xff, 0x8c, 0x00);
        t.glyph_datum = Color32::from_rgb(0x5a, 0xd8, 0xaa);
        t.pick_selected = Color32::from_rgb(0xff, 0x3f, 0xb4);
        t.pick_hover = Color32::from_rgb(0xff, 0xa6, 0xdc);
        t
    }
    /// The contour colour scale for the renderer, low to high.
    pub fn contour_stops(&self) -> [Rgb; 5] {
        self.contour.map(rgb)
    }

    /// The colours of the analysis setup's symbols; their names in the overlay text colour, the
    /// selected definition's in the selection colour.
    pub fn glyph_colors(&self) -> GlyphColors {
        GlyphColors {
            constraint: rgb(self.glyph_constraint),
            load: rgb(self.glyph_load),
            interaction: rgb(self.glyph_interaction),
            contact: rgb(self.glyph_contact),
            highlight: rgb(self.pick_selected),
            thermal: rgb(self.glyph_thermal),
            label: rgb(self.overlay_text),
            datum: rgb(self.glyph_datum),
        }
    }

    /// The viewport's colours for the scene builder ([`fieldcraft_engine::view::build`]).
    pub fn palette(&self) -> Palette {
        Palette {
            background: (rgb(self.viewport_top), rgb(self.viewport_bottom)),
            hex: rgb(self.hex),
            tet: rgb(self.tet),
            wedge: rgb(self.wedge),
            pyramid: rgb(self.pyramid),
            shell: rgb(self.shell),
            line_element: rgb(self.line_element),
            mesh_edge: rgb(self.mesh_edge),
            streamline: rgb(self.streamline),
            node: rgb(self.node),
            geometry: rgb(self.geometry),
            geometry_edge: rgb(self.geometry_edge),
            bodies: self.bodies.map(rgb),
            no_result: rgb(self.no_result),
            contour: self.contour_stops(),
            contour_over: rgb(self.contour_over),
            marked: rgb(self.error),
            merged: rgb(self.warm),
            glyphs: self.glyph_colors(),
        }
    }

    /// The colour elements of this kind are drawn in.
    pub fn element(&self, kind: ElementKind) -> Color32 {
        match kind {
            ElementKind::Hex8 | ElementKind::Hex20 => self.hex,
            ElementKind::Tet4 | ElementKind::Tet10 => self.tet,
            ElementKind::Wedge6 => self.wedge,
            ElementKind::Pyramid5 => self.pyramid,
            ElementKind::Tri3 | ElementKind::Tri6 | ElementKind::Quad4 | ElementKind::Quad8 => self.shell,
            ElementKind::Rod2 | ElementKind::Rod3 | ElementKind::Beam2 => self.line_element,
        }
    }
}

/// A colour for the renderer.
pub fn rgb(c: Color32) -> Rgb {
    Rgb(c.r(), c.g(), c.b())
}

/// The label font at `size`: the kit's semibold once the context has it (fonts set with [`apply`]
/// take effect from the next pass), the default proportional font until then.
pub fn label_font(ctx: &egui::Context, size: f32) -> egui::FontId {
    let font = type_face::semibold(size);
    if ctx.fonts(|f| f.families().contains(&font.family)) { font } else { egui::FontId::proportional(size) }
}

/// The icons' base accent ramp (`active-1` palest to `active-3` strongest): SolidCraft's blue
/// re-hued to FieldCraft's teal at the same OKLab lightness and chroma, so icons carry the app's
/// accent as SolidCraft's carry its own.
const ICON_RAMP_LIGHT: [[u8; 3]; 3] = [[0xb0, 0xe4, 0xdd], [0x06, 0x9a, 0x90], [0x05, 0x7b, 0x73]];
const ICON_RAMP_DARK: [[u8; 3]; 3] = [[0x02, 0xb0, 0xa5], [0x0c, 0x96, 0x8d], [0x09, 0x84, 0x7b]];

/// The icon palette for a theme: SolidCraft's, with the base ramp in teal.
fn icon_palette(dark: bool) -> kit_icons::Palette {
    let mut palette = kit_icons::Palette::default();
    palette.accents[0] = Some(if dark { ICON_RAMP_DARK } else { ICON_RAMP_LIGHT });
    palette
}

/// FieldCraft's English text for the kit's message ids (its tooltips, history controls, view
/// compass): the kit reads every one through `solidcraft_ui_text`, and an id with no catalogue
/// reads as itself.
/// The catalogue's text: `locales/en.json`, with the workspaces' lines and tab lists from the
/// capability catalogue (`crate::toolbar::workspace_ids`).
static ENGLISH_SOURCE: std::sync::LazyLock<String> = std::sync::LazyLock::new(|| {
    let mut messages: serde_json::Map<String, serde_json::Value> = serde_json::from_str(include_str!("../locales/en.json")).unwrap_or_default();
    for w in fieldcraft_engine::catalogue::Workspace::ALL {
        let (about, tabs) = crate::toolbar::workspace_ids(w);
        messages.insert(about.into(), w.about().into());
        messages.insert(tabs.into(), w.tabs().iter().map(|t| t.label()).collect::<Vec<_>>().join(" · ").into());
    }
    serde_json::Value::Object(messages).to_string()
});
static ENGLISH: std::sync::LazyLock<solidcraft_ui_text::Catalog> =
    std::sync::LazyLock::new(|| solidcraft_ui_text::Catalog::new("en", "English", ENGLISH_SOURCE.as_str()));
static MESSAGES: std::sync::LazyLock<[&'static solidcraft_ui_text::Catalog; 1]> = std::sync::LazyLock::new(|| [&*ENGLISH]);

/// SolidCraft's look installed into the kit, once per process: FieldCraft's icon set (SolidCraft's
/// with ours, [`crate::icons::SET`]), SolidCraft's code-drawn glyphs and its built-in themes, and
/// FieldCraft's text for the kit's messages.
pub fn install_look() {
    solidcraft_ui_text::install(&*MESSAGES);
    crate::icons::install();
    solidcraft_look::glyphs::install();
    solidcraft_look::themes::install();
}

/// Style the context for the dark or light theme: SolidCraft's look installed in the kit, with
/// FieldCraft's tokens (its teal accent on SolidCraft's chrome) and teal icon accents; the kit's
/// visuals from them, and SolidCraft's interface faces (IBM Plex). egui draws selected labels'
/// text and links in the accent's ink, which reads on the selection fill and on panels.
pub fn apply(ctx: &egui::Context, dark: bool) {
    install_look();
    kit_theme::set_dark(dark);
    kit_theme::set_custom(Some(kit_theme::Custom { dark, tokens: kit(dark) }));
    kit_icons::set_theme_palette(Some((dark, icon_palette(dark))));
    solidcraft_look::fonts::install(ctx);
    kit_theme::apply(ctx);
    let (k, ink) = (kit(dark), Tokens::get(dark).accent_ink);
    ctx.global_style_mut(|s| {
        s.visuals.selection.stroke.color = ink;
        s.visuals.hyperlink_color = ink;
        // A plain egui field outside the kit's forms keeps an outline at rest, as the kit's do.
        s.visuals.widgets.inactive.bg_stroke = k.control_border_stroke(k.control_border);
    });
}

#[cfg(test)]
#[path = "theme_colour_tests.rs"]
mod colour_tests;

#[cfg(test)]
mod tests {
    use super::*;

    /// The visuals [`apply`] gives a context in a theme.
    pub(super) fn visuals(dark: bool) -> egui::Visuals {
        let ctx = egui::Context::default();
        apply(&ctx, dark);
        ctx.global_style().visuals.clone()
    }

    /// Empty, unfocused fields still draw an outline in either theme.
    #[test]
    fn empty_text_fields_have_a_visible_outline() {
        for dark in [false, true] {
            let ctx = egui::Context::default();
            apply(&ctx, dark);
            let mut field_rect = egui::Rect::NOTHING;
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                field_rect = ui.text_edit_singleline(&mut String::new()).rect;
            });
            out.textures_delta.clear();
            assert!(
                out.shapes.iter().any(|shape| {
                    matches!(&shape.shape, egui::Shape::Rect(rect)
                    if rect.rect == field_rect
                        && rect.stroke.width >= 1.0
                        && rect.stroke.color != Color32::TRANSPARENT
                        && rect.stroke.color != Tokens::get(dark).panel)
                }),
                "empty field has no visible outline: dark={dark}"
            );
        }
    }

    /// Title and section bands stand out from the form behind them, the title more than its
    /// sections.
    #[test]
    fn heading_bands_stand_out_in_order() {
        let luma = |c: Color32| i32::from(c.r()) * 3 + i32::from(c.g()) * 6 + i32::from(c.b());
        for t in [Tokens::light(), Tokens::dark()] {
            let off = |c: Color32| (luma(c) - luma(t.dialog_bg)).abs();
            assert!(off(t.heading_bg) > off(t.subheading_bg), "{t:?}");
            assert!(off(t.subheading_bg) > 60, "{t:?}");
        }
    }

    /// The kit's messages read as FieldCraft's text, not as their ids (the browser's eye and fold
    /// tooltips, the history controls).
    #[test]
    fn the_kits_messages_read_as_text() {
        install_look();
        for (id, text) in [("ui.hide", "Hide"), ("ui.expand", "Expand"), ("ui.history_back", "One step back (undo)")] {
            assert_eq!(solidcraft_ui_text::tr(id), text);
        }
        assert_eq!(solidcraft_ui_text::format("ui.viewport_status_clear", &[("name", "Mesh")]), "Clear Mesh");
    }

    /// The label font is used only once the context has it; before that, the default font.
    #[test]
    fn the_label_font_falls_back_until_it_is_installed() {
        let ctx = egui::Context::default();
        let run = |ctx: &egui::Context| {
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                let f = label_font(ui.ctx(), 12.0);
                // Laying text out in it never panics.
                let _ = ui.painter().layout_no_wrap("PREPARE".into(), f, Color32::WHITE);
            });
            out.textures_delta.clear();
        };
        run(&ctx);
        assert_eq!(label_font(&ctx, 12.0).family, egui::FontFamily::Proportional);
        apply(&ctx, true);
        run(&ctx);
        assert_eq!(label_font(&ctx, 12.0), type_face::semibold(12.0));
        run(&ctx);
    }

    /// The kit draws with FieldCraft's colours once they are applied: its tokens are [`kit`]'s,
    /// in the theme applied, and the chrome FieldCraft draws itself is the same.
    #[test]
    fn the_kit_draws_with_fieldcrafts_tokens() {
        let ctx = egui::Context::default();
        for dark in [false, true] {
            apply(&ctx, dark);
            let (k, t) = (kit_theme::Tokens::get(), Tokens::get(dark));
            assert_eq!(k, kit(dark), "dark={dark}");
            assert!(kit_theme::is_dark() == dark);
            for (a, b) in [(k.accent, t.accent), (k.panel, t.panel), (k.text, t.text), (k.dialog_body, t.dialog_bg), (k.selected_ink, t.accent_ink)] {
                assert_eq!(a, b, "dark={dark}");
            }
            assert_eq!(ctx.global_style().visuals.panel_fill, t.panel, "dark={dark}");
        }
    }

    /// Setup symbols differ from every element colour and from each other, and stand out from
    /// the viewport background.
    #[test]
    fn symbol_colours_stand_out() {
        let luma = |c: Color32| i32::from(c.r()) * 3 + i32::from(c.g()) * 6 + i32::from(c.b());
        let far = |a: Color32, b: Color32| {
            (i32::from(a.r()) - i32::from(b.r())).abs() + (i32::from(a.g()) - i32::from(b.g())).abs() + (i32::from(a.b()) - i32::from(b.b())).abs()
        };
        for dark in [true, false] {
            let t = Tokens::get(dark);
            let symbols = [t.glyph_constraint, t.glyph_load, t.glyph_interaction, t.glyph_contact, t.glyph_thermal];
            for (i, s) in symbols.iter().enumerate() {
                assert!((luma(*s) - luma(t.viewport_bottom)).abs() > 400 || far(*s, t.viewport_bottom) > 250, "{s:?} dark={dark}");
                for k in ElementKind::ALL {
                    assert!(far(*s, t.element(k)) > 90, "{s:?} vs {k:?} dark={dark}");
                }
                for o in symbols.iter().skip(i + 1) {
                    assert!(far(*s, *o) > 90, "{s:?} vs {o:?} dark={dark}");
                }
            }
        }
    }

    /// Streamlines in one colour (the legend showing pressure) stand out from the viewport
    /// background in both themes (#193 review: they were near-black on the dark one).
    #[test]
    fn streamlines_stand_out_from_the_background() {
        let luma = |c: Color32| i32::from(c.r()) * 3 + i32::from(c.g()) * 6 + i32::from(c.b());
        for dark in [true, false] {
            let t = Tokens::get(dark);
            for bg in [t.viewport_top, t.viewport_bottom] {
                assert!((luma(t.streamline) - luma(bg)).abs() > 1200, "dark={dark}");
            }
            assert_eq!(t.palette().streamline, rgb(t.streamline));
        }
    }

    /// The selection and hover colours differ from each other, from every element, geometry,
    /// node and symbol colour, and from the background.
    #[test]
    fn selection_colours_stand_out() {
        let far = |a: Color32, b: Color32| {
            (i32::from(a.r()) - i32::from(b.r())).abs() + (i32::from(a.g()) - i32::from(b.g())).abs() + (i32::from(a.b()) - i32::from(b.b())).abs()
        };
        for dark in [true, false] {
            let t = Tokens::get(dark);
            assert!(far(t.pick_selected, t.pick_hover) > 90, "dark={dark}");
            for c in [t.pick_selected, t.pick_hover] {
                let others = [
                    t.geometry,
                    t.node,
                    t.glyph_constraint,
                    t.glyph_load,
                    t.glyph_interaction,
                    t.glyph_contact,
                    t.glyph_thermal,
                    t.viewport_bottom,
                    t.mesh_edge,
                ];
                for o in others.into_iter().chain(ElementKind::ALL.map(|k| t.element(k))) {
                    assert!(far(c, o) > 100, "{c:?} vs {o:?} dark={dark}");
                }
            }
        }
    }

    /// Every element kind has a colour, and edges stand out from the faces they lie on.
    #[test]
    fn element_colours_contrast_with_edges() {
        let luma = |c: Color32| u32::from(c.r()) * 3 + u32::from(c.g()) * 6 + u32::from(c.b());
        for dark in [true, false] {
            let t = Tokens::get(dark);
            for k in ElementKind::ALL {
                assert!(luma(t.element(k)) > luma(t.mesh_edge) + 600, "{k:?} dark={dark}");
            }
        }
    }
}
