//! Colours and widget style: the Craft family's dark theme (the default, blue-grey chrome around a
//! #2e3440 viewport) and its light one, with the same chrome as SolidCraft's tokens so the apps
//! look alike side by side, and FieldCraft's own violet accent so they are told apart. Every
//! colour the UI uses comes from [`Tokens`].
//!
//! Adapted from SolidCraft's `crates/ui-egui/src/theme.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.

use egui::{Color32, Stroke};
use fieldcraft_engine::view::Palette;
use fieldcraft_model::ElementKind;
use fieldcraft_render::Rgb;
use fieldcraft_render::glyphs::GlyphColors;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tokens {
    /// The title bar (dark in both themes) and its text.
    pub app_bar: Color32,
    pub app_bar_text: Color32,
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
    /// FieldCraft's accent (violet): the primary action, the current tab, focus and links. Violet
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
    /// Icon line work and fills.
    pub icon: Color32,
    pub icon_fill: Color32,
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

/// The accent in the light theme: #6F47C2.
const LIGHT_ACCENT: Color32 = Color32::from_rgb(0x6f, 0x47, 0xc2);
/// The accent in the dark theme: #8861CF, the chosen #8A63D2 darkened by 0.008 in OKLCH lightness
/// (hue and chroma kept) so white on it reaches 4.5:1 (#8A63D2 gives 4.38:1).
const DARK_ACCENT: Color32 = Color32::from_rgb(0x88, 0x61, 0xcf);
/// How far a hovered accent fill moves toward black: darker keeps white on it legible.
const ACCENT_HOVER: f32 = 0.12;
const LIGHT_PANEL: Color32 = Color32::from_rgb(250, 251, 252);
const DARK_PANEL: Color32 = Color32::from_rgb(39, 43, 52);
const DARK_PANEL_HEADER: Color32 = Color32::from_rgb(47, 52, 62);

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

impl Tokens {
    /// Both themes, evaluated at compile time.
    const LIGHT: Tokens = Tokens::light();
    const DARK: Tokens = Tokens::dark();

    pub fn get(dark: bool) -> Tokens {
        if dark { Tokens::DARK } else { Tokens::LIGHT }
    }

    pub const fn light() -> Tokens {
        Tokens {
            app_bar: Color32::from_rgb(52, 56, 64),
            app_bar_text: Color32::from_rgb(230, 234, 242),
            toolbar: Color32::from_rgb(246, 247, 249),
            panel: LIGHT_PANEL,
            panel_header: Color32::from_rgb(236, 238, 242),
            dialog_bg: Color32::from_rgb(250, 251, 252),
            heading_bg: Color32::from_rgb(218, 223, 231),
            subheading_bg: Color32::from_rgb(234, 237, 242),
            border: Color32::from_rgb(208, 212, 220),
            field: Color32::from_rgb(232, 236, 242),
            text: Color32::from_rgb(38, 42, 50),
            text_dim: Color32::from_rgb(110, 116, 128),
            accent: LIGHT_ACCENT,
            accent_ink: LIGHT_ACCENT,
            on_accent: Color32::WHITE,
            accent_hover: mix(LIGHT_ACCENT, Color32::BLACK, ACCENT_HOVER),
            accent_soft: mix(LIGHT_PANEL, LIGHT_ACCENT, 0.18),
            hover: mix(LIGHT_PANEL, LIGHT_ACCENT, 0.09),
            error: Color32::from_rgb(214, 54, 54),
            ok: Color32::from_rgb(44, 122, 58),
            warm: Color32::from_rgb(222, 150, 30),
            icon: Color32::from_rgb(70, 78, 92),
            icon_fill: Color32::from_rgb(178, 196, 222),
            timeline: Color32::from_rgb(232, 235, 240),
            timeline_item: Color32::WHITE,
            timeline_marker: Color32::from_rgb(70, 76, 88),
            viewport_top: Color32::from_rgb(240, 243, 247),
            viewport_bottom: Color32::from_rgb(203, 210, 221),
            overlay: Color32::from_rgba_unmultiplied_const(255, 255, 255, 235),
            overlay_text: Color32::from_rgb(38, 42, 50),
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
    }

    /// Dark blue-grey: the viewport #2e3440, the chrome a little darker, light text, the violet
    /// accent.
    pub const fn dark() -> Tokens {
        let mut t = Tokens::light();
        t.app_bar = Color32::from_rgb(28, 31, 38);
        t.toolbar = Color32::from_rgb(37, 41, 49);
        t.panel = DARK_PANEL;
        t.panel_header = DARK_PANEL_HEADER;
        t.dialog_bg = Color32::from_rgb(59, 68, 83);
        t.heading_bg = Color32::from_rgb(38, 44, 55);
        t.subheading_bg = Color32::from_rgb(48, 55, 68);
        t.border = Color32::from_rgb(63, 69, 82);
        t.field = Color32::from_rgb(30, 34, 41);
        t.text = Color32::from_rgb(222, 226, 232);
        t.text_dim = Color32::from_rgb(148, 155, 168);
        t.accent = DARK_ACCENT;
        // The accent lightened in OKLCH (hue kept) until it reads at 4.5:1 on every dark ground
        // text sits on, the form background (#3b4453) being the closest.
        t.accent_ink = Color32::from_rgb(0xc3, 0x9e, 0xff);
        t.accent_hover = mix(DARK_ACCENT, Color32::BLACK, ACCENT_HOVER);
        t.accent_soft = mix(DARK_PANEL, DARK_ACCENT, 0.30);
        t.hover = mix(DARK_PANEL_HEADER, DARK_ACCENT, 0.10);
        t.error = Color32::from_rgb(240, 100, 96);
        t.ok = Color32::from_rgb(0x8c, 0xc8, 0x8c);
        t.warm = Color32::from_rgb(236, 170, 70);
        t.icon = Color32::from_rgb(208, 214, 226);
        t.icon_fill = Color32::from_rgb(86, 112, 150);
        t.timeline = Color32::from_rgb(33, 37, 45);
        t.timeline_item = Color32::from_rgb(52, 58, 70);
        t.timeline_marker = Color32::from_rgb(170, 178, 192);
        t.viewport_top = Color32::from_rgb(46, 52, 64);
        t.viewport_bottom = Color32::from_rgb(46, 52, 64);
        t.overlay = Color32::from_rgba_unmultiplied_const(34, 38, 46, 240);
        t.overlay_text = Color32::from_rgb(222, 226, 232);
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

fn visuals(dark: bool) -> egui::Visuals {
    let t = Tokens::get(dark);
    let mut v = if dark { egui::Visuals::dark() } else { egui::Visuals::light() };
    v.panel_fill = t.panel;
    v.window_fill = t.panel;
    v.window_stroke = Stroke::new(1.0, t.border);
    v.extreme_bg_color = t.field;
    v.faint_bg_color = t.panel_header;
    v.override_text_color = Some(t.text);
    let w = &mut v.widgets;
    let noninteractive = w.noninteractive.bg_fill;
    for (state, fill) in [
        (&mut w.noninteractive, noninteractive),
        (&mut w.inactive, t.panel_header),
        (&mut w.hovered, t.hover),
        (&mut w.active, t.accent_soft),
        (&mut w.open, t.panel_header),
    ] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.corner_radius = egui::CornerRadius::same(3);
        state.fg_stroke = Stroke::new(1.0, t.text);
    }
    v.selection.bg_fill = t.accent_soft;
    // egui draws selected labels' text in the selection stroke's colour, on its fill.
    v.selection.stroke = Stroke::new(1.0, t.accent_ink);
    v.hyperlink_color = t.accent_ink;
    v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, t.border);
    v.widgets.inactive.bg_stroke = Stroke::new(1.0, t.border);
    v.window_corner_radius = egui::CornerRadius::same(6);
    v.menu_corner_radius = egui::CornerRadius::same(6);
    v
}

/// The font family of the interface's caps labels (tabs, panel labels, workspaces, headings):
/// Inter SemiBold, then egui's own proportional fonts for what it lacks. egui has no bold weight
/// of its default font, so labels need a face of their own to read at small sizes.
pub const LABEL_FAMILY: &str = "label";
/// Inter SemiBold 4.1 (SIL OFL 1.1, unmodified; fonts/OFL.txt, ATTRIBUTION.md).
static INTER_SEMIBOLD: &[u8] = include_bytes!("../fonts/Inter-SemiBold.ttf");

/// The font definitions: egui's defaults plus the label family.
pub fn fonts() -> egui::FontDefinitions {
    let mut f = egui::FontDefinitions::default();
    f.font_data.insert("inter-semibold".into(), std::sync::Arc::new(egui::FontData::from_static(INTER_SEMIBOLD)));
    let mut family = vec!["inter-semibold".to_string()];
    family.extend(f.families.get(&egui::FontFamily::Proportional).cloned().unwrap_or_default());
    f.families.insert(egui::FontFamily::Name(LABEL_FAMILY.into()), family);
    f
}

/// The caps-label font at `size`: the label family once the context has it (fonts set with
/// [`apply`] take effect from the next pass), the default proportional font until then.
pub fn label_font(ctx: &egui::Context, size: f32) -> egui::FontId {
    let family = egui::FontFamily::Name(LABEL_FAMILY.into());
    if ctx.fonts(|f| f.families().contains(&family)) { egui::FontId::new(size, family) } else { egui::FontId::proportional(size) }
}

/// Style the context for the dark or light theme, with the interface's fonts.
pub fn apply(ctx: &egui::Context, dark: bool) {
    ctx.set_fonts(fonts());
    ctx.set_theme(if dark { egui::Theme::Dark } else { egui::Theme::Light });
    ctx.set_visuals(visuals(dark));
    ctx.global_style_mut(|s| {
        s.spacing.item_spacing = egui::vec2(6.0, 4.0);
        s.spacing.button_padding = egui::vec2(6.0, 3.0);
    });
}

#[cfg(test)]
#[path = "theme_colour_tests.rs"]
mod colour_tests;

#[cfg(test)]
mod tests {
    use super::*;

    /// Empty, unfocused fields still draw an outline in either theme.
    #[test]
    fn empty_text_fields_have_a_visible_outline() {
        for dark in [false, true] {
            let ctx = egui::Context::default();
            ctx.set_visuals(visuals(dark));
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

    /// Title and section bands stand out from a panel, the title more than its sections.
    #[test]
    fn heading_bands_are_darker_than_panels_in_order() {
        let luma = |c: Color32| u32::from(c.r()) * 3 + u32::from(c.g()) * 6 + u32::from(c.b());
        for t in [Tokens::light(), Tokens::dark()] {
            assert!(luma(t.heading_bg) < luma(t.subheading_bg), "{t:?}");
            assert!(luma(t.subheading_bg) < luma(t.dialog_bg), "{t:?}");
        }
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
        assert_eq!(label_font(&ctx, 12.0).family, egui::FontFamily::Name(LABEL_FAMILY.into()));
        run(&ctx);
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
