//! The navigation bar, and the mouse navigation presets: which mouse button, with which keys
//! held, orbits, pans or zooms.
//!
//! The bar holds SolidCraft's default slots where FieldCraft has them (SolidCraft's
//! `crates/ui-egui/src/navigation_slots.rs` `default_bar` and `navigation_bar_view.rs`; MIT OR
//! Apache-2.0, see ATTRIBUTION.md): orbit, pan and zoom; fit and home; the visual style and the
//! display menu, which ends the bar. The kit paints the bar and its slots
//! (`solidcraft_ui_kit::nav_bar`); the display menu is FieldCraft's own (its switches and slider
//! are more than the kit's slot menus hold).
//!
//! The presets after FieldCraft's own follow SolidCraft's `crates/ui-egui/src/prefs.rs`
//! (`NavPreset`, `NavBinding`, `NAV_PRESETS`; MIT OR Apache-2.0), see ATTRIBUTION.md. In every
//! preset a left drag the preset does not bind does what the navigation bar's mode says (orbit
//! by default), and a Ctrl (Cmd) + left drag that starts on empty space orbits, unless the
//! preset binds Ctrl + left itself (SolidCraft #471).

use egui::{Rect, pos2, vec2};
use solidcraft_ui_kit::{self as kit, BackingTooltipExt as _};

use crate::{ModelApp, NavMode, Style};

/// The bar's slots: (key, kind, name, tooltip, shortcut), in SolidCraft's default bar's order
/// where FieldCraft has the slot (it has no look-at-a-face, between fit and home).
const BAR: [(&str, kit::NavigationKind, &str, &str, Option<&str>); 9] = [
    ("orbit", kit::NavigationKind::Tool, "Orbit", "Orbit: left drag turns the view", None),
    ("pan", kit::NavigationKind::Tool, "Pan", "Pan: left drag moves the view (or right or middle drag)", None),
    ("zoom", kit::NavigationKind::Tool, "Zoom", "Zoom: left drag up and down (or the wheel)", None),
    ("separator", kit::NavigationKind::Separator, "", "", None),
    ("fit", kit::NavigationKind::Action, "Fit", "Fit the model (or a double click)", Some("F")),
    ("home", kit::NavigationKind::Action, "Home view", "Home view", None),
    ("separator", kit::NavigationKind::Separator, "", "", None),
    ("visual_style", kit::NavigationKind::Menu, "Visual style", "Visual style", None),
    ("display", kit::NavigationKind::Menu, "Display", "Display: what the view shows, the projection, the panels and the mouse", None),
];

/// A slot's icon (SolidCraft's artwork): the visual style's is the current style's, as
/// SolidCraft's `navigation_slots::visual_style_icon`.
fn icon(key: &str, style: Style) -> &'static str {
    match key {
        "orbit" => "ui.orbit",
        "pan" => "ui.pan",
        "zoom" => "ui.zoom",
        "fit" => "ui.fit",
        "home" => "ui.home",
        "visual_style" => style_icon(style),
        "display" => "display",
        _ => "",
    }
}

/// Each visual style's icon: SolidCraft's for the same style.
fn style_icon(style: Style) -> &'static str {
    match style {
        Style::ShadedEdges => "ui.style",
        Style::Shaded => "style_shaded",
        Style::Wireframe => "style_wireframe",
    }
}

/// Where a slot is published in the window's layout (`ui.inspect`'s `layout`).
fn region(key: &str) -> Option<&'static str> {
    Some(match key {
        "orbit" => "nav_bar:orbit",
        "pan" => "nav_bar:pan",
        "zoom" => "nav_bar:zoom",
        "fit" => "nav_bar:fit",
        "home" => "nav_bar:home",
        "visual_style" => "nav_bar:visual_style",
        "display" => "nav_bar:display",
        _ => return None,
    })
}

/// The menu key ([`crate::UiState::menu`]) a menu slot opens.
fn menu_key(key: &str) -> Option<&'static str> {
    match key {
        "visual_style" => Some("visual_style"),
        "display" => Some("view"),
        _ => None,
    }
}

/// How wide a slot of `kind` is on the bar.
fn slot_width(kind: kit::NavigationKind) -> f32 {
    let m = kit::NavigationMetrics::for_size(kit::NavigationSize::Regular);
    match kind {
        kit::NavigationKind::Separator => kit::NAV_SEPARATOR,
        kit::NavigationKind::Menu => m.icon + kit::NAV_MENU_LEFT + kit::NAV_MENU_RIGHT + kit::NAV_CHEVRON + kit::NAV_MENU_GAP,
        _ => m.button,
    }
}

/// The bar's rect, bottom centre of the view `rect`.
pub(crate) fn bar_rect(rect: Rect) -> Rect {
    let m = kit::NavigationMetrics::for_size(kit::NavigationSize::Regular);
    let inner = BAR.iter().map(|s| slot_width(s.1)).sum::<f32>() + kit::NAV_GAP * (BAR.len() - 1) as f32;
    let width = inner + 2.0 * (kit::NAV_BORDER + kit::NAV_PADDING);
    Rect::from_min_size(pos2(rect.center().x - width / 2.0, rect.bottom() - m.height() - kit::NAV_BOTTOM), vec2(width, m.height()))
}

/// The navigation bar at the bottom centre of the view `rect`, and its two menus.
pub(crate) fn bar(app: &mut ModelApp, ui: &mut egui::Ui, rect: Rect) {
    let bar = bar_rect(rect);
    if !rect.contains_rect(bar) {
        return;
    }
    app.place("nav_bar", bar);
    let m = kit::NavigationMetrics::for_size(kit::NavigationSize::Regular);
    let mut x = bar.left() + kit::NAV_BORDER + kit::NAV_PADDING;
    let slots: Vec<kit::nav_bar::SlotView> = BAR
        .iter()
        .enumerate()
        .map(|(index, &(key, kind, name, tip, shortcut))| {
            let rect = Rect::from_min_size(pos2(x, bar.top() + kit::NAV_BORDER + kit::NAV_PADDING), vec2(slot_width(kind), m.button));
            x = rect.right() + kit::NAV_GAP;
            let selected = match key {
                "orbit" => app.nav == NavMode::Orbit,
                "pan" => app.nav == NavMode::Pan,
                "zoom" => app.nav == NavMode::Zoom,
                _ => menu_key(key).is_some_and(|k| app.ui.menu.as_deref() == Some(k)),
            };
            let tip = if key == "visual_style" { format!("{tip}: {}", app.ui.style.label()) } else { tip.into() };
            kit::nav_bar::SlotView {
                index,
                id: key,
                kind,
                icon: icon(key, app.ui.style),
                text: name.into(),
                name: name.into(),
                tip,
                key: shortcut.map(str::to_string),
                selected,
                enabled: true,
                hidden: false,
                show_text: false,
                rect,
                alias: None,
            }
        })
        .collect();
    let responses = kit::nav_bar::paint_slots(ui, bar, kit::NavigationSize::Regular, kit::NavigationEdge::Bottom, &slots, None);
    for (slot, resp) in slots.iter().zip(&responses) {
        if let Some(region) = region(slot.id) {
            app.place(region, slot.rect);
        }
        if resp.clicked() {
            match slot.id {
                "orbit" => app.nav = NavMode::Orbit,
                "pan" => app.nav = NavMode::Pan,
                "zoom" => app.nav = NavMode::Zoom,
                "fit" => app.fit_view(),
                "home" => app.home_view(),
                key => {
                    if let Some(k) = menu_key(key) {
                        app.toggle_menu(k);
                    }
                }
            }
        }
    }
    let at = |key: &str| slots.iter().find(|s| s.id == key).map(|s| s.rect);
    if let Some(anchor) = at("visual_style") {
        crate::frame::settings_menu(app, ui, "visual_style", anchor, style_menu);
    }
    if let Some(anchor) = at("display") {
        crate::frame::settings_menu(app, ui, "view", anchor, display_menu);
    }
}

/// The visual style menu: the styles, the current one checked (SolidCraft's `style_menu`).
fn style_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    for s in Style::ALL {
        let mut row = kit::MenuItem::new(s.label(), s.label(), style_icon(s));
        row.checked = Some(app.ui.style == s);
        let r = row.show(ui);
        if s == Style::ALL[0] {
            app.place("view_menu:style", r.rect);
        }
        if r.clicked() {
            app.ui.style = s;
            app.ui.menu = None;
        }
    }
}

/// The display menu (SolidCraft's View menu at the bar's end): what the view shows and the
/// symbols' size, the projection, the panels, and the mouse preset.
fn display_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    crate::frame::menu_heading(ui, "Show");
    // Checked menu rows: one text inset for the whole menu.
    let mut first = None;
    for (label, on) in [
        ("Geometry", &mut app.ui.show_geometry),
        ("Mesh", &mut app.ui.show_mesh),
        ("Nodes", &mut app.ui.show_nodes),
        ("Loads & BCs", &mut app.ui.show_loads),
        ("Interactions", &mut app.ui.show_interactions),
    ] {
        let r = toggle_row(ui, label, on);
        first.get_or_insert(r.rect);
    }
    if let Some(r) = first {
        app.place("view_menu:show", r);
    }
    kit::scope(ui, |ui| {
        let (lo, hi) = fieldcraft_render::glyphs::SIZE_RANGE;
        kit::row(ui, "Symbol size", |ui| {
            kit::slider(ui, &mut app.ui.symbol_size, lo..=hi, kit::Slider { step: Some(1.0), format: |v| format!("{v:.0} pt"), ..Default::default() })
        })
        .backing_on_hover_text("The size of boundary condition, load and interaction symbols on screen, at any zoom");
    });
    crate::frame::menu_separator(ui);
    crate::frame::menu_heading(ui, "Projection");
    for (perspective, label) in [(false, "Orthographic"), (true, "Perspective")] {
        let mut row = kit::MenuItem::new(label, label, "ui.perspective");
        row.checked = Some(app.ui.perspective == perspective);
        let r = row.show(ui);
        if perspective {
            app.place("view_menu:perspective", r.rect);
        }
        if r.clicked() {
            app.ui.perspective = perspective;
        }
    }
    crate::frame::menu_separator(ui);
    toggle_row(ui, "Browser", &mut app.ui.show_tree);
    toggle_row(ui, "Timeline", &mut app.ui.show_timeline);
    crate::frame::menu_separator(ui);
    crate::frame::menu_heading(ui, "Navigation");
    kit::scope(ui, |ui| {
        let current = preset(&app.ui.nav_preset);
        let labels: Vec<&str> = NAV_PRESETS.iter().map(|p| p.label).collect();
        let mut selected = NAV_PRESETS.iter().position(|p| p.id == current.id).unwrap_or(0);
        kit::choice_row(ui, "Preset", "fc_nav_preset", &mut selected, &labels).backing_on_hover_text(current.help);
        if let Some(p) = NAV_PRESETS.get(selected).filter(|p| p.id != current.id) {
            app.ui.nav_preset = p.id.into();
            app.ui.reverse_zoom = p.reverse_zoom;
        }
        kit::checkbox(ui, "Reverse zoom", &mut app.ui.reverse_zoom).backing_on_hover_text("The wheel zooms the other way");
    });
}

/// A menu row that shows or hides something: a checked menu row, flipped by a click.
fn toggle_row(ui: &mut egui::Ui, label: &str, on: &mut bool) -> egui::Response {
    let mut row = kit::MenuItem::new(label, label, "");
    row.checked = Some(*on);
    let r = row.show(ui);
    if r.clicked() {
        *on = !*on;
    }
    r
}

/// A mouse button that can drag the camera.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NavButton {
    Left,
    Middle,
    Right,
}

/// Modifier keys a binding needs held (others may be held too).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct NavKeys {
    pub shift: bool,
    pub ctrl: bool,
    pub alt: bool,
}

impl NavKeys {
    const NONE: NavKeys = NavKeys { shift: false, ctrl: false, alt: false };
    const SHIFT: NavKeys = NavKeys { shift: true, ctrl: false, alt: false };
    const CTRL: NavKeys = NavKeys { shift: false, ctrl: true, alt: false };
    const SHIFT_ALT: NavKeys = NavKeys { shift: true, ctrl: false, alt: true };

    /// All of these keys are held (Ctrl is Cmd on macOS).
    fn held(self, m: egui::Modifiers) -> bool {
        (!self.shift || m.shift) && (!self.ctrl || m.ctrl || m.command) && (!self.alt || m.alt)
    }
}

/// One camera binding: dragging `button` with `keys` held does `mode`.
#[derive(Clone, Copy, Debug)]
pub struct NavBinding {
    pub button: NavButton,
    pub keys: NavKeys,
    pub mode: NavMode,
}

const fn bind(button: NavButton, keys: NavKeys, mode: NavMode) -> NavBinding {
    NavBinding { button, keys, mode }
}

/// A mouse mapping: its bindings, most specific first (the first one whose button is dragged
/// with all its keys held wins).
pub struct NavPreset {
    pub id: &'static str,
    pub label: &'static str,
    /// What the mouse does, for the settings menu.
    pub help: &'static str,
    /// Choosing the preset turns the wheel the other way.
    pub reverse_zoom: bool,
    pub bindings: &'static [NavBinding],
}

impl NavPreset {
    /// What dragging `button` with `mods` held does, if the preset binds it.
    pub fn mode(&self, button: NavButton, mods: egui::Modifiers) -> Option<NavMode> {
        self.bindings.iter().find(|b| b.button == button && b.keys.held(mods)).map(|b| b.mode)
    }
}

use NavButton::{Left, Middle, Right};
use NavMode::{Orbit, Pan, Zoom};

/// The id of the default preset.
pub const DEFAULT_PRESET: &str = "fieldcraft";

/// The mouse mappings. The first is the default.
pub const NAV_PRESETS: [NavPreset; 7] = [
    NavPreset {
        id: DEFAULT_PRESET,
        label: "FieldCraft",
        help: "Left drag does the navigation bar's mode (orbit, pan or zoom); right, middle or Shift+left drag pans",
        reverse_zoom: false,
        bindings: &[bind(Left, NavKeys::SHIFT, Pan), bind(Right, NavKeys::NONE, Pan), bind(Middle, NavKeys::NONE, Pan)],
    },
    NavPreset {
        id: "fusion",
        label: "Fusion",
        help: "Middle drag pans, Shift+middle orbits, right drag orbits",
        reverse_zoom: false,
        bindings: &[bind(Middle, NavKeys::SHIFT, Orbit), bind(Middle, NavKeys::NONE, Pan), bind(Right, NavKeys::NONE, Orbit)],
    },
    NavPreset {
        id: "solidworks",
        label: "SolidWorks",
        help: "Middle drag orbits, Ctrl+middle pans, Shift+middle zooms, right drag orbits",
        reverse_zoom: true,
        bindings: &[
            bind(Middle, NavKeys::CTRL, Pan),
            bind(Middle, NavKeys::SHIFT, Zoom),
            bind(Middle, NavKeys::NONE, Orbit),
            bind(Right, NavKeys::NONE, Orbit),
        ],
    },
    NavPreset {
        id: "inventor",
        label: "Inventor",
        help: "Middle drag pans, Shift+middle orbits, Ctrl+middle zooms, right drag orbits",
        reverse_zoom: false,
        bindings: &[
            bind(Middle, NavKeys::SHIFT, Orbit),
            bind(Middle, NavKeys::CTRL, Zoom),
            bind(Middle, NavKeys::NONE, Pan),
            bind(Right, NavKeys::NONE, Orbit),
        ],
    },
    NavPreset {
        id: "alias",
        label: "Alias",
        help: "Shift+Alt+left drag orbits, Shift+Alt+middle pans, Shift+Alt+right zooms; middle drag also pans and right drag orbits",
        reverse_zoom: false,
        bindings: &[
            bind(Left, NavKeys::SHIFT_ALT, Orbit),
            bind(Middle, NavKeys::SHIFT_ALT, Pan),
            bind(Right, NavKeys::SHIFT_ALT, Zoom),
            bind(Middle, NavKeys::CTRL, Zoom),
            bind(Middle, NavKeys::NONE, Pan),
            bind(Right, NavKeys::NONE, Orbit),
        ],
    },
    NavPreset {
        id: "tinkercad",
        label: "Tinkercad",
        help: "Right drag orbits, middle or Shift+right drag pans",
        reverse_zoom: false,
        bindings: &[bind(Right, NavKeys::SHIFT, Pan), bind(Right, NavKeys::NONE, Orbit), bind(Middle, NavKeys::NONE, Pan)],
    },
    NavPreset {
        id: "openscad",
        label: "OpenSCAD",
        help: "Left drag orbits; right or middle drag pans, Shift+right or Shift+middle drag zooms",
        reverse_zoom: false,
        bindings: &[
            bind(Left, NavKeys::NONE, Orbit),
            bind(Right, NavKeys::SHIFT, Zoom),
            bind(Middle, NavKeys::SHIFT, Zoom),
            bind(Right, NavKeys::NONE, Pan),
            bind(Middle, NavKeys::NONE, Pan),
        ],
    },
];

/// The preset with this id (the default for an unknown one).
pub fn preset(id: &str) -> &'static NavPreset {
    NAV_PRESETS.iter().find(|p| p.id == id).unwrap_or(&NAV_PRESETS[0])
}

/// What a drag of `button` with `mods` held does under `preset`: the preset's binding; else, for
/// the left button, an orbit when only Ctrl (Cmd) is held and the drag starts off the model,
/// nothing in box selection (the drag draws the box), and the navigation bar's mode `bar`
/// otherwise. `None`: the drag does not move the camera.
pub fn drag_mode(preset: &NavPreset, button: NavButton, mods: egui::Modifiers, bar: NavMode, boxing: bool, on_model: bool) -> Option<NavMode> {
    if let Some(m) = preset.mode(button, mods) {
        return Some(m);
    }
    if button != Left {
        return None;
    }
    let ctrl_only = NavKeys::CTRL.held(mods) && !mods.shift && !mods.alt;
    if ctrl_only && !on_model {
        Some(Orbit)
    } else if boxing {
        None
    } else {
        Some(bar)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui::Modifiers;

    /// The bar holds SolidCraft's default slots where FieldCraft has them, in its order: orbit,
    /// pan, zoom, fit, home, the visual style, and the display menu at the end. The display menu
    /// holds the projection, as SolidCraft's does.
    #[test]
    fn the_bar_holds_solidcrafts_slots_ending_with_display() {
        let mut app = ModelApp::new(fieldcraft_engine::Session::default());
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let frame = |app: &mut ModelApp, events: Vec<egui::Event>| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(pos2(0.0, 0.0), vec2(1400.0, 900.0))), events, ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        };
        let draw = |app: &mut ModelApp| (0..3).for_each(|_| frame(app, Vec::new()));
        draw(&mut app);
        let at = |app: &ModelApp, k: &str| app.layout.get(k).copied().unwrap_or_else(|| panic!("{k} not drawn"));
        let order = ["orbit", "pan", "zoom", "fit", "home", "visual_style", "display"];
        let slots: Vec<Rect> = order.iter().map(|k| at(&app, region(k).unwrap_or_default())).collect();
        for (pair, keys) in slots.windows(2).zip(order.windows(2)) {
            assert!(pair[0].right() <= pair[1].left(), "{keys:?}: {pair:?}");
        }
        let bar = at(&app, "nav_bar");
        let last = slots.last().copied().unwrap_or(Rect::NOTHING);
        assert!(bar.right() - last.right() < 2.0 * (kit::NAV_BORDER + kit::NAV_PADDING) + 0.5, "display ends the bar: {last:?} in {bar:?}");
        assert_eq!(BAR.iter().filter(|s| s.1 != kit::NavigationKind::Separator).count(), order.len(), "no other slot");
        app.ui.menu = Some("view".into());
        draw(&mut app);
        let perspective = at(&app, "view_menu:perspective");
        app.ui.perspective = false;
        let press = |pressed| egui::Event::PointerButton {
            pos: perspective.center(),
            button: egui::PointerButton::Primary,
            pressed,
            modifiers: Modifiers::default(),
        };
        frame(&mut app, vec![egui::Event::PointerMoved(perspective.center())]);
        frame(&mut app, vec![press(true)]);
        frame(&mut app, vec![press(false)]);
        assert!(app.ui.perspective, "the display menu's Perspective row");
    }

    const SHIFT: Modifiers = Modifiers::SHIFT;
    const CTRL: Modifiers = Modifiers::CTRL;
    const NONE: Modifiers = Modifiers::NONE;

    /// Each preset maps its buttons and keys as its help says.
    #[test]
    fn presets_bind_as_listed() {
        let shift_alt = Modifiers { shift: true, alt: true, ..NONE };
        let cases: [(&str, &[(NavButton, Modifiers, Option<NavMode>)]); 7] = [
            ("fieldcraft", &[(Left, SHIFT, Some(Pan)), (Right, NONE, Some(Pan)), (Middle, NONE, Some(Pan)), (Left, NONE, None)]),
            ("fusion", &[(Middle, NONE, Some(Pan)), (Middle, SHIFT, Some(Orbit)), (Right, NONE, Some(Orbit)), (Left, NONE, None)]),
            ("solidworks", &[(Middle, NONE, Some(Orbit)), (Middle, CTRL, Some(Pan)), (Middle, SHIFT, Some(Zoom)), (Right, NONE, Some(Orbit))]),
            ("inventor", &[(Middle, NONE, Some(Pan)), (Middle, SHIFT, Some(Orbit)), (Middle, CTRL, Some(Zoom)), (Right, NONE, Some(Orbit))]),
            ("alias", &[(Left, shift_alt, Some(Orbit)), (Middle, shift_alt, Some(Pan)), (Right, shift_alt, Some(Zoom)), (Middle, NONE, Some(Pan))]),
            ("tinkercad", &[(Right, NONE, Some(Orbit)), (Right, SHIFT, Some(Pan)), (Middle, NONE, Some(Pan)), (Left, NONE, None)]),
            ("openscad", &[(Left, NONE, Some(Orbit)), (Right, NONE, Some(Pan)), (Right, SHIFT, Some(Zoom)), (Middle, SHIFT, Some(Zoom))]),
        ];
        for (id, binds) in cases {
            let p = preset(id);
            assert_eq!(p.id, id);
            for &(button, mods, want) in binds {
                assert_eq!(p.mode(button, mods), want, "{id}: {button:?} {mods:?}");
            }
        }
        assert_eq!(preset("no such preset").id, DEFAULT_PRESET);
        assert!(preset("solidworks").reverse_zoom);
    }

    /// An unbound left drag does the bar's mode; Ctrl+left off the model orbits in every preset
    /// (in box selection too); box selection keeps the plain left drag for the box.
    #[test]
    fn left_drag_follows_the_bar_and_ctrl_orbits_from_empty_space() {
        for p in &NAV_PRESETS {
            assert_eq!(drag_mode(p, Left, CTRL, Pan, true, false), Some(Orbit), "{}", p.id);
            assert_eq!(drag_mode(p, Left, CTRL, Pan, false, true), Some(p.mode(Left, CTRL).unwrap_or(Pan)), "{}", p.id);
        }
        let fc = preset(DEFAULT_PRESET);
        assert_eq!(drag_mode(fc, Left, NONE, Zoom, false, true), Some(Zoom));
        assert_eq!(drag_mode(fc, Left, NONE, Orbit, true, false), None);
        assert_eq!(drag_mode(fc, Left, SHIFT, Orbit, false, false), Some(Pan));
        assert_eq!(drag_mode(preset("fusion"), Right, CTRL, Orbit, false, false), Some(Orbit));
    }
}
