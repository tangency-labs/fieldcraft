//! The one place windows and menus are made, so every prompt and menu shares the same frames:
//! the UI kit's windows, confirmations and menus (`solidcraft_ui_kit`).

use egui::Ui;
use solidcraft_ui_kit as kit;

/// How wide a window is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Width {
    /// Prompts and small panels.
    Normal,
    /// Windows holding a table or a list.
    Wide,
}

impl Width {
    /// The window's width in the kit's form cards.
    fn cards(self) -> f32 {
        match self {
            Width::Normal => 1.75,
            Width::Wide => 2.9,
        }
    }
}

/// A system window (the kit's): `icon` (an icon id) and `title` in its header with a close button,
/// centred on the screen, `width` wide. Footer actions are set with [`kit::Window::actions`].
pub fn window<'a>(id: egui::Id, icon: &'a str, title: &'a str, width: Width) -> kit::Window<'a> {
    kit::Window::new(id, icon, title).cards(width.cards())
}

/// Where a window or confirmation `id` was drawn last frame.
pub fn drawn(ctx: &egui::Context, id: egui::Id) -> Option<egui::Rect> {
    ctx.memory(|m| m.area_rect(id))
}

/// A menu hanging from `anchor` (its button), open while [`crate::UiState::menu`] is `key` (so
/// `ui.set {menu}` and offscreen captures open it too): the kit's popup and list surface, [`MENU_W`]
/// wide. A click on one of its rows, a click outside or Escape closes it. Returns where it was
/// drawn.
pub fn menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    list(app, ui, key, anchor, (MENU_W, Close::OnClick), add)
}

/// As [`menu`], [`WIDE_MENU_W`] wide, for a menu with notes and cards (Help).
pub fn wide_menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    list(app, ui, key, anchor, (WIDE_MENU_W, Close::OnClick), add)
}

/// As [`menu`], but clicks inside keep it open (settings with checkboxes).
pub fn settings_menu<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    list(app, ui, key, anchor, (MENU_W, Close::OnClickOutside), add)
}

fn list<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    (width, close): (f32, Close),
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::Rect> {
    hang(app, ui, key, anchor, close, |app, ui| kit::Menu::show(ui, width, |ui| add(app, ui))).map(|r| r.response.rect)
}

/// As [`menu`], for content that draws its own surface (the kit's workspace menu, the ribbon's
/// drop-downs); clicks inside keep it open, so its rows close it themselves (or `ui.close()`).
pub fn surface_popup<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::InnerResponse<R>> {
    hang(app, ui, key, anchor, Close::OnClickOutside, add)
}

/// What closes a menu besides Escape and a click outside it.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Close {
    /// Any click inside it too (a list of actions).
    OnClick,
    /// Only what its rows ask for ([`Ui::close`]).
    OnClickOutside,
}

/// The menu `key` as SolidCraft hangs its ribbon drop-downs and app bar menus: a kit
/// `Surface::Popup` area on the foreground whose top-left corner is `anchor`'s bottom-left
/// (moved to stay on screen), `add` drawing its surface. It closes on Escape, a click outside it
/// (but not the click on `anchor` that opened it: its button toggles it), a click inside for
/// [`Close::OnClick`], or [`Ui::close`] from a row.
fn hang<R>(
    app: &mut crate::ModelApp,
    ui: &Ui,
    key: &str,
    anchor: egui::Rect,
    close: Close,
    add: impl FnOnce(&mut crate::ModelApp, &mut Ui) -> R,
) -> Option<egui::InnerResponse<R>> {
    if app.ui.menu.as_deref() != Some(key) {
        return None;
    }
    let ctx = ui.ctx().clone();
    let shown = kit::Container::new(&ctx, kit::Surface::Popup)
        .area(egui::Id::unique(("fc_menu", key)))
        .fixed_pos(anchor.left_bottom())
        .constrain(true)
        .order(egui::Order::Foreground)
        .show(&ctx, |ui| {
            ui.scope_builder(egui::UiBuilder::new().closable(), |ui| {
                let inner = add(app, ui);
                (inner, ui.should_close())
            })
            .inner
        });
    let (inner, asked) = shown.inner;
    let (click, on_anchor) =
        ctx.input(|i| (i.pointer.any_click(), i.pointer.press_origin().or(i.pointer.interact_pos()).is_some_and(|p| anchor.contains(p))));
    let clicked = click && !on_anchor && (close == Close::OnClick || shown.response.clicked_elsewhere());
    if (asked || clicked || ctx.input(|i| i.key_pressed(egui::Key::Escape))) && app.ui.menu.as_deref() == Some(key) {
        app.ui.menu = None;
    }
    Some(egui::InnerResponse { inner, response: shown.response })
}

/// The menu a right click on `response` opens: the kit's popup, [`MENU_W`] wide, for `add`'s rows.
pub fn context_menu<R>(response: &egui::Response, add: impl FnOnce(&mut Ui) -> R) -> Option<R> {
    kit::Menu::context(response, |ui| {
        let t = solidcraft_theme::Tokens::get();
        ui.set_width((MENU_W - t.space_6 * 2.0 - 2.0).max(0.0));
        ui.spacing_mut().item_spacing.y = 0.0;
        add(ui)
    })
    .map(|r| r.inner)
}

/// A menu's width.
pub const MENU_W: f32 = kit::LIST_W;
/// A wide menu's width: Help, whose rows carry notes and cards. FieldCraft's own size; the kit
/// has list menus only ([`kit::LIST_W`]).
pub const WIDE_MENU_W: f32 = 360.0;

/// A menu row: its label, and its shortcut (or a note) at the right; dimmed and inert when
/// disabled. The kit's menu row; `id` names it in the widget registry (`menu:<id>`).
pub fn menu_row(ui: &mut Ui, id: &str, label: &str, right: Option<&str>, enabled: bool) -> egui::Response {
    let mut item = kit::MenuItem::new(id, label, "");
    item.shortcut = right.unwrap_or("");
    item.enabled = enabled;
    item.show(ui)
}

/// A menu section's small heading.
pub fn menu_heading(ui: &mut Ui, text: &str) {
    kit::Menu::heading(ui, text);
}

/// A rule between a menu's groups.
pub fn menu_separator(ui: &mut Ui) {
    kit::Menu::separator(ui);
}

#[cfg(test)]
mod tests {
    use crate::ModelApp;
    use egui::{Rect, pos2, vec2};
    use fieldcraft_engine::Session;

    const SCREEN: Rect = Rect::from_min_max(pos2(0.0, 0.0), pos2(1400.0, 900.0));
    const ANCHOR: Rect = Rect::from_min_max(pos2(300.0, 40.0), pos2(380.0, 60.0));

    /// One frame of `draw`, with these input events.
    fn frame(ctx: &egui::Context, events: Vec<egui::Event>, draw: impl FnMut(&mut egui::Ui)) {
        let raw = egui::RawInput { screen_rect: Some(SCREEN), events, ..Default::default() };
        let mut out = ctx.run_ui(raw, draw);
        out.textures_delta.clear();
    }

    fn press(at: egui::Pos2, pressed: bool) -> Vec<egui::Event> {
        let button = egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: egui::Modifiers::default() };
        vec![egui::Event::PointerMoved(at), button]
    }

    /// A list menu "k" under [`ANCHOR`] (a settings menu when `settings`): where it was drawn.
    fn show(app: &mut ModelApp, ui: &egui::Ui, settings: bool) -> Option<Rect> {
        let rows = |_: &mut ModelApp, ui: &mut egui::Ui| {
            super::menu_row(ui, "a", "Row A", None, true);
        };
        if settings { super::settings_menu(app, ui, "k", ANCHOR, rows) } else { super::menu(app, ui, "k", ANCHOR, rows) }
    }

    /// A menu hangs from its button, as SolidCraft's drop-downs do: the kit's popup surface with
    /// its top-left corner on the button's bottom-left corner, no gap.
    #[test]
    fn a_menu_hangs_from_its_anchor() {
        let mut app = ModelApp::new(Session::default());
        app.ui.menu = Some("k".into());
        let ctx = egui::Context::default();
        let mut drawn = None;
        for _ in 0..2 {
            frame(&ctx, Vec::new(), |ui| drawn = show(&mut app, ui, false));
        }
        let menu = drawn.unwrap_or_else(|| panic!("menu not drawn"));
        assert_eq!(menu.left_top(), ANCHOR.left_bottom(), "{menu:?}");
        // On the foreground, above the panels it covers.
        let layer = ctx.layer_id_at(menu.center());
        assert_eq!(layer.map(|l| l.order), Some(egui::Order::Foreground), "{layer:?}");
    }

    /// A menu near the window's right edge is moved left to stay on screen.
    #[test]
    fn a_menu_stays_on_screen() {
        let mut app = ModelApp::new(Session::default());
        app.ui.menu = Some("w".into());
        let ctx = egui::Context::default();
        let at = Rect::from_min_size(pos2(SCREEN.right() - 30.0, 10.0), vec2(24.0, 24.0));
        let mut drawn = None;
        for _ in 0..2 {
            frame(&ctx, Vec::new(), |ui| {
                drawn = super::wide_menu(&mut app, ui, "w", at, |_, ui| {
                    super::menu_row(ui, "a", "Row A", None, true);
                })
            });
        }
        let menu = drawn.unwrap_or_else(|| panic!("menu not drawn"));
        assert!(menu.right() <= SCREEN.right() + 0.5 && (menu.top() - at.bottom()).abs() < 0.5, "{menu:?}");
    }

    /// The click on its button that opens a menu leaves it open; a click elsewhere or Escape
    /// closes it; a click on a row closes a menu but not a settings menu.
    #[test]
    fn menus_close_as_menus_do() {
        let ctx = egui::Context::default();
        let open = |settings: bool| {
            let mut app = ModelApp::new(Session::default());
            for _ in 0..2 {
                frame(&ctx, Vec::new(), |ui| {
                    show(&mut app, ui, settings);
                });
            }
            app
        };
        let click = |app: &mut ModelApp, at: egui::Pos2, settings: bool, toggles: bool| {
            for down in [true, false] {
                frame(&ctx, press(at, down), |ui| {
                    if !down && toggles {
                        app.toggle_menu("k");
                    }
                    show(app, ui, settings);
                });
            }
            frame(&ctx, Vec::new(), |ui| {
                show(app, ui, settings);
            });
        };
        let row = ANCHOR.left_bottom() + vec2(40.0, 20.0);
        for settings in [false, true] {
            // Opened by its button's click: still open after that click.
            let mut app = open(settings);
            click(&mut app, ANCHOR.center(), settings, true);
            assert_eq!(app.ui.menu.as_deref(), Some("k"), "the opening click closed it");
            // A click elsewhere closes it.
            click(&mut app, pos2(900.0, 600.0), settings, false);
            assert_eq!(app.ui.menu, None, "a click elsewhere");
            // Escape closes it.
            app.ui.menu = Some("k".into());
            frame(&ctx, Vec::new(), |ui| {
                show(&mut app, ui, settings);
            });
            let escape =
                egui::Event::Key { key: egui::Key::Escape, physical_key: None, pressed: true, repeat: false, modifiers: egui::Modifiers::default() };
            frame(&ctx, vec![escape], |ui| {
                show(&mut app, ui, settings);
            });
            assert_eq!(app.ui.menu, None, "Escape");
            // A row's click.
            app.ui.menu = Some("k".into());
            click(&mut app, row, settings, false);
            assert_eq!(app.ui.menu.is_some(), settings, "a row's click, settings {settings}");
        }
    }

    /// Windows, popups and outer hosts are made by the kit only, so none can skip its frames:
    /// no raw egui window, modal or popup anywhere (menus are the kit's `Surface::Popup`
    /// containers, made here in frame.rs), and no raw egui panel, central panel or area anywhere
    /// in the app's code: every outer host is a `kit::Container` with a semantic `Surface`
    /// (SolidCraft docs/ui-container-audit.md).
    #[test]
    fn windows_and_hosts_come_from_the_kit() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let windows = [concat!("egui::Window", "::new("), concat!("Popup", "::new("), concat!("egui::Popup", "::"), concat!("egui::Modal", "::new(")];
        let hosts = [
            concat!("egui::Panel", "::"),
            concat!("egui::Central", "Panel"),
            concat!("egui::Area", "::new"),
            concat!("egui::Side", "Panel"),
            concat!("egui::TopBottom", "Panel"),
        ];
        let mut bad = Vec::new();
        for e in std::fs::read_dir(&dir).unwrap() {
            let p = e.unwrap().path();
            if p.extension().is_none_or(|x| x != "rs") {
                continue;
            }
            let text = std::fs::read_to_string(&p).unwrap();
            // The app's code, not the tests that name what they look for.
            let code = text.split("#[cfg(test)]").next().unwrap_or_default();
            for raw in windows.iter().filter(|raw| code.contains(*raw)) {
                bad.push(format!("{}: {raw}", p.display()));
            }
            bad.extend(hosts.iter().filter(|raw| code.contains(*raw)).map(|raw| format!("{}: {raw}", p.display())));
        }
        assert!(bad.is_empty(), "use crate::frame (the kit's windows and menus) or a kit::Container instead of {bad:?}");
    }
}
