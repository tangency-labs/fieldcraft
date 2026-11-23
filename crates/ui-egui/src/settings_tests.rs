use fieldcraft_engine::Session;
use fieldcraft_engine::usage::{self, TEST_LOCK};
use serde_json::json;
use solidcraft_ui_registry as registry;

use super::*;
use crate::files::{PrefsStore, Services};
use crate::privacy::UsageHost;

/// The window a test draws.
const SCREEN: egui::Vec2 = egui::vec2(1400.0, 1000.0);

/// A widget the frame drew: its accessible label, its hover text and its centre.
struct Drawn {
    label: String,
    hover: Option<String>,
    at: egui::Pos2,
}

/// A temporary preferences file, its folder removed when dropped.
struct Temp(std::path::PathBuf);

impl Temp {
    fn new(name: &str) -> Temp {
        let dir = std::env::temp_dir().join(format!("fieldcraft-settings-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Temp(dir.join("prefs.json"))
    }

    fn store(&self) -> PrefsStore {
        PrefsStore::file(self.0.clone())
    }

    fn saved(&self) -> Prefs {
        Prefs::load(Some(&self.store()))
    }
}

impl Drop for Temp {
    fn drop(&mut self) {
        if let Some(dir) = self.0.parent() {
            let _ = std::fs::remove_dir_all(dir);
        }
    }
}

/// An app whose preferences live in `prefs`, that may share usage and check for updates.
fn app(prefs: &Temp) -> ModelApp {
    let (_tx, events) = crate::inbox::channel();
    ModelApp::new(Session::default())
        .with_services(Services { prefs: Some(prefs.store()), ..Default::default() })
        .with_usage(UsageHost {
            platform: "desktop",
            os: "test",
            arch: "test",
            new_id: Box::new(|| Some("11111111-2222-4333-8444-555555555555".into())),
        })
        .with_updates(crate::updates::UpdateHost {
            unsupported: None,
            automatic: false,
            start: Box::new(|_| {}),
            events,
            now: Box::new(|| 0),
            open_file: Box::new(|_| Ok(())),
        })
}

/// One frame of the whole window with `events`; what it drew.
fn frame(app: &mut ModelApp, ctx: &egui::Context, events: Vec<egui::Event>) -> Vec<Drawn> {
    ctx.enable_accesskit();
    let raw = egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, SCREEN)), focused: true, events, ..Default::default() };
    let mut out = ctx.run_ui(raw, |ui| {
        app.logic(ui.ctx());
        app.ui(ui);
    });
    out.textures_delta.clear();
    let ppp = f64::from(ctx.pixels_per_point());
    out.platform_output
        .accesskit_update
        .iter()
        .flat_map(|u| u.nodes.iter().map(|(_, n)| n))
        .filter_map(|n| {
            let b = n.bounds()?;
            Some(Drawn {
                label: n.label().or(n.value()).unwrap_or_default().to_string(),
                hover: n.description().map(str::to_string),
                at: egui::pos2(((b.x0 + b.x1) / 2.0 / ppp) as f32, ((b.y0 + b.y1) / 2.0 / ppp) as f32),
            })
        })
        .collect()
}

fn frames(app: &mut ModelApp, ctx: &egui::Context, n: usize) -> Vec<Drawn> {
    (0..n).map(|_| frame(app, ctx, Vec::new())).last().unwrap_or_default()
}

/// A click at `at`, then frames for it to settle; what the last drew.
fn click_at(app: &mut ModelApp, ctx: &egui::Context, at: egui::Pos2) -> Vec<Drawn> {
    let press = |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: egui::Modifiers::NONE };
    frame(app, ctx, vec![egui::Event::PointerMoved(at), press(true)]);
    frame(app, ctx, vec![press(false)]);
    frames(app, ctx, 2)
}

/// Click what the kit published under `name`.
fn click(app: &mut ModelApp, ctx: &egui::Context, name: &str) -> Vec<Drawn> {
    let at = registry::rect_of(name).map(|r| r.center()).or_else(|| registry::handle_at(name)).unwrap_or_else(|| panic!("{name} is published"));
    click_at(app, ctx, at)
}

/// Click the widget labelled `label` inside the window.
fn click_label(app: &mut ModelApp, ctx: &egui::Context, drawn: &[Drawn], label: &str) -> Vec<Drawn> {
    let window = app.layout.get("settings").copied().expect("Settings is drawn");
    let at =
        drawn.iter().find(|d| d.label == label && window.contains(d.at)).map(|d| d.at).unwrap_or_else(|| panic!("{label} is not drawn in Settings"));
    click_at(app, ctx, at)
}

fn type_text(app: &mut ModelApp, ctx: &egui::Context, text: &str) {
    frame(app, ctx, vec![egui::Event::Text(text.into())]);
    frames(app, ctx, 2);
}

fn key(app: &mut ModelApp, ctx: &egui::Context, key: egui::Key) {
    let event = |pressed| egui::Event::Key { key, physical_key: None, pressed, repeat: false, modifiers: egui::Modifiers::NONE };
    frame(app, ctx, vec![event(true), event(false)]);
    frames(app, ctx, 2);
}

/// Every page is in the list and shows when its tab is clicked; a word from each page's settings
/// finds it from another page, outlines the match, and Enter opens the page with the row outlined.
#[test]
fn every_page_is_reachable_from_the_list_and_through_search() {
    let prefs = Temp::new("reach");
    let ctx = egui::Context::default();
    let mut app = app(&prefs);
    app.open_settings();
    frames(&mut app, &ctx, 3);
    for (index, page) in PAGES.iter().enumerate().rev() {
        let Page::Host(page) = page;
        click(&mut app, &ctx, page.handle);
        assert_eq!(app.settings.section, index, "{} shows its page", page.handle);
    }
    for (word, index, row) in
        [("units", 0, "Display units"), ("anonymous", 1, "Share anonymous usage"), ("automatically", 2, "Check for updates automatically")]
    {
        let other = (index + 1) % PAGES.len();
        app.settings.section = other;
        frames(&mut app, &ctx, 2);
        type_text(&mut app, &ctx, word);
        let (query, _, hits) = kit::settings_search::nav(&ctx);
        assert_eq!(query, word);
        let found: Vec<(usize, &str)> = hits.iter().map(|h| (h.page, h.label.as_str())).collect();
        assert!(found.iter().all(|(page, _)| *page == index), "{word} finds only its page's settings: {found:?}");
        let at = found.iter().position(|hit| *hit == (index, row)).unwrap_or_else(|| panic!("{word} finds {row}: {found:?}"));
        assert!(registry::handle_at("prefs:search:current").is_some(), "{word}: the match is outlined");
        for _ in 0..at {
            key(&mut app, &ctx, egui::Key::ArrowDown);
        }
        assert_eq!(kit::settings_search::nav(&ctx).1, at, "{word}: the arrows reach {row}");
        key(&mut app, &ctx, egui::Key::Enter);
        assert_eq!((app.settings.section, app.settings.filter.as_str()), (index, ""), "{word}: Enter opens the page");
        assert!(registry::handle_at("prefs:revealed").is_some(), "{word}: the opened row is outlined");
    }
}

/// Cancel puts back what the window found, through every setter: the units, History's dock and a
/// usage question never answered, changed by clicks and by `ui.settings`; on screen, in the
/// session, in recording and on disk.
#[test]
fn cancel_restores_the_snapshot_down_to_an_unanswered_question() {
    let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    usage::disable();
    let prefs = Temp::new("cancel");
    let ctx = egui::Context::default();
    let mut app = app(&prefs);
    assert_eq!(app.prefs.usage_sharing, None, "never asked");
    let before = Snapshot::of(&app);
    app.open_settings();
    frames(&mut app, &ctx, 3);
    let drawn = click(&mut app, &ctx, "prefs:page:Privacy");
    click_label(&mut app, &ctx, &drawn, "Share anonymous usage");
    assert_eq!(app.prefs.usage_sharing, Some(true), "the click shares");
    assert!(usage::enabled());
    from_control(&mut app, &json!({"set": {"units": "m-kg-s-N-Pa", "history_dock": "top", "update_check": false}})).unwrap();
    frames(&mut app, &ctx, 2);
    assert_ne!(Snapshot::of(&app), before);
    assert_eq!(prefs.saved().history_dock, HistoryDock::Top, "a change is kept as it is made");

    click(&mut app, &ctx, "prefs:cancel");
    assert!(!app.settings.open && registry::rect_of("prefs:window").is_none(), "Cancel closes");
    assert_eq!(Snapshot::of(&app), before, "Cancel puts everything back");
    assert_eq!(app.session.default_units, before.units, "the session shows the old units again");
    assert!(!usage::enabled(), "nothing is recorded again");
    let saved = prefs.saved();
    assert_eq!((saved.usage_sharing, saved.units, saved.history_dock, saved.update_check), (None, before.units, before.history_dock, true));
    app.ask_consent();
    assert_eq!(app.prompt, Some(crate::files::Prompt::Consent), "the question is still to be asked");
    usage::disable();
}

/// Apply keeps the changes and the window; the saved preferences, read back, equal the window's.
/// OK closes keeping them: changes are kept as they are made (the kit's OK does not call Apply).
#[test]
fn apply_and_ok_keep_the_changes_on_disk() {
    let prefs = Temp::new("apply");
    let ctx = egui::Context::default();
    let mut app = app(&prefs);
    app.open_settings();
    frames(&mut app, &ctx, 3);
    from_control(&mut app, &json!({"set": {"units": "mm-kg-ms-kN-GPa", "history_dock": "top", "update_check": false, "update_skip": "9.9.9"}}))
        .unwrap();
    frames(&mut app, &ctx, 2);
    assert!(window::apply_pending(&app));
    click(&mut app, &ctx, "prefs:apply");
    assert!(app.settings.open && !window::apply_pending(&app), "Apply keeps the window with nothing left to apply");
    let saved = prefs.saved();
    assert_eq!(saved, app.prefs, "what is saved reads back as the window has it");
    assert_eq!(
        (saved.units.preset(), saved.history_dock, saved.update_check, saved.update_skip.as_str()),
        (Preset::MmKgMs, HistoryDock::Top, false, "9.9.9")
    );

    let drawn = click(&mut app, &ctx, "prefs:page:General");
    click_label(&mut app, &ctx, &drawn, "Translucent browser");
    click(&mut app, &ctx, "prefs:ok");
    assert!(!app.settings.open);
    let saved = prefs.saved();
    assert!(saved.browser_float, "OK keeps a change clicked after the last Apply");
    let again = ModelApp::new(Session::default()).with_services(Services { prefs: Some(prefs.store()), ..Default::default() });
    assert_eq!(Snapshot::of(&again), Snapshot::of(&app), "the next run starts with them");
}

/// Restore Defaults puts back only the page shown.
#[test]
fn restore_defaults_resets_the_page_shown() {
    let prefs = Temp::new("defaults");
    let mut app = app(&prefs);
    from_control(&mut app, &json!({"page": "general", "set": {"units": "in-lbf-s-psi", "status_bar": true, "update_check": false}})).unwrap();
    from_control(&mut app, &json!({"press": "defaults"})).unwrap();
    let d = Prefs::default();
    assert_eq!((app.prefs.units.clone(), app.ui.status_bar, app.prefs.update_check), (d.units, d.status_bar, false));
}

/// With Settings open, typing goes to its search: the window's shortcuts (S, the command search)
/// do not act on the model behind it.
#[test]
fn typing_in_settings_searches_and_leaves_the_model_alone() {
    let prefs = Temp::new("typing");
    let ctx = egui::Context::default();
    let mut app = app(&prefs);
    app.open_settings();
    frames(&mut app, &ctx, 3);
    let s = |pressed| egui::Event::Key { key: egui::Key::S, physical_key: None, pressed, repeat: false, modifiers: egui::Modifiers::NONE };
    frame(&mut app, &ctx, vec![s(true), egui::Event::Text("s".into()), s(false)]);
    frames(&mut app, &ctx, 2);
    assert!(!app.ui.palette, "S does not open the command search over Settings");
    assert_eq!(app.settings.filter, "s");
}

/// No page shows JSON, a command id or a message id the catalogue lacks: every label and hover on
/// every page, and in a search over them all.
#[test]
fn no_page_shows_json_or_ids() {
    let ids: Vec<&str> = fieldcraft_engine::commands().iter().map(|c| c.id).collect();
    let id_like = |t: &str| t.contains('.') && t.split('.').all(|w| !w.is_empty() && w.chars().all(|c| c.is_ascii_lowercase() || c == '_'));
    let prefs = Temp::new("words");
    let ctx = egui::Context::default();
    let mut app = app(&prefs);
    app.prefs.update_skip = "9.9.9".into();
    app.open_settings();
    frames(&mut app, &ctx, 3);
    let mut texts = Vec::new();
    for Page::Host(page) in PAGES {
        let drawn = click(&mut app, &ctx, page.handle);
        let window = app.layout.get("settings").copied().expect("Settings is drawn");
        texts.extend(drawn.into_iter().filter(|d| window.contains(d.at)).flat_map(|d| std::iter::once(d.label).chain(d.hover)));
    }
    type_text(&mut app, &ctx, "e");
    let window = app.layout.get("settings").copied().expect("Settings is drawn");
    texts.extend(frames(&mut app, &ctx, 1).into_iter().filter(|d| window.contains(d.at)).flat_map(|d| std::iter::once(d.label).chain(d.hover)));
    for reached in ["Display units", "Translucent browser", "Share anonymous usage", "Check for updates automatically", "FieldCraft 9.9.9"] {
        assert!(texts.iter().any(|t| t == reached), "{reached} not drawn: {texts:?}");
    }
    let bad: Vec<&String> =
        texts.iter().filter(|t| t.contains('{') || t.contains('"') || id_like(t) || ids.iter().any(|id| t.contains(id))).collect();
    assert!(bad.is_empty(), "{bad:?}");
}

/// The window's gear and every page's icon have artwork in FieldCraft's set.
#[test]
fn the_window_and_its_pages_have_icons() {
    let names = std::iter::once("settings").chain(PAGES.iter().map(|Page::Host(p)| match p.icon {
        kit::Icon::Kit(name) => name,
        _ => "",
    }));
    for name in names {
        assert!(!crate::icons::SET.resolve(name).placeholder, "{name} has no artwork");
    }
}

/// `ui.settings` checks everything before it changes anything: an unknown page, button or
/// setting, a skipped version that is none, and usage in a run that cannot share are refused.
#[test]
fn hostile_settings_requests_change_nothing() {
    let prefs = Temp::new("hostile");
    let mut app = ModelApp::new(Session::default()).with_services(Services { prefs: Some(prefs.store()), ..Default::default() });
    assert!(from_control(&mut app, &json!({"press": "ok"})).is_err(), "nothing to press while closed");
    from_control(&mut app, &json!({})).unwrap();
    let before = Snapshot::of(&app);
    for bad in [
        json!({"page": "appearance"}),
        json!({"page": 3}),
        json!({"press": "close"}),
        json!({"set": {"units": "furlongs"}}),
        json!({"set": {"history_dock": "left", "status_bar": true}}),
        json!({"set": {"theme": "dark"}}),
        json!({"set": {"update_skip": "x".repeat(1 << 16)}}),
        json!({"set": {"usage_sharing": true, "status_bar": true}}),
    ] {
        assert!(from_control(&mut app, &bad).is_err(), "{bad} is refused");
        assert_eq!(Snapshot::of(&app), before, "{bad} changed nothing");
    }
    let reply = from_control(&mut app, &json!({"page": "updates"})).unwrap();
    assert_eq!((reply["open"].clone(), reply["page"].clone()), (json!(true), json!("updates")));
}

/// Settings opened from the start page closes on Esc and leaves the page up: the page's own Esc
/// (Back) waits while Settings is over it.
#[test]
fn esc_closes_settings_over_the_start_page_and_keeps_the_page() {
    let prefs = Temp::new("esc-start");
    let mut app = app(&prefs);
    let ctx = egui::Context::default();
    app.open_start();
    frames(&mut app, &ctx, 2);
    app.open_settings();
    frames(&mut app, &ctx, 2);
    key(&mut app, &ctx, egui::Key::Escape);
    frames(&mut app, &ctx, 2);
    assert!(!app.settings.open, "Esc closes Settings");
    assert!(app.start_shown(), "and the start page stays");
}

/// Restore Defaults never answers the first run's question: on Privacy, an answered question goes
/// back to unanswered and an unanswered one stays so; on Updates, the update choices go back to
/// theirs.
#[test]
fn restore_defaults_never_answers_the_usage_question() {
    let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    usage::disable();
    let prefs = Temp::new("defaults-usage");
    let mut app = app(&prefs);
    assert_eq!(app.prefs.usage_sharing, None, "never asked");
    from_control(&mut app, &json!({"page": "privacy"})).unwrap();
    from_control(&mut app, &json!({"press": "defaults"})).unwrap();
    assert_eq!(app.prefs.usage_sharing, None, "an unanswered question stays unanswered");
    from_control(&mut app, &json!({"set": {"usage_sharing": true, "usage_automation": false}})).unwrap();
    from_control(&mut app, &json!({"press": "defaults"})).unwrap();
    let d = Prefs::default();
    assert_eq!((app.prefs.usage_sharing, app.prefs.usage_automation), (None, d.usage_automation), "Defaults unanswer, never answer");
    assert_eq!(prefs.saved().usage_sharing, None, "and keep it so");
    from_control(&mut app, &json!({"page": "updates", "set": {"update_check": false, "update_skip": "9.9.9"}})).unwrap();
    from_control(&mut app, &json!({"press": "defaults"})).unwrap();
    assert_eq!((app.prefs.update_check, app.prefs.update_skip.as_str()), (d.update_check, d.update_skip.as_str()), "Updates' defaults");
    usage::disable();
}

/// Changing only what is counted (`usage_automation`) leaves the first run's question unanswered.
#[test]
fn counting_scripts_does_not_answer_the_usage_question() {
    let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    usage::disable();
    let prefs = Temp::new("automation-only");
    let mut app = app(&prefs);
    assert_eq!(app.prefs.usage_sharing, None, "never asked");
    from_control(&mut app, &json!({"page": "privacy"})).unwrap();
    from_control(&mut app, &json!({"set": {"usage_automation": false}})).unwrap();
    assert_eq!((app.prefs.usage_sharing, app.prefs.usage_automation), (None, false), "the question stays open");
    assert!(!usage::enabled(), "nothing is recorded");
    usage::disable();
}

/// Selecting Hebrew (or legacy `iw`) in Settings persists `"he"`, renders Hebrew headings,
/// restores `"en"` on Cancel, and rejects unknown language tags without mutating state.
#[test]
fn language_preference_switches_to_hebrew_and_restores_on_cancel() {
    let _g = TEST_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let prefs = Temp::new("language");
    let ctx = egui::Context::default();
    let mut app = ModelApp::new(Session::default()).with_services(Services { prefs: Some(prefs.store()), ..Default::default() });
    assert_eq!(app.prefs.language, "en");
    let before = Snapshot::of(&app);

    app.open_settings();
    frames(&mut app, &ctx, 2);
    assert!(from_control(&mut app, &json!({"set": {"language": "xx"}})).is_err());
    assert_eq!(Snapshot::of(&app), before);

    from_control(&mut app, &json!({"set": {"language": "iw"}})).unwrap();
    assert_eq!(app.prefs.language, "he");
    assert_eq!(prefs.saved().language, "he");
    let drawn = frames(&mut app, &ctx, 2);
    let he_heading = {
        let _lang = crate::i18n::scope("he");
        solidcraft_ui_text::tr("settings.general").to_owned()
    };
    assert!(drawn.iter().any(|d| d.label == he_heading), "expected Hebrew heading {he_heading:?} in drawn labels");

    click(&mut app, &ctx, "prefs:cancel");
    assert_eq!(app.prefs.language, "en");
    assert_eq!(prefs.saved().language, "en");
}
