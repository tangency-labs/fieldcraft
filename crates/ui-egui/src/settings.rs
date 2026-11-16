//! Settings: SolidCraft's Preferences window from the kit ([`kit::prefs_window`]) over
//! FieldCraft's own pages. The kit draws the window (page list, search with outline, Restore
//! Defaults, OK / Cancel / Apply over a snapshot); [`ModelApp`] is its [`Host`]. The pages are
//! drawn only with the kit's `settings_*` helpers, in FieldCraft's words (`locales/en.json`,
//! `settings.*`), and are grouped as SolidCraft's window groups its pages: General first on its
//! own, then a group of FieldCraft's own (Privacy, Updates). SolidCraft's Appearance, Layout,
//! Graphics and Navigation pages join from the kit once their bodies and words are kit code.
//!
//! A change takes effect and is kept as it is made, through the same setters the rest of the
//! window uses ([`ModelApp::set_default_units`], [`ModelApp::set_usage`],
//! [`ModelApp::set_update_prefs`], and the preferences [`ModelApp::save_prefs`] keeps); Cancel
//! puts the snapshot back through them.

use fieldcraft_model::DisplaySystem;
use fieldcraft_model::units::Preset;
use serde_json::{Value, json};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::BackingTooltipExt as _;
use solidcraft_ui_kit::prefs_window::{self as window, Button, Host, HostPage, Labels, Options, Page, State};
use solidcraft_ui_text::tr;

use crate::ModelApp;
use crate::files::Prefs;
use crate::shell::HistoryDock;

/// The prefix of every name the window publishes (`prefs:window`, `prefs:page:General`, …): the
/// kit's settings helpers publish theirs under `prefs:` too.
pub const PREFIX: &str = "prefs";

/// The pages, in the list's order. A page's control-channel name is its handle's last part,
/// lower-cased (`general`, `privacy`, `updates`).
pub const PAGES: &[Page<ModelApp>] = &[
    Page::Host(HostPage {
        icon: kit::Icon::Kit("prefs_general"),
        label: "settings.general",
        group: None,
        handle: "prefs:page:General",
        scroll: true,
        draw: general,
        defaults: general_defaults,
    }),
    Page::Host(HostPage {
        icon: kit::Icon::Kit("shield"),
        label: "settings.privacy",
        group: Some("settings.group_app"),
        handle: "prefs:page:Privacy",
        scroll: true,
        draw: privacy,
        defaults: privacy_defaults,
    }),
    Page::Host(HostPage {
        icon: kit::Icon::Kit("update"),
        label: "settings.updates",
        group: None,
        handle: "prefs:page:Updates",
        scroll: true,
        draw: updates,
        defaults: updates_defaults,
    }),
];

/// Everything the pages change, as Cancel puts it back (and `ui.settings` reports it).
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct Snapshot {
    units: DisplaySystem,
    history_dock: HistoryDock,
    browser_float: bool,
    browser_eye_left: bool,
    status_bar: bool,
    /// `None` until the first run's question is answered.
    usage_sharing: Option<bool>,
    usage_automation: bool,
    update_check: bool,
    update_skip: String,
}

impl Snapshot {
    fn of(app: &ModelApp) -> Snapshot {
        Snapshot {
            units: app.prefs.units.clone(),
            history_dock: app.ui.history_dock,
            browser_float: app.ui.browser_float,
            browser_eye_left: app.ui.browser_eye_left,
            status_bar: app.ui.status_bar,
            usage_sharing: app.prefs.usage_sharing,
            usage_automation: app.prefs.usage_automation,
            update_check: app.prefs.update_check,
            update_skip: app.prefs.update_skip.clone(),
        }
    }

    /// Put these settings back, each through its setter, and keep them.
    fn put_back(&self, app: &mut ModelApp) {
        if app.prefs.units != self.units {
            app.set_default_units(self.units.clone());
        }
        app.ui.history_dock = self.history_dock;
        app.ui.browser_float = self.browser_float;
        app.ui.browser_eye_left = self.browser_eye_left;
        app.ui.status_bar = self.status_bar;
        if (app.prefs.usage_sharing, app.prefs.usage_automation) != (self.usage_sharing, self.usage_automation) {
            app.put_usage(self.usage_sharing, self.usage_automation);
        }
        if (app.prefs.update_check, app.prefs.update_skip.as_str()) != (self.update_check, self.update_skip.as_str()) {
            app.set_update_prefs(self.update_check, Some(self.update_skip.clone()));
        }
        app.save_prefs();
    }
}

impl Host for ModelApp {
    type Snapshot = Snapshot;

    fn window(&mut self) -> &mut State<Snapshot> {
        &mut self.settings
    }

    fn window_ref(&self) -> &State<Snapshot> {
        &self.settings
    }

    fn take_snapshot(&self) -> Snapshot {
        Snapshot::of(self)
    }

    fn matches(&self, snapshot: &Snapshot) -> bool {
        Snapshot::of(self) == *snapshot
    }

    fn restore(&mut self, snapshot: &Snapshot) -> Result<(), String> {
        snapshot.put_back(self);
        Ok(())
    }

    /// Every change is kept as it is made; Apply makes sure the window's own preferences are too.
    fn apply(&mut self) -> Result<(), String> {
        self.save_prefs();
        Ok(())
    }
}

impl ModelApp {
    /// Open Settings at the page it last showed (the app bar's gear, Help › Settings…).
    pub fn open_settings(&mut self) {
        let page = self.settings.section;
        self.settings.open_at(page, PAGES.len());
        self.ui.menu = None;
    }
}

/// The window's frame and words.
fn options() -> Options<'static, ModelApp> {
    Options {
        id: egui::Id::unique("fc_settings"),
        prefix: PREFIX,
        labels: Labels {
            title: tr("settings.title"),
            icon: "settings",
            ok: tr("ui.ok"),
            cancel: tr("ui.cancel"),
            apply: tr("ui.apply"),
            search: tr("settings.search"),
            no_results: tr("settings.no_results"),
            restore_defaults: tr("settings.restore_defaults"),
            restore_defaults_hint: tr("settings.restore_defaults_hint"),
        },
        sidebar_footer: None,
        above_page: None,
    }
}

/// Draw Settings while it is open; where it stands is `ui.inspect`'s `settings` region.
pub fn show(app: &mut ModelApp, ctx: &egui::Context) {
    window::show(ctx, app, PAGES, &options(), |_| {});
    if app.settings.open
        && let Some(rect) = solidcraft_ui_registry::rect_of(&format!("{PREFIX}:window"))
    {
        app.place("settings", rect);
    }
}

fn page(index: usize) -> Option<&'static HostPage<ModelApp>> {
    PAGES.get(index).map(|Page::Host(page)| page)
}

/// A page's control-channel name: its handle's last part, lower-cased.
fn page_name(index: usize) -> Option<String> {
    page(index).and_then(|p| p.handle.rsplit(':').next()).map(str::to_lowercase)
}

/// Settings for `ui.inspect` and `ui.settings`.
pub fn describe(app: &ModelApp) -> Value {
    json!({
        "open": app.settings.open,
        "page": page_name(app.settings.section),
        "apply_pending": window::apply_pending(app),
        "usage_available": app.usage_available(),
        "updates_available": app.updates_available(),
        "values": Snapshot::of(app),
    })
}

/// A change `ui.settings` makes, as a page's controls make it.
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct Change {
    units: Option<DisplaySystem>,
    history_dock: Option<HistoryDock>,
    browser_float: Option<bool>,
    browser_eye_left: Option<bool>,
    status_bar: Option<bool>,
    usage_sharing: Option<bool>,
    usage_automation: Option<bool>,
    update_check: Option<bool>,
    update_skip: Option<String>,
}

/// `ui.settings`: open Settings (at `page`), change what its pages change (`set`), or press a
/// footer button (`press`). Everything is checked before anything changes.
pub fn from_control(app: &mut ModelApp, p: &Value) -> Result<Value, String> {
    let at = match p.get("page") {
        None => None,
        Some(name) => {
            let name = name.as_str().ok_or("`page` must be a page name: general, privacy or updates")?;
            Some(
                (0..PAGES.len())
                    .find(|i| page_name(*i).as_deref() == Some(name))
                    .ok_or_else(|| format!("no Settings page {name:?} (general, privacy, updates)"))?,
            )
        }
    };
    let change = match p.get("set") {
        None => None,
        Some(set) => Some(serde_json::from_value::<Change>(set.clone()).map_err(|e| format!("`set`: {e}"))?),
    };
    let press = match p.get("press").map(|b| b.as_str()) {
        None => None,
        Some(Some("ok")) => Some(Button::Ok),
        Some(Some("cancel")) => Some(Button::Cancel),
        Some(Some("apply")) => Some(Button::Apply),
        Some(Some("defaults")) => Some(Button::Defaults),
        Some(_) => return Err("`press` must be ok, cancel, apply or defaults".into()),
    };
    if let Some(change) = &change {
        let usage = change.usage_sharing.is_some() || change.usage_automation.is_some();
        if usage && !app.usage_available() {
            return Err("this run cannot share usage (no analytics key, or FIELDCRAFT_TELEMETRY=0, CI or automation)".into());
        }
        if let Some(skip) = &change.update_skip
            && !skip.is_empty()
            && fieldcraft_update::Version::parse(skip).is_none()
        {
            return Err(format!("`update_skip` must be a version such as 1.2.3, or \"\" for none: {skip:?}"));
        }
    }
    if (change.is_some() || press.is_some()) && !app.settings.open && at.is_none() {
        return Err("Settings is not open (`ui.settings {}` opens it)".into());
    }
    if let Some(index) = at {
        app.settings.clear_search();
        app.settings.open_at(index, PAGES.len());
        app.ui.menu = None;
    } else if change.is_none() && press.is_none() {
        app.open_settings();
    }
    // The snapshot Cancel puts back is what the window found, before this request's changes.
    window::hold(app);
    if let Some(change) = change {
        change_to(app, change);
    }
    if let Some(button) = press {
        window::press(app, PAGES, button)?;
    }
    Ok(describe(app))
}

fn change_to(app: &mut ModelApp, c: Change) {
    if let Some(units) = c.units {
        app.set_default_units(units);
    }
    if let Some(dock) = c.history_dock {
        app.ui.history_dock = dock;
    }
    if let Some(on) = c.browser_float {
        app.ui.browser_float = on;
    }
    if let Some(on) = c.browser_eye_left {
        app.ui.browser_eye_left = on;
    }
    if let Some(on) = c.status_bar {
        app.ui.status_bar = on;
    }
    // Only an explicit `usage_sharing` answers the first run's question; changing what is counted
    // leaves the question as it stands.
    match (c.usage_sharing, c.usage_automation) {
        (Some(sharing), automation) => app.set_usage(sharing, automation.unwrap_or(app.prefs.usage_automation)),
        (None, Some(automation)) => app.put_usage(app.prefs.usage_sharing, automation),
        (None, None) => {}
    }
    if c.update_check.is_some() || c.update_skip.is_some() {
        app.set_update_prefs(c.update_check.unwrap_or(app.prefs.update_check), c.update_skip);
    }
    app.save_prefs();
}

/// A checkbox with its hover; true when it was clicked.
fn check(ui: &mut egui::Ui, label: &str, hint: &str, on: &mut bool) -> bool {
    kit::checkbox(ui, tr(label), on).backing_on_hover_text(tr(hint)).changed()
}

/// A dim paragraph under a page's settings.
fn note(ui: &mut egui::Ui, id: &str) {
    kit::text(ui, tr(id), kit::Text::Note(id));
}

/// General: the display units, where History stands, and the window's own cards and bars.
fn general(ui: &mut egui::Ui, app: &mut ModelApp) {
    kit::settings_page(ui, |ui| {
        kit::settings_heading(ui, tr("settings.general"), tr("settings.general_note"));
        kit::settings_section(ui, "units", tr("settings.units"), |ui| {
            let labels: Vec<&str> = Preset::ALL.iter().map(|p| p.label()).collect();
            let current = app.prefs.units.preset();
            let mut selected = Preset::ALL.iter().position(|p| *p == current).unwrap_or(0);
            kit::row(ui, tr("settings.units_display"), |ui| {
                kit::dropdown(ui, "fc_settings_units", &mut selected, &labels).backing_on_hover_text(tr("settings.units_display_hint"))
            });
            if let Some(preset) = Preset::ALL.get(selected).copied().filter(|p| *p != current) {
                app.set_default_units(DisplaySystem::new(preset));
            }
        });
        kit::settings_section(ui, "history", tr("settings.history"), |ui| {
            const DOCKS: [HistoryDock; 2] = [HistoryDock::Bottom, HistoryDock::Top];
            let options = [tr("settings.history_bottom"), tr("settings.history_top")];
            let mut selected = DOCKS.iter().position(|d| *d == app.ui.history_dock).unwrap_or(0);
            let label = tr("settings.history_dock");
            let changed = kit::row(ui, label, |ui| kit::settings_segments(ui, label, &mut selected, &options));
            if changed && let Some(dock) = DOCKS.get(selected) {
                app.ui.history_dock = *dock;
            }
            note(ui, "settings.history_dock_hint");
        });
        kit::settings_section(ui, "window", tr("settings.window"), |ui| {
            check(ui, "settings.browser_float", "settings.browser_float_hint", &mut app.ui.browser_float);
            check(ui, "settings.browser_eye_left", "settings.browser_eye_left_hint", &mut app.ui.browser_eye_left);
            check(ui, "settings.status_bar", "settings.status_bar_hint", &mut app.ui.status_bar);
        });
    });
}

fn general_defaults(app: &mut ModelApp) {
    let d = Prefs::default();
    if app.prefs.units != d.units {
        app.set_default_units(d.units);
    }
    app.ui.history_dock = d.history_dock;
    app.ui.browser_float = d.browser_float;
    app.ui.browser_eye_left = d.browser_eye_left;
    app.ui.status_bar = d.status_bar;
    app.save_prefs();
}

/// Privacy: anonymous usage sharing (docs/privacy.md).
fn privacy(ui: &mut egui::Ui, app: &mut ModelApp) {
    kit::settings_page(ui, |ui| {
        kit::settings_heading(ui, tr("settings.privacy"), tr("settings.privacy_note"));
        kit::settings_section(ui, "usage", tr("settings.usage"), |ui| {
            let available = app.usage_available();
            let (mut sharing, mut automation) = (app.prefs.usage_sharing == Some(true), app.prefs.usage_automation);
            let mut changed = false;
            ui.add_enabled_ui(available, |ui| {
                changed |= check(ui, "settings.usage_share", "settings.usage_share_hint", &mut sharing);
            });
            ui.add_enabled_ui(available && sharing, |ui| {
                changed |= check(ui, "settings.usage_automation", "settings.usage_automation_hint", &mut automation);
            });
            if changed && available {
                app.set_usage(sharing, automation);
            }
            if !available {
                note(ui, "settings.usage_unavailable");
            }
            note(ui, "settings.usage_detail");
            note(ui, "settings.usage_local");
            let link =
                kit::button(ui, tr("settings.usage_what_is_sent"), kit::Button::Ghost).backing_on_hover_text(tr("settings.usage_what_is_sent_hint"));
            if link.clicked() {
                app.open_link(ui.ctx(), crate::privacy::PRIVACY_URL);
            }
        });
    });
}

/// The defaults never answer the first run's question for the user: sharing goes back to
/// unanswered (the question is asked again at the next start), what is counted to its default.
/// A run that cannot share leaves the choice alone.
fn privacy_defaults(app: &mut ModelApp) {
    if app.usage_available() {
        let d = Prefs::default();
        app.put_usage(d.usage_sharing, d.usage_automation);
    }
}

/// Updates: the automatic check and the skipped version (docs/updates.md).
fn updates(ui: &mut egui::Ui, app: &mut ModelApp) {
    kit::settings_page(ui, |ui| {
        kit::settings_heading(ui, tr("settings.updates"), tr("settings.updates_note"));
        let available = app.updates_available();
        kit::settings_section(ui, "check", tr("settings.updates_check"), |ui| {
            let mut on = app.prefs.update_check;
            let changed =
                ui.add_enabled_ui(available, |ui| check(ui, "settings.updates_automatic", "settings.updates_automatic_hint", &mut on)).inner;
            if changed {
                app.set_update_prefs(on, None);
            }
            note(ui, if available { "settings.updates_request" } else { "settings.updates_unavailable" });
            let now = ui.add_enabled_ui(available, |ui| {
                kit::button(ui, tr("settings.updates_now"), kit::Button::Normal).backing_on_hover_text(tr("settings.updates_now_hint"))
            });
            if now.inner.clicked() {
                app.check_for_updates();
            }
        });
        kit::settings_section(ui, "skipped", tr("settings.updates_skipped"), |ui| {
            if app.prefs.update_skip.is_empty() {
                note(ui, "settings.updates_skipped_none");
                return;
            }
            let version = solidcraft_ui_text::format("settings.updates_skipped_version", &[("version", app.prefs.update_skip.as_str())]);
            let again = kit::row(ui, tr("settings.updates_skipped"), |ui| {
                kit::text(ui, &version, kit::Text::Inline).backing_on_hover_text(tr("settings.updates_skipped_hint"));
                kit::button(ui, tr("settings.updates_offer_again"), kit::Button::Normal).clicked()
            });
            if again {
                app.set_update_prefs(app.prefs.update_check, Some(String::new()));
            }
        });
    });
}

fn updates_defaults(app: &mut ModelApp) {
    let d = Prefs::default();
    app.set_update_prefs(d.update_check, Some(d.update_skip));
}

#[cfg(test)]
#[path = "settings_tests.rs"]
mod tests;
