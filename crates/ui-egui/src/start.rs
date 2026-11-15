//! The start page: the kit's start page (`solidcraft_ui_kit::start_page`, as SolidCraft's Home
//! and SurfaceCraft's start page use it) in place of the workspace: a rail (back to the model,
//! New, Open…, Recover, Settings, Crash Reports, Report a Bug) beside the action tiles, the recent
//! files and the samples with their thumbnails. Shown at launch when nothing opens and from the
//! application bar's Home; Back (or Esc) returns to the model.
//!
//! Every entry runs what the File menu runs (`file.new`, the Open chooser or `file.open`'s form,
//! `file.sample`, `file.recover`, `deck.read`'s form…). The page stands over the document it was
//! opened on and goes once that document is replaced or changed, or a form opens: whatever its
//! entry started (a pick, an unsaved-changes question, a background open) has then landed in the
//! workspace. An open that is cancelled or fails leaves the document as it was, and the page up.

use std::borrow::Cow;

use serde_json::{Value, json};
use solidcraft_ui_kit::start_page::{self, Back, Body, Icon, Page, RailEntry, RailItem, Recent, Review, Reviews, Start, Tile, Tiles};

use crate::ModelApp;
use crate::files::{self, Kind};

/// The page's state: whether it is up, its section, and what it stands over.
#[derive(Default)]
pub struct StartPage {
    /// The document the page was opened over, while it is up.
    over: Option<Over>,
    section: Section,
    /// The autosaved models the Recover section lists, read from disk once it is shown.
    recoverable: Option<Vec<Value>>,
    /// Something of the window's was over the page when the last frame ended (a question, a
    /// menu, the command search): an Esc pressed now is that one's, not the page's.
    covered: bool,
}

/// What the page stands over: the open document ([`fieldcraft_engine::Session::opened`]), its
/// state ([`fieldcraft_engine::Session::revision`]) and the forms opened so far. Any of them
/// moving means the workspace has something new to show.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Over {
    document: u64,
    revision: u64,
    forms: u64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Section {
    #[default]
    Start,
    Recover,
}

/// What an entry of the page stands for.
#[derive(Clone, Debug, PartialEq)]
enum Action {
    Back,
    Section(Section),
    New,
    /// The Open chooser for files of this kind (any FieldCraft opens for `None`).
    Open(Option<Kind>),
    ReadDeck,
    Recent(String),
    Sample(&'static str),
    Recover(String),
    Discard(String),
    Settings,
    CrashReports,
    ReportBug,
}

impl ModelApp {
    /// Show the start page over the open model (the application bar's Home).
    pub fn open_start(&mut self) {
        self.start_page.over = Some(self.over());
        self.start_page.section = Section::Start;
        self.start_page.recoverable = None;
    }

    /// The start page at launch: shown unless the launch is opening something (a sample or a
    /// file named on the command line or in the page's address).
    pub fn launch_start(&mut self) {
        if self.loading.is_none() && self.queued.is_empty() {
            self.open_start();
        }
    }

    /// Back to the model.
    pub fn close_start(&mut self) {
        self.start_page = StartPage::default();
    }

    /// Whether the start page is up.
    pub fn start_shown(&self) -> bool {
        self.start_page.over.is_some()
    }

    fn over(&self) -> Over {
        Over { document: self.session.opened(), revision: self.session.revision(), forms: self.forms_opened }
    }

    /// Nothing the page started is still on its way: no chooser, question or background job.
    fn settled(&self) -> bool {
        self.loading.is_none() && self.queued.is_empty() && self.picker.is_none() && self.prompt.is_none()
    }
}

/// The page's rail rows and tiles, by the names the kit publishes them under, also kept in
/// `ui.inspect`'s layout (docs/control-protocol.md).
const REGIONS: &[&str] = &[
    "start:back",
    "start:start",
    "start:new",
    "start:open",
    "start:recover",
    "start:settings",
    "start:crash_reports",
    "start:report_bug",
    "start:tile:new",
    "start:tile:open",
    "start:tile:design",
    "start:tile:mesh",
    "start:tile:deck",
];

/// Draw the start page in `ui` when it is up; `false` when it is not, and the caller draws the
/// workspace.
pub fn page(app: &mut ModelApp, ui: &mut egui::Ui) -> bool {
    let Some(over) = app.start_page.over else { return false };
    if app.settled() && app.over() != over {
        app.close_start();
        return false;
    }
    let interactive = app.busy().is_none() && app.picker.is_none();
    let rect = ui.available_rect_before_wrap();
    app.place("start_page", rect);
    if app.start_page.section == Section::Recover && app.start_page.recoverable.is_none() {
        app.start_page.recoverable = Some(app.recoverable());
    }
    let mut child = ui.new_child(egui::UiBuilder::new().id_salt("fc_start_page").max_rect(rect));
    let chosen = start_page::show(&mut child, contents(app, ui.ctx(), interactive));
    for region in REGIONS {
        if let Some(r) = solidcraft_ui_registry::rect_of(region) {
            app.place(region, r);
        }
    }
    // Esc chooses Back. A stop-gap until the kit's start page binds Esc itself
    // (solidcraft#17391); delete it then, or it takes the key first. It runs after the kit has
    // drawn, so the kit's open rail menu (a narrow window) takes Esc first, and only when nothing
    // of the window's (a prompt, a menu, the command search, Settings) was over the page as the
    // last frame ended: a menu drawn before the page this frame has already closed on this Esc.
    let back = interactive && !app.start_page.covered && ui.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape));
    if let Some(action) = chosen.or(back.then_some(Action::Back)) {
        run(app, action);
    }
    app.start_page.covered = app.prompt.is_some() || app.ui.menu.is_some() || app.ui.palette || app.settings.open;
    true
}

/// The page as the kit draws it.
fn contents<'a>(app: &'a ModelApp, ctx: &egui::Context, interactive: bool) -> Page<'a, Action> {
    let key = crate::appbar::shortcut;
    let entry = |id: &'static str, label: &'static str, shortcut: String, selected: bool, chosen: Action| {
        RailItem::Entry(RailEntry { id: id.into(), label: label.into(), shortcut: shortcut.into(), selected, enabled: true, chosen })
    };
    let section = app.start_page.section;
    let rail = vec![
        entry("start:start", "Start", String::new(), section == Section::Start, Action::Section(Section::Start)),
        entry("start:new", "New", key("N"), false, Action::New),
        entry("start:open", "Open…", key("O"), false, Action::Open(None)),
        RailItem::Separator,
        entry("start:recover", "Recover", String::new(), section == Section::Recover, Action::Section(Section::Recover)),
        RailItem::Separator,
        entry("start:settings", "Settings…", String::new(), false, Action::Settings),
        entry("start:crash_reports", "Crash Reports…", String::new(), false, Action::CrashReports),
        entry("start:report_bug", "Report a Bug…", String::new(), false, Action::ReportBug),
    ];
    // An error (a file that would not open) stays under the hint until the next message.
    let error = app.status.as_ref().filter(|s| s.error).map(|s| Cow::from(s.text.as_str()));
    let (title, hint, body) = match section {
        Section::Start => (
            "Start",
            "Open a model, a SolidCraft design or a tagged mesh, or drop one here. The samples are ready to solve.",
            Body::Start(Box::new(start(app, ctx))),
        ),
        Section::Recover => (
            "Recover",
            "Models autosaved by a FieldCraft that crashed or closed with unsaved changes.",
            Body::Reviews(recover(app.start_page.recoverable.as_deref().unwrap_or_default())),
        ),
    };
    Page {
        id: "start",
        name: "Start".into(),
        interactive,
        back: Back { template: "Back to {name}".into(), document: app.model_name().into(), shortcut: "Esc".into(), chosen: Action::Back },
        rail,
        title: title.into(),
        hint: hint.into(),
        error,
        body,
    }
}

/// The start section: the ways to begin, the recent files and the samples.
fn start<'a>(app: &'a ModelApp, ctx: &egui::Context) -> Start<'a, Action> {
    let tile = |id: &'static str, icon: &'static str, title: &'static str, hint: String, primary: bool, chosen: Action| Tile {
        id: id.into(),
        icon: Icon::Kit(icon),
        title: title.into(),
        hint: hint.into(),
        primary,
        chosen,
    };
    let open_key = crate::appbar::shortcut("O");
    let actions = vec![
        tile("start:tile:new", "new", "New model", "Empty: add geometry, a mesh and the setup".into(), true, Action::New),
        tile("start:tile:open", "open", "Open…", format!("A model, a design or a mesh ({open_key})"), false, Action::Open(None)),
        tile("start:tile:design", "open", "Open design…", "A SolidCraft design, to mesh and set up".into(), false, Action::Open(Some(Kind::Design))),
        tile("start:tile:mesh", "import", "Import mesh…", "A face-tagged 3MF surface mesh".into(), false, Action::Open(Some(Kind::Mesh))),
        tile("start:tile:deck", "import", "Read deck…", "An Abaqus, Nastran or LS-DYNA input deck".into(), false, Action::ReadDeck),
    ];
    let recent = app.recent.paths.iter().map(|path| recent_row(path)).collect();
    let samples = crate::samples::ordered()
        .into_iter()
        .map(|s| start_page::Sample {
            title: s.label.into(),
            description: s.about.into(),
            icon: Icon::Kit(crate::samples::ICON),
            thumbnail: crate::samples::thumbnail(ctx, s.name).map(|t| (t.id(), t.size_vec2())),
            chosen: Action::Sample(s.name),
        })
        .collect();
    Start {
        banner: None,
        actions_title: "New".into(),
        actions: Tiles { grid: "fc_start_new", tiles: actions },
        recent_title: "Recent".into(),
        recent_empty: "Models, designs and meshes you open or save are listed here.".into(),
        recent,
        samples_title: "Samples".into(),
        samples,
    }
}

/// A recent file's row: its name, its folder and what it is.
fn recent_row(path: &str) -> Recent<'static, Action> {
    let folder = std::path::Path::new(path).parent().map(|p| p.display().to_string()).unwrap_or_default();
    let kind = match Kind::of(path) {
        Some(Kind::Model) => "Model",
        Some(Kind::Design) => "SolidCraft design",
        Some(Kind::Mesh) => "Tagged mesh",
        Some(Kind::Orphan) => "Orphan mesh",
        None => "",
    };
    Recent {
        title: files::file_name(path).into(),
        path: folder.into(),
        when: "".into(),
        kind: kind.into(),
        // FieldCraft reads no signing provenance from the files it lists.
        security: None,
        thumbnail: None,
        thumbnail_failed: false,
        chosen: Action::Recent(path.to_string()),
        menu: Vec::new(),
    }
}

/// The Recover section: each autosaved model with Recover and Discard.
fn recover(entries: &[Value]) -> Reviews<'static, Action> {
    let items = entries
        .iter()
        .filter_map(|e| {
            let id = e.get("id").and_then(Value::as_str)?.to_string();
            let n = |k: &str| e.get(k).and_then(Value::as_u64).unwrap_or(0);
            let title = e.get("path").and_then(Value::as_str).map_or_else(|| "Untitled".to_string(), files::file_name);
            let detail = format!("Autosaved {} · {} nodes, {} elements", files::ago(n("age")), n("nodes"), n("elements"));
            Some(Review {
                title: title.into(),
                detail: detail.into(),
                recover: ("Recover".into(), Action::Recover(id.clone())),
                discard: ("Discard".into(), Action::Discard(id)),
            })
        })
        .collect();
    Reviews { empty: "No autosaved models to recover.".into(), items }
}

/// Run the entry chosen. What opens a document or a form leaves the page up until it lands
/// ([`page`]).
fn run(app: &mut ModelApp, action: Action) {
    match action {
        Action::Back => app.close_start(),
        Action::Section(section) => {
            app.start_page.section = section;
            app.start_page.recoverable = None;
        }
        Action::New => app.new_model(),
        Action::Open(kind) => app.open_dialog(kind),
        Action::ReadDeck => {
            // A refusal says why in the status, which the page shows.
            let _ = app.start("deck.read");
        }
        Action::Recent(path) => app.open_path(&path),
        Action::Sample(name) => app.open_sample(name),
        Action::Recover(id) => app.recover(&id),
        Action::Discard(id) => {
            let _ = app.run("file.recovery.discard", json!({"id": id}));
            app.start_page.recoverable = None;
        }
        Action::Settings => app.open_settings(),
        Action::CrashReports => app.open_crash_reports(),
        Action::ReportBug => app.report_bug(None),
    }
}

#[cfg(test)]
mod tests {
    use fieldcraft_engine::Session;
    use serde_json::json;

    use super::*;

    const W: f32 = 1280.0;
    const H: f32 = 800.0;

    /// Draws `app` frame by frame as the window does (its logic, then its layout), in a window
    /// `size` points.
    struct Drive {
        ctx: egui::Context,
        time: f64,
        size: egui::Vec2,
    }

    impl Drive {
        fn new() -> Drive {
            Drive::sized(W, H)
        }

        fn sized(w: f32, h: f32) -> Drive {
            Drive { ctx: egui::Context::default(), time: 0.0, size: egui::vec2(w, h) }
        }

        fn frame(&mut self, app: &mut ModelApp, events: Vec<egui::Event>) {
            self.time += 0.1;
            let raw = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, self.size)),
                time: Some(self.time),
                events,
                ..Default::default()
            };
            let mut out = self.ctx.run_ui(raw, |ui| {
                app.logic(ui.ctx());
                app.ui(ui);
            });
            out.textures_delta.clear();
        }

        /// A key pressed and released with `modifiers`, then a frame to settle.
        fn key(&mut self, app: &mut ModelApp, key: egui::Key, modifiers: egui::Modifiers) {
            let event = |pressed| egui::Event::Key { key, physical_key: None, pressed, repeat: false, modifiers };
            self.frame(app, vec![egui::Event::ModifiersChanged(modifiers), event(true), event(false)]);
            self.frame(app, vec![egui::Event::ModifiersChanged(egui::Modifiers::NONE)]);
        }

        /// A click at `at`, as [`Drive::click`] makes it.
        fn click_at(&mut self, app: &mut ModelApp, at: egui::Pos2) {
            let button =
                |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
            self.frame(app, vec![egui::Event::PointerMoved(at)]);
            self.frame(app, vec![button(true)]);
            self.frame(app, vec![button(false)]);
        }

        /// A click where the kit published `name` last frame, as a pointer would: move, press,
        /// release, each a frame.
        fn click(&mut self, app: &mut ModelApp, name: &str) {
            let at = solidcraft_ui_registry::handle_at(name).unwrap_or_else(|| panic!("{name} was not drawn"));
            self.click_at(app, at);
        }

        /// Wheel the page's body until what the kit published as `name` is inside the page;
        /// whether it had to.
        fn scroll_to(&mut self, app: &mut ModelApp, name: &str) -> bool {
            for wheeled in 0..60 {
                let page = solidcraft_ui_registry::rect_of("start:page");
                let rect = solidcraft_ui_registry::rect_of(name);
                if let (Some(page), Some(rect)) = (page, rect)
                    && page.contains_rect(rect)
                {
                    return wheeled > 0;
                }
                let over = page.map_or(egui::Pos2::ZERO, |p| p.center() + egui::vec2(p.width() / 4.0, 0.0));
                let wheel = egui::Event::MouseWheel {
                    unit: egui::MouseWheelUnit::Point,
                    delta: egui::vec2(0.0, -80.0),
                    phase: egui::TouchPhase::Move,
                    modifiers: Default::default(),
                    source: egui::MouseWheelSource::Wheel,
                };
                self.frame(app, vec![egui::Event::PointerMoved(over), wheel]);
                self.frame(app, vec![]);
            }
            panic!("{name} never scrolls into the page");
        }
    }

    fn shown(app: &ModelApp) -> bool {
        app.layout.contains_key("start_page")
    }

    /// At launch with nothing to open, the page is the window below the application bar; the
    /// workspace (its viewport) is not drawn. A launch that opens a sample shows the sample.
    #[test]
    fn the_page_is_up_at_launch_unless_the_launch_opens_something() {
        let mut drive = Drive::new();
        let mut app = ModelApp::new(Session::default());
        app.launch_start();
        drive.frame(&mut app, vec![]);
        drive.frame(&mut app, vec![]);
        assert!(shown(&app) && !app.layout.contains_key("viewport"), "{:?}", app.layout.keys());
        let page = app.layout["start_page"];
        let bar = app.layout["app_bar"];
        assert!((page.top() - bar.bottom()).abs() < 0.5 && (page.width() - W).abs() < 0.5, "{page:?} under {bar:?}");

        for region in REGIONS {
            assert!(app.layout.get(region).is_some_and(|r| page.contains_rect(*r)), "{region} on the page");
        }

        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.launch_start();
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app) && app.layout.contains_key("viewport"));
    }

    /// Home brings the page up over an open model; Back returns to it as it was.
    #[test]
    fn home_opens_the_page_over_a_model_and_back_returns_to_it() {
        let mut drive = Drive::new();
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        drive.frame(&mut app, vec![]);
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app));
        let (document, revision) = (app.session.opened(), app.session.revision());
        drive.click(&mut app, "app_bar:home");
        drive.frame(&mut app, vec![]);
        assert!(shown(&app), "Home shows the start page");
        drive.click(&mut app, "start:back");
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app) && app.layout.contains_key("viewport"), "Back returns to the model");
        assert_eq!((app.session.opened(), app.session.revision()), (document, revision), "the model is as it was");
        // Esc does as Back does.
        app.open_start();
        drive.frame(&mut app, vec![]);
        let esc = |pressed| egui::Event::Key { key: egui::Key::Escape, physical_key: None, pressed, repeat: false, modifiers: Default::default() };
        drive.frame(&mut app, vec![esc(true), esc(false)]);
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app));
    }

    /// Esc goes to what is over the page first: with the File menu open, it closes the menu and
    /// the page stays; the next Esc returns to the model.
    #[test]
    fn esc_closes_a_menu_over_the_page_before_the_page() {
        let mut drive = Drive::new();
        let mut app = ModelApp::new(Session::default());
        app.open_start();
        drive.frame(&mut app, vec![]);
        drive.click(&mut app, "app_bar:file");
        drive.frame(&mut app, vec![]);
        assert_eq!(app.ui.menu.as_deref(), Some("file"));
        let esc = |pressed| egui::Event::Key { key: egui::Key::Escape, physical_key: None, pressed, repeat: false, modifiers: Default::default() };
        drive.frame(&mut app, vec![esc(true), esc(false)]);
        drive.frame(&mut app, vec![]);
        assert!(app.ui.menu.is_none() && shown(&app), "the menu closes, the page stays");
        drive.frame(&mut app, vec![esc(true), esc(false)]);
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app), "then Esc returns to the model");
    }

    /// Every sample's card is on the page and a click opens that sample, at 1280 × 800 and in a
    /// small window where the page must scroll to reach them (#245): with unsaved work (a
    /// material) it asks first, naming the sample, and the page stays up behind the question.
    #[test]
    fn every_sample_card_opens_its_sample() {
        for (w, h) in [(W, H), (900.0, 560.0)] {
            let mut scrolled = false;
            for (i, sample) in crate::samples::ordered().into_iter().enumerate() {
                let mut drive = Drive::sized(w, h);
                let mut app = ModelApp::new(Session::default());
                app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}))
                    .unwrap();
                app.open_start();
                drive.frame(&mut app, vec![]);
                drive.frame(&mut app, vec![]);
                let card = format!("start:sample:{i}");
                scrolled |= drive.scroll_to(&mut app, &card);
                drive.click(&mut app, &card);
                assert_eq!(app.prompt, Some(files::Prompt::Unsaved(files::Replace::Sample(sample.name))), "{w}×{h}: clicking {}", sample.name);
                drive.frame(&mut app, vec![]);
                assert!(shown(&app), "{w}×{h}: {}: the page waits for the answer", sample.name);
            }
            // The small window cannot show them all at once.
            if w < W {
                assert!(scrolled, "{w}×{h}");
            }
        }
    }

    /// In a narrow window (the web app) the kit folds the rail into a menu; Esc with that menu
    /// open closes the menu, and the page stays.
    #[test]
    fn esc_closes_the_narrow_rail_menu_before_the_page() {
        // Narrower than the kit's rail beside one tile column with its insets (start_page::show).
        let k = solidcraft_theme::Tokens::get();
        let narrow = k.file_nav_width + solidcraft_ui_kit::FileView::tile_min_width() + k.file_content_inset_x * 2.0 - k.space_24;
        let mut drive = Drive::sized(narrow, H);
        let mut app = ModelApp::new(Session::default());
        app.open_start();
        drive.frame(&mut app, vec![]);
        drive.frame(&mut app, vec![]);
        assert!(solidcraft_ui_registry::handle_at("start:back").is_none(), "the rail is folded into its menu");
        // The menu's button heads the page: the kit's rail frame (`FileView::navigation`) is
        // inset `space_14` from the page's top, and the button is a control high.
        let page = app.layout["start_page"];
        drive.click_at(&mut app, egui::pos2(page.center().x, page.top() + k.space_14 + k.control_height / 2.0));
        drive.frame(&mut app, vec![]);
        assert!(solidcraft_ui_registry::handle_at("start:back").is_some(), "the rail menu is open");
        drive.key(&mut app, egui::Key::Escape, egui::Modifiers::NONE);
        assert!(shown(&app), "Esc closes the menu, not the page");
        assert!(solidcraft_ui_registry::handle_at("start:back").is_none(), "the menu is closed");
        drive.key(&mut app, egui::Key::Escape, egui::Modifiers::NONE);
        assert!(!shown(&app), "the next Esc returns to the model");
    }

    /// A sample opened from the page replaces it once it has loaded; while it loads the page
    /// says so and takes no clicks.
    #[test]
    fn a_sample_opened_from_the_page_replaces_it_once_loaded() {
        let mut drive = Drive::new();
        let mut app = ModelApp::new(Session::default());
        app.open_start();
        drive.frame(&mut app, vec![]);
        let i = crate::samples::ordered().iter().position(|s| s.name == "cantilever").unwrap();
        let card = format!("start:sample:{i}");
        drive.scroll_to(&mut app, &card);
        drive.click(&mut app, &card);
        // Both at once: a later frame's logic may already have finished the load.
        assert!(app.busy().is_some(), "the sample opens in the background");
        assert!(shown(&app), "the page stays up while the sample loads");
        app.wait_loaded();
        drive.frame(&mut app, vec![]);
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app) && app.layout.contains_key("viewport"));
        let mut cantilever = Session::default();
        cantilever.execute("file.sample", &json!({"name": "cantilever"})).unwrap();
        assert_eq!(app.session.model.steps(), cantilever.model.steps(), "the cantilever is open");
    }

    /// New starts an empty model and shows it; Open… without a file chooser opens the open form
    /// in the workspace.
    #[test]
    fn new_and_open_land_in_the_workspace() {
        let mut drive = Drive::new();
        let mut app = ModelApp::new(Session::default());
        app.launch_start();
        drive.frame(&mut app, vec![]);
        drive.click(&mut app, "start:tile:new");
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app), "New shows the new model");

        app.launch_start();
        drive.frame(&mut app, vec![]);
        drive.click(&mut app, "start:tile:open");
        drive.frame(&mut app, vec![]);
        assert!(!shown(&app) && app.dialog.is_some(), "the open form is in the workspace");
    }

    /// A recent file that is no longer there leaves the page up, saying why.
    #[test]
    fn a_recent_file_that_fails_to_open_keeps_the_page_with_the_error() {
        let mut drive = Drive::new();
        let mut app = ModelApp::new(Session::default());
        let gone = std::env::temp_dir().join(format!("fieldcraft-start-{}-gone.fieldcraft", std::process::id()));
        app.recent.paths = vec![gone.display().to_string()];
        app.launch_start();
        drive.frame(&mut app, vec![]);
        drive.click(&mut app, "start:recent:0");
        app.wait_loaded();
        drive.frame(&mut app, vec![]);
        drive.frame(&mut app, vec![]);
        assert!(shown(&app), "the page stays");
        assert!(app.status.as_ref().is_some_and(|s| s.error), "{:?}", app.status);
    }
}
