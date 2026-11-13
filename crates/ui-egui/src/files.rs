//! Files the user opens and saves: the host's native file pickers ([`Services`]), opening off the
//! window thread, drag and drop, the recent-files list, and the prompts around them (unsaved
//! changes, autosaved models to recover).
//!
//! Only user actions reach a picker: a click (toolbar, title bar, File menu), a shortcut or a drop.
//! Programmatic calls (`engine.execute`, the control channel's `ui.start`) run commands with the
//! paths they are given and never open a native dialog. Without a picker service (tests, the
//! offscreen capture, the web), the same entries open their command's parameter form instead.
//!
//! A file opens on a copy of the session in a worker thread; the window keeps drawing and shows
//! the file's progress with Cancel. The finished copy replaces the session only when the command
//! succeeded, so a failed or cancelled open leaves the model as it was. Structure follows
//! SolidCraft's `crates/ui-egui/src/file_load.rs` and `documents.rs` (MIT OR Apache-2.0), see
//! ATTRIBUTION.md.

use std::future::Future;
use std::path::{Path, PathBuf};
use std::pin::Pin;
use std::sync::Arc;
use std::sync::mpsc::{Sender, TryRecvError};
use std::task::{Context, Poll, Waker};

use fieldcraft_model::progress::Progress;

use fieldcraft_engine::Session;
use fieldcraft_engine::catalogue::Workspace;
use serde_json::{Value, json};
use solidcraft_ui_kit as kit;

use crate::inbox::Inbox;
use crate::{ModelApp, Prebuilt};
use solidcraft_ui_kit::BackingTooltipExt as _;

/// A file type a picker offers: a name and its extensions (without the dot).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Filter {
    pub name: &'static str,
    pub extensions: &'static [&'static str],
}

/// The host picks a file to open from the filters (the first is selected), or `None` when the
/// user cancelled.
pub type PickOpen = Box<dyn Fn(&[Filter]) -> Option<String>>;
/// The host picks where to save, starting from a file name; `None` when the user cancelled.
pub type PickSave = Box<dyn Fn(&str, &[Filter]) -> Option<String>>;
/// An owned request for a host's asynchronous file chooser.
#[derive(Debug, PartialEq, Eq)]
pub enum FilePicker {
    Open { filters: Vec<Filter> },
    Save { name: String, filters: Vec<Filter> },
}
/// Hosts return promptly; the future is polled on the UI thread without waiting.
pub type PickFile = Box<dyn Fn(FilePicker) -> Pin<Box<dyn Future<Output = Option<String>>>>>;
/// The host opens a web link in the user's browser.
pub type OpenUrl = Box<dyn Fn(&str) -> Result<(), String>>;

/// What the host provides: native file pickers and where the recent-files list is kept. The
/// desktop app installs them; tests answer the pickers with closures.
#[derive(Default)]
pub struct Services {
    /// Preferred over the synchronous services (which remain useful for the browser and tests).
    pub pick_file: Option<PickFile>,
    pub pick_open: Option<PickOpen>,
    pub pick_save: Option<PickSave>,
    /// The recent-files list's file (JSON); `None` keeps the list for this run only.
    pub recent_path: Option<PathBuf>,
    /// Where the user's preferences ([`Prefs`]) are kept: a file on the desktop, the browser's
    /// storage on the web; `None` keeps them for this run only.
    pub prefs: Option<PrefsStore>,
    /// Opens a link in the browser (Help › Report a Bug); without it, egui's own
    /// `open_url` (the browser build opens the link itself).
    pub open_url: Option<OpenUrl>,
}

enum PickerAction {
    Open,
    Save {
        params: Value,
        next: Option<Replace>,
    },
    Export {
        command: String,
        params: Value,
        extension: &'static str,
        alternatives: &'static [&'static str],
    },
    /// File › Run Script…: run the chosen JSON script.
    Script,
}

/// The status when the native file chooser fails inside the host ([`ModelApp::poll_picker`]).
const CHOOSER_FAILED: &str = "The file chooser could not open; try again";

pub(crate) struct PendingPicker {
    future: Pin<Box<dyn Future<Output = Option<String>>>>,
    action: PickerAction,
    revision: (u64, u64, Option<String>),
    document: u64,
}

const MODEL: Filter = Filter { name: "FieldCraft model", extensions: &["fieldcraft"] };
const DESIGN: Filter = Filter { name: "SolidCraft design", extensions: &["solidcraft"] };
const MESH: Filter = Filter { name: "Face-tagged mesh (3MF)", extensions: &["3mf", "3MF"] };
const ORPHAN: Filter = Filter { name: "Orphan mesh (Gmsh MSH 4.1)", extensions: &["msh", "MSH"] };
/// The Delete confirmation's icon: the kit's bin, as its delete buttons draw.
pub(crate) const DELETE_ICON: &str = "delete";
const SCRIPT: Filter = Filter { name: "FieldCraft script (JSON)", extensions: &["json", "JSON"] };
/// Largest script File › Run Script… reads (bytes), as `fieldcraft-cli run`.
const MAX_SCRIPT_BYTES: u64 = 256 * 1024 * 1024;
const ANY: Filter = Filter { name: "Models, designs and meshes", extensions: &["fieldcraft", "solidcraft", "3mf", "3MF", "msh", "MSH"] };

/// The Open picker's filters: everything FieldCraft opens, then each type.
pub const OPEN_FILTERS: &[Filter] = &[ANY, MODEL, DESIGN, MESH, ORPHAN];

/// What opening a file does with it, from its extension.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// A `.fieldcraft` model replaces the session (`file.open`).
    Model,
    /// A `.solidcraft` design becomes the geometry (`geometry.open`).
    Design,
    /// A `.3mf` face-tagged mesh is added as a surface (`geometry.import.mesh`).
    Mesh,
    /// A `.msh` orphan mesh replaces the mesh (`mesh.import.orphan`).
    Orphan,
}

impl Kind {
    /// The kind of a file by its extension (any case); `None` for a file FieldCraft does not open.
    pub fn of(path: &str) -> Option<Kind> {
        let ext = Path::new(path).extension()?.to_str()?.to_ascii_lowercase();
        match ext.as_str() {
            "fieldcraft" => Some(Kind::Model),
            "solidcraft" => Some(Kind::Design),
            "3mf" => Some(Kind::Mesh),
            "msh" => Some(Kind::Orphan),
            _ => None,
        }
    }

    pub fn command(self) -> &'static str {
        match self {
            Kind::Model => "file.open",
            Kind::Design => "geometry.open",
            Kind::Mesh => "geometry.import.mesh",
            Kind::Orphan => "mesh.import.orphan",
        }
    }

    /// The picker's filters when the user asked for this kind (its own type first).
    fn filters(self) -> &'static [Filter] {
        match self {
            Kind::Model => &[MODEL, ANY],
            Kind::Design => &[DESIGN, ANY],
            Kind::Mesh => &[MESH, ANY],
            Kind::Orphan => &[ORPHAN, ANY],
        }
    }

    /// What the window says while it opens and once it has.
    fn verbs(self) -> (&'static str, &'static str) {
        match self {
            Kind::Model | Kind::Design => ("Opening", "Opened"),
            Kind::Mesh | Kind::Orphan => ("Importing", "Imported"),
        }
    }
}

/// The file a command opens, if it is one the File menu handles with a picker.
fn kind_of_command(command: &str) -> Option<Kind> {
    [Kind::Model, Kind::Design, Kind::Mesh, Kind::Orphan].into_iter().find(|k| k.command() == command)
}

/// Whether `command` may run long: a solve, a mesher or an open (a sample meshes, a file can be
/// large), run as a job off the window thread.
pub fn is_long(command: &str) -> bool {
    BACKGROUND.contains(&command)
        || opens(command)
        || matches!(
            command,
            "solve.run" | "solve.static" | "solve.nonlinear" | "solve.heat" | "solve.frequency" | "solve.buckle" | "solve.explicit" | "solve.fluid"
        )
}

/// Whether `command` opens a model, a sample, a design or a mesh.
fn opens(command: &str) -> bool {
    fresh(command) || kind_of_command(command).is_some()
}

/// Whether `command` replaces the whole model: its job starts from an empty session rather than
/// a copy of the window's (copying a large model would stall the window).
fn fresh(command: &str) -> bool {
    matches!(command, "file.open" | "file.recover" | "file.sample" | "deck.read")
}

/// Commands a form runs in the background: they can take a while on a large part.
pub const BACKGROUND: &[&str] = &[
    "mesh.import.orphan",
    "deck.read",
    "mesh.tet.generate",
    "mesh.hex.generate",
    "mesh.hex.sweep",
    "mesh.hex.mapped",
    "mesh.hex.blocks",
    "mesh.hex.octree",
    "mesh.hex.project",
    "mesh.node.relax",
    "mesh.plan.apply",
    "mesh.cells",
    "geometry.partition",
    "results.path",
];

/// Most entries the recent-files list keeps, and the longest path it keeps.
pub const MAX_RECENT: usize = 10;
const MAX_RECENT_PATH: usize = 4096;
/// Largest recent-files file read.
const MAX_RECENT_BYTES: u64 = 256 * 1024;

/// Recently opened and saved files, newest first.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Recent {
    pub paths: Vec<String>,
    file: Option<PathBuf>,
}

impl Recent {
    /// The list kept in `file` (an unreadable or malformed file gives an empty list).
    pub fn load(file: Option<PathBuf>) -> Recent {
        let paths = file
            .as_ref()
            .filter(|f| std::fs::metadata(f).is_ok_and(|m| m.is_file() && m.len() <= MAX_RECENT_BYTES))
            .and_then(|f| std::fs::read(f).ok())
            .and_then(|b| serde_json::from_slice::<Vec<String>>(&b).ok())
            .unwrap_or_default()
            .into_iter()
            .filter(|p| !p.is_empty() && p.len() <= MAX_RECENT_PATH && Kind::of(p).is_some())
            .take(MAX_RECENT)
            .collect();
        Recent { paths, file }
    }

    /// Put `path` first and save the list.
    pub fn add(&mut self, path: &str) {
        if path.is_empty() || path.len() > MAX_RECENT_PATH {
            return;
        }
        self.paths.retain(|p| p != path);
        self.paths.insert(0, path.to_string());
        self.paths.truncate(MAX_RECENT);
        self.save();
    }

    fn save(&self) {
        let Some(file) = &self.file else { return };
        let write = || -> std::io::Result<()> {
            if let Some(dir) = file.parent() {
                std::fs::create_dir_all(dir)?;
            }
            std::fs::write(file, serde_json::to_vec_pretty(&self.paths).unwrap_or_default())
        };
        if let Err(e) = write() {
            log::warn!("recent files {}: {e}", file.display());
        }
    }
}

/// Largest preferences file read.
const MAX_PREFS_BYTES: u64 = 64 * 1024;

/// Where the preferences are kept: reads their JSON (`None` when there is none) and writes it.
pub struct PrefsStore {
    pub load: Box<dyn Fn() -> Option<Vec<u8>>>,
    pub save: Box<dyn Fn(&[u8]) -> Result<(), String>>,
}

impl PrefsStore {
    /// The preferences in `file` (an unreadable or oversized file reads as none).
    pub fn file(file: PathBuf) -> PrefsStore {
        let read = file.clone();
        PrefsStore {
            load: Box::new(move || {
                std::fs::metadata(&read).is_ok_and(|m| m.is_file() && m.len() <= MAX_PREFS_BYTES).then(|| std::fs::read(&read).ok()).flatten()
            }),
            save: Box::new(move |bytes| {
                if let Some(dir) = file.parent() {
                    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
                }
                std::fs::write(&file, bytes).map_err(|e| format!("{}: {e}", file.display()))
            }),
        }
    }
}

/// What the window remembers for the user between runs (`prefs.json` beside `recent.json`).
/// Keys it no longer knows (`show_setup`, from before the setup card went) are ignored.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Prefs {
    /// The mouse navigation preset ([`crate::navigation::NAV_PRESETS`] id).
    pub nav_preset: String,
    /// The wheel zooms the other way.
    pub reverse_zoom: bool,
    /// Share anonymous usage (docs/privacy.md); `None` until the user is asked.
    pub usage_sharing: Option<bool>,
    /// Count commands run by scripts and AI agents (MCP) too.
    pub usage_automation: bool,
    /// The random id usage is sent under, made when the user first agrees; empty before.
    pub install_id: String,
    /// Check for updates at start-up, at most once a day ([`crate::updates`]).
    pub update_check: bool,
    /// When the last automatic check was made (seconds since the Unix epoch; 0: never).
    pub update_checked_at: u64,
    /// The version the user skipped (the automatic check does not offer it); empty for none.
    pub update_skip: String,
    /// The units numbers are shown and typed in, for a model without its own (docs/units.md).
    pub units: fieldcraft_model::DisplaySystem,
}

impl Default for Prefs {
    fn default() -> Self {
        Prefs {
            nav_preset: crate::navigation::DEFAULT_PRESET.into(),
            reverse_zoom: false,
            usage_sharing: None,
            usage_automation: true,
            install_id: String::new(),
            update_check: true,
            update_checked_at: 0,
            update_skip: String::new(),
            units: fieldcraft_model::DisplaySystem::default(),
        }
    }
}

impl Prefs {
    /// The preferences kept in `store` (unreadable or malformed ones give the defaults).
    pub fn load(store: Option<&PrefsStore>) -> Prefs {
        store.and_then(|s| (s.load)()).and_then(|b| serde_json::from_slice(&b).ok()).unwrap_or_default()
    }

    /// Write them to `store`; a failure is logged (the preferences stay for this run).
    pub fn save(&self, store: Option<&PrefsStore>) {
        let Some(store) = store else { return };
        if let Err(e) = (store.save)(&serde_json::to_vec_pretty(self).unwrap_or_default()) {
            log::warn!("preferences: {e}");
        }
    }
}

/// A command running in the background on a copy of the session: a file opening, a solve, a
/// mesher started from its form.
pub struct Loading {
    /// What the window shows while it runs ("Opening beam.fieldcraft").
    pub title: String,
    /// The status when it succeeded ("Opened beam.fieldcraft").
    done: String,
    /// The file it opens, for the recent list.
    recent: Option<String>,
    /// The command form that started it: closed when it succeeds, showing the error when not.
    form: Option<&'static str>,
    /// The finished session copy and the command's reply.
    rx: Inbox<Finished>,
    /// The work, while it waits to run on the window thread (the web, which has no threads).
    pending: Option<Work>,
    /// Whether the window has drawn the job's progress: on the web the work runs (and the page
    /// stops drawing) only once the user has seen that it runs.
    shown: bool,
    /// The command it runs.
    pub command: &'static str,
    /// Where the command is, reported from the worker ([`fieldcraft_model::progress`]).
    pub progress: Arc<Progress>,
    /// When it started, on a monotonic clock: the window's frame time stands still while the
    /// window is hidden or covered (no frames are drawn), which made a 6 s solve "done in 0.0 s"
    /// (#590).
    pub started: web_time::Instant,
    /// Control-channel requests waiting for its reply (`engine.execute` of a long command,
    /// `job.wait`): each gets the command's result, or why it has none.
    pub waiters: Vec<Sender<Value>>,
}

impl Loading {
    /// The job for `ui.inspect` and `job.status`: what runs, where it is and for how long.
    pub fn to_json(&self) -> Value {
        let p = self.progress.snapshot();
        json!({
            "title": self.title,
            "command": self.command,
            "elapsed": self.started.elapsed().as_secs_f64(),
            "progress": p,
            "fraction": p.fraction(),
            "line": p.line(),
        })
    }

    /// The overlay's line: `Solving · Step 1/2 Load · Load factor 0.3 of 1 · 12 s`.
    pub fn line(&self) -> String {
        let p = self.progress.snapshot().line();
        let elapsed = elapsed(self.started.elapsed().as_secs_f64());
        if p.is_empty() { format!("{}… · {elapsed}", self.title) } else { format!("{}… · {p} · {elapsed}", self.title) }
    }
}

/// What the status bar says when job `command` finished: its `done` text with the time it took,
/// and for a mesher the elements made, their worst quality and, for the octree mesher, how deep
/// it refined and for which wall (lengths in `units`, the reply's).
fn finished(command: &str, done: &str, took: &str, reply: &Value, units: &fieldcraft_model::DisplaySystem) -> String {
    let r = &reply["result"];
    let bodies = || r["bodies"].as_array().into_iter().flatten();
    // The element kind: the reply's, else the one its bodies share.
    let kind = r["kind"].as_str().or_else(|| {
        let mut kinds = bodies().filter_map(|b| b["kind"].as_str());
        let first = kinds.next()?;
        kinds.all(|k| k == first).then_some(first)
    });
    match (r["elements"].as_u64(), kind) {
        (Some(n), Some(kind)) if BACKGROUND.contains(&command) => {
            let least = |f: &dyn Fn(&Value) -> Option<f64>| bodies().filter_map(f).reduce(f64::min);
            let mut q = match (least(&|b| b["quality"]["mesh"]["min_dihedral"].as_f64()), least(&|b| b["quality"]["min_scaled_jacobian"].as_f64())) {
                (Some(d), _) => format!(", smallest dihedral angle {d:.1}°"),
                (None, Some(j)) => format!(", smallest scaled Jacobian {}", fieldcraft_model::units::format_number(sig3(j))),
                _ => String::new(),
            };
            let levels = bodies().filter_map(|b| b["grid"]["levels"].as_u64()).max();
            if let Some(levels) = levels {
                let wall = least(&|b| b["grid"]["thinnest_wall"].as_f64())
                    .map(|w| {
                        format!(
                            " (thinnest wall {} {})",
                            fieldcraft_model::units::format_number(sig3(w)),
                            units.symbol(fieldcraft_model::units::Quantity::Length)
                        )
                    })
                    .unwrap_or_default();
                q.push_str(&format!(", {levels} octree levels{wall}"));
            }
            let n = fieldcraft_engine::thousands(usize::try_from(n).unwrap_or(usize::MAX));
            format!("{done}: {n} {kind}{q} in {took} — Mesh › Quality checks them")
        }
        _ if is_long(command) && !opens(command) => format!("{done} in {took}"),
        _ => done.to_string(),
    }
}

/// `v` to three significant digits.
fn sig3(v: f64) -> f64 {
    if v == 0.0 || !v.is_finite() {
        return v;
    }
    let scale = 10f64.powi(2 - v.abs().log10().floor() as i32);
    (v * scale).round() / scale
}

/// What the window says while a sample opens, and once it has.
fn sample_job(sample: &fieldcraft_engine::sample::Sample) -> (String, String) {
    (format!("Opening the {} sample", sample.label), format!("Opened the {} sample", sample.label))
}

/// A duration as `0.4 s`, `12 s` or `3 min 05 s`.
pub fn elapsed(seconds: f64) -> String {
    match seconds {
        s if s < 10.0 => format!("{s:.1} s"),
        s if s < 60.0 => format!("{s:.0} s"),
        s => format!("{} min {:02} s", (s / 60.0) as u64, (s % 60.0) as u64),
    }
}

/// A finished background job whose session has not replaced the window's yet.
pub struct Opened {
    /// The command the job ran.
    command: &'static str,
    session: Box<Session>,
    /// The session's display data, built with it.
    prebuilt: Box<Prebuilt>,
    reply: Value,
    done: String,
    recent: Option<String>,
    form: Option<&'static str>,
}

/// What waits for an answer about unsaved changes.
#[derive(Clone, Debug, PartialEq)]
pub enum Replace {
    /// A `.fieldcraft` file.
    File(String),
    /// An Abaqus model-data import with its explicit input-unit choice.
    Deck(Value),
    /// A sample (`file.sample`).
    Sample(&'static str),
    /// Restarting to finish installing an update ([`crate::updates`]).
    Update,
    /// An empty model (`file.new`).
    New,
}

impl Replace {
    fn name(&self) -> String {
        match self {
            Replace::File(p) => file_name(p),
            Replace::Deck(p) => p.get("path").and_then(Value::as_str).map(file_name).unwrap_or_else(|| "an Abaqus deck".into()),
            Replace::Sample(n) => fieldcraft_engine::sample::find(n).map_or_else(|| (*n).to_string(), |s| format!("the {} sample", s.label)),
            Replace::Update => "the update".into(),
            Replace::New => "a new model".into(),
        }
    }
}

/// A question the window asks before going on.
#[derive(Clone, Debug, PartialEq)]
pub enum Prompt {
    /// Opening a model or a sample discards the current one's unsaved changes.
    Unsaved(Replace),
    /// Autosaved models to recover (`file.recovery.list`).
    Recover,
    /// The samples to open (`file.sample`).
    Samples,
    /// A mesher's form would replace the mesh of bodies that have one: Replace runs it, Cancel
    /// leaves the form open.
    ReplaceMesh(fieldcraft_engine::Replacement),
    /// A design (its path) has none of the bodies the model's mesh and sets are on: start a new
    /// model with it, or keep the model. The opened session waits in [`ModelApp::held`].
    Design(String),
    /// Delete a definition (the item menu's Delete…, the Delete key): what uses it is listed, and
    /// the deletion refused or (a datum's partitions) done with it.
    Delete(crate::menu::DeleteAsk),
    /// The first run asks once whether to share anonymous usage ([`crate::privacy`]).
    Consent,
    /// Help › Settings: the privacy and update choices ([`crate::privacy`]).
    Settings,
    /// An update is downloaded and verified: restart to finish it, or later ([`crate::updates`]).
    Update,
    /// The last run crashed (its new crash reports, newest first): report it, or not now.
    Crashed(Vec<fieldcraft_engine::crash::Summary>),
    /// Help › Report a Bug ([`crate::report`]; the form is [`ModelApp::report`]).
    ReportBug,
    /// Help › Crash Reports: the saved reports, newest first, each to report, copy or delete
    /// ([`crate::report`]).
    CrashReports(Vec<fieldcraft_engine::crash::Summary>),
}

/// The answer to [`Prompt::Design`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DesignAnswer {
    /// Open the design into a new model (`geometry.open` with `"model": "new"`).
    New,
    /// Open the design and keep the model's mesh and definitions.
    Keep,
    /// Do not open the design.
    Cancel,
}

/// The answer to [`Prompt::Unsaved`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Answer {
    Save,
    Discard,
    Cancel,
}

/// The file name a model is saved under by default.
/// The document's name: its file's, else the sample it was opened as ("Cantilever (sample)").
fn model_name(session: &Session) -> String {
    match session.path.as_deref().and_then(|p| Path::new(p).file_stem()).and_then(|s| s.to_str()) {
        Some(stem) => stem.to_string(),
        None => session.sample.map_or_else(|| "Untitled".to_string(), |label| format!("{label} (sample)")),
    }
}

/// `path` with `ext` added when it has no extension of that name (some pickers leave it off).
fn with_extension(path: String, ext: &str) -> String {
    let has = Path::new(&path).extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case(ext));
    if has { path } else { format!("{path}.{ext}") }
}

/// Run `work` off the window thread (inline on the web, which has no threads).
pub(crate) fn off_thread(name: String, work: impl FnOnce() + Send + 'static) {
    if let Some(work) = start(name, Box::new(work)) {
        work();
    }
}

/// A background job's work.
type Work = Box<dyn FnOnce() + Send>;

/// What a background job sends back: the command's reply and what it made, or why it failed.
type Finished = Result<(Value, Product), String>;

/// What a successful background job hands the window.
enum Product {
    /// The session the job ran on, which replaces the window's, with its display data.
    Session(Box<Session>, Box<Prebuilt>),
    /// A read-only query ([`reads_only`]) leaves the model as it was: only the path plot and its
    /// search index come back, so the window keeps its session and its display.
    PathPlot(Option<Arc<fieldcraft_engine::path::Plot>>, fieldcraft_engine::path::Cache),
}

/// Whether a background `command` only reads the model: its job returns the view state it made
/// ([`Product::PathPlot`]) instead of a whole session.
fn reads_only(command: &str) -> bool {
    command == "results.path"
}

/// Start `work` on a worker thread. The web has no threads: the work comes back, to run on the
/// window thread once the window has shown that it runs ([`Loading::shown`]).
fn start(name: String, work: Work) -> Option<Work> {
    #[cfg(not(target_arch = "wasm32"))]
    {
        if let Err(e) = std::thread::Builder::new().name(name.clone()).spawn(work) {
            // The work went with the closure: a job then reports that it could not start.
            log::warn!("{name}: no worker thread: {e}");
        }
        None
    }
    #[cfg(target_arch = "wasm32")]
    {
        let _ = name;
        Some(work)
    }
}

/// Run `command` on `session` off the window thread: the reply's receiver, and the work when it
/// waits to run on the window thread (the web).
/// The display data is built with it, the design tessellated unless it is `drawn` (the window's).
fn spawn(
    session: Session,
    command: &'static str,
    params: Value,
    progress: Arc<Progress>,
    drawn: Option<Arc<fieldcraft_geometry::Design>>,
) -> (Inbox<Finished>, Option<Work>) {
    let (tx, rx) = crate::inbox::channel();
    let work = start(
        format!("fieldcraft {command}"),
        Box::new(move || {
            let mut s = session;
            #[allow(clippy::disallowed_methods, reason = "the job's worker thread")]
            let r = fieldcraft_model::progress::scoped(Some(progress), || s.execute(command, &params)).map_err(|e| e.to_string());
            let r = r.map(|v| {
                if reads_only(command) {
                    // The session copy is dropped here, on the worker, not in a frame.
                    (v, Product::PathPlot(s.path_plot.take(), std::mem::take(&mut s.path_index)))
                } else {
                    // What the window draws of the new session is built here too, not in a frame.
                    let prebuilt = Box::new(Prebuilt::of(&s, &drawn));
                    (v, Product::Session(Box::new(s), prebuilt))
                }
            });
            let _ = tx.send(r);
        }),
    );
    (rx, work)
}

impl ModelApp {
    fn has_open_picker(&self) -> bool {
        self.services.pick_file.is_some() || self.services.pick_open.is_some()
    }

    fn has_save_picker(&self) -> bool {
        self.services.pick_file.is_some() || self.services.pick_save.is_some()
    }

    /// Only user actions call this. A pending chooser never blocks commands or the event loop.
    fn pick_file(&mut self, request: FilePicker, action: PickerAction) -> bool {
        if self.picker.is_some() {
            self.set_status("Finish or cancel the open file chooser first", false);
            return false;
        }
        if self.ensure_idle().is_err() {
            return false;
        }
        if let Some(pick) = &self.services.pick_file {
            // A host that builds its native chooser here, and panics doing it, is caught as in
            // `poll_picker`.
            let Ok(future) = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| pick(request))) else {
                self.set_status(CHOOSER_FAILED, true);
                return false;
            };
            self.picker = Some(PendingPicker {
                future,
                action,
                revision: (self.revision, self.session.revision(), self.session.path.clone()),
                document: self.session.opened(),
            });
            return false;
        }
        let path = match request {
            FilePicker::Open { filters } => self.services.pick_open.as_ref().and_then(|p| p(&filters)),
            FilePicker::Save { name, filters } => self.services.pick_save.as_ref().and_then(|p| p(&name, &filters)),
        };
        path.is_some_and(|path| self.picked_file(action, path))
    }

    /// Poll only while a chooser is pending. No receive, wait or nested event loop runs here.
    pub(crate) fn poll_picker(&mut self) -> bool {
        let Some(mut pending) = self.picker.take() else { return false };
        // The native chooser runs in the host's code, outside any command's guard: AppKit can
        // return no save panel, which its binding turns into a panic (#780). That panic is caught
        // here and the chooser given up, so the window and the unsaved work stay; the crash hook
        // still records it. The broken future is dropped under the same guard.
        let polled = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| pending.future.as_mut().poll(&mut Context::from_waker(Waker::noop()))));
        let Ok(polled) = polled else {
            let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || drop(pending)));
            self.set_status(CHOOSER_FAILED, true);
            return false;
        };
        match polled {
            Poll::Pending => {
                self.picker = Some(pending);
                true
            }
            Poll::Ready(None) => false,
            Poll::Ready(Some(path)) => {
                // A control command, drop, new/open, save elsewhere or completed background job
                // may have changed what this choice was for while the native chooser was open.
                if self.loading.is_some()
                    || pending.document != self.session.opened()
                    || pending.revision != (self.revision, self.session.revision(), self.session.path.clone())
                {
                    self.set_status("The model changed while the file chooser was open; choose the file again", true);
                } else {
                    self.picked_file(pending.action, path);
                }
                false
            }
        }
    }

    fn picked_file(&mut self, action: PickerAction, path: String) -> bool {
        match action {
            PickerAction::Open => {
                self.open_path(&path);
                false
            }
            PickerAction::Script => self.run_script_file(&path).is_ok(),
            PickerAction::Save { mut params, next } => {
                params["path"] = json!(with_extension(path, "fieldcraft"));
                let saved = self.run("file.save_as", params).is_ok();
                self.saved_to_recent(saved);
                if saved && let Some(next) = next {
                    self.replace_after_save(next);
                }
                saved
            }
            PickerAction::Export { command, mut params, extension, alternatives } => {
                let path = if alternatives
                    .iter()
                    .any(|e| Path::new(&path).extension().and_then(|x| x.to_str()).is_some_and(|x| x.eq_ignore_ascii_case(e)))
                {
                    path
                } else {
                    with_extension(path, extension)
                };
                params["path"] = json!(path);
                if is_long(&command) {
                    let started = self.start_job(&command, params);
                    if let Err(e) = &started {
                        self.set_status(e.clone(), true);
                    }
                    return started.is_ok();
                }
                let written = self.run(&command, params).is_ok();
                if written && command != "deck.write" {
                    self.set_status(format!("Wrote {}", file_name(&path)), false);
                }
                written
            }
        }
    }

    /// The model's name: its file's stem, the sample it was opened as, or "Untitled".
    pub fn model_name(&self) -> String {
        model_name(&self.session)
    }

    /// What is running in the background ("Opening beam.fieldcraft", "Solving"), if anything.
    pub fn busy(&self) -> Option<&str> {
        self.loading.as_ref().map(|l| l.title.as_str())
    }

    /// `Err` (and the status says so) while a background job runs: commands wait for it.
    pub fn ensure_idle(&mut self) -> Result<(), String> {
        let Some(busy) = self.busy() else { return Ok(()) };
        let e = format!("{busy}: wait for it or cancel it");
        self.set_status(e.clone(), true);
        Err(e)
    }

    /// Run `command` in the background; the session it ran on replaces this one when it
    /// succeeds. A command that brings a whole model (`fresh`: opening a model, recovering one)
    /// starts from an empty session, the others from a copy of this one. One at a time.
    fn begin_job(&mut self, command: &'static str, mut params: Value, title: String, done: String, recent: Option<String>) -> bool {
        if self.ensure_idle().is_err() {
            return false;
        }
        // A fresh session is the app's own model too (its commands go into the crash journal),
        // with the same work caps and file access, as the engine keeps them when a command
        // replaces the session, and the window's count of camera requests.
        self.freeze_deck_units(command, &mut params);
        let base = if fresh(command) { self.session.empty_like() } else { self.session.clone() };
        let progress = Progress::new();
        let (rx, pending) = spawn(base, command, params, progress.clone(), self.bodies.0.clone());
        self.set_status(format!("{title}…"), false);
        let started = web_time::Instant::now();
        self.loading = Some(Loading { title, done, recent, form: None, rx, pending, shown: false, command, progress, started, waiters: Vec::new() });
        true
    }

    /// The form's OK: run its command, in the background for meshers and deck imports
    /// ([`BACKGROUND`]) so the window keeps drawing with Cancel; the form closes when the
    /// command succeeds and shows the error when it fails.
    /// A mesher that would replace a body's mesh asks first ([`Prompt::ReplaceMesh`]).
    pub fn submit_dialog(&mut self) {
        let Ok((command, mut params)) = self.dialog_params() else { return };
        if command == "deck.read" && self.has_unsaved_work() {
            self.freeze_deck_units(command, &mut params);
            self.prompt = Some(Prompt::Unsaved(Replace::Deck(params)));
            return;
        }
        if let Some(r) = fieldcraft_engine::replacement(&self.session, command, &params) {
            self.prompt = Some(Prompt::ReplaceMesh(r));
            return;
        }
        self.submit_now();
    }

    /// Answer [`Prompt::ReplaceMesh`]: run the form (replacing the mesh), or leave it open.
    pub fn answer_replace_mesh(&mut self, replace: bool) {
        if !matches!(self.prompt, Some(Prompt::ReplaceMesh(_))) {
            return;
        }
        self.prompt = None;
        if replace {
            self.submit_now();
        }
    }

    /// The form's OK once confirmed.
    fn submit_now(&mut self) {
        let Some(command) = self.dialog.as_ref().map(|d| d.command) else { return };
        if !BACKGROUND.contains(&command) {
            let _ = self.run_dialog();
            return;
        }
        let Ok((command, params)) = self.dialog_params() else { return };
        self.submit_background(command, params);
    }

    /// Preserve the input preset before a fresh job drops the current model. Display overrides
    /// are intentionally ignored: an Abaqus deck needs one consistent unit system.
    fn freeze_deck_units(&self, command: &str, params: &mut Value) {
        if command == "deck.read"
            && params.get("units").is_none_or(Value::is_null)
            && let Some(object) = params.as_object_mut()
        {
            object.insert("units".into(), json!(self.session.units().preset().label()));
        }
    }

    fn submit_background(&mut self, command: &'static str, params: Value) {
        let label = crate::command_label(command).to_string();
        // A file opened from its form says so as a dropped one does ("Opened box.fieldcraft"),
        // and is a recent file too (A4-6).
        let path = kind_of_command(command).and(params.get("path")).and_then(Value::as_str).map(str::to_string);
        let (title, done) = match (kind_of_command(command), &path) {
            (Some(kind), Some(path)) => {
                let (doing, done) = kind.verbs();
                (format!("{doing} {}", file_name(path)), format!("{done} {}", file_name(path)))
            }
            _ => (label.clone(), format!("{label} done")),
        };
        if self.begin_job(command, params, title, done, path)
            && let Some(l) = self.loading.as_mut()
        {
            l.form = Some(command);
        }
    }

    /// File > Open: pick a model, design or tagged mesh and open it. `None` asks for any of them.
    /// Without a picker service the command's form opens instead.
    pub fn open_dialog(&mut self, kind: Option<Kind>) {
        if !self.has_open_picker() {
            self.open_form(kind.unwrap_or(Kind::Model).command(), json!({}));
            return;
        }
        self.pick_file(FilePicker::Open { filters: kind.map_or(OPEN_FILTERS, Kind::filters).to_vec() }, PickerAction::Open);
    }

    /// File › Run Script…: pick a JSON script (an array of `{"command", "params"}`, or a journal)
    /// and run it, as `fieldcraft-cli run` does.
    pub fn run_script_dialog(&mut self) {
        if !self.has_open_picker() {
            self.set_status("Run Script needs a file chooser; scripts also run with fieldcraft-cli run", true);
            return;
        }
        self.pick_file(FilePicker::Open { filters: vec![SCRIPT] }, PickerAction::Script);
    }

    /// Read the script at `path` (bounded in size) and run it; the status says what happened.
    pub fn run_script_file(&mut self, path: &str) -> Result<(), String> {
        let name = file_name(path);
        let script = std::fs::metadata(path)
            .map_err(|e| e.to_string())
            .and_then(|m| if m.len() > MAX_SCRIPT_BYTES { Err(format!("larger than {} MB", MAX_SCRIPT_BYTES / (1024 * 1024))) } else { Ok(()) })
            .and_then(|()| std::fs::read_to_string(path).map_err(|e| e.to_string()))
            .and_then(|text| serde_json::from_str::<Value>(&text).map_err(|e| format!("not a JSON script ({e})")));
        match script {
            Ok(script) => self.run_script(&script),
            Err(e) => {
                let e = format!("{name}: {e}");
                self.set_status(e.clone(), true);
                Err(e)
            }
        }
    }

    /// Whether replacing the model would lose work: edits, or a still-unsaved imported deck
    /// whose history was intentionally cleared. A different document cannot inherit the marker.
    pub(crate) fn has_unsaved_work(&self) -> bool {
        self.session.is_dirty() && (!self.session.journal.is_empty() || self.imported_document == Some(self.session.opened()))
    }

    /// Open a file the user chose (picked, dropped or from the recent list) by its type. A model
    /// replacing unsaved changes asks first.
    pub fn open_path(&mut self, path: &str) {
        let Some(kind) = Kind::of(path) else {
            self.set_status(
                format!("{path}: FieldCraft opens .fieldcraft models, .solidcraft designs, .3mf tagged meshes and .msh orphan meshes"),
                true,
            );
            return;
        };
        if kind == Kind::Model && self.has_unsaved_work() {
            self.prompt = Some(Prompt::Unsaved(Replace::File(path.to_string())));
            return;
        }
        self.load(kind, path.to_string());
    }

    /// Open `path` in the background; files dropped while another opens wait their turn.
    fn load(&mut self, kind: Kind, path: String) {
        if self.loading.is_some() {
            self.queued.push_back(path);
            return;
        }
        let name = file_name(&path);
        let (doing, done) = kind.verbs();
        self.begin_job(kind.command(), json!({"path": path}), format!("{doing} {name}"), format!("{done} {name}"), Some(path));
    }

    /// Help › Samples (and the empty viewport): open a sample model, asking first when that
    /// would discard unsaved changes.
    pub fn open_sample(&mut self, name: &str) {
        let Some(sample) = fieldcraft_engine::sample::find(name) else {
            // An unknown name (a typo, an old `?sample=` link): the samples to choose from (#200).
            self.set_status(format!("There is no sample `{name}`: choose one of these"), true);
            self.prompt = Some(Prompt::Samples);
            return;
        };
        if self.has_unsaved_work() {
            self.prompt = Some(Prompt::Unsaved(Replace::Sample(sample.name)));
            return;
        }
        self.load_sample(sample.name);
    }

    /// File › New (`file.new`): start an empty model, asking first when that would discard
    /// unsaved changes.
    pub fn new_model(&mut self) {
        if self.has_unsaved_work() {
            self.prompt = Some(Prompt::Unsaved(Replace::New));
            return;
        }
        self.load_new();
    }

    fn load_new(&mut self) {
        if self.run("file.new", json!({})).is_ok() {
            self.model_changed();
            self.set_status("New model", false);
        }
    }

    /// Open a sample in the background: building its design and meshing it can take seconds.
    fn load_sample(&mut self, name: &str) {
        let Some(sample) = fieldcraft_engine::sample::find(name) else { return };
        let (title, done) = sample_job(sample);
        self.begin_job("file.sample", json!({"name": sample.name}), title, done, None);
    }

    /// Solve every step (`solve.run`) in the background; the results come up when it is done.
    pub fn solve_all(&mut self) {
        if let Some(u) = self.unavailable("solve.run") {
            self.refuse(u);
            return;
        }
        let n = self.session.model.steps().len();
        let done = format!("Solved {n} step{}", if n == 1 { "" } else { "s" });
        self.begin_job("solve.run", json!({}), "Solving".into(), done, None);
    }

    /// Run `command` as a background job for the control channel: the job's number (its
    /// opening), or why it cannot start (another job runs, the command cannot run yet).
    pub fn start_job(&mut self, command: &str, params: Value) -> Result<u64, String> {
        let spec = fieldcraft_engine::find_command(command).ok_or_else(|| format!("unknown command `{command}`"))?;
        if let Err(u) = self.session.availability(spec) {
            return Err(u.to_string());
        }
        let label = crate::command_label(spec.id).to_string();
        let title = if spec.id.starts_with("solve.") { "Solving".to_string() } else { label.clone() };
        let steps = self.session.model.steps().len();
        let sample = || params.get("name").and_then(Value::as_str).and_then(fieldcraft_engine::sample::find);
        let (title, done) = match spec.id {
            "file.sample" if let Some(sample) = sample() => sample_job(sample),
            "solve.run" => (title, format!("Solved {steps} step{}", if steps == 1 { "" } else { "s" })),
            _ => (title, format!("{label} done")),
        };
        if !self.begin_job(spec.id, params, title, done, None) {
            return Err(self.busy().map_or_else(|| "busy".to_string(), |b| format!("{b}: wait for it or cancel it")));
        }
        self.jobs = self.jobs.wrapping_add(1);
        Ok(self.jobs)
    }

    /// Recover an autosaved model (`file.recover`) in the background.
    pub fn recover(&mut self, id: &str) {
        self.prompt = None;
        self.begin_job("file.recover", json!({"id": id}), "Recovering the autosaved model".into(), "Recovered the autosaved model".into(), None);
    }

    /// Cancel the background job: its command stops at its next check
    /// ([`fieldcraft_model::progress::check`]) and what it made is dropped; the model stays as
    /// it was.
    pub fn cancel_loading(&mut self) {
        if let Some(l) = self.loading.take() {
            l.progress.cancel();
            self.queued.clear();
            for w in &l.waiters {
                let _ = w.send(json!({"ok": false, "error": format!("cancelled: {}", l.title), "cancelled": true}));
            }
            let line = if l.title == "Solving" { "Solve cancelled".to_string() } else { format!("{} cancelled", l.title) };
            self.set_status(line, false);
        }
    }

    /// Take in a finished background job (called every frame); `true` while one is running.
    pub fn poll_loading(&mut self) -> bool {
        let Some(l) = &mut self.loading else { return false };
        if l.pending.is_some() && !l.shown {
            return true;
        }
        if let Some(work) = l.pending.take() {
            work();
        }
        let outcome = match l.rx.try_recv() {
            Ok(done) => Some(done),
            Err(TryRecvError::Empty) => return true,
            Err(TryRecvError::Disconnected) => None,
        };
        let Some(l) = self.loading.take() else { return false };
        let reply = match &outcome {
            Some(Ok((v, _))) => json!({"ok": true, "result": v}),
            Some(Err(e)) => json!({"ok": false, "error": e}),
            None => json!({"ok": false, "error": format!("{}: stopped unexpectedly", l.title)}),
        };
        for w in &l.waiters {
            let _ = w.send(reply.clone());
        }
        let took = elapsed(l.started.elapsed().as_secs_f64());
        let l = Loading { done: finished(l.command, &l.done, &took, &reply, self.session.units()), ..l };
        match outcome {
            Some(Ok((reply, Product::PathPlot(plot, index)))) => {
                self.session.path_plot = plot;
                self.session.path_index = index;
                self.restore_path_plot();
                if self.dialog.as_ref().is_some_and(|d| l.form == Some(d.command)) {
                    self.form_done();
                }
                let done = reply.get("path").and_then(Value::as_str).map_or(l.done, |p| format!("Wrote {}", file_name(p)));
                self.set_status(done, false);
            }
            Some(Ok((reply, Product::Session(session, prebuilt)))) => {
                let opened = Opened { command: l.command, session, prebuilt, reply, done: l.done, recent: l.recent, form: l.form };
                // A design of another part over a meshed or defined model: ask what to do with
                // the model before it changes (the job ran with the default, keep).
                let unrelated = opened.reply.get("related").and_then(Value::as_bool) == Some(false)
                    && opened.reply.get("model").and_then(Value::as_str) == Some("kept");
                match opened.recent.clone() {
                    Some(path) if unrelated && Kind::of(&path) == Some(Kind::Design) => {
                        self.held = Some(opened);
                        self.set_status(format!("{} has none of the bodies the model is on", file_name(&path)), false);
                        self.prompt = Some(Prompt::Design(path));
                        return false;
                    }
                    _ => self.take_opened(opened),
                }
            }
            Some(Err(e)) => {
                let what = l.recent.as_deref().map_or_else(|| l.title.clone(), file_name);
                self.set_status(format!("{what}: {e}"), true);
                if let Some(d) = self.dialog.as_mut().filter(|d| l.form == Some(d.command)) {
                    d.error = Some(e);
                }
            }
            None => self.set_status(format!("{}: stopped unexpectedly", l.title), true),
        }
        self.open_queued();
        self.loading.is_some()
    }

    /// Open the next dropped file waiting its turn.
    fn open_queued(&mut self) {
        if let Some(next) = self.queued.pop_front() {
            self.open_path(&next);
        }
    }

    /// A finished job's session replaces the window's.
    fn take_opened(&mut self, o: Opened) {
        // The session it replaces may be large: it is freed off the window thread.
        let mut session = *o.session;
        // Its commands come from where the window's do now (the job may have started from the
        // command search or the control channel).
        session.origin = self.session.origin;
        let old = std::mem::replace(&mut self.session, session);
        off_thread("fieldcraft drop".into(), move || drop(old));
        self.model_changed();
        // A model without results has nothing for the Results workspace (every tool there would
        // be greyed): it comes up where its work is, the geometry and mesh or the definitions.
        if self.ui.workspace == Workspace::Results && self.session.result_steps().next().is_none() {
            self.ui.select_workspace(if self.session.model.elements().is_empty() { Workspace::Prepare } else { Workspace::Model });
        }
        self.prebuilt = Some(*o.prebuilt);
        if let Some(p) = &o.recent {
            self.recent.add(p);
        }
        if self.dialog.as_ref().is_some_and(|d| o.form == Some(d.command)) {
            self.form_done();
        }
        let warnings: Vec<&str> = o.reply.get("warnings").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str).collect();
        let mut done = o.done;
        if o.command == "deck.read" {
            self.imported_document = Some(self.session.opened());
            done = format!("Read Abaqus Deck: {} unread keywords (details in the panel)", o.reply["unknown_count"]);
            self.reply = Some(crate::reply::Reply {
                command: o.command.into(),
                label: crate::command_label(o.command).into(),
                params: json!({"format":o.reply["format"], "path":o.reply["path"], "units":o.reply["units"]}),
                value: o.reply.clone(),
            });
        }
        if o.command == "file.sample" && !self.session.model.steps().is_empty() {
            done.push_str(" — Solve shows its results");
        }
        if warnings.is_empty() {
            self.set_status(done, false);
        } else {
            self.set_status(format!("{done} — {}", warnings.join("; ")), true);
        }
    }

    /// Answer [`Prompt::Design`]: a new model with the design (opened again, in the background),
    /// the design over the kept model, or neither.
    pub fn answer_design(&mut self, answer: DesignAnswer) {
        let Some(Prompt::Design(path)) = self.prompt.clone() else { return };
        self.prompt = None;
        let Some(opened) = self.held.take() else { return };
        match answer {
            DesignAnswer::Keep => {
                let stale = opened.reply.get("stale").map(stale_note).unwrap_or_default();
                let done = format!("{}{stale}", opened.done);
                self.take_opened(Opened { done, ..opened });
            }
            DesignAnswer::New => {
                let name = file_name(&path);
                self.begin_job(
                    "geometry.open",
                    json!({"path": path, "model": "new"}),
                    format!("Opening {name}"),
                    format!("Opened {name} as a new model"),
                    Some(path),
                );
            }
            DesignAnswer::Cancel => self.set_status(format!("Did not open {}", file_name(&path)), false),
        }
        if self.loading.is_none() {
            self.open_queued();
        }
    }

    /// Block until the background job finishes (tests: a frame never waits for a job).
    #[cfg(test)]
    pub fn wait_loaded(&mut self) {
        for _ in 0..6000 {
            if !self.poll_loading() {
                return;
            }
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
    }

    /// File > Save: to the model's file, or Save As when it has none. Returns whether it saved.
    pub fn save(&mut self) -> bool {
        if self.session.path.is_none() {
            return self.save_as();
        }
        let saved = self.run("file.save", json!({})).is_ok();
        self.saved_to_recent(saved)
    }

    /// After a save: the model's file goes first in the recent list.
    fn saved_to_recent(&mut self, saved: bool) -> bool {
        if saved && let Some(p) = self.session.path.clone() {
            self.recent.add(&p);
        }
        saved
    }

    /// File > Save As: pick a new file for the model. Without a picker the command's form opens.
    pub fn save_as(&mut self) -> bool {
        self.save_as_with(json!({}))
    }

    /// File > Save With Geometry Embedded: Save As, keeping a copy of the design and meshes in
    /// the file (`file.save_as`'s `geometry: "embed"`).
    pub fn save_packed(&mut self) -> bool {
        self.save_as_with(json!({"geometry": "embed"}))
    }

    /// Save As with `params` besides the path.
    fn save_as_with(&mut self, params: Value) -> bool {
        self.save_as_then(params, None)
    }

    fn save_as_then(&mut self, params: Value, next: Option<Replace>) -> bool {
        if !self.has_save_picker() {
            self.open_form("file.save_as", params);
            return false;
        }
        self.pick_file(
            FilePicker::Save { name: format!("{}.fieldcraft", self.model_name()), filters: vec![MODEL] },
            PickerAction::Save { params, next },
        )
    }

    /// Write Deck: pick where to write the model as a `format` deck. Without a picker the
    /// command's form opens with the format filled in.
    pub fn write_deck(&mut self, format: &str) {
        let Some(f) = fieldcraft_decks::format(format) else {
            self.set_status(format!("no {format} deck writer yet"), true);
            return;
        };
        if !self.has_save_picker() {
            self.open_form("deck.write", json!({"format": format}));
            return;
        }
        let exts: &'static [&'static str] = match f.extension {
            "inp" => &["inp"],
            "bdf" => &["bdf", "dat", "nas"],
            "k" => &["k", "key", "dyn"],
            _ => &[],
        };
        let filter = Filter { name: "Solver deck", extensions: exts };
        self.pick_file(
            FilePicker::Save { name: format!("{}.{}", self.model_name(), f.extension), filters: vec![filter] },
            PickerAction::Export { command: "deck.write".into(), params: json!({"format": format}), extension: f.extension, alternatives: exts },
        );
    }

    /// Results › Export CSV… (`results.export.csv` with `params`, its `what` and `step` among
    /// them): asks the save picker for `<model>-<step>-<what>.csv`; without a picker, the
    /// command's form opens instead.
    pub fn export_csv(&mut self, params: Value) {
        self.export_csv_with("results.export.csv", params);
    }

    /// [`Self::export_csv`] through `command`, which takes a `path` (results.export.csv,
    /// results.xy).
    pub fn export_csv_with(&mut self, command: &str, params: Value) {
        if !self.has_save_picker() {
            self.open_form(command, params);
            return;
        }
        let part = |k: &str| params.get(k).and_then(Value::as_str).unwrap_or_default().to_string();
        let what = match command {
            "results.xy" => "xy".to_string(),
            "results.path" => "path".to_string(),
            _ => part("what"),
        };
        let name = format!("{}-{}-{what}.csv", self.model_name(), part("step"));
        self.pick_file(
            FilePicker::Save { name, filters: vec![Filter { name: "CSV", extensions: &["csv"] }] },
            PickerAction::Export { command: command.into(), params, extension: "csv", alternatives: &["csv"] },
        );
    }

    /// The autosaved models that can be recovered (`file.recovery.list`), newest first.
    pub fn recoverable(&mut self) -> Vec<Value> {
        #[allow(clippy::disallowed_methods, reason = "lists the autosave folder: a few small files")]
        let listed = self.session.execute("file.recovery.list", &json!({}));
        listed.ok().and_then(|v| v.as_array().cloned()).unwrap_or_default()
    }

    /// Answer the unsaved-changes prompt.
    pub fn answer_unsaved(&mut self, answer: Answer) {
        let Some(Prompt::Unsaved(next)) = self.prompt.clone() else { return };
        self.prompt = None;
        if answer == Answer::Save && self.session.path.is_none() {
            self.save_as_then(json!({}), Some(next));
            return;
        }
        let go = match answer {
            Answer::Save => self.save(),
            Answer::Discard => true,
            Answer::Cancel => false,
        };
        if !go {
            return;
        }
        self.replace_after_save(next);
    }

    fn replace_after_save(&mut self, next: Replace) {
        match next {
            Replace::File(path) => {
                if let Some(kind) = Kind::of(&path) {
                    self.load(kind, path);
                }
            }
            Replace::Deck(params) => self.submit_background("deck.read", params),
            Replace::Sample(name) => self.load_sample(name),
            Replace::New => self.load_new(),
            Replace::Update => self.close_to_update(),
        }
    }

    /// What a click on a file entry does, when it is one: open, save, write a deck or recover
    /// (through the host's pickers when `pickers`), choose a sample, or solve in the background.
    /// `false` for every other command, and for picker entries without pickers (their form
    /// opens instead).
    pub(crate) fn file_action(&mut self, command: &str, preset: &Value, pickers: bool) -> bool {
        let can_save = pickers && self.has_save_picker();
        match command {
            "file.save" if self.session.path.is_some() || can_save => {
                self.save();
            }
            "file.save_as" if can_save => {
                self.save_as();
            }
            "file.recover" | "file.recovery.list" => self.prompt = Some(Prompt::Recover),
            "file.sample" => self.prompt = Some(Prompt::Samples),
            "file.new" => self.new_model(),
            "diagnostics.bug_report" => self.report_bug(None),
            "diagnostics.crash_reports" | "diagnostics.crash_report" | "diagnostics.crash_reports.delete" => self.open_crash_reports(),
            "solve.run" => self.solve_all(),
            "results.xy" => self.open_xy(),
            "results.path" => self.open_path_plot(),
            "deck.write" if can_save => match preset.get("format").and_then(Value::as_str) {
                Some(f) => self.write_deck(f),
                None => return false,
            },
            other => match kind_of_command(other) {
                Some(kind) if pickers && self.has_open_picker() => self.open_dialog(Some(kind)),
                _ => return false,
            },
        }
        true
    }

    /// Files dropped on the window (and by the control channel's `ui.drop`) open in turn; a
    /// hovered drop shows where it goes.
    pub(crate) fn dropped_files(&mut self, ctx: &egui::Context) {
        let mut dropped: Vec<String> = ctx
            .input(|i| i.raw.dropped_files.iter().filter_map(|f| Some(f.path().to_string_lossy().to_string()).filter(|p| !p.is_empty())).collect());
        dropped.append(&mut self.synthetic_drops);
        // A model first: it replaces the session, then designs and meshes go into it.
        let (models, rest): (Vec<String>, Vec<String>) = dropped.into_iter().partition(|p| Kind::of(p) == Some(Kind::Model));
        for p in models.into_iter().take(1).chain(rest) {
            self.open_path(&p);
        }
    }
}

/// What a kept model has that the design does not (`geometry.open`'s `stale`), for the status:
/// empty when everything resolves.
fn stale_note(stale: &Value) -> String {
    let names = |k: &str| stale.get(k).and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str).collect::<Vec<_>>().join(", ");
    let parts: Vec<String> = [("mesh of", names("mesh_bodies")), ("sets", names("sets")), ("seeds on", names("seeds"))]
        .into_iter()
        .filter(|(_, n)| !n.is_empty())
        .map(|(what, n)| format!("{what} {n}"))
        .collect();
    if parts.is_empty() { String::new() } else { format!(" — not on this design: {}", parts.join("; ")) }
}

/// The last part of a path, for messages.
pub fn file_name(path: &str) -> String {
    Path::new(path).file_name().and_then(|n| n.to_str()).map_or_else(|| path.to_string(), str::to_string)
}

/// How long ago, from an age in seconds: "just now", "5 min ago", "3 h ago", "2 days ago".
pub fn ago(seconds: u64) -> String {
    match seconds {
        0..60 => "just now".into(),
        60..3600 => format!("{} min ago", seconds / 60),
        3600..86_400 => format!("{} h ago", seconds / 3600),
        86_400..172_800 => "1 day ago".into(),
        _ => format!("{} days ago", seconds / 86_400),
    }
}

/// A confirmation's or window's message, wrapping in its column.
pub fn message(ui: &mut egui::Ui, text: &str) {
    ui.add(egui::Label::new(egui::RichText::new(text).color(solidcraft_theme::Tokens::get().form_ink)).wrap());
}

/// A dim note under a message, wrapping in its column, at the kit's label size.
pub fn note(ui: &mut egui::Ui, text: &str) {
    kit::status_line(ui, text, kit::StatusTone::Dim);
}

/// Whether a window or confirmation was dismissed: its close button, or its cancel action (Esc).
pub fn dismissed(closed: bool, action: kit::Action) -> bool {
    closed || action == kit::Action::Cancel
}

/// What a confirmation answered: `ok`, `left` (its left action) or, dismissed, `cancel`.
fn answer<T>(closed: bool, action: kit::Action, ok: T, left: T, cancel: T) -> Option<T> {
    match action {
        kit::Action::Ok => Some(ok),
        kit::Action::Apply => Some(left),
        _ if dismissed(closed, action) => Some(cancel),
        _ => None,
    }
}

/// The prompts (unsaved changes, recover), the loading overlay and the drop hint.
pub fn show(app: &mut ModelApp, ctx: &egui::Context) {
    let t = crate::theme::Tokens::get(app.ui.dark);
    let hovering = ctx.input(|i| !i.raw.hovered_files.is_empty());
    if hovering || app.loading.is_some() {
        let fraction = app.loading.as_ref().and_then(|l| l.progress.snapshot().fraction());
        let text = match app.loading.as_ref() {
            Some(l) => l.line(),
            None => "Drop to open: .fieldcraft models, .solidcraft designs, .3mf tagged meshes, .msh orphan meshes".into(),
        };
        let mut cancel = false;
        // A notification surface (the kit's Toast), independent of the workspace layout, with the
        // kit's progress meter while the fraction is known.
        let k = solidcraft_theme::Tokens::get();
        let host = kit::Container::new(ctx, kit::Surface::Toast);
        let loading = app.loading.is_some();
        host.area(egui::Id::unique("fc_loading")).order(egui::Order::Foreground).anchor(egui::Align2::CENTER_CENTER, egui::Vec2::ZERO).show(
            ctx,
            |ui| {
                let frame =
                    host.frame(ctx, k.panel, kit::frame::Side::Left).inner_margin(egui::Margin::symmetric(k.space_16 as i8, k.space_12 as i8));
                host.show_frame(ui, frame, ctx.content_rect(), |ui| match fraction.filter(|_| loading) {
                    Some(f) => cancel = kit::progress_action(ui, "fc_loading", &text, "", f as f32, "Cancel").clicked(),
                    None => {
                        ui.horizontal(|ui| {
                            if loading {
                                ui.add(egui::Spinner::new().size(k.space_14).color(k.text));
                            }
                            kit::text(ui, &text, kit::Text::Inline);
                            if loading {
                                cancel = kit::button(ui, "Cancel", kit::Button::Normal).clicked();
                            }
                        });
                    }
                });
            },
        );
        if let Some(l) = app.loading.as_mut() {
            l.shown = true;
        }
        if cancel {
            app.cancel_loading();
        }
        if app.loading.is_some() {
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
    }
    match app.prompt.clone() {
        Some(Prompt::Unsaved(next)) => {
            let question = match next {
                Replace::Update => format!("Save the changes to {} before restarting to update?", app.model_name()),
                _ => format!("Save the changes to {} before opening {}?", app.model_name(), next.name()),
            };
            let (closed, action) = kit::Confirm::new(egui::Id::unique("fc_unsaved"), "file.save", "Unsaved changes")
                .tone(kit::Tone::Warning)
                .actions("Save", "Cancel", Some("Don't Save"))
                .secondary_left()
                .show(ctx, true, |ui| {
                    message(ui, &question);
                    note(ui, "Without saving, they are lost.");
                });
            if let Some(a) = answer(closed, action, Answer::Save, Answer::Discard, Answer::Cancel) {
                app.answer_unsaved(a);
            }
        }
        Some(Prompt::Design(path)) => {
            let question = format!("{} has none of the bodies {}'s mesh and sets are on.", file_name(&path), app.model_name());
            let (closed, action) = kit::Confirm::new(egui::Id::unique("fc_design"), "geometry.open", "Open design")
                .actions("Start a New Model", "Cancel", Some("Keep the Current Model"))
                .secondary_left()
                .show(ctx, true, |ui| {
                    message(ui, &question);
                    note(ui, "A new model starts empty with this design; keeping the model keeps its mesh and definitions, which then refer to bodies the design does not have.");
                });
            if let Some(a) = answer(closed, action, DesignAnswer::New, DesignAnswer::Keep, DesignAnswer::Cancel) {
                app.answer_design(a);
            }
        }
        Some(Prompt::Delete(ask)) => {
            let allowed = ask.allowed();
            let (closed, action) = kit::Confirm::new(egui::Id::unique("fc_delete"), DELETE_ICON, "Delete")
                .tone(kit::Tone::Danger)
                .actions(if allowed { "Delete" } else { "" }, if allowed { "Cancel" } else { "Close" }, None)
                .show(ctx, allowed, |ui| {
                    message(ui, &ask.question());
                    if allowed && !ask.users.is_empty() {
                        note(ui, "They are deleted too; Undo brings everything back.");
                    }
                });
            if let Some(delete) = answer(closed, action, allowed, false, false) {
                app.answer_delete(delete);
            }
        }
        Some(Prompt::ReplaceMesh(r)) => {
            let (closed, action) = kit::Confirm::new(egui::Id::unique("fc_replace_mesh"), "mesh.tet.generate", "Replace mesh")
                .tone(kit::Tone::Warning)
                .actions("Replace", "Cancel", None)
                .show(ctx, true, |ui| {
                    message(ui, &r.question);
                    for w in &r.warnings {
                        note(ui, w);
                    }
                });
            if let Some(replace) = answer(closed, action, true, false, false) {
                app.answer_replace_mesh(replace);
            }
        }
        Some(Prompt::Recover) => {
            let entries = app.recover_list(ctx);
            let (mut recover, mut discard) = (None, None);
            let (closed, action) = crate::frame::window(egui::Id::unique("fc_recover"), "file.recover", "Recover", crate::frame::Width::Wide)
                .actions("", "Close", (entries.len() > 1).then_some("Discard All"))
                .show(ctx, true, |ui| {
                    if entries.is_empty() {
                        note(ui, "No autosaved models to recover.");
                    }
                    // Names, ages and sizes whole: the kit's window body truncates text by default.
                    ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
                    egui::Grid::new("fc_recover_grid").num_columns(4).spacing(egui::vec2(10.0, 4.0)).show(ui, |ui| {
                        for e in &entries {
                            let id = e.get("id").and_then(Value::as_str).unwrap_or_default().to_string();
                            let path = e.get("path").and_then(Value::as_str);
                            let n = |k: &str| e.get(k).and_then(Value::as_u64).unwrap_or(0);
                            let name = ui.label(path.map_or_else(|| "Untitled".to_string(), file_name));
                            if let Some(path) = path {
                                name.backing_on_hover_text(path);
                            }
                            ui.label(egui::RichText::new(ago(n("age"))).color(t.text_dim));
                            ui.label(egui::RichText::new(format!("{} nodes, {} elements", n("nodes"), n("elements"))).color(t.text_dim));
                            ui.horizontal(|ui| {
                                if kit::button(ui, "Recover", kit::Button::Small).clicked() {
                                    recover = Some(id.clone());
                                }
                                if kit::button(ui, "Discard", kit::Button::Small).clicked() {
                                    discard = Some(id.clone());
                                }
                            });
                            ui.end_row();
                        }
                    });
                });
            let discard_all = action == kit::Action::Apply;
            if let Some(id) = recover {
                app.recover(&id);
            } else if discard.is_some() || discard_all {
                let params = discard.map_or_else(|| json!({"all": true}), |id| json!({"id": id}));
                let _ = app.run("file.recovery.discard", params);
                ctx.data_mut(|d| d.remove::<Vec<Value>>(egui::Id::unique("fc_recover_entries")));
            } else if dismissed(closed, action) {
                app.prompt = None;
                ctx.data_mut(|d| d.remove::<Vec<Value>>(egui::Id::unique("fc_recover_entries")));
            }
        }
        Some(Prompt::Samples) => {
            let mut open = None;
            let (closed, action) = crate::frame::window(egui::Id::unique("fc_samples"), "file.sample", "Samples", crate::frame::Width::Wide)
                .actions("", "Close", None)
                .show(ctx, true, |ui| {
                    crate::samples::grid(ui, Some((ctx.content_rect().height() * 0.6).max(160.0)), &mut open);
                });
            if let Some(name) = open {
                app.prompt = None;
                app.open_sample(name);
            } else if dismissed(closed, action) {
                app.prompt = None;
            }
        }
        // Drawn by `crate::report::show` and `crate::privacy::show`.
        Some(Prompt::Crashed(_) | Prompt::CrashReports(_) | Prompt::ReportBug | Prompt::Consent | Prompt::Settings | Prompt::Update) | None => {}
    }
}

impl ModelApp {
    /// The recover window's entries, read from disk once while it is open.
    fn recover_list(&mut self, ctx: &egui::Context) -> Vec<Value> {
        let key = egui::Id::unique("fc_recover_entries");
        if let Some(v) = ctx.data(|d| d.get_temp::<Vec<Value>>(key)) {
            return v;
        }
        let v = self.recoverable();
        ctx.data_mut(|d| d.insert_temp(key, v.clone()));
        v
    }
}

/// The File menu's rows: open, recent files, save, decks, samples, recover. The title bar's
/// File menu.
pub fn file_menu(app: &mut ModelApp, ui: &mut egui::Ui) {
    let key = |k: &str| crate::appbar::shortcut(k);
    let row = |ui: &mut egui::Ui, id: &str, label: &str, right: Option<&str>| crate::frame::menu_row(ui, id, label, right, true).clicked();
    let mut act: Option<Box<dyn FnOnce(&mut ModelApp)>> = None;
    if row(ui, "file.new", "New", Some(&key("N"))) {
        act = Some(Box::new(|app| app.new_model()));
    }
    if row(ui, "file.open", "Open…", Some(&key("O"))) {
        act = Some(Box::new(|app| app.open_dialog(None)));
    }
    if row(ui, "file.sample", "Samples…", None) {
        act = Some(Box::new(|app| app.prompt = Some(Prompt::Samples)));
    }
    if !app.recent.paths.is_empty() {
        crate::frame::menu_heading(ui, "Recent");
        for (i, p) in app.recent.paths.iter().take(6).enumerate() {
            let name = file_name(p);
            if crate::frame::menu_row(ui, &format!("file.recent.{i}"), &name, None, true).backing_on_hover_text(p).clicked() {
                let p = p.clone();
                act = Some(Box::new(move |app| app.open_path(&p)));
            }
        }
    }
    crate::frame::menu_separator(ui);
    if row(ui, "file.save", "Save", Some(&key("S"))) {
        act = Some(Box::new(|app| {
            app.save();
        }));
    }
    if row(ui, "file.save_as", "Save As…", Some(&key("Shift+S"))) {
        act = Some(Box::new(|app| {
            app.save_as();
        }));
    }
    let packed = crate::frame::menu_row(
        ui,
        "file.save_packed",
        "Save With Geometry Embedded…",
        None,
        app.session.design.is_some() || !app.session.surfaces.is_empty(),
    )
    .backing_on_hover_text(
        "Save As, with a copy of the design and imported meshes inside the file, so it opens anywhere: whoever you send it to gets your geometry too",
    );
    if packed.clicked() {
        act = Some(Box::new(|app| {
            app.save_packed();
        }));
    }
    crate::frame::menu_separator(ui);
    if row(ui, "deck.read", "Read Deck…", None) {
        act = Some(Box::new(|app| {
            let _ = app.start("deck.read");
        }));
    }
    if row(ui, "file.script", "Run Script…", None) {
        act = Some(Box::new(|app| app.run_script_dialog()));
    }
    // Each group behind a rule, as the others.
    crate::frame::menu_separator(ui);
    crate::frame::menu_heading(ui, "Write Deck");
    for c in fieldcraft_engine::catalogue::CATALOGUE.iter().filter(|c| c.panel == "Decks" && c.id != "deck.read") {
        let format = c.preset().get("format").and_then(Value::as_str).map(str::to_string);
        let soon = c.milestone().map(|m| format!("soon · {m}"));
        let enabled = format.is_some() && soon.is_none();
        if crate::frame::menu_row(ui, c.id, &format!("{}…", c.label), soon.as_deref(), enabled).backing_on_hover_text(c.tip).clicked()
            && let Some(f) = format
        {
            act = Some(Box::new(move |app| app.write_deck(&f)));
        }
    }
    crate::frame::menu_separator(ui);
    if crate::frame::menu_row(ui, "file.recover", "Recover…", None, true)
        .backing_on_hover_text("Models autosaved by an app that crashed or closed with unsaved changes")
        .clicked()
    {
        act = Some(Box::new(|app| app.prompt = Some(Prompt::Recover)));
    }
    if let Some(f) = act {
        app.ui.menu = None;
        ui.close();
        f(app);
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::cell::RefCell;
    use std::rc::Rc;

    use super::*;

    fn temp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("fieldcraft-files-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[derive(Default)]
    struct Chooser {
        answer: RefCell<Option<Option<String>>>,
        requests: RefCell<Vec<FilePicker>>,
    }

    fn deferred(app: &mut ModelApp) -> Rc<Chooser> {
        let chooser = Rc::new(Chooser::default());
        let host = chooser.clone();
        app.services.pick_file = Some(Box::new(move |request| {
            host.requests.borrow_mut().push(request);
            let host = host.clone();
            Box::pin(std::future::poll_fn(move |_| match host.answer.borrow_mut().take() {
                Some(answer) => Poll::Ready(answer),
                None => Poll::Pending,
            }))
        }));
        chooser
    }

    /// A native chooser that panics, when it is built or while it runs (AppKit returning no save
    /// panel, #780), is given up with a status; the app keeps its model and can choose again.
    #[test]
    fn a_chooser_that_panics_is_given_up_without_taking_the_app_down() {
        let mut app = ModelApp::new(Session::default());
        app.services.pick_file =
            Some(Box::new(|_| Box::pin(std::future::poll_fn(|_| -> Poll<Option<String>> { panic!("unexpected NULL from the save panel") }))));
        app.save_as();
        assert!(app.picker.is_some());
        assert!(!app.poll_picker());
        assert!(app.picker.is_none() && app.session.path.is_none());
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text == CHOOSER_FAILED));
        app.services.pick_file = Some(Box::new(|_| panic!("the panel could not be built")));
        app.save_as();
        assert!(app.picker.is_none());
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text == CHOOSER_FAILED));
        // A working chooser afterwards opens as usual.
        let chooser = deferred(&mut app);
        app.save_as();
        assert_eq!(chooser.requests.borrow().len(), 1);
    }

    #[test]
    fn a_pending_chooser_keeps_frames_running_and_only_opens_once() {
        let mut app = ModelApp::new(Session::default());
        let chooser = deferred(&mut app);
        app.open_dialog(None);
        app.save_as();
        app.open_dialog(Some(Kind::Design));
        assert_eq!(chooser.requests.borrow().len(), 1);
        assert_eq!(chooser.requests.borrow()[0], FilePicker::Open { filters: OPEN_FILTERS.to_vec() });
        let ctx = egui::Context::default();
        for _ in 0..8 {
            let mut output = ctx.run_ui(Default::default(), |ui| app.logic(ui.ctx()));
            output.textures_delta.clear();
            assert!(app.picker.is_some() && app.loading.is_none());
        }
        *chooser.answer.borrow_mut() = Some(None);
        assert!(!app.poll_picker());
        assert!(app.picker.is_none() && app.session.path.is_none());
        app.save_as();
        assert_eq!(chooser.requests.borrow().len(), 2, "cancelling releases the chooser guard");
    }

    #[test]
    fn a_chooser_started_by_layout_schedules_its_first_poll_without_more_input() {
        let mut app = ModelApp::new(Session::default());
        let _chooser = deferred(&mut app);
        let ctx = egui::Context::default();
        let screen_rect = Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0)));
        // Settle fonts/layout/animations so an initial-frame repaint cannot hide the bug.
        for frame in 0..10 {
            let input = egui::RawInput { time: Some(f64::from(frame)), screen_rect, ..Default::default() };
            let mut output = ctx.run_ui(input, |ui| {
                app.logic(ui.ctx());
                app.ui(ui);
            });
            output.textures_delta.clear();
        }
        let input = egui::RawInput { time: Some(10.0), screen_rect, ..Default::default() };
        let mut output = ctx.run_ui(input, |ui| {
            app.logic(ui.ctx());
            // A menu/button action arrives during layout, after logic has already polled.
            app.open_dialog(None);
            app.ui(ui);
        });
        let repaint = output.viewport_output.get(&egui::ViewportId::ROOT).map(|v| v.repaint_delay);
        output.textures_delta.clear();
        assert!(app.picker.is_some());
        assert!(repaint.is_some_and(|delay| delay <= std::time::Duration::from_millis(50)));
    }

    #[test]
    fn asynchronous_save_then_open_and_import_use_the_existing_paths() {
        let dir = temp("async-save-open");
        let path = dir.join("kept.fieldcraft").to_string_lossy().into_owned();
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let chooser = deferred(&mut app);
        assert!(!app.save_packed(), "a pending save has not saved yet");
        assert!(!Path::new(&path).exists());
        *chooser.answer.borrow_mut() = Some(Some(path.clone()));
        assert!(!app.poll_picker());
        assert_eq!(app.session.path.as_deref(), Some(path.as_str()));
        assert!(app.session.pack_geometry && !app.session.is_dirty());
        assert_eq!(app.recent.paths.as_slice(), std::slice::from_ref(&path));
        app.run("file.new", json!({})).unwrap();
        app.open_dialog(Some(Kind::Model));
        *chooser.answer.borrow_mut() = Some(Some(path));
        app.poll_picker();
        app.wait_loaded();
        assert!(app.session.model.node(1).is_some());
        app.open_dialog(Some(Kind::Mesh));
        assert_eq!(chooser.requests.borrow().last(), Some(&FilePicker::Open { filters: Kind::Mesh.filters().to_vec() }));
        // A bad import is still handled by the normal error path without replacing the model.
        *chooser.answer.borrow_mut() = Some(Some(dir.join("missing.3mf").to_string_lossy().into_owned()));
        app.poll_picker();
        app.wait_loaded();
        assert!(app.session.model.node(1).is_some());
        assert!(app.status.as_ref().is_some_and(|s| s.error));
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn chooser_results_for_a_changed_or_replaced_model_are_discarded() {
        let dir = temp("async-stale");
        for command in ["mesh.node.add", "file.new"] {
            let target = dir.join(format!("{command}.fieldcraft"));
            let mut app = ModelApp::new(Session::default());
            let chooser = deferred(&mut app);
            app.save_as();
            let params = if command == "file.new" { json!({}) } else { json!({"label": 42, "xyz": [0, 0, 0]}) };
            app.run(command, params).unwrap();
            *chooser.answer.borrow_mut() = Some(Some(target.to_string_lossy().into_owned()));
            app.poll_picker();
            assert!(!target.exists());
            assert!(app.session.path.is_none());
            assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.contains("model changed")));
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn a_chooser_cannot_apply_a_selection_while_a_background_job_runs() {
        let dir = temp("async-background");
        let target = dir.join("stale.fieldcraft");
        let mut app = ModelApp::new(Session::default());
        let chooser = deferred(&mut app);
        app.save_as();
        app.start_job("file.new", json!({})).unwrap();
        *chooser.answer.borrow_mut() = Some(Some(target.to_string_lossy().into_owned()));
        app.poll_picker();
        assert!(!target.exists() && app.session.path.is_none());
        app.wait_loaded();
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn save_before_replacing_waits_for_success_and_never_replaces_after_cancel_or_failure() {
        let dir = temp("async-before-new");
        for outcome in ["success", "cancel", "failure"] {
            let mut app = ModelApp::new(Session::default());
            app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
            let chooser = deferred(&mut app);
            app.new_model();
            app.answer_unsaved(Answer::Save);
            assert!(app.session.model.node(1).is_some() && app.picker.is_some());
            let answer = match outcome {
                "cancel" => None,
                "failure" => {
                    let blocker = dir.join("not-a-directory");
                    std::fs::write(&blocker, "blocked").unwrap();
                    Some(blocker.join("model.fieldcraft").to_string_lossy().into_owned())
                }
                _ => Some(dir.join("kept.fieldcraft").to_string_lossy().into_owned()),
            };
            *chooser.answer.borrow_mut() = Some(answer);
            app.poll_picker();
            assert_eq!(app.session.model.node(1).is_some(), outcome != "success", "{outcome}");
            assert!(app.picker.is_none());
        }
        let mut kept = Session::default();
        kept.execute("file.open", &json!({"path": dir.join("kept.fieldcraft")})).unwrap();
        assert!(kept.model.node(1).is_some(), "replacement only follows a real save");
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn asynchronous_deck_and_csv_exports_keep_names_parameters_and_extensions() {
        let dir = temp("async-exports");
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("solve.static", json!({})).unwrap();
        let chooser = deferred(&mut app);
        app.write_deck("abaqus");
        let deck = dir.join("beam.inp");
        *chooser.answer.borrow_mut() = Some(Some(dir.join("beam").to_string_lossy().into_owned()));
        app.poll_picker();
        assert!(std::fs::read_to_string(deck).unwrap().to_ascii_uppercase().contains("*NODE"));
        assert_eq!(
            chooser.requests.borrow()[0],
            FilePicker::Save { name: "Cantilever (sample).inp".into(), filters: vec![Filter { name: "Solver deck", extensions: &["inp"] }] }
        );
        app.export_csv(json!({"what": "reactions", "step": "Load"}));
        let csv = dir.join("reactions.CSV");
        *chooser.answer.borrow_mut() = Some(Some(csv.to_string_lossy().into_owned()));
        app.poll_picker();
        let text = std::fs::read_to_string(csv).unwrap();
        assert!(text.starts_with("set,nodes,fx [N]") && text.contains("\r\ntotal,"));
        assert_eq!(
            chooser.requests.borrow()[1],
            FilePicker::Save { name: "Cantilever (sample)-Load-reactions.csv".into(), filters: vec![Filter { name: "CSV", extensions: &["csv"] }] }
        );
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn save_before_opening_keeps_both_models_and_rejects_old_choices_after_open() {
        let dir = temp("async-before-open");
        let other = dir.join("other.fieldcraft");
        let kept = dir.join("kept.fieldcraft");
        let mut saved = Session::default();
        saved.execute("mesh.node.add", &json!({"label": 2, "xyz": [0, 0, 0]})).unwrap();
        saved.execute("file.save_as", &json!({"path": other})).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let chooser = deferred(&mut app);
        app.open_path(&other.to_string_lossy());
        app.answer_unsaved(Answer::Save);
        *chooser.answer.borrow_mut() = Some(Some(kept.to_string_lossy().into_owned()));
        app.poll_picker();
        app.wait_loaded();
        assert!(app.session.model.node(2).is_some() && app.session.model.node(1).is_none());
        saved.execute("file.open", &json!({"path": kept})).unwrap();
        assert!(saved.model.node(1).is_some());

        app.save_as();
        app.run("file.open", json!({"path": kept})).unwrap();
        let stale = dir.join("stale.fieldcraft");
        *chooser.answer.borrow_mut() = Some(Some(stale.to_string_lossy().into_owned()));
        app.poll_picker();
        assert!(!stale.exists() && app.session.model.node(1).is_some());
        let _ = std::fs::remove_dir_all(dir);
    }

    #[test]
    fn save_before_update_only_closes_after_a_completed_save() {
        let dir = temp("async-before-update");
        let ctx = egui::Context::default();
        for cancelled in [true, false] {
            let (_tx, events) = crate::inbox::channel();
            let mut app = ModelApp::new(Session::default()).with_updates(crate::updates::UpdateHost {
                unsupported: None,
                automatic: false,
                start: Box::new(|_| {}),
                events,
                now: Box::new(|| 0),
                open_file: Box::new(|_| Ok(())),
            });
            app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
            let chooser = deferred(&mut app);
            app.answer_update(true);
            app.answer_unsaved(Answer::Save);
            app.poll_updates(&ctx);
            assert!(!app.quit_requested);
            *chooser.answer.borrow_mut() = Some((!cancelled).then(|| dir.join("kept.fieldcraft").to_string_lossy().into_owned()));
            app.poll_picker();
            app.poll_updates(&ctx);
            assert_eq!(app.quit_requested, !cancelled);
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    /// The navigation preset and the wheel's direction are remembered in the preferences file;
    /// an older file with keys no longer used still loads, a malformed one gives the defaults.
    #[test]
    fn navigation_preferences_are_remembered() {
        let file = temp("prefs").join("prefs.json");
        std::fs::write(&file, br#"{"show_setup": false}"#).unwrap();
        let store = PrefsStore::file(file.clone());
        assert_eq!(Prefs::load(Some(&store)), Prefs::default());
        let services = || Services { prefs: Some(PrefsStore::file(file.clone())), ..Default::default() };
        let mut app = ModelApp::new(Session::default()).with_services(services());
        app.ui.nav_preset = "solidworks".into();
        app.ui.reverse_zoom = true;
        app.save_prefs();
        let again = ModelApp::new(Session::default()).with_services(services());
        assert_eq!((again.ui.nav_preset.as_str(), again.ui.reverse_zoom), ("solidworks", true));
        // The display units too: the session shows and takes numbers in them.
        let mut app = again;
        let si = fieldcraft_model::DisplaySystem::new(fieldcraft_model::units::Preset::MKgS);
        app.set_default_units(si.clone());
        let again = ModelApp::new(Session::default()).with_services(services());
        assert_eq!((&again.prefs().units, again.session.units()), (&si, &si));
        std::fs::write(&file, b"{nope").unwrap();
        assert_eq!(Prefs::load(Some(&store)), Prefs::default());
    }

    /// Save With Geometry Embedded saves through the picker with `geometry: "embed"`, the file then
    /// holding a copy of the design; without a picker the form opens with it chosen.
    #[test]
    fn save_with_geometry_embedded_packs_the_design() {
        let dir = temp("packed");
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("geometry.open", json!({"path": design})).unwrap();
        app.save_packed();
        let d = app.dialog.as_ref().unwrap();
        assert_eq!((d.command, serde_json::from_str::<Value>(&d.text).unwrap()["geometry"].as_str()), ("file.save_as", Some("embed")));
        let target = dir.join("packed.fieldcraft").to_string_lossy().into_owned();
        let answer = target.clone();
        let mut app =
            app.with_services(Services { pick_save: Some(Box::new(move |_: &str, _: &[Filter]| Some(answer.clone()))), ..Default::default() });
        assert!(app.save_packed());
        let mut opened = Session::default();
        opened.execute("file.open", &json!({"path": target})).unwrap();
        assert!(opened.pack_geometry && opened.embedded.len() == 1, "the file holds a copy of the design");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A picker that answers `path` and records the filters it was offered.
    fn answering(path: Option<String>, offered: Rc<RefCell<Vec<Vec<Filter>>>>) -> PickOpen {
        Box::new(move |f: &[Filter]| {
            offered.borrow_mut().push(f.to_vec());
            path.clone()
        })
    }

    /// Original, dimensionless model data; no commercial deck or history is reused.
    fn input_deck(dir: &Path) -> String {
        let path = dir.join("original.inp");
        std::fs::write(&path, "*NODE\n1,1,0,0\n2,2,0,0\n3,1,1,0\n4,1,0,1\n*ELEMENT,TYPE=C3D4,ELSET=Solid\n10,1,2,3,4\n*MATERIAL,NAME=Steel\n*ELASTIC\n1000,0.3\n*SOLID SECTION,ELSET=Solid,MATERIAL=Steel\n*STEP\n*STATIC\n*END STEP\n").unwrap();
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn deck_form_confirms_replacement_and_keeps_the_import_report() {
        let dir = temp("deck-form");
        let path = input_deck(&dir);
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label":99,"xyz":[0,0,0]})).unwrap();
        app.open_form("deck.read", json!({"format":"abaqus","path":path,"units":"mm-t-s-N-MPa"}));
        app.submit_dialog();
        assert!(matches!(app.prompt, Some(Prompt::Unsaved(Replace::Deck(_)))));
        assert!(app.loading.is_none());
        app.answer_unsaved(Answer::Cancel);
        assert!(app.session.model.node(99).is_some() && app.dialog.is_some());
        app.submit_dialog();
        app.answer_unsaved(Answer::Discard);
        assert!(app.loading.is_some());
        app.wait_loaded();
        assert!(app.session.model.node(99).is_none());
        assert_eq!(app.session.model.elements().len(), 1);
        assert!(app.session.is_dirty() && app.session.path.is_none() && app.dialog.is_none());
        let reply = app.reply.as_ref().unwrap();
        assert_eq!(reply.command, "deck.read");
        assert_eq!(reply.value["history_imported"], false);
        assert_eq!(
            reply.value["unknown"],
            json!([
                {"keyword":"STEP","line":12}, {"keyword":"STATIC","line":13}, {"keyword":"ENDSTEP","line":14}
            ])
        );
        assert!(app.status.as_ref().unwrap().text.contains("3 unread keywords"));
        assert!(app.session.journal.is_empty());
        app.new_model();
        assert_eq!(app.prompt, Some(Prompt::Unsaved(Replace::New)), "an untouched import is still unsaved work");
        app.answer_unsaved(Answer::Cancel);
        assert_eq!(app.session.model.elements().len(), 1);
    }

    #[test]
    fn fresh_deck_jobs_freeze_the_input_preset_without_prompting_programmatic_calls() {
        let dir = temp("deck-units");
        let path = input_deck(&dir);
        let mut app = ModelApp::new(Session::default());
        app.run("model.units", json!({"system":"in-lbf-s-psi","overrides":{"length":"mm"}})).unwrap();
        app.run("mesh.node.add", json!({"label":99,"xyz":[0,0,0]})).unwrap();
        app.start_job("deck.read", json!({"format":"abaqus","path":path,"units":null})).unwrap();
        assert!(app.prompt.is_none());
        app.wait_loaded();
        let inch = fieldcraft_model::units::Quantity::Length.parse("1 in", &fieldcraft_model::DisplaySystem::default()).unwrap();
        assert_eq!(app.session.model.node(1).unwrap().xyz, [inch, 0., 0.]);
        assert!(app.session.model.node(99).is_none());
        assert_eq!(app.reply.as_ref().unwrap().value["units"], "in-lbf-s-psi");
    }

    #[test]
    fn deck_import_protection_survives_failed_or_cancelled_jobs_and_stops_after_saving() {
        let dir = temp("deck-marker");
        let path = input_deck(&dir);
        let mut app = ModelApp::new(Session::default());
        app.start_job("deck.read", json!({"format":"abaqus","path":path})).unwrap();
        app.wait_loaded();
        let imported = app.session.opened();
        assert_eq!(app.imported_document, Some(imported));
        assert!(app.has_unsaved_work());
        app.start_job("deck.read", json!({"format":"abaqus","path":dir.join("missing.inp")})).unwrap();
        app.wait_loaded();
        assert_eq!(app.session.opened(), imported);
        assert!(app.has_unsaved_work());
        app.start_job("deck.read", json!({"format":"abaqus","path":path})).unwrap();
        app.cancel_loading();
        assert_eq!(app.session.opened(), imported);
        assert!(app.has_unsaved_work());
        app.new_model();
        assert_eq!(app.prompt, Some(Prompt::Unsaved(Replace::New)));
        app.answer_unsaved(Answer::Cancel);
        app.run("file.save_as", json!({"path":dir.join("imported.fieldcraft")})).unwrap();
        assert!(!app.has_unsaved_work());
        app.new_model();
        assert!(app.prompt.is_none() && app.session.model.nodes().is_empty());
        assert_ne!(app.session.opened(), imported);
        assert!(!app.has_unsaved_work());
    }

    #[test]
    fn failed_deck_job_keeps_the_model_and_the_form_for_correction() {
        let dir = temp("deck-failure");
        let path = dir.join("broken.inp");
        std::fs::write(&path, "*NODE\n1,0,0,0\n*ELEMENT,TYPE=C3D4\n10,1,2,3,4\n").unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label":99,"xyz":[0,0,0]})).unwrap();
        let before = app.session.model.clone();
        let undo = app.session.undo_depth();
        app.open_form("deck.read", json!({"format":"abaqus","path":path}));
        app.submit_dialog();
        app.answer_unsaved(Answer::Discard);
        app.wait_loaded();
        assert_eq!(app.session.model, before);
        assert_eq!(app.session.undo_depth(), undo);
        assert!(app.session.is_dirty() && app.status.as_ref().unwrap().error);
        assert!(app.dialog.as_ref().unwrap().error.is_some());
        app.run("edit.undo", json!({})).unwrap();
        assert!(app.session.model.node(99).is_none());
    }

    #[test]
    fn file_types_dispatch_to_their_command() {
        assert_eq!(Kind::of("a/b.fieldcraft"), Some(Kind::Model));
        assert_eq!(Kind::of("B.SOLIDCRAFT"), Some(Kind::Design));
        assert_eq!(Kind::of("part.3MF").map(Kind::command), Some("geometry.import.mesh"));
        assert_eq!(Kind::of("part.MSH").map(Kind::command), Some("mesh.import.orphan"));
        assert!(is_long("mesh.import.orphan"));
        assert_eq!(Kind::of("part.step"), None);
        assert_eq!(Kind::of("noext"), None);
        assert_eq!(with_extension("a/b".into(), "fieldcraft"), "a/b.fieldcraft");
        assert_eq!(with_extension("a/b.FIELDCRAFT".into(), "fieldcraft"), "a/b.FIELDCRAFT");
    }

    #[test]
    fn orphan_picker_imports_in_background_and_failed_replacement_keeps_the_model() {
        let dir = temp("orphan");
        let path = dir.join("input.msh").to_string_lossy().into_owned();
        // Original one-tet ASCII fixture, written directly from the public format grammar.
        std::fs::write(&path, "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$Nodes\n1 4 1 4\n3 1 0 4\n1\n2\n3\n4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n$EndNodes\n$Elements\n1 1 7 7\n3 1 4 1\n7 1 2 3 4\n$EndElements\n").unwrap();
        let offered = Rc::new(RefCell::new(Vec::new()));
        let mut app = ModelApp::new(Session::default());
        app.services.pick_open = Some(answering(Some(path.clone()), offered.clone()));
        app.start("mesh.import.orphan").unwrap();
        assert!(app.dialog.is_none());
        assert_eq!(offered.borrow()[0][0], ORPHAN);
        assert_eq!(app.busy(), Some("Importing input.msh"));
        app.wait_loaded();
        assert_eq!(app.session.model.elements().len(), 1);
        assert!(app.session.model.element(7).is_some());
        app.run("edit.undo", serde_json::json!({})).unwrap();
        assert!(app.session.model.nodes().is_empty(), "background import remains one undo step");
        app.run("edit.redo", serde_json::json!({})).unwrap();
        app.run("model.set.create", serde_json::json!({"name":"selection", "nodes":[1]})).unwrap();
        let before = app.session.model.clone();
        let journal = app.session.journal.clone();
        app.open_path(&path);
        app.wait_loaded();
        assert_eq!(app.session.model, before);
        assert_eq!(app.session.journal, journal);
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.contains("label-based set `selection`")));
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_dir(&dir);
    }

    /// Open (picker) → background open → the model, the recent list; Save As (picker) and Save.
    #[test]
    fn the_pickers_open_and_save_models() {
        let dir = temp("models");
        let saved = dir.join("beam.fieldcraft").to_string_lossy().to_string();
        let mut first = ModelApp::new(Session::default());
        first.run_script(&fieldcraft_engine::sample::script()).unwrap();
        // Save As through the save picker: the extension is added, the file becomes the model's.
        let bare = dir.join("beam").to_string_lossy().to_string();
        first.services.pick_save = Some(Box::new(move |name: &str, f: &[Filter]| {
            assert_eq!((name, f.first().map(|f| f.extensions)), ("Untitled.fieldcraft", Some(&["fieldcraft"][..])));
            Some(bare.clone())
        }));
        first.start("file.save_as").unwrap();
        assert_eq!(first.session.path.as_deref(), Some(saved.as_str()));
        assert!(!first.session.is_dirty());
        assert_eq!(first.recent.paths, vec![saved.clone()]);

        let offered = Rc::new(RefCell::new(Vec::new()));
        let mut app = ModelApp::new(Session::default());
        app.services.recent_path = None;
        app.services.pick_open = Some(answering(Some(saved.clone()), offered.clone()));
        app.start("file.open").unwrap();
        assert!(app.dialog.is_none(), "a picker, not the form");
        assert_eq!(offered.borrow()[0][0], MODEL, "Open Model offers models first");
        assert_eq!(app.busy(), Some("Opening beam.fieldcraft"));
        app.wait_loaded();
        assert_eq!(app.session.model.elements().len(), first.session.model.elements().len());
        assert_eq!(app.session.path.as_deref(), Some(saved.as_str()));
        assert_eq!(app.recent.paths, vec![saved.clone()]);
        assert_eq!(app.status.as_ref().map(|s| (s.text.as_str(), s.error)), Some(("Opened beam.fieldcraft", false)));
        // Save goes to the model's own file without asking.
        app.services.pick_save = Some(Box::new(|_: &str, _: &[Filter]| -> Option<String> { panic!("Save asked for a file") }));
        app.run("mesh.node.add", serde_json::json!({"label": 99999, "xyz": [0, 0, 0]})).unwrap();
        assert!(app.save());
        assert!(!app.session.is_dirty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The icons a few frames of `draw` ask the kit for, with whether each has artwork.
    pub(crate) fn icons_drawn(ctx: &egui::Context, mut draw: impl FnMut(&egui::Context)) -> Vec<(String, bool)> {
        crate::icons::install();
        solidcraft_ui_kit::icons::audit::start();
        for _ in 0..3 {
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| draw(ui.ctx()));
            out.textures_delta.clear();
        }
        let drawn = solidcraft_ui_kit::icons::audit::take().into_iter().map(|r| (r.name, !r.placeholder)).collect();
        solidcraft_ui_kit::icons::audit::stop();
        drawn
    }

    /// The Delete confirmation draws the kit's delete artwork, not an empty icon (K11).
    #[test]
    fn the_delete_prompt_draws_its_own_icon() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        let load = app.session.model.loads().keys().next().cloned().unwrap();
        app.ask_delete(crate::browser::Def::Load, &load);
        let ctx = egui::Context::default();
        crate::theme::apply(&ctx, true);
        let drawn = icons_drawn(&ctx, |ctx| show(&mut app, ctx));
        assert!(drawn.contains(&(DELETE_ICON.to_string(), true)), "{drawn:?}");
    }

    /// File › Run Script… runs a JSON script the picker chose, as `fieldcraft-cli run` and the
    /// control channel's `engine.script` do; a file that is not a script says why and changes
    /// nothing.
    #[test]
    fn run_script_runs_the_picked_json_script() {
        let dir = temp("script");
        let script = dir.join("beam.json");
        std::fs::write(&script, serde_json::to_string(&fieldcraft_engine::sample::script()).unwrap()).unwrap();
        let bad = dir.join("bad.json");
        std::fs::write(&bad, "not json").unwrap();
        let offered = Rc::new(RefCell::new(Vec::new()));
        let mut app = ModelApp::new(Session::default());
        app.services.pick_open = Some(answering(Some(script.to_string_lossy().into_owned()), offered.clone()));
        app.run_script_dialog();
        assert_eq!(offered.borrow()[0][0], SCRIPT, "a script picker offers JSON scripts");
        assert!(!app.session.model.elements().is_empty(), "the script built the model");
        assert_eq!(app.status.as_ref().map(|s| (s.text.as_str(), s.error)), Some(("Script done", false)));
        let journal = app.session.journal.len();
        app.services.pick_open = Some(answering(Some(bad.to_string_lossy().into_owned()), offered));
        app.run_script_dialog();
        assert_eq!(app.session.journal.len(), journal);
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.starts_with("bad.json:")), "{:?}", app.status);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn cancelled_pickers_and_failed_opens_leave_the_model() {
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", serde_json::json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        let offered = Rc::new(RefCell::new(Vec::new()));
        app.services.pick_open = Some(answering(None, offered.clone()));
        app.open_dialog(None);
        assert_eq!(offered.borrow()[0], OPEN_FILTERS.to_vec(), "File > Open offers every type");
        assert!(app.loading.is_none());
        // A file that does not open: the model stays, the error is shown.
        let dir = temp("broken");
        let broken = dir.join("broken.solidcraft");
        std::fs::write(&broken, b"not a design").unwrap();
        app.open_path(&broken.to_string_lossy());
        app.wait_loaded();
        assert!(app.session.model.node(1).is_some());
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.starts_with("broken.solidcraft:")), "{:?}", app.status);
        assert!(app.recent.paths.is_empty());
        app.open_path("drawing.dxf");
        assert!(app.status.as_ref().is_some_and(|s| s.error));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn opening_a_model_over_unsaved_changes_asks_first() {
        let dir = temp("unsaved");
        let path = dir.join("other.fieldcraft").to_string_lossy().to_string();
        let mut other = ModelApp::new(Session::default());
        other.run("mesh.node.add", serde_json::json!({"label": 5, "xyz": [0, 0, 0]})).unwrap();
        other.run("file.save_as", serde_json::json!({"path": path})).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", serde_json::json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        app.open_path(&path);
        assert_eq!(app.prompt, Some(Prompt::Unsaved(Replace::File(path.clone()))));
        assert!(app.loading.is_none());
        app.answer_unsaved(Answer::Cancel);
        assert!(app.prompt.is_none() && app.loading.is_none() && app.session.model.node(1).is_some());
        app.open_path(&path);
        app.answer_unsaved(Answer::Discard);
        app.wait_loaded();
        assert!(app.session.model.node(5).is_some() && app.session.model.node(1).is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A design of another part over a meshed model asks; each answer does what it says, and the
    /// same design opens again without asking.
    #[test]
    fn opening_a_design_of_another_part_asks_what_to_do_with_the_model() {
        let dir = temp("design");
        let path = dir.join("plate.solidcraft").to_string_lossy().to_string();
        std::fs::write(&path, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let cube = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        let meshed = || {
            let mut app = ModelApp::new(Session::default());
            app.run("mesh.block", serde_json::json!({"body": "Seal", "kind": "hex8", "corners": cube, "divisions": [1, 1, 1]})).unwrap();
            app.run("model.set.create", serde_json::json!({"name": "seal", "geometry": [{"body": "Seal"}]})).unwrap();
            app
        };
        let opened = |app: &mut ModelApp| {
            app.open_path(&path);
            app.wait_loaded();
        };

        // Cancel: nothing changes.
        let mut app = meshed();
        opened(&mut app);
        assert_eq!(app.prompt, Some(Prompt::Design(path.clone())));
        assert!(app.session.design.is_none(), "the model waits for the answer");
        app.answer_design(DesignAnswer::Cancel);
        assert!(app.prompt.is_none() && app.held.is_none() && app.session.design.is_none() && app.session.model.elements().len() == 1);

        // Keep: the design over the model, and the status says what does not resolve.
        let mut app = meshed();
        opened(&mut app);
        app.answer_design(DesignAnswer::Keep);
        assert!(app.session.design.is_some() && app.session.model.elements().len() == 1);
        assert!(app.status.as_ref().is_some_and(|s| s.text.contains("not on this design: mesh of Seal; sets seal")), "{:?}", app.status);
        assert_eq!(app.recent.paths, vec![path.clone()]);

        // New: an empty model with the design, one undo step from the old one.
        let mut app = meshed();
        opened(&mut app);
        app.answer_design(DesignAnswer::New);
        app.wait_loaded();
        assert!(app.prompt.is_none() && app.session.design.is_some());
        assert!(app.session.model.elements().is_empty() && app.session.model.sets().is_empty());
        app.run("edit.undo", serde_json::json!({})).unwrap();
        assert!(app.session.design.is_none() && app.session.model.elements().len() == 1);

        // The same design again (its bodies are the model's): no question.
        app.run("edit.redo", serde_json::json!({})).unwrap();
        app.run("model.set.create", serde_json::json!({"name": "plate", "geometry": [{"body": "Body1"}]})).unwrap();
        opened(&mut app);
        assert!(app.prompt.is_none() && app.session.model.set("plate").is_some());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A mesher's status says what it made: the octree mesher its count, worst scaled Jacobian,
    /// and how deep it refined for the thinnest wall, as the tet mesher says its count (A4-14).
    #[test]
    fn a_finished_mesher_says_what_it_made() {
        let units = fieldcraft_model::DisplaySystem::default();
        let octree = serde_json::json!({"ok": true, "result": {"elements": 154710, "nodes": 170000, "bodies": [
            {"body": "Case", "kind": "hex8", "elements": 154710, "quality": {"min_scaled_jacobian": 0.3812},
             "grid": {"levels": 4, "thinnest_wall": 2.0}}]}});
        assert_eq!(
            finished("mesh.hex.octree", "Octree Hex Mesh done", "1 min 10 s", &octree, &units),
            "Octree Hex Mesh done: 154,710 hex8, smallest scaled Jacobian 0.381, 4 octree levels (thinnest wall 2 mm) in 1 min 10 s — Mesh › Quality checks them"
        );
        let tet =
            serde_json::json!({"ok": true, "result": {"elements": 3125, "kind": "tet4", "bodies": [{"quality": {"mesh": {"min_dihedral": 36.54}}}]}});
        assert_eq!(
            finished("mesh.tet.generate", "Tet Mesh done", "2.2 s", &tet, &units),
            "Tet Mesh done: 3,125 tet4, smallest dihedral angle 36.5° in 2.2 s — Mesh › Quality checks them"
        );
    }

    /// A file opened from its form reads as a dropped one does: "Opened plate.solidcraft" (A4-6).
    #[test]
    fn a_file_opened_from_its_form_says_which() {
        let dir = temp("form-open");
        let path = dir.join("plate.solidcraft").to_string_lossy().to_string();
        std::fs::write(&path, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.submit_background("geometry.open", serde_json::json!({"path": path}));
        app.wait_loaded();
        assert_eq!(app.status.as_ref().map(|s| s.text.as_str()), Some("Opened plate.solidcraft"));
        assert_eq!(app.recent.paths, vec![path.clone()]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Saving and writing a deck finish with their sentence in the status bar, not the reply's
    /// fields (A4-6).
    #[test]
    fn a_save_and_a_deck_say_what_they_wrote() {
        let dir = temp("sentences");
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let status = |app: &ModelApp| app.status.as_ref().map(|s| s.text.clone()).unwrap_or_default();
        let model = dir.join("beam.fieldcraft").to_string_lossy().to_string();
        app.run("file.save_as", serde_json::json!({"path": model})).unwrap();
        assert_eq!(status(&app), "Saved beam.fieldcraft");
        let deck = dir.join("beam.inp").to_string_lossy().to_string();
        app.run("deck.write", serde_json::json!({"format": "abaqus", "path": deck})).unwrap();
        assert!(status(&app).starts_with("Wrote beam.inp (Abaqus"), "{}", status(&app));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A model opened without results leaves the Results workspace, where every tool would be
    /// greyed: a design opened as a new model comes up in Prepare (A4-15).
    #[test]
    fn a_model_without_results_does_not_open_in_the_results_workspace() {
        let dir = temp("workspace");
        let path = dir.join("plate.solidcraft").to_string_lossy().to_string();
        std::fs::write(&path, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.ui.select_workspace(Workspace::Results);
        app.start_job("geometry.open", serde_json::json!({"path": path, "model": "new"})).unwrap();
        app.wait_loaded();
        assert!(app.session.design.is_some() && app.session.model.elements().is_empty());
        assert_eq!(app.ui.workspace, Workspace::Prepare);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn decks_are_written_where_the_picker_says() {
        let dir = temp("decks");
        let target = dir.join("cantilever").to_string_lossy().to_string();
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let t = target.clone();
        app.services.pick_save = Some(Box::new(move |name: &str, f: &[Filter]| {
            assert_eq!((name, f.first().map(|f| f.extensions)), ("Untitled.inp", Some(&["inp"][..])));
            Some(t.clone())
        }));
        app.start("deck.write.abaqus").unwrap();
        assert!(app.dialog.is_none());
        let text = std::fs::read_to_string(format!("{target}.inp")).unwrap();
        assert!(text.contains("*NODE") || text.contains("*Node"), "{}", &text[..text.len().min(200)]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Export CSV… asks the picker for `<model>-<step>-reactions.csv` and writes it; without a
    /// picker the command's form opens (#92).
    #[test]
    fn reactions_export_where_the_picker_says() {
        let dir = temp("csv");
        let target = dir.join("sums").to_string_lossy().to_string();
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("solve.static", json!({})).unwrap();
        let params = json!({"what": "reactions", "step": "Load"});
        app.export_csv(params.clone());
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("results.export.csv"), "no picker: the form");
        app.dialog = None;
        let t = target.clone();
        app.services.pick_save = Some(Box::new(move |name: &str, f: &[Filter]| {
            assert_eq!((name, f.first().map(|f| f.extensions)), ("Cantilever (sample)-Load-reactions.csv", Some(&["csv"][..])));
            Some(t.clone())
        }));
        app.export_csv(params);
        let text = std::fs::read_to_string(format!("{target}.csv")).unwrap();
        // The headers state the units (the display units, here mm-t-s-N-MPa).
        assert!(text.starts_with("set,nodes,fx [N],fy [N],fz [N],mx [N·mm],my [N·mm],mz [N·mm]\r\n") && text.contains("\r\ntotal,"), "{text}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn recovered_entries_say_how_long_ago() {
        let cases = [
            (0, "just now"),
            (59, "just now"),
            (60, "1 min ago"),
            (3599, "59 min ago"),
            (7200, "2 h ago"),
            (90_000, "1 day ago"),
            (13 * 86_400, "13 days ago"),
        ];
        for (s, want) in cases {
            assert_eq!(ago(s), want, "{s}");
        }
    }

    #[test]
    fn the_recent_list_is_kept_newest_first_and_survives_bad_files() {
        let dir = temp("recent");
        let file = dir.join("recent.json");
        let mut r = Recent::load(Some(file.clone()));
        assert!(r.paths.is_empty());
        for i in 0..12 {
            r.add(&format!("/m/{i}.fieldcraft"));
        }
        r.add("/m/3.fieldcraft");
        let again = Recent::load(Some(file.clone()));
        assert_eq!(again.paths.len(), MAX_RECENT);
        assert_eq!(again.paths.first().map(String::as_str), Some("/m/3.fieldcraft"));
        assert_eq!(again.paths.iter().filter(|p| *p == "/m/3.fieldcraft").count(), 1);
        std::fs::write(&file, b"{not json").unwrap();
        assert!(Recent::load(Some(file.clone())).paths.is_empty());
        std::fs::write(&file, br#"["/x/a.fieldcraft", "/x/b.exe", ""]"#).unwrap();
        assert_eq!(Recent::load(Some(file)).paths, vec!["/x/a.fieldcraft".to_string()]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A sample opens in the background with nothing to save (building and meshing one takes
    /// seconds: the window keeps drawing, #552); Solve runs its steps in the background and the
    /// Results tab comes up with the first step's result.
    #[test]
    fn a_sample_solves_into_the_results_view() {
        let mut app = ModelApp::new(Session::default());
        app.start("file.sample").unwrap();
        assert_eq!(app.prompt, Some(Prompt::Samples));
        app.prompt = None;
        app.open_sample("cantilever");
        assert_eq!(app.busy(), Some("Opening the Cantilever sample"));
        assert!(app.session.model.nodes().is_empty(), "the window's model stays until the job is done");
        app.wait_loaded();
        assert!(!app.session.is_dirty() && !app.session.model.steps().is_empty());
        assert_eq!(app.status.as_ref().map(|s| s.text.as_str()), Some("Opened the Cantilever sample — Solve shows its results"));
        app.start("solve.run").unwrap();
        assert_eq!(app.busy(), Some("Solving"));
        assert!(app.run("model.info", serde_json::json!({})).is_err(), "commands wait for the solve");
        app.wait_loaded();
        assert!(app.status.as_ref().is_some_and(|s| s.text.starts_with("Solved 2 steps in ")), "{:?}", app.status);
        app.logic(&egui::Context::default());
        assert_eq!(app.ui.tab, fieldcraft_engine::catalogue::Tab::Results);
        assert!(app.plot().is_some());
        // Over unsaved changes, opening a sample asks first.
        app.run("mesh.node.add", serde_json::json!({"label": 900000, "xyz": [0, 0, 0]})).unwrap();
        app.open_sample("impact");
        assert_eq!(app.prompt, Some(Prompt::Unsaved(Replace::Sample("impact"))));
        app.answer_unsaved(Answer::Discard);
        assert_eq!(app.busy(), Some("Opening the Bar Impact sample"));
        app.wait_loaded();
        assert!(app.session.model.step("Impact").is_some());
        // An unknown sample name opens the samples to choose from, and says why.
        app.open_sample("crank-armm");
        assert_eq!(app.prompt, Some(Prompt::Samples));
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.contains("crank-armm")), "{:?}", app.status);
        app.prompt = None;
        // File › New over unsaved changes asks first too, then empties the model.
        app.run("mesh.node.add", serde_json::json!({"label": 900000, "xyz": [0, 0, 0]})).unwrap();
        app.start("file.new").unwrap();
        assert_eq!(app.prompt, Some(Prompt::Unsaved(Replace::New)));
        app.answer_unsaved(Answer::Discard);
        assert!(app.session.model.nodes().is_empty() && !app.session.is_dirty());
        // A model without steps says so instead of solving nothing.
        let mut empty = ModelApp::new(Session::default());
        empty.solve_all();
        assert!(empty.busy().is_none() && empty.status.as_ref().is_some_and(|s| s.error));
    }

    /// A cancelled or refused solve leaves the window in the workspace it was in, and a job
    /// waiting over the control channel hears that it was cancelled (#324, #367).
    #[test]
    fn a_cancelled_or_refused_solve_stays_in_its_workspace() {
        use fieldcraft_engine::catalogue::Workspace;
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.ui.select_workspace(Workspace::Model);
        let job = app.start_job("solve.run", serde_json::json!({})).unwrap();
        assert_eq!(job, 1);
        let (tx, rx) = std::sync::mpsc::channel();
        app.loading.as_mut().unwrap().waiters.push(tx);
        assert_eq!(app.loading.as_ref().unwrap().to_json()["command"], "solve.run");
        app.cancel_loading();
        app.logic(&egui::Context::default());
        assert_eq!(app.ui.workspace, Workspace::Model);
        assert_eq!(rx.recv().unwrap()["cancelled"], true);
        let mut empty = ModelApp::new(Session::default());
        empty.ui.select_workspace(Workspace::Model);
        assert!(empty.start_job("solve.run", serde_json::json!({})).is_err());
        empty.solve_all();
        empty.logic(&egui::Context::default());
        assert_eq!(empty.ui.workspace, Workspace::Model);
    }

    /// A job's time comes from a clock, not the window's frame time, which stands still while no
    /// frames are drawn (a hidden or covered window): a solve that took 6 s said "done in 0.0 s"
    /// (#590). Here no frame is drawn at all (`now` stays 0) while the job runs 0.3 s or more.
    #[test]
    fn a_job_reports_its_time_without_frames() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        let frozen = app.now;
        app.start_job("solve.run", serde_json::json!({})).unwrap();
        std::thread::sleep(std::time::Duration::from_millis(300));
        assert!(app.loading.as_ref().unwrap().to_json()["elapsed"].as_f64().unwrap() >= 0.3);
        app.wait_loaded();
        assert_eq!(app.now, frozen);
        let status = app.status.clone().unwrap();
        let took: f64 = status.text.rsplit(" in ").next().and_then(|t| t.trim_end_matches(" s").parse().ok()).unwrap_or(0.0);
        assert!(status.text.starts_with("Solved ") && took >= 0.3, "{}", status.text);
    }

    /// A first-time user's path: Open picks a SolidCraft design, the Tet Mesh button opens its
    /// form, OK with the defaults meshes it in the background, and the mesh is there. A mesher
    /// that does not fit says why in the form, which stays open.
    #[test]
    fn a_picked_design_meshes_from_the_ribbon_with_the_defaults() {
        let dir = temp("path");
        let design = dir.join("plate.solidcraft").to_string_lossy().to_string();
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let offered = Rc::new(RefCell::new(Vec::new()));
        let mut app = ModelApp::new(Session::default());
        app.services.pick_open = Some(answering(Some(design.clone()), offered));
        app.open_dialog(None);
        app.wait_loaded();
        assert_eq!(app.session.design.as_ref().map(|d| d.bodies.len()), Some(1), "{:?}", app.status);
        assert!(!app.display().bounds.is_empty());

        app.start("mesh.tet.generate").unwrap();
        assert_eq!(app.dialog.as_ref().map(|d| (d.command, d.text.as_str())), Some(("mesh.tet.generate", r#"{"order":"linear","quality":true}"#)));
        app.submit_dialog();
        assert_eq!(app.busy(), Some("Generate Tet Mesh"), "meshing runs off the window thread");
        app.wait_loaded();
        assert!(app.dialog.is_none(), "the form closes when it succeeded: {:?}", app.status);
        // The status sums up the mesh made: elements, worst angle, time, where to check it.
        let status = app.status.clone().unwrap();
        assert!(
            !status.error && status.text.starts_with("Generate Tet Mesh done: ") && status.text.contains(" tet4, smallest dihedral angle "),
            "{}",
            status.text
        );
        let counts = app.session.model.element_counts();
        assert!(matches!(counts.as_slice(), [(fieldcraft_model::ElementKind::Tet4, n)] if *n > 0), "{counts:?}");
        assert_eq!(app.session.journal.len(), 2, "open and mesh, each an undo step");

        // Meshing the meshed body again asks first (#267); Cancel keeps the form and the mesh.
        app.start("mesh.tet.generate").unwrap();
        app.submit_dialog();
        let Some(Prompt::ReplaceMesh(r)) = app.prompt.clone() else { panic!("no question: {:?}", app.prompt) };
        let tets = counts.first().map_or(0, |c| c.1);
        assert!(r.question.starts_with("Body1 already has a mesh (") && r.question.ends_with(" tet4). Replace it?"), "{r:?}");
        assert!(tets > 999 && r.question.contains(','), "thousands are separated: {r:?}");
        assert!(app.busy().is_none());
        app.answer_replace_mesh(false);
        assert!(app.prompt.is_none() && app.dialog.is_some(), "Cancel leaves the form open");
        assert_eq!(app.session.journal.len(), 2);
        // Replace meshes it again, in place of the old mesh: as many tets, one undo step.
        app.submit_dialog();
        app.answer_replace_mesh(true);
        app.wait_loaded();
        assert!(app.dialog.is_none(), "{:?}", app.status);
        assert_eq!(app.session.model.element_counts(), counts, "replaced, not stacked");
        // The Mesh node shows the run that made it, once though it ran twice (#343).
        let made = app.mesh_history();
        assert_eq!(made.len(), 1, "{made:?}");
        assert_eq!(made[0].0, "mesh.tet.generate");
        assert_eq!(made[0].1["order"], "linear");

        // A plate with a hole is not one mapped block: the form stays with the reason.
        app.start("mesh.hex.mapped").unwrap();
        app.submit_dialog();
        app.answer_replace_mesh(true);
        app.wait_loaded();
        let error = app.dialog.as_ref().and_then(|d| d.error.clone()).unwrap_or_default();
        assert!(error.contains("not a logical cube"), "{error}");
        assert!(app.status.as_ref().is_some_and(|s| s.error));
        assert_eq!(app.session.model.element_counts(), counts, "a failed mesher changes nothing");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The partition, the cells mesher and the plan's apply run as jobs (#545): the window
    /// stays free while they run, the model they leave is the one the command makes run
    /// directly, bit for bit, in one undo step, and a cancelled one leaves the model as it was.
    #[test]
    fn partition_cells_and_plan_apply_run_as_jobs_like_the_command() {
        let dir = temp("jobs-545");
        let design = dir.join("plate.solidcraft").to_string_lossy().to_string();
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.open_path(&design);
        app.wait_loaded();
        app.run("geometry.datum", json!({"name": "X30", "by": {"type": "offset", "plane": "yz", "offset": 30}})).unwrap();
        let plan = app.run("mesh.plan", json!({"body": "Body1"})).unwrap();
        let steps = [
            ("mesh.plan.apply", json!({"plan": plan["name"]})),
            ("geometry.partition", json!({"name": "Half", "datum": "X30"})),
            ("mesh.cells", json!({"bodies": ["Body1/a", "Body1/b"], "size": 2.5})),
        ];
        for (command, params) in steps {
            assert!(is_long(command), "{command}");
            let before = app.session.clone();
            let mut direct = app.session.clone();
            let expected = direct.execute(command, &params).unwrap_or_else(|e| panic!("{command}: {e}"));

            // Cancelled: nothing changes.
            app.start_job(command, params.clone()).unwrap();
            app.cancel_loading();
            assert!(app.busy().is_none());
            assert!(app.session.model == before.model && app.session.journal.len() == before.journal.len(), "{command}: cancelled");

            app.start_job(command, params.clone()).unwrap();
            assert_eq!(app.busy(), Some(crate::command_label(command)), "{command} runs off the window thread");
            app.wait_loaded();
            assert!(app.status.as_ref().is_some_and(|s| !s.error), "{command}: {:?}", app.status);
            assert!(app.session.model == direct.model, "{command}: the job's model is the command's");
            let bodies = |s: &Session| s.design.as_ref().map(|d| d.bodies.iter().map(|b| b.name.clone()).collect::<Vec<_>>());
            assert_eq!(bodies(&app.session), bodies(&direct), "{command}");
            assert_eq!(app.session.journal.len(), before.journal.len() + 1, "{command}: one undo step");
            let elements = expected["elements"].as_u64().map(|n| n as usize);
            assert!(elements.is_none_or(|n| n == app.session.model.elements().len()), "{command}: {expected}");

            // Undo takes the whole command back; redo brings it again.
            app.run("edit.undo", json!({})).unwrap();
            assert!(app.session.model == before.model, "{command}: undo");
            app.run("edit.redo", json!({})).unwrap();
            assert!(app.session.model == direct.model, "{command}: redo");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The control channel's `ui.start` never opens a picker: it opens the form.
    #[test]
    fn programmatic_starts_never_open_a_picker() {
        let mut app = ModelApp::new(Session::default());
        app.services.pick_open = Some(Box::new(|_: &[Filter]| -> Option<String> { panic!("a picker opened") }));
        app.services.pick_save = Some(Box::new(|_: &str, _: &[Filter]| -> Option<String> { panic!("a picker opened") }));
        app.services.pick_file = Some(Box::new(|_| panic!("an asynchronous picker opened")));
        app.start_programmatic("file.open").unwrap();
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("file.open"));
        // A deck needs a model to write.
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.start_programmatic("deck.write.abaqus").unwrap();
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("deck.write"));
        // In-window choices are not pickers: they behave as a click.
        app.start_programmatic("file.sample").unwrap();
        assert_eq!(app.prompt, Some(Prompt::Samples));
        // Explicit paths retain the ordinary command result contract and never ask the host.
        let dir = temp("programmatic-save");
        let path = dir.join("model.fieldcraft");
        app.run("file.save_as", json!({"path": path})).unwrap();
        assert!(path.exists() && app.picker.is_none());
        let _ = std::fs::remove_dir_all(dir);
    }
}
