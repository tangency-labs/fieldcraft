//! The FieldCraft egui front end.
//!
//! A thin shell over [`fieldcraft_engine::Session`]: panels read engine state and act through
//! commands ([`ModelApp::run`]); nothing here owns model data. The window follows the Craft
//! family's layout (SolidCraft, VectorCraft):
//!
//! - the title bar ([`appbar`]): the app, quick actions (File menu, save, undo, redo, Solve,
//!   search, theme, help) and the model's document tab;
//! - the ribbon ([`toolbar`]): the workspace switcher, its tabs, and panels of icon buttons whose
//!   labelled drop-downs list the whole catalogue (coming-soon entries greyed with their
//!   milestone);
//! - the browser on the left ([`browser`]): the model's geometry, mesh and definitions;
//! - the viewport ([`viewport`]) with its navigation bar, orientation gizmo and selection bar
//!   (picking faces, edges, vertices, bodies, nodes and elements: [`selection`]);
//! - the context panel on the right ([`context`]), only while there is context: the open command
//!   form, the results controls, the selected definition's properties;
//! - the command line and status ([`palette`]), and the journal timeline at the bottom
//!   ([`timeline`]). Every command is also one search away (S or Cmd/Ctrl+K).
//!
//! The control channel's requests are handled here ([`control`]); the desktop app owns the
//! socket. Structure and conventions follow SolidCraft's `crates/ui-egui` (MIT OR Apache-2.0),
//! see ATTRIBUTION.md.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]
#![forbid(unsafe_code)]
// Slow work stays off the window thread (clippy.toml's `disallowed-methods`, #552); tests drive the
// session and wait for jobs directly.
#![cfg_attr(test, allow(clippy::disallowed_methods, clippy::disallowed_types))]

pub mod appbar;
pub mod browser;
pub mod capture;
pub mod context;
pub mod control;
pub mod describe;
pub mod elide;
pub mod feature_graph;
pub mod fields;
pub mod files;
pub mod form;
pub mod frame;
pub mod framing;
pub mod gpu;
pub mod icons;
pub mod inbox;
pub mod menu;
pub mod navigation;
pub mod node_move;
pub mod palette;
pub mod privacy;
pub mod reply;
pub mod report;
pub mod results;
pub mod results_panel;
pub mod samples;
pub mod scene;
pub mod selection;
pub mod theme;
pub mod timeline;
pub mod titlebar;
pub mod toolbar;
pub mod updates;
pub mod viewport;

use std::hash::{Hash, Hasher};
use std::sync::Arc;
use std::sync::mpsc::Sender;

use fieldcraft_engine::catalogue::{self, Tab, Workspace};
use fieldcraft_engine::view::display_tolerance;
use fieldcraft_engine::{EngineError, Session};
use fieldcraft_model::DisplaySystem;
use fieldcraft_render::glyphs::{self, Definition, Drawn, Glyphs};
use fieldcraft_render::{Aabb3, Camera, DisplayMesh, Scene, StandardView, Vec3};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

pub use control::ControlRequest;
pub use files::{Filter, Services};
pub use scene::Style;

/// The viewport height symbols are sized for before the viewport is first drawn (points).
const DEFAULT_VIEW_HEIGHT: f64 = 800.0;

/// The version shown in Help and over the control channel.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// View settings: what is shown and how, and which menu or search is open. `ui.set` changes any
/// of them by name.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct UiState {
    pub style: Style,
    pub show_nodes: bool,
    /// The open design and imported surfaces.
    pub show_geometry: bool,
    /// Elements.
    pub show_mesh: bool,
    /// Boundary condition and load symbols.
    pub show_loads: bool,
    /// Reference point, coupling and rigid wall symbols.
    pub show_interactions: bool,
    /// The setup's symbols over a results plot (those of the step shown); off by default, so
    /// they do not cover the contours.
    pub result_symbols: bool,
    /// Markers at a contour's largest and smallest value.
    pub result_extremes: bool,
    /// The setup symbols' size on screen, in points (held to
    /// [`fieldcraft_render::glyphs::SIZE_RANGE`]).
    pub symbol_size: f64,
    pub perspective: bool,
    pub dark: bool,
    /// The browser (left).
    pub show_tree: bool,
    /// The journal timeline (bottom).
    pub show_timeline: bool,
    /// The mouse navigation preset ([`navigation::NAV_PRESETS`] id; remembered per user).
    pub nav_preset: String,
    /// The wheel zooms the other way (remembered per user).
    pub reverse_zoom: bool,
    /// The workspace shown, and its toolbar tab.
    pub workspace: Workspace,
    pub tab: Tab,
    /// The command search (S, Cmd/Ctrl+K).
    pub palette: bool,
    /// The open menu: `file`, `help`, `workspace`, `view`, or a toolbar panel's drop-down as
    /// `<tab>/<panel>` (`mesh/Tet`).
    pub menu: Option<String>,
}

impl Default for UiState {
    fn default() -> Self {
        UiState {
            style: Style::ShadedEdges,
            show_nodes: false,
            show_geometry: true,
            show_mesh: true,
            show_loads: true,
            show_interactions: true,
            result_symbols: false,
            result_extremes: true,
            symbol_size: glyphs::DEFAULT_SIZE,
            perspective: false,
            dark: true,
            show_tree: true,
            show_timeline: true,
            nav_preset: navigation::DEFAULT_PRESET.into(),
            reverse_zoom: false,
            workspace: Workspace::default(),
            tab: Tab::default(),
            palette: false,
            menu: None,
        }
    }
}

impl UiState {
    /// Show `tab`, on the current workspace when it has it, else on the first that does.
    pub fn select_tab(&mut self, tab: Tab) {
        self.tab = tab;
        self.workspace = Workspace::showing(tab, self.workspace);
    }

    /// Show workspace `w` on its first tab.
    pub fn select_workspace(&mut self, w: Workspace) {
        self.workspace = w;
        self.tab = w.tabs().first().copied().unwrap_or(Tab::Tools);
    }

    /// What the viewport shows.
    pub fn show(&self) -> scene::Show {
        scene::Show {
            geometry: self.show_geometry,
            mesh: self.show_mesh,
            nodes: self.show_nodes,
            loads: self.show_loads,
            interactions: self.show_interactions,
        }
    }
}

/// What a left drag in the viewport does (the navigation bar's modes).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum NavMode {
    #[default]
    Orbit,
    Pan,
    Zoom,
}

/// The parameter form a command opens when started from the Commands panel.
#[derive(Clone, Debug, PartialEq)]
pub struct ParamsDialog {
    pub command: &'static str,
    /// The parameters as JSON text.
    pub text: String,
    pub error: Option<String>,
    /// The command refused one parameter: its path and the problem phrased for the form, shown
    /// beside its field.
    pub refused: Option<(String, String)>,
    /// The paths of the fields edited: their problems show; the others wait for OK.
    pub touched: std::collections::BTreeSet<String>,
    /// OK was tried with something missing: every problem shows.
    pub attempted: bool,
    /// The JSON text is shown instead of the typed fields.
    pub advanced: bool,
    /// Which opening of a form this is (counts up): a form opened over another one, even of the
    /// same command, gets fresh widget ids, so nothing typed into the one it replaces carries over.
    pub opening: u64,
    /// Opened from another form's "Create …" button: the form to go back to.
    pub back: Option<Box<Back>>,
    /// The set Use selection made of the selection (`model.set.create`'s parameters): created on
    /// OK when the form still names it, so Cancel leaves the model as it was.
    pub pending_set: Option<Value>,
}

/// A command that failed: its message, and the parameter it refused (path, form message).
struct Failure {
    message: String,
    refused: Option<(String, String)>,
}

/// The form a "Create …" button left, and where the item made goes in it.
#[derive(Clone, Debug, PartialEq)]
pub struct Back {
    pub dialog: ParamsDialog,
    pub collection: fieldcraft_engine::schema::Collection,
    /// The reference field's path in its parameters (`a/b`).
    pub path: String,
    /// The collection's names before: the one added is the new item.
    pub before: Vec<String>,
}

/// A status bar message.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Status {
    pub text: String,
    pub error: bool,
}

/// Selection by picking in the viewport ([`selection`]): the filter, what is selected and what
/// is under the pointer, and what they are drawn with.
#[derive(Default)]
pub struct PickState {
    pub filter: selection::Filter,
    pub selection: selection::Selection,
    /// What a click would pick under the pointer.
    pub hover: Option<selection::Item>,
    /// The faces of the feature whose row in the Feature Graph is under the pointer.
    feature_hover: Vec<selection::Item>,
    /// A left drag draws a selection box (nodes, element faces, elements) instead of orbiting.
    pub box_mode: bool,
    /// The box being dragged (window points).
    pub dragging: Option<(egui::Pos2, egui::Pos2)>,
    /// Where the pointer was when the hover was worked out (window points).
    pub hover_at: Option<egui::Pos2>,
    /// The hover is a body whose row in a command form is under the pointer, not the viewport's.
    pub form_hover: bool,
    /// Bumped whenever the selection changes, and whenever the hover does.
    generation: u64,
    hover_generation: u64,
    /// The geometry tags behind the drawn mesh and the model revision they were read at.
    index: (Option<u64>, Arc<selection::Index>),
    /// The names of the drawn geometry surfaces, in the order they were added to the display.
    names: Vec<selection::BodyNames>,
    /// The selection's highlight and what it was built for (kept while only the hover changes).
    selected_overlay: (Option<u64>, Arc<Scene>),
    /// The filter an open form switched to: the form's opening, the filter it set and the one
    /// to restore when it closes.
    form_filter: Option<(u64, selection::Filter, selection::Filter)>,
    /// The highlight overlay (the selection's and the hover's) and what it was built for.
    overlay: (Option<u64>, Arc<Scene>),
    /// What the hover was last worked out for: the pointer, the camera, the scene, the filter
    /// and the viewport.
    pub(crate) hover_key: Option<(Option<egui::Pos2>, Camera, u64, selection::Filter, Option<egui::Rect>)>,
    /// The model revision the selection's mesh labels were last checked at.
    checked: Option<u64>,
}

/// The viewport's renderer state.
#[derive(Default)]
pub struct Viewport {
    /// Where the viewport was drawn last frame (points).
    pub rect: Option<egui::Rect>,
    /// The GPU target when drawing with wgpu; `None` draws on the CPU.
    pub gpu: Option<gpu::GpuTarget>,
    pub slot: gpu::SceneSlot,
    /// The scene key last handed to the GPU.
    pub uploaded: Option<u64>,
    /// The selection highlight for the GPU, and the overlay key last handed to it.
    pub overlay_slot: gpu::SceneSlot,
    pub overlay_uploaded: Option<u64>,
    /// The setup's symbols for the GPU, and the symbols key last handed to it.
    pub symbols_slot: gpu::SceneSlot,
    pub symbols_uploaded: Option<u64>,
    /// CPU fallback: the texture and what it shows (scene, symbols and overlay keys, camera,
    /// size).
    cpu: Option<(egui::TextureHandle, (u64, u64, u64), Camera, [usize; 2])>,
    /// CPU fallback: the scene key the renderer last refused.
    cpu_failed: Option<u64>,
    /// GPU errors reported by the device (see [`gpu::GpuErrors`]), shown in the status bar.
    pub gpu_errors: gpu::GpuErrors,
    /// The overlays drawn over the viewport last frame (legend, selection bar, gizmo, navigation
    /// bar): a fit leaves them free ([`framing::free_part`]).
    pub overlays: Vec<egui::Rect>,
    /// What the current drag does to the camera, chosen when it starts.
    pub drag: Option<NavMode>,
    /// The point the current orbit turns about (drawn as a marker while orbiting).
    pub orbit_pivot: Option<Vec3>,
    /// The item the open context menu was opened on (a datum or a partition), if any.
    pub menu_target: Option<menu::Target>,
}

/// The XY plot in the results panel (results.xy): the quantity picked, the node it is read at,
/// and the series the command replied with for `key` (step, result identity, `y`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct XyPlot {
    /// `results.xy`'s `y`; a nodal field's `node` is filled from [`XyPlot::node`].
    pub y: Option<Value>,
    /// `results.xy`'s `x` when it is a travel (`None`: the step's own axis).
    pub x: Option<Value>,
    pub node: Option<u32>,
    pub key: Option<(String, usize, Value)>,
    pub series: Option<Value>,
}

/// A probed node: its values as (what, value) rows and as lines of text, for the step, mode
/// and frame they were read from.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Probe {
    pub node: u32,
    pub step: String,
    pub mode: usize,
    /// The explicit field frame (`None`: the last).
    pub frame: Option<usize>,
    pub rows: Vec<(String, String)>,
    pub lines: Vec<String>,
}

/// What a results plot was built for.
#[derive(Clone, Debug, PartialEq)]
struct PlotKey {
    revision: u64,
    /// The result's identity: every solve stores a new one.
    result: usize,
    current: bool,
    view: fieldcraft_engine::results::ResultsView,
    /// The mode and field frame drawn (an animated explicit step changes frames).
    pick: fieldcraft_engine::results::Pick,
    frame: Option<u64>,
    /// The display units its numbers (legend, scale) are in.
    units: String,
}

/// A datum form's datum and what it was worked out for (revision, command, form text).
type DatumPreview = (Option<(u64, &'static str, String)>, Option<fieldcraft_model::Datum>);
/// The nodes the open Merge Nodes form would merge (both of each pair), and the revision and form
/// text they are for.
type MergePreview = (Option<(u64, String)>, Arc<Vec<u32>>);
/// The open Relax Nodes form's preview, by model revision and form text.
type RelaxPreview = (Option<(u64, String)>, Option<Arc<fieldcraft_engine::cmd::RelaxPreview>>);
/// The open Solve Fluid form's memory estimate in a line and why the run would be refused (or
/// why there is no estimate), by model revision and form text.
type FluidNote = (Option<(u64, String)>, Option<(String, Option<String>)>);

/// A screenshot asked for over the control channel.
struct Shot {
    token: u64,
    path: Option<String>,
    reply: Sender<Value>,
    /// Frames waited so far: the request's own changes get drawn before the picture is taken.
    frames: u32,
    requested: bool,
    /// App time after which the request is answered with an error.
    deadline: f64,
}

pub struct ModelApp {
    pub session: Session,
    pub ui: UiState,
    pub cam: Camera,
    pub status: Option<Status>,
    pub dialog: Option<ParamsDialog>,
    /// The datum the open form defines, and the revision, command and form text it is for.
    preview: std::cell::RefCell<DatumPreview>,
    merge_preview: std::cell::RefCell<MergePreview>,
    relax_preview: std::cell::RefCell<RelaxPreview>,
    fluid_note: std::cell::RefCell<FluidNote>,
    /// The open Move Nodes form's move, worked out for the view ([`node_move`]).
    node_move: std::cell::RefCell<node_move::Cache>,
    /// The Move Nodes drag in progress.
    pub node_drag: Option<node_move::Drag>,
    /// Forms opened so far (numbers each [`ParamsDialog::opening`]).
    forms_opened: u64,
    pub viewport: Viewport,
    pub icons: icons::IconCache,
    /// The toolbar's buttons drawn last frame, by catalogue id (for the control channel).
    pub toolbar_rects: Vec<(&'static str, egui::Rect, toolbar::Fit)>,
    /// The GPU adapter, for the status bar and diagnostics.
    pub gpu_info: Option<String>,
    /// Draws `ui.render`'s pictures on the window's device with the viewport's pipelines; the
    /// CPU rasterizer draws them without it (no GPU, the browser, or a failed set-up, kept in
    /// [`App::render_png`]'s report).
    #[cfg(not(target_arch = "wasm32"))]
    capturer: Result<Option<fieldcraft_gpu::capture::Capturer>, String>,
    /// Bumped whenever the model may have changed.
    revision: u64,
    display: (Option<u64>, DisplayMesh),
    /// The analysis setup's symbols, built with the display data.
    glyphs: Glyphs,
    /// The results plot and what it was built for.
    plot: (Option<PlotKey>, Option<Arc<results::Plot>>),
    /// The node probed last (a click on a results plot).
    pub probe: Option<Probe>,
    /// The display units the probe, the reaction sums and the XY series were read in.
    shown_units: String,
    /// The reaction sums asked for last in the results panel.
    pub reactions: Option<Value>,
    /// The XY plot, while it is open in the results panel.
    pub xy: Option<XyPlot>,
    /// The last command reply with more to it than its status line, shown in the context panel.
    pub reply: Option<reply::Reply>,
    /// The Feature Graph rows opened to their details: (body, feature id).
    pub features_open: std::collections::BTreeSet<(String, u64)>,
    /// Whether results were shown last frame (the Results tab comes up when they appear).
    results_shown: bool,
    /// The framing the viewport does when it next draws ([`framing::Reframe`]).
    reframe: Option<framing::Reframe>,
    /// The session's camera requests seen ([`Session::camera_requests`]: a saved view shown).
    camera_seen: u64,
    /// The open document ([`Session::opened`]) the view was framed for.
    framed_document: Option<u64>,
    /// Tessellated design bodies with their names, kept while the same design is open.
    bodies: (Option<Arc<fieldcraft_geometry::Design>>, BodyMeshes),
    /// The display data a finished job built for the session it brought, adopted by the next
    /// [`Self::display`]; any later change to the model drops it.
    prebuilt: Option<Prebuilt>,
    /// The frame-independent part of the results plot (scale, contour) and what it was built for.
    prepared: (Option<PlotKey>, results::Prepared),
    /// The built scene, its key and the box around its datums (they reach beyond the model).
    scene: (Option<u64>, Scene, Aabb3),
    /// The setup's symbols, a layer of their own drawn at the camera's scale, its key and what
    /// it drew.
    symbols: (Option<u64>, Scene, Drawn),
    /// The symbols of the step a results plot shows, and the model revision and step they are of.
    step_glyphs: (Option<(u64, String)>, Glyphs),
    /// The theme last applied to the context.
    styled: Option<bool>,
    /// Input injected by the control channel, fed to egui one pointer event per frame.
    pub synthetic: Vec<egui::Event>,
    synthetic_mods: egui::Modifiers,
    /// Files dropped on the window by the control channel (`ui.drop`): the next frame takes them
    /// with the window system's drops.
    pub synthetic_drops: Vec<String>,
    control_rx: Option<inbox::Inbox<ControlRequest>>,
    shots: Vec<Shot>,
    shot_token: u64,
    /// The offscreen context `ui.capture` draws the window with (made on first use).
    offscreen: Option<Box<capture::Offscreen>>,
    /// The app time (seconds) of the current frame.
    pub now: f64,
    pub quit_requested: bool,
    /// Host services: native file pickers, the recent-files file.
    pub services: files::Services,
    /// Recently opened and saved files.
    pub recent: files::Recent,
    /// The user's preferences as last saved.
    prefs: files::Prefs,
    /// The host's go-ahead to share usage, when this run may ([`privacy`]).
    usage: Option<privacy::UsageHost>,
    /// The consent question was asked in this run.
    consent_asked: bool,
    /// Checking for and installing updates, when this run may ([`updates`]).
    updates: Option<updates::State>,
    /// The web app: a new deployment is out, and this reloads the page into it.
    reload: Option<Box<dyn Fn()>>,
    /// The setup progress the browser marks, and what it was worked out for.
    setup: (Option<(u64, u64, bool, bool)>, std::sync::Arc<Vec<fieldcraft_engine::progress::StageState>>),
    /// The file opening in the background, and files waiting to open after it.
    pub loading: Option<files::Loading>,
    /// The last command refused for what the model misses, listed in the context panel until
    /// it can run or is dismissed.
    pub refused: Option<fieldcraft_engine::requirement::Unavailable>,
    /// Jobs started over the control channel (numbers each).
    pub jobs: u64,
    pub queued: std::collections::VecDeque<String>,
    /// The question the window is asking (unsaved changes, recover).
    pub prompt: Option<files::Prompt>,
    /// Help › Report a Bug's form, while open ([`files::Prompt::ReportBug`]).
    pub report: Option<report::ReportForm>,
    /// The last error the status bar showed (what Report a Bug starts from).
    pub last_error: Option<String>,
    /// A design opened over a model of another part, waiting for [`files::Prompt::Design`]'s answer.
    pub held: Option<files::Opened>,
    /// An imported deck's document identity: unsaved model data has no command history.
    imported_document: Option<u64>,
    /// The window draws its own title bar and caption buttons (Windows, Linux).
    pub custom_titlebar: bool,
    /// The title bar runs under the macOS traffic lights.
    pub integrated_titlebar: bool,
    /// The command search's query and the command line's text.
    pub palette_text: String,
    pub command_line: String,
    /// Browser sections that are unfolded.
    pub tree_open: std::collections::BTreeSet<String>,
    /// The definition selected in the browser: (kind, name), its properties on the right.
    pub selected: Option<(browser::Def, String)>,
    /// The edge seed selected in the browser (by its edge), its properties on the right; never
    /// together with [`ModelApp::selected`].
    pub selected_seed: Option<fieldcraft_model::GeometryRef>,
    /// The Mesh node is selected in the browser: how the current mesh was made shows on the
    /// right ([`ModelApp::mesh_history`]); never together with a definition or a seed.
    pub selected_mesh: bool,
    /// What is picked in the viewport.
    pub pick: PickState,
    /// Set while [`ModelApp::capture`] draws the window offscreen: no pointer there is the
    /// user's, so the hover stays as it is.
    capturing: bool,
    /// What a left drag in the viewport does.
    pub nav: NavMode,
    /// Where the window's regions were drawn last frame (`ui.inspect`'s `layout`).
    pub layout: std::collections::BTreeMap<&'static str, egui::Rect>,
    /// The toolbar's panel labels drawn last frame, by `<tab>/<panel>`.
    pub panel_rects: Vec<(String, egui::Rect)>,
    /// The timeline's items and the journal state they were built for.
    timeline: Option<((u64, usize, (usize, usize)), std::sync::Arc<Vec<timeline::Item>>)>,
    /// The browser's element counts and the model revision they were counted at.
    element_counts: Option<(u64, Vec<(fieldcraft_model::ElementKind, usize)>)>,
    /// The selected definition's property rows and what they were read for.
    /// The selected definition's rows, by (revision, display units, kind, name).
    properties: Option<((u64, String, browser::Def, String), std::sync::Arc<Vec<(String, String)>>)>,
    /// [`ModelApp::mesh_history`] with each command's rows, by model revision.
    mesh_made: Option<((u64, String), std::sync::Arc<Vec<(String, Value, Vec<(String, String)>)>>)>,
    /// The search's query and its results.
    search: Option<(String, std::sync::Arc<Vec<palette::Hit>>)>,
}

/// Register the engine's picture drawer ([`fieldcraft_engine::view::register_drawer`]) for the
/// pictures commands draw themselves (`results.export.animation`): the viewport's GPU pipeline on
/// the process's shared headless device under `FIELDCRAFT_GPU`, the CPU rasterizer without one.
#[cfg(not(target_arch = "wasm32"))]
pub fn register_picture_drawer() {
    fieldcraft_engine::view::register_drawer(Box::new(|f, w, h| fieldcraft_gpu::capture::rgba_shared(&f.scene, &[&f.symbols], &f.cam, w, h)));
}

impl ModelApp {
    pub fn new(mut session: Session) -> ModelApp {
        // The window's own model: its commands count as clicks unless said otherwise.
        session.origin.get_or_insert(fieldcraft_engine::usage::Origin::Ui);
        ModelApp {
            session,
            ui: UiState::default(),
            cam: Camera::default(),
            status: None,
            dialog: None,
            preview: Default::default(),
            merge_preview: Default::default(),
            relax_preview: Default::default(),
            fluid_note: Default::default(),
            node_move: Default::default(),
            node_drag: None,
            forms_opened: 0,
            viewport: Viewport::default(),
            icons: icons::IconCache::default(),
            toolbar_rects: Vec::new(),
            gpu_info: None,
            #[cfg(not(target_arch = "wasm32"))]
            capturer: Ok(None),
            revision: 0,
            display: (None, DisplayMesh::default()),
            glyphs: Glyphs::default(),
            plot: (None, None),
            probe: None,
            shown_units: String::new(),
            reactions: None,
            xy: None,
            reply: None,
            features_open: Default::default(),
            results_shown: false,
            reframe: None,
            camera_seen: 0,
            framed_document: None,
            bodies: (None, Vec::new()),
            prebuilt: None,
            prepared: (None, results::Prepared::default()),
            scene: (None, Scene::default(), Aabb3::EMPTY),
            symbols: (None, Scene::default(), Drawn::default()),
            step_glyphs: (None, Glyphs::default()),
            styled: None,
            synthetic: Vec::new(),
            synthetic_mods: egui::Modifiers::default(),
            synthetic_drops: Vec::new(),
            control_rx: None,
            shots: Vec::new(),
            shot_token: 0,
            offscreen: None,
            now: 0.0,
            quit_requested: false,
            services: files::Services::default(),
            recent: files::Recent::default(),
            prefs: files::Prefs::default(),
            usage: None,
            consent_asked: false,
            updates: None,
            reload: None,
            setup: (None, std::sync::Arc::default()),
            loading: None,
            refused: None,
            jobs: 0,
            queued: std::collections::VecDeque::new(),
            prompt: None,
            report: None,
            last_error: None,
            held: None,
            imported_document: None,
            custom_titlebar: false,
            integrated_titlebar: false,
            palette_text: String::new(),
            command_line: String::new(),
            tree_open: ["mesh", "results"].into_iter().map(str::to_string).collect(),
            selected: None,
            selected_seed: None,
            selected_mesh: false,
            pick: PickState::default(),
            capturing: false,
            nav: NavMode::default(),
            layout: std::collections::BTreeMap::new(),
            panel_rects: Vec::new(),
            timeline: None,
            element_counts: None,
            search: None,
            properties: None,
            mesh_made: None,
        }
    }

    /// Serve control-channel requests (see [`control`]).
    pub fn with_control(mut self, rx: impl Into<inbox::Inbox<ControlRequest>>) -> Self {
        self.control_rx = Some(rx.into());
        self
    }

    /// Use the host's services (file pickers, the recent-files file) and load the recent list.
    pub fn with_services(mut self, services: files::Services) -> Self {
        self.recent = files::Recent::load(services.recent_path.clone());
        self.prefs = files::Prefs::load(services.prefs.as_ref());
        self.session.default_units = self.prefs.units.clone();
        self.ui.nav_preset = self.prefs.nav_preset.clone();
        self.ui.reverse_zoom = self.prefs.reverse_zoom;
        self.services = services;
        self
    }

    /// Draw the viewport on the GPU (eframe's wgpu render state, with the app's depth buffer
    /// bits and MSAA sample count). Without it the viewport renders on the CPU. Device errors
    /// show in the status bar.
    pub fn set_wgpu(&mut self, rs: &egui_wgpu::RenderState, depth_bits: u8, samples: u32, ctx: &egui::Context) {
        self.viewport.gpu = Some(gpu::install(rs, depth_bits, samples, self.viewport.gpu_errors.clone(), ctx.clone()));
        let info = rs.adapter.get_info();
        self.gpu_info = Some(format!("{} ({:?})", info.name.trim(), info.backend));
        #[cfg(not(target_arch = "wasm32"))]
        {
            let tier = fieldcraft_gpu::select::rank(info.device_type).0;
            self.capturer =
                fieldcraft_gpu::capture::Capturer::on_device(&rs.adapter, rs.device.clone(), rs.queue.clone(), tier).map(Some).map_err(|e| {
                    log::warn!("ui.render draws on the CPU: {e}");
                    e.to_string()
                });
        }
    }

    /// Run a command with JSON parameters: what every panel, the control channel and scripts
    /// do. The result or the error goes to the status bar.
    pub fn run(&mut self, id: &str, params: Value) -> Result<Value, String> {
        self.run_form(id, params, false).map_err(|f| f.message)
    }

    /// [`Self::run`], from a form when `form`: a refused parameter is then named by its field's
    /// label (the engine's message names its key, which agents need).
    #[allow(clippy::disallowed_methods, reason = "the window's quick commands; long ones (`files::is_long`) start as jobs")]
    fn run_form(&mut self, id: &str, params: Value, form: bool) -> Result<Value, Failure> {
        self.ensure_idle().map_err(|message| Failure { message, refused: None })?;
        let (revision, view) = (self.session.revision(), self.session.view.clone());
        let label = command_label(id);
        let r = self.session.execute(id, &params).map_err(|e| match e {
            EngineError::Param(p) => {
                let message = if form { p.for_form() } else { p.for_agents() };
                let shown = p.for_form();
                Failure { refused: Some((p.path, shown)), message }
            }
            e => Failure { message: e.to_string(), refused: None },
        });
        let changed = self.session.revision() != revision;
        // The definition selected in the browser keeps its selection under its new name.
        if id == "model.rename"
            && let Ok(v) = &r
            && let Some((def, name)) = &mut self.selected
            && def.rename_kind() == v["kind"].as_str()
            && v["from"].as_str() == Some(name.as_str())
            && let Some(to) = v["name"].as_str()
        {
            *name = to.to_string();
        }
        // A view change (hide, show) draws the model again, as an edit does.
        if changed || self.session.view != view {
            self.model_changed();
        }
        // The results display commands reply with what is shown, which the viewport shows.
        let display = (id.starts_with("results.") && id != "results.reactions") || id == "solve.run";
        // The reaction sums show in the results panel, however they were asked for.
        if id == "results.reactions"
            && let Ok(v) = &r
        {
            self.reactions = Some(v.clone());
            self.shown_units = self.session.units().label();
        }
        // The XY plot's series, however it was asked for; the panel opens on it.
        if id == "results.xy" && (r.is_ok() || params.get("path").is_none()) {
            let xy = self.xy.get_or_insert_with(XyPlot::default);
            xy.series = r.as_ref().ok().cloned();
            // The picker follows a series asked for elsewhere (the control channel, a script).
            if let Some(y) = params.get("y").filter(|_| r.is_ok()) {
                let mut template = y.clone();
                if let Some(node) = template.as_object_mut().and_then(|o| o.remove("node")) {
                    xy.node = node.as_u64().and_then(|n| u32::try_from(n).ok());
                }
                xy.y = Some(template);
                xy.x = params.get("x").filter(|x| x.is_object()).cloned();
            }
        }
        // What the duplicate check found is selected, so the viewport shows it (#234).
        if id == "mesh.check.duplicates"
            && let Ok(v) = &r
        {
            self.select_duplicates(v);
        }
        if r.is_ok() && id == "deck.read" {
            self.imported_document = Some(self.session.opened());
        }
        if r.is_ok() && (changed || display) {
            self.follow_stage(id);
        }
        match &r {
            Ok(v) if v.is_null() || (changed && id != "deck.read") || display => self.set_status(reply::done(id, label, &params, v), false),
            // A reply reads as one line; the whole of a rich one shows in the context panel (the
            // reaction sums have their own place in the results panel).
            Ok(v) => {
                let mut line = format!("{label}: {}", reply::summary_in(id, v, Some(self.session.units())));
                if reply::is_rich(v) && !(id == "results.reactions" && self.plot().is_some()) {
                    line.push_str(" (details in the panel)");
                    // Kept with their units: a live reply runs again in whatever units are shown then.
                    let params = self.session.params_with_units(id, &params);
                    self.reply = Some(reply::Reply { command: id.to_string(), label: label.to_string(), params, value: v.clone() });
                }
                self.set_status(line, false)
            }
            Err(e) => self.set_status(format!("{label}: {}", e.message), true),
        }
        r
    }

    /// Run a script (a JSON array of `{command, params}`), stopping at the first error. The
    /// ribbon follows its last command that has a stage.
    pub fn run_script(&mut self, script: &Value) -> Result<(), String> {
        // A script waits for the background job, as a command does: the job's session would
        // replace what it changed.
        self.ensure_idle()?;
        #[allow(clippy::disallowed_methods, reason = "scripts are not jobs yet")]
        let r = self.session.run_script(script).map(|_| ()).map_err(|e| e.to_string());
        self.model_changed();
        if r.is_ok() {
            let commands = script.as_array().into_iter().flatten().filter_map(|step| step.get("command").and_then(Value::as_str));
            if let Some(id) = commands.rev().find(|id| catalogue::stage_tab(id).is_some()) {
                self.follow_stage(id);
            }
        }
        match &r {
            Ok(()) => self.set_status("Script done", false),
            Err(e) => self.set_status(format!("Script: {e}"), true),
        }
        r
    }

    /// The ribbon turns to the tab where `id` lives ([`catalogue::stage_tab`]), as a person's click
    /// there would, so commands run by scripts, the control channel and MCP show their stage.
    /// Called after a command changed the model or the results display; queries leave the ribbon
    /// alone, and a ribbon click never moves it.
    fn follow_stage(&mut self, id: &str) {
        if let Some(tab) = catalogue::stage_tab(id) {
            self.ui.select_tab(tab);
        }
    }

    /// What a click on a toolbar entry, a menu row or a palette result does: a command without
    /// parameters runs, one with parameters opens its form, and the file entries (open, save,
    /// write a deck, recover) go through the host's file pickers. `id` is a command or a
    /// catalogue entry; an entry that is coming soon does nothing and says when it comes.
    /// Programmatic calls use [`ModelApp::run`] or [`ModelApp::start_programmatic`].
    pub fn start(&mut self, id: &str) -> Result<(), String> {
        self.begin(id, true)
    }

    /// `ui.start` from the control channel: as a click, but never a native file picker (entries
    /// that need one open their command's form).
    pub fn start_programmatic(&mut self, id: &str) -> Result<(), String> {
        self.begin(id, false)
    }

    fn begin(&mut self, id: &str, interactive: bool) -> Result<(), String> {
        let (command, preset) = match catalogue::find(id) {
            Some(c) => match c.command() {
                Some(command) => (command, c.preset()),
                None => {
                    let e = format!("{} is coming soon ({})", c.label, c.milestone().unwrap_or("later"));
                    self.set_status(e.clone(), false);
                    return Err(e);
                }
            },
            None => (id, json!({})),
        };
        let spec = fieldcraft_engine::find_command(command).ok_or_else(|| format!("unknown command `{command}`"))?;
        // A command missing what it needs opens no form: the status names what is missing.
        if let Err(u) = self.session.availability(spec) {
            self.record_ui(fieldcraft_engine::usage::Event::UnavailableClicked {
                command: u.command,
                missing: u.missing.iter().map(|m| m.0).collect(),
            });
            let e = u.to_string();
            self.refuse(u);
            return Err(e);
        }
        if self.file_action(spec.id, &preset, interactive) {
            self.dialog = None;
            return Ok(());
        }
        if spec.params.trim() == "{}" {
            self.dialog = None;
            self.run(spec.id, json!({})).map(|_| ())
        } else {
            // The form starts with the entry's preset; the rest (a path…) is filled in there.
            self.open_form(spec.id, preset);
            Ok(())
        }
    }

    /// Open a command's parameter form with `params` filled in.
    pub fn open_form(&mut self, command: &str, params: Value) {
        if let Some(spec) = fieldcraft_engine::find_command(command) {
            // A mesher starts from how the mesh it would replace was made (size, order…), so a
            // remesh changes only what is typed (#343).
            let params = match params.as_object() {
                Some(o) if o.is_empty() => self.last_mesher_params(spec.id).unwrap_or(params),
                _ => params,
            };
            // The schema's defaults are filled in, so the form starts as the command would run.
            let mut params = spec.form.map_or(params.clone(), |form| fieldcraft_engine::schema::prefill(form, &params));
            // A required reference with one obvious candidate starts picked.
            if let (Some(form), Some(obj)) = (spec.form, params.as_object_mut()) {
                let selected = self.selected.clone();
                let session = &self.session;
                let pick = |c| match &selected {
                    Some((kind, name)) if kind.collection() == c => Some(name.clone()),
                    _ => Some(session.names(c)).filter(|n| n.len() == 1).and_then(|n| n.into_iter().next()),
                };
                form::preselect(form, obj, &pick);
                // A new datum or partition starts with the next free name (Datum-1, …), editable.
                if let Some(name) = self.default_name(spec.id)
                    && obj.get("name").and_then(Value::as_str).is_none_or(str::is_empty)
                {
                    obj.insert("name".into(), json!(name));
                }
            }
            self.forms_opened = self.forms_opened.wrapping_add(1);
            self.dialog = Some(ParamsDialog {
                command: spec.id,
                text: params.to_string(),
                error: None,
                refused: None,
                touched: Default::default(),
                attempted: false,
                advanced: spec.form.is_none(),
                opening: self.forms_opened,
                back: None,
                pending_set: None,
            });
        }
    }

    /// The parameters mesher `id` last ran with in this model (its journal), while the mesh it
    /// made is there; `None` for other commands.
    fn last_mesher_params(&self, id: &str) -> Option<Value> {
        fieldcraft_engine::MESHERS.iter().find(|(m, _)| *m == id)?;
        let (_, params) = self.session.journal.iter().rev().find(|(c, _)| c == id)?;
        let targets = fieldcraft_engine::MESHERS.iter().find(|(m, _)| *m == id)?.1;
        let meshed = targets(&self.session, params).ok()?.iter().any(|b| !self.session.model.body_mesh(b).elements.is_empty());
        meshed.then(|| params.clone())
    }

    /// A command that cannot run yet: the status bar says so in short, and the context panel
    /// lists each missing piece with the command that provides it ([`ModelApp::refused`]).
    pub fn refuse(&mut self, u: fieldcraft_engine::requirement::Unavailable) {
        let label = command_label(u.command);
        let text = match u.missing.as_slice() {
            [one] => format!("{label} cannot run yet. {one}"),
            many => format!("{label} cannot run yet: {} things missing (listed in the panel)", many.len()),
        };
        self.set_status(text, true);
        self.refused = Some(u);
    }

    /// The name a new item of command `id` starts with: `<Kind>-<n>`, the first not taken.
    fn default_name(&self, id: &str) -> Option<String> {
        use fieldcraft_engine::schema::Collection;
        let (base, collection) = match id {
            "geometry.datum" | "geometry.datum.axis" | "geometry.datum.point" => ("Datum", Collection::Datum),
            "geometry.partition" => ("Partition", Collection::Partition),
            _ => return None,
        };
        let taken = self.session.names(collection);
        (1..=fieldcraft_engine::schema::MAX_NAMES).map(|n| format!("{base}-{n}")).find(|n| !taken.contains(n))
    }

    /// A reference field's "Create …": open the command that makes an item of `collection` over
    /// the open form, which comes back (with the new item at `path`) when that command runs or
    /// is cancelled. A provider that opens no form (a file picker, a missing requirement) leaves
    /// the form as it was.
    pub fn create_from_form(&mut self, collection: fieldcraft_engine::schema::Collection, path: String) {
        let Some(parent) = self.dialog.take() else { return };
        let before = self.session.names(collection);
        let provider = collection.requirement().provider();
        let _ = self.start(provider);
        match self.dialog.as_mut() {
            Some(d) if d.command == provider && d.back.is_none() => d.back = Some(Box::new(Back { dialog: parent, collection, path, before })),
            _ => self.dialog = Some(parent),
        }
    }

    /// Close the open form: back to the form it was opened from, if any.
    pub fn close_dialog(&mut self) {
        if let Some(d) = &self.dialog {
            self.record_ui(fieldcraft_engine::usage::Event::FormCancelled { command: d.command });
        }
        self.dialog = self.dialog.take().and_then(|d| d.back).map(|b| b.dialog);
    }

    /// The open form's command ran: back to the form it was opened from, with the item it made
    /// in the reference field.
    fn dialog_done(&mut self) {
        let Some(back) = self.dialog.take().and_then(|d| d.back) else { return };
        let Back { mut dialog, collection, path, before } = *back;
        let made = self.session.names(collection).into_iter().find(|n| !before.contains(n));
        if let (Some(name), Ok(mut params)) = (made, serde_json::from_str::<Value>(&dialog.text)) {
            let (parent, key) = path.rsplit_once('/').unwrap_or(("", &path));
            let pointer = if parent.is_empty() { String::new() } else { format!("/{parent}") };
            if let Some(Value::Object(obj)) = params.pointer_mut(&pointer) {
                obj.insert(key.into(), json!(name));
            }
            dialog.text = params.to_string();
        }
        self.dialog = Some(dialog);
    }

    /// The open form's command succeeded (here or as a job): a definition it made is the
    /// selection, its properties in the panel, and the form closes (back to the form it was
    /// opened from, if any).
    pub(crate) fn form_done(&mut self) {
        if let Some(d) = &self.dialog
            && let Some(kind) = browser::Def::ALL.into_iter().find(|k| k.command("create") == d.command)
            && let Ok(params) = serde_json::from_str::<Value>(&d.text)
            && let Some(name) = params.get("name").and_then(Value::as_str)
            && kind.exists(&self.session.model, name)
        {
            self.selected = Some((kind, name.to_string()));
        }
        self.dialog_done();
    }

    /// Run the open parameter form; it closes when the command succeeds.
    #[allow(clippy::disallowed_methods, reason = "creates the form's new set, and takes it back when the form's command fails")]
    pub fn run_dialog(&mut self) -> Result<Value, String> {
        let (command, params) = self.dialog_params()?;
        // The set Use selection made, when the form still names it: created first, and taken back
        // if the command then fails.
        let pending = self.dialog.as_ref().and_then(|d| d.pending_set.clone()).filter(|set| {
            let name = set.get("name").and_then(Value::as_str);
            name.is_some() && [params.pointer("/target/set"), params.get("set")].into_iter().flatten().any(|v| v.as_str() == name)
        });
        if let Some(set) = &pending
            && let Err(e) = self.session.execute("model.set.create", set)
        {
            let message = format!("the selection's set: {e}");
            if let Some(open) = self.dialog.as_mut() {
                open.error = Some(message.clone());
            }
            return Err(message);
        }
        match self.run_form(command, params, true) {
            Ok(v) => {
                self.form_done();
                Ok(v)
            }
            Err(Failure { message, refused }) => {
                if pending.is_some() {
                    let _ = self.session.execute("edit.undo", &json!({}));
                    self.model_changed();
                }
                if let Some(open) = self.dialog.as_mut() {
                    open.error = Some(message.clone());
                    open.refused = refused;
                }
                Err(message)
            }
        }
    }

    /// The open form's command and parameters; parameters that are not JSON are the form's error.
    pub(crate) fn dialog_params(&mut self) -> Result<(&'static str, Value), String> {
        let Some(d) = self.dialog.as_mut() else { return Err("no command form is open".into()) };
        match serde_json::from_str::<Value>(&d.text) {
            Ok(params) => {
                d.error = None;
                Ok((d.command, params))
            }
            Err(e) => {
                let e = format!("the parameters are not JSON: {e}");
                d.error = Some(e.clone());
                Err(e)
            }
        }
    }

    /// Probe a node of the results shown (`solve.probe`) and keep its values for the viewport.
    pub fn probe_node(&mut self, node: u32) -> Result<(), String> {
        let Some(view) = self.session.results_view.clone() else { return Err("no results are shown".into()) };
        let pick = self.plot().map_or(view.pick(), |p| p.pick);
        let mut params = json!({"step": view.step, "nodes": [node], "mode": pick.mode});
        if let Some(frame) = pick.frame {
            params["frame"] = json!(frame);
        }
        let reply = match self.run("solve.probe", params) {
            Ok(v) => v,
            Err(e) => {
                self.probe = None;
                return Err(e);
            }
        };
        // The probe shows in the results panel and beside the node, not as a reply.
        if self.reply.as_ref().is_some_and(|r| r.command == "solve.probe") {
            self.reply = None;
        }
        let probe = self.read_probe(&reply, node, view.step, pick);
        self.set_status(probe.lines.join("   "), false);
        self.probe = Some(probe);
        self.shown_units = self.session.units().label();
        Ok(())
    }

    /// The probe of `node` from `solve.probe`'s reply, in the display units.
    fn read_probe(&mut self, reply: &Value, node: u32, step: String, pick: fieldcraft_engine::results::Pick) -> Probe {
        let plot = self.plot();
        let rows = probe_rows(reply, plot.as_deref(), self.session.units(), node, pick.mode);
        let lines: Vec<String> = std::iter::once(format!("Node {node}")).chain(rows.iter().map(|(what, value)| format!("{what} {value}"))).collect();
        Probe { node, step, mode: pick.mode, frame: pick.frame, rows, lines }
    }

    pub fn set_status(&mut self, text: impl Into<String>, error: bool) {
        let text = text.into();
        if error {
            log::warn!("{text}");
            self.last_error = Some(text.clone());
        }
        self.status = Some(Status { text, error });
    }

    /// Move an error the graphics device reported into the status bar.
    fn show_gpu_errors(&mut self) {
        let error = self.viewport.gpu_errors.lock().ok().and_then(|mut e| e.take());
        if let Some(e) = error {
            self.status = Some(Status { text: format!("Graphics: {e}"), error: true });
        }
    }

    /// The model may have changed: display data and the scene are rebuilt when next drawn.
    pub fn model_changed(&mut self) {
        self.revision = self.revision.wrapping_add(1);
        // Display data built for the model before this change is out of date (and may be large).
        if let Some(stale) = self.prebuilt.take() {
            files::off_thread("fieldcraft drop".into(), move || drop(stale));
        }
        self.frame_new_document();
        self.follow_units();
        // A reply about the model is read again; any other one is out of date.
        self.reply = self.reply.take().and_then(|mut r| {
            reply::LIVE.contains(&r.command.as_str()).then_some(())?;
            #[allow(clippy::disallowed_methods, reason = "reads the model again for its reply (`reply::LIVE`: listings and checks)")]
            let value = self.session.execute(&r.command, &r.params);
            r.value = value.ok()?;
            Some(r)
        });
    }

    /// Read the probe, the reaction sums and the XY series again when the display units changed
    /// (`model.units`, Settings › Units): they are kept as shown, in the units they were read in.
    #[allow(clippy::disallowed_methods, reason = "reads one probed node and one step's reaction sums again")]
    fn follow_units(&mut self) {
        let units = self.session.units().label();
        if units == self.shown_units {
            return;
        }
        self.shown_units = units;
        if let Some(old) = self.probe.take() {
            let mut params = json!({"step": old.step, "nodes": [old.node], "mode": old.mode});
            if let Some(frame) = old.frame {
                params["frame"] = json!(frame);
            }
            // A failure (a stale result) leaves no probe, as a click would.
            if let Ok(reply) = self.session.execute("solve.probe", &params) {
                let pick = fieldcraft_engine::results::Pick { mode: old.mode, frame: old.frame };
                self.probe = Some(self.read_probe(&reply, old.node, old.step, pick));
            }
        }
        if let Some(step) = self.reactions.as_ref().and_then(|r| r["step"].as_str()).map(str::to_string) {
            self.reactions = self.session.execute("results.reactions", &json!({"step": step})).ok();
        }
        if let Some(xy) = self.xy.as_mut() {
            (xy.key, xy.series) = (None, None);
        }
    }

    /// The model's display data with its geometry (the open design's bodies, tessellated once
    /// per design, and the imported surfaces), rebuilt when the model changed; what a finished
    /// job built off the window thread ([`Prebuilt`]) is taken as it is. The first content to
    /// appear is framed.
    pub fn display(&mut self) -> &DisplayMesh {
        if self.display.0 != Some(self.revision) {
            let prebuilt = self.prebuilt.take().filter(|p| p.revision == self.session.revision());
            let (mut d, prebodies) = match prebuilt {
                Some(p) => (p.model, p.bodies),
                None => (draw(&self.session), None),
            };
            let design = self.session.design.clone();
            if !same_design(&design, &self.bodies.0) {
                let (meshes, failed) = prebodies.unwrap_or_else(|| design.as_deref().map(tessellate).unwrap_or_default());
                self.bodies = (design, meshes);
                if !failed.is_empty() {
                    self.set_status(format!("Some bodies cannot be drawn: {}", failed.join("; ")), true);
                }
            }
            self.pick.names.clear();
            // Hidden bodies and surfaces are neither drawn nor picked.
            for (names, mesh) in self.bodies.1.iter().filter(|(n, _)| self.session.shows(fieldcraft_engine::ItemKind::Body, &n.name)) {
                d.add_curves(&mesh.edges);
                d.add_surface(&mesh.positions, &mesh.triangles, &mesh.tri_face);
                self.pick.names.push(names.clone());
            }
            // An imported surface's faces are named by their tags.
            for s in self.session.surfaces.iter().filter(|s| self.session.shows(fieldcraft_engine::ItemKind::Surface, &s.name)) {
                let positions: Vec<[f64; 3]> = s.mesh.positions.iter().map(|p| fieldcraft_engine::units::geometry_point(*p)).collect();
                d.add_surface(&positions, &s.mesh.triangles, &s.mesh.tri_face);
                self.pick.names.push(selection::BodyNames {
                    name: s.name.clone(),
                    faces: (0..s.faces).map(|f| f.to_string()).collect(),
                    edges: Vec::new(),
                    vertices: Vec::new(),
                });
            }
            let was_empty = self.display.1.bounds.is_empty();
            self.display = (Some(self.revision), d);
            self.glyphs = self.session.glyphs();
            if was_empty && !self.display.1.bounds.is_empty() {
                self.reframe = framing::Reframe::merge(self.reframe, framing::Reframe::Fit);
            }
        }
        &self.display.1
    }

    /// Definitions whose symbols are a sample (`Clamp: 400 of 12000 nodes`), as drawn last.
    pub fn glyph_notes(&self) -> &[String] {
        &self.symbols.2.notes
    }

    /// The key of the current scene: the model revision, every setting the scene depends on and
    /// the results plot shown.
    pub fn scene_key(&self) -> u64 {
        let mut h = std::hash::DefaultHasher::new();
        let plot = self.plot.1.as_ref().map(|p| Arc::as_ptr(p) as usize);
        let marks = fieldcraft_engine::view::geometry_marks(&self.session).map(|m| m as *const _ as usize);
        (self.revision, self.ui.style, self.ui.show(), self.ui.dark, plot, self.selected_datum(), marks).hash(&mut h);
        format!("{:?}", self.datum_preview()).hash(&mut h);
        h.finish()
    }

    /// What the results plot is built for: the model, the result, its view and the animation
    /// frame; `None` when no results are shown.
    fn plot_key(&self) -> Option<PlotKey> {
        let (view, result, current) = results::shown(&self.session)?;
        let frame = results::frame(view, result, self.now);
        let pick = results::pick(view, result, frame);
        Some(PlotKey {
            revision: self.revision,
            result: Arc::as_ptr(result) as usize,
            current,
            view: view.clone(),
            pick,
            frame,
            units: self.session.units().label(),
        })
    }

    /// The results plot the viewport draws, rebuilt when its key changes; `None` shows the
    /// model. Only the deformation follows the animation frame; the scale and contour values are
    /// kept while the rest of the key (the mode and field frame drawn among it) stays.
    pub fn plot(&mut self) -> Option<Arc<results::Plot>> {
        let key = self.plot_key();
        if self.plot.0 != key {
            self.display();
            let still = key.as_ref().map(|k| PlotKey { frame: None, ..k.clone() });
            if self.prepared.0 != still {
                let pick = key.as_ref().map(|k| k.pick).unwrap_or_default();
                let bodies: Vec<&str> = self.pick.names.iter().map(|n| n.name.as_str()).collect();
                self.prepared = (still, results::prepare(&self.session, &self.display.1, pick, &bodies).unwrap_or_default());
            }
            let frame = key.as_ref().and_then(|k| k.frame);
            let p = key.as_ref().and_then(|_| results::plot(&self.session, &self.display.1, frame, &self.prepared.1));
            self.plot = (key, p.map(Arc::new));
        }
        self.plot.1.clone()
    }

    /// The scene the viewport draws, rebuilt when its key changes; its radius fits the camera.
    pub fn scene(&mut self) -> (u64, &Scene) {
        self.display();
        self.plot();
        let key = self.scene_key();
        if self.scene.0 != Some(key) {
            let t = theme::Tokens::get(self.ui.dark);
            let (style, show) = (self.ui.style, self.ui.show());
            let mut s = match &self.plot.1 {
                Some(p) => results::scene(p, style, show, &t),
                None => scene::build(&self.display.1, style, show, &t),
            };
            // The elements Mark Bad Elements marked, in the error colour over the element faces.
            if self.plot.1.is_none() && show.mesh && style != scene::Style::Wireframe {
                let marked = fieldcraft_engine::view::marked(&self.session);
                fieldcraft_engine::view::mark(&mut s, &self.display.1, &marked, theme::rgb(t.error));
            }
            // Faces merged for meshing, each one face in the warm colour where geometry faces are
            // drawn; then what Check Geometry found, in the error colour: its faces, short edges
            // and a dot at each finding over everything.
            let faces = style != scene::Style::Wireframe && !(show.mesh && !self.display.1.tris.is_empty());
            let bodies: Vec<&str> = self.pick.names.iter().map(|n| n.name.as_str()).collect();
            if show.geometry && self.plot.1.is_none() && faces {
                fieldcraft_engine::view::draw_merged_faces(&mut s, &self.display.1, &self.session, &bodies, theme::rgb(t.warm));
            }
            if let Some(marks) = fieldcraft_engine::view::geometry_marks(&self.session).filter(|_| show.geometry && self.plot.1.is_none()) {
                fieldcraft_engine::view::draw_geometry_marks(&mut s, &self.display.1, &bodies, marks, faces, theme::rgb(t.error));
            }
            let mut datums = Aabb3::EMPTY;
            if show.geometry && self.plot.1.is_none() {
                datums = scene::datums(&mut s, &self.session, &self.display.1.bounds, self.selected_datum().as_deref(), &t);
                if let Some(d) = self.datum_preview() {
                    datums.union(scene::datum_preview(&mut s, d, &self.display.1.bounds, &t));
                }
                scene::seeds(&mut s, &self.session.model, &self.pick.names, &self.display.1, &t);
                scene::sizes(&mut s, &self.session.model, &self.pick.names, &self.display.1, &t);
            }
            self.scene = (Some(key), s, datums);
        }
        // The clipping radius holds the symbols and the datums too (framing leaves them out: the
        // symbols' size in the model follows the zoom, the datums are drawn at the model).
        let mut bounds = self.drawn_bounds();
        bounds.union(self.symbols.2.bounds);
        bounds.union(self.scene.2);
        self.scene.1.radius = scene::radius(&bounds, &self.cam);
        (key, &self.scene.1)
    }

    /// The datum drawn highlighted, shown or hidden: the one the open Partition form cuts by,
    /// else the one selected in the browser.
    fn selected_datum(&self) -> Option<String> {
        if let Some(d) = self.dialog.as_ref().filter(|d| d.command == "geometry.partition") {
            return serde_json::from_str::<Value>(&d.text).ok()?.get("datum").and_then(Value::as_str).map(str::to_string);
        }
        match &self.selected {
            Some((browser::Def::Datum, name)) => Some(name.clone()),
            _ => None,
        }
    }

    /// The bodies the open Partition form would cut: those it lists, or every drawn body whose
    /// surface lies on both sides of its datum plane. Empty without such a form.
    pub fn partition_preview(&self) -> Vec<String> {
        let Some(d) = self.dialog.as_ref().filter(|d| d.command == "geometry.partition") else { return Vec::new() };
        let Ok(params) = serde_json::from_str::<Value>(&d.text) else { return Vec::new() };
        if let Some(listed) = params.get("bodies").and_then(Value::as_array).filter(|l| !l.is_empty()) {
            return listed.iter().filter_map(Value::as_str).map(str::to_string).collect();
        }
        let Some((point, normal)) = params.get("datum").and_then(Value::as_str).and_then(|n| self.session.model.datum(n)).and_then(|d| d.plane())
        else {
            return Vec::new();
        };
        let (point, normal) = (Vec3::from(point), Vec3::from(normal));
        let display = &self.display.1;
        let tol = SECTION_SIDE * display.bounds.diagonal().max(MIN_SECTION_DIAGONAL);
        // Per body: whether it has surface points below and above the plane.
        let mut sides = vec![(false, false); self.pick.names.len()];
        for t in &display.surface_tris {
            let Some(side) = sides.get_mut(t.body as usize) else { continue };
            for p in t.p {
                let h = (p - point).dot(normal);
                side.0 |= h < -tol;
                side.1 |= h > tol;
            }
        }
        self.pick.names.iter().zip(sides).filter(|(_, s)| s.0 && s.1).map(|(b, _)| b.name.clone()).collect()
    }

    /// The nodes the open Merge Nodes form would merge, both of each pair, drawn before OK.
    /// Worked out once per form text and model revision (the search covers the whole mesh).
    pub fn merge_preview(&self) -> Arc<Vec<u32>> {
        let Some(d) = self.dialog.as_ref().filter(|d| d.command == "mesh.node.merge") else { return Arc::default() };
        let key = (self.revision, d.text.clone());
        let mut cache = self.merge_preview.borrow_mut();
        if cache.0.as_ref() != Some(&key) {
            let pairs = serde_json::from_str::<Value>(&d.text).ok().and_then(|p| fieldcraft_engine::cmd::node_merge_preview(&self.session, &p));
            let mut nodes: Vec<u32> = pairs.unwrap_or_default().into_iter().flatten().collect();
            nodes.sort_unstable();
            nodes.dedup();
            *cache = (Some(key), Arc::new(nodes));
        }
        cache.1.clone()
    }

    /// What the open Relax Nodes form would do (the nodes' new positions and the region's
    /// quality), worked out once per form text and model revision.
    pub fn relax_preview(&self) -> Option<Arc<fieldcraft_engine::cmd::RelaxPreview>> {
        let d = self.dialog.as_ref().filter(|d| d.command == "mesh.node.relax")?;
        let key = (self.revision, d.text.clone());
        let mut cache = self.relax_preview.borrow_mut();
        if cache.0.as_ref() != Some(&key) {
            let preview = serde_json::from_str::<Value>(&d.text).ok().and_then(|p| fieldcraft_engine::cmd::relax_preview(&self.session, &p));
            *cache = (Some(key), preview.map(Arc::new));
        }
        cache.1.clone()
    }

    /// The open Solve Fluid form's memory estimate in a line, and why the run would be refused
    /// before it starts (or why there is no estimate), worked out once per form text and model
    /// revision.
    pub fn fluid_note(&self) -> Option<(String, Option<String>)> {
        let d = self.dialog.as_ref().filter(|d| d.command == "solve.fluid")?;
        let mut cache = self.fluid_note.borrow_mut();
        if !cache.0.as_ref().is_some_and(|(revision, text)| *revision == self.revision && *text == d.text) {
            let note = match serde_json::from_str::<Value>(&d.text)
                .map_err(|e| e.to_string())
                .and_then(|p| fieldcraft_engine::cmd::fluid_estimate(&self.session, &p).map_err(|e| e.to_string()))
            {
                Ok(e) => (e.summary(), e.problem),
                Err(e) => ("Memory: no estimate".into(), Some(e)),
            };
            *cache = (Some((self.revision, d.text.clone())), Some(note));
        }
        cache.1.clone()
    }

    /// The open Relax Nodes form's preview in a line: nodes moved and the region's quality.
    pub fn relax_note(&self) -> Option<String> {
        let p = self.relax_preview()?;
        let r = &p.relaxed;
        Some(format!(
            "Preview: {} nodes move (largest {:.3}); scaled Jacobian min {:.3} to {:.3}, edge ratio max {:.2} to {:.2}, invalid {} to {}",
            r.moved,
            r.max_displacement,
            r.before.min_quality,
            r.after.min_quality,
            r.before.max_edge_ratio,
            r.after.max_edge_ratio,
            r.before.invalid,
            r.after.invalid
        ))
    }

    /// Whether the Relax Nodes form is open.
    pub fn relaxing(&self) -> bool {
        self.dialog.as_ref().is_some_and(|d| d.command == "mesh.node.relax")
    }

    /// Whether the open Relax Nodes form relaxes listed faces: clicks in the view then pick faces
    /// instead of placing a sphere.
    pub fn relax_by_faces(&self) -> bool {
        self.dialog
            .as_ref()
            .filter(|d| d.command == "mesh.node.relax")
            .and_then(|d| serde_json::from_str::<Value>(&d.text).ok())
            .is_some_and(|p| p.get("faces").and_then(Value::as_array).is_some_and(|f| !f.is_empty()))
    }

    /// A length of the scene (metres) in the display units, as a form takes it.
    pub(crate) fn shown_length(&self, m: f64) -> f64 {
        fieldcraft_engine::results::shown(fieldcraft_model::units::Quantity::Length, m, self.session.units())
    }

    /// A point of the scene (metres) in the display units, as a form takes it.
    pub(crate) fn shown_point(&self, p: Vec3) -> [f64; 3] {
        [p.x, p.y, p.z].map(|c| self.shown_length(c))
    }

    /// A point given in the display units (an agent's) as a point of the scene, in metres.
    pub(crate) fn scene_point(&self, p: [f64; 3]) -> Vec3 {
        let [x, y, z] = p.map(|c| self.session.units().to_si(fieldcraft_model::units::Quantity::Length, c).unwrap_or(f64::NAN));
        Vec3::new(x, y, z)
    }

    /// Put the open Relax Nodes form's sphere at `at` (a click on the mesh); a sphere without a
    /// radius yet gets a twentieth of the mesh's size.
    pub fn relax_place(&mut self, at: Vec3) -> bool {
        let mut bounds = Aabb3::EMPTY;
        for n in self.session.model.nodes().values() {
            bounds.add(Vec3::from(n.xyz));
        }
        let size = bounds.diagonal();
        // The form is in the display units; the scene in metres.
        let (center, radius) = (self.shown_point(at), self.shown_length(size / 20.0));
        let Some(d) = self.dialog.as_mut().filter(|d| d.command == "mesh.node.relax") else { return false };
        let mut params = serde_json::from_str::<Value>(&d.text).ok().filter(Value::is_object).unwrap_or_else(|| json!({}));
        params["center"] = json!(center);
        if let Some(o) = params.as_object_mut() {
            for key in ["around_element", "nodes", "elements", "faces"] {
                o.remove(key);
            }
        }
        if params.get("radius").and_then(Value::as_f64).is_none_or(|r| !(r > 0.0)) && size > 0.0 && size.is_finite() {
            params["radius"] = json!(radius);
        }
        d.text = params.to_string();
        d.touched.insert("center".into());
        true
    }

    /// Scale the open Relax Nodes form's radius by `factor` (Alt + scroll over the view).
    pub fn relax_scale_radius(&mut self, factor: f64) -> bool {
        let Some(sphere) = self
            .dialog
            .as_ref()
            .filter(|d| d.command == "mesh.node.relax")
            .and_then(|d| serde_json::from_str::<Value>(&d.text).ok().and_then(|p| fieldcraft_engine::cmd::relax_sphere(&self.session, &p)))
        else {
            return false;
        };
        let radius = self.shown_length(sphere.1 * factor);
        let Some(d) = self.dialog.as_mut() else { return false };
        let mut params = serde_json::from_str::<Value>(&d.text).ok().filter(Value::is_object).unwrap_or_else(|| json!({}));
        params["radius"] = json!(radius);
        d.text = params.to_string();
        true
    }

    /// The datum the open datum form defines so far, drawn before OK.
    /// Worked out once per form text and model revision (a face centre costs a design query).
    fn datum_preview(&self) -> Option<fieldcraft_model::Datum> {
        let d = self.dialog.as_ref()?;
        let key = (self.revision, d.command, d.text.clone());
        let mut cache = self.preview.borrow_mut();
        if cache.0.as_ref() != Some(&key) {
            let params = serde_json::from_str::<Value>(&d.text).ok();
            let datum = params.and_then(|p| fieldcraft_engine::cmd::datum_preview(&self.session, d.command, &p));
            *cache = (Some(key), datum);
        }
        cache.1
    }

    /// The definition selected in the browser, when it has symbols.
    fn highlighted(&self) -> Option<(Definition, &str)> {
        let (def, name) = self.selected.as_ref()?;
        let def = def.symbols()?;
        Some((def, name.as_str()))
    }

    /// The setup's symbols for the viewport as drawn last frame (its height in points).
    pub fn symbols(&mut self) -> (u64, &Scene) {
        let height = self.viewport.rect.map_or(DEFAULT_VIEW_HEIGHT, |r| f64::from(r.height()));
        self.symbols_for(height, 1.0)
    }

    /// The setup's symbols at the size setting on screen (in points), for a view `height` pixels
    /// high at `pixels_per_point` ([`glyphs::unit_for`]) with the current camera: a layer of their own, rebuilt when the model, what
    /// is shown, the plot, the highlighted definition or the size changes, and while zooming
    /// only when the symbol unit moves to another quarter of a doubling (the mesh's scene is
    /// never rebuilt for it). Their nodes follow the drawn shape (a results plot's deformed one).
    pub fn symbols_for(&mut self, height: f64, pixels_per_point: f64) -> (u64, &Scene) {
        self.display();
        self.plot();
        let unit = glyphs::unit_for(&self.cam, height, self.ui.symbol_size, pixels_per_point).and_then(glyphs::quantise_unit);
        let mut h = std::hash::DefaultHasher::new();
        let plot = self.plot.1.as_ref().map(|p| Arc::as_ptr(p) as usize);
        let show = self.ui.show();
        (self.revision, show.loads, show.interactions, self.ui.result_symbols, self.ui.dark, plot, unit.map(|u| u.0), self.highlighted())
            .hash(&mut h);
        let key = h.finish();
        if self.symbols.0 != Some(key) {
            // Over a results plot: only when asked for, and only the shown step's.
            let step = self.plot.1.as_ref().map(|p| p.view.step.clone());
            if let Some(step) = step.as_ref().filter(|_| self.ui.result_symbols) {
                let for_plot = Some((self.revision, step.clone()));
                if self.step_glyphs.0 != for_plot {
                    self.step_glyphs = (for_plot, self.session.step_glyphs(step));
                }
            }
            let glyphs = if step.is_some() { &self.step_glyphs.1 } else { &self.glyphs };
            let t = theme::Tokens::get(self.ui.dark);
            let drawn = self.plot.1.as_ref().map_or(&self.display.1, |p| &p.display);
            let mut s = Scene::default();
            let unit = unit.map_or(f64::NAN, |u| u.1);
            let d = if step.is_some() && !self.ui.result_symbols {
                Drawn::default()
            } else {
                scene::symbols(&mut s, glyphs, |n| drawn.node_position(n), show, step.is_some(), unit, self.highlighted(), &t)
            };
            self.symbols = (Some(key), s, d);
        }
        (key, &self.symbols.1)
    }

    /// What the scene built last draws, as framing sees it: the model, a results plot's
    /// (deformed) shape and the reference points. Datums are never framed (#295): they are drawn
    /// at the model's extent, so framing them would only shrink the part.
    fn drawn_bounds(&self) -> Aabb3 {
        let mut bounds = self.display.1.bounds;
        if let Some(p) = &self.plot.1 {
            bounds.union(p.display.bounds);
            // A fluid step's slice and streamlines fill its domain.
            if let Some(f) = &p.flow {
                bounds.union(f.bounds());
            }
        }
        // Reference points are places in the model, framed with it (their symbols are not).
        for p in self.session.model.reference_points().values() {
            let p = Vec3::from(p.xyz);
            if p.is_finite() {
                bounds.add(p);
            }
        }
        bounds
    }

    /// A render of the viewport's scene with the current camera and the selection, as PNG
    /// bytes, `w` × `h` pixels at `pixels_per_point`: symbols, lines and nodes are as large as
    /// in a window of that scale (#285). Drawn with the viewport's pipelines on the window's
    /// device when there is one, else (or when the GPU fails) by the CPU rasterizer, with why;
    /// the second value says which.
    pub fn render_png(&mut self, w: usize, h: usize, pixels_per_point: f64) -> Result<(Vec<u8>, fieldcraft_gpu::Renderer), String> {
        let cam = self.cam;
        let (_, overlay) = self.overlay();
        self.symbols_for(h as f64, pixels_per_point);
        self.scene();
        let scaled;
        let (scene, symbols, overlay) = if pixels_per_point == 1.0 {
            (&self.scene.1, &self.symbols.1, &*overlay)
        } else {
            let (mut scene, mut symbols, mut overlay) = (self.scene.1.clone(), self.symbols.1.clone(), (*overlay).clone());
            for s in [&mut scene, &mut symbols, &mut overlay] {
                s.scale_marks(pixels_per_point);
            }
            scaled = (scene, symbols, overlay);
            (&scaled.0, &scaled.1, &scaled.2)
        };
        #[cfg(not(target_arch = "wasm32"))]
        let why = match &self.capturer {
            Ok(Some(c)) => {
                let (pw, ph) = (u32::try_from(w).unwrap_or(u32::MAX), u32::try_from(h).unwrap_or(u32::MAX));
                match c.rgba(scene, &[symbols, overlay], &cam, pw, ph) {
                    Ok(px) => {
                        let png = fieldcraft_render::png(w, h, &px).ok_or("the picture could not be encoded")?;
                        return Ok((png, fieldcraft_gpu::Renderer::Gpu { adapter: c.adapter.clone(), tier: c.tier, samples: c.samples() }));
                    }
                    Err(e @ fieldcraft_gpu::GpuError::Size(_)) => return Err(e.to_string()),
                    Err(e) => Some(e.to_string()),
                }
            }
            Ok(None) => None,
            Err(e) => Some(e.clone()),
        };
        #[cfg(target_arch = "wasm32")]
        let why = None;
        let png = fieldcraft_render::render_layers_png(scene, &[symbols, overlay], &cam, w, h).map_err(|e| e.to_string())?;
        Ok((png, fieldcraft_gpu::Renderer::Cpu { why }))
    }

    /// The geometry tags behind the drawn mesh, read again when the model changed.
    fn pick_index(&mut self) -> Arc<selection::Index> {
        self.display();
        if self.pick.index.0 != Some(self.revision) {
            self.pick.index = (Some(self.revision), Arc::new(selection::Index::new(&self.session.model, &self.display.1)));
        }
        self.pick.index.1.clone()
    }

    /// Run `f` with the picker of what the viewport draws: a results plot's elements, or the
    /// model's elements and its geometry where no element faces cover it (as [`scene::build`]).
    pub fn with_picker<R>(&mut self, f: impl FnOnce(&selection::Picker) -> R) -> R {
        let index = self.pick_index();
        let plot = self.plot();
        let (show_mesh, show_geometry) = (self.ui.show_mesh, self.ui.show_geometry);
        let display = plot.as_ref().map_or(&self.display.1, |p| &p.display);
        let mesh_faces = show_mesh && !display.tris.is_empty();
        let picker = selection::Picker {
            display,
            names: &self.pick.names,
            model: &self.session.model,
            index: &index,
            elements: show_mesh,
            geometry: show_geometry && !mesh_faces,
            curves: show_geometry,
        };
        f(&picker)
    }

    /// The viewport's projection as drawn last frame (`None` before it was drawn).
    pub fn projector(&mut self) -> Option<fieldcraft_render::pick::Projector> {
        let rect = self.viewport.rect?;
        let radius = self.scene().1.radius;
        Some(fieldcraft_render::pick::Projector::new(&self.cam, f64::from(rect.width()), f64::from(rect.height()), radius))
    }

    /// The nearest drawn point under window point `at`: the pick ray's first hit on the results
    /// plot when one shows, else on the shown elements' faces or the design's surfaces.
    pub fn surface_point(&mut self, at: egui::Pos2) -> Option<Vec3> {
        let rect = self.viewport.rect.filter(|r| r.contains(at))?;
        let (o, d) = self.cam.ray(f64::from(at.x - rect.left()), f64::from(at.y - rect.top()), f64::from(rect.width()), f64::from(rect.height()));
        let t = if let Some(plot) = self.plot() {
            plot.display.hit_element_face(o, d).map(|h| h.t)
        } else {
            let (mesh, geometry) = (self.ui.show_mesh, self.ui.show_geometry);
            let display = self.display();
            let faces = mesh.then(|| display.hit_element_face(o, d)).flatten().map(|h| h.t);
            let surface = geometry.then(|| display.hit_surface(o, d)).flatten().map(|h| h.t);
            faces.into_iter().chain(surface).min_by(f64::total_cmp)
        }?;
        let p = o + d * t;
        p.is_finite().then_some(p)
    }

    /// What an orbit that starts on empty space turns about: the centre of what is drawn (the
    /// target when nothing is).
    pub fn drawn_centre(&mut self) -> Vec3 {
        self.scene();
        let b = self.drawn_bounds();
        if b.is_empty() || !b.center().is_finite() { self.cam.target } else { b.center() }
    }

    /// What a click at window point `at` would pick with the current filter.
    pub fn pick_at(&mut self, at: egui::Pos2) -> Option<selection::Item> {
        let rect = self.viewport.rect.filter(|r| r.contains(at))?;
        let proj = self.projector()?;
        let filter = self.pick.filter;
        let local = [f64::from(at.x - rect.left()), f64::from(at.y - rect.top())];
        self.with_picker(|p| p.pick(&proj, local, filter))
    }

    /// Change the selection by `items` (see [`selection::Mode`]); the status bar says what is
    /// selected.
    pub fn select(&mut self, items: impl IntoIterator<Item = selection::Item>, mode: selection::Mode) {
        if self.pick.selection.apply(items, mode) {
            self.pick.generation = self.pick.generation.wrapping_add(1);
        }
        // The filter label follows what the selection holds (a seed picked in the browser selects
        // its edge whatever the filter was, #341).
        if let Some(f) = self.pick.selection.filter() {
            self.pick.filter = f;
        }
        let mut text = format!("Selected: {}", self.pick.selection.summary());
        if self.pick.selection.capped() {
            text.push_str(&format!(" (at most {} items)", selection::MAX_SELECTION));
        }
        self.set_status(text, false);
    }

    pub fn clear_selection(&mut self) {
        if !self.pick.selection.is_empty() {
            self.pick.selection.clear();
            self.pick.generation = self.pick.generation.wrapping_add(1);
            self.set_status("Selection cleared", false);
        }
    }

    /// Pick at window point `at` and change the selection by it; a plain click on nothing clears
    /// the selection. Returns what was picked.
    pub fn select_at(&mut self, at: egui::Pos2, mode: selection::Mode) -> Option<selection::Item> {
        // A click where the hover was just worked out picks what is highlighted.
        let current =
            self.pick.hover_key.as_ref().is_some_and(|k| k.0 == Some(at) && k.1 == self.cam && k.2 == self.scene_key() && k.3 == self.pick.filter);
        let item = if current { self.pick.hover.clone() } else { self.pick_at(at) };
        match &item {
            Some(i) => self.select([i.clone()], mode),
            None if mode == selection::Mode::Replace => self.clear_selection(),
            None => {}
        }
        item
    }

    /// Select what is visible inside the box between window points `a` and `b` (nodes, element
    /// faces or elements, with the filter); returns how many items it held.
    pub fn select_box(&mut self, a: egui::Pos2, b: egui::Pos2, mode: selection::Mode) -> usize {
        let (Some(rect), Some(proj)) = (self.viewport.rect, self.projector()) else { return 0 };
        let filter = self.pick.filter;
        let r = [a.x - rect.left(), a.y - rect.top(), b.x - rect.left(), b.y - rect.top()].map(f64::from);
        let items = self.with_picker(|p| p.pick_box(&proj, r, filter));
        let n = items.len();
        self.select(items, mode);
        n
    }

    /// Set what is under the pointer (`None` when nothing, or the pointer is elsewhere).
    pub fn set_hover(&mut self, item: Option<selection::Item>) {
        if self.pick.hover != item {
            self.pick.hover = item;
            self.pick.hover_generation = self.pick.hover_generation.wrapping_add(1);
        }
    }

    /// Highlight `items` as hovered (the faces of a Feature Graph row under the pointer; empty:
    /// none).
    pub fn set_feature_hover(&mut self, items: Vec<selection::Item>) {
        if self.pick.feature_hover != items {
            self.pick.feature_hover = items;
            self.pick.hover_generation = self.pick.hover_generation.wrapping_add(1);
        }
    }

    /// The selection and the hover drawn over the scene in the theme's highlight colours,
    /// rebuilt when either or the scene changes; and its key.
    pub fn overlay(&mut self) -> (u64, Arc<Scene>) {
        let hash = |v: &dyn Fn(&mut std::hash::DefaultHasher)| {
            let mut h = std::hash::DefaultHasher::new();
            v(&mut h);
            h.finish()
        };
        let scene = self.scene_key();
        let set = self.selected_set();
        let cut = self.partition_preview();
        let merge = self.merge_preview();
        let relax_key = self.dialog.as_ref().filter(|d| d.command == "mesh.node.relax").map(|d| d.text.clone());
        let move_key = self.dialog.as_ref().filter(|d| d.command == node_move::COMMAND).map(|d| (d.text.clone(), self.node_move_pending()));
        let selected = hash(&|h| (scene, self.pick.generation, &set, &cut, &*merge, &relax_key, &move_key).hash(h));
        let key = hash(&|h| (selected, self.pick.hover_generation).hash(h));
        if self.pick.overlay.0 == Some(key) {
            return (key, self.pick.overlay.1.clone());
        }
        let t = theme::Tokens::get(self.ui.dark);
        // The selection's part is kept while only the hover changes.
        if self.pick.selected_overlay.0 != Some(selected) {
            let sel = std::mem::take(&mut self.pick.selection);
            let mut s = self.with_picker(|p| p.overlay(&sel, None, theme::rgb(t.pick_selected), theme::rgb(t.pick_hover)));
            self.pick.selection = sel;
            // The set selected in the browser: its members in the highlight of selected symbols.
            if let Some(items) = set.as_ref().and_then(|name| self.set_items(name)) {
                let highlight = t.glyph_colors().highlight;
                let m = self.with_picker(|p| p.items_overlay(&items, highlight));
                s.tris.extend(m.tris);
                s.lines.extend(m.lines);
                s.points.extend(m.points);
            }
            // The bodies the open Partition form would cut (#347).
            if !cut.is_empty() {
                let items: Vec<selection::Item> = cut
                    .iter()
                    .filter_map(|b| fieldcraft_model::GeometryRef::new(b.as_str(), fieldcraft_model::Entity::Body).ok())
                    .map(selection::Item::Geometry)
                    .collect();
                let m = self.with_picker(|p| p.items_overlay(&items, theme::rgb(t.pick_hover)));
                s.tris.extend(m.tris);
                s.lines.extend(m.lines);
            }
            // The nodes the open Merge Nodes form would merge (#234).
            if !merge.is_empty() {
                let items: Vec<selection::Item> = merge.iter().map(|&n| selection::Item::Node(n)).collect();
                let m = self.with_picker(|p| p.items_overlay(&items, theme::rgb(t.pick_hover)));
                s.points.extend(m.points);
            }
            // The Relax Nodes sphere and its preview (#398).
            if relax_key.is_some() {
                self.relax_overlay(&mut s, &t);
            }
            // The Move Nodes preview: the ghost, the moved elements and the sphere.
            if move_key.is_some() {
                self.node_move_overlay(&mut s, &t);
            }
            self.pick.selected_overlay = (Some(selected), Arc::new(s));
        }
        let mut s = (*self.pick.selected_overlay.1).clone();
        if let Some(hover) = self.pick.hover.clone() {
            // Over result contours the hover is an outline, so it never hides the result (#207).
            let fill = self.plot.1.is_none();
            let h = self.with_picker(|p| p.hover(&hover, theme::rgb(t.pick_hover), fill));
            s.tris.extend(h.tris);
            s.lines.extend(h.lines);
            s.points.extend(h.points);
        }
        if !self.pick.feature_hover.is_empty() {
            let items = self.pick.feature_hover.clone();
            let h = self.with_picker(|p| p.items_overlay(&items, theme::rgb(t.pick_hover)));
            s.tris.extend(h.tris);
            s.lines.extend(h.lines);
        }
        self.pick.overlay = (Some(key), Arc::new(s));
        (key, self.pick.overlay.1.clone())
    }

    /// The open Relax Nodes form drawn over the scene: its sphere of influence, translucent, and
    /// the edges of the elements it changes at their relaxed positions.
    fn relax_overlay(&self, s: &mut Scene, t: &theme::Tokens) {
        let sphere = self
            .dialog
            .as_ref()
            .and_then(|d| serde_json::from_str::<Value>(&d.text).ok())
            .and_then(|p| fieldcraft_engine::cmd::relax_sphere(&self.session, &p));
        if let Some((c, r)) = sphere {
            scene::sphere(s, Vec3::from(c), r, theme::rgb(t.pick_hover));
        }
        let Some(p) = self.relax_preview() else { return };
        let model = &self.session.model;
        let at = |n: u32| p.relaxed.positions.get(&n).copied().or_else(|| model.node(n).map(|x| x.xyz)).map(Vec3::from);
        let color = theme::rgb(t.pick_selected);
        for e in model.elements().values().filter(|e| e.nodes.iter().any(|n| p.relaxed.positions.contains_key(n))) {
            for edge in e.kind.edges() {
                let ends = (edge.first().and_then(|&i| e.nodes.get(i)), edge.get(1).and_then(|&i| e.nodes.get(i)));
                if let (Some(&a), Some(&b)) = ends
                    && let (Some(a), Some(b)) = (at(a), at(b))
                {
                    s.lines.push(fieldcraft_render::SceneLine { a, b, color, width: 1.5 });
                }
            }
        }
        for x in p.relaxed.positions.values() {
            s.points.push(fieldcraft_render::ScenePoint { p: Vec3::from(*x), color, size: 4.0 });
        }
    }

    /// The seed on `edge`, if it has one.
    pub fn seed(&self, edge: &fieldcraft_model::GeometryRef) -> Option<&fieldcraft_model::Seed> {
        self.session.model.seeds().iter().find(|s| s.edge == *edge)
    }

    /// Select (or with `None` deselect) a definition in the browser, as a click on its row does.
    pub fn select_definition(&mut self, selected: Option<(browser::Def, String)>) {
        if selected.is_some() {
            self.selected_seed = None;
            self.selected_mesh = false;
        }
        self.selected = selected;
    }

    /// Select what `mesh.check.duplicates` reported: the coincident nodes, the repeated elements
    /// and the faces of cracks.
    fn select_duplicates(&mut self, v: &Value) {
        let label = |x: &Value| x.as_u64().and_then(|l| u32::try_from(l).ok());
        let mut items: Vec<selection::Item> = Vec::new();
        for pair in v["coincident"].as_array().into_iter().flatten() {
            items.extend(pair.as_array().into_iter().flatten().take(2).filter_map(label).map(selection::Item::Node));
        }
        for group in v["duplicate_elements"].as_array().into_iter().flatten() {
            items.extend(group.as_array().into_iter().flatten().filter_map(label).map(selection::Item::Element));
        }
        items.extend(v["cracks"].as_array().into_iter().flatten().filter_map(|f| serde_json::from_value(f.clone()).ok()).map(selection::Item::Face));
        if items.is_empty() {
            self.clear_selection();
        } else {
            self.select(items, selection::Mode::Replace);
        }
    }

    /// Select the Mesh node in the browser: how the current mesh was made, on the right.
    pub fn select_mesh(&mut self) {
        self.selected = None;
        self.selected_seed = None;
        self.selected_mesh = true;
    }

    /// The meshing commands that made the current mesh with their parameters, oldest first: from
    /// the journal, the latest run for each body that still has elements (a remesh replaces an
    /// earlier run, a mirror adds its image). Empty for a mesh no mesher made (imported, built
    /// node by node).
    pub fn mesh_history(&self) -> Vec<(String, Value)> {
        // The journal's numbers are in its own units; the forms show the display units.
        let mut covered: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        let mut out = Vec::new();
        for (command, params) in self.session.journal.iter().rev() {
            let Some((_, targets)) = fieldcraft_engine::MESHERS.iter().find(|(m, _)| m == command) else { continue };
            let Ok(bodies) = targets(&self.session, params) else { continue };
            let fresh: Vec<String> =
                bodies.into_iter().filter(|b| !covered.contains(b) && !self.session.model.body_mesh(b).elements.is_empty()).collect();
            if fresh.is_empty() {
                continue;
            }
            covered.extend(fresh);
            out.push((command.clone(), self.session.journal_params_shown(command, params)));
        }
        out.reverse();
        out
    }

    /// Select the seed on `edge` in the browser: its properties on the right, its edge picked in
    /// the viewport.
    pub fn select_seed(&mut self, edge: fieldcraft_model::GeometryRef) {
        self.selected = None;
        self.selected_mesh = false;
        self.select([selection::Item::Geometry(edge.clone())], selection::Mode::Replace);
        self.selected_seed = Some(edge);
    }

    /// The set selected in the browser, if a set is.
    fn selected_set(&self) -> Option<String> {
        match &self.selected {
            Some((browser::Def::Set, name)) => Some(name.clone()),
            _ => None,
        }
    }

    /// Set `name`'s members as viewport items (design entities by name, mesh items by label).
    fn set_items(&self, name: &str) -> Option<Vec<selection::Item>> {
        let set = self.session.model.sets().get(name)?;
        let items = set.geometry.iter().cloned().map(selection::Item::Geometry);
        let items = items.chain(set.nodes.iter().map(|n| selection::Item::Node(*n)));
        let items = items.chain(set.elements.iter().map(|e| selection::Item::Element(*e)));
        Some(items.chain(set.faces.iter().map(|f| selection::Item::Face(*f))).collect())
    }

    /// Drop selected nodes, elements and element faces the model no longer has (after an edit,
    /// an undo or a remesh); geometry stays selected by name.
    fn check_selection(&mut self) {
        if self.pick.checked == Some(self.revision) {
            return;
        }
        self.pick.checked = Some(self.revision);
        let model = &self.session.model;
        let before = self.pick.selection.len();
        self.pick.selection.retain(|item| match item {
            selection::Item::Geometry(_) => true,
            selection::Item::Node(n) => model.node(*n).is_some(),
            selection::Item::Element(e) => model.element(*e).is_some(),
            selection::Item::Face(f) => model.face_nodes(*f).is_some(),
        });
        if self.pick.selection.len() != before {
            self.pick.generation = self.pick.generation.wrapping_add(1);
        }
    }

    /// Open `model.set.create`'s form with the selection as the new set's members.
    pub fn create_set_from_selection(&mut self) -> Result<(), String> {
        let params = self.selection_set(None)?;
        self.open_form("model.set.create", params);
        Ok(())
    }

    /// `model.set.create`'s parameters for a set of the selection, named after `base` (`Set-1`
    /// by default; the first name not taken).
    fn selection_set(&self, base: Option<&str>) -> Result<Value, String> {
        if self.pick.selection.is_empty() {
            return Err("nothing is selected".into());
        }
        let name = selection::fresh_set_name(&self.session.model, base.unwrap_or("Set-1"));
        Ok(selection::with_members(&json!({"name": name}), &self.pick.selection))
    }

    /// Put the selection into the open command form: a set's members, or a new set made from it
    /// (`model.set.create`, named after the definition) as the target of a boundary condition,
    /// load, section or coupling.
    pub fn use_selection_in_dialog(&mut self) -> Result<(), String> {
        let Some(d) = self.dialog.as_ref() else { return Err("no command form is open".into()) };
        let command = d.command;
        if self.pick.selection.is_empty() {
            return Err("nothing is selected".into());
        }
        let params: Value = serde_json::from_str(&d.text).map_err(|e| format!("the parameters are not JSON: {e}"))?;
        if !params.is_object() {
            return Err("the parameters must be a JSON object".into());
        }
        let Some(target) = selection::target(command, &params) else { return Err(format!("{command} takes no selection")) };
        let filled = match target {
            selection::Target::Members => selection::with_members(&params, &self.pick.selection),
            selection::Target::Edges => {
                selection::with_edges(&params, &self.pick.selection).ok_or("no edge is selected: pick edges with Select: Edges")?
            }
            selection::Target::Position | selection::Target::DatumPoint => {
                // A query that fills the form: it changes nothing and has no reply panel of its own.
                #[allow(clippy::disallowed_methods, reason = "the selection's centre: a sum over the picked items")]
                let centre = self.session.execute("geometry.centre", &self.pick.selection.params()).map_err(|e| e.to_string())?;
                let mut out = params.as_object().cloned().unwrap_or_default();
                if target == selection::Target::Position {
                    out.insert("xyz".into(), centre.get("point").cloned().unwrap_or(Value::Null));
                } else {
                    out.insert("by".into(), selection::datum_point_by(&self.pick.selection, &centre));
                }
                Value::Object(out)
            }
            selection::Target::ExtrudeFaces => selection::with_extrude_faces(&params, &self.pick.selection)
                .ok_or("no face is selected: pick the source face, then the target surface (Select: Faces)")?,
            selection::Target::Faces => {
                selection::with_faces(&params, &self.pick.selection).ok_or("no face is selected: pick the faces (Select: Faces)")?
            }
            selection::Target::Bodies => {
                let bodies = selection::bodies(&self.pick.selection).ok_or("no geometry is selected: pick cells with Select: Bodies")?;
                let mut out = params.as_object().cloned().unwrap_or_default();
                out.insert("bodies".into(), json!(bodies));
                Value::Object(out)
            }
            selection::Target::Set | selection::Target::TargetSet => {
                let base = params.get("name").and_then(Value::as_str).map(|n| format!("{n}-set"));
                let set = self.selection_set(base.as_deref())?;
                // A set with exactly these members is used as it is; otherwise the new one is made
                // on OK ([`ParamsDialog::pending_set`]).
                let (name, pending) = match selection::same_set(&self.session.model, &set) {
                    Some(existing) => (existing, None),
                    None => (set.get("name").and_then(Value::as_str).unwrap_or_default().to_string(), Some(set)),
                };
                if let Some(open) = self.dialog.as_mut() {
                    open.pending_set = pending;
                }
                selection::with_set(&params, target, &name)
            }
        };
        if let Some(open) = self.dialog.as_mut() {
            open.text = serde_json::to_string_pretty(&filled).unwrap_or_else(|_| filled.to_string());
            open.error = None;
        }
        Ok(())
    }

    /// Frame what is drawn in the part of the viewport the overlays leave free (see
    /// [`framing`]), keeping the direction of view.
    pub fn fit_view(&mut self) {
        // A view asked for replaces one still waiting to be done.
        self.reframe = None;
        self.scene();
        let b = self.drawn_bounds();
        let (aspect, free) = self.free_part();
        self.cam.fit_in(&b, aspect, free);
    }

    /// Select elements `labels` in the viewport and frame them.
    pub fn show_elements(&mut self, labels: &[u32]) {
        self.select(labels.iter().map(|l| selection::Item::Element(*l)), selection::Mode::Replace);
        self.frame_elements(labels);
    }

    /// Frame elements `labels` (as drawn undeformed), keeping the direction of view; a tiny
    /// element is framed with some of the model around it.
    pub fn frame_elements(&mut self, labels: &[u32]) {
        let model = &self.session.model;
        let mut b = Aabb3::EMPTY;
        for e in labels.iter().filter_map(|l| model.element(*l)) {
            for n in e.nodes.iter().filter_map(|n| model.node(*n)) {
                b.add(Vec3::from(n.xyz));
            }
        }
        if b.is_empty() {
            return;
        }
        self.reframe = None;
        self.scene();
        // At least an eighth of the model across.
        let pad = (self.drawn_bounds().diagonal() / 8.0 - b.diagonal()).max(0.0) / 2.0;
        let d = Vec3::new(pad, pad, pad);
        let (lo, hi) = (b.min - d, b.max + d);
        b.add(lo);
        b.add(hi);
        let (aspect, free) = self.free_part();
        self.cam.fit_in(&b, aspect, free);
    }

    /// The viewport's aspect and the part of it the overlays leave free, as fractions (the whole
    /// of a square view before the viewport is first drawn).
    fn free_part(&self) -> (f64, [f64; 4]) {
        match self.viewport.rect.filter(|r| r.width() > 1.0 && r.height() > 1.0) {
            Some(r) => (f64::from(r.width()) / f64::from(r.height()), framing::free_part(r, &self.viewport.overlays)),
            None => (1.0, fieldcraft_render::FULL_VIEW),
        }
    }

    /// Whether every corner of what is drawn shows in the free part of the viewport.
    fn framed(&mut self) -> bool {
        let Some(r) = self.viewport.rect.filter(|r| r.width() > 1.0 && r.height() > 1.0) else { return true };
        self.scene();
        let b = self.drawn_bounds();
        if b.is_empty() {
            return true;
        }
        let [x0, y0, x1, y1] = framing::free_part(r, &self.viewport.overlays);
        let (w, h) = (f64::from(r.width()), f64::from(r.height()));
        let p = fieldcraft_render::pick::Projector::new(&self.cam, w, h, self.scene.1.radius);
        b.corners().into_iter().all(|c| p.screen(c).is_some_and(|[x, y, _]| x >= x0 * w && x <= x1 * w && y >= y0 * h && y <= y1 * h))
    }

    /// Another model, sample or design opened (or the window's first one): a fresh start's home
    /// view when the viewport next draws, and nothing of the last document left selected or
    /// reported in the context panel (a definition of the same name is another one).
    fn frame_new_document(&mut self) {
        if self.framed_document != Some(self.session.opened()) {
            if self.framed_document.is_some() {
                self.selected = None;
                self.selected_seed = None;
                self.selected_mesh = false;
                self.reply = None;
                // Picks are by name, and the names of another document are other things.
                if !self.pick.selection.is_empty() {
                    self.pick.selection.clear();
                    self.pick.generation = self.pick.generation.wrapping_add(1);
                }
            }
            self.framed_document = Some(self.session.opened());
            // The probe and the reaction sums were of the old document's results.
            self.probe = None;
            self.reactions = None;
            self.xy = None;
            self.reframe = framing::Reframe::merge(self.reframe, framing::Reframe::Home);
        }
    }

    /// The home view: isometric, the whole model framed.
    pub fn home_view(&mut self) {
        self.set_view(StandardView::Iso);
    }

    /// Look along a standard view, the whole model framed.
    pub fn set_view(&mut self, v: StandardView) {
        self.cam.set_view(v);
        self.fit_view();
    }

    /// Aim the camera as `aim` says ([`fieldcraft_engine::view::camera`], the control channel's
    /// `ui.camera`): turned, and with a target moved to look at it from the view's free part, as
    /// a standard view moves it. A framing still waiting is dropped.
    pub fn aim_camera(&mut self, aim: &fieldcraft_engine::view::Aim) -> Result<(), String> {
        let (aspect, _) = self.free_part();
        self.cam = fieldcraft_engine::view::aim(&self.session, aim, self.cam, aspect).map_err(|e| e.to_string())?;
        self.reframe = None;
        Ok(())
    }

    /// Aim at the camera a saved view asks for (`view.restore`, a model opened showing one), once
    /// per request; then keep the session's camera the viewport's, so `view.save` keeps it.
    fn follow_camera(&mut self) {
        if self.session.camera_requests() != self.camera_seen {
            self.camera_seen = self.session.camera_requests();
            if let Some(aim) = self.session.camera.clone()
                && let Err(e) = self.aim_camera(&aim)
            {
                self.set_status(format!("Saved view: {e}"), true);
            }
        }
        self.session.camera = Some(fieldcraft_engine::view::aim_of(&self.cam, &self.session));
    }

    /// Do the framing asked for ([`framing::Reframe`]) once there is something to frame; the
    /// viewport calls it when it knows its size and its overlays, before drawing.
    pub(crate) fn apply_reframe(&mut self) {
        let Some(r) = self.reframe else { return };
        if self.display().bounds.is_empty() {
            return;
        }
        self.reframe = None;
        match r {
            framing::Reframe::Home => self.home_view(),
            framing::Reframe::Fit => self.fit_view(),
            framing::Reframe::IfOffscreen => {
                if !self.framed() {
                    self.fit_view();
                }
            }
        }
    }

    /// Per-frame logic before layout: theme, camera projection, control requests, screenshots.
    pub fn logic(&mut self, ctx: &egui::Context) {
        self.now = ctx.input(|i| i.time);
        if self.styled != Some(self.ui.dark) {
            theme::apply(ctx, self.ui.dark);
            self.styled = Some(self.ui.dark);
        }
        self.drain_control(ctx);
        if self.poll_loading() {
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
        // A new document's home framing first, so the saved view it opens showing turns the
        // camera after it (and drops that framing) rather than being framed over.
        self.frame_new_document();
        self.follow_camera();
        self.dropped_files(ctx);
        self.poll_updates(ctx);
        self.shortcuts(ctx);
        self.cam.fov = if self.ui.perspective { PERSPECTIVE_FOV } else { 0.0 };
        if !self.synthetic.is_empty() {
            ctx.request_repaint();
        }
        // Results coming up (a solve, a results command) bring the Results tab.
        let shown = self.session.results_view.is_some();
        if shown && !self.results_shown {
            self.ui.select_tab(Tab::Results);
            self.reframe = framing::Reframe::merge(self.reframe, framing::Reframe::IfOffscreen);
        }
        self.results_shown = shown;
        self.check_selection();
        self.follow_form_filter();
        // A refusal goes once the command can run.
        if self.refused.as_ref().is_some_and(|u| self.unavailable(u.command).is_none()) {
            self.refused = None;
        }
        // A selected definition that is gone (undo, delete) is no longer selected.
        if self.selected.as_ref().is_some_and(|(kind, name)| !kind.exists(&self.session.model, name)) {
            self.selected = None;
        }
        if self.selected_seed.as_ref().is_some_and(|edge| self.seed(edge).is_none()) {
            self.selected_seed = None;
        }
        if self.session.results_view.as_ref().is_some_and(|v| v.animate) {
            ctx.request_repaint_after(std::time::Duration::from_secs_f64(results::PERIOD / results::FRAMES as f64));
        }
        self.screenshots(ctx);
    }

    /// A form that picks one kind of entity (Seed Edges: edges) switches the selection filter to
    /// it while open, and the filter goes back when the form closes. A filter the user changes
    /// while the form is open stays until the form wants another kind.
    fn follow_form_filter(&mut self) {
        let want = self
            .dialog
            .as_ref()
            .and_then(|d| serde_json::from_str::<Value>(&d.text).ok().and_then(|p| selection::form_filter(d.command, &p)).map(|f| (d.opening, f)));
        match (want, self.pick.form_filter) {
            (Some((opening, f)), Some((o, set, _))) if opening == o && f == set => {}
            (Some((opening, f)), previous) => {
                let restore = previous.map_or(self.pick.filter, |(_, _, r)| r);
                self.pick.form_filter = Some((opening, f, restore));
                self.pick.filter = f;
            }
            (None, Some((_, _, restore))) => {
                self.pick.form_filter = None;
                self.pick.filter = restore;
            }
            (None, None) => {}
        }
    }

    /// Undo (Ctrl/Cmd+Z), redo (Ctrl/Cmd+Shift+Z or Ctrl/Cmd+Y), open (Ctrl/Cmd+O), save
    /// (Ctrl/Cmd+S), save as (Ctrl/Cmd+Shift+S) and delete the browser's selection (Delete),
    /// unless a text field has the keyboard (it undoes its own typing).
    fn shortcuts(&mut self, ctx: &egui::Context) {
        if ctx.memory(|m| m.focused().is_some()) {
            return;
        }
        // Escape cancels a Move Nodes drag first.
        if self.node_drag.is_some() && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape)) {
            self.node_move_cancel();
        }
        // Escape clears the viewport's selection first (a second one cancels an open form).
        if !self.pick.selection.is_empty()
            && self.ui.menu.is_none()
            && !self.ui.palette
            && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape))
        {
            self.clear_selection();
        }
        let shortcut = |modifiers, key| egui::KeyboardShortcut::new(modifiers, key);
        let shift = egui::Modifiers::COMMAND | egui::Modifiers::SHIFT;
        // Shift first: the plain shortcut matches it too.
        let (redo, undo, save_as, save, open, new) = ctx.input_mut(|i| {
            let redo = i.consume_shortcut(&shortcut(shift, egui::Key::Z)) || i.consume_shortcut(&shortcut(egui::Modifiers::COMMAND, egui::Key::Y));
            let undo = i.consume_shortcut(&shortcut(egui::Modifiers::COMMAND, egui::Key::Z));
            let save_as = i.consume_shortcut(&shortcut(shift, egui::Key::S));
            let save = i.consume_shortcut(&shortcut(egui::Modifiers::COMMAND, egui::Key::S));
            let open = i.consume_shortcut(&shortcut(egui::Modifiers::COMMAND, egui::Key::O));
            (redo, undo, save_as, save, open, i.consume_shortcut(&shortcut(egui::Modifiers::COMMAND, egui::Key::N)))
        });
        if redo {
            let _ = self.run("edit.redo", json!({}));
        } else if undo {
            let _ = self.run("edit.undo", json!({}));
        }
        if save_as {
            self.save_as();
        } else if save {
            self.save();
        }
        if open {
            self.open_dialog(None);
        } else if new {
            self.new_model();
        }
        // Delete (or Backspace) deletes the definition selected in the browser, after asking.
        if self.prompt.is_none()
            && let Some((kind, name)) = self.selected.clone()
            && ctx
                .input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Delete) || i.consume_key(egui::Modifiers::NONE, egui::Key::Backspace))
        {
            self.ask_delete(kind, &name);
        }
        // The command search: S (as in SolidCraft) or Cmd/Ctrl+K.
        let search = ctx.input_mut(|i| {
            i.consume_shortcut(&shortcut(egui::Modifiers::COMMAND, egui::Key::K)) || i.consume_key(egui::Modifiers::NONE, egui::Key::S)
        });
        if search {
            self.ui.palette = true;
            self.ui.menu = None;
        }
    }

    /// Inject synthetic events: one pointer event per frame (egui reads a press and release in
    /// the same frame as a click, but a drag needs frames between its moves), keys up to their
    /// release.
    pub fn raw_input_hook(&mut self, raw: &mut egui::RawInput) {
        let Some(first) = self.synthetic.first() else { return };
        let n = match first {
            egui::Event::PointerMoved(_) | egui::Event::PointerButton { .. } | egui::Event::MouseWheel { .. } => 1,
            _ => self.synthetic.iter().position(|e| matches!(e, egui::Event::Key { pressed: false, .. })).map_or(self.synthetic.len(), |i| i + 1),
        };
        if let egui::Event::PointerMoved(p) | egui::Event::PointerButton { pos: p, .. } = first {
            raw.events.push(egui::Event::PointerMoved(*p));
        }
        // egui keeps the modifier keys as input state, changed by their own event.
        match first {
            egui::Event::PointerButton { modifiers, pressed, .. } => {
                // A drag keeps its press's modifiers until the release.
                self.synthetic_mods = if *pressed { *modifiers } else { egui::Modifiers::default() };
                raw.events.push(egui::Event::ModifiersChanged(*modifiers));
            }
            egui::Event::Key { modifiers, .. } => raw.events.push(egui::Event::ModifiersChanged(*modifiers)),
            egui::Event::PointerMoved(_) => raw.events.push(egui::Event::ModifiersChanged(self.synthetic_mods)),
            _ => {}
        }
        let n = n.min(self.synthetic.len());
        raw.events.extend(self.synthetic.drain(..n));
    }

    /// Lay out the whole window.
    pub fn ui(&mut self, ui: &mut egui::Ui) {
        self.show_gpu_errors();
        self.layout.clear();
        appbar::app_bar(self, ui);
        toolbar::toolbar(self, ui);
        if self.ui.show_timeline {
            timeline::timeline(self, ui);
        }
        palette::command_line(self, ui);
        if self.ui.show_tree {
            browser::browser(self, ui);
        }
        context::panel(self, ui);
        egui::CentralPanel::default().frame(egui::Frame::NONE).show(ui, |ui| viewport::show(self, ui));
        palette::popup(self, ui.ctx());
        files::show(self, ui.ctx());
        report::show(self, ui.ctx());
        self.ask_consent();
        privacy::show(self, ui.ctx());
        updates::show(self, ui.ctx());
        if self.custom_titlebar {
            titlebar::resize_zones(ui);
        }
        self.save_prefs();
    }

    /// The setup progress ([`fieldcraft_engine::progress`]) the browser marks, worked out again
    /// only when the model or its results change.
    pub fn setup_stages(&mut self) -> std::sync::Arc<Vec<fieldcraft_engine::progress::StageState>> {
        use fieldcraft_engine::requirement::Requirement;
        let key = (self.revision, self.session.revision(), Requirement::Result.met(&self.session), Requirement::ResultsShown.met(&self.session));
        if self.setup.0 != Some(key) {
            self.setup = (Some(key), std::sync::Arc::new(self.session.progress()));
        }
        self.setup.1.clone()
    }

    /// Save the preferences the window changed (a View menu toggle, `ui.set`).
    /// Show and take numbers in `units` from now on, in models without their own, and remember
    /// it (Help › Settings). Nothing in the model changes: it is stored in SI.
    pub fn set_default_units(&mut self, units: fieldcraft_model::DisplaySystem) {
        self.session.default_units = units.clone();
        self.prefs.units = units;
        self.prefs.save(self.services.prefs.as_ref());
        self.follow_units();
    }

    pub(crate) fn save_prefs(&mut self) {
        if self.ui.nav_preset != self.prefs.nav_preset || self.ui.reverse_zoom != self.prefs.reverse_zoom {
            self.prefs.nav_preset = self.ui.nav_preset.clone();
            self.prefs.reverse_zoom = self.ui.reverse_zoom;
            self.prefs.save(self.services.prefs.as_ref());
        }
    }

    /// Remember where a region of the window was drawn (`ui.inspect`'s `layout`).
    pub fn place(&mut self, region: &'static str, rect: egui::Rect) {
        self.layout.insert(region, rect);
    }

    /// Open, or close when it is open, a menu (see [`UiState::menu`]).
    pub fn toggle_menu(&mut self, key: &str) {
        self.ui.menu = if self.ui.menu.as_deref() == Some(key) { None } else { Some(key.to_string()) };
    }

    fn drain_control(&mut self, ctx: &egui::Context) {
        /// Requests run on this thread, so a flood would starve drawing; the rest wait for the
        /// next frame.
        const PER_FRAME: usize = 64;
        let Some(rx) = self.control_rx.take() else { return };
        for _ in 0..PER_FRAME {
            let Ok(req) = rx.try_recv() else { break };
            // The MCP server says so; anything else on the channel is a script.
            let origin = match req.params.get("client").and_then(Value::as_str) {
                Some("mcp") => fieldcraft_engine::usage::Origin::Mcp,
                _ => fieldcraft_engine::usage::Origin::Script,
            };
            match self.with_origin(origin, |app| control::handle(app, ctx, &req)) {
                control::Outcome::Done(v) => {
                    let _ = req.reply.send(v);
                }
                control::Outcome::Pending => {
                    let _ = req.reply.send(json!({control::PENDING: true}));
                }
                control::Outcome::Screenshot { path } => {
                    self.shot_token = self.shot_token.wrapping_add(1);
                    let deadline = self.now + SHOT_TIMEOUT;
                    self.shots.push(Shot { token: self.shot_token, path, reply: req.reply.clone(), frames: 0, requested: false, deadline });
                }
            }
            ctx.request_repaint();
        }
        self.control_rx = Some(rx);
    }

    /// Ask for pending screenshots once their request's changes are drawn, and save the ones
    /// that arrived.
    fn screenshots(&mut self, ctx: &egui::Context) {
        if self.shots.is_empty() {
            return;
        }
        let arrived: Vec<(u64, std::sync::Arc<egui::ColorImage>)> = ctx.input(|i| {
            i.raw
                .events
                .iter()
                .filter_map(|e| match e {
                    egui::Event::Screenshot { user_data, image, .. } => {
                        user_data.data.as_ref().and_then(|d| d.downcast_ref::<u64>()).map(|t| (*t, image.clone()))
                    }
                    _ => None,
                })
                .collect()
        });
        for (token, image) in arrived {
            if let Some(i) = self.shots.iter().position(|s| s.token == token) {
                let shot = self.shots.remove(i);
                let _ = shot.reply.send(control::save_screenshot(&image, shot.path.as_deref()));
            }
        }
        let now = self.now;
        // No frame was presented in time (a hidden window, a locked screen): draw the window
        // offscreen instead.
        let (late, waiting): (Vec<Shot>, Vec<Shot>) = std::mem::take(&mut self.shots).into_iter().partition(|s| now > s.deadline);
        self.shots = waiting;
        for s in late {
            let _ = s.reply.send(control::capture_to(self, control::window_size(ctx), s.path.as_deref()));
        }
        self.shots.retain_mut(|s| {
            s.frames += 1;
            if !s.requested && s.frames >= 3 {
                ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::UserData::new(s.token)));
                s.requested = true;
            }
            true
        });
        if !self.shots.is_empty() {
            ctx.request_repaint_after(std::time::Duration::from_millis(16));
        }
    }
}

/// A design's bodies tessellated for drawing and picking, each with its names.
type BodyMeshes = Vec<(selection::BodyNames, fieldcraft_geometry::SurfaceMesh)>;

/// Whether `a` and `b` are the same open design (or both none).
fn same_design(a: &Option<Arc<fieldcraft_geometry::Design>>, b: &Option<Arc<fieldcraft_geometry::Design>>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => Arc::ptr_eq(a, b),
        (None, None) => true,
        _ => false,
    }
}

/// Tessellate `design`'s bodies (seconds for a large part): the bodies drawn, and why the
/// others cannot be.
#[allow(clippy::disallowed_methods, reason = "on a job's worker (`Prebuilt`), or in the window when a command run there changed the design")]
fn tessellate(design: &fieldcraft_geometry::Design) -> (BodyMeshes, Vec<String>) {
    let mut failed = Vec::new();
    let meshes = design
        .bodies
        .iter()
        .filter_map(|body| {
            let tolerance = body.measure().map(|m| display_tolerance(&Aabb3 { min: m.bbox_min.into(), max: m.bbox_max.into() }));
            let mut mesh = tolerance.and_then(|t| body.tessellate(t)).map_err(|e| failed.push(format!("{}: {e}", body.name))).ok()?;
            // The geometry is in millimetres, the scene in metres.
            let m = fieldcraft_engine::units::geometry_point;
            let names = selection::BodyNames {
                name: body.name.clone(),
                faces: body.faces().iter().map(|f| f.name.clone()).collect(),
                edges: body.edge_names(&mesh),
                vertices: body.vertex_names(&mesh).into_iter().map(|v| (v.name, m(v.point).into())).collect(),
            };
            mesh.positions.iter_mut().chain(mesh.edges.iter_mut().flatten()).for_each(|p| *p = m(*p));
            Some((names, mesh))
        })
        .collect();
    (meshes, failed)
}

/// A session's display data, built by the job that made the session on its worker so that the
/// window only adopts it ([`ModelApp::display`]): drawing a large model and tessellating its
/// design take long enough to stall a frame (#552).
pub struct Prebuilt {
    /// The session revision it was built for.
    revision: u64,
    /// The model's elements as drawn ([`draw`]).
    model: DisplayMesh,
    /// The design's bodies tessellated and why the others cannot be drawn; `None` when the
    /// window has drawn this design already.
    bodies: Option<(BodyMeshes, Vec<String>)>,
}

impl Prebuilt {
    /// Build `session`'s display data (on a job's worker); its design is tessellated unless it
    /// is `drawn`, the one the window draws.
    pub fn of(session: &Session, drawn: &Option<Arc<fieldcraft_geometry::Design>>) -> Prebuilt {
        let bodies = (!same_design(&session.design, drawn)).then(|| session.design.as_deref().map(tessellate).unwrap_or_default());
        Prebuilt { revision: session.revision(), model: draw(session), bodies }
    }
}

/// The model's elements as drawn, without the hidden bodies'.
#[allow(clippy::disallowed_methods, reason = "on a job's worker (`Prebuilt`), or in the window after a change made there")]
fn draw(session: &Session) -> DisplayMesh {
    DisplayMesh::from_model_hiding(&session.model, &session.hidden_elements())
}

/// Seconds a screenshot waits for a presented frame.
const SHOT_TIMEOUT: f64 = 8.0;
/// Vertical field of view of the perspective camera (radians).
const PERSPECTIVE_FOV: f64 = 0.6;

/// A body is cut by a section plane when it has surface points farther from the plane than this
/// fraction of the model's diagonal, on both sides.
const SECTION_SIDE: f64 = 1e-6; // unit-factors: relative to the model's diagonal
/// The smallest diagonal [`SECTION_SIDE`] is taken of (metres), so an empty model still has a side.
const MIN_SECTION_DIAGONAL: f64 = 1e-9; // unit-factors: a floor in metres (SI)

impl ModelApp {
    /// Why a toolbar entry or a command cannot start now: what its command needs that the
    /// session lacks (`fieldcraft_engine::requirement`). `None` when it can start, and for
    /// coming-soon entries and unknown ids (which say so themselves).
    pub fn unavailable(&self, id: &str) -> Option<fieldcraft_engine::requirement::Unavailable> {
        let command = catalogue::find(id).map_or(Some(id), |c| c.command())?;
        let spec = fieldcraft_engine::find_command(command)?;
        self.session.availability(spec).err()
    }
}

/// A command's name as people see it (its label).
pub fn command_label(id: &str) -> &str {
    fieldcraft_engine::find_command(id).map_or(id, |c| c.label)
}

/// What a probe of `node` read (`solve.probe`'s reply), as (what, value) rows with units from
/// `units` (the reply's): displacement, reaction, stresses, temperature, velocity or a mode's shape, and the
/// contour shown when it is none of those. Round-off beside the result's largest value is 0.
fn probe_rows(reply: &Value, plot: Option<&results::Plot>, units: &DisplaySystem, node: u32, mode: usize) -> Vec<(String, String)> {
    use fieldcraft_engine::results::{Field, Pick, Quantity, displacements, magnitude, unit_of, with_unit as titled};
    let at = reply.pointer("/nodes/0").cloned().unwrap_or(Value::Null);
    let vector = |key: &str| at.get(key).and_then(|v| serde_json::from_value::<[f64; 3]>(v.clone()).ok());
    let unit = |q: Quantity| unit_of(q, units);
    // The largest displacement of the result shown, in the display units: the scale round-off
    // is read against.
    let largest = plot.map_or(0.0, |p| {
        let si = displacements(&p.result, Pick { mode, frame: p.pick.frame })
            .map_or(0.0, |u| u.iter().map(|v| magnitude(*v)).filter(|m| m.is_finite()).fold(0.0, f64::max));
        Field::Displacement.shown(&p.result, si, units)
    });
    let mut rows = Vec::new();
    let vector_row = |rows: &mut Vec<(String, String)>, name: &str, v: [f64; 3], q: Quantity, scale: f64| {
        let u = unit(q);
        rows.push((titled(name, u.as_deref()), results::components(&v, scale)));
        rows.push((format!("|{name}|"), results::with_unit(results::number_in(magnitude(v), scale), u.as_deref())));
    };
    if let Some(u) = vector("displacement") {
        vector_row(&mut rows, "U", u, Quantity::Length, largest);
    }
    if let Some(rf) = vector("reaction").filter(|r| r.iter().any(|c| *c != 0.0)) {
        rows.push((titled("RF", unit(Quantity::Force).as_deref()), results::components(&rf, 0.0)));
    }
    let stress = at.get("stress").and_then(|v| serde_json::from_value::<[f64; 6]>(v.clone()).ok());
    let mises = at.get("von_mises").and_then(Value::as_f64);
    let stress_unit = unit(Quantity::Pressure);
    if let Some(m) = mises {
        rows.push(("Mises".into(), results::with_unit(results::number(m), stress_unit.as_deref())));
    }
    if let Some(t) = stress {
        let scale = mises.unwrap_or(0.0);
        rows.push((titled("S xx, yy, zz", stress_unit.as_deref()), results::components(&t[..3], scale)));
        rows.push((titled("S xy, xz, yz", stress_unit.as_deref()), results::components(&t[3..], scale)));
        let p = fieldcraft_engine::results::principal(&t);
        rows.push((titled("Principal", stress_unit.as_deref()), results::components(&p, scale)));
    }
    if let Some(t) = at.get("temperature").and_then(Value::as_f64) {
        rows.push(("T".into(), results::with_unit(results::number(t), unit(Quantity::Temperature).as_deref())));
    }
    if let Some(q) = vector("heat_flux") {
        rows.push((titled("q", unit(Quantity::HeatFlux).as_deref()), results::components(&q, 0.0)));
    }
    if let Some(v) = vector("velocity") {
        vector_row(&mut rows, "V", v, Quantity::Velocity, 0.0);
    }
    if let Some(shape) = vector("shape") {
        rows.push((format!("Mode {mode} shape (normalised)"), results::components(&shape, largest)));
    }
    let contour = plot.and_then(|p| p.contour.as_ref().map(|c| (c.field, c.at(p.result.nodes(), node), &p.result)));
    if let Some((field, v, result)) =
        contour.filter(|(f, ..)| !matches!(f, Field::VonMises | Field::Temperature) && !f.label().starts_with("Displacement"))
    {
        rows.push((field.label().to_string(), results::with_unit(results::number(v), field.unit(result, units).as_deref())));
    }
    rows
}

/// The tests that open a GPU device take turns, for as long as the device lives: several devices
/// open at once in one process hung NVIDIA's driver on a GB10 (#614). A test that failed while
/// holding the turn leaves it to the next.
#[cfg(test)]
pub(crate) fn gpu_turn() -> std::sync::MutexGuard<'static, ()> {
    static TURN: std::sync::Mutex<()> = std::sync::Mutex::new(());
    TURN.lock().unwrap_or_else(|e| e.into_inner())
}

#[cfg(test)]
mod tests {
    use fieldcraft_render::Vec3;

    use super::*;

    #[test]
    fn commands_run_through_the_engine_and_report_in_the_status_bar() {
        let mut app = ModelApp::new(Session::default());
        let key = app.scene_key();
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        assert_eq!(app.status.as_ref().map(|s| s.text.as_str()), Some("Add Node done"));
        assert_ne!(app.scene_key(), key, "a model change rebuilds the scene");
        let err = app.run("mesh.node.add", json!({"label": 0, "xyz": [0, 0, 0]})).unwrap_err();
        assert!(app.status.as_ref().is_some_and(|s| s.error && s.text.starts_with("Add Node:")), "{err}");
        app.run("model.info", json!({})).unwrap();
        assert!(app.status.as_ref().is_some_and(|s| s.text.starts_with("Model Info: 1 node, 0 elements")), "{:?}", app.status);
    }

    /// Save View keeps the viewport's camera; showing the view turns the viewport back to it.
    #[test]
    fn a_saved_view_turns_the_viewport_back() {
        let mut app = ModelApp::new(Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.cam.yaw = 0.5;
        app.cam.pitch = 0.25;
        app.follow_camera();
        app.run("view.save", json!({"name": "Side"})).unwrap();
        assert!((app.session.view.views[0].camera.as_ref().unwrap()["yaw"].as_f64().unwrap() - 0.5f64.to_degrees()).abs() < 1e-9);
        app.cam.yaw = -1.0;
        app.run("view.restore", json!({"name": "Side"})).unwrap();
        app.follow_camera();
        assert!((app.cam.yaw - 0.5).abs() < 1e-9 && (app.cam.pitch - 0.25).abs() < 1e-9, "{:?}", app.cam);
        // Once per request: a later turn stays.
        app.cam.yaw = 1.0;
        app.follow_camera();
        assert_eq!(app.cam.yaw, 1.0);
    }

    /// A model saved showing a view opens with that view's camera, not the home view (#646).
    #[test]
    fn a_model_opens_with_the_camera_of_the_view_it_was_saved_showing() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-view-open-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("side.fieldcraft").to_string_lossy().into_owned();
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(Session::default());
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.run("solve.run", json!({})).unwrap();
        app.run("view.save", json!({"name": "Side", "camera": {"yaw": -60, "pitch": 25}})).unwrap();
        app.run("file.save_as", json!({"path": path})).unwrap();
        // Opened as the window opens files: in the background, from an empty session.
        for _ in 0..2 {
            app.run("file.new", json!({})).unwrap();
            draw(&mut app, &ctx, 1280.0, 800.0);
            app.start_job("file.open", json!({"path": path})).unwrap();
            while app.loading.is_some() {
                draw(&mut app, &ctx, 1280.0, 800.0);
                std::thread::sleep(std::time::Duration::from_millis(5));
            }
            draw(&mut app, &ctx, 1280.0, 800.0);
            assert!((app.cam.yaw - (-60f64).to_radians()).abs() < 1e-9 && (app.cam.pitch - 25f64.to_radians()).abs() < 1e-9, "{:?}", app.cam);
        }
        assert!((app.cam.yaw - (-60f64).to_radians()).abs() < 1e-9 && (app.cam.pitch - 25f64.to_radians()).abs() < 1e-9, "{:?}", app.cam);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Status lines say what changed, in words, whatever started the command (#333), and the
    /// document is named after the sample opened (#334).
    #[test]
    fn status_lines_say_what_changed() {
        let mut app = ModelApp::new(Session::default());
        assert_eq!(app.model_name(), "Untitled");
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        assert_eq!(app.model_name(), "Cantilever (sample)");
        let status = |app: &ModelApp| app.status.as_ref().map(|s| s.text.clone()).unwrap_or_default();
        assert_eq!(status(&app), "Opened the Cantilever sample");
        app.run("solve.run", json!({})).unwrap();
        assert_eq!(status(&app), "Solved 2 steps");
        app.run("results.deformed", json!({"step": "Load", "show": true})).unwrap();
        assert_eq!(status(&app), "Showing Load deformed");
        app.run("results.deformed", json!({"step": "Load", "scale": 2})).unwrap();
        assert_eq!(status(&app), "Showing Load deformed ×2");
        app.run("results.contour", json!({"step": "Load", "field": "von_mises"})).unwrap();
        assert_eq!(status(&app), "Showing Von Mises stress on Load");
        app.run("results.modes", json!({"step": "Modes", "mode": 2})).unwrap();
        assert!(status(&app).starts_with("Mode 2: ") && status(&app).ends_with(" Hz"), "{}", status(&app));
        app.run("results.clear", json!({})).unwrap();
        assert_eq!(status(&app), "Back to the model");
        app.run("model.material.create", json!({"name": "Al", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        assert_eq!(status(&app), "Added material Al");
        app.run("model.material.delete", json!({"name": "Al"})).unwrap();
        assert_eq!(status(&app), "Deleted material Al");
    }

    /// Every command run from the UI reports one readable line, never its JSON, in every state
    /// tried (empty, a solved sample); a reply with more to it opens in the context panel.
    #[test]
    fn every_reply_reads_as_one_line() {
        let empty = Session::default();
        let mut solved = ModelApp::new(Session::default());
        solved.open_sample("cantilever");
        solved.wait_loaded();
        solved.run("solve.run", json!({})).unwrap();
        let mut ran = 0;
        for session in [empty, solved.session.clone()] {
            // Every command bare, and the queries that need parameters with ones the sample has.
            let bare = fieldcraft_engine::commands().iter().filter(|c| !c.id.starts_with("file.")).map(|c| (c.id, json!({})));
            let queries = [
                ("model.set.resolve", json!({"name": "root"})),
                ("solve.probe", json!({"step": "Load", "nodes": [1, 2]})),
                ("solve.probe", json!({"step": "Modes", "set": "tip", "mode": 1})),
                ("results.reactions", json!({"step": "Load"})),
            ];
            for (id, params) in bare.chain(queries) {
                let mut app = ModelApp::new(session.clone());
                let r = app.run(id, params);
                let Some(status) = app.status.clone() else { continue };
                let text = status.text;
                assert!(!text.contains('\n'), "{id}: {text}");
                if r.is_ok() {
                    ran += 1;
                    assert!(!text.contains("{\"") && !text.contains("[{") && text.chars().count() <= 200, "{id}: {text}");
                }
            }
        }
        assert!(ran >= 20, "only {ran} commands ran");
    }

    /// Check Model and the Step Table open their readable view in the context panel, which
    /// follows the model as it changes.
    #[test]
    fn check_model_opens_its_diagnostics_in_the_panel() {
        let mut app = ModelApp::new(Session::default());
        app.run("model.check", json!({})).unwrap();
        let status = app.status.clone().unwrap();
        assert_eq!(status.text, "Check Model: 2 errors, 0 warnings (details in the panel)");
        assert!(app.reply.as_ref().is_some_and(|r| r.command == "model.check"));
        assert!(context::has_context(&mut app));
        let corners = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
        app.run("mesh.block", json!({"body": "A", "kind": "hex8", "corners": corners, "divisions": [1, 1, 1]})).unwrap();
        let messages = app.reply.as_ref().map(|r| r.value["diagnostics"].to_string()).unwrap_or_default();
        assert!(!messages.contains("the mesh has no elements"), "the check follows the model: {messages}");
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("model.activity", json!({})).unwrap();
        assert_eq!(app.status.unwrap().text, "Step Table: 2 steps, 1 boundary condition, 1 load (details in the panel)");
    }

    /// Commands run from outside the ribbon bring their tab; queries and Tools commands leave it.
    #[test]
    fn the_ribbon_follows_commands_to_their_stage() {
        let mut app = ModelApp::new(Session::default());
        app.ui.select_tab(Tab::Results);
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        assert_eq!((app.ui.workspace, app.ui.tab), (Workspace::Prepare, Tab::Mesh));
        app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}})).unwrap();
        assert_eq!((app.ui.workspace, app.ui.tab), (Workspace::Model, Tab::Properties));
        // A query and a failed command leave it; so does undo (Tools).
        app.run("mesh.quality", json!({})).ok();
        let _ = app.run("model.load.create", json!({"name": "bad"}));
        app.run("edit.undo", json!({})).unwrap();
        assert_eq!(app.ui.tab, Tab::Properties);
        // A script follows its last command with a stage.
        app.run_script(&json!([{"command": "mesh.node.add", "params": {"label": 2, "xyz": [1, 0, 0]}}, {"command": "model.info", "params": {}}]))
            .unwrap();
        assert_eq!(app.ui.tab, Tab::Mesh);
        // The ribbon's own entries never move it.
        app.ui.select_tab(Tab::Properties);
        app.start("model.units").unwrap();
        assert_eq!(app.ui.tab, Tab::Properties);
    }

    /// Two blocks meshed apart: the Merge Nodes form shows the nodes it would merge before OK,
    /// and Check Duplicates selects what it finds (#234).
    #[test]
    fn node_merge_previews_and_the_duplicate_check_selects() {
        let mut app = ModelApp::new(Session::default());
        for (body, x0) in [("A", 0.0), ("B", 1.0)] {
            let x1 = x0 + 1.0;
            let c = [[x0, 0., 0.], [x1, 0., 0.], [x1, 1., 0.], [x0, 1., 0.], [x0, 0., 1.], [x1, 0., 1.], [x1, 1., 1.], [x0, 1., 1.]];
            app.run("mesh.block", json!({"body": body, "kind": "hex8", "corners": c, "divisions": [1, 1, 1], "merge": 0})).unwrap();
        }
        assert!(app.merge_preview().is_empty(), "no form open");
        app.start("mesh.node.merge").unwrap();
        assert_eq!(app.merge_preview().len(), 8, "four pairs on x = 1");
        let (_, overlay) = app.overlay();
        assert!(!overlay.points.is_empty());
        app.dialog = None;
        app.run("mesh.check.duplicates", json!({})).unwrap();
        assert_eq!(app.pick.selection.len(), 8 + 2, "the coincident nodes and the two crack faces");
        app.run("mesh.node.merge", json!({})).unwrap();
        assert_eq!(app.session.model.nodes().len(), 12);
        app.run("mesh.check.duplicates", json!({})).unwrap();
        assert!(app.pick.selection.is_empty());
    }

    /// #650: the Solve Fluid form shows the run's memory against this machine's bound before
    /// OK, takes caps raised above their defaults, and says what to change when a cap is short.
    #[test]
    fn the_solve_fluid_form_estimates_the_memory_and_takes_raised_caps() {
        let mut app = ModelApp::new(Session::default());
        app.session.limits.memory = 32_000_000_000;
        for (id, p) in [
            ("model.material.create", json!({"name": "Oil", "density": 1e-9, "fluid": {"type": "newtonian", "viscosity": 1e-9}})),
            ("model.fluid_domain.create", json!({"name": "Box", "material": "Oil", "min": [0, 0, 0], "max": [40, 10, 10], "spacing": 1})),
            ("model.step.create", json!({"name": "Flow", "procedure": {"type": "fluid", "domain": "Box", "end_time": 1, "frames": 24}})),
        ] {
            app.run(id, p).unwrap();
        }
        assert!(app.fluid_note().is_none(), "no form open");
        app.start("solve.fluid").unwrap();
        let (note, problem) = app.fluid_note().unwrap();
        assert!(note.starts_with("Memory: 25 frames × 4,000 cells: 2 MB of frames"), "{note}");
        assert!(note.contains("Runs with raised caps may take up to 32.0 GB on this machine"), "{note}");
        assert!(problem.is_none(), "{problem:?}");
        let fields = fieldcraft_engine::find_command("solve.fluid").and_then(|c| c.form).unwrap();
        let set = |app: &mut ModelApp, p: Value| {
            if let Some(d) = app.dialog.as_mut() {
                d.text = p.to_string();
            }
        };
        // A raised frame cap is a valid entry and leaves the run within the bound.
        let raised = json!({"max_frame_bytes": 8_000_000_000u64, "max_cells": 20_000_000});
        assert!(crate::form::messages(fields, raised.as_object().unwrap()).is_empty());
        set(&mut app, raised);
        assert_eq!(app.fluid_note().unwrap().1, None);
        // A short cap: the note says how many frames fit and what to raise the cap to.
        set(&mut app, json!({"max_frame_bytes": 1_000_000}));
        let problem = app.fluid_note().unwrap().1.unwrap();
        assert!(problem.contains("15 frames fit: set the step's `frames` to at most 14, or raise `max_frame_bytes` to 1600000"), "{problem}");
        // Above what this machine allows.
        set(&mut app, json!({"max_frame_bytes": 40_000_000_000u64}));
        let problem = app.fluid_note().unwrap().1.unwrap();
        assert!(problem.contains("this machine allows"), "{problem}");
    }

    /// Relax Nodes: a click places the sphere, Alt + scroll sizes it, the form previews the moved
    /// nodes and the region's quality, and the sphere is drawn translucent (#398).
    #[test]
    fn relax_nodes_places_its_sphere_and_previews() {
        let mut app = ModelApp::new(Session::default());
        let c = [[0., 0., 0.], [4., 0., 0.], [4., 4., 0.], [0., 4., 0.], [0., 0., 4.], [4., 0., 4.], [4., 4., 4.], [0., 4., 4.]];
        app.run("mesh.block", json!({"body": "A", "kind": "hex8", "corners": c, "divisions": [4, 4, 4]})).unwrap();
        // A spike inside.
        let n = node_at(&app, [2.0, 2.0, 2.0]);
        let spike = fieldcraft_engine::units::geometry_point([2.7, 2.6, 2.5]);
        app.session.model.set_node(n, spike).unwrap();
        app.model_changed();
        assert!(!app.relaxing() && app.relax_preview().is_none());
        app.start("mesh.node.relax").unwrap();
        assert!(app.relaxing());
        assert!(app.relax_place(results::tests::mm([2.0, 2.0, 2.0])));
        assert!(app.relax_scale_radius(10.0));
        let note = app.relax_note().unwrap();
        assert!(note.starts_with("Preview:"), "{note}");
        let p = app.relax_preview().unwrap();
        assert!(p.relaxed.positions.contains_key(&n) && p.relaxed.after.min_quality > p.relaxed.before.min_quality);
        let (_, overlay) = app.overlay();
        assert!(!overlay.fills.is_empty() && !overlay.lines.is_empty(), "the sphere and the preview are drawn");
        // OK runs it as the form says.
        let text = app.dialog.as_ref().unwrap().text.clone();
        app.dialog = None;
        app.run("mesh.node.relax", serde_json::from_str(&text).unwrap()).unwrap();
        assert_ne!(app.session.model.node(n).unwrap().xyz, spike);
    }

    /// Move Nodes: the form previews the move (the region's quality before and now, the ghost and
    /// the coloured elements drawn), Alt + scroll sizes the radius, a drag's release runs the
    /// command in one undo step and keeps the form for the next drag, Esc cancels a drag.
    #[test]
    fn move_nodes_previews_and_runs_on_release() {
        let mut app = ModelApp::new(Session::default());
        let c = [[0., 0., 0.], [4., 0., 0.], [4., 4., 0.], [0., 4., 0.], [0., 0., 4.], [4., 0., 4.], [4., 4., 4.], [0., 4., 4.]];
        app.run("mesh.block", json!({"body": "A", "kind": "hex8", "corners": c, "divisions": [4, 4, 4]})).unwrap();
        let n = node_at(&app, [2.0, 2.0, 2.0]);
        assert!(!app.moving_nodes() && app.node_move_preview().is_none());
        app.start(node_move::COMMAND).unwrap();
        assert!(app.moving_nodes());
        if let Some(d) = app.dialog.as_mut() {
            d.text = json!({"nodes": [n], "by": [0.4, 0.2, 0.0], "on_geometry": true, "guard": "no_inversion"}).to_string();
        }
        assert!(app.node_move_scale_radius(1.0), "a radius to start from");
        assert!(app.node_move_scale_radius(5.0), "Alt + scroll grows it");
        let note = app.node_move_note().unwrap();
        assert!(note.starts_with("Preview:") && note.contains("worse"), "{note}");
        let p = app.node_move_preview().unwrap();
        let m = p.moved.as_ref().unwrap();
        assert!(m.positions.contains_key(&n) && m.moved > 1 && p.sphere.is_some());
        let (_, overlay) = app.overlay();
        assert!(
            !overlay.fills.is_empty() && !overlay.lines.is_empty() && !overlay.on_top.is_empty(),
            "sphere, ghost and moved edges, the grabbed node"
        );
        // A changed radius waits out the debounce (the sphere drawn meanwhile), then prepares.
        assert!(app.node_move_scale_radius(1.5));
        let waiting = app.node_move_preview().unwrap();
        assert!(app.node_move_pending() && waiting.moved.is_none() && waiting.sphere.is_some());
        app.now += node_move::DEBOUNCE * 2.0;
        assert!(app.node_move_preview().unwrap().moved.is_some() && !app.node_move_pending());
        // Esc during a drag: nothing moves, the form stays.
        let start = app.session.model.clone();
        app.node_drag = Some(node_move::Drag { start: results::tests::mm([2.0, 2.0, 2.0]), normal: Vec3::Z });
        app.node_move_cancel();
        assert!(app.node_drag.is_none() && app.moving_nodes() && app.session.model == start);
        assert!(app.dialog.as_ref().unwrap().text.contains("\"by\":[0.0,0.0,0.0]"));
        // A drag's release runs the move.
        if let Some(d) = app.dialog.as_mut() {
            let mut p: Value = serde_json::from_str(&d.text).unwrap();
            p["by"] = json!([0.4, 0.2, 0.0]);
            d.text = p.to_string();
        }
        app.node_drag = Some(node_move::Drag { start: results::tests::mm([2.0, 2.0, 2.0]), normal: Vec3::Z });
        app.node_move_release().unwrap();
        // Moved by the form's 0.4 and 0.2 mm (the model is in metres).
        let moved = app.session.model.node(n).unwrap().xyz;
        assert!(Vec3::from(moved).dist(results::tests::mm([2.4, 2.2, 2.0])) < 1e-15, "{moved:?}");
        assert!(app.moving_nodes() && app.node_drag.is_none(), "the form stays open for the next drag");
        app.run("edit.undo", json!({})).unwrap();
        assert_eq!(app.session.model, start);
    }

    /// Relax Nodes by face (#471): Use selection puts the picked design faces in the form in
    /// place of a sphere, clicks then pick faces, and OK relaxes them.
    #[test]
    fn relax_nodes_takes_picked_faces() {
        let mut app = ModelApp::new(Session::default());
        let dir = std::env::temp_dir().join(format!("fieldcraft-relax-faces-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        std::fs::write(&design, fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap()).unwrap();
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        app.run("mesh.tet.generate", json!({"size": 4.0})).unwrap();
        let body = app.session.design.as_ref().unwrap().bodies[0].clone();
        let face = body.faces()[0].name.clone();
        app.start("mesh.node.relax").unwrap();
        assert!(app.relax_place(Vec3 { x: 1.0, y: 1.0, z: 1.0 }) && !app.relax_by_faces());
        let picked = fieldcraft_model::geometry::GeometryRef::face(body.name.clone(), face.clone()).unwrap();
        app.select([selection::Item::Geometry(picked)], selection::Mode::Replace);
        app.use_selection_in_dialog().unwrap();
        let params: Value = serde_json::from_str(&app.dialog.as_ref().unwrap().text).unwrap();
        assert_eq!(params["faces"], json!([{"body": body.name, "face": face}]), "{params}");
        assert!(params.get("center").is_none() && params.get("radius").is_none(), "{params}");
        assert!(app.relax_by_faces() && app.relax_preview().is_some());
        app.dialog = None;
        app.run("mesh.node.relax", params).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Mesh Quality reads as a sentence in the status bar, not as its JSON report.
    #[test]
    fn mesh_quality_reports_in_a_line() {
        let mut app = ModelApp::new(Session::default());
        // Nothing to measure: the command needs a mesh and says so.
        assert!(app.run("mesh.quality", json!({})).is_err());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        app.run("mesh.quality", json!({})).unwrap();
        let status = app.status.clone().unwrap();
        assert!(status.text.starts_with("Mesh Quality: 226 elements, all valid; scaled Jacobian min "), "{}", status.text);
        assert!(!status.error && !status.text.contains('{'));
    }

    /// Mark Bad Elements draws the marked elements' faces in the error colour, and only while
    /// they are marked.
    #[test]
    fn marked_elements_draw_in_the_error_colour() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        let error = theme::rgb(theme::Tokens::get(app.ui.dark).error);
        let red = |app: &mut ModelApp| app.scene().1.tris.iter().filter(|t| t.color == error).count();
        assert_eq!(red(&mut app), 0);
        let some = app.session.model.elements().keys().take(3).copied().collect::<Vec<_>>();
        app.run("mesh.quality.mark", json!({"elements": some, "limits": {"scaled_jacobian": 1.5}})).unwrap();
        let marked = red(&mut app);
        assert!(marked > 0 && marked < app.scene().1.tris.len(), "{marked}");
        app.run("mesh.quality.mark", json!({})).unwrap();
        assert_eq!(red(&mut app), 0);
    }

    /// Check Geometry's findings draw in the error colour (faces, and a dot at each over
    /// everything) until a check with `mark: false`.
    #[test]
    fn geometry_check_findings_draw_in_the_error_colour() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-check-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        let bytes = fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap();
        std::fs::write(&design, bytes).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        let error = theme::rgb(theme::Tokens::get(app.ui.dark).error);
        let red = |app: &mut ModelApp| {
            let s = app.scene().1;
            (s.tris.iter().filter(|t| t.color == error).count(), s.on_top.iter().filter(|l| l.color == error && l.a == l.b).count())
        };
        assert_eq!(red(&mut app), (0, 0));
        // Against a 1 m element every face of the plate is small.
        let r = app.run("geometry.check", json!({"size": 1000})).unwrap();
        let (faces, dots) = red(&mut app);
        assert!(faces > 0, "{r}");
        assert_eq!(dots, r["findings"].as_array().unwrap().len());
        app.run("geometry.check", json!({"size": 1000, "mark": false})).unwrap();
        assert_eq!(red(&mut app), (0, 0));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn starting_a_command_opens_its_form_and_running_it_closes_the_form() {
        let mut app = ModelApp::new(Session::default());
        app.start("model.info").unwrap();
        assert!(app.dialog.is_none());
        app.start("mesh.node.add").unwrap();
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("mesh.node.add"));
        app.dialog.as_mut().unwrap().text = "not json".into();
        assert!(app.run_dialog().is_err());
        assert!(app.dialog.as_ref().is_some_and(|d| d.error.is_some()));
        app.dialog.as_mut().unwrap().text = r#"{"label": 7, "xyz": [1, 2, 3]}"#.into();
        app.run_dialog().unwrap();
        assert!(app.dialog.is_none());
        assert!(app.session.model.node(7).is_some());
        assert!(app.start("no.such").is_err());
    }

    /// A parameter the command refuses is named by its field's label in the form and the status
    /// bar; the command itself (agents) still names the key.
    #[test]
    fn a_refused_parameter_is_named_by_its_label() {
        let mut app = ModelApp::new(Session::default());
        app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}})).unwrap();
        app.start("model.material.plasticity").unwrap();
        app.dialog.as_mut().unwrap().text = r#"{"name": 7}"#.into();
        let e = app.run_dialog().unwrap_err();
        assert_eq!(e, "Material: must be a name");
        let d = app.dialog.as_ref().unwrap();
        assert_eq!(d.refused, Some(("name".to_string(), "Material: must be a name".to_string())));
        assert!(app.status.as_ref().is_some_and(|s| s.text == "Plasticity: Material: must be a name"));
        let raw = app.session.execute("model.material.plasticity", &json!({"name": 7})).unwrap_err();
        assert_eq!(raw.to_string(), "`name` must be a name");
    }

    /// A required reference with one candidate starts picked (the browser's selection wins); with
    /// none, "Create …" opens the providing command over the form, which comes back with the new
    /// item picked (or as it was, when that command is cancelled).
    #[test]
    fn reference_fields_pick_the_only_candidate_and_create_a_missing_one() {
        use fieldcraft_engine::schema::Collection;
        let params = |app: &ModelApp| serde_json::from_str::<Value>(&app.dialog.as_ref().unwrap().text).unwrap();
        let steel = json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}});
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.node.add", json!({"label": 1, "xyz": [0, 0, 0]})).unwrap();
        app.run("model.set.create", json!({"name": "all", "nodes": [1]})).unwrap();
        // (Opened directly: from the toolbar, Section waits for a material.)
        app.open_form("model.section.create", json!({}));
        let opened = params(&app);
        assert_eq!(opened["set"], "all");
        assert!(opened.get("material").is_none());

        // Cancelled: back to the section form as it was.
        app.create_from_form(Collection::Material, "material".into());
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.material.create"));
        app.close_dialog();
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.section.create"));
        assert_eq!(params(&app), opened);

        // Made: back with it picked.
        app.create_from_form(Collection::Material, "material".into());
        app.dialog.as_mut().unwrap().text = steel.to_string();
        app.run_dialog().unwrap();
        let d = app.dialog.as_ref().unwrap();
        assert_eq!(d.command, "model.section.create");
        assert!(d.back.is_none());
        assert_eq!(params(&app)["material"], "Steel");
        app.close_dialog();
        assert!(app.dialog.is_none());
        // What a form made is selected (#360); the rest of the test starts with nothing selected.
        assert_eq!(app.selected, Some((browser::Def::Material, "Steel".into())));
        app.selected = None;

        // One material: Plasticity starts with it.
        app.start("model.material.plasticity").unwrap();
        assert_eq!(params(&app)["name"], "Steel");
        // Two: the one selected in the browser, else nothing.
        app.run("model.material.create", json!({"name": "Al", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        app.start("model.material.plasticity").unwrap();
        assert!(params(&app).get("name").is_none());
        app.selected = Some((browser::Def::Material, "Steel".into()));
        app.start("model.material.plasticity").unwrap();
        assert_eq!(params(&app)["name"], "Steel");
    }

    #[test]
    fn toolbar_entries_start_their_command_and_planned_ones_do_nothing() {
        let mut app = ModelApp::new(Session::default());
        app.start("model.material.create").unwrap();
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.material.create"));
        app.dialog = None;
        let journal = app.session.journal.len();
        let err = app.start("mesh.order.quadratic").unwrap_err();
        assert_eq!(err, "Quadratic is coming soon (M2)");
        assert!(app.dialog.is_none());
        assert_eq!(app.session.journal.len(), journal);
        assert_eq!(app.status.as_ref().map(|s| s.error), Some(false));
        // An entry whose command misses what it needs opens nothing and says what is missing.
        let err = app.start("deck.write.abaqus").unwrap_err();
        assert_eq!(err, "deck.write cannot run yet. Needs a mesh — Mesh › Tet Mesh");
        assert!(app.dialog.is_none());
        let status = "Write Solver Deck cannot run yet. Needs a mesh — Mesh › Tet Mesh";
        assert_eq!(app.status.as_ref().map(|s| (s.error, s.text.as_str())), Some((true, status)));
        // The panel lists what is missing (#323), until the command can run.
        assert_eq!(app.refused.as_ref().map(|u| u.command), Some("deck.write"));
        // With a model, a deck entry opens deck.write's form with its format filled in.
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        app.start("deck.write.abaqus").unwrap();
        app.logic(&egui::Context::default());
        assert!(app.refused.is_none());
        let d = app.dialog.as_ref().unwrap();
        assert_eq!((d.command, d.text.as_str()), ("deck.write", r#"{"format":"abaqus"}"#));
    }

    /// The ribbon, the command search and Solve grey what the model cannot run yet, and say why.
    #[test]
    fn entries_missing_what_they_need_are_greyed_with_the_reason() {
        let mut app = ModelApp::new(Session::default());
        let entry =
            |app: &ModelApp, id: &str| crate::toolbar::describe(app)["entries"].as_array().unwrap().iter().find(|e| e["id"] == id).unwrap().clone();
        let plasticity = entry(&app, "model.material.plasticity");
        assert_eq!(plasticity["enabled"], false);
        assert_eq!(plasticity["unavailable"]["reason"], "Needs a material — Properties › Material");
        assert_eq!(plasticity["unavailable"]["missing"][0]["provider"], "model.material.create");
        assert_eq!(entry(&app, "model.material.create")["enabled"], true);
        assert_eq!(entry(&app, "mesh.order.quadratic")["enabled"], false, "coming soon");
        let solve = app.unavailable("solve.run").unwrap();
        assert_eq!(solve.lines(), "Needs a mesh — Mesh › Tet Mesh\nNeeds a step — Steps › Step");

        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        assert_eq!(entry(&app, "model.material.plasticity")["enabled"], true);
        assert!(app.unavailable("solve.run").is_none());
        assert_eq!(entry(&app, "results.contour")["unavailable"]["reason"], "Needs a solved step — Solve › Solve");
    }

    #[test]
    fn the_sample_fills_the_view_and_renders() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        let d = app.display().clone();
        assert!(!d.tris.is_empty() && !d.edges.is_empty());
        // Framed when the viewport first draws.
        draw(&mut app, &egui::Context::default(), 1280.0, 800.0);
        assert!(app.framed(), "{:?}", app.cam);
        let png = app.render_png(64, 48, 1.0).unwrap().0;
        assert!(png.starts_with(&[0x89, b'P', b'N', b'G']));
    }

    /// The window and the headless pictures build one scene: the same session gives the same
    /// faces, lines and points from the app's viewport and from `engine::view` (the light theme
    /// is `Palette::LIGHT`), meshed or not.
    #[test]
    fn the_window_and_headless_pictures_build_one_scene() {
        use fieldcraft_engine::view::{self, Palette};
        assert_eq!(theme::Tokens::light().palette(), Palette::LIGHT);
        let mut app = ModelApp::new(Session::default());
        app.ui.dark = false;
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        draw(&mut app, &egui::Context::default(), 1280.0, 800.0);
        let headless =
            view::frame(&app.session, &fieldcraft_render::Orientation::of(fieldcraft_render::StandardView::Iso), 320, 200, &Palette::LIGHT, None);
        let window = app.scene().1.clone();
        assert!(!window.tris.is_empty() && !window.lines.is_empty());
        assert_eq!(window.tris, headless.scene.tris);
        assert_eq!(window.lines, headless.scene.lines);
        assert_eq!(window.points, headless.scene.points);
        assert_eq!(window.fills, headless.scene.fills);
        assert_eq!(window.background, headless.scene.background);
    }

    /// `ui.render` draws with the viewport's pipelines on the window's device when there is one
    /// (here a device the `FIELDCRAFT_GPU` policy opens, as the app's), says so, and draws on
    /// the CPU with the reason when the GPU could not be set up.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn renders_on_the_windows_gpu() {
        use fieldcraft_gpu::{Renderer, capture::Capturer, select};
        let _turn = crate::gpu_turn();
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::script()).unwrap();
        draw(&mut app, &egui::Context::default(), 1280.0, 800.0);
        let size = |png: &[u8]| (u32::from_be_bytes(png[16..20].try_into().unwrap()), u32::from_be_bytes(png[20..24].try_into().unwrap()));
        match Capturer::open(select::Policy::from_env().unwrap()) {
            Ok(Some(c)) => {
                app.capturer = Ok(Some(c));
                let (png, drawn) = app.render_png(160, 100, 2.0).unwrap();
                assert!(matches!(drawn, Renderer::Gpu { .. }), "{drawn:?}");
                assert_eq!(size(&png), (160, 100));
                let (gpu, _) = app.render_png(160, 100, 1.0).unwrap();
                assert_ne!(png, gpu, "the scale changes the marks");
                assert!(app.render_png(0, 100, 1.0).is_err());
            }
            Ok(None) | Err(_) if select::gpu_required() => panic!("FIELDCRAFT_REQUIRE_GPU=1 but no GPU"),
            Ok(None) | Err(_) => eprintln!("no GPU adapter: the GPU half skipped"),
        }
        app.capturer = Err("no device".into());
        let (png, drawn) = app.render_png(160, 100, 1.0).unwrap();
        assert_eq!(drawn, Renderer::Cpu { why: Some("no device".into()) });
        assert_eq!(size(&png), (160, 100));
    }

    #[test]
    fn the_open_design_and_imported_surfaces_are_drawn() {
        let dir = std::env::temp_dir().join(format!("fieldcraft-ui-design-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let design = dir.join("plate.solidcraft");
        let bytes = fieldcraft_geometry::sample::plate_design(&fieldcraft_geometry::sample::Plate::default()).unwrap();
        std::fs::write(&design, bytes).unwrap();
        let mut app = ModelApp::new(Session::default());
        app.run("geometry.open", json!({"path": design.to_string_lossy()})).unwrap();
        let d = app.display().clone();
        assert!(!d.surface_tris.is_empty() && !d.surface_edges.is_empty(), "the design's bodies are drawn");
        assert!(d.tris.is_empty());
        // Without elements, geometry draws as faces with their boundaries; the camera frames it.
        let t = theme::Tokens::dark();
        let s = scene::build(&d, Style::ShadedEdges, scene::Show::DEFAULT, &t);
        assert_eq!(s.tris.len(), d.surface_tris.len());
        assert!(s.tris.iter().all(|tri| tri.color == theme::rgb(t.geometry)));
        assert!(!d.bounds.is_empty());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The sample cantilever in an app, solved (static step `Load`, frequency step `Modes`).
    fn solved_app() -> ModelApp {
        ModelApp::new(results::tests::solved_cantilever())
    }

    fn node_at(app: &ModelApp, xyz: [f64; 3]) -> u32 {
        results::tests::node_at(&app.session, xyz)
    }

    #[test]
    fn setup_symbols_show_with_their_toggles_and_fade_over_results() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        let t = theme::Tokens::get(app.ui.dark);
        let count = |app: &mut ModelApp, c: egui::Color32| app.symbols().1.lines.iter().filter(|x| x.color == theme::rgb(c)).count();
        assert!(count(&mut app, t.glyph_constraint) > 0, "the clamp's triangles");
        assert!(count(&mut app, t.glyph_load) > 0, "the push's arrows");
        // Wireframes, in a layer of their own.
        assert!(app.symbols().1.tris.is_empty());
        assert!(!app.scene().1.lines.iter().any(|l| l.color == theme::rgb(t.glyph_constraint)));
        // The symbols reach beyond the model: the clipping radius holds them.
        let radius = app.scene().1.radius;
        app.ui.show_loads = false;
        assert_eq!(count(&mut app, t.glyph_constraint), 0);
        assert!(app.scene().1.radius <= radius);
        app.ui.show_loads = true;
        // Selecting the clamp in the browser highlights its symbols.
        app.selected = Some((browser::Def::BoundaryCondition, "Clamp".into()));
        assert_eq!(count(&mut app, t.glyph_constraint), 0);
        assert!(count(&mut app, t.pick_selected) > 0);
        assert!(count(&mut app, t.glyph_load) > 0);
        app.selected = None;
        // Over a results plot they are off until asked for (#327), then muted.
        app.run("solve.static", json!({})).unwrap();
        app.run("results.deformed", json!({"step": "Load"})).unwrap();
        assert!(app.symbols().1.lines.is_empty());
        app.ui.result_symbols = true;
        assert_eq!(count(&mut app, t.glyph_constraint), 0);
        let muted = t.glyph_colors().over_results(theme::rgb(t.viewport_bottom));
        assert!(app.symbols().1.lines.iter().any(|x| x.color == muted.constraint));
        assert!(app.glyph_notes().is_empty());
    }

    /// Zooming resizes the symbols (a layer of their own) and leaves the mesh's scene alone; the
    /// symbols measure the size setting on screen at any zoom.
    #[test]
    fn zooming_rebuilds_the_symbols_only_and_keeps_their_size_on_screen() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        app.set_view(StandardView::Front);
        let h = 600.0;
        // The length of the longest constraint line on screen: a triangle's side, about its
        // length (a symbol unit) long.
        let measure = |app: &mut ModelApp| {
            let cam = app.cam;
            let lines = app.symbols_for(h, 1.0).1.lines.clone();
            let t = theme::Tokens::get(app.ui.dark);
            let at = |p: Vec3| cam.to_screen(p, 800.0, h, 1e6).unwrap();
            lines
                .iter()
                .filter(|l| l.color == theme::rgb(t.glyph_constraint))
                .map(|l| {
                    let (a, b) = (at(l.a), at(l.b));
                    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
                })
                .fold(0.0, f64::max)
        };
        let scene = app.scene_key();
        let (symbols, _) = app.symbols_for(h, 1.0);
        let mut sizes = Vec::new();
        for distance in [40.0, 150.0, 900.0, 5000.0] {
            app.cam.distance = distance;
            sizes.push(measure(&mut app));
        }
        assert_eq!(app.scene_key(), scene, "the camera is not in the scene's key");
        assert_ne!(app.symbols_for(h, 1.0).0, symbols, "zooming resizes the symbols");
        let (lo, hi) = sizes.iter().fold((f64::MAX, 0.0f64), |(lo, hi), s| (lo.min(*s), hi.max(*s)));
        assert!(lo > 10.0 && hi / lo < 1.2, "{sizes:?}");
        // The size setting scales them.
        let before = measure(&mut app);
        app.ui.symbol_size = 2.0 * glyphs::DEFAULT_SIZE;
        let after = measure(&mut app);
        assert!((after / before - 2.0).abs() < 0.2, "{before} -> {after}");
        assert_eq!(app.scene_key(), scene);
    }

    #[test]
    fn solving_shows_the_result_on_the_results_tab() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::cantilever()).unwrap();
        app.start("solve.static").unwrap();
        app.run_dialog().unwrap();
        let view = app.session.results_view.clone().unwrap();
        // A static step first shows its contours on the undeformed shape (#77).
        assert_eq!((view.step.as_str(), view.contour, view.deformed), ("Load", Some(fieldcraft_engine::results::Field::VonMises), false));
        app.logic(&egui::Context::default());
        assert_eq!(app.ui.tab, Tab::Results);
        app.run("solve.frequency", json!({})).unwrap();
        app.run("results.deformed", json!({"step": "Modes"})).unwrap();
        let view = app.session.results_view.clone().unwrap();
        assert_eq!((view.step.as_str(), view.mode, view.contour), ("Modes", 1, Some(fieldcraft_engine::results::Field::Displacement)));
        // The plot draws contour bands: more triangles than the model's faces.
        let faces = app.display().tris.len();
        assert!(app.scene().1.tris.len() > faces);
        app.run("results.clear", json!({})).unwrap();
        assert!(app.plot().is_none());
        // The model's faces again (without the setup's symbols).
        app.ui.show_loads = false;
        assert_eq!(app.scene().1.tris.len(), faces);
    }

    /// One frame of the window at `w` × `h` points.
    fn draw(app: &mut ModelApp, ctx: &egui::Context, w: f32, h: f32) {
        let raw = egui::RawInput { screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(w, h))), ..Default::default() };
        let mut out = ctx.run_ui(raw, |ui| {
            app.logic(ui.ctx());
            app.ui(ui);
        });
        out.textures_delta.clear();
    }

    /// XY Plot from the ribbon opens the plot in the results panel, which runs results.xy for
    /// its quantity and keeps the series; a nodal field waits for a probed node (#94).
    #[test]
    fn the_xy_plot_runs_results_xy() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(Session::default());
        app.open_sample("impact");
        app.wait_loaded();
        app.run("solve.run", json!({})).unwrap();
        app.start("results.xy").unwrap();
        assert!(app.dialog.is_none(), "the plot opens in the panel, not a form");
        draw(&mut app, &ctx, 1280.0, 800.0);
        let s = app.xy.as_ref().and_then(|x| x.series.clone()).expect("the panel ran results.xy");
        assert_eq!((s["step"].as_str(), s["x"]["axis"].as_str()), (Some("Impact"), Some("time")));
        assert!(s["points"].as_u64().unwrap_or(0) > 1, "{s}");
        // A nodal field needs a node: the probed one.
        app.xy.as_mut().unwrap().y = Some(json!({"field": "velocity"}));
        draw(&mut app, &ctx, 1280.0, 800.0);
        assert!(app.xy.as_ref().unwrap().series.is_none());
        app.probe_node(1).unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        let s = app.xy.as_ref().and_then(|x| x.series.clone()).expect("read at the probed node");
        assert_eq!(s["frames"].as_array().map(Vec::len), Some(11), "{s}");
        // A linear static step offers nothing; a frequency step its frequencies.
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("solve.run", json!({})).unwrap();
        app.run("results.modes", json!({"step": "Modes"})).unwrap();
        app.start("results.xy").unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        let s = app.xy.as_ref().and_then(|x| x.series.clone()).expect("the frequencies");
        assert_eq!(s["x"]["axis"], "mode");
    }

    /// A nonlinear static step's XY plot offers the force–travel curve of each boundary-condition
    /// set (#204): picked, the panel runs results.xy with the set's travel as x.
    #[test]
    fn the_xy_plot_offers_force_travel_curves() {
        let ctx = egui::Context::default();
        let mut app = ModelApp::new(Session::default());
        app.open_sample("rubber-seal");
        app.wait_loaded();
        app.run("solve.run", json!({})).unwrap();
        app.start("results.xy").unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        let result = app.session.latest_result("Squeeze").unwrap().0.clone();
        let options = crate::results_panel::xy_options(&app.session, &result);
        let (_, y, x) = options.iter().find(|(l, ..)| l == "Force–travel Z: top").expect("the top's force–travel curve").clone();
        let xy = app.xy.as_mut().unwrap();
        (xy.y, xy.x) = (Some(y), x);
        draw(&mut app, &ctx, 1280.0, 800.0);
        let s = app.xy.as_ref().and_then(|x| x.series.clone()).expect("the panel ran results.xy");
        assert_eq!(s["x"]["axis"], "displacement", "{s}");
        assert!(s["points"].as_u64().unwrap_or(0) > 5, "{s}");
    }

    fn is_iso(cam: &Camera) -> bool {
        let (yaw, pitch) = StandardView::Iso.angles();
        (cam.yaw - yaw).abs() < 1e-9 && (cam.pitch - pitch).abs() < 1e-9
    }

    /// Every registered sample is on the Start card and opens with a click, at 1280 × 800 and in a
    /// small window where the grid scrolls (#245). Unsaved work (a material, no geometry) keeps the
    /// Start card up and turns the click into the question naming that sample, so nothing loads.
    #[test]
    fn every_sample_is_reachable_from_the_start_card() {
        for (w, h) in [(1280.0, 800.0), (900.0, 560.0)] {
            let mut scrolled = false;
            for sample in fieldcraft_engine::sample::SAMPLES {
                let mut app = ModelApp::new(Session::default());
                app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210000, "poisson": 0.3}}))
                    .unwrap();
                let ctx = egui::Context::default();
                let mut time = 0.0;
                let mut frame = |app: &mut ModelApp, events: Vec<egui::Event>| {
                    time += 0.1;
                    let raw = egui::RawInput {
                        screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(w, h))),
                        time: Some(time),
                        events,
                        ..Default::default()
                    };
                    let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
                    out.textures_delta.clear();
                };
                frame(&mut app, vec![]);
                frame(&mut app, vec![]);
                let card = |ctx: &egui::Context| ctx.read_response(samples::card_id(sample.name)).map(|r| (r.rect, r.interact_rect));
                // Wheel down over the grid until the card is in view: wholly, or (when the grid is
                // shorter than a card) most of it.
                let in_view = |c: Option<(egui::Rect, egui::Rect)>| {
                    c.is_some_and(|(r, i)| i.is_positive() && i.height() >= (r.height() * 0.6).min(r.height() - 0.5))
                };
                // Over the first card, which is in view before any scrolling.
                let first = ctx
                    .read_response(samples::card_id(samples::ordered()[0].name))
                    .map_or(egui::pos2(w / 2.0, h / 2.0), |r| r.interact_rect.center());
                for _ in 0..40 {
                    if in_view(card(&ctx)) {
                        break;
                    }
                    let wheel = egui::Event::MouseWheel {
                        unit: egui::MouseWheelUnit::Point,
                        delta: egui::vec2(0.0, -60.0),
                        phase: egui::TouchPhase::Move,
                        modifiers: Default::default(),
                    };
                    frame(&mut app, vec![egui::Event::PointerMoved(first), wheel]);
                    frame(&mut app, vec![]);
                    scrolled = true;
                }
                assert!(in_view(card(&ctx)), "{w}×{h}: {} is never in view: {:?}", sample.name, card(&ctx));
                let at = card(&ctx).map_or(egui::Pos2::ZERO, |(_, i)| i.center());
                let button =
                    |pressed| egui::Event::PointerButton { pos: at, button: egui::PointerButton::Primary, pressed, modifiers: Default::default() };
                frame(&mut app, vec![egui::Event::PointerMoved(at)]);
                frame(&mut app, vec![button(true)]);
                frame(&mut app, vec![button(false)]);
                assert_eq!(app.prompt, Some(files::Prompt::Unsaved(files::Replace::Sample(sample.name))), "{w}×{h}: clicking {}", sample.name);
            }
            // The small window cannot show them all at once.
            if w < 1000.0 {
                assert!(scrolled, "{w}×{h}");
            }
        }
    }

    /// Fit and home never frame datums (#295); a hidden datum selected in the browser is drawn,
    /// highlighted, until the selection moves on.
    #[test]
    fn datums_are_not_framed_and_a_selected_one_shows() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        app.run("file.sample", json!({"name": "cantilever"})).unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        app.home_view();
        let home = app.cam;
        app.run("geometry.datum", json!({"name": "Far", "by": {"type": "offset", "plane": "xy", "offset": 5000}})).unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        app.home_view();
        // What framing sees is the model alone; the view stays as close as before (the overlays
        // may differ a little after a command).
        assert_eq!(app.drawn_bounds(), app.display.1.bounds, "a datum far off is not framed");
        assert!(app.cam.distance < home.distance * 1.2, "{:?} vs {:?}", app.cam, home);
        app.fit_view();
        assert!(app.cam.distance < home.distance * 1.2);
        let highlight = theme::Tokens::get(app.ui.dark).glyph_colors().highlight;
        let lit = |app: &mut ModelApp| app.scene().1.on_top.iter().filter(|l| l.color == highlight).count();
        app.run("view.hide", json!({"items": [{"kind": "datum", "name": "Far"}]})).unwrap();
        assert_eq!(lit(&mut app), 0);
        app.selected = Some((browser::Def::Datum, "Far".into()));
        assert!(lit(&mut app) > 0, "the selected datum is drawn though hidden");
        app.selected = None;
        assert_eq!(lit(&mut app), 0);
    }

    /// Opening a sample, then another, frames each in the home view; an edit keeps the camera
    /// (#233, #223).
    #[test]
    fn opening_another_sample_frames_it_and_an_edit_keeps_the_camera() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        app.run("file.sample", json!({"name": "crank-arm"})).unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        assert!(is_iso(&app.cam) && app.framed(), "{:?}", app.cam);
        // Looked at from elsewhere, an edit keeps that camera.
        app.cam.orbit(120.0, -40.0);
        app.cam.zoom_at(0.2, None);
        let turned = app.cam;
        app.run("model.material.create", json!({"name": "Al", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}})).unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        assert_eq!(app.cam, turned);
        // Another sample: a fresh start's home view, all of it in sight.
        app.run("file.sample", json!({"name": "bracket-hex"})).unwrap();
        draw(&mut app, &ctx, 1280.0, 800.0);
        assert!(is_iso(&app.cam) && app.framed(), "{:?}", app.cam);
    }

    /// After a solve the model is in sight, clear of the legend, at the window sizes people use;
    /// a fit leaves the legend free with the deformed shape too (#223, #76).
    #[test]
    fn results_come_up_framed_clear_of_the_legend() {
        for (w, h) in [(1280.0, 800.0), (700.0, 1000.0)] {
            let mut app = ModelApp::new(Session::default());
            let ctx = egui::Context::default();
            app.run("file.sample", json!({"name": "crank-arm"})).unwrap();
            draw(&mut app, &ctx, w, h);
            // Zoomed in on a corner before solving: the results bring it all back.
            let corner = app.display().bounds.max;
            app.cam.zoom_at(0.3, Some(corner));
            app.run("solve.run", json!({})).unwrap();
            draw(&mut app, &ctx, w, h);
            draw(&mut app, &ctx, w, h);
            let legend = *app.layout.get("legend").expect("the legend is placed");
            assert!(app.framed(), "{w}×{h}: {:?}", app.cam);
            app.run("results.deformed", json!({"show": true, "scale": 200})).unwrap();
            draw(&mut app, &ctx, w, h);
            app.fit_view();
            assert!(app.framed(), "{w}×{h} deformed");
            let view = app.viewport.rect.unwrap();
            let b = app.drawn_bounds();
            for (k, c) in b.corners().into_iter().enumerate() {
                let at = crate::viewport::to_screen(&app, c).unwrap();
                assert!(!legend.contains(at) && view.contains(at), "{w}×{h}: corner {k} at {at:?}, legend {legend:?}");
            }
        }
    }

    /// At the smallest window the legend goes below the selection bar instead of under it; a wide
    /// one keeps it at the top (#559).
    #[test]
    fn the_legend_stays_clear_of_the_selection_bar() {
        for (w, h, at_top) in [(960.0, 600.0, false), (1400.0, 900.0, true)] {
            let mut app = ModelApp::new(Session::default());
            let ctx = egui::Context::default();
            app.run("file.sample", json!({"name": "cantilever"})).unwrap();
            app.run("solve.run", json!({})).unwrap();
            for _ in 0..3 {
                draw(&mut app, &ctx, w, h);
            }
            let (legend, bar) = (region(&app, "legend"), region(&app, "selection_bar"));
            assert!(!legend.intersects(bar), "{w}×{h}: legend {legend:?}, selection bar {bar:?}");
            let view = app.viewport.rect.unwrap();
            assert!(view.contains_rect(legend), "{w}×{h}: legend {legend:?} outside the viewport {view:?}");
            assert_eq!(legend.top() == view.top() + 12.0, at_top, "{w}×{h}: legend {legend:?}, viewport {view:?}");
        }
    }

    /// Every standard view frames the whole model: the cantilever's front view keeps both ends
    /// (#147), in orthographic and perspective.
    #[test]
    fn standard_views_frame_the_whole_model() {
        for perspective in [false, true] {
            let mut app = ModelApp::new(Session::default());
            app.ui.perspective = perspective;
            let ctx = egui::Context::default();
            app.run("file.sample", json!({"name": "cantilever"})).unwrap();
            draw(&mut app, &ctx, 1280.0, 800.0);
            for v in StandardView::ALL {
                app.set_view(v);
                assert!(app.framed(), "{} view, perspective {perspective}: {:?}", v.label(), app.cam);
            }
        }
    }

    /// Rendered on the CPU from the front, undeformed, with Z displacement contours: the free end
    /// (lowest band) is blue, the clamped end (highest band) red.
    #[test]
    fn contours_render_with_known_colours_at_known_nodes() {
        let mut app = solved_app();
        app.run("results.contour", json!({"step": "Load", "field": "displacement_z", "bands": 6})).unwrap();
        app.run("results.deformed", json!({"show": false})).unwrap();
        app.ui.style = Style::Shaded;
        app.set_view(StandardView::Front);
        app.fit_view();
        let (w, h) = (400, 200);
        let cam = app.cam;
        let radius = app.scene().1.radius;
        let canvas = fieldcraft_render::render(app.scene().1, &cam, w, h).unwrap();
        let pixel = |p: Vec3| {
            let (x, y, _) = cam.to_screen(p, w as f64, h as f64, radius).unwrap();
            let i = (y as usize * w + x as usize) * 4;
            [canvas.rgba[i], canvas.rgba[i + 1], canvas.rgba[i + 2]]
        };
        // Cell centres on the front face (y = 0), one at each end.
        let tip = pixel(results::tests::mm([97.5, 0.0, 12.5]));
        let root = pixel(results::tests::mm([2.5, 0.0, 12.5]));
        assert!(tip[2] > tip[0] + 60, "the free end is blue: {tip:?}");
        assert!(root[0] > root[2] + 60, "the clamped end is red: {root:?}");
    }

    #[test]
    fn probing_a_node_reads_its_values_and_stale_results_say_so() {
        let mut app = solved_app();
        app.run("results.deformed", json!({"step": "Load"})).unwrap();
        // The window has framed its document, as it does on its first frame.
        app.model_changed();
        let tip = node_at(&app, [100.0, 0.0, 0.0]);
        app.probe_node(tip).unwrap();
        let probe = app.probe.clone().unwrap();
        assert_eq!((probe.node, probe.step.as_str()), (tip, "Load"));
        assert!(probe.lines.iter().any(|l| l.starts_with("U ")) && probe.lines.iter().any(|l| l.starts_with("Mises")), "{probe:?}");
        // Units, and round-off shown as 0 (#329).
        assert!(probe.lines.iter().any(|l| l.starts_with("U [mm] ")) && probe.lines.iter().any(|l| l.ends_with(" MPa")), "{probe:?}");
        // Other display units read the probe again in them.
        app.run("model.units", json!({"system": "m-kg-s-N-Pa"})).unwrap();
        let si = app.probe.clone().unwrap();
        assert!(si.node == tip && si.lines.iter().any(|l| l.starts_with("U [m] ")) && si.lines.iter().any(|l| l.ends_with(" Pa")), "{si:?}");
        app.run("model.units", json!({"system": "mm-t-s-N-MPa"})).unwrap();
        assert_eq!(app.probe.as_ref().map(|p| &p.lines), Some(&probe.lines), "and back");
        // A root node carries a reaction.
        app.probe_node(node_at(&app, [0.0, 0.0, 0.0])).unwrap();
        assert!(app.probe.as_ref().unwrap().lines.iter().any(|l| l.starts_with("RF ")));
        // The plot under a click: a ray onto the tip's end face picks a node there.
        let plot = app.plot().unwrap();
        let (node, _) = plot.display.pick(results::tests::mm([1000.0, 0.5, 0.5]), Vec3::new(-1.0, 0.0, 0.0)).unwrap();
        assert_eq!(app.session.model.node(node).unwrap().xyz[0], 0.1);
        // An edit leaves the plot shown but stale; probing it is refused.
        app.run(
            "model.material.edit",
            json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 70000, "poisson": 0.33}, "density": 2.7e-9}),
        )
        .unwrap();
        assert!(app.plot().is_some_and(|p| p.stale));
        let e = app.probe_node(tip).unwrap_err();
        assert!(e.contains("stale") && app.probe.is_none(), "{e}");
    }

    #[test]
    fn explicit_frames_are_drawn_and_animated() {
        let mut app = ModelApp::new(Session::default());
        app.run_script(&fieldcraft_engine::sample::impact()).unwrap();
        app.run("solve.explicit", json!({})).unwrap();
        // The last frame at true scale: the bar has moved back off the wall's side (x < 0 is
        // the wall, so nothing passes it) and carries stress contours.
        let last = app.plot().unwrap();
        assert_eq!(last.frame().map(|(k, _, n)| (k, n)), Some((10, 11)));
        assert!(last.contour.as_ref().is_some_and(|c| c.field == fieldcraft_engine::results::Field::VonMises));
        assert!(last.display.bounds.min.x > -1e-9, "{:?}", last.display.bounds);
        // Frame 0 is the undeformed bar at rest on its initial velocity.
        app.run("results.frame", json!({"frame": 0})).unwrap();
        let first = app.plot().unwrap();
        assert_eq!(first.display.tris, app.display().tris);
        // Animated, the frames run through once per cycle.
        app.run("results.animate", json!({})).unwrap();
        app.now = results::PERIOD * 0.5;
        assert_eq!(app.plot().unwrap().pick.frame, Some(5));
        let image = app.capture(capture::CaptureSize { width: 1200.0, height: 800.0, scale: 1.0 }).unwrap();
        assert_eq!(image.size, [1200, 800]);
        // Probing reads the frame drawn.
        app.run("results.animate", json!({"on": false})).unwrap();
        app.run("results.frame", json!({"frame": 4})).unwrap();
        app.probe_node(node_at(&app, [100.0, 0.0, 0.0])).unwrap();
        let probe = app.probe.clone().unwrap();
        assert_eq!(probe.frame, Some(4));
        assert!(probe.lines.iter().any(|l| l.starts_with("U ")), "{probe:?}");
    }

    #[test]
    fn the_capture_shows_the_legend_and_panel_over_the_plot() {
        let mut app = solved_app();
        app.run("results.deformed", json!({"step": "Load"})).unwrap();
        let image = app.capture(capture::CaptureSize { width: 1200.0, height: 800.0, scale: 1.0 }).unwrap();
        let t = theme::Tokens::dark();
        // The legend's top swatch (the highest band, red) is drawn left of the middle.
        let red = t.contour[4];
        let near = |c: egui::Color32, d: egui::Color32| {
            (i32::from(c.r()) - i32::from(d.r())).abs() < 30
                && (i32::from(c.g()) - i32::from(d.g())).abs() < 40
                && (i32::from(c.b()) - i32::from(d.b())).abs() < 30
        };
        let found = (0..600).any(|x| (150..400).any(|y| near(image.pixels[y * 1200 + x], red)));
        assert!(found, "a red legend swatch");
    }

    fn region(app: &ModelApp, name: &str) -> egui::Rect {
        app.layout.get(name).copied().unwrap_or_else(|| panic!("no {name} in {:?}", app.layout.keys()))
    }

    /// The family layout: title bar, ribbon, browser left, viewport, command line and timeline at
    /// the bottom; the context panel on the right only while there is context.
    #[test]
    fn the_window_is_laid_out_like_the_other_crafts() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        let size = capture::CaptureSize { width: 1400.0, height: 900.0, scale: 1.0 };
        // `layout` is the live window's; a capture leaves it, so draw a frame of the UI directly.
        let frame = |app: &mut ModelApp| {
            let ctx = egui::Context::default();
            for _ in 0..2 {
                let raw = egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(size.width, size.height))),
                    ..Default::default()
                };
                let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
                out.textures_delta.clear();
            }
        };
        frame(&mut app);
        let (bar, ribbon, browser, view) = (region(&app, "app_bar"), region(&app, "ribbon"), region(&app, "browser"), region(&app, "viewport"));
        let (line, timeline) = (region(&app, "command_line"), region(&app, "timeline"));
        assert_eq!((bar.top(), bar.height()), (0.0, appbar::HEIGHT));
        assert!((ribbon.top() - bar.bottom()).abs() < 1.0 && (ribbon.height() - toolbar::HEIGHT).abs() < 1.0, "{ribbon:?}");
        // The browser's content, inside its margin.
        assert!(browser.left() < 8.0 && (browser.width() - browser::WIDTH).abs() <= 14.0, "{browser:?}");
        assert!(view.left() >= browser.right() - 1.0 && view.top() >= ribbon.bottom() - 1.0 && view.bottom() <= line.top() + 1.0);
        assert!((timeline.bottom() - 900.0).abs() < 1.0 && line.bottom() <= timeline.top() + 1.0);
        assert!((view.right() - 1400.0).abs() < 1.0, "no context: the viewport reaches the right edge");
        assert!(!app.layout.contains_key("context"));
        for chrome in ["nav_bar", "axis_view", "workspace_switcher", "document_tab", "app_bar:solve"] {
            let r = region(&app, chrome);
            assert!(r.is_positive(), "{chrome}");
        }
        assert!(view.contains_rect(region(&app, "nav_bar")) && view.contains_rect(region(&app, "axis_view")));
        // A form brings the context panel; the viewport makes room.
        app.start("model.material.create").unwrap();
        frame(&mut app);
        let context = region(&app, "context");
        assert!(context.right() > 1390.0 && context.width() > 280.0, "{context:?}");
        assert!(region(&app, "viewport").right() <= context.left() + 1.0);
        // The search and a drop-down are drawn when their state says so.
        app.dialog = None;
        app.ui.palette = true;
        app.ui.menu = Some("mesh/Tet".into());
        app.ui.select_tab(Tab::Mesh);
        frame(&mut app);
        assert!(app.layout.contains_key("palette") && app.layout.contains_key("toolbar_dropdown"), "{:?}", app.layout.keys());
        assert_eq!(app.ui.workspace, Workspace::Prepare);
    }

    /// No entry is lost: on every tab, at a wide and a narrow window, every catalogue entry is a
    /// button, a collapsed panel or a » panel, and a row of its panel's drop-down once it is open.
    #[test]
    fn every_catalogue_entry_is_reachable_from_its_tab() {
        for width in [1400.0, 960.0] {
            for tab in Tab::ALL {
                let mut app = ModelApp::new(Session::default());
                app.ui.select_tab(tab);
                let ctx = egui::Context::default();
                let run = |app: &mut ModelApp| {
                    let raw = egui::RawInput {
                        screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(width, 800.0))),
                        ..Default::default()
                    };
                    let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
                    out.textures_delta.clear();
                };
                run(&mut app);
                for c in catalogue::tab_entries(tab) {
                    assert!(app.toolbar_rects.iter().any(|(id, ..)| *id == c.id), "{width}: {} is not on the {tab:?} tab", c.id);
                }
                for p in catalogue::panels(tab) {
                    let Some(panel) = p.first().map(|c| c.panel) else { continue };
                    let collapsed_or_more = app.toolbar_rects.iter().any(|(id, _, f)| *id == p[0].id && matches!(f, toolbar::Fit::More));
                    if collapsed_or_more {
                        continue;
                    }
                    app.ui.menu = Some(toolbar::panel_key(tab, panel));
                    run(&mut app);
                    run(&mut app);
                    for c in p {
                        let row = app.toolbar_rects.iter().any(|(id, r, f)| *id == c.id && *f == toolbar::Fit::Menu && r.is_positive());
                        assert!(row, "{width}: {} has no row in the {panel} drop-down", c.id);
                    }
                    app.ui.menu = None;
                }
            }
        }
    }

    /// A face picked by name stays selected through a remesh and highlights the new mesh's faces;
    /// the set made from it resolves on the new mesh.
    #[test]
    fn a_geometry_pick_survives_a_remesh() {
        let block = |n: u32| {
            let corners = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]];
            json!({"body": "Bar", "kind": "hex8", "corners": corners, "divisions": [n, n, n], "faces": {"4": "tip"}})
        };
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.block", block(2)).unwrap();
        app.viewport.rect = Some(egui::Rect::from_min_size(egui::pos2(100.0, 50.0), egui::vec2(800.0, 600.0)));
        app.set_view(StandardView::Right);
        app.fit_view();
        let at = viewport::to_screen(&app, results::tests::mm([10.0, 5.5, 4.5])).unwrap();
        let tip = selection::Item::Geometry(fieldcraft_model::GeometryRef::face("Bar", "tip").unwrap());
        assert_eq!(app.select_at(at, selection::Mode::Replace), Some(tip.clone()));
        assert_eq!(app.overlay().1.tris.len(), 2 * 2 * 2);
        // A click beside the model clears the selection; a Shift click adds back.
        assert_eq!(app.select_at(egui::pos2(105.0, 55.0), selection::Mode::Replace), None);
        assert!(app.pick.selection.is_empty());
        app.select_at(at, selection::Mode::Add);
        // Remeshed finer: the tip is still selected and highlights the new faces.
        app.run("edit.undo", json!({})).unwrap();
        app.run("mesh.block", block(4)).unwrap();
        app.logic(&egui::Context::default());
        assert_eq!(app.pick.selection.items().iter().collect::<Vec<_>>(), vec![&tip]);
        assert_eq!(app.overlay().1.tris.len(), 4 * 4 * 2);
        app.create_set_from_selection().unwrap();
        assert_eq!(app.dialog.as_ref().map(|d| d.command), Some("model.set.create"));
        app.run_dialog().unwrap();
        assert_eq!(app.session.model.resolve_set("Set-1").unwrap().faces.len(), 16);
        // Mesh labels that are gone leave the selection.
        app.select([selection::Item::Node(1), selection::Item::Element(9999)], selection::Mode::Replace);
        app.run("mesh.node.add", json!({"label": 100_000, "xyz": [50, 50, 50]})).unwrap();
        app.logic(&egui::Context::default());
        assert_eq!(app.pick.selection.items().iter().collect::<Vec<_>>(), vec![&selection::Item::Node(1)]);
    }

    /// A job builds the display data of the session it brings on its worker (the model drawn,
    /// a new design tessellated), and the window adopts it as it is: the same as building it in a
    /// frame, which a large part would stall (#552). A change before it is drawn drops it.
    #[test]
    fn a_job_brings_its_display_data() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("ibeam");
        app.wait_loaded();
        assert!(app.prebuilt.is_some(), "built with the session");
        let adopted = app.display().clone();
        assert!(app.prebuilt.is_none() && !app.bodies.1.is_empty());
        let mut framed = ModelApp::new(app.session.clone());
        assert_eq!(&adopted, framed.display(), "the same as built in a frame");
        assert_eq!(app.bodies.1, framed.bodies.1);
        // A job on the same design (a solve, a mesher) does not tessellate it again.
        assert!(Prebuilt::of(&app.session, &app.bodies.0).bodies.is_none());
        assert!(Prebuilt::of(&app.session, &None).bodies.is_some_and(|(b, _)| b == app.bodies.1));
        app.open_sample("ibeam");
        app.wait_loaded();
        app.run("mesh.node.add", json!({"label": 900000, "xyz": [0, 0, 0]})).unwrap();
        assert!(app.prebuilt.is_none(), "the model changed after the job");
        assert!(app.display().nodes.len() > adopted.nodes.len());
    }

    #[test]
    fn display_tolerance_follows_the_size() {
        let b = |max: [f64; 3]| Aabb3 { min: Vec3::ZERO, max: max.into() };
        assert_eq!(display_tolerance(&b([100.0, 0.0, 0.0])), 0.1);
        assert_eq!(display_tolerance(&b([0.0; 3])), 1e-3);
        assert_eq!(display_tolerance(&b([1e9, 0.0, 0.0])), 1.0);
        assert_eq!(display_tolerance(&Aabb3::EMPTY), 1e-3);
    }

    #[test]
    fn graphics_device_errors_reach_the_status_bar() {
        let mut app = ModelApp::new(Session::default());
        if let Ok(mut e) = app.viewport.gpu_errors.lock() {
            *e = Some("Out of Memory".into());
        }
        app.show_gpu_errors();
        assert_eq!(app.status, Some(Status { text: "Graphics: Out of Memory".into(), error: true }));
        assert!(app.viewport.gpu_errors.lock().unwrap().is_none(), "each error is shown once");
    }

    /// An orbit pivots on the model point under the press (the pick ray's first hit, on the
    /// face toward the eye) and keeps it still on screen; from empty space it pivots on the
    /// centre of what is drawn (#262).
    #[test]
    fn orbit_pivots_on_the_point_under_the_cursor_or_the_drawn_centre() {
        let corners = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0], [0, 0, 10], [10, 0, 10], [10, 10, 10], [0, 10, 10]];
        let mut app = ModelApp::new(Session::default());
        app.run("mesh.block", json!({"body": "Bar", "kind": "hex8", "corners": corners, "divisions": [2, 2, 2]})).unwrap();
        app.viewport.rect = Some(egui::Rect::from_min_size(egui::pos2(100.0, 50.0), egui::vec2(800.0, 600.0)));
        app.set_view(StandardView::Right);
        app.fit_view();
        // The right view looks along −X: the face under the cursor is x = 10.
        let on = results::tests::mm([10.0, 7.0, 3.0]);
        let at = viewport::to_screen(&app, on).unwrap();
        let pivot = app.surface_point(at).unwrap();
        assert!(pivot.dist(on) < 1e-9, "{pivot:?}");
        app.cam.orbit_about(pivot, 60.0, -30.0);
        let after = viewport::to_screen(&app, pivot).unwrap();
        assert!((after - at).length() < 1e-3, "{at:?} -> {after:?}");
        // Empty space: the centre of the drawn bounds.
        let pivot = app.drawn_centre();
        assert!(pivot.dist(results::tests::mm([5.0, 5.0, 5.0])) < 1e-9, "{pivot:?}");
        assert!(app.surface_point(egui::pos2(105.0, 55.0)).is_none());
        // Nothing drawn: the target.
        let mut empty = ModelApp::new(Session::default());
        empty.cam.target = Vec3::new(1.0, 2.0, 3.0);
        assert_eq!(empty.drawn_centre(), Vec3::new(1.0, 2.0, 3.0));
    }
}
