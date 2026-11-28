//! Shared Craft-family panel host (`solidcraft_ui_kit::panel_runtime`).
//!
//! Owns the 4-panel registry (`browser`, `history`, `results`, `inspector`), `PanelRuntime`,
//! `WorkspacePanelHost`, `FieldCraftLayoutStore` persistence, read-through `UiState` sync,
//! `HistoryDock::Top` strip fallback, and `browser_float` translucency scoping.

pub mod content;
#[cfg(not(target_arch = "wasm32"))]
pub mod native;
pub mod store;

use std::cell::RefCell;
use std::collections::BTreeSet;
use std::rc::Rc;

use egui::{Rect, pos2, vec2};
use serde_json::{Value, json};
use solidcraft_panel_layout::{
    Column, ColumnId, Direction, DocumentId, EdgeRule, HeaderActionSpec, HostCapabilities, HostId, InstanceId, Layout, LayoutAction, LayoutService,
    LogicalRect, MAX_LOGICAL, MAX_SECTIONS, Mutation, PanelFeatures, PanelInstance, PanelKind, PanelRegistry, PanelResult, PanelSpec, PanelSurface,
    PlacementMemory, Point, Rule, Section, SectionId, SectionLocation, SideHeightPolicy, Size, ViewportId, WorkspaceId, Zone, column_bounds, error,
    normalize_column_widths, tab_bounds, validate,
};
use solidcraft_ui_kit as kit;
use solidcraft_ui_kit::panel_runtime::{
    ContentServices, ControllerBinding, HostEvent, HostOptions, HostPointerTarget, HostResizeEdge, InAppNativeAdapter, Motion, NativeHostAdapter,
    Ownership, PanelRuntime, RuntimeEnvironment, WorkspacePanelHost,
};
use solidcraft_ui_text::tr;

use self::content::{
    FieldCraftCommands, FieldCraftFrameContent, FieldCraftLocal, FieldCraftModels, FieldCraftPanelController, PanelCommand, PanelSharedState,
};
use self::store::FieldCraftLayoutStore;
use crate::ModelApp;
use crate::files::Prefs;
use crate::shell::HistoryDock;

pub const WORKSPACE_ID: &str = "fieldcraft.main";
pub const PANEL_BROWSER: &str = "browser";
pub const PANEL_HISTORY: &str = "history";
pub const PANEL_RESULTS: &str = "results";
pub const PANEL_INSPECTOR: &str = "inspector";
pub const ALL_PANEL_IDS: [&str; 4] = [PANEL_BROWSER, PANEL_HISTORY, PANEL_RESULTS, PANEL_INSPECTOR];

pub const KIND_BROWSER: &str = "fieldcraft.browser";
pub const KIND_HISTORY: &str = "fieldcraft.history";
pub const KIND_RESULTS: &str = "fieldcraft.results";
pub const KIND_INSPECTOR: &str = "fieldcraft.inspector";

pub const MIN_PANEL_SIZE: f32 = 180.0;
pub const DEFAULT_BROWSER_HEIGHT: f32 = 520.0;
pub const DEFAULT_HISTORY_WIDTH: f32 = 640.0;
pub const DEFAULT_RESULTS_HEIGHT: f32 = 420.0;
pub const DEFAULT_INSPECTOR_HEIGHT: f32 = 360.0;
pub const RIGHT_INSPECTOR_WEIGHT: f32 = 0.45;
pub const RIGHT_RESULTS_WEIGHT: f32 = 0.55;

pub(crate) fn panel_kind(id: &str) -> &'static str {
    match id {
        PANEL_BROWSER => KIND_BROWSER,
        PANEL_HISTORY => KIND_HISTORY,
        PANEL_RESULTS => KIND_RESULTS,
        PANEL_INSPECTOR => KIND_INSPECTOR,
        _ => "",
    }
}

fn logical_rect(r: Rect) -> LogicalRect {
    let valid = r.is_finite() && r.width() >= 0.0 && r.height() >= 0.0;
    if valid {
        LogicalRect { origin: Point { x: r.left(), y: r.top() }, size: Size { width: r.width().max(1.0), height: r.height().max(1.0) } }
    } else {
        LogicalRect { origin: Point { x: 0.0, y: 0.0 }, size: Size { width: 1.0, height: 1.0 } }
    }
}

fn to_egui_rect(r: LogicalRect) -> Rect {
    Rect::from_min_size(pos2(r.origin.x, r.origin.y), vec2(r.size.width, r.size.height))
}

/// Register FieldCraft's four persistent panels.
pub(crate) fn build_registry() -> PanelResult<PanelRegistry> {
    let mut registry = PanelRegistry::default();
    let side_zones: BTreeSet<Zone> = [Zone::Left, Zone::Right, Zone::Floating, Zone::Window].into_iter().collect();
    let history_zones: BTreeSet<Zone> = [Zone::Bottom, Zone::Floating, Zone::Window].into_iter().collect();

    let browser_kind = PanelKind::new(KIND_BROWSER)?;
    registry.register(PanelSpec {
        kind: browser_kind.clone(),
        title_key: "panel.browser".into(),
        icon: "list".into(),
        allowed_zones: side_zones.clone(),
        default_zone: Zone::Left,
        min_size: Size { width: crate::browser::MIN_WIDTH, height: MIN_PANEL_SIZE },
        default_size: Size { width: crate::browser::WIDTH, height: DEFAULT_BROWSER_HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Translucent,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: false, popout: true, tabbable: true, ..Default::default() },
        side_height_policy: SideHeightPolicy::Fill,
        fixed_bottom_height: None,
        header_actions: Vec::new(),
        text_heavy: false,
    })?;
    registry.set_side_widths(&browser_kind, vec![(crate::browser::MIN_WIDTH, crate::browser::MAX_WIDTH)])?;

    let history_actions = [
        ("first", "ui.history_first", "timeline_first"),
        ("back", "ui.history_back", "timeline_back"),
        ("forward", "ui.history_forward", "timeline_forward"),
        ("last", "ui.history_last", "timeline_last"),
    ]
    .into_iter()
    .map(|(id, label_key, icon)| {
        Ok(HeaderActionSpec { id: solidcraft_panel_layout::ControlId::new(id)?, label_key: label_key.into(), icon: icon.into() })
    })
    .collect::<PanelResult<Vec<_>>>()?;

    registry.register(PanelSpec {
        kind: PanelKind::new(KIND_HISTORY)?,
        title_key: "panel.history".into(),
        icon: "ui.views".into(),
        allowed_zones: history_zones,
        default_zone: Zone::Bottom,
        min_size: Size { width: 240.0, height: crate::timeline::HEIGHT },
        default_size: Size { width: DEFAULT_HISTORY_WIDTH, height: crate::timeline::HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Opaque,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: false, popout: true, tabbable: false, ..Default::default() },
        side_height_policy: SideHeightPolicy::Content,
        fixed_bottom_height: Some(crate::timeline::HEIGHT),
        header_actions: history_actions,
        text_heavy: false,
    })?;

    registry.register(PanelSpec {
        kind: PanelKind::new(KIND_INSPECTOR)?,
        title_key: "panel.inspector".into(),
        icon: "list".into(),
        allowed_zones: side_zones.clone(),
        default_zone: Zone::Right,
        min_size: Size { width: MIN_PANEL_SIZE, height: MIN_PANEL_SIZE },
        default_size: Size { width: crate::context::WIDTH, height: DEFAULT_INSPECTOR_HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Opaque,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: true, popout: true, tabbable: true, ..Default::default() },
        side_height_policy: SideHeightPolicy::Fill,
        fixed_bottom_height: None,
        header_actions: Vec::new(),
        text_heavy: false,
    })?;

    registry.register(PanelSpec {
        kind: PanelKind::new(KIND_RESULTS)?,
        title_key: "panel.results".into(),
        icon: "results.deformed".into(),
        allowed_zones: side_zones,
        default_zone: Zone::Right,
        min_size: Size { width: MIN_PANEL_SIZE, height: MIN_PANEL_SIZE },
        default_size: Size { width: crate::context::WIDTH, height: DEFAULT_RESULTS_HEIGHT },
        max_size: Size { width: MAX_LOGICAL, height: MAX_LOGICAL },
        edge_rule: EdgeRule::None,
        surface: PanelSurface::Opaque,
        features: PanelFeatures { collapsible: true, closable: true, hot_closable: true, popout: true, tabbable: true, ..Default::default() },
        side_height_policy: SideHeightPolicy::Fill,
        fixed_bottom_height: None,
        header_actions: Vec::new(),
        text_heavy: false,
    })?;

    Ok(registry)
}

/// Seed the initial layout from legacy `Prefs` (`browser_width`, `browser_collapsed`).
/// `history_dock == HistoryDock::Top` is tracked alongside `Zone::Bottom` via `PanelHost`.
pub(crate) fn seeded_default_layout(registry: &PanelRegistry, prefs: &Prefs) -> PanelResult<Layout> {
    let workspace = WorkspaceId::new(WORKSPACE_ID)?;
    let mut layout = Layout::empty(workspace);

    let left_col = ColumnId::new("col-left")?;
    let right_col = ColumnId::new("col-right")?;
    let browser_w = crate::browser::width(prefs.browser_width);
    let right_w = crate::context::WIDTH;

    layout.columns.push(Column { id: left_col.clone(), zone: Zone::Left, width: browser_w });
    layout.columns.push(Column { id: right_col.clone(), zone: Zone::Right, width: right_w });

    struct Entry {
        id: &'static str,
        section: &'static str,
        zone: Zone,
        column: Option<ColumnId>,
        size: Size,
        collapsed: bool,
        weight: f32,
    }

    let entries = [
        Entry {
            id: PANEL_BROWSER,
            section: "sec-browser",
            zone: Zone::Left,
            column: Some(left_col),
            size: Size { width: browser_w, height: DEFAULT_BROWSER_HEIGHT },
            collapsed: prefs.browser_collapsed,
            weight: 1.0,
        },
        Entry {
            id: PANEL_INSPECTOR,
            section: "sec-inspector",
            zone: Zone::Right,
            column: Some(right_col.clone()),
            size: Size { width: right_w, height: DEFAULT_INSPECTOR_HEIGHT },
            collapsed: false,
            weight: RIGHT_INSPECTOR_WEIGHT,
        },
        Entry {
            id: PANEL_RESULTS,
            section: "sec-results",
            zone: Zone::Right,
            column: Some(right_col),
            size: Size { width: right_w, height: DEFAULT_RESULTS_HEIGHT },
            collapsed: false,
            weight: RIGHT_RESULTS_WEIGHT,
        },
        Entry {
            id: PANEL_HISTORY,
            section: "sec-history",
            zone: Zone::Bottom,
            column: None,
            size: Size { width: DEFAULT_HISTORY_WIDTH, height: crate::timeline::HEIGHT },
            collapsed: false,
            weight: 1.0,
        },
    ];

    for (section_index, entry) in entries.into_iter().enumerate() {
        let id = InstanceId::new(entry.id)?;
        let kind = PanelKind::new(panel_kind(entry.id))?;
        let section = SectionId::new(entry.section)?;
        let location = SectionLocation::Dock { zone: entry.zone };
        let column_width = entry.column.as_ref().map(|_| entry.size.width);

        layout.sections.push(Section {
            id: section.clone(),
            location: location.clone(),
            tabs: vec![id.clone()],
            active: id.clone(),
            size: entry.size,
            collapsed: entry.collapsed,
            column: entry.column.clone(),
            weight: entry.weight,
            preferred_extent: None,
            extent_priority: 0,
        });
        layout.instances.push(PanelInstance {
            id,
            kind,
            requested_open: true,
            memory: PlacementMemory {
                section,
                location,
                section_index,
                tab_index: 0,
                size: entry.size,
                collapsed: entry.collapsed,
                column: entry.column,
                column_index: 0,
                column_width,
                weight: entry.weight,
                preferred_extent: None,
                extent_priority: 0,
            },
        });
    }

    normalize_column_widths(&mut layout, registry)?;
    validate(&layout, registry)?;
    Ok(layout)
}

/// Saved History placement captured before Settings edits `history_dock`, so Cancel restores a
/// pre-edit floating position as well as `HistoryDock`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct HistoryPlacement {
    pub(crate) dock: HistoryDock,
    pub(crate) section: Option<(SectionId, SectionLocation, Size, bool)>,
    pub(crate) memory: PlacementMemory,
}

#[cfg(not(target_arch = "wasm32"))]
type ActiveNativeAdapter = native::NativeAdapter;
#[cfg(target_arch = "wasm32")]
type ActiveNativeAdapter = InAppNativeAdapter;

/// Host wrapper owning `PanelRuntime`, `WorkspacePanelHost`, `FieldCraftLayoutStore`, and
/// read/write adapters for FieldCraft's window shell.
pub(crate) struct PanelHost {
    runtime: PanelRuntime,
    host: WorkspacePanelHost,
    native: ActiveNativeAdapter,
    #[cfg(not(target_arch = "wasm32"))]
    coordinator: native::NativeViewportCoordinator,
    store: FieldCraftLayoutStore,
    shared: Rc<RefCell<PanelSharedState>>,
    local: FieldCraftLocal,
    environment: RuntimeEnvironment,
    native_windows: bool,
    history_dock: HistoryDock,
    pending_browser_width: Option<f32>,
    pending_browser_collapsed: Option<bool>,
    surface_alpha: kit::blur::State,
    omitted_native: BTreeSet<String>,
}

impl PanelHost {
    pub(crate) fn from_prefs(prefs: &Prefs, native_windows: bool) -> PanelResult<Self> {
        let registry = build_registry()?;
        let (layout, store) = FieldCraftLayoutStore::load_from_prefs(prefs, &registry, native_windows)?;
        Self::from_loaded(layout, store, registry, prefs.history_dock, native_windows)
    }

    /// Build a fresh read-only in-window host seeded from the live layout snapshot for `ui.capture`.
    /// Popped-out native OS windows (`Zone::Window`) are omitted from the offscreen capture.
    pub(crate) fn headless_snapshot(&self, prefs: &Prefs) -> PanelResult<Self> {
        let registry = self.runtime.layout().registry().clone();
        let mut layout = self.runtime.layout().layout().clone();
        let mut omitted = BTreeSet::new();
        for section in &layout.sections {
            if section.location.zone() == Zone::Window {
                for id in &section.tabs {
                    omitted.insert(id.0.clone());
                }
            }
        }
        layout.sections.retain(|s| s.location.zone() != Zone::Window);
        for inst in &mut layout.instances {
            if omitted.contains(&inst.id.0) {
                inst.requested_open = false;
                if matches!(inst.memory.location, SectionLocation::Window { .. }) {
                    inst.memory.location = SectionLocation::Floating { origin: Point { x: 0.0, y: 0.0 } };
                    inst.memory.column = None;
                    inst.memory.column_index = 0;
                    inst.memory.column_width = None;
                }
            }
        }
        validate(&layout, &registry)?;
        let mut snap = prefs.clone();
        snap.history_dock = self.history_dock;
        snap.panel_layout = serde_json::to_string(&solidcraft_panel_layout::LayoutEnvelope::new(layout))
            .ok()
            .and_then(|s| crate::files::RawJson::from_json_string(s).ok());
        let (sanitized, store) = FieldCraftLayoutStore::load_from_prefs(&snap, &registry, false)?;
        let mut host = Self::from_loaded(sanitized, store, registry, self.history_dock, false)?;
        host.omitted_native = omitted;
        Ok(host)
    }

    fn build_runtime(
        layout: Layout,
        registry: PanelRegistry,
        environment: RuntimeEnvironment,
        shared: &Rc<RefCell<PanelSharedState>>,
        native: &mut ActiveNativeAdapter,
        local: &mut FieldCraftLocal,
        native_windows: bool,
    ) -> PanelResult<PanelRuntime> {
        let capabilities = HostCapabilities { persistence: true, native_windows, max_native_hosts: if native_windows { MAX_SECTIONS } else { 0 } };
        let has_native = layout.sections.iter().any(|s| s.location.zone() == Zone::Window);
        let initial_layout = if native_windows && has_native { seeded_default_layout(&registry, &Prefs::default())? } else { layout.clone() };
        let service = LayoutService::new(registry, initial_layout, capabilities)?;
        let mut runtime = PanelRuntime::new(service, environment)?;
        for id in ALL_PANEL_IDS {
            let shared_ref = shared.clone();
            runtime.register(
                PanelKind::new(panel_kind(id))?,
                ControllerBinding {
                    ownership: Ownership::Workspace,
                    create: Box::new(move |_inst, _owner| Ok(Box::new(FieldCraftPanelController::new(id, shared_ref.clone())))),
                },
            )?;
        }
        if native_windows && has_native {
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local };
            runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout }, native, &mut services)?;
        }
        Ok(runtime)
    }

    fn from_loaded(
        layout: Layout,
        store: FieldCraftLayoutStore,
        registry: PanelRegistry,
        history_dock: HistoryDock,
        native_windows: bool,
    ) -> PanelResult<Self> {
        let environment = RuntimeEnvironment {
            workspace: WorkspaceId::new(WORKSPACE_ID)?,
            document: Some(DocumentId::new("fieldcraft.active")?),
            window: HostId::new("main")?,
            viewport: ViewportId::new("main")?,
            direction: Direction::Ltr,
            motion: Motion::Animated,
        };
        let shared = Rc::new(RefCell::new(PanelSharedState {
            history_available: true,
            results_available: false,
            inspector_available: false,
            inspector_heading: None,
            commands: Default::default(),
        }));
        let mut local = FieldCraftLocal::default();
        #[cfg(not(target_arch = "wasm32"))]
        let mut native = if native_windows {
            native::NativeAdapter::Desktop(solidcraft_ui_kit::panel_runtime::DesktopNativeAdapter::default())
        } else {
            native::NativeAdapter::InApp(InAppNativeAdapter)
        };
        #[cfg(target_arch = "wasm32")]
        let mut native = InAppNativeAdapter;

        let runtime = Self::build_runtime(layout, registry, environment.clone(), &shared, &mut native, &mut local, native_windows)?;

        let mut options = HostOptions { content_header_height: 0.0, ..Default::default() };
        options.geometry.default_opaque = false;
        options.metrics.footer_height = 0.0;
        options.metrics.gap = kit::SHELL_GAP;
        let mut host = WorkspacePanelHost::new(options)?;
        host.note_native_capabilities(native.capabilities());

        let mut this = Self {
            runtime,
            host,
            native,
            #[cfg(not(target_arch = "wasm32"))]
            coordinator: native::NativeViewportCoordinator::default(),
            store,
            shared,
            local,
            environment,
            native_windows,
            history_dock,
            pending_browser_width: None,
            pending_browser_collapsed: None,
            surface_alpha: kit::blur::State::default(),
            omitted_native: BTreeSet::new(),
        };
        let _ = this.reconcile_runtime();
        Ok(this)
    }

    fn reconcile_runtime(&mut self) -> PanelResult<()> {
        let models = FieldCraftModels;
        let mut commands = FieldCraftCommands;
        let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
        self.runtime.reconcile_context(self.environment.clone(), &mut services)
    }

    pub(crate) fn is_open(&self, id: &str) -> bool {
        InstanceId::new(id).ok().is_some_and(|inst| self.runtime.layout().layout().instance(&inst).is_ok_and(|p| p.requested_open))
    }

    fn history_in_bottom_zone(&self) -> bool {
        InstanceId::new(PANEL_HISTORY)
            .ok()
            .is_some_and(|inst| self.runtime.layout().layout().containing(&inst).is_ok_and(|s| s.location.zone() == Zone::Bottom))
    }

    pub(crate) fn apply_mutations(&mut self, mutations: Vec<Mutation>, mark_dirty: bool) -> PanelResult<()> {
        if mutations.is_empty() {
            return Ok(());
        }
        let models = FieldCraftModels;
        let mut commands = FieldCraftCommands;
        let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
        let action = LayoutAction::Apply { expected_revision: None, mutations };
        self.runtime.dispatch(action, &mut self.native, &mut services)?;
        if mark_dirty {
            self.store.replace_live(self.runtime.layout().layout(), self.runtime.layout().registry(), true)?;
        }
        Ok(())
    }

    fn ensure_open(&mut self, id: &str, open: bool, mark_dirty: bool) -> PanelResult<()> {
        if self.is_open(id) == open {
            return Ok(());
        }
        let inst = InstanceId::new(id)?;
        let kind = PanelKind::new(panel_kind(id))?;
        let mutation = if open { Mutation::Open { kind, instance: inst.clone(), destination: None } } else { Mutation::Close { instance: inst } };
        self.apply_mutations(vec![mutation], mark_dirty)?;
        if open && id == PANEL_BROWSER {
            let inst = InstanceId::new(PANEL_BROWSER)?;
            if let Some(w) = self.pending_browser_width.take() {
                const WIDTH_EPSILON_PT: f32 = 1e-3; // unit-factors: logical-point UI width tolerance
                let needs_width = self.runtime.layout().layout().containing(&inst).ok().is_some_and(|sec| match &sec.column {
                    Some(col_id) => self
                        .runtime
                        .layout()
                        .layout()
                        .columns
                        .iter()
                        .find(|c| c.id == *col_id)
                        .is_some_and(|c| (c.width - w).abs() > WIDTH_EPSILON_PT),
                    None => (sec.size.width - w).abs() > WIDTH_EPSILON_PT,
                });
                if needs_width {
                    self.set_browser_width(w)?;
                }
            }
            if let Some(c) = self.pending_browser_collapsed.take()
                && self.runtime.layout().layout().containing(&inst).is_ok_and(|sec| sec.collapsed != c)
            {
                self.set_browser_collapsed(c)?;
            }
        }
        Ok(())
    }

    /// Synchronize availability, `HistoryDock::Top`, and `browser_float` translucency, prepare
    /// `WorkspacePanelHost` geometry with `SHELL_GAP` outer insets, and populate `app.shell`
    /// before `context::panel` and `viewport::show` lay out against it.
    pub(crate) fn prepare(&mut self, app: &mut ModelApp, ui: &mut egui::Ui, render_domain: Rect) -> bool {
        if !render_domain.is_finite() || render_domain.width() <= 0.0 || render_domain.height() <= 0.0 {
            return false;
        }
        self.history_dock = app.ui.history_dock;
        let show_tree = app.ui.show_tree;
        let show_timeline = app.ui.show_timeline;

        let results_active = crate::results::shown(&app.session).is_some();
        let inspector_active = crate::context::inspector_active(app);
        let inspector_heading = crate::context::inspector_heading(app);

        // Publish current content availability before calling `ensure_open`, because
        // `PanelRuntime` rejects `Open` for a panel whose current availability is `false`.
        {
            let mut shared = self.shared.borrow_mut();
            shared.history_available = true;
            shared.results_available = results_active;
            shared.inspector_available = inspector_active;
            shared.inspector_heading = inspector_heading.clone();
        }
        if let Err(err) = self.reconcile_runtime() {
            log::warn!("panel runtime reconcile failed: {err}");
            return false;
        }

        // Sync visibility changes initiated outside `panel_runtime` (View menu, `ui.set`).
        if !self.omitted_native.contains(PANEL_BROWSER) {
            let _ = self.ensure_open(PANEL_BROWSER, show_tree, true);
        }
        if !self.omitted_native.contains(PANEL_HISTORY) {
            let _ = self.ensure_open(PANEL_HISTORY, show_timeline, true);
        }
        if !self.omitted_native.contains(PANEL_RESULTS) && results_active && !self.is_open(PANEL_RESULTS) {
            let _ = self.ensure_open(PANEL_RESULTS, true, false);
        }
        if !self.omitted_native.contains(PANEL_INSPECTOR) && inspector_active && !self.is_open(PANEL_INSPECTOR) {
            let _ = self.ensure_open(PANEL_INSPECTOR, true, false);
        }

        let top_strip =
            !self.omitted_native.contains(PANEL_HISTORY) && show_timeline && self.history_dock == HistoryDock::Top && self.history_in_bottom_zone();
        {
            let mut shared = self.shared.borrow_mut();
            shared.history_available = !top_strip && !self.omitted_native.contains(PANEL_HISTORY);
            shared.results_available = results_active && !self.omitted_native.contains(PANEL_RESULTS);
            shared.inspector_available = inspector_active && !self.omitted_native.contains(PANEL_INSPECTOR);
            shared.inspector_heading = inspector_heading;
        }

        let top_rect = if top_strip {
            let cards = app.shell.cards.clone();
            let rect = crate::timeline::top_strip_rect(render_domain, &cards);
            if let Some(r) = rect {
                app.shell.record_top_dock(r);
            }
            rect
        } else {
            None
        };

        if let Err(err) = self.reconcile_runtime() {
            log::warn!("panel runtime reconcile failed: {err}");
            return false;
        }

        let gap = kit::Container::new(ui.ctx(), kit::Surface::Panel).gap();
        self.host.note_native_capabilities(self.native.capabilities());
        let _ = self.host.set_requested_opacity(None);
        let _ = self.host.set_gap(gap);
        let float_alpha = app.ui.browser_float.then_some((crate::browser::FLOAT_ALPHA * 255.0).round() as u8);
        self.surface_alpha.set_surface_alpha(float_alpha);
        let surface_scope = self.surface_alpha.begin(ui.ctx());
        self.host.note_egui_context(ui.ctx(), !app.ui.browser_float);
        drop(surface_scope);

        if let Err(err) = self.host.capture_metrics(&self.runtime, &FieldCraftModels) {
            log::warn!("panel capture_metrics failed: {err}");
            return false;
        }
        if let Err(err) = self.host.capture_reservations(&self.runtime, &self.store) {
            log::warn!("panel capture_reservations failed: {err}");
            return false;
        }

        let mut placement = app.shell.panel_space(render_domain);
        placement.min.x += gap.min(placement.width() / 2.0);
        placement.max.x -= gap.min(placement.width());
        placement.min.y += gap.min(placement.height() / 2.0);
        placement.max.y -= gap.min(placement.height());

        let chrome: Vec<LogicalRect> = app
            .shell
            .cards
            .iter()
            .copied()
            .chain(top_rect)
            .map(|r| r.intersect(render_domain))
            .filter(|r| r.is_positive() && r.is_finite())
            .map(logical_rect)
            .collect();

        if let Err(err) = self.host.prepare_in_domain(&self.runtime, logical_rect(render_domain), logical_rect(placement), &chrome) {
            log::warn!("panel prepare_in_domain failed: {err}");
            return false;
        }

        if let Some(frame) = self.host.frame() {
            for section in &frame.layout.sections {
                app.shell.set_panel_geometry(section.zone, to_egui_rect(section.outer));
            }
        }
        true
    }

    /// Paint prepared `WorkspacePanelHost` surfaces above `viewport::show` and publish layout
    /// regions and registry handles for automation/tests.
    pub(crate) fn paint_and_finish(&mut self, app: &mut ModelApp, ui: &mut egui::Ui, render_domain: Rect) {
        let top_strip = app.ui.show_timeline && self.history_dock == HistoryDock::Top && self.history_in_bottom_zone();
        if top_strip {
            let cards = app.shell.cards.clone();
            crate::timeline::top_strip_contents(app, ui, render_domain, &cards);
        }

        {
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
            let mut content = FieldCraftFrameContent { app };
            if let Err(err) = self.host.show_deferred_with_content(ui, &mut self.runtime, &mut self.native, &mut services, &mut content) {
                log::warn!("panel show_deferred_with_content failed: {err}");
                return;
            }
        }

        {
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
            if let Err(err) = self.host.finish(ui, &mut self.runtime, &mut self.native, &mut services, &mut self.store) {
                log::warn!("panel host finish failed: {err}");
            }
        }

        // Publish zone geometry and legacy layout/test handles from the completed frame.
        if let Some(frame) = self.host.frame() {
            let committed = frame.events.iter().any(|e| matches!(e, HostEvent::Committed { .. } | HostEvent::PreparedCommitted { .. }));
            let active_layout = frame.next_layout.as_ref().unwrap_or(&frame.layout);
            for section in &active_layout.sections {
                let outer = to_egui_rect(section.outer);
                app.shell.set_panel_geometry(section.zone, outer);
                match section.instance.0.as_str() {
                    PANEL_BROWSER => {
                        app.place("browser", outer);
                        app.place("browser:header", to_egui_rect(section.header));
                        if let Some(c) = section.collapse {
                            app.place("browser:collapse", to_egui_rect(c));
                        }
                        if let Some(rail) = active_layout.controls.iter().find(|ctrl| {
                            matches!(
                                &ctrl.target,
                                HostPointerTarget::Resize { section: sid, edge }
                                    if *sid == section.section && matches!(edge, HostResizeEdge::Right | HostResizeEdge::Left)
                            )
                        }) {
                            app.place("browser:rail", to_egui_rect(rail.bounds));
                        }
                    }
                    PANEL_HISTORY => {
                        app.place("timeline", outer);
                        if let Some(title) = section.title {
                            solidcraft_ui_registry::publish_handle("timeline:grip", to_egui_rect(title).center());
                        }
                        for (action_id, rect) in &section.actions {
                            solidcraft_ui_registry::publish_handle(&format!("timeline:{}", action_id.0), to_egui_rect(*rect).center());
                        }
                    }
                    PANEL_RESULTS => {
                        app.place("context:results", outer);
                        let merged = app.layout.get("context").copied().map_or(outer, |prev| prev.union(outer));
                        app.place("context", merged);
                    }
                    PANEL_INSPECTOR => {
                        if let Some(key) = crate::context::inspector_region_key(app) {
                            app.place(key, outer);
                            let merged = app.layout.get("context").copied().map_or(outer, |prev| prev.union(outer));
                            app.place("context", merged);
                        }
                    }
                    _ => {}
                }
            }
            if committed {
                self.sync_read_through_to_app(app);
            }
        }
    }

    pub(crate) fn drain_commands(&mut self) -> Vec<PanelCommand> {
        self.shared.borrow_mut().drain_commands()
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn downgrade_to_in_window(&mut self) {
        let registry = self.runtime.layout().registry().clone();
        if let Ok(projected) = self.store.project_live_to_floating(self.runtime.layout().layout(), &registry) {
            let mut next_native = native::NativeAdapter::InApp(InAppNativeAdapter);
            if let Ok(runtime) =
                Self::build_runtime(projected, registry, self.environment.clone(), &self.shared, &mut next_native, &mut self.local, false)
            {
                self.runtime = runtime;
                self.native = next_native;
                self.native_windows = false;
                self.coordinator.clear_active_state();
                self.host.note_native_capabilities(self.native.capabilities());
                let _ = self.reconcile_runtime();
                log::info!("downgraded panel host to in-window mode");
            }
        }
    }

    /// Update desktop monitor topology and dynamically upgrade/downgrade `native_windows` capability.
    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn sync_desktop_monitors(
        &mut self,
        _prefs: &Prefs,
        monitors: impl FnOnce() -> Vec<solidcraft_ui_kit::panel_runtime::NativeMonitor>,
        global_coords_available: bool,
    ) -> bool {
        let supported = self.coordinator.update_monitors(self.native.desktop(), monitors, global_coords_available);
        if supported && !self.native_windows {
            let registry = self.runtime.layout().registry().clone();
            let layout = if self.coordinator.has_pending_close_persist() {
                self.runtime.layout().layout().clone()
            } else {
                self.store.take_preserved_native().map(|e| e.layout).unwrap_or_else(|| self.runtime.layout().layout().clone())
            };
            let mut next_native = native::NativeAdapter::Desktop(solidcraft_ui_kit::panel_runtime::DesktopNativeAdapter::default());
            match Self::build_runtime(layout, registry, self.environment.clone(), &self.shared, &mut next_native, &mut self.local, true) {
                Ok(runtime) => {
                    if let Some(desktop) = next_native.desktop() {
                        self.coordinator.push_to_desktop(desktop);
                    }
                    self.runtime = runtime;
                    self.native = next_native;
                    self.native_windows = true;
                    self.host.note_native_capabilities(self.native.capabilities());
                    let _ = self.reconcile_runtime();
                    log::info!("upgraded panel host to desktop native windows");
                }
                Err(err) => {
                    log::warn!("failed to upgrade panel host to desktop native windows: {err}");
                }
            }
        } else if !supported && self.native_windows {
            self.downgrade_to_in_window();
        }
        self.native_windows
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn mark_exiting(&mut self) {
        self.coordinator.mark_exiting();
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn show_desktop_windows(&mut self, app: &mut ModelApp, ctx: &egui::Context) {
        let mut driver = native::EguiViewportDriver;
        if let Err(err) = self.show_desktop_windows_with(&mut driver, app, ctx) {
            log::warn!("show desktop panel windows failed: {err}");
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn show_desktop_windows_with<D: native::ViewportDriver>(
        &mut self,
        driver: &mut D,
        app: &mut ModelApp,
        ctx: &egui::Context,
    ) -> PanelResult<()> {
        // Keep shared availability/headings current even when Start page is shown.
        let top_strip = app.ui.show_timeline && self.history_dock == HistoryDock::Top && self.history_in_bottom_zone();
        let results_active = crate::results::shown(&app.session).is_some();
        let inspector_active = crate::context::inspector_active(app);
        let inspector_heading = crate::context::inspector_heading(app);
        {
            let mut shared = self.shared.borrow_mut();
            shared.history_available = !top_strip;
            shared.results_available = results_active;
            shared.inspector_available = inspector_active;
            shared.inspector_heading = inspector_heading;
        }
        let _ = self.reconcile_runtime();
        let mut coordinator = std::mem::take(&mut self.coordinator);
        let res =
            coordinator.show_desktop_windows(driver, app, ctx, &mut self.runtime, &mut self.host, &mut self.native, &mut self.store, &mut self.local);
        self.coordinator = coordinator;
        if self.coordinator.is_embedded_downgraded() && self.native_windows {
            self.downgrade_to_in_window();
        } else {
            let _ = self.reconcile_runtime();
        }
        self.host.note_native_capabilities(self.native.capabilities());
        res
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn finish_cross_viewport_sync(
        &mut self,
        app: &ModelApp,
        ctx: &egui::Context,
        had_dialog: bool,
        had_refused: bool,
        had_drained_commands: bool,
    ) {
        let mut driver = native::EguiViewportDriver;
        self.finish_cross_viewport_sync_with(&mut driver, app, ctx, had_dialog, had_refused, had_drained_commands);
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(crate) fn finish_cross_viewport_sync_with<D: native::ViewportDriver>(
        &mut self,
        driver: &mut D,
        app: &ModelApp,
        ctx: &egui::Context,
        had_dialog: bool,
        had_refused: bool,
        had_drained_commands: bool,
    ) {
        self.coordinator.sync_post_step7(driver, app, ctx, had_dialog, had_refused, had_drained_commands);
    }

    fn sync_read_through_to_app(&mut self, app: &mut ModelApp) {
        app.ui.browser_width = self.browser_width();
        app.ui.browser_collapsed = self.browser_collapsed();
        if !self.history_in_bottom_zone() && app.ui.history_dock == HistoryDock::Top {
            self.history_dock = HistoryDock::Bottom;
            app.ui.history_dock = HistoryDock::Bottom;
        }
    }

    pub(crate) fn browser_width(&self) -> f32 {
        let Ok(inst) = InstanceId::new(PANEL_BROWSER) else { return crate::browser::WIDTH };
        let layout = self.runtime.layout().layout();
        if let Ok(sec) = layout.containing(&inst) {
            if let Some(col_id) = &sec.column
                && let Some(col) = layout.columns.iter().find(|c| c.id == *col_id)
            {
                return crate::browser::width(col.width);
            }
            return crate::browser::width(sec.size.width);
        }
        layout.instance(&inst).map_or(crate::browser::WIDTH, |i| crate::browser::width(i.memory.column_width.unwrap_or(i.memory.size.width)))
    }

    pub(crate) fn set_browser_width(&mut self, width: f32) -> PanelResult<()> {
        let clamped = crate::browser::width(width);
        let inst = InstanceId::new(PANEL_BROWSER)?;
        let layout = self.runtime.layout().layout().clone();
        if let Ok(sec) = layout.containing(&inst) {
            self.pending_browser_width = None;
            let mutation = if let Some(col_id) = &sec.column {
                Mutation::SetColumnWidth { column: col_id.clone(), width: clamped }
            } else {
                Mutation::SetSize { section: sec.id.clone(), size: Size { width: clamped, height: sec.size.height } }
            };
            self.apply_mutations(vec![mutation], true)
        } else {
            let registry = self.runtime.layout().registry();
            if let Ok(entry) = layout.instance(&inst)
                && let Some(sec) = layout.sections.iter().find(|s| s.id == entry.memory.section)
            {
                if let Some(col_id) = &sec.column {
                    let (min, max) = column_bounds(&layout, registry, col_id)?;
                    if clamped < min || clamped > max {
                        return Err(error(Rule::Bounds, "browser width outside shared column bounds"));
                    }
                } else {
                    let (min, max) = tab_bounds(&layout, registry, sec)?;
                    if clamped < min.width || clamped > max.width {
                        return Err(error(Rule::Bounds, "browser width outside shared section bounds"));
                    }
                }
            }
            let mut updated = layout;
            if let Some(entry) = updated.instances.iter_mut().find(|i| i.id == inst) {
                entry.memory.size.width = clamped;
                if entry.memory.column.is_some() {
                    entry.memory.column_width = Some(clamped);
                }
                if let Some(col_id) = &entry.memory.column
                    && let Some(col) = updated.columns.iter_mut().find(|c| c.id == *col_id)
                {
                    col.width = clamped;
                }
            }
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
            self.runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout: updated.clone() }, &mut self.native, &mut services)?;
            self.store.replace_live(&updated, self.runtime.layout().registry(), true)?;
            self.pending_browser_width = Some(clamped);
            Ok(())
        }
    }

    pub(crate) fn browser_collapsed(&self) -> bool {
        let Ok(inst) = InstanceId::new(PANEL_BROWSER) else { return false };
        let layout = self.runtime.layout().layout();
        if let Ok(sec) = layout.containing(&inst) {
            return sec.collapsed;
        }
        layout.instance(&inst).is_ok_and(|i| i.memory.collapsed)
    }

    pub(crate) fn set_browser_collapsed(&mut self, collapsed: bool) -> PanelResult<()> {
        let inst = InstanceId::new(PANEL_BROWSER)?;
        let layout = self.runtime.layout().layout().clone();
        if let Ok(sec) = layout.containing(&inst) {
            self.pending_browser_collapsed = None;
            if sec.collapsed == collapsed {
                return Ok(());
            }
            self.apply_mutations(vec![Mutation::Collapse { section: sec.id.clone(), collapsed }], true)
        } else {
            let mut updated = layout;
            if let Some(entry) = updated.instances.iter_mut().find(|i| i.id == inst) {
                entry.memory.collapsed = collapsed;
            }
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
            self.runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout: updated.clone() }, &mut self.native, &mut services)?;
            self.store.replace_live(&updated, self.runtime.layout().registry(), true)?;
            self.pending_browser_collapsed = Some(collapsed);
            Ok(())
        }
    }

    pub(crate) fn set_panel_visible(&mut self, id: &str, open: bool) -> PanelResult<()> {
        self.ensure_open(id, open, true)
    }

    pub(crate) fn history_dock(&self) -> HistoryDock {
        self.history_dock
    }

    pub(crate) fn history_placement_snapshot(&self) -> Option<HistoryPlacement> {
        let inst = InstanceId::new(PANEL_HISTORY).ok()?;
        let layout = self.runtime.layout().layout();
        let entry = layout.instance(&inst).ok()?;
        let section = layout.containing(&inst).ok().map(|s| (s.id.clone(), s.location.clone(), s.size, s.collapsed));
        Some(HistoryPlacement { dock: self.history_dock, section, memory: entry.memory.clone() })
    }

    pub(crate) fn restore_history_placement(&mut self, saved: &HistoryPlacement) -> PanelResult<()> {
        self.history_dock = saved.dock;
        let inst = InstanceId::new(PANEL_HISTORY)?;
        let mut updated = self.runtime.layout().layout().clone();
        if let Some(entry) = updated.instances.iter_mut().find(|i| i.id == inst) {
            entry.memory = saved.memory.clone();
        }
        if let Some((sec_id, loc, size, collapsed)) = &saved.section
            && let Some(idx) = updated.sections.iter().position(|s| s.tabs.contains(&inst))
            && updated.sections[idx].tabs.len() == 1
        {
            let id_conflict = *sec_id != updated.sections[idx].id && updated.sections.iter().any(|s| s.id == *sec_id);
            let cur_sec = &mut updated.sections[idx];
            if !id_conflict {
                cur_sec.id = sec_id.clone();
            } else if let Some(entry) = updated.instances.iter_mut().find(|i| i.id == inst) {
                entry.memory.section = cur_sec.id.clone();
            }
            cur_sec.location = loc.clone();
            cur_sec.size = *size;
            cur_sec.collapsed = *collapsed;
            cur_sec.column = None;
        }
        let models = FieldCraftModels;
        let mut commands = FieldCraftCommands;
        let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
        self.runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout: updated.clone() }, &mut self.native, &mut services)?;
        self.store.replace_live(&updated, self.runtime.layout().registry(), true)?;
        let _ = self.reconcile_runtime();
        Ok(())
    }

    pub(crate) fn set_history_dock(&mut self, dock: HistoryDock) -> PanelResult<()> {
        self.history_dock = dock;
        let inst = InstanceId::new(PANEL_HISTORY)?;
        let layout = self.runtime.layout().layout().clone();
        if let Ok(sec) = layout.containing(&inst) {
            if sec.location.zone() != Zone::Bottom {
                let section_id = SectionId::new("sec-history")?;
                let size = Size { width: DEFAULT_HISTORY_WIDTH, height: crate::timeline::HEIGHT };
                self.apply_mutations(
                    vec![Mutation::Move {
                        instance: inst,
                        destination: solidcraft_panel_layout::Destination::NewSection {
                            section: section_id,
                            index: 0,
                            location: SectionLocation::Dock { zone: Zone::Bottom },
                            size,
                        },
                    }],
                    true,
                )?;
            }
        } else if layout.instance(&inst).is_ok_and(|i| i.memory.location.zone() != Zone::Bottom) {
            let section_id = SectionId::new("sec-history")?;
            let size = Size { width: DEFAULT_HISTORY_WIDTH, height: crate::timeline::HEIGHT };
            let mut updated = layout;
            if let Some(entry) = updated.instances.iter_mut().find(|i| i.id == inst) {
                entry.memory.section = section_id;
                entry.memory.location = SectionLocation::Dock { zone: Zone::Bottom };
                entry.memory.section_index = 0;
                entry.memory.tab_index = 0;
                entry.memory.size = size;
                entry.memory.column = None;
                entry.memory.column_index = 0;
                entry.memory.column_width = None;
            }
            let models = FieldCraftModels;
            let mut commands = FieldCraftCommands;
            let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
            self.runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout: updated.clone() }, &mut self.native, &mut services)?;
            self.store.replace_live(&updated, self.runtime.layout().registry(), true)?;
        }
        Ok(())
    }

    pub(crate) fn reset_layout(&mut self) -> PanelResult<()> {
        let defaults = Prefs::default();
        let registry = self.runtime.layout().registry().clone();
        let fresh = seeded_default_layout(&registry, &defaults)?;
        let models = FieldCraftModels;
        let mut commands = FieldCraftCommands;
        let mut services = ContentServices { models: &models, commands: &mut commands, local: &mut self.local };
        self.runtime.dispatch(LayoutAction::Replace { expected_revision: None, layout: fresh.clone() }, &mut self.native, &mut services)?;
        self.history_dock = HistoryDock::Bottom;
        self.pending_browser_width = None;
        self.pending_browser_collapsed = None;
        self.store.replace_live(&fresh, &registry, true)?;
        let _ = self.reconcile_runtime();
        Ok(())
    }

    pub(crate) fn is_dirty(&self) -> bool {
        self.store.is_dirty()
    }

    pub(crate) fn ack_dirty(&mut self) {
        self.store.clear_dirty();
    }

    pub(crate) fn encode_prefs_layout(&self) -> store::EncodedLayout {
        self.store.encode_prefs_layout(self.runtime.layout().layout(), self.runtime.layout().registry())
    }

    #[cfg(test)]
    pub(crate) fn to_prefs_value(&self) -> Option<crate::files::RawJson> {
        self.store.to_prefs_value(self.runtime.layout().layout(), self.runtime.layout().registry())
    }

    /// Deterministic `ui.inspect.panels` array sorted by `id` (`browser`, `history`, `inspector`, `results`).
    pub(crate) fn inspect_panels(&self, app: &ModelApp) -> Vec<Value> {
        let layout = self.runtime.layout().layout();
        let frame = self.host.frame();
        let frame_sections = frame.map(|f| f.next_layout.as_ref().unwrap_or(&f.layout).sections.as_slice()).unwrap_or(&[]);
        let start_shown = app.start_shown();
        let mut ids = ALL_PANEL_IDS;
        ids.sort_unstable();
        ids.into_iter()
            .map(|id| {
                let inst_id = InstanceId::new(id).ok();
                let requested_open = inst_id.as_ref().is_some_and(|i| layout.instance(i).is_ok_and(|p| p.requested_open));
                let has_content = match id {
                    PANEL_BROWSER => true,
                    PANEL_HISTORY => true,
                    PANEL_RESULTS => crate::results::shown(&app.session).is_some(),
                    PANEL_INSPECTOR => crate::context::inspector_active(app),
                    _ => false,
                };
                let section = inst_id.as_ref().and_then(|i| layout.containing(i).ok());
                let zone =
                    if id == PANEL_HISTORY && self.history_dock == HistoryDock::Top && section.is_some_and(|s| s.location.zone() == Zone::Bottom) {
                        "top"
                    } else if let Some(s) = section {
                        match s.location.zone() {
                            Zone::Left => "left",
                            Zone::Right => "right",
                            Zone::Bottom => "bottom",
                            Zone::Floating => "floating",
                            Zone::Window => "native",
                        }
                    } else {
                        "hidden"
                    };
                let is_top_history = id == PANEL_HISTORY && zone == "top";
                let is_native = section.is_some_and(|s| s.location.zone() == Zone::Window);
                let presented = frame_sections.iter().find(|s| s.instance.0 == id);
                let visible = if start_shown {
                    false
                } else if is_top_history {
                    requested_open && has_content
                } else if is_native {
                    requested_open && has_content && section.is_some_and(|s| s.active.0 == id)
                } else if frame.is_some() {
                    presented.is_some()
                } else {
                    requested_open && has_content && section.is_some_and(|s| s.active.0 == id)
                };
                let host = if !visible {
                    "hidden"
                } else if is_native {
                    "native"
                } else if section.is_some_and(|s| s.location.zone() == Zone::Floating) {
                    "floating"
                } else {
                    "docked"
                };
                let rect = if !visible || is_native {
                    None
                } else if is_top_history {
                    app.layout.get("timeline").copied().map(|r| json!([r.left(), r.top(), r.width(), r.height()]))
                } else {
                    presented.map(|s| json!([s.outer.origin.x, s.outer.origin.y, s.outer.size.width, s.outer.size.height]))
                };
                let active_title = match id {
                    PANEL_BROWSER => tr("panel.browser").to_string(),
                    PANEL_HISTORY => tr("panel.history").to_string(),
                    PANEL_RESULTS => tr("panel.results").to_string(),
                    PANEL_INSPECTOR => crate::context::inspector_heading(app).map_or_else(|| tr("panel.inspector").to_string(), |(_, t)| t),
                    _ => id.to_string(),
                };
                json!({
                    "id": id,
                    "visible": visible,
                    "has_content": has_content,
                    "host": host,
                    "zone": zone,
                    "rect": rect,
                    "active_title": active_title,
                    "native_supported": self.native_windows,
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fieldcraft_engine::Session;
    use solidcraft_panel_layout::Destination;

    #[test]
    fn default_layout_has_expected_zones_and_weights() {
        let host = PanelHost::from_prefs(&Prefs::default(), false).unwrap();
        let layout = host.runtime.layout().layout();
        let zone_of = |id: &str| layout.containing(&InstanceId::new(id).unwrap()).unwrap().location.zone();
        assert_eq!(zone_of(PANEL_BROWSER), Zone::Left);
        assert_eq!(zone_of(PANEL_HISTORY), Zone::Bottom);
        assert_eq!(zone_of(PANEL_INSPECTOR), Zone::Right);
        assert_eq!(zone_of(PANEL_RESULTS), Zone::Right);
    }

    #[test]
    fn history_rejects_left_and_right_docks() {
        let mut host = PanelHost::from_prefs(&Prefs::default(), false).unwrap();
        let history = InstanceId::new(PANEL_HISTORY).unwrap();
        for forbidden in [Zone::Left, Zone::Right] {
            let err = host.apply_mutations(
                vec![Mutation::Move {
                    instance: history.clone(),
                    destination: Destination::NewSection {
                        section: SectionId::new("sec-bad-hist").unwrap(),
                        index: 0,
                        location: SectionLocation::Dock { zone: forbidden },
                        size: Size { width: 280.0, height: 320.0 },
                    },
                }],
                true,
            );
            assert!(err.is_err(), "history must reject {forbidden:?}");
        }
    }

    #[test]
    fn browser_supports_dock_float_resize_and_inspector_results_tab_and_split() {
        let mut host = PanelHost::from_prefs(&Prefs::default(), false).unwrap();
        let browser = InstanceId::new(PANEL_BROWSER).unwrap();
        let inspector = InstanceId::new(PANEL_INSPECTOR).unwrap();
        let results = InstanceId::new(PANEL_RESULTS).unwrap();

        // Float browser.
        host.apply_mutations(
            vec![Mutation::Move {
                instance: browser.clone(),
                destination: Destination::NewSection {
                    section: SectionId::new("sec-browser-float").unwrap(),
                    index: 0,
                    location: SectionLocation::Floating { origin: Point { x: 40.0, y: 120.0 } },
                    size: Size { width: 280.0, height: 480.0 },
                },
            }],
            true,
        )
        .unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&browser).unwrap().location.zone(), Zone::Floating);

        // Redock browser to Right.
        host.apply_mutations(
            vec![Mutation::Move {
                instance: browser.clone(),
                destination: Destination::NewSection {
                    section: SectionId::new("sec-browser-right").unwrap(),
                    index: 0,
                    location: SectionLocation::Dock { zone: Zone::Right },
                    size: Size { width: 290.0, height: 480.0 },
                },
            }],
            true,
        )
        .unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&browser).unwrap().location.zone(), Zone::Right);

        // Join inspector and results as tabs in the same section.
        host.shared.borrow_mut().results_available = true;
        host.shared.borrow_mut().inspector_available = true;
        host.reconcile_runtime().unwrap();
        let results_sec = host.runtime.layout().layout().containing(&results).unwrap().id.clone();
        host.apply_mutations(
            vec![Mutation::Move { instance: inspector.clone(), destination: Destination::JoinTab { section: results_sec.clone(), index: 0 } }],
            true,
        )
        .unwrap();
        assert_eq!(host.runtime.layout().layout().containing(&inspector).unwrap().tabs.len(), 2);

        // Reset restores defaults and dirties store.
        host.reset_layout().unwrap();
        assert!(host.store.is_dirty());
        assert_eq!(host.runtime.layout().layout().containing(&browser).unwrap().location.zone(), Zone::Left);
    }

    #[test]
    fn clear_inspector_top_pops_one_item_at_a_time_and_keeps_lower_items() {
        let mut app = ModelApp::new(Session::default());
        app.run("model.material.create", json!({"name": "Steel", "elastic": {"type": "isotropic", "young": 210e9, "poisson": 0.3}})).unwrap();
        app.select_definition(Some((crate::browser::Def::Material, "Steel".into())));
        app.selected_mesh = true;
        app.reply =
            Some(crate::reply::Reply { command: "model.info".into(), label: "Model Info".into(), params: json!({}), value: json!({"nodes": 1}) });
        assert!(crate::context::inspector_active(&app));
        assert_eq!(crate::context::inspector_region_key(&app), Some("context:reply"));
        crate::context::clear_inspector_top(&mut app);
        assert!(app.reply.is_none());
        assert_eq!(crate::context::inspector_region_key(&app), Some("context:properties"));
        crate::context::clear_inspector_top(&mut app);
        assert!(app.selected.is_none());
        assert!(app.selected_mesh);
        assert_eq!(crate::context::inspector_region_key(&app), Some("context:mesh"));
        crate::context::clear_inspector_top(&mut app);
        assert!(!app.selected_mesh);
        assert!(!crate::context::inspector_active(&app));
    }

    #[test]
    fn native_windows_false_hides_popout_capability() {
        let mut app = ModelApp::new(Session::default());
        let host = app.panels.as_ref().unwrap();
        let panels = host.inspect_panels(&app);
        assert_eq!(panels.len(), 4);
        assert!(panels.iter().all(|p| p["native_supported"] == json!(false)));
        // Force an unavailable panel open check does not crash.
        app.selected_mesh = true;
        let panels = host.inspect_panels(&app);
        let insp = panels.iter().find(|p| p["id"] == "inspector").unwrap();
        assert_eq!(insp["has_content"], json!(true));
    }

    #[test]
    fn prepare_refreshes_results_and_inspector_availability_before_reopening() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        let ctx = egui::Context::default();
        let draw = |app: &mut ModelApp| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1280.0, 800.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        };

        // Initial frame has no results and no inspector item.
        draw(&mut app);
        app.with_panels(|_, panels| {
            let _ = panels.set_panel_visible(PANEL_RESULTS, false);
            let _ = panels.set_panel_visible(PANEL_INSPECTOR, false);
        });
        assert!(!app.panels.as_ref().unwrap().is_open(PANEL_RESULTS));
        assert!(!app.panels.as_ref().unwrap().is_open(PANEL_INSPECTOR));

        // While content is still absent, prepare leaves both closed.
        draw(&mut app);
        assert!(!app.panels.as_ref().unwrap().is_open(PANEL_RESULTS));
        assert!(!app.panels.as_ref().unwrap().is_open(PANEL_INSPECTOR));

        // Supply results and an inspector selection: prepare reopens both on that same frame.
        app.run("solve.static", json!({"step": "Load"})).unwrap();
        app.selected_mesh = true;
        draw(&mut app);
        assert!(app.panels.as_ref().unwrap().is_open(PANEL_RESULTS));
        assert!(app.panels.as_ref().unwrap().is_open(PANEL_INSPECTOR));
    }

    #[test]
    fn set_history_dock_updates_placement_while_history_is_hidden() {
        let mut host = PanelHost::from_prefs(&Prefs::default(), false).unwrap();
        let history = InstanceId::new(PANEL_HISTORY).unwrap();
        host.apply_mutations(
            vec![Mutation::Move {
                instance: history,
                destination: Destination::NewSection {
                    section: SectionId::new("sec-history-float").unwrap(),
                    index: 0,
                    location: SectionLocation::Floating { origin: Point { x: 80.0, y: 120.0 } },
                    size: Size { width: 520.0, height: 180.0 },
                },
            }],
            true,
        )
        .unwrap();
        assert!(!host.history_in_bottom_zone());

        host.set_panel_visible(PANEL_HISTORY, false).unwrap();
        assert!(!host.is_open(PANEL_HISTORY));
        host.set_history_dock(HistoryDock::Top).unwrap();
        host.set_panel_visible(PANEL_HISTORY, true).unwrap();

        assert!(host.history_in_bottom_zone());
        assert_eq!(host.history_dock(), HistoryDock::Top);
    }

    #[test]
    fn start_page_reports_panels_hidden_with_no_rect_without_losing_requested_open() {
        let mut app = ModelApp::new(Session::default());
        let ctx = egui::Context::default();
        let draw = |app: &mut ModelApp| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1280.0, 800.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        };

        draw(&mut app);
        let panels = app.panels.as_ref().unwrap().inspect_panels(&app);
        let browser = panels.iter().find(|p| p["id"] == "browser").unwrap();
        assert_eq!(browser["visible"], json!(true));
        assert!(browser["rect"].is_array());

        app.open_start();
        draw(&mut app);
        let panels = app.panels.as_ref().unwrap().inspect_panels(&app);
        let browser = panels.iter().find(|p| p["id"] == "browser").unwrap();
        assert_eq!(browser["visible"], json!(false));
        assert_eq!(browser["host"], json!("hidden"));
        assert_eq!(browser["rect"], Value::Null);
        assert!(app.panels.as_ref().unwrap().is_open(PANEL_BROWSER));

        app.close_start();
        draw(&mut app);
        let panels = app.panels.as_ref().unwrap().inspect_panels(&app);
        let browser = panels.iter().find(|p| p["id"] == "browser").unwrap();
        assert_eq!(browser["visible"], json!(true));
        assert!(browser["rect"].is_array());
    }

    #[test]
    fn view_menu_visibility_transitions_persist_across_reload() {
        let dir = std::env::temp_dir().join(format!("fc-view-menu-vis-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("prefs.json");
        let services = || crate::files::Services { prefs: Some(crate::files::PrefsStore::file(file.clone())), ..Default::default() };

        let mut app = ModelApp::new(Session::default()).with_services(services());
        let ctx = egui::Context::default();
        let draw = |app: &mut ModelApp| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1280.0, 800.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        };
        draw(&mut app);

        // Simulate View menu toggle directly changing UiState.
        app.ui.show_tree = false;
        app.ui.show_timeline = false;
        draw(&mut app);

        let reloaded = ModelApp::new(Session::default()).with_services(services());
        assert!(!reloaded.ui.show_tree);
        assert!(!reloaded.ui.show_timeline);
        assert!(!reloaded.panels.as_ref().unwrap().is_open(PANEL_BROWSER));
        assert!(!reloaded.panels.as_ref().unwrap().is_open(PANEL_HISTORY));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn hidden_browser_edits_apply_when_rejoining_shared_floating_section() {
        let mut host = PanelHost::from_prefs(&Prefs::default(), false).unwrap();
        let browser = InstanceId::new(PANEL_BROWSER).unwrap();
        let inspector = InstanceId::new(PANEL_INSPECTOR).unwrap();

        host.shared.borrow_mut().inspector_available = true;
        host.reconcile_runtime().unwrap();

        // Put Inspector in a 310-pt expanded floating section, then join Browser as a tab.
        host.apply_mutations(
            vec![Mutation::Move {
                instance: inspector.clone(),
                destination: Destination::NewSection {
                    section: SectionId::new("sec-shared-float").unwrap(),
                    index: 0,
                    location: SectionLocation::Floating { origin: Point { x: 60.0, y: 100.0 } },
                    size: Size { width: 310.0, height: 420.0 },
                },
            }],
            true,
        )
        .unwrap();
        let shared_sec = host.runtime.layout().layout().containing(&inspector).unwrap().id.clone();
        host.apply_mutations(
            vec![Mutation::Move { instance: browser.clone(), destination: Destination::JoinTab { section: shared_sec.clone(), index: 0 } }],
            true,
        )
        .unwrap();

        // Hide and reopen without pending edits: existing shared section geometry stays unchanged.
        host.set_panel_visible(PANEL_BROWSER, false).unwrap();
        host.set_panel_visible(PANEL_BROWSER, true).unwrap();
        let sec = host.runtime.layout().layout().containing(&browser).unwrap();
        assert!((sec.size.width - 310.0).abs() < 1e-6);
        assert!(!sec.collapsed);

        // Hide Browser, edit width/collapsed (rejecting 170 since Inspector requires 180), then reopen.
        host.set_panel_visible(PANEL_BROWSER, false).unwrap();
        assert!(host.set_browser_width(170.0).is_err());
        host.set_browser_width(400.0).unwrap();
        host.set_browser_collapsed(true).unwrap();
        host.set_panel_visible(PANEL_BROWSER, true).unwrap();

        let sec = host.runtime.layout().layout().containing(&browser).unwrap();
        assert_eq!(sec.id, shared_sec);
        assert!((sec.size.width - 400.0).abs() < 1e-6);
        assert!(sec.collapsed);
    }

    #[test]
    fn inspect_panels_reports_only_presented_tab_visible_when_results_and_inspector_share_section() {
        let mut app = ModelApp::new(Session::default());
        app.open_sample("cantilever");
        app.wait_loaded();
        app.run("solve.static", json!({"step": "Load"})).unwrap();
        app.selected_mesh = true;

        let ctx = egui::Context::default();
        let draw = |app: &mut ModelApp| {
            let raw = egui::RawInput { screen_rect: Some(Rect::from_min_size(egui::Pos2::ZERO, vec2(1280.0, 800.0))), ..Default::default() };
            let mut out = ctx.run_ui(raw, |ui| app.ui(ui));
            out.textures_delta.clear();
        };
        draw(&mut app);

        let results = InstanceId::new(PANEL_RESULTS).unwrap();
        let inspector = InstanceId::new(PANEL_INSPECTOR).unwrap();
        app.with_panels(|_, panels| {
            let insp_sec = panels.runtime.layout().layout().containing(&inspector).unwrap().id.clone();
            panels
                .apply_mutations(
                    vec![
                        Mutation::Move { instance: results.clone(), destination: Destination::JoinTab { section: insp_sec, index: 1 } },
                        Mutation::Activate { instance: inspector.clone() },
                    ],
                    true,
                )
                .unwrap();
        });

        draw(&mut app);

        let panels = app.panels.as_ref().unwrap().inspect_panels(&app);
        let insp_entry = panels.iter().find(|p| p["id"] == "inspector").unwrap();
        let res_entry = panels.iter().find(|p| p["id"] == "results").unwrap();
        assert_eq!(insp_entry["visible"], json!(true));
        assert_eq!(insp_entry["has_content"], json!(true));
        assert_eq!(insp_entry["host"], json!("docked"));
        assert!(insp_entry["rect"].is_array());

        assert_eq!(res_entry["visible"], json!(false));
        assert_eq!(res_entry["has_content"], json!(true));
        assert_eq!(res_entry["host"], json!("hidden"));
        assert_eq!(res_entry["rect"], Value::Null);
    }
}
